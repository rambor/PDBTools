require "PDB/version"
require "PDB/atom"
require "PDB/molecule"
require "gsl"

module PDB
  # Model contains many molecules derived from PDB File
  # define global variables
  $twoPI = Math::PI*2.0
  $letters = %{A B C D E F G H I J K L M N O P Q R S T U V W X Y Z}

  class Model

    # active_set array of atoms that are in use
    #
    # :molecule => holds a hash, key is chain and value is Molecule class
    #
    attr_reader :filename, :molecule, :active_set, :molecules, :dmax, :extrema, :centering_coordinates, :random_extrema

    def initialize(file, options={})
      @filename = file
      @active_set =[]
      @molecules = Hash.new # use chain as molecule key for the hash
      raise ArgumentError.new("File does not exists in directory #{file}") unless File.exists?(file)
      @extrema=Array.new(14)
      @abcdefgh=Array.new(8)
      @centering_coordinates=Array.new(3)

      reset_extremes
      openPDBFile(file, options)
    end


    # Open PDB file and create Molecule
    def openPDBFile(filename, flags={})

      puts "OPENING #{filename}"
      PDB::report_log("OPENING #{filename}")
      pdbLines=[]

      open(filename){|x| pdbLines = x.readlines}

      pdbLines.select!{ |x| x =~ /^ATOM/}
      pdbLines.collect!{|x| PDB::Atom.new(x) }

      # group atoms by chains organized by collection of residues
      current_chain = pdbLines[0].chain
      @molecules[current_chain.to_sym] = PDB::Molecule.new

      pdbLines.each do |atom|

        if (atom.chain != current_chain)
          current_chain = atom.chain
          @molecules[current_chain.to_sym] = PDB::Molecule.new
        end

        # strip out waters
        if atom.residue == "HOH"
          @molecules[:HOH] ||= PDB::Molecule.new
          atom.chain = "HOH" # set chain to WAT
          @molecules[:HOH].addAtom(atom)
        else
          @molecules[current_chain.to_sym].addAtom(atom)
        end
      end


      # create working set based on selection (flags)
      if flags.size == 0
        fillActiveSet
      else
        @active_set.clear
        # :CA => false
        # :CA => true
        selectAtomsByAttribute(flags)
      end

    end



    # Select Atoms by Attribute making a copy of the atom, only selects as AND statements (or intersections)
    #
    # :atom_type => {:CA => true}
    # usage : selectAtomsByAttributes( :atom_type => {:CA => true, :CB => true})
    # usage : selectAtomsByAttributes({:atom_type => {:CA => true, :CB => true}, :chain => {:A => true}})
    # usage : selectAtomsByAttributes({:atom_type => {:CA => true, :CB => true}, :chain => {:B => false}})
    # selectAtomsByAttributes(
    #   {:atom_type => {:CA => true, :CB => true}, :chain => {:B => false}},
    #   {:atom_type => {:O => true}, :residue => {:HOH => true}},
    # )
    # Different Hashes at top level are AND (intersection)
    #
    # selectAtomsByAttribute(:atom_type => {:CB => true}, :residue => {:ASP => true})
    #
    # To select as a compound AND statement, use an array of hashes, each element of array is a hash
    #
    # usage : selectAtomsByAttribute([{:atom_type => {:CB => true}, :residue => {:ASP => true}},{:residue => {:HOH => true}}])
    #
    def selectAtomsByAttribute(selection=[])

      @active_set.clear
      #
      # atom_type => {:CA => false, :CB => true}, :chain => {:A => false}
      # apply false to all selected items ???
      # each element of the array is an AND so, we concatenate the output
      if selection.class == Array
        selection.each do |hash_selection|
          @active_set.concat(select_from_molecules_true(hash_selection))
        end
      else
        @active_set.concat(select_from_molecules_true(selection))
      end

      remove_duplicate_atom_selections(@active_set)

#
      # if active_set is empty, assume we select all
      if @active_set.size == 0
        puts "Filling? #{@active_set.size}"
        @active_set.clear
        fillActiveSet
      end

      #
      # remove any unwanted attributes => false flags
      #
      if selection.class == Array
        # hash_selection is something like:
        # atom_type => {:CA => false, :CB => true}
        # selection = [{:atom_type => {:CA => true, :CB => true}}, {:residue => {:HOH => true}}]
        selection.each do |hash_selection|
          # attribute = :atom_type
          select_from_active_set_false(@active_set, hash_selection)
        end
      else
        select_from_active_set_false(@active_set, selection)
      end

      setExtrema
    end


    def select_from_active_set_false(array_of_atoms, selection_hash)
      #
      # attribute == :atom_type, :chain, etc
      # keys must be attributes in Atom class
      #
      temp=[]
      selection_hash.each_pair do |attribute, parameters|
        # iterate over each molecule, concatenate for each time it is true
        raise ArgumentError.new("Property not found in Atom Class #{attribute}") unless PDB::Atom.method_defined? attribute
          #
          # type_symbol = :CA, :CB, :CHAIN_IDENTIFIER
          # value = true or false
          parameters.each_pair do |type_symbol, value|

            type = (type_symbol.is_a? Symbol) ? type_symbol.to_s.upcase : type_symbol.upcase

            if !value
              PDB::report_log("SELECTOR : REMOVING #{attribute} => #{type}")
              #temp.concat( array_of_atoms.select{|atom| atom.send(attribute) != type} )
              temp.concat( array_of_atoms.delete_if{|atom| atom.send(attribute) == type} )
            end
          end
      end
    end



    # Example : selectAtomsBy(:atom_type => {:CA => true, :CB => false}, :chain => {:B => false})
    def select_from_molecules_true(selection_hash)
      keys = selection_hash.keys
      attribute = keys[0]

      raise ArgumentError.new("Property not found in Atom Class #{attribute}") unless PDB::Atom.method_defined? attribute

      temp=Array.new
      # create initial selection by selecting over all molecules in Model

      selection_hash[attribute].each_pair do |type_symbol, value|
        @molecules.each_value do |mol|
          type = type_symbol.to_sym
          if value
            PDB::report_log("SELECTOR : ACCEPTING #{attribute} => #{type} : CHAIN => #{mol.chain}")
            temp.concat( mol.selectAtomsBy(attribute, type) )
          end
        end
      end

      # all subsequent selections are applied to temp
      selected_set=[]
      if keys.size > 1
        for i in 1...keys.size
          attribute = keys[i]
          # type_symbol = :CA, :CB, etc
          selection_hash[attribute].each_pair do |type_symbol, value|
            type = type_symbol.to_sym
            puts "SELECTING #{attribute} #{type} <=> #{value} #{selected_set.size}"
            if value
              PDB::report_log("SELECTOR : ACCEPTING #{attribute} => #{type}")
              selected_set.concat(temp.select{|atom| atom.send(attribute) == (type.to_s).upcase})
            end
          end
        end
      else
        selected_set = temp
      end

      # remove duplicates
      remove_duplicate_atom_selections(selected_set)
      return selected_set
    end


    #
    # use this method for removing duplicate entries that may occur during a selection
    # :chain => {:A => true}
    # :atom_type => {:CA = true}
    #
    def remove_duplicate_atom_selections(array_of_atoms)
      array_of_atoms.uniq!{|atom| atom.atom_number || atom.chain }
    end


    #
    # Add all atoms to active_set
    #
    def fillActiveSet
      @active_set.clear
      @molecules.each_value do |mol|
        mol.residues.each do |res|
          res.atoms.each do |atom|
            @active_set << atom.dup
          end
        end
      end
      # reset extreme values
      setExtrema
    end



    def checkIfExtrema(atom)

      if !atom.atom_type.include?("H")

        x = atom.xpos
        y = atom.ypos
        z = atom.zpos

        if (x > @extrema[0].xpos)
          @extrema[0] = atom
        elsif (x < @extrema[1].xpos)
          @extrema[1] = atom
        elsif (y > @extrema[2].ypos)
          @extrema[2] = atom
        elsif (y < @extrema[3].ypos)
          @extrema[3] = atom
        elsif (z > @extrema[4].zpos)
          @extrema[4] = atom
        elsif (z < @extrema[5].zpos)
          @extrema[5] = atom
        end

        # maximizes
        if (x-y+z) > @abcdefgh[0]
          @extrema[6] = atom
          @abcdefgh[0] = (x-y+z)
        end

        if (x-y-z) > @abcdefgh[1]
          @extrema[7] = atom
          @abcdefgh[1] = (x-y-z)
        end

        if (x+y+z) > @abcdefgh[2]
          @extrema[8] = atom
          @abcdefgh[2] = (x+y+z)
        end

        if (x+y-z) > @abcdefgh[3]
          @extrema[9] = atom
          @abcdefgh[3] = (x+y-z)
        end

        # minimizes
        if (x-y+z) < @abcdefgh[4]
          @extrema[10] = atom
          @abcdefgh[4] = (x-y+z)
        end

        if (x-y-z) < @abcdefgh[5]
          @extrema[11] = atom
          @abcdefgh[5] = (x-y-z)
        end

        if (x+y+z) < @abcdefgh[6]
          @extrema[12] = atom
          @abcdefgh[6] = (x+y+z)
        end

        if (x+y-z) < @abcdefgh[7]
          @extrema[13] = atom
          @abcdefgh[7] = (x+y-z)
        end

      end
    end


    def setExtrema()
      # output extreme values
      # only consider non-hydrgen atoms
      reset_extremes
      @active_set.each do |atom|
        checkIfExtrema(atom)
      end
    end


    def setDmax
      max = 0
      for i in 0...@total
        atom1 = @active_set[i]
        nextAtom = i + 1
        for j in nextAtom...@total
          atom2 = @active_set[j]
          dx = atom1.xpos - atom2.xpos
          dy = atom1.ypos - atom2.ypos
          dz = atom1.zpos - atom2.zpos

          dis2 = dx*dx + dy*dy + dz*dz
          if (dis2 > max)
            max = dis2
            @extrema[6] = @active_set[i]
            @extrema[7] = @active_set[j]
          end
        end
      end

      @dmax = Math::sqrt(max)
      PDB::report_log("DMAX => #{@dmax}")
    end



    def writeModelToFile(filename)
      writeToFile(@active_set, filename)
    end



    def writeExtremaToFile(filename)
      writeToFile(@extrema, filename)
    end


    def writeRandomExtremaToFile(filename)
      writeToFile(@molecule.random_extrema, filename)
    end


    # Write PDB file from array of PDBAtom objects
    def writeToFile(atoms, filename)

      if filename.split(/\./).last != "pdb"
        filename = filename + ".pdb"
      end

      filename.tr!(" ", "_")

      newLines =[]
      count = 1
      current_chain = atoms[0].chain

      atoms.each do |atom|

        # if next chain is different, add ter
        if (current_chain != atom.chain)
          newLines << "TER\n"
          current_chain = atom.chain
        end

        if atom.residue.size == 1
          newLines << sprintf("ATOM  %5s %-4s %3s %1s%4s    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n", count, atom.atom_type, convert_to_3_letter(atom.residue, atom.atom_type), atom.chain, atom.resid, atom.xpos, atom.ypos, atom.zpos, atom.occ, atom.temp, atom.atom)
        elsif atom.residue.size == 3
          newLines << sprintf("ATOM  %5s %-4s %3s %1s%4s    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n", count, atom.atom_type, atom.residue, atom.chain, atom.resid, atom.xpos, atom.ypos, atom.zpos, atom.occ, atom.temp, atom.atom)
        else
          abort("Improper residue in PDB file -> can't write to file")
        end

        count += 1
      end

      newLines << "END\n"
      open(filename, 'w'){|x| newLines.each {|xx| x << xx}}
    end


    #
    # geometric methods for translating and rotating
    #
    def translate(delx, dely, delz)

      @active_set.each do |atom|
        atom.xpos = atom.xpos + delx
        atom.ypos = atom.ypos + dely
        atom.zpos = atom.zpos + delz
      end
    end


    # Rotate Molecule by specific angles in radians
    def rotate(euler_x, euler_y, euler_z)
      rotate_x(euler_x)
      rotate_y(euler_y)
      rotate_z(euler_z)
    end


    # Rotates molecule randomly
    def random_rotate()
      PDB::report_log("RANDOM ROTATION")
      angle = rand()*$twoPI
      rotate_x(angle)

      angle = rand()*$twoPI
      rotate_y(angle)

      angle = rand()*$twoPI
      rotate_z(angle)
      PDB::report_log("RANDOM ROTATION ENDED")
    end


    # rotate along x
    def rotate_x(angle)
      newCoords = GSL::Vector.alloc(3)
      xrotation = GSL::Matrix.alloc(3,3)
      cos = Math::cos(angle)
      sin = Math::sin(angle)

      xrotation[0,0] = 1
      xrotation[1,1] = cos
      xrotation[1,2] = -1.0*sin
      xrotation[2,1] = sin
      xrotation[2,2] = cos

      @active_set.each do |atom|
        newCoords[0] = atom.xpos
        newCoords[1] = atom.ypos
        newCoords[2] = atom.zpos
        newCoords = xrotation*newCoords
        atom.xpos = newCoords[0]
        atom.ypos = newCoords[1]
        atom.zpos = newCoords[2]
      end
      PDB::report_log("ROTATED X-AXIS => #{angle}")
    end


    # rotate along y
    def rotate_y(angle)
      newCoords = GSL::Vector.alloc(3)
      yrotation = GSL::Matrix.alloc(3,3)
      cos = Math::cos(angle)
      sin = Math::sin(angle)

      yrotation[0,0] = cos
      yrotation[0,2] = sin
      yrotation[1,1] = 1
      yrotation[2,0] = -sin
      yrotation[2,2] = cos

      @active_set.each do |atom|
        newCoords[0] = atom.xpos
        newCoords[1] = atom.ypos
        newCoords[2] = atom.zpos
        newCoords = yrotation*newCoords
        atom.xpos = newCoords[0]
        atom.ypos = newCoords[1]
        atom.zpos = newCoords[2]
      end
      PDB::report_log("ROTATED Y-AXIS => #{angle}")
    end


    # rotate along z
    def rotate_z(angle)
      newCoords = GSL::Vector.alloc(3)
      zrotation = GSL::Matrix.alloc(3,3)
      cos = Math::cos(angle)
      sin = Math::sin(angle)

      zrotation[0,0] = cos
      zrotation[0,2] = sin
      zrotation[1,1] = 1
      zrotation[2,0] = -sin
      zrotation[2,2] = cos

      @active_set.each do |atom|
        newCoords[0] = atom.xpos
        newCoords[1] = atom.ypos
        newCoords[2] = atom.zpos
        newCoords = zrotation*newCoords
        atom.xpos = newCoords[0]
        atom.ypos = newCoords[1]
        atom.zpos = newCoords[2]
      end
      PDB::report_log("ROTATED Z-AXIS => #{angle}")
    end



    def center_molecule
      calculateCenteringCoordinates

      @active_set.each do |atom|
        atom.xpos = atom.xpos - @centering_coordinates[0]
        atom.ypos = atom.ypos - @centering_coordinates[1]
        atom.zpos = atom.zpos - @centering_coordinates[2]
      end
    end


    #
    # select random extrema points
    # first three points maximize area of a triangle
    # fourth point maximizes distance from a point ot the triangular plane
    # additional points are taken at random
    def select_random_extrema(total_to_select)

      begin
        if total_to_select < 3 || total_to_select > 12
          raise "Total number of random points must be > 3 and < 12"
        end
      rescue => error
        PDB::report_error("#{error.class} and #{error.message} : #{total_to_select} ")
      end

      @extrema.shuffle!

      points=Array.new(total_to_select)
      points[0] = @extrema[0].dup
      # select longest distance
      prior = 0
      for i in 1...@extrema.size

        anchor = @extrema[i]

        delx = points[0].xpos - anchor.xpos
        dely = points[0].ypos - anchor.ypos
        delz = points[0].zpos - anchor.zpos

        dis2 = delx*delx + dely*dely + delz*delz
        if dis2 > prior
          points[1] = @extrema[i].dup
          prior = dis2
        end
      end

      prior = 0
      # find the next point that maximizes distance between first two select points
      for i in 1...@extrema.size

        anchor = @extrema[i]

        delx = points[0].xpos - anchor.xpos
        dely = points[0].ypos - anchor.ypos
        delz = points[0].zpos - anchor.zpos

        delx2 = points[1].xpos - anchor.xpos
        dely2 = points[1].ypos - anchor.ypos
        delz2 = points[1].zpos - anchor.zpos

        dis2 = Math::sqrt(delx*delx + dely*dely + delz*delz) + Math::sqrt(delx2*delx2 + dely2*dely2 + delz2*delz2)

        if dis2 > prior
          points[2] = @extrema[i].dup
          prior = dis2
        end
      end

      # find fourth point as a non-planar point (distances between point and plane is maximized)
      vec1 = GSL::Vector[points[0].xpos - points[2].xpos, points[0].ypos - points[2].ypos, points[0].zpos - points[2].zpos]
      vec2 = GSL::Vector[points[1].xpos - points[2].xpos, points[1].ypos - points[2].ypos, points[1].zpos - points[2].zpos]

      # cross-product - get parameters for a plane
      a_coeff = vec1[1]*vec2[2] - vec1[2]*vec2[1]
      b_coeff = vec1[2]*vec2[0] - vec1[0]*vec2[2]
      c_coeff = vec1[0]*vec2[1] - vec1[1]*vec2[0]
      d_coeff = a_coeff*points[2].xpos + b_coeff*points[2].ypos + c_coeff*points[2].zpos

      inv_coeff = 1.0/Math::sqrt(a_coeff*a_coeff + b_coeff*b_coeff + c_coeff*c_coeff)

      prior = 0
      # find the next point that maximizes distance between first two select points
      for i in 1...@extrema.size

        anchor = @extrema[i]
        dis = (anchor.xpos*a_coeff + anchor.ypos*b_coeff + anchor.zpos*c_coeff - d_coeff).abs * inv_coeff

        if (dis > prior)
          points[3] = @extrema[i].dup
          prior = dis
        end
        #end
      end

      # add additional points from random selection
      if total_to_select > 4
        for i in 4...total_to_select
            for j in 0...14
              temp = points.select{|p| p.xpos == @extrema[j].xpos}
              if temp.size == 0
                @points[i] = @extrema[j].dup
                break
              end
            end
        end
      end

      @random_extrema = points
    end


    def reset_extremes
      temp = "ATOM      1  N   ASP A   1       0.000   0.000   0.000  0.00  0.00           N"
      @extrema[0] = Atom.new(temp) # holds max X
      @extrema[0].xpos = -100000
      @extrema[1] = Atom.new(temp) # holds min X
      @extrema[1].xpos = 100000
      @extrema[2] = Atom.new(temp) # holds max Y
      @extrema[2].ypos = -100000
      @extrema[3] = Atom.new(temp) # holds min Y
      @extrema[3].ypos = 100000
      @extrema[4] = Atom.new(temp) # holds max Z
      @extrema[4].zpos = -100000
      @extrema[5] = Atom.new(temp) # holds min Z
      @extrema[5].zpos = 100000


      @abcdefgh[0] = -1000000  # 6
      @abcdefgh[1] = -1000000  # 7
      @abcdefgh[2] = -1000000  # 8
      @abcdefgh[3] = -1000000  # 9
      @abcdefgh[4] = 1000000   # 10
      @abcdefgh[5] = 1000000   # 11
      @abcdefgh[6] = 1000000   # 12
      @abcdefgh[7] = 1000000   # 13
      PDB::report_log("EXTREMA POINTS RESET")
    end
    private :reset_extremes


    def calculateCenteringCoordinates
      x=0
      y=0
      z=0
      invCount = 1.0/@total

      @active_set.each do |atom|
        x += atom.xpos
        y += atom.ypos
        z += atom.zpos
      end

      @centering_coordinates[0] = x*invCount
      @centering_coordinates[1] = y*invCount
      @centering_coordinates[2] = z*invCount
    end
    private :calculateCenteringCoordinates

  end # end of Model Class


  #
  # == PDB Module methods defined below
  #

  # Extract sequence of residues from array of atoms.
  # Order of the sequence is set by the order of the atoms in the array.
  #
  # ==== Attributes
  #
  # * +atoms+ - Array of Atoms
  #
  def extract_sequence(atoms)

    sequence={}
    first = atoms[0].resid
    sequence[atoms[0].chain] =[]
    sequence[atoms[0].chain] << atoms[0].residue

    atoms.each do |atom|
      #if first != atom.resid
      #  @sequence << atom.resid
      #  first = atom.resid
      #end
      if (sequence.key?(atom.chain))

        if first != atom.resid
          sequence[atom.chain] << atom.residue
          first = atom.resid
        end
      else
        sequence[atom.chain] =[]
        first = atom.resid
        sequence[atom.chain] << atom.residue
      end
    end
  end
  module_function :extract_sequence


  # Convert 3-letter residue to 1-letter
  # search a hash for the 3 letter input string
  #
  # ==== Attributes
  #
  # * +residue+ - 3 letter string
  #
  # @param residue string
  def convert_to_one_letter(residue)

    res = {
        :GUA => "G",
        :ADE => "A",
        :CYT => "C",
        :URI => "U",
        :THY => "T",
        :Gr  => "G",
        :Ar  => "A",
        :Cr  => "C",
        :Ur  => "U",
        :ALA => "A",
        :ARG => "R",
        :ASN => "N",
        :ASP => "D",
        :ASX => "B",
        :CYS => "C",
        :GLU => "E",
        :GLN => "Q",
        :GLX => "Z",
        :GLY => "G",
        :HIS => "H",
        :ILE => "I",
        :LEU => "L",
        :LYS => "K",
        :MET => "M",
        :PHE => "F",
        :PRO => "P",
        :SER => "S",
        :THR => "T",
        :TRP => "W",
        :TYR => "Y",
        :VAL => "V",
        :SEC => "U",
        :PCA => "J"
    }

    residue = residue.upcase.to_sym
    return res[residue]
  end
  module_function :convert_to_one_letter


  def report_log(message)
    (Thread.current[:messages] ||= []) << "#{Time.now}\t#{message}"

    File.open('PDBTools.log', 'a') do |file|
      file.puts "#{Time.now}  #{message}"
    end
  end
  module_function :report_log


  def report_error(error_message)
    (Thread.current[:errors] ||= []) << "#{Time.now}\t#{error_message}"

    File.open('PDBTools_errors.txt', 'a') do |file|
      file.puts  "#{Time.now}  #{error_message}"
     # (Thread.current[:errors] ||= []).each do |error|
     #   file.puts error
     # end
    end
  end
  module_function :report_error

end



