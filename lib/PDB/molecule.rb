
require "PDB/atom"
require "PDB/residue"
require 'linked-list'
require 'PDB/fasta'

module PDB
  #
  # == Molecule
  # Molecule Class holds the collection of Atoms that were derived from PDB file
  # Extreme values of the coordinates are determined for each atom added or after selection
  # Methods are available for various atom selections attributes of the Atom class
  #
  class Molecule


    attr_reader :index_of_last_residue_added, :chain_breaks, :mass, :total, :extrema, :dmax, :centering_coordindates, :sequence, :resids, :random_extrema, :residues, :chain

    # The +new+ class method initializes the class.
    # === Parameters
    # * _non required_
    #
    def initialize()
      @total=0
      @mass=0
      @sequence=[] # use chain ID as key and sequence array as value
      @resids=[]   # use chain ID as key and sequence array as value
      @chain_breaks=[]

      @residues = LinkedList::List.new
      # set the dummy extreme atom
      @extrema=Array.new(14)
      @abcdefgh=Array.new(8)

      reset_extremes
      @centering_coordinates=Array.new(3)

    end



    def addAtom(atom)
      # handle exception
      begin

        @chain ||= atom.chain
        raise ("CHAIN MISMATCH FOR NEW ATOM #{@chain} != #{atom.chain} FOR ATOM #{atom.index}") unless atom.chain == @chain

        # get resid
        # if resid is not found, make new resid
        # molecule is managed as a collection of residues
        if !@residues.last.nil? && @residues.last.resid == atom.resid
          @residues.last.add_atom(atom)
        else
          current_residue = checkForResid(atom.resid)
          if current_residue.nil?
            @residues.push(Residue.new(atom))
          else
            current_residue.add_atom(atom)
          end
        end


        checkIfExtrema(atom)
        @total += 1
      rescue => error
        PDB::report_error("#{error.class} and #{error.message} : #{atom.inspect} ")
      end
    end



    def checkForResid(resid)
      if @residues.size == 0
        return nil
      else
        @residues.each do |residue|
          if residue.resid == resid
            return residue
          end
        end

        return nil
      end
    end



    def change_chain(new_chain)
      @residues.each do |residue|
        residue.atoms.each do |atom|
          atom.chain = new_chain
        end
      end
    end


    # Extract sequence from input PDB after loading atoms, can't be used in initializer
    # Sequences are extracted as an array of strings for each chain
    # Resids are also collected for each chain
    def extractSequence


      @index_of_last_residue_added = @residues.first.resid;
      #@sequence << first
      @residues.each do |res|
        @sequence << res.resname
        @resids << res.resid

        # detect discontinuities in chain
        if ( (res.resid - @index_of_last_residue_added) > 1)
          @chain_breaks << @index_of_last_residue_added
        end

        @index_of_last_residue_added = res.resid

        end
    end


    def getTotalChainBreaks
      return @chain_breaks.size
    end


    def getMass
      @mass = 0
      @residues.each do |res|
        res.atoms.each do |atom|
          @mass += atom.mass
        end
      end
      @mass
    end

    # read in a fasta file and compare sequence to sequence
    # PDB could be more or less than input FASTA
    def compareToFasta(filename)

        fasta = Fasta.new(filename)

      # residues may not start in same place
      # can FASTA sequence be smaller than PDB sequence?

    end


    # Select atoms using a hash container
    #
    # * *Args*
    #   - +hash_key+ -> must be an attribute in atom.rb (Atom Class)
    #   - +selection+ -> hash using symbols of the attribute
    #
    # @param hash_key [:symbol]
    # Example : selectAtomsBy(:atom_type, :CA)
    #

    def selectAtomsBy(attribute, type)

      temp=[]

      @residues.each do |res|
        res.atoms.each do |atom|
          if atom.send(attribute) == (type.to_s).upcase
            temp << atom.dup
          end
        end
      end

      #setExtrema
      temp
    end


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
      points[0] = @extrema[0]
      # select longest distance
      prior = 0
      for i in 1...@extrema.size

        anchor = @extrema[i]

        delx = points[0].xpos - anchor.xpos
        dely = points[0].ypos - anchor.ypos
        delz = points[0].zpos - anchor.zpos

        dis2 = delx*delx + dely*dely + delz*delz
        if dis2 > prior
          points[1] = @extrema[i]
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
            points[2] = @extrema[i]
            prior = dis2
          end
      end

      # find fourth point
      vec1 = GSL::Vector[points[0].xpos - points[2].xpos, points[0].ypos - points[2].ypos, points[0].zpos - points[2].zpos]
      vec2 = GSL::Vector[points[1].xpos - points[2].xpos, points[1].ypos - points[2].ypos, points[1].zpos - points[2].zpos]

      # cross-product
      a_coeff = vec1[1]*vec2[2] - vec1[2]*vec2[1]
      b_coeff = vec1[2]*vec2[0] - vec1[0]*vec2[2]
      c_coeff = vec1[0]*vec2[1] - vec1[1]*vec2[0]
      d_coeff = a_coeff*points[2].xpos + b_coeff*points[2].ypos + c_coeff*points[2].zpos

      inv_coeff = 1.0/Math::sqrt(a_coeff*a_coeff + b_coeff*b_coeff + c_coeff*c_coeff)

      prior = 0
      # find the next point that maximizes distance between first two select points
      for i in 1...@extrema.size
        anchor = @extrema[i]
        #if anchor.atom_number != points[0].atom_number && anchor.atom_number != points[1].atom_number && anchor.atom_number != points[2].atom_number

          dis = (anchor.xpos*a_coeff + anchor.ypos*b_coeff + anchor.zpos*c_coeff - d_coeff).abs * inv_coeff

          if (dis > prior)
            points[3] = @extrema[i]
            prior = dis
          end
        #end
      end

      # add additional points from random selection
      @random_extrema = points
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


  end # end of class definition


end
