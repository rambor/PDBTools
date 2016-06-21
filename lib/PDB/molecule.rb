
require "PDB/atom"

module PDB
  #
  # == Molecule
  # Molecule Class holds the collection of Atoms that were derived from PDB file
  # Extreme values of the coordinates are determined for each atom added or after selection
  # Methods are available for various atom selections attributes of the Atom class
  #
  class Molecule

    attr_reader :atoms, :total, :extrema, :dmax, :centering_coordindates, :sequence, :resids

    # The +new+ class method initializes the class.
    # === Parameters
    # * _non required_
    #
    def initialize()
      @total=0
      @atoms=[]
      @sequence={} # use chain ID as key and sequence array as value
      @resids={}   # use chain ID as key and sequence array as value

      # set the dummy extreme atom
      @extrema=Array.new(8)
      temp = "ATOM      1  N   ASP A   1       0.000   0.000   0.000  0.00  0.00           N"
      @extrema[0] = Atom.new(temp)
      @extrema[0].xpos = -10000 # maxX
      @extrema[1] = Atom.new(temp)
      @extrema[1].xpos = 10000  # minX
      @extrema[2] = Atom.new(temp)
      @extrema[2].ypos = -10000 # maxY
      @extrema[3] = Atom.new(temp)
      @extrema[3].ypos = 10000  # minY
      @extrema[4] = Atom.new(temp)
      @extrema[4].zpos = -10000 # maxZ
      @extrema[5] = Atom.new(temp)
      @extrema[5].zpos = 10000  # minZ

      @centering_coordinates=Array.new(3)
    end


    def addAtom(atom)
      # handle exception
      begin
        @atoms << atom
        checkIfExtrema(@atoms.last)
        @total += 1
      rescue => error
        PDB::report_error("#{error.class} and #{error.message} : #{atom.inspect} ")
      end
    end

    # Extract sequence from input PDB
    # Sequences are extracted as an array of strings for each chain
    # Resids are also collected for each chain
    def extractSequence

      first = @atoms[0].resid
      @sequence[@atoms[0].chain] =[]
      @resids[@atoms[0].chain] =[]
      @sequence[@atoms[0].chain] << @atoms[0].residue
      @resids[@atoms[0].chain] << first
      #@sequence << first
      atoms.each do |atom|
        #if first != atom.resid
        #  @sequence << atom.resid
        #  first = atom.resid
        #end
        if (@sequence.key?(atom.chain))

          if first != atom.resid
            @sequence[atom.chain] << atom.residue
            first = atom.resid
            @resids[atom.chain] << first
          end
        else

          @sequence[atom.chain] =[]
          first = atom.resid
          @sequence[atom.chain] << atom.residue
          @resids[atom.chain] << first
        end

      end
    end


    def getAtomByIndex(index)
      return @atoms[index]
    end


    def removeAtomsByType(atom_type)
      @atoms.select!{|atom| atom.atom_type != atom_type}
      @total = @atoms.size
    end



    def removeWaters
      @atoms.select!{|atom| atom.residue != "HOH" }
    end



    # Select atoms using a hash container
    #
    # * *Args*
    #   - +hash_key+ -> must be an attribute in atom.rb (Atom Class)
    #   - +selection+ -> hash using symbols of the attribute
    #
    # @param hash_key [:symbol]
    def selectAtomsByType(hash_key, selection={})

      temp=[]

      selection.each_pair do |type, value|
        if (value)
          temp = temp + @atoms.select{|atom| atom.send(hash_key) == (type.to_s).upcase }
        end
      end

      @atoms = temp

      selection.each_pair do |type, value|
        if (value == false)
          @atoms.select!{|atom| atom.send(hash_key) != (type.to_s).upcase }
        end
      end

      setExtrema
      @total = atoms.size
    end



    def checkIfExtrema(atom)

      if !atom.atom_type.include?("H")

        if (atom.xpos > @extrema[0].xpos)
          @extrema[0] = atom
        elsif (atom.xpos < @extrema[1].xpos)
          @extrema[1] = atom
        elsif (atom.ypos > @extrema[2].ypos)
          @extrema[2] = atom
        elsif (atom.ypos < @extrema[3].ypos)
          @extrema[3] = atom
        elsif (atom.zpos > @extrema[4].zpos)
          @extrema[4] = atom
        elsif (atom.zpos < @extrema[5].zpos)
          @extrema[5] = atom
        end
      end
    end



    def setExtrema()
      # output extreme values
      # only consider non-hydrgen atoms
      @atoms.each do |atom|
        checkIfExtrema(atom)
      end
    end



    def setDmax
      max = 0
      for i in 0...@total
        atom1 = @atoms[i]
        nextAtom = i + 1
        for j in nextAtom...@total
          atom2 = @atoms[j]
          dx = atom1.xpos - atom2.xpos
          dy = atom1.ypos - atom2.ypos
          dz = atom1.zpos - atom2.zpos

          dis2 = dx*dx + dy*dy + dz*dz
          if (dis2 > max)
            max = dis2
            @extrema[6] = @atoms[i]
            @extrema[7] = @atoms[j]
          end
        end
      end

      @dmax = Math::sqrt(max)
    end



    def translate(delx, dely, delz)

      @atoms.each do |atom|
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

      angle = rand()*$twoPI
      rotate_x(angle)

      angle = rand()*$twoPI
      rotate_y(angle)

      angle = rand()*$twoPI
      rotate_z(angle)
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

      @atoms.each do |atom|
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

      @atoms.each do |atom|
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

      @atoms.each do |atom|
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

      @atoms.each do |atom|
        atom.xpos = atom.xpos - @centering_coordinates[0]
        atom.ypos = atom.ypos - @centering_coordinates[1]
        atom.zpos = atom.zpos - @centering_coordinates[2]
      end
    end



    def calculateCenteringCoordinates
      x=0
      y=0
      z=0
      invCount = 1.0/@total

      @atoms.each do |atom|
        x += atom.xpos
        y += atom.ypos
        z += atom.zpos
      end

      @centering_coordinates[0] = x*invCount
      @centering_coordinates[1] = y*invCount
      @centering_coordinates[2] = z*invCount
    end
    private :calculateCenteringCoordinates


  end # end of class definition


end
