
require "PDB/atom"

module PDB
  class Molecule

    attr_reader :atoms, :total, :extrema, :dmax, :centering_coordindates, :sequence

    def initialize()
      @total=0
      @atoms=[]
      @sequence=[]

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


    def extractSequence

      first = @atoms[0].resid
      @sequence << first
      atoms.each do |atom|
        if first != atom.resid
          @sequence << atom.resid
          first = atom.resid
        end
      end
    end


    def getAtomByIndex(index)
      return @atoms[index]
    end


    def removeAtomsByType(atom_type)
      @atoms.select!{|atom| atom.atomType != atom_type}
      @total = @atoms.size
    end


    def selectAtomsByType(atom_types={})

      temp=[]
      atom_types.each_pair do |key, value|
        if (value)
          temp = temp + @atoms.select{|atom| atom.atomType == (key.to_s).upcase}
        end
      end

      @atoms = temp
      self.setExtrema
      @total = atoms.size
    end


    def selectAtomsByChains(chains={})

      temp=[]
      chains.each_pair do |key, value|
        if (value)
          temp = temp + @atoms.select{|atom| atom.chain == (key.to_s).upcase}
        end
      end

      @atoms = temp
      self.setExtrema
      @total = atoms.size
    end



    def checkIfExtrema(atom)

      if !atom.atomType.include?("H")

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


    def center_molecule
      self.calculateCenteringCoordinates

      @atoms.each do |atom|
        atom.xpos = atom.xpos - @centering_coordinates[0]
        atom.ypos = atom.ypos - @centering_coordinates[1]
        atom.zpos = atom.zpos - @centering_coordinates[2]
      end
    end

    def self.calculateCenteringCoordinates
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
    private_class_method :calculateCenteringCoordinates

  end # end of class definition




end
