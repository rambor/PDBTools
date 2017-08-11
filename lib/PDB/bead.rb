module PDB
  class Bead

    attr_reader :radius, :xpos, :ypos, :zpos, :index

    def initialize(radius, x, y, z, index)
      @radius = radius
      @xpos = x
      @ypos = y
      @zpos = z
      @index = index
    end

  end
end