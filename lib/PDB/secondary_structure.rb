

module PDB

  # From Stride Documentation
  # H alpha helix
  # G 3-10
  # I PI helix
  # Extended conformation
  # T turn - many types
  # C coil
  class SecondaryStructure

    attr_reader :structure, :subtype_index

    def initialize(type_of, subclass_index)
      @structure = type_of
      @subtype_index = subclass_index
    end

  end


end
