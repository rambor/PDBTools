

module PDB

  class SecondaryStructure

    attr_reader :structure, :subtype_index

    def initialize(type_of, subclass_index)
      @structure = type_of
      @subtype_index = subclass_index
    end

  end
end
