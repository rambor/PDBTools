module PDB

  class Residue

    attr_accessor :atoms, :resid, :resname, :chain

    def initialize(atom)
      @atoms=[]
      @atoms << atom
      @resid = atom.resid
      @chain = atom.chain
      @resname = atom.residue
    end


    def add_atom(atom)
      raise ArgumentError, 'CHAIN does not match' unless atom.chain == @chain
      raise ArgumentError, 'RESID does not match' unless atom.resid == @resid
      @atoms << atom
    end

  end

end
