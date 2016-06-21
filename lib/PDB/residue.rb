module PDB

  class Residue

    attr_accessor :atoms, :resid, :chain, :resname

    def initialize(resname, resid, chain)
      @atoms=[]
      @resid = resid
      @chain = chain
      @resname = resname
    end

    # store reference to the atom
    def add_atom(atom)
      @atoms << atom
    end


    def convert_to_single


    end

  end

end
