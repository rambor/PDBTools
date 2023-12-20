require "PDB/secondary_structure"
require 'set'

module PDB

  class Residue

    # atom definition excludes terminal Oxygen
    AMINO_ACID_SIZES=Hash.new
    AMINO_ACID_SIZES["ALA"]=5
    AMINO_ACID_SIZES["VAL"]=7
    AMINO_ACID_SIZES["ILE"]=8
    AMINO_ACID_SIZES["LEU"]=8
    AMINO_ACID_SIZES["MET"]=8
    AMINO_ACID_SIZES["PHE"]=11
    AMINO_ACID_SIZES["TYR"]=12
    AMINO_ACID_SIZES["TRP"]=14
    AMINO_ACID_SIZES["ARG"]=11
    AMINO_ACID_SIZES["HIS"]=10
    AMINO_ACID_SIZES["LYS"]=9
    AMINO_ACID_SIZES["ASP"]=8
    AMINO_ACID_SIZES["GLU"]=9
    AMINO_ACID_SIZES["SER"]=6
    AMINO_ACID_SIZES["THR"]=7
    AMINO_ACID_SIZES["ASN"]=8
    AMINO_ACID_SIZES["GLN"]=9
    AMINO_ACID_SIZES["CYS"]=6
    AMINO_ACID_SIZES["GLY"]=4
    AMINO_ACID_SIZES["PRO"]=7
    AMINO_ACID_SIZES["PCA"]=7
    AMINO_ACID_SIZES["UNK"]=4

    attr_reader :valid, :resid, :resname, :chain, :ss, :atoms, :alt, :alternates, :icode, :residue_type

    def initialize(atom)
      @atoms=[]
      @atoms << atom
      @resid = atom.resid
      @chain = atom.chain
      @resname = atom.residue
      @icode = atom.icode
      @valid = false

      # @alt=Set(atom.alt)
      @alt=atom.alt
      @alternates = false
    end


    def add_atom(atom)
      raise ArgumentError, 'CHAIN does not match' unless atom.chain == @chain
      raise ArgumentError, 'RESID does not match' unless atom.resid == @resid

      if (@alt.size == 0 && atom.alt.size > 0) # ignore alternates, use first identifier in alt to select rotomer
        @alt = atom.alt
      end

      if (atom.alt.size == 0 || @alt == atom.alt)
        @atoms << atom
      end
      # @alt.add(atom.alt)
      # @alternates = @alt.size > 1 ? true : false
    end

    # determines if the number of non-hydrogen atoms matches required number for residue
    def isValid
      # check that the amino acid is complete
      nonHydrogenAtoms = self.atoms.select{|atom| atom.atom != "H"}
      if (nonHydrogenAtoms.size == AMINO_ACID_SIZES[self.resname])
        @isvalid = true
      end
      puts "ISVALID : #{self.isvalid} #{nonHydrogenAtoms.size} "
    end

    # @param residueToCompare
    def same_conformation(residueToCompare)

      raise ArgumentError, 'Comparing two residues with different number of atoms' unless self.atoms.size == residueToCompare.atoms.size
      raise ArgumentError, ''
      # create distance matrix (linear form) of reference (self) residue
      totalAtoms = self.atoms.size

      refDistances=[]
      for i in 0...totalAtoms
        nextAtom = i + 1
        atom1 = self.atoms[i]
        for j in nextAtom...totalAtoms
          atom2 = atoms[j]
          refDistances << (atom1.xpos - atom2.xpos)**2 + (atom1.ypos - atom2.ypos)**2 + (atom1.zpos - atom2.zpos)**2
        end
      end
      # ADD O atom from residue (k-1) of reference to refDistances

      # populate target insure same order
      targetAtoms = residueToCompare.atoms
      sortedTargetAtoms=[]
      for i in 0...totalAtoms
        atomtype = self.atoms[i].atom_type
        if (atomtype == targetAtoms[i].atom_type)
          sortedTargetAtoms << targetAtoms[i]
        else
          indexIt = targetAtoms.index{|x| x.atom_type == atomtype}
          raise ArgumentError, 'AtomType not found for resid #{self.resid} => #{atomtype}' unless !indexIt.nil?
          sortedTargetAtoms << targetAtoms[indexIt]
        end
      end
      # ADD O atom from residue (k-1) of target to sortedTargetAtoms

      tarDistances=[]
      for i in 0...totalAtoms
        nextAtom = i + 1
        atom1 = sortedTargetAtoms[i]
        for j in nextAtom...totalAtoms
          atom2 = sortedTargetAtoms[j]
          tarDistances << (atom1.xpos - atom2.xpos)**2 + (atom1.ypos - atom2.ypos)**2 + (atom1.zpos - atom2.zpos)**2
        end
      end

      for i in 0...tarDistances.size
        puts "#{i} #{refDistances[i]} <=> #{tarDistances[i]}"
      end

    end

    def printResidue

      @atoms.each do |atom|
        puts "#{atom.atom_number} #{atom.resid} #{atom.atom_type} alt #{atom.alt} icode #{atom.icode}"

      end
    end


    # Accepts string that is the secondary structure assignement of the residue
    # Note: not all residues will have a secondary structure assignment
    #
    # @param value String (HELIX, SHEET, TURN)
    # @param subtype (index)
    def setSecondaryStructure(ssString, subtype)
      @ss = SecondaryStructure.new(ssString, subtype)
    end

    def validate

      temp = @atoms.select{|at| at.atom != "H"}
      if (!AMINO_ACID_SIZES[@resname].nil? && temp.size >= AMINO_ACID_SIZES[@resname])
        @valid = true
      end
    end

    def updateResID(newResId)
      @resid = newResId
      @atoms.each do |atom|
        atom.set_resid(newResId)
      end
    end

    def setResidueType(type)
      @residue_type = type
      @atoms.each do |atom|
        atom.setResidueType(type)
      end
    end

  end
end
