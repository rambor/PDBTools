module PDB
  class Atom

    attr_accessor :atomNumber, :atomType, :residue, :residueType, :resid, :chain, :xpos, :ypos, :zpos, :atom, :occ, :temp, :alt, :mass

    def initialize(line)

      if !(line =~ /ATOM/)
        raise ("Not a PDB formatted line")
      end

      @atomNumber = line[6,5].to_i
      @atomType = line[12,4].strip

      @residue = line[17,3].strip

      @residueType = "" # protein, nucleic acid, etc

      @alt = line[16].strip
      @resid = line[22,4].to_i
      @chain = line[21].strip


      @xpos = line[30,8].to_f
      @ypos = line[38,8].to_f
      @zpos = line[46,8].to_f

      @atom = line[76,2].strip

      if self.atom == "O1P"
        self.atom == "OP1"
      end

      if self.atom == "O2P"
        self.atom == "OP2"
      end

      @occ = line[54,6].to_f
      @temp = line[60,6].to_f

      self.setMass()

    end

    def self.setResidueType(type)
      @residueType = type
    end

    # Convert 3-letter residue to 1-letter
    # input is a string
    def convert_to_one_letter(residue)
      res = {
          "GUA" => "G",
          "ADE" => "A",
          "CYT" => "C",
          "URI" => "U",
          "THY" => "T",
          "Gr" => "G",
          "Ar" => "A",
          "Cr" => "C",
          "Ur" => "U",
          "ALA" => "A",
          "ARG" => "R",
          "ASN" => "N",
          "ASP" => "D",
          "CYS" => "C",
          "GLU" => "E",
          "GLN" => "Q",
          "GLY" => "G",
          "HIS" => "H",
          "ILE" => "I",
          "LEU" => "L",
          "LYS" => "K",
          "MET" => "M",
          "PHE" => "F",
          "PRO" => "P",
          "SER" => "S",
          "THR" => "T",
          "TRP" => "W",
          "TYR" => "Y",
          "VAL" => "V",
          "SEC" => "U",
          "PCA" => "J"
      }

      return res[residue]
    end

    def setMass()
      if (self.atom == "")


      end
      @mass = 12
    end


  end
end




