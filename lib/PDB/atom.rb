module PDB
  class Atom

    attr_accessor :atom_number, :atom_type, :residue, :residue_type, :resid, :chain, :xpos, :ypos, :zpos, :atom, :occ, :temp, :alt, :mass

    @@masses = Hash.new

    @@masses["H"] = 1.008
    @@masses["LI"] = 6.941
    @@masses["C"] = 12.011
    @@masses["N"] = 14.007
    @@masses["O"] = 15.999
    @@masses["NA"] = 22.990
    @@masses["SI"] = 28.086
    @@masses["P"] = 30.974
    @@masses["S"] = 32.066
    @@masses["CL"] = 35.453
    @@masses["K"] = 39.098
    @@masses["CA"] = 40.078
    @@masses["SE"] = 78.971


    def initialize(line)

      if !(line =~ /ATOM/)
        raise ("Not a PDB formatted line")
      end

      @atom_number = line[6,5].to_i
      @atom_type = line[12,4].strip

      @residue = line[17,3].strip

      @residue_type = "" # protein, nucleic acid, etc

      @alt = line[16].strip
      @resid = line[22,4].to_i
      @chain = line[21].strip


      @xpos = line[30,8].to_f
      @ypos = line[38,8].to_f
      @zpos = line[46,8].to_f

      @atom = line[76,2].strip

      if @atom == "O1P"
        @atom == "OP1"
      end

      if @atom == "O2P"
        @atom == "OP2"
      end

      @occ = line[54,6].to_f
      @temp = line[60,6].to_f

      setMass()

    end

    def setResidueType(type)
      @residue_type = type
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
      carbon = ["C", "CA", "CB", "CG", "CD", "CE", "CZ", "CG1", "CG2", "CE1", "CD2", "CD1"]
      nitrogen = ["N", "NE1", "NE2", "NH", "NZ", "ND1"]
      oxygen = ["O", "OD1", "OD2", "OE1", "OE2"]


      if (@atom == "") # not specified, try to infer from atom_type

        if !carbon.find_index(@atom_type).nil?
          @mass = @@masses["C"]
          @atom = "C"
        elsif !nitrogen.find_index(@atom_type).nil?
          @mass = @@masses["N"]
          @atom = "N"
        elsif !oxygen.find_index(@atom_type).nil?
          @mass = @@masses["O"]
          @atom = "O"
        end

      elsif @atom.length > 0
        @mass = @@masses[@atom]
      else
        @mass = 12
      end

    end


  end
end




