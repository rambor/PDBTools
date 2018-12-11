module PDB
  class Atom

    attr_accessor :atom_number, :chain, :xpos, :ypos, :zpos, :occ, :temp, :alt
    attr_reader :atom_type, :residue, :residue_type, :resid, :icode, :mass, :atom

    @@masses = Hash.new

    @@masses["H"]  = 1.008
    @@masses["LI"] = 6.941
    @@masses["C"]  = 12.011
    @@masses["N"]  = 14.007
    @@masses["O"]  = 15.999
    @@masses["NA"] = 22.990
    @@masses["SI"] = 28.086
    @@masses["P"]  = 30.974
    @@masses["S"]  = 32.066
    @@masses["CL"] = 35.453
    @@masses["K"]  = 39.098
    @@masses["CA"] = 40.078
    @@masses["SE"] = 78.971

    def initialize(line)

      if !(line =~ /ATOM/)
        raise ("Not a PDB formatted line")
      end

      @atom_number = line[6,5].to_i
      @atom_type = line[12,4].strip

      @residue = line[17,3].strip # such as GLY, ADE

      @residue_type = "" # protein, nucleic acid, etc

      @alt = line[16].strip # alternate location

      @resid = line[22,4].to_i
      @icode = line[26].strip # insertion
      @chain = line[21].strip


      @xpos = line[30,8].to_f
      @ypos = line[38,8].to_f
      @zpos = line[46,8].to_f

      @atom = line[76,2].strip
      @old = line

      #if (@atom =~ /[0-9]+/)
        @atom = ""
      #end

      if @atom_type == "O1P"
        @atom_type = "OP1"
        @atom = "P"
      end

      if @atom_type == "O2P"
        @atom_type = "OP2"
        @atom = "P"
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

    def centerOn(x,y,z)
      @xpos -= x
      @ypos -= y
      @zpos -= z
    end

    def setMass()

      carbon = ["C", "CA", "CB", "CG", "CD", "CE", "CZ", "CG1", "CG2", "CE1", "CE2", "CE3", "CZ2", "CZ3", "CD2", "CD1", "CH2", "CH3"]
      nitrogen = ["N", "NE", "NE1", "NE2", "NH", "NH1", "NH2", "NZ", "ND1", "ND2"]
      oxygen = ["O", "OD1", "OD2", "OE1", "OE2", "OG", "OXT"]

      if (@atom == "") # not specified, try to infer from atom_type

        if !carbon.find_index(@atom_type).nil?
          @mass = @@masses["C"]
          @atom = "C"
        elsif !nitrogen.find_index(@atom_type).nil?
          @mass = @@masses["N"]
          @atom = "N"
        elsif !oxygen.find_index(@atom_type).nil? || @atom_type[0] == "O"
          @mass = @@masses["O"]
          @atom = "O"
        elsif @atom_type[0] == "H"
          @mass = @@masses["H"]
          @atom = "H"
        elsif @atom_type[0] == "S"
          @mass = @@masses["S"]
          @atom = "S"
        elsif @atom_type.size > 1 && @atom_type[0] =~ /[1-9]/ && @atom_type[1] == "H"
          @mass = @@masses["H"]
          @atom = "H"
        else
          @atom = "C"
          @mass = 12
        end

      elsif @atom.length > 0
        @mass = @@masses[@atom]
      else
        @atom = "C"
        @mass = 12
        puts "not found"
        exit
      end

      if (@mass.nil? || @mass == 0)
        puts "NOT SET => #{@residue} #{@resid} #{@atom_type} #{@atom_type.size} | #{@mass} | #{@atom} \n#{@old}"
        exit
      end

    end


  end
end




