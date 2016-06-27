require "PDB/version"
require "PDB/atom"
require "PDB/molecule"
require "gsl"

module PDB
  # Model contains many molecules derived from PDB File
  # define global variables
  $twoPI = Math::PI*2.0
  $letters = %{A B C D E F G H I J K L M N O P Q R S T U V W X Y Z}

  class Model

    # active_set array of atoms that are in use
    #
    # :molecule => holds a hash, key is chain and value is Molecule class
    #
    attr_reader :filename, :molecule, :active_set, :molecules

    def initialize(file, options={})
      @filename = file
      @active_set =[]
      @molecules = Hash.new # use chain as molecule key for the hash
      raise ArgumentError.new("File does not exists in directory #{file}") unless File.exists?(file)

      openPDBFile(file, options)
    end



    # Open PDB file and create Molecule
    def openPDBFile(filename, flags={})

      puts "OPENING #{filename}"
      pdbLines=[]

      open(filename){|x| pdbLines = x.readlines}

      pdbLines.select!{ |x| x =~ /^ATOM/}
      pdbLines.collect!{|x| PDB::Atom.new(x) }

      # group atoms by chains organized by collection of residues
      current_chain = pdbLines[0].chain
      @molecules[current_chain.to_sym] = PDB::Molecule.new

      pdbLines.each do |atom|

        if (atom.chain != current_chain)
          current_chain = atom.chain
          @molecules[current_chain.to_sym] = PDB::Molecule.new
        end

        # strip out waters
        if atom.residue == "HOH"
          @molecules[:HOH] ||= PDB::Molecule.new
          atom.chain = "HOH" # set chain to WAT
          @molecules[:HOH].addAtom(atom)
        else
          @molecules[current_chain.to_sym].addAtom(atom)
        end
      end


      # create working set based on selection (flags)
      if flags.size == 0
        fillActiveSet
      else
        @active_set.clear
        # :CA => false
        # :CA => true
        selectAtomsByAttribute(flags)
      end

    end



    # Select Atoms by Attribute making a copy of the atom, only selects as AND statements (or intersections)
    #
    # :atom_type => {:CA => true}
    # usage : selectAtomsByAttributes( :atom_type => {:CA => true, :CB => true})
    # usage : selectAtomsByAttributes({:atom_type => {:CA => true, :CB => true}, :chain => {:A => true}})
    # usage : selectAtomsByAttributes({:atom_type => {:CA => true, :CB => true}, :chain => {:B => false}})
    # selectAtomsByAttributes(
    #   {:atom_type => {:CA => true, :CB => true}, :chain => {:B => false}},
    #   {:atom_type => {:O => true}, :residue => {:HOH => true}},
    # )
    # Different Hashes at top level are AND (intersection)
    #
    # selectAtomsByAttribute(:atom_type => {:CB => true}, :residue => {:ASP => true})
    #
    def selectAtomsByAttribute(selection=[])

      @active_set.clear
      #
      # atom_type => {:CA => false, :CB => true}, :chain => {:A => false}
      # apply false to all selected items ???
      # each element of the array is an AND so, we concatenate the output
      if selection.class == Array
        selection.each do |hash_selection|
          @active_set.concat(select_from_molecules_true(hash_selection))
        end
      else
        @active_set.concat(select_from_molecules_true(selection))
      end

      remove_duplicate_atom_selections(@active_set)

#
      # if active_set is empty, assume we select all
      if @active_set.size == 0
        puts "Filling? #{@active_set.size}"
        @active_set.clear
        fillActiveSet
      end

      #
      # remove any unwanted attributes => false flags
      #
      if selection.class == Array
        # hash_selection is something like:
        # atom_type => {:CA => false, :CB => true}
        # selection = [{:atom_type => {:CA => true, :CB => true}}, {:residue => {:HOH => true}}]
        selection.each do |hash_selection|
          # attribute = :atom_type
          select_from_active_set_false(@active_set, hash_selection)
        end
      else
        select_from_active_set_false(@active_set, selection)
      end

    end



    def select_from_active_set_false(array_of_atoms, selection_hash)
      #
      # attribute == :atom_type, :chain, etc
      # keys must be attributes in Atom class
      #
      temp=[]
      selection_hash.each_pair do |attribute, parameters|
        # iterate over each molecule, concatenate for each time it is true
        raise ArgumentError.new("Property not found in Atom Class #{attribute}") unless PDB::Atom.method_defined? attribute
          #
          # type_symbol = :CA, :CB, :CHAIN_IDENTIFIER
          # value = true or false
          parameters.each_pair do |type_symbol, value|

            type = (type_symbol.is_a? Symbol) ? type_symbol.to_s.upcase : type_symbol.upcase

            if !value
              puts "FALSE :#{attribute} => #{value} #{array_of_atoms.size} #{type} SEND => #{array_of_atoms[0].send(attribute)}"
              #temp.concat( array_of_atoms.select{|atom| atom.send(attribute) != type} )
              temp.concat( array_of_atoms.delete_if{|atom| atom.send(attribute) == type} )
            end
          end
      end
    end



    # Example : selectAtomsBy(:atom_type => {:CA => true, :CB => false}, :chain => {:B => false})
    def select_from_molecules_true(selection_hash)
      keys = selection_hash.keys
      attribute = keys[0]

      raise ArgumentError.new("Property not found in Atom Class #{attribute}") unless PDB::Atom.method_defined? attribute

      temp=Array.new
      # create initial selection by selecting over all molecules in Model
      selection_hash[attribute].each_pair do |type_symbol, value|
        @molecules.each_value do |mol|
          type = type_symbol.to_sym
          if value
            temp.concat( mol.selectAtomsBy(attribute, type) )
          end
        end
      end

      # all subsequent selections are applied to temp
      selected_set=[]
      if keys.size > 1
        for i in 1...keys.size
          attribute = keys[i]
          # type_symbol = :CA, :CB, etc
          selection_hash[attribute].each_pair do |type_symbol, value|
            type = type_symbol.to_sym
            puts "SELECTING #{attribute} #{type} <=> #{value} #{selected_set.size}"
            if value
              selected_set.concat(temp.select{|atom| atom.send(attribute) == (type.to_s).upcase})
            end
          end
        end
      else
        selected_set = temp
      end

      # remove duplicates
      remove_duplicate_atom_selections(selected_set)
      return selected_set
    end


    #
    # use this method for removing duplicate entries that may occur during a selection
    # :chain => {:A => true}
    # :atom_type => {:CA = true}
    #
    def remove_duplicate_atom_selections(array_of_atoms)
      array_of_atoms.uniq!{|atom| atom.atom_number || atom.chain }
    end


    def fillActiveSet
      @molecules.each_value do |mol|
        mol.residues.each do |res|
          res.atoms.each do |atom|
            @active_set << atom.dup
          end
        end
      end
    end


    #
    # determine extreme points from active_set
    #
    def setExtrema



    end



    def writeModelToFile(filename)
      writeToFile(@molecule.atoms, filename)
    end



    def writeExtremaToFile(filename)
      writeToFile(@molecule.extrema, filename)
    end



    def writeRandomExtremaToFile(filename)
      writeToFile(@molecule.random_extrema, filename)
    end



    # Write PDB file from array of PDBAtom objects
    def writeToFile(atoms, filename)

      if filename.split(/\./).last != "pdb"
        filename = filename + ".pdb"
      end

      newLines =[]
      count = 1
      current_chain = atoms[0].chain

      atoms.each do |atom|

        # if next chain is different, add ter
        if (current_chain != atom.chain)
          newLines << "TER\n"
          current_chain = atom.chain
        end

        if atom.residue.size == 1
          newLines << sprintf("ATOM  %5s %-4s %3s %1s%4s    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n", count, atom.atom_type, convert_to_3_letter(atom.residue, atom.atom_type), atom.chain, atom.resid, atom.xpos, atom.ypos, atom.zpos, atom.occ, atom.temp, atom.atom)
        elsif atom.residue.size == 3
          newLines << sprintf("ATOM  %5s %-4s %3s %1s%4s    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n", count, atom.atom_type, atom.residue, atom.chain, atom.resid, atom.xpos, atom.ypos, atom.zpos, atom.occ, atom.temp, atom.atom)
        else
          abort("Improper residue in PDB file -> can't write to file")
        end

        count += 1
      end

      newLines << "END\n"
      open(filename, 'w'){|x| newLines.each {|xx| x << xx}}
    end
  end


  # Convert 3-letter residue to 1-letter
  # input is a string
  def convert_to_one_letter(residue)
    res = {
        :GUA => "G",
        :ADE => "A",
        :CYT => "C",
        :URI => "U",
        :THY => "T",
        :Gr  => "G",
        :Ar  => "A",
        :Cr  => "C",
        :Ur  => "U",
        :ALA => "A",
        :ARG => "R",
        :ASN => "N",
        :ASP => "D",
        :ASX => "B",
        :CYS => "C",
        :GLU => "E",
        :GLN => "Q",
        :GLX => "Z",
        :GLY => "G",
        :HIS => "H",
        :ILE => "I",
        :LEU => "L",
        :LYS => "K",
        :MET => "M",
        :PHE => "F",
        :PRO => "P",
        :SER => "S",
        :THR => "T",
        :TRP => "W",
        :TYR => "Y",
        :VAL => "V",
        :SEC => "U",
        :PCA => "J"
    }

    return res[residue]
  end
  module_function :convert_to_one_letter



  def report_log(message)
    (Thread.current[:messages] ||= []) << "#{message}"
    puts message

    #File.open('PDB_errors.txt', 'a') do |file|
    #  (Thread.current[:messages] ||= []).each do |error|
    #    file.puts error
    #  end
    #end
  end
  module_function :report_log


  def report_error(error_message)
    (Thread.current[:errors] ||= []) << "#{error_message}"
    puts error_message

    File.open('PDB_errors.txt', 'a') do |file|
      (Thread.current[:errors] ||= []).each do |error|
        file.puts error
      end
    end
  end
  module_function :report_error


end



