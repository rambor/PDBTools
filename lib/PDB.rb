require "PDB/version"
require "PDB/atom"
require "PDB/molecule"
require "gsl"

module PDB
  # Model contains many molecules derived from PDB File
  $twoPI = Math::PI*2.0

  class Model

    attr_reader :filename, :molecule

    def initialize(file, options={})
      @filename = file
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

      @molecule = PDB::Molecule.new()

      # fill Molecule
      pdbLines.each{|atom| @molecule.addAtom(atom) }

      # :CA => false
      # :CA => true
      flags.each_pair do |key, value|
        if PDB::Atom.method_defined? key.downcase
          @molecule.selectAtomsByType(key.downcase, value)
        end
      end

    end


    def writeModelToFile(filename)
      writeToFile(@molecule.atoms, filename)
    end


    def writeExtremaToFile(filename)
      writeToFile(@molecule.extrema, filename)
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



