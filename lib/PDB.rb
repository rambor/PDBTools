require "PDB/version"
require "PDB/atom"
require "PDB/molecule"
require "gsl"

module PDB
  # Model contains many molecules derived from PDB File
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



