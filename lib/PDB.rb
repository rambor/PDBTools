require "PDB/version"
require "PDB/atom"
require "PDB/molecule"
require "gsl"

module PDB
  # Model contains many molecules derived from PDB File
  class Model

    attr_reader :filename, :chains, :molecule

    def initialize(file, options={})
      @filename = file
      raise ArgumentError.new("File does not exists in directory #{file}") unless File.exists?(file)
      @chains = Hash.new

      @molecule = PDB::Molecule.new()
      openPDBFile(file, options)
    end

    # Open PDB file and return array of PDBAtom objects
    def openPDBFile(filename, flags={})

      puts "Opening #{filename}"
      pdbLines=[]

      open(filename){|x| pdbLines = x.readlines}

      puts "BEFORE : #{pdbLines.size}"
      pdbLines.select!{ |x| x =~ /^ATOM/}
      puts " AFTER : #{pdbLines.size}"

      pdbLines.collect!{|x| x = PDB::Atom.new(x) }
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


  def randomRotate(atoms)

    newPositions = atoms.clone

    newCoords = GSL::Vector.alloc(3)
    twoPI = Math::PI*2.0

    # rotate along x
    angle = rand()*twoPI
    xrotation = GSL::Matrix.alloc(3,3)
    cos = Math::cos(angle)
    sin = Math::sin(angle)

    xrotation[0,0] = 1
    xrotation[1,1] = cos
    xrotation[1,2] = -1.0*sin
    xrotation[2,1] = sin
    xrotation[2,2] = cos

    newPositions.each do |atom|
      newCoords[0] = atom.xpos
      newCoords[1] = atom.ypos
      newCoords[2] = atom.zpos
      newCoords = xrotation*newCoords
      atom.xpos = newCoords[0]
      atom.ypos = newCoords[1]
      atom.zpos = newCoords[2]
    end

    # rotate along y
    angle = rand()*twoPI
    yrotation = GSL::Matrix.alloc(3,3)
    cos = Math::cos(angle)
    sin = Math::sin(angle)

    yrotation[0,0] = cos
    yrotation[0,2] = sin
    yrotation[1,1] = 1
    yrotation[2,0] = -sin
    yrotation[2,2] = cos

    atoms.each do |atom|
      newCoords[0] = atom.xpos
      newCoords[1] = atom.ypos
      newCoords[2] = atom.zpos
      newCoords = yrotation*newCoords
      atom.xpos = newCoords[0]
      atom.ypos = newCoords[1]
      atom.zpos = newCoords[2]
    end

    # rotate along z
    angle = rand()*twoPI
    zrotation = GSL::Matrix.alloc(3,3)
    cos = Math::cos(angle)
    sin = Math::sin(angle)

    zrotation[0,0] = cos
    zrotation[0,2] = sin
    zrotation[1,1] = 1
    zrotation[2,0] = -sin
    zrotation[2,2] = cos

    atoms.each do |atom|
      newCoords[0] = atom.xpos
      newCoords[1] = atom.ypos
      newCoords[2] = atom.zpos
      newCoords = zrotation*newCoords
      atom.xpos = newCoords[0]
      atom.ypos = newCoords[1]
      atom.zpos = newCoords[2]
    end

    return newPositions

  end
  module_function :randomRotate

end



