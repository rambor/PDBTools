# coding: utf-8
lib = File.expand_path('../lib', __FILE__)
$LOAD_PATH.unshift(lib) unless $LOAD_PATH.include?(lib)
require 'PDB/version'

Gem::Specification.new do |spec|
  spec.name          = "PDB"
  spec.version       = PDB::VERSION
  spec.authors       = ["rambor"]
  spec.email         = ["robert_p_rambo@hotmail.com"]
  spec.summary       = %q{Ruby class for managing PDB files}
  spec.description   = %q{Write a longer description. Optional.}
  spec.homepage      = ""
  spec.license       = "MIT"

  #spec.files         = `git ls-files -z`.split("\x0")
  spec.files         = ["lib/PDB.rb", "lib/PDB/version.rb", "lib/PDB/atom.rb", "lib/PDB/molecule.rb"]
  spec.executables   = spec.files.grep(%r{^bin/}) { |f| File.basename(f) }
  spec.test_files    = spec.files.grep(%r{^(test|spec|features)/})
  spec.require_paths = ["lib"]

  spec.add_development_dependency "bundler", "~> 1.6"
  spec.add_development_dependency "rake"
  spec.add_runtime_dependency "gsl", "~> 2.1"
end
