# PDB

Library for reading in a PDB (Protein DataBank File) and manipulating coordinates
An ATOM line from the PDB file is used to create an instance of the Atom Class
Some attribute mappings are described below:

      :atom_number => line[6,5]
        :atom_type => line[12,4]
          :residue => line[17,3]
            :resid => line[22,4]
            :chain => line[21]
             :xpos => line[30,8]
             :ypos => line[38,8]
             :zpos => line[46,8]

Additional attribute assignments can be found in the Atom class (atom.rb).
PDBTools includes limited selection scheme based on attributes of the Atom Class.
The library is meant to be used in a executable scripts that perform operations on the PDB coordinates.

## Installation

Add this line to your application's Gemfile:

    gem 'PDB'

And then execute:

    $ bundle

Or install it yourself as:

    $ gem install PDB

## Usage

TODO: Write usage instructions here

## Contributing

1. Fork it ( https://github.com/[my-github-username]/PDB/fork )
2. Create your feature branch (`git checkout -b my-new-feature`)
3. Commit your changes (`git commit -am 'Add some feature'`)
4. Push to the branch (`git push origin my-new-feature`)
5. Create a new Pull Request
