module PDB
  class Fasta

    attr_reader :residues, :total_residues

    def initialize(filename)

      puts "OPENING #{filename}"
      PDB::report_log("OPENING #{filename}")
      lines=[]
      open(filename){|x| lines = x.readlines}
      @residues=[]

      begin

        lines.each do |line|
          if ( (line =~ /^>/).nil? )

            tempresidues = line.chomp.strip.split("")

            tempresidues.each do |res|

              if !((res =~ /[GACTRNDEQHILKMFPSWYVUJ]/).nil?)
                @residues << res
              else
                raise ("RESIDUE NOT FOUND IN LIB #{res}")
              end

            end
          end
        end

      rescue => error
        PDB::report_error("#{error.class} and #{error.message} : #{lines.inspect} ")
      end

      @total_residues = @residues.size
      puts @residues.inspect

    end

  end
end