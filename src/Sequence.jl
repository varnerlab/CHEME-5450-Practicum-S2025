
# --- PRIVATE METHODS BELOW HERE ---------------------------------------------------------------------------------- #
function _parse_gene_seq(gene_seq::String)::Dict{String,Int}

    # initialize -
    nucleotide_dictionary = Dict{String,Int}(); # We will return a dictionary w/nucleotide keys and the number of nucleotides as values -
    nucleotide_dictionary["a"] = 0;
    nucleotide_dictionary["t"] = 0;
    nucleotide_dictionary["g"] = 0;
    nucleotide_dictionary["c"] = 0;
    number_of_nucleotides = length(gene_seq); # What is the length of the gene sequence -
  
    # main loop -
    for nucleotide_index ∈ 1:number_of_nucleotides
      
        # get the test nucleotide -
        test_nucleotide = gene_seq[nucleotide_index]; # this gives back a character?
        
        # check for a -
        if (test_nucleotide == 'a')
            value = nucleotide_dictionary["a"];
            nucleotide_dictionary["a"] = value + 1;
        end
  
        # check for t -
        if (test_nucleotide == 't')
            value = nucleotide_dictionary["t"];
            nucleotide_dictionary["t"] = value + 1;
        end
  
        # check for g -
        if (test_nucleotide == 'g')
            value = nucleotide_dictionary["g"];
            nucleotide_dictionary["g"] = value + 1;
        end
  
        # check for c -
        if (test_nucleotide == 'c')
            value = nucleotide_dictionary["c"];
            nucleotide_dictionary["c"] = value + 1;
        end
    end
  
    # return -
    return nucleotide_dictionary;
end
# -- PRIVATE METHODS ABOVE HERE ---------------------------------------------------------------------------------- #

# -- PUBLIC METHODS BELOW HERE ----------------------------------------------------------------------------------- #
function transcription(sequence::String; configuration::Dict{String,Any} = nothing)::Array{String,1}

    # initialize -
    MRNA_type = configuration_dictionary["mRNA_type_prefix"]
    GENE_type = configuration_dictionary["gene_type_prefix"]
    nucleotide_dictionary = _parse_gene_seq(gene_seq);  # generate the sequence dictionary -
  

end

function translation()::Array{String,1}
end
# -- PUBLIC METHODS ABOVE HERE ----------------------------------------------------------------------------------- #