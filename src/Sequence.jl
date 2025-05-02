
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
function transcription(sequence::String; mrnaprefix::String = "mRNA", geneprefix::String = "gene", 
    genename::String = "test", rnapsymbol::String = "RNAP")::Vector{String}

    # initialize -
    reactions = Vector{String}(); # reactions to hold the output 
    total_ntp = 0; # total number of nucleotides -
    nucleotide_dictionary = _parse_gene_seq(sequence);  # generate the sequence dictionary -
  
    # write the RNAP binding step -
    "transcriptional_initiation_$(genename),$(geneprefix)_$(genename)+$(rnapsymbol),$(rnapsymbol)_OPEN_$(geneprefix)_$(genename),false" |> s-> push!(reactions, s);
    
    # go through by dictionary, and get the base count -
    buffer= ""; # buffer to hold the output
    buffer*="transcription_$(genename),$(rnapsymbol)_OPEN_$(geneprefix)_$(genename)" 
    for (key,value) in nucleotide_dictionary
  
      if (key == "a")
  
        # write the M_atp_c line -
        buffer*="+$(value)*M_atp_c"
  
        # How many a's do we have?
        total_ntp += value;
  
      elseif (key == "t")
  
        # write the M_utp_c line -
        buffer*="+$(value)*M_utp_c"
  
        # How many u's do we have?
        total_ntp += value;
  
      elseif (key == "g")
  
        # write the M_gtp_c line -
        buffer*="+$(value)*M_gtp_c"
  
        # How many g's do we have?
        total_ntp += value;
  
      else
  
        # write the M_gtp_c line -
        buffer*="+$(value)*M_ctp_c"
  
        # How many c's do we have?
        total_ntp += value;
      end
    end
  
    # mRNA+GENE+RNAP+1320*M_pi_c,0,inf;
    buffer*="+$(total_ntp)*M_h2o_c,$(mrnaprefix)_$(genename)+$(geneprefix)_$(genename)+$(rnapsymbol)+$(total_ntp)*M_ppi_c,false"
    buffer |> s-> push!(reactions, s);
  
    # mRNA_decay degradation reaction -
    # mRNA_decay,[],mRNA,144*M_cmp_c+151*M_gmp_c+189*M_ump_c+176*M_amp_c,0,inf;
    buffer = ""; # reset the buffer
    buffer*="$(mrnaprefix)_degradation_$(genename),$(mrnaprefix)_$(genename),"
    local_buffer = "";
    for (key,value) in nucleotide_dictionary
  
      if (key == "a")
  
        # write the M_atp_c line -
        local_buffer*="+$(value)*M_amp_c"
  
      elseif (key == "t")
  
        # write the M_utp_c line -
        local_buffer*="+$(value)*M_ump_c"
  
      elseif (key == "g")
  
        # write the M_gtp_c line -
        local_buffer*="+$(value)*M_gmp_c"
  
      else
  
        # write the M_gtp_c line -
        local_buffer*="+$(value)*M_cmp_c"
  
      end
    end
  
    buffer*="$(local_buffer[2:end]),false"
    buffer |> s-> push!(reactions, s);
  
    # return the buffer -
    return reactions;
end

function translation(sequence::String; mrnaprefix::String = "mRNA", proteinprefix::String = "prot", 
    trnasymbol::String = "tRNA", proteinname::String = "test", ribosomesymbol::String = "RIBOSOME",)::Array{String,1}

    # initialize -
    reactions = Vector{String}(); # reactions to hold the output
    path_to_mapping_file = joinpath(_PATH_TO_CONFIGURATION,"aa_map.csv")
    map_array = readdlm(path_to_mapping_file,','); #metabolite 1, one letter 2
    protein_aa_dictionary = Dict();
  
    # Create a mapping dictionary -
    symbol_metabolite_map = Dict();
    for map_index in collect(1:20)
      one_letter_aa_symbol = map_array[map_index,2];
      metabolite_symbol = map_array[map_index,1];
      symbol_metabolite_map[one_letter_aa_symbol] = metabolite_symbol;
      protein_aa_dictionary[metabolite_symbol] = 0.0;
    end
  
    # Parse the protein seq -
    number_aa_residues = length(sequence);
    local_counter = 0;
    for aa_index in collect(1:number_aa_residues)
  
      # What AA do we have?
      aa_value = string(sequence[aa_index]);
      if (aa_value != "\n" && aa_value != " ")
  
        key = symbol_metabolite_map[aa_value];
  
        # Update the dictionary -
        quantity_aa = protein_aa_dictionary[key];
        protein_aa_dictionary[key] = quantity_aa + 1;
        local_counter+=1;
      end
    end
  
    # Ok, we have the protein sequence , build the reaction string buffer -
    "translation_initiation_$(proteinname),$(mrnaprefix)_$(proteinname)+$(ribosomesymbol),$(ribosomesymbol)_START_$(proteinname),false" |> s-> push!(reactions, s);
    
    # translation reaction -
    buffer="";
    buffer*="translation_$(proteinname),$(ribosomesymbol)_START_$(proteinname)+$(2*local_counter)*M_gtp_c+$(2*local_counter)*M_h2o_c";
    for aa_index in collect(1:20)
  
      # Get charged tRNA -
      metabolite_symbol = map_array[aa_index,1];
  
      # number of this AA -
      value = protein_aa_dictionary[metabolite_symbol];
  
      # Add charged tRNA species to buffer -
      buffer*="+$(value)*$(metabolite_symbol)_$(trnasymbol)";
    end
  
    # products -
    buffer*=",$(ribosomesymbol)+$(mrnaprefix)_$(proteinname)+$(proteinprefix)_$(proteinname)+$(2*local_counter)*M_gdp_c+$(2*local_counter)*M_pi_c+$(local_counter)*$(trnasymbol),false"
    buffer |> s-> push!(reactions, s);
  
    # Write the reactions for charing the tRNA -
    for aa_index in collect(1:20)
  
        buffer = ""; # reset the buffer
        metabolite_symbol = map_array[aa_index,1];
        value = protein_aa_dictionary[metabolite_symbol];
  
        # Add charged tRNA species to buffer -
        buffer*="$(trnasymbol)_charging_$(metabolite_symbol)_$(proteinname),$(value)*$(metabolite_symbol)+$(value)*M_atp_c+$(value)*$(trnasymbol)+$(value)*M_h2o_c,";
        buffer*="$(value)*$(metabolite_symbol)_$(trnasymbol)+$(value)*M_amp_c+$(value)*M_ppi_c,false";
        buffer |> s-> push!(reactions, s);
    end

  
    return reactions;
end

function exchangereactions(reactions::Vector{String})::Array{String,1}
    
    # initialize -
    exchange_reactions = Vector{String}(); # reactions that hold the exchange reactions
    global_species_set = Set{String}(); # global species set -
    

    # process the reactions -
    reaction_species_set = nothing;
    for reaction ∈ reactions
        
      # split the reaction -
      component_array = split(reaction,',');
      reactants = component_array[2] |> String; # reactants
      products = component_array[3] |> String; # products

      # process the reaction, and product sets -
      reactant_set = _extract_species_symbols_set(reactants);
      product_set = _extract_species_symbols_set(products);
      reaction_species_set = reactant_set ∪ product_set; # union of the two sets
      global_species_set = global_species_set ∪ reaction_species_set;
    end

    # convert the set to an array -
    global_species_array = global_species_set |> collect |> sort;

    # build the exchange reactions -
    for species ∈ global_species_array
  
      # build the exchange reaction -
      "exchange_$(species),[],$(species),true" |> s-> push!(exchange_reactions, s);
    end
  
    # return -
    return exchange_reactions;
end

# -- PUBLIC METHODS ABOVE HERE ----------------------------------------------------------------------------------- #