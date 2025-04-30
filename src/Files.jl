function load_protein_sequence_from_file(path_to_protein_file::AbstractString)

    # Load the file, read the sequence into a string and return
    file = open(path_to_protein_file)
    protein_seq = read(file, String)
    close(file)
    return protein_seq
end
  
function load_gene_sequence_from_file(path_to_gene_file::AbstractString)
  
    # Load the file, read the sequence into a string and return
    file = open(path_to_gene_file)
    gene_seq = read(file, String)
    close(file)
    return gene_seq
end