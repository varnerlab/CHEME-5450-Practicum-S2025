import Base.+

function +(buffer::Array{String,1}, line::String)
    push!(buffer, line)
end

"""
    read_reaction_file(path_to_file::String) -> Array{String,1}

Read in a VFF formatted reaction file and return the reaction lines as an array of strings.

### Arguments
- `path_to_file::String`: The path to the VFF formatted reaction file.

### Returns
- `Array{String,1}`: An array of strings containing the reaction lines from the file.
"""
function read_reaction_file(path_to_file::String)::Array{String,1}

    # initialize -
    vff_file_buffer = String[]
    vff_reaction_array = Array{String,1}()

    # Read in the file -
    open("$(path_to_file)", "r") do file
        for line in eachline(file)
            +(vff_file_buffer, line)
        end
    end

    # process -
    for reaction_line ∈ vff_file_buffer
        
        # skip comments and empty lines -
        if (occursin("//", reaction_line) == false && 
            isempty(reaction_line) == false)
        
            # grab -
            push!(vff_reaction_array,reaction_line)
        end
    end

    # return -
    return vff_reaction_array
end