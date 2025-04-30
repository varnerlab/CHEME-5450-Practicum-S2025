# setup paths -
const _ROOT = @__DIR__
const _PATH_TO_SRC = joinpath(_ROOT, "src");
const _PATH_TO_DATA = joinpath(_ROOT, "data");
const _PATH_TO_CONFIGURATION = joinpath(_ROOT, "configuration");

# flag to check if the include file was called -
const _DID_INCLUDE_FILE_GET_CALLED = true;

# check: do we have a manifest file?
using Pkg
if (isfile(joinpath(_ROOT, "Manifest.toml")) == false) # have manifest file, we are good. Otherwise, we need to instantiate the environment
    Pkg.activate("."); Pkg.resolve(); Pkg.instantiate(); Pkg.update();
end

# load external packages
using DelimitedFiles
using JSON
using DataFrames
using JLD2
using CSV
using FileIO
using LinearAlgebra
using Statistics
using PrettyTables
using JuMP
using GLPK
using Test

# load my codes -
include(joinpath(_PATH_TO_SRC, "Types.jl"));
include(joinpath(_PATH_TO_SRC, "Factory.jl"));
include(joinpath(_PATH_TO_SRC, "Parser.jl"));
include(joinpath(_PATH_TO_SRC, "Flux.jl"));
include(joinpath(_PATH_TO_SRC, "Files.jl"));
include(joinpath(_PATH_TO_SRC, "Sequence.jl"));