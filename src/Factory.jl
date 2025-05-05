function build(modeltype::Type{MyPrimalFluxBalanceAnalysisCalculationModel}, 
    data::NamedTuple)::MyPrimalFluxBalanceAnalysisCalculationModel

    # get data -
    S = data.S;
    fluxbounds = data.fluxbounds;
    objective = data.objective;
    species = data.species;
    reactions = data.reactions;
 
    # build an empty model -
    model = modeltype();
 
    # add data to the model -
    model.S = S;
    model.fluxbounds = fluxbounds;
    model.objective = objective;
    model.species = species;
    model.reactions = reactions;
     
    # return the model -
    return model;
end