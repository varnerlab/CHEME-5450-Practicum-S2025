function _flux(problem::MyPrimalFluxBalanceAnalysisCalculationModel)

    # initialize -
    results = Dict{String,Any}()
    c = problem.objective; # objective function coefficients

    # bounds -
    fluxbounds = problem.fluxbounds;
    lb = fluxbounds[:, 1]; # lower bounds
    ub = fluxbounds[:, 2]; # upper bounds
        
    # species constraints -
    A = problem.S; # constraint matrix

    # how many variables do we have?
    d = length(c);

    # Setup the problem -
    model = Model(GLPK.Optimizer)
    @variable(model, lb[i,1] <= x[i=1:d] <= ub[i,1], start=0.0) # we have d variables
    
    # set objective function -   
    @objective(model, Max, transpose(c)*x);
    @constraints(model, 
        begin
            A*x == 0 # my material balance constraints 
        end
    );

    # run the optimization -
    optimize!(model)

    # check: was the optimization successful?
    @assert is_solved_and_feasible(model)

    # populate -
    x_opt = value.(x);
    results["argmax"] = x_opt
    results["objective_value"] = objective_value(model);

    # return -
    return results
end


"""
    solve(model::AbstractFluxCalculationModel)
"""
function solve(model::AbstractFluxCalculationModel)
    return _flux(model);
end