function blank_model(; silent=false, threads=nothing)
    model = Model(() -> SCIP.Optimizer())
    set_optimizer_attribute(model, MOI.Silent(), silent)
    # set_optimizer_attribute(model, MOI.NumberOfThreads(), threads)
    return model
end