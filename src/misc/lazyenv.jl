_current_env = nothing

function current_env()
    global _current_env
    # if isnothing(_current_env)
    #     _current_env = Gurobi.Env()
    # end
    return _current_env
end