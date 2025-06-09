follower_c(prob::PricingProblem) = base_costs(prob)

function ct(t, prob::PricingProblem)
    c = follower_c(prob)
    t_full = expand_t(t, prob)
    return c + t_full
end
