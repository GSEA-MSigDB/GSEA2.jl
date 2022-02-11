function select_set(se_fe_, mi, ma)

    println("Before selecting set:")

    summarize(se_fe_, n_pr = 0)

    se_fe_ = Dict(se => fe_ for (se, fe_) in se_fe_ if mi <= length(fe_) <= ma)

    println("After:")

    summarize(se_fe_, n_pr = 0)

    se_fe_

end
