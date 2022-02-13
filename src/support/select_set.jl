function select_set(se_fe_, mi, ma; in_ = [])

    println("Before selecting set:")

    summarize(se_fe_, n_pr = 0)

    if !isempty(in_)

        println("Intersecting")

        se_fe_ = Dict(se => intersect(fe_, in_) for (se, fe_) in se_fe_)

    end

    se_fe_ = Dict(se => fe_ for (se, fe_) in se_fe_ if mi <= length(fe_) <= ma)

    println("After:")

    summarize(se_fe_, n_pr = 0)

    se_fe_

end
