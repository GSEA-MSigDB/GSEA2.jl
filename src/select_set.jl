function select_set(
    se_fe_::Dict{String, Vector{String}},
    mi::Int64,
    ma::Int64,
)::Dict{String, Vector{String}}

    println("Before selecting set:")

    summarize(se_fe_; n_pr = 0)

    se_fe_ = Dict(se => fe_ for (se, fe_) in se_fe_ if mi <= length(fe_) <= ma)

    println("After:")

    summarize(se_fe_; n_pr = 0)

    return se_fe_

end

export select_set
