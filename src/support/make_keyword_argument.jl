function make_keyword_argument(ke_ar)

    Dict(
        Symbol(ke1) => ke_ar[ke2] for
        (ke1, ke2) in [["al", "algorithm"], ["we", "weight"], ["n_jo", "number_of_jobs"]]
    )

end
