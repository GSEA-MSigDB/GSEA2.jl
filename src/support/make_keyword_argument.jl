function make_keyword_argument(ke_ar)

    Dict(
        Symbol(kes) => ke_ar[ke] for
        (kes, ke) in (("ex", "exponent"), ("al", "algorithm"), ("n_jo", "number_of_jobs"))
    )

end
