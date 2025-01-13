using Comonicon: @cast, @main

using Random: seed!, shuffle!

using StatsBase: mean, sample, std

function plot(
    di,
    al,
    fe_,
    fe_x_sa_x_sc,
    fe1___,
    nac,
    se_,
    sa_,
    se_x_sa_x_en;
    ex = 1,
    n_pl = 4,
)

    Omics.Error.error_missing(di)

    Omics.Plot.plot_heat_map(
        joinpath(di, "set_x_$(Omics.Path.clean(nac))_x_enrichment.html"),
        se_x_sa_x_en;
        y = se_,
        x = sa_,
        layout = Dict(
            "title" => Dict("text" => "Enrichment using $(string(al)[6:(end - 2)])"),
            "yaxis" => Dict("title" => Dict("text" => "Set ($(lastindex(se_)))")),
            "xaxis" => Dict("title" => Dict("text" => "$nac ($(lastindex(sa_)))")),
        ),
    )

    noe = .!isnan.(se_x_sa_x_en)

    no_ = BitVector(undef, lastindex(fe_))

    for id_ in
        view(view(CartesianIndices(se_x_sa_x_en), noe), sortperm(view(se_x_sa_x_en, noe)))[Omics.Rank.get_extreme(
        sum(noe),
        n_pl,
    )]

        id1, id2 = Tuple(id_)

        sc_ = view(fe_x_sa_x_sc, :, id2)

        no_ .= .!isnan.(sc_)

        scn_ = view(sc_, no_)

        so_ = sortperm(scn_; rev = true)

        title_text = "$(sa_[id2]) Enriching $(se_[id1])"

        plot(
            joinpath(di, "$(Omics.Path.clean(title_text)).html"),
            al,
            view(fe_, no_)[so_],
            scn_[so_],
            fe1___[id1];
            ex,
            title_text,
        )

    end

end

function _normalize!(ma, di, st)

    if isone(di)

        fu = eachrow

    elseif di == 2

        fu = eachcol

    else

        error("Dimension is not 1 or 2.")

    end

    @info "Normalizing dimension $di using standard-deviation $st"

    foreach(Omics.Normalization.normalize_with_0!, fu(ma))

    clamp!(ma, -st, st)

end

function _read_set(js, fe_, mi, ma, fr)

    se_fe1_ = Omics.Dict.read(js)

    se_ = collect(keys(se_fe1_))

    fe1___ = collect(values(se_fe1_))

    n = lastindex(se_)

    @info "Selecting ($mi <= N <= $ma && $fr <= %) from $(Omics.String.count(n, "set"))"

    ke_ = BitVector(undef, n)

    for (id, fe1_) in enumerate(fe1___)

        n1 = lastindex(fe1_)

        intersect!(fe1_, fe_)

        n2 = lastindex(fe1_)

        ke_[id] = mi <= n2 <= ma && fr <= n2 / n1

    end

    n_ke = sum(ke_)

    me = "Selected $(Omics.String.count(n_ke, "set"))."

    if iszero(n_ke)

        error(me)

    end

    @info me

    se_[ke_], fe1___[ke_]

end

function _set_algorithm(al)

    if al == "ks"

        KS()

    elseif al == "ksa"

        KSa()

    elseif al == "kli1"

        KLi1()

    elseif al == "kli"

        KLi()

    elseif al == "kliom"

        KLioM()

    elseif al == "kliop"

        KLioP()

    else

        error("`$al` is not \"ks\", \"ksa\", \"kli1\", \"kli\", \"kliom\", or \"kliop\".")

    end

end

function _error(fe_, an_)

    Omics.Error.error_duplicate(fe_)

    Omics.Error.error_bad(!isfinite, an_)

end

"""
Convert `.cls` and `.gct` to `.tsv`s.

# Arguments

  - `target_x_sample_x_number_tsv`: Output target `.tsv`.
  - `feature_x_sample_x_score_tsv`: Output feature `.tsv`.
  - `cls`: Input `.cls`.
  - `gct`: Input `.gct`.
"""
@cast function convert_cls_gct(
    target_x_sample_x_number_tsv,
    feature_x_sample_x_score_tsv,
    cls,
    gct,
)

    _nat, ta_, _sa_, ta_x_sa_x_nu = Omics.DataFrame.separate(Omics.CLS.read(cls))

    _naf, fe_, sa_, fe_x_sa_x_nu = Omics.DataFrame.separate(Omics.GCT.read(gct))

    n_sac = lastindex(_sa_)

    n_sag = lastindex(sa_)

    if n_sac != n_sag

        error("Numbers of samples differ. $n_sac (`.cls`) != $n_sag (`.gct`).")

    end

    _error(fe_, fe_x_sa_x_nu)

    Omics.DataFrame.write(
        target_x_sample_x_number_tsv,
        "Target",
        ta_,
        sa_,
        ta_x_sa_x_nu .- 1,
    )

    Omics.DataFrame.write(feature_x_sample_x_score_tsv, "Feature", fe_, sa_, fe_x_sa_x_nu)

    target_x_sample_x_number_tsv, feature_x_sample_x_score_tsv

end

"""
Convert one or more `.gmt`s to a `.json`.

# Arguments

  - `set_features_json`: Output `.json`.
  - `gmt_`: Input `.gmt`s.
"""
@cast function convert_gmt(set_features_json, gmt_...)

    Omics.Dict.write(set_features_json, merge!((Omics.GMT.read(gmt) for gmt in gmt_)...))

end

"""
Run data-rank (single-sample) GSEA.

# Arguments

  - `output_directory`:
  - `feature_x_sample_x_score_tsv`:
  - `set_features_json`:

# Options

  - `--normalization-dimension`: = 0. 0 (not normalizing) | 1 | 2.
  - `--normalization-standard-deviation`: = 4.
  - `--minimum-set-size`: = 15.
  - `--maximum-set-size`: = 500.
  - `--minimum-set-fraction`: = 0.0.
  - `--algorithm`: = "ks". "ks" | "ksa" | "kli1" | "kli" | "kliom" | "kliop".
  - `--post-skip-minimum-set-size`: = 1.
  - `--exponent`: = 1.0.

# Flags

  - `--skip-0`: = false. Set this to true for single-cell or other sparse data.
"""
@cast function data_rank(
    output_directory,
    feature_x_sample_x_score_tsv,
    set_features_json;
    skip_0::Bool = false,
    normalization_dimension::Int = 0,
    normalization_standard_deviation::Float64 = 4.0,
    minimum_set_size::Int = 15,
    maximum_set_size::Int = 500,
    minimum_set_fraction::Float64 = 0.0,
    algorithm = "ks",
    post_skip_minimum_set_size::Int = 1,
    exponent::Float64 = 1.0,
)

    Omics.Error.error_missing(output_directory)

    _naf, fe_, sa_, fe_x_sa_x_sc = Omics.DataFrame.separate(feature_x_sample_x_score_tsv)

    _error(fe_, fe_x_sa_x_sc)

    if skip_0

        replace!(fe_x_sa_x_sc, 0 => NaN)

    end

    if !iszero(normalization_dimension)

        _normalize!(fe_x_sa_x_sc, normalization_dimension, normalization_standard_deviation)

    end

    al = _set_algorithm(algorithm)

    se_, fe1___ = _read_set(
        set_features_json,
        fe_,
        minimum_set_size,
        maximum_set_size,
        minimum_set_fraction,
    )

    se_x_sa_x_en = enrich(
        al,
        fe_,
        fe_x_sa_x_sc,
        fe1___;
        mi = post_skip_minimum_set_size,
        ex = exponent,
    )

    Omics.DataFrame.write(
        joinpath(output_directory, "set_x_sample_x_enrichment.tsv"),
        "Set",
        se_,
        sa_,
        se_x_sa_x_en,
    )

    plot(
        output_directory,
        al,
        fe_,
        fe_x_sa_x_sc,
        fe1___,
        "Sample",
        se_,
        sa_,
        se_x_sa_x_en;
        ex = exponent,
    )

end

function _normalize_enrichment(nu, nem, pom)

    Omics.Number.is_negative(nu) ? -nu / nem : nu / pom

end

function _normalize_enrichment(nu, nem, pom, nes, pos)

    Omics.Number.is_negative(nu) ? -1 + (nu - nem) / 3nes : 1 + (nu - pom) / 3pos

end

function _normalize_enrichment!(::Union{KS, KSa}, en_, se_x_id_x_ra)

    enn_ = Vector{Float64}(undef, lastindex(en_))

    for (id, (en, ra_)) in enumerate(zip(en_, eachrow(se_x_id_x_ra)))

        ne_, po_ = Omics.Number.separate(ra_)

        nem = mean(ne_)

        pom = mean(po_)

        enn_[id] = _normalize_enrichment(en, nem, pom)

        ra_ .= _normalize_enrichment.(ra_, nem, pom)

    end

    enn_

end

function _normalize_enrichment!(::Union{KLi1, KLi, KLioM, KLioP}, en_, se_x_id_x_ra)

    enn_ = Vector{Float64}(undef, lastindex(en_))

    for (id, (en, ra_)) in enumerate(zip(en_, eachrow(se_x_id_x_ra)))

        ne_, po_ = Omics.Number.separate(ra_)

        nem = mean(ne_)

        pom = mean(po_)

        nes = std(ne_)

        pos = std(po_)

        enn_[id] = _normalize_enrichment(en, nem, pom, nes, pos)

        ra_ .= _normalize_enrichment.(ra_, nem, pom, nes, pos)

    end

    enn_

end

function _write(
    ou,
    wr,
    se_,
    en_,
    se_x_id_x_ra,
    n_pl,
    pl_,
    al,
    fe_,
    sc_,
    fe1___,
    ex,
    naf,
    nas,
    nal,
    nah,
)

    n_se = lastindex(se_)

    se_x_st_x_nu = fill(NaN, n_se, 4)

    so_ = sortperm(en_)

    en_ = en_[so_]

    se_ = se_[so_]

    fe1___ = view(fe1___, so_)

    se_x_id_x_ra = se_x_id_x_ra[so_, :]

    if wr

        Omics.DataFrame.write(
            joinpath(ou, "set_x_index_x_random.tsv"),
            "Set",
            se_,
            1:size(se_x_id_x_ra, 2),
            se_x_id_x_ra,
        )

    end

    se_x_st_x_nu[:, 1] = en_

    enn_ = _normalize_enrichment!(al, en_, se_x_id_x_ra)

    se_x_st_x_nu[:, 2] = enn_

    idl = findlast(Omics.Number.is_negative, en_)

    nei_ = 1:idl

    poi_ = (idl + 1):n_se

    nep_, nea_, pop_, poa_ = Omics.Statistics.get_p_value(enn_, nei_, poi_, se_x_id_x_ra)

    se_x_st_x_nu[nei_, 3] = nep_

    se_x_st_x_nu[poi_, 3] = pop_

    se_x_st_x_nu[nei_, 4] = nea_

    se_x_st_x_nu[poi_, 4] = poa_

    Omics.DataFrame.write(
        joinpath(ou, "set_x_statistic_x_number.tsv"),
        "Set",
        se_,
        ["Enrichment", "Normalized Enrichment", "P-Value", "Adjusted P-Value"],
        se_x_st_x_nu,
    )

    for id in unique!(vcat(Omics.Rank.get_extreme(enn_, n_pl), indexin(pl_, se_)))

        if isnothing(id)

            continue

        end

        title_text = "$id $(se_[id])"

        plot(
            joinpath(ou, "$(Omics.Path.clean(title_text)).html"),
            al,
            fe_,
            sc_,
            fe1___[id];
            ex,
            title_text,
            naf,
            nas,
            nal,
            nah,
        )

    end

end

function _permute_set(n, se, al, fe_, sc_, fe1___, ex)

    se_x_id_x_ra = Matrix{Float64}(undef, lastindex(fe1___), n)

    if 0 < n

        @info "Calculating significance by permuting sets"

        le_ = lastindex.(fe1___)

        seed!(se)

        @showprogress for id in 1:n

            se_x_id_x_ra[:, id] =
                enrich(al, fe_, sc_, (le -> sample(fe_, le; replace = false)).(le_); ex)

        end

    end

    se_x_id_x_ra

end

function _use_permutation(permutation, al, fe_, fe1___, ex)

    _nar, ro_, id_, ro_x_id_x_ra = Omics.DataFrame.separate(permutation)

    fe_x_id_x_ra = view(ro_x_id_x_ra, indexin(fe_, ro_), :)

    n = lastindex(id_)

    se_x_id_x_ra = Matrix{Float64}(undef, lastindex(fe1___), n)

    @info "Calculating significance using predefined $n random scores"

    @showprogress for (id, ra_) in enumerate(eachcol(fe_x_id_x_ra))

        so_ = sortperm(ra_; rev = true)

        se_x_id_x_ra[:, id] = enrich(al, fe_[so_], ra_[so_], fe1___; ex)

    end

    se_x_id_x_ra

end

"""
Run user-rank (pre-rank) GSEA.

# Arguments

  - `output_directory`:
  - `feature_x_metric_x_score_tsv`:
  - `set_features_json`:

# Options

  - `--minimum-set-size`: = 15.
  - `--maximum-set-size`: = 500.
  - `--minimum-set-fraction`: = 0.0.
  - `--algorithm`: = "ks". "ks" | "ksa" | "kli1" | "kli" | "kliom" | "kliop".
  - `--exponent`: = 1.0.
  - `--permutation`: = "set". "set" | feature_x_index_x_random.tsv.
  - `--number-of-permutations`: = 100.
  - `--random-seed`: = 20150603.
  - `--number-of-sets-to-plot`: = 4.
  - `--more-sets-to-plot`: = "". Space-separated set names.
  - `--feature-name`: = "Gene".
  - `--score-name`: = "User-Defined Score".
  - `--low-text`: = "Low Side".
  - `--high-text`: = "High Side".

# Flags

  - `--write-set-x-index-x-random-tsv`: = false.
"""
@cast function user_rank(
    output_directory,
    feature_x_metric_x_score_tsv,
    set_features_json;
    minimum_set_size::Int = 15,
    maximum_set_size::Int = 500,
    minimum_set_fraction::Float64 = 0.0,
    algorithm = "ks",
    exponent::Float64 = 1.0,
    permutation = "set",
    number_of_permutations::Int = 100,
    random_seed::Int = 20150603,
    write_set_x_index_x_random_tsv::Bool = false,
    number_of_sets_to_plot::Int = 4,
    more_sets_to_plot = "",
    feature_name = "Gene",
    score_name = "User-Defined Score",
    low_text = "Low Side",
    high_text = "High Side",
)

    Omics.Error.error_missing(output_directory)

    al = _set_algorithm(algorithm)

    feature_x_metric_x_score = Omics.DataFrame.read(feature_x_metric_x_score_tsv)

    fe_ = feature_x_metric_x_score[!, 1]

    sc_ = feature_x_metric_x_score[!, 2]

    _error(fe_, sc_)

    so_ = sortperm(sc_; rev = true)

    sc_ = sc_[so_]

    fe_ = fe_[so_]

    se_, fe1___ = _read_set(
        set_features_json,
        fe_,
        minimum_set_size,
        maximum_set_size,
        minimum_set_fraction,
    )

    if permutation == "set"

        se_x_id_x_ra = _permute_set(
            number_of_permutations,
            random_seed,
            al,
            fe_,
            sc_,
            fe1___,
            exponent,
        )

    elseif isfile(permutation)

        se_x_id_x_ra = _use_permutation(permutation, al, fe_, fe1___, exponent)

    else

        error("`$permutation` is not \"set\" or feature_x_index_x_random.tsv.")

    end

    _write(
        output_directory,
        write_set_x_index_x_random_tsv,
        se_,
        enrich(al, fe_, sc_, fe1___; ex = exponent),
        se_x_id_x_ra,
        number_of_sets_to_plot,
        split(more_sets_to_plot),
        al,
        fe_,
        sc_,
        fe1___,
        exponent,
        feature_name,
        score_name,
        low_text,
        high_text,
    )

end

function _get_mean_difference(nu1_, nu2_)

    mean(nu1_) - mean(nu2_)

end

function _get_standard_deviation(nu_, me)

    max(0.2abs(me), std(nu_; corrected = true))

end

function _get_signal_to_noise_ratio(nu1_, nu2_)

    me1 = mean(nu1_)

    me2 = mean(nu2_)

    (me1 - me2) / (_get_standard_deviation(nu1_, me1) + _get_standard_deviation(nu2_, me2))

end

function _target_sort(fu, is_, fe_x_sa_x_sc, fe_)

    sc_ = fu.(eachrow(fe_x_sa_x_sc[:, .!is_]), eachrow(fe_x_sa_x_sc[:, is_]))

    so_ = sortperm(sc_; rev = true)

    sc_[so_], fe_[so_]

end

"""
Run metric-rank (standard) GSEA.

# Arguments

  - `output_directory`:
  - `target_x_sample_x_number_tsv`:
  - `feature_x_sample_x_score_tsv`:
  - `set_features_json`:

# Options

  - `--minimum-set-size`: = 15.
  - `--maximum-set-size`: = 500.
  - `--minimum-set-fraction`: = 0.0.
  - `--normalization-dimension`: = 0. 0 (not normalizing) | 1 | 2.
  - `--normalization-standard-deviation`: = 4.
  - `--metric`: = "signal-to-noise-ratio". "mean-difference" | "signal-to-noise-ratio".
  - `--algorithm`: = "ks". "ks" | "ksa" | "kli1" | "kli" | "kliom" | "kliop".
  - `--exponent`: = 1.0.
  - `--permutation`: = "sample". "sample" | "set" | feature_x_index_x_random.tsv.
  - `--number-of-permutations`: = 100.
  - `--random-seed`: = 20150603.
  - `--number-of-sets-to-plot`: = 4.
  - `--more-sets-to-plot`: = "". Space-separated set names.
  - `--feature-name`: = "Gene".
  - `--score-name`: = "Signal-to-Noise Ratio".
  - `--low-text`: = "Low Side".
  - `--high-text`: = "High Side".

# Flags

  - `--write-set-x-index-x-random-tsv`: = false.
"""
@cast function metric_rank(
    output_directory,
    target_x_sample_x_number_tsv,
    feature_x_sample_x_score_tsv,
    set_features_json;
    minimum_set_size::Int = 15,
    maximum_set_size::Int = 500,
    minimum_set_fraction::Float64 = 0.0,
    normalization_dimension::Int = 0,
    normalization_standard_deviation::Float64 = 4.0,
    metric = "signal-to-noise-ratio",
    algorithm = "ks",
    exponent::Float64 = 1.0,
    permutation = "sample",
    number_of_permutations::Int = 100,
    random_seed::Int = 20150603,
    write_set_x_index_x_random_tsv::Bool = false,
    number_of_sets_to_plot::Int = 4,
    more_sets_to_plot = "",
    feature_name = "Gene",
    score_name = "Signal-to-Noise Ratio",
    low_text = "Low Side",
    high_text = "High Side",
)

    Omics.Error.error_missing(output_directory)

    _nat, ta_, sat_, ta_x_sa_x_nu = Omics.DataFrame.separate(target_x_sample_x_number_tsv)

    _error(ta_, ta_x_sa_x_nu)

    un_ = Set(ta_x_sa_x_nu)

    if un_ != Set((0, 1))

        error("Target numbers are not all 0 or 1. $un_.")

    end

    _naf, fe_, saf_, fe_x_sa_x_sc = Omics.DataFrame.separate(feature_x_sample_x_score_tsv)

    fe_x_sa_x_sc = fe_x_sa_x_sc[:, indexin(sat_, saf_)]

    _error(fe_, fe_x_sa_x_sc)

    if !iszero(normalization_dimension)

        _normalize!(fe_x_sa_x_sc, normalization_dimension, normalization_standard_deviation)

    end

    if metric == "mean-difference"

        fu = _get_mean_difference

    elseif metric == "signal-to-noise-ratio"

        fu = _get_signal_to_noise_ratio

    else

        error("`$metric` is not \"mean-difference\" or \"signal-to-noise-ratio\".")

    end

    is_ = convert(BitVector, view(ta_x_sa_x_nu, 1, :))

    sc_, fe_ = _target_sort(fu, is_, fe_x_sa_x_sc, fe_)

    Omics.DataFrame.write(
        joinpath(output_directory, "feature_x_metric_x_score.tsv"),
        "Feature",
        fe_,
        [metric],
        reshape(sc_, :, 1),
    )

    se_, fe1___ = _read_set(
        set_features_json,
        fe_,
        minimum_set_size,
        maximum_set_size,
        minimum_set_fraction,
    )

    al = _set_algorithm(algorithm)

    if permutation == "sample"

        se_x_id_x_ra = Matrix{Float64}(undef, lastindex(fe1___), number_of_permutations)

        if 0 < number_of_permutations

            @info "Calculating significance by permuting samples"

            seed!(random_seed)

            @showprogress for id in 1:number_of_permutations

                ra_, fer_ = _target_sort(fu, shuffle!(is_), fe_x_sa_x_sc, fe_)

                se_x_id_x_ra[:, id] = enrich(al, fer_, ra_, fe1___; ex = exponent)

            end

        end

    elseif permutation == "set"

        se_x_id_x_ra = _permute_set(
            number_of_permutations,
            random_seed,
            al,
            fe_,
            sc_,
            fe1___,
            exponent,
        )

    elseif isfile(permutation)

        se_x_id_x_ra = _use_permutation(permutation, al, fe_, fe1___, exponent)

    else

        error("`$permutation` is not \"sample\", \"set\", or feature_x_index_x_random.tsv.")

    end

    _write(
        output_directory,
        write_set_x_index_x_random_tsv,
        se_,
        enrich(al, fe_, sc_, fe1___; ex = exponent),
        se_x_id_x_ra,
        number_of_sets_to_plot,
        split(more_sets_to_plot),
        al,
        fe_,
        sc_,
        fe1___,
        exponent,
        feature_name,
        score_name,
        low_text,
        high_text,
    )

end

"""
The official command-line program for gene-set-enrichment analysis (GSEA). Learn more at https://github.com/KwatMDPhD/GSEA.jl.
"""
@main
