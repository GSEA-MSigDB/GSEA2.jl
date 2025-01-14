using Comonicon: @cast, @main

using ProgressMeter: @showprogress

using Random: seed!, shuffle!

using StatsBase: mean, sample, std

function plot(ou, al, fe_, sc, me___, ns, se_, sa_, en; ex = 1.0, up = 4)

    Omics.Error.error_missing(ou)

    Omics.Plot.plot_heat_map(
        joinpath(ou, "set_x_$(ns)_x_enrichment.html"),
        en;
        y = se_,
        x = sa_,
        layout = Dict(
            "title" => Dict("text" => "Enrichment using $(string(al)[6:(end - 2)])"),
            "yaxis" =>
                Dict("title" => Dict("text" => Omics.Strin.coun(lastindex(se_), "Set"))),
            "xaxis" =>
                Dict("title" => Dict("text" => Omics.Strin.coun(lastindex(sa_), ns))),
        ),
    )

    ge_ = .!isnan.(en)

    gs_ = BitVector(undef, lastindex(fe_))

    for id_ in
        view(view(CartesianIndices(en), ge_), sortperm(view(en, ge_)))[Omics.Rank.get_extreme(
        sum(ge_),
        up,
    )]

        ie, ia = Tuple(id_)

        sc_ = view(sc, :, ia)

        gs_ .= .!isnan.(sc_)

        so_ = view(sc_, gs_)

        id_ = sortperm(so_; rev = true)

        ti = "$(sa_[ia]) Enriching $(se_[ie])"

        plot(joinpath(ou, "$ti.html"), al, view(fe_, gs_)[id_], so_[id_], me___[ie]; ex, ti)

    end

end

function _normalize!(ma, di, st)

    if isone(di)

        ea = eachrow

    elseif di == 2

        ea = eachcol

    else

        error("dimension is not 1 or 2.")

    end

    @info "Normalizing dimension $di using standard-deviation $st"

    foreach(Omics.Normalization.normalize_with_0!, ea(ma))

    clamp!(ma, -st, st)

end

function _read_set(js, fe_, mi, ma, fr)

    se_me_ = Omics.Dic.read(js)

    se_ = collect(keys(se_me_))

    me___ = collect(values(se_me_))

    ue = lastindex(se_)

    @info "Selecting ($mi <= N <= $ma && $fr <= %) from $(Omics.Strin.coun(ue, "set"))"

    ke_ = BitVector(undef, ue)

    for (id, me_) in enumerate(me___)

        ub = lastindex(me_)

        intersect!(me_, fe_)

        ua = lastindex(me_)

        ke_[id] = mi <= ua <= ma && fr <= ua / ub

    end

    uk = sum(ke_)

    me = "Selected $(Omics.String.count(uk, "set"))."

    if iszero(uk)

        error(me)

    end

    @info me

    se_[ke_], me___[ke_]

end

function _set_algorithm(al)

    if al == "ks"

        KS()

    elseif al == "ksa"

        KSa()

    elseif al == "kliom"

        KLioM()

    elseif al == "kliop"

        KLioP()

    elseif al == "kli"

        KLi()

    elseif al == "kli1"

        KLi1()

    else

        error("`$al` is not \"ks\", \"ksa\", \"kliom\", \"kliop\", \"kli\", or \"kli1\".")

    end

end

function _error(fe_, nu_)

    Omics.Error.error_duplicate(fe_)

    Omics.Error.error_bad(!isfinite, nu_)

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

    _, ta_, _sa_, nu = Omics.DataFrame.separate(Omics.CLS.read(cls))

    _, fe_, sa_, sc = Omics.DataFrame.separate(Omics.GCT.read(gct))

    ut = lastindex(_sa_)

    uf = lastindex(sa_)

    if ut != uf

        error("numbers of samples differ. $ut (`.cls`) != $uf (`.gct`).")

    end

    _error(fe_, sc)

    Omics.DataFrame.write(target_x_sample_x_number_tsv, "Target", ta_, sa_, nu .- 1.0)

    Omics.DataFrame.write(feature_x_sample_x_score_tsv, "Feature", fe_, sa_, sc)

    target_x_sample_x_number_tsv, feature_x_sample_x_score_tsv

end

"""
Convert one or more `.gmt`s to a `.json`.

# Arguments

  - `set_features_json`: Output `.json`.
  - `gmt_`: Input `.gmt`s.
"""
@cast function convert_gmt(set_features_json, gmt_...)

    Omics.Dic.write(set_features_json, merge!((Omics.GMT.read(gmt) for gmt in gmt_)...))

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

    _, fe_, sa_, sc = Omics.DataFrame.separate(feature_x_sample_x_score_tsv)

    _error(fe_, sc)

    if skip_0

        replace!(sc, 0.0 => NaN)

    end

    if !iszero(normalization_dimension)

        _normalize!(sc, normalization_dimension, normalization_standard_deviation)

    end

    al = _set_algorithm(algorithm)

    se_, me___ = _read_set(
        set_features_json,
        fe_,
        minimum_set_size,
        maximum_set_size,
        minimum_set_fraction,
    )

    en = enrich(al, fe_, sc, me___; um = post_skip_minimum_set_size, ex = exponent)

    Omics.DataFrame.write(
        joinpath(output_directory, "set_x_sample_x_enrichment.tsv"),
        "Set",
        se_,
        sa_,
        en,
    )

    plot(output_directory, al, fe_, sc, me___, "Sample", se_, sa_, en; ex = exponent)

end

function _normalize_enrichment(nu, mn, mp)

    Omics.Number.is_negative(nu) ? -nu / mn : nu / mp

end

function _normalize_enrichment(nu, mn, mp, sn, sp)

    Omics.Number.is_negative(nu) ? -1.0 + (nu - mn) / 3.0 * sn : 1.0 + (nu - mp) / 3.0 * sp

end

function _normalize_enrichment!(::Union{KS, KSa}, en_, ra)

    er_ = Vector{Float64}(undef, lastindex(en_))

    for (id, (en, ra_)) in enumerate(zip(en_, eachrow(ra)))

        ne_, po_ = Omics.Number.separate(ra_)

        mn = mean(ne_)

        mp = mean(po_)

        er_[id] = _normalize_enrichment(en, mn, mp)

        ra_ .= _normalize_enrichment.(ra_, mn, mp)

    end

    er_

end

function _normalize_enrichment!(::Union{KLi1, KLi, KLioM, KLioP}, en_, ra)

    er_ = Vector{Float64}(undef, lastindex(en_))

    for (id, (en, ra_)) in enumerate(zip(en_, eachrow(ra)))

        ne_, po_ = Omics.Number.separate(ra_)

        mn = mean(ne_)

        mp = mean(po_)

        sn = std(ne_)

        sp = std(po_)

        er_[id] = _normalize_enrichment(en, mn, mp, sn, sp)

        ra_ .= _normalize_enrichment.(ra_, mn, mp, sn, sp)

    end

    er_

end

function _write(ou, wr, se_, en_, ra, up, pl_, al, fe_, sc_, me___, ex, nf, ns, nl, nh)

    ue = lastindex(se_)

    re = fill(NaN, ue, 4)

    id_ = sortperm(en_)

    en_ = en_[id_]

    se_ = se_[id_]

    me___ = view(me___, id_)

    ra = ra[id_, :]

    if wr

        Omics.DataFrame.write(
            joinpath(ou, "set_x_index_x_random.tsv"),
            "Set",
            se_,
            1:size(ra, 2),
            ra,
        )

    end

    re[:, 1] = en_

    er_ = _normalize_enrichment!(al, en_, ra)

    re[:, 2] = er_

    id = findlast(Omics.Number.is_negative, en_)

    il_ = 1:id

    ig_ = (id + 1):ue

    re[il_, 3], re[il_, 4], re[ig_, 3], re[ig_, 4] =
        Omics.Statistics.get(ra_, er_, il_, ig_)

    Omics.DataFrame.write(
        joinpath(ou, "set_x_statistic_x_number.tsv"),
        "Set",
        se_,
        ["Enrichment", "Normalized Enrichment", "P-Value", "Q-Value"],
        re,
    )

    for id in unique!(vcat(Omics.Rank.get_extreme(er_, up), indexin(pl_, se_)))

        if isnothing(id)

            continue

        end

        ti = "$id $(se_[id])"

        plot(joinpath(ou, "$ti.html"), al, fe_, sc_, me___[id]; ex, ti, nf, ns, nl, nh)

    end

end

function _permute_set(ur, se, al, fe_, sc_, me___, ex)

    ra = Matrix{Float64}(undef, lastindex(me___), ur)

    if 0 < ur

        @info "Calculating significance by permuting sets"

        le_ = lastindex.(me___)

        seed!(se)

        @showprogress for id in 1:ur

            ra[:, id] =
                enrich(al, fe_, sc_, (le -> sample(fe_, le; replace = false)).(le_); ex)

        end

    end

    ra

end

function _use_permutation(permutation, al, fe_, me___, ex)

    _, ro_, id_, rn = Omics.DataFrame.separate(permutation)

    fe_x_id_x_ra = view(rn, indexin(fe_, ro_), :)

    ur = lastindex(id_)

    ra = Matrix{Float64}(undef, lastindex(me___), ur)

    @info "Calculating significance using predefined $ur random scores"

    @showprogress for (id, ra_) in enumerate(eachcol(fe_x_id_x_ra))

        id_ = sortperm(ra_; rev = true)

        ra[:, id] = enrich(al, fe_[id_], ra_[id_], me___; ex)

    end

    ra

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

    id_ = sortperm(sc_; rev = true)

    sc_ = sc_[id_]

    fe_ = fe_[id_]

    se_, me___ = _read_set(
        set_features_json,
        fe_,
        minimum_set_size,
        maximum_set_size,
        minimum_set_fraction,
    )

    if permutation == "set"

        ra =
            _permute_set(number_of_permutations, random_seed, al, fe_, sc_, me___, exponent)

    elseif isfile(permutation)

        ra = _use_permutation(permutation, al, fe_, me___, exponent)

    else

        error("`$permutation` is not \"set\" or feature_x_index_x_random.tsv.")

    end

    _write(
        output_directory,
        write_set_x_index_x_random_tsv,
        se_,
        enrich(al, fe_, sc_, me___; ex = exponent),
        ra,
        number_of_sets_to_plot,
        split(more_sets_to_plot),
        al,
        fe_,
        sc_,
        me___,
        exponent,
        feature_name,
        score_name,
        low_text,
        high_text,
    )

end

function _get_mean_difference(n1_, n2_)

    mean(n1_) - mean(n2_)

end

function _get_standard_deviation(nu_, me)

    max(0.2 * abs(me), std(nu_; corrected = true))

end

function _get_signal_to_noise_ratio(n1_, n2_)

    m1 = mean(n1_)

    m2 = mean(n2_)

    (m1 - m2) / (_get_standard_deviation(n1_, m1) + _get_standard_deviation(n2_, m2))

end

function _target_sort(fu, is_, sc, fe_)

    sc_ = fu.(eachrow(sc[:, .!is_]), eachrow(sc[:, is_]))

    id_ = sortperm(sc_; rev = true)

    sc_[id_], fe_[id_]

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

    _, ta_, st_, nu = Omics.DataFrame.separate(target_x_sample_x_number_tsv)

    _error(ta_, nu)

    un_ = Set(nu)

    if un_ != Set((0, 1))

        error("target numbers are not all 0 or 1. $un_.")

    end

    _, fe_, sf_, sc = Omics.DataFrame.separate(feature_x_sample_x_score_tsv)

    sc = sc[:, indexin(st_, sf_)]

    _error(fe_, sc)

    if !iszero(normalization_dimension)

        _normalize!(sc, normalization_dimension, normalization_standard_deviation)

    end

    if metric == "mean-difference"

        fu = _get_mean_difference

    elseif metric == "signal-to-noise-ratio"

        fu = _get_signal_to_noise_ratio

    else

        error("`$metric` is not \"mean-difference\" or \"signal-to-noise-ratio\".")

    end

    is_ = convert(BitVector, view(nu, 1, :))

    sc_, fe_ = _target_sort(fu, is_, sc, fe_)

    Omics.DataFrame.write(
        joinpath(output_directory, "feature_x_metric_x_score.tsv"),
        "Feature",
        fe_,
        [metric],
        reshape(sc_, :, 1),
    )

    se_, me___ = _read_set(
        set_features_json,
        fe_,
        minimum_set_size,
        maximum_set_size,
        minimum_set_fraction,
    )

    al = _set_algorithm(algorithm)

    if permutation == "sample"

        ra = Matrix{Float64}(undef, lastindex(me___), number_of_permutations)

        if 0 < number_of_permutations

            @info "Calculating significance by permuting samples"

            seed!(random_seed)

            @showprogress for id in 1:number_of_permutations

                ra_, fer_ = _target_sort(fu, shuffle!(is_), sc, fe_)

                ra[:, id] = enrich(al, fer_, ra_, me___; ex = exponent)

            end

        end

    elseif permutation == "set"

        ra =
            _permute_set(number_of_permutations, random_seed, al, fe_, sc_, me___, exponent)

    elseif isfile(permutation)

        ra = _use_permutation(permutation, al, fe_, me___, exponent)

    else

        error("`$permutation` is not \"sample\", \"set\", or feature_x_index_x_random.tsv.")

    end

    _write(
        output_directory,
        write_set_x_index_x_random_tsv,
        se_,
        enrich(al, fe_, sc_, me___; ex = exponent),
        ra,
        number_of_sets_to_plot,
        split(more_sets_to_plot),
        al,
        fe_,
        sc_,
        me___,
        exponent,
        feature_name,
        score_name,
        low_text,
        high_text,
    )

end

"""
# Gene set enrichment analysis 🏔️
"""
@main
