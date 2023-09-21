module GSEA

using Comonicon: @cast, @main

using ProgressMeter: @showprogress

using Random: seed!, shuffle!

using StatsBase: mean, sample, std

using BioLab

include("FeatureSetEnrichment.jl")

function _normalize_with_0_clamp!(ma, di, st)

    if isone(di)

        fu = eachrow

    elseif di == 2

        fu = eachcol

    else

        error("Dimension $di is not 1 or 2.")

    end

    @info "Normalizing dimension $di using standard deviation $st"

    foreach(BioLab.Normalization.normalize_with_0!, fu(ma))

    clamp!(ma, -st, st)

end

function _read_set(js, fe_, mi, ma, mif)

    se_fe1_ = BioLab.Dict.read(js, Dict{String, Vector{String}})

    se_ = collect(keys(se_fe1_))

    fe1___ = collect(values(se_fe1_))

    n = length(se_)

    @info "Selecting ($mi <= N <= $ma && $mif <= %) from $(BioLab.String.count(n, "set"))"

    ke_ = BitVector(undef, n)

    for (id, fe1_) in enumerate(fe1___)

        n1 = length(fe1_)

        intersect!(fe1_, fe_)

        n2 = length(fe1_)

        ke_[id] = mi <= n2 <= ma && mif <= n2 / n1

    end

    n_ke = sum(ke_)

    if iszero(n_ke)

        error("Selected 0 set.")

    end

    @info "Selected $(BioLab.String.count(n_ke, "set"))."

    se_[ke_], fe1___[ke_]

end

function _set_algorithm(al)

    if al == "ks"

        FeatureSetEnrichment.KS()

    elseif al == "ksa"

        FeatureSetEnrichment.KSa()

    elseif al == "kli1"

        FeatureSetEnrichment.KLi1()

    elseif al == "kli"

        FeatureSetEnrichment.KLi()

    elseif al == "kliom"

        FeatureSetEnrichment.KLioM()

    elseif al == "kliop"

        FeatureSetEnrichment.KLioP()

    else

        error(
            "Algorithm \"$al\" is not \"ks\", \"ksa\", \"kli1\", \"kli\", \"kliom\", or \"kliop\".",
        )

    end

end

"""
Convert `.cls` and `.gct` to `.tsv`s.

# Arguments

  - `target_x_sample_x_number_tsv`: Output target `.tsv`.
  - `feature_x_sample_x_score_tsv`: Output feature `.tsv`.
  - `cls`: Input `.cls`.
  - `gct`: Input `.gct`.
"""
#@cast function convert_cls_gct(
function convert_cls_gct(target_x_sample_x_number_tsv, feature_x_sample_x_score_tsv, cls, gct)

    _nat, ta_, _sa_, ta_x_sa_x_nu = BioLab.DataFrame.separate(BioLab.CLS.read(cls))

    BioLab.Error.error_bad(ta_)

    BioLab.Error.error_bad(ta_x_sa_x_nu)

    _naf, fe_, sa_, fe_x_sa_x_nu = BioLab.DataFrame.separate(BioLab.GCT.read(gct))

    if length(_sa_) != length(sa_)

        error("Sample (column) lengths differ.")

    end

    BioLab.Error.error_bad(fe_)

    BioLab.Error.error_duplicate(fe_)

    BioLab.Error.error_bad(fe_x_sa_x_nu)

    BioLab.DataFrame.write(
        target_x_sample_x_number_tsv,
        BioLab.DataFrame.make("Target", ta_, sa_, ta_x_sa_x_nu .- 1),
    )

    BioLab.DataFrame.write(
        feature_x_sample_x_score_tsv,
        BioLab.DataFrame.make("Feature", fe_, sa_, fe_x_sa_x_nu),
    )

    target_x_sample_x_number_tsv, feature_x_sample_x_score_tsv

end

"""
Convert one or more `.gmt`s to a `.json`.

# Arguments

  - `set_features_json`: Output `.json`.
  - `gmt_`: Input `.gmt`s.
"""
#@cast function convert_gmt(set_features_json, gmt_...)
function convert_gmt(set_features_json, gmt_...)

    BioLab.Dict.write(set_features_json, merge((BioLab.GMT.read(gmt) for gmt in gmt_)...))

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

  - `--skip-0`: = false.
"""
#@cast function data_rank(
function data_rank(
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

    BioLab.Error.error_missing(output_directory)

    _naf, fe_, sa_, fe_x_sa_x_sc =
        BioLab.DataFrame.separate(BioLab.DataFrame.read(feature_x_sample_x_score_tsv))

    BioLab.Error.error_bad(fe_)

    BioLab.Error.error_duplicate(fe_)

    BioLab.Error.error_bad(fe_x_sa_x_sc)

    if skip_0

        replace!(fe_x_sa_x_sc, 0 => NaN)

    end

    if !iszero(normalization_dimension)

        _normalize_with_0_clamp!(
            fe_x_sa_x_sc,
            normalization_dimension,
            normalization_standard_deviation,
        )

    end

    al = _set_algorithm(algorithm)

    se_, fe1___ =
        _read_set(set_features_json, fe_, minimum_set_size, maximum_set_size, minimum_set_fraction)

    se_x_sa_x_en = FeatureSetEnrichment.enrich(
        al,
        fe_,
        fe_x_sa_x_sc,
        fe1___;
        mi = post_skip_minimum_set_size,
        ex = exponent,
    )

    BioLab.DataFrame.write(
        joinpath(output_directory, "set_x_sample_x_enrichment.tsv"),
        BioLab.DataFrame.make("Set", se_, sa_, se_x_sa_x_en),
    )

    FeatureSetEnrichment.plot(
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

    output_directory

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

    n_se = length(se_)

    se_x_st_x_nu = fill(NaN, n_se, 4)

    so_ = sortperm(en_)

    en_ = view(en_, so_)

    se_ = view(se_, so_)

    fe1___ = view(fe1___, so_)

    se_x_id_x_ra = view(se_x_id_x_ra, so_, :)

    se_x_st_x_nu[:, 1] = en_

    if wr

        BioLab.DataFrame.write(
            joinpath(ou, "set_x_index_x_random.tsv"),
            BioLab.DataFrame.make("Set", se_, collect(1:size(se_x_id_x_ra, 2)), se_x_id_x_ra),
        )

    end

    nef_ = Vector{Float64}(undef, n_se)

    pof_ = Vector{Float64}(undef, n_se)

    for id in 1:n_se

        ne_ = Vector{Float64}()

        po_ = Vector{Float64}()

        for ra in view(se_x_id_x_ra, id, :)

            if ra < 0

                push!(ne_, ra)

            else

                push!(po_, ra)

            end

        end

        nef_[id] = -1 / mean(ne_)

        pof_[id] = 1 / mean(po_)

    end

    idl = findlast(<(0), en_)

    enn_ = Vector{Float64}(undef, n_se)

    for id in eachindex(enn_)

        if id <= idl

            fa_ = nef_

        else

            fa_ = pof_

        end

        enn_[id] = en_[id] * fa_[id]

    end

    se_x_st_x_nu[:, 2] = enn_

    nei_ = 1:idl

    poi_ = (idl + 1):n_se

    npv_, nad_, ppv_, pad_ =
        BioLab.Statistics.get_p_value(enn_, nei_, poi_, se_x_id_x_ra; nef_, pof_)

    se_x_st_x_nu[nei_, 3] = npv_

    se_x_st_x_nu[poi_, 3] = ppv_

    se_x_st_x_nu[nei_, 4] = nad_

    se_x_st_x_nu[poi_, 4] = pad_

    BioLab.DataFrame.write(
        joinpath(ou, "set_x_statistic_x_number.tsv"),
        BioLab.DataFrame.make(
            "Set",
            se_,
            ["Enrichment", "Normalized Enrichment", "P-Value", "Adjusted P-Value"],
            se_x_st_x_nu,
        ),
    )

    for id in unique(vcat(BioLab.Rank.get_extreme(en_, n_pl), indexin(pl_, se_)))

        if isnothing(id)

            continue

        end

        title_text = "$id $(se_[id])"

        FeatureSetEnrichment.plot(
            joinpath(ou, "$(BioLab.Path.clean(title_text)).html"),
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

    ou

end

function _permute_set(n_pe, se, al, fe_, sc_, fe1___, ex)

    se_x_id_x_ra = Matrix{Float64}(undef, length(fe1___), n_pe)

    if 0 < n_pe

        @info "Permuting sets to compute significance"

        le_ = length.(fe1___)

        seed!(se)

        @showprogress for id in 1:n_pe

            se_x_id_x_ra[:, id] = FeatureSetEnrichment.enrich(
                al,
                fe_,
                sc_,
                [sample(fe_, le; replace = false) for le in le_];
                ex,
            )

        end

    end

    se_x_id_x_ra

end

function _use_permutation(permutation, al, fe_, fe1___, ex)

    feature_x_index_x_random = BioLab.DataFrame.read(permutation)

    fe_x_id_x_ra = view(
        Matrix(feature_x_index_x_random[!, 2:end]),
        indexin(fe_, feature_x_index_x_random[!, 1]),
        :,
    )

    n_pe = size(fe_x_id_x_ra, 2)

    se_x_id_x_ra = Matrix{Float64}(undef, length(fe1___), n_pe)

    @info "Using predefined $n_pe random scores to compute significance"

    @showprogress for id in 1:n_pe

        ra_ = fe_x_id_x_ra[:, id]

        so_ = sortperm(ra_; rev = true)

        se_x_id_x_ra[:, id] =
            FeatureSetEnrichment.enrich(al, view(fe_, so_), view(ra_, so_), fe1___; ex)

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
#@cast function user_rank(
function user_rank(
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

    BioLab.Error.error_missing(output_directory)

    feature_x_metric_x_score = BioLab.DataFrame.read(feature_x_metric_x_score_tsv)

    fe_ = feature_x_metric_x_score[!, 1]

    BioLab.Error.error_bad(fe_)

    BioLab.Error.error_duplicate(fe_)

    sc_ = feature_x_metric_x_score[!, 2]

    BioLab.Error.error_bad(sc_)

    so_ = sortperm(sc_; rev = true)

    sc_ = view(sc_, so_)

    fe_ = view(fe_, so_)

    se_, fe1___ =
        _read_set(set_features_json, fe_, minimum_set_size, maximum_set_size, minimum_set_fraction)

    al = _set_algorithm(algorithm)

    if permutation == "set"

        se_x_id_x_ra =
            _permute_set(number_of_permutations, random_seed, al, fe_, sc_, fe1___, exponent)

    elseif isfile(permutation)

        se_x_id_x_ra = _use_permutation(permutation, al, fe_, fe1___, exponent)

    else

        error("Permutation \"$permutation\" is not \"set\" or feature_x_index_x_random.tsv.")

    end

    _write(
        output_directory,
        write_set_x_index_x_random_tsv,
        se_,
        FeatureSetEnrichment.enrich(al, fe_, sc_, fe1___; ex = exponent),
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

function _get_standard_deviation(nu_, me)

    fr = 0.2

    #if iszero(me)

    #    return fr

    #end

    max(abs(me) * fr, std(nu_; corrected = true))

end

function _get_signal_to_noise_ratio(nu1_, nu2_)

    me1 = mean(nu1_)

    me2 = mean(nu2_)

    (me1 - me2) / (_get_standard_deviation(nu1_, me1) + _get_standard_deviation(nu2_, me2))

end

# TODO: Benchmark.
function _target_sort(fu, is_, fe_x_sa_x_sc, fe_)

    # TODO: Try `view`.
    sc_ = fu.(eachrow(fe_x_sa_x_sc[:, .!is_]), eachrow(fe_x_sa_x_sc[:, is_]))

    so_ = sortperm(sc_; rev = true)

    view(sc_, so_), view(fe_, so_)

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
  - `--metric`: = "signal-to-noise-ratio". "signal-to-noise-ratio" | (coming soon).
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
#@cast function metric_rank(
function metric_rank(
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

    BioLab.Error.error_missing(output_directory)

    _nat, ta_, sat_, ta_x_sa_x_nu =
        BioLab.DataFrame.separate(BioLab.DataFrame.read(target_x_sample_x_number_tsv))

    BioLab.Error.error_bad(ta_)

    BioLab.Error.error_duplicate(ta_)

    BioLab.Error.error_bad(ta_x_sa_x_nu)

    un_ = Set(ta_x_sa_x_nu)

    if un_ != Set((0, 1))

        error("Target has a number other than 0 and 1. $un_.")

    end

    _naf, fe_, saf_, fe_x_sa_x_sc =
        BioLab.DataFrame.separate(BioLab.DataFrame.read(feature_x_sample_x_score_tsv))

    BioLab.Error.error_bad(fe_)

    BioLab.Error.error_duplicate(fe_)

    BioLab.Error.error_bad(fe_x_sa_x_sc)

    if !iszero(normalization_dimension)

        _normalize_with_0_clamp!(
            fe_x_sa_x_sc,
            normalization_dimension,
            normalization_standard_deviation,
        )

    end

    fe_x_sa_x_sc = view(fe_x_sa_x_sc, :, indexin(sat_, saf_))

    if metric == "signal-to-noise-ratio"

        fu = _get_signal_to_noise_ratio

    else

        error("Metric \"$metric\" is not \"signal-to-noise-ratio\".")

    end

    is_ = convert(BitVector, view(ta_x_sa_x_nu, 1, :))

    sc_, fe_ = _target_sort(fu, is_, fe_x_sa_x_sc, fe_)

    BioLab.DataFrame.write(
        joinpath(output_directory, "feature_x_metric_x_score.tsv"),
        BioLab.DataFrame.make("Feature", fe_, [metric], reshape(sc_, :, 1)),
    )

    se_, fe1___ =
        _read_set(set_features_json, fe_, minimum_set_size, maximum_set_size, minimum_set_fraction)

    al = _set_algorithm(algorithm)

    if permutation == "sample"

        se_x_id_x_ra = Matrix{Float64}(undef, length(fe1___), number_of_permutations)

        if 0 < number_of_permutations

            @info "Permuting samples to compute significance"

            seed!(random_seed)

            @showprogress for id in 1:number_of_permutations

                ra_, fe2_ = _target_sort(fu, shuffle!(is_), fe_x_sa_x_sc, fe_)

                se_x_id_x_ra[:, id] =
                    FeatureSetEnrichment.enrich(al, fe2_, ra_, fe1___; ex = exponent)

            end

        end

    elseif permutation == "set"

        se_x_id_x_ra =
            _permute_set(number_of_permutations, random_seed, al, fe_, sc_, fe1___, exponent)

    elseif isfile(permutation)

        se_x_id_x_ra = _use_permutation(permutation, al, fe_, fe1___, exponent)

    else

        error(
            "Permutation \"$permutation\" is not \"sample\", \"set\", or feature_x_index_x_random.tsv.",
        )

    end

    _write(
        output_directory,
        write_set_x_index_x_random_tsv,
        se_,
        FeatureSetEnrichment.enrich(al, fe_, sc_, fe1___; ex = exponent),
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
#@main

end
