module GSEA

using Comonicon: @cast, @main

using DataFrames: DataFrame, insertcols!

using ProgressMeter: @showprogress

using Random: seed!, shuffle!

using StatsBase: sample

using BioLab

function _filter_set!(se_fe_, it, it_, mi, ma)

    println("⛵️ Before filtering sets")

    BioLab.Dict.print(se_fe_; n = 0)

    if it

        println("🐡 Removing non-intersecting genes")

        for (se, fe_) in se_fe_

            se_fe_[se] = intersect(fe_, it_)

        end

    end

    println("🎣 Removing sets whose size is not between $mi and $ma")

    for (se, fe_) in se_fe_

        if !(mi <= length(fe_) <= ma)

            pop!(se_fe_, se)

        end

    end

    println("🍣 After")

    BioLab.Dict.print(se_fe_; n = 0)

    return nothing

end

function _use_algorithm(al)

    if al == "ks"

        return BioLab.FeatureSetEnrichment.KS()

    elseif al == "ksa"

        return BioLab.FeatureSetEnrichment.KSa()

    elseif al == "kli"

        return BioLab.FeatureSetEnrichment.KLi()

    elseif al == "kliop"

        return BioLab.FeatureSetEnrichment.KLioP()

    elseif al == "kliom"

        return BioLab.FeatureSetEnrichment.KLioM()

    else

        error()

    end

end

# TODO: Benchmark if the lower-level functions are getting the correct type.

"""
Run data-rank (single-sample) GSEA.

# Arguments

  - `setting_json`:
  - `gene_x_sample_x_score_tsv`:
  - `set_genes_json`:
  - `output_directory`:
"""
@cast function data_rank(setting_json, gene_x_sample_x_score_tsv, set_genes_json, output_directory)

    ke_ar = BioLab.Dict.read(setting_json)

    fe_x_sa_x_sc = BioLab.Table.read(gene_x_sample_x_score_tsv)

    se_fe_ = convert(Dict{String, Vector{String}}, BioLab.Dict.read(set_genes_json))

    _filter_set!(
        se_fe_,
        ke_ar["remove_gene_set_genes"],
        fe_x_sa_x_sc[!, 1],
        ke_ar["minimum_gene_set_size"],
        ke_ar["maximum_gene_set_size"],
    )

    se_x_sa_x_en = BioLab.FeatureSetEnrichment.score_set(
        _use_algorithm(ke_ar["algorithm"]),
        fe_x_sa_x_sc,
        se_fe_;
        ex = ke_ar["exponent"],
        n_jo = ke_ar["number_of_jobs"],
    )

    BioLab.Table.write(
        joinpath(mkpath(output_directory), "set_x_sample_x_enrichment.tsv"),
        se_x_sa_x_en,
    )

    return nothing

end

function _tabulate_statistic(se_en, se_ra__, ou)

    se_ = collect(keys(se_en))

    en_ = collect(values(se_en))

    n = length(se_)

    mkpath(ou)

    if isempty(se_ra__)

        gl_ = fill(NaN, n)

        gla_ = fill(NaN, n)

    else

        ra__ = [collect(values(se_ra)) for se_ra in se_ra__]

        gl_, gla_ = BioLab.Significance.get_p_value_and_adjust(en_, vcat(ra__...))

        se_x_ra_x_en = DataFrame("Set" => se_)

        insertcols!(se_x_ra_x_en, (string(id) => ra_ for (id, ra_) in enumerate(ra__))...)

        BioLab.Table.write(joinpath(ou, "set_x_random_x_enrichment.tsv"), se_x_ra_x_en)

    end

    se_x_st_x_nu = sort(
        DataFrame(
            "Set" => se_,
            "Enrichment" => en_,
            "Global p value" => gl_,
            "Adjusted global p value" => gla_,
        ),
        "Enrichment",
    )

    BioLab.Table.write(joinpath(ou, "set_x_statistic_x_number.tsv"), se_x_st_x_nu)

    return se_x_st_x_nu

end

function _plot_mountain(se_x_st_x_nu, fe, sc, lo, hi, n_ex, pl_, al, fe_, sc_, se_fe_, ex, di)

    n_se = size(se_x_st_x_nu, 1)

    n_ex = min(n_ex, n_se)

    # TODO: Try `1:2`.
    co_ = [1, 2]

    for ro in 1:n_ex

        se, en = se_x_st_x_nu[ro, co_]

        if en <= 0 && !(se in pl_)

            push!(pl_, se)

        end

    end

    for ro in n_se:-1:(n_se - n_ex + 1)

        se, en = se_x_st_x_nu[ro, co_]

        if 0 <= en && !(se in pl_)

            push!(pl_, se)

        end

    end

    pl = mkpath(joinpath(di, "plot"))

    for se in pl_

        BioLab.FeatureSetEnrichment.score_set(
            al,
            fe_,
            sc_,
            se_fe_[se];
            ex,
            title_text = se,
            fe,
            sc,
            lo,
            hi,
            ht = joinpath(pl, "$(BioLab.Path.clean(se)).html"),
        )

    end

    return nothing

end

function user_rank(al, fe_, sc_, se_fe_, fe, sc, lo, hi, ex, ra, n_pe, n_ex, pl_, ou)

    se_en = BioLab.FeatureSetEnrichment.score_set(al, fe_, sc_, se_fe_; ex)

    if 0 < n_pe

        println("Permuting sets to compute significance")

        se_si = Dict(se => length(fe_) for (se, fe_) in se_fe_)

        seed!(ra)

        se_ra__ = @showprogress [
            BioLab.FeatureSetEnrichment.score_set(
                al,
                fe_,
                sc_,
                Dict(se => sample(fe_, si; replace = false) for (se, si) in se_si);
                ex,
            ) for _ in 1:n_pe
        ]

    else

        se_ra__ = []

    end

    se_x_st_x_nu = _tabulate_statistic(se_en, se_ra__, ou)

    _plot_mountain(se_x_st_x_nu, fe, sc, lo, hi, n_ex, pl_, al, fe_, sc_, se_fe_, ex, ou)

    return se_x_st_x_nu

end

"""
Run user-rank (pre-rank) GSEA.

# Arguments

  - `setting_json`:
  - `gene_x_metric_x_score_tsv`:
  - `set_genes_json`:
  - `output_directory`:
"""
@cast function user_rank(setting_json, gene_x_metric_x_score_tsv, set_genes_json, output_directory)

    ke_ar = BioLab.Dict.read(setting_json)

    _nag, fe_, _me, fe_x_me_x_sc =
        BioLab.DataFrame.separate(BioLab.Table.read(gene_x_metric_x_score_tsv))

    BioLab.Array.error_duplicate(fe_)

    BioLab.Matrix.error_bad(fe_x_me_x_sc, Real)

    sc_ = fe_x_me_x_sc[:, 1]

    sc_, fe_ = BioLab.Collection.sort_like((sc_, fe_); ic = false)

    se_fe_ = convert(Dict{String, Vector{String}}, BioLab.Dict.read(set_genes_json))

    _filter_set!(
        se_fe_,
        ke_ar["remove_gene_set_genes"],
        fe_,
        ke_ar["minimum_gene_set_size"],
        ke_ar["maximum_gene_set_size"],
    )

    user_rank(
        _use_algorithm(ke_ar["algorithm"]),
        fe_,
        sc_,
        se_fe_,
        ke_ar["feature_name"],
        ke_ar["score_name"],
        ke_ar["low_text"],
        ke_ar["high_text"],
        ke_ar["exponent"],
        ke_ar["random_seed"],
        ke_ar["number_of_permutations"],
        ke_ar["number_of_extreme_gene_sets_to_plot"],
        convert(Vector{String}, ke_ar["gene_sets_to_plot"]),
        output_directory,
    )

    return nothing

end

function _compare_and_sort(fu, bo_, fe_x_sa_x_sc, fe_)

    sc_, fes_ = BioLab.Collection.sort_like(
        (BioLab.FeatureXSample.target(fu, bo_, fe_x_sa_x_sc), fe_);
        ic = false,
    )

    return fes_, sc_

end

"""
Run metric-rank (standard) GSEA.

# Arguments

  - `setting_json`:
  - `target_x_sample_x_number_tsv`:
  - `gene_x_sample_x_score_tsv`:
  - `set_genes_json`:
  - `output_directory`:
"""
@cast function metric_rank(
    setting_json,
    target_x_sample_x_number_tsv,
    gene_x_sample_x_score_tsv,
    set_genes_json,
    output_directory,
)

    ke_ar = BioLab.Dict.read(setting_json)

    _nat, ta_, sat_, ta_x_sa_x_nu =
        BioLab.DataFrame.separate(BioLab.Table.read(target_x_sample_x_number_tsv))

    BioLab.Array.error_duplicate(ta_)

    BioLab.Matrix.error_bad(ta_x_sa_x_nu, Real)

    _nag, fe_, saf_, fe_x_sa_x_sc =
        BioLab.DataFrame.separate(BioLab.Table.read(gene_x_sample_x_score_tsv))

    BioLab.Array.error_duplicate(fe_)

    BioLab.Matrix.error_bad(fe_x_sa_x_sc, Real)

    fe_x_sa_x_sc = fe_x_sa_x_sc[:, indexin(sat_, saf_)]

    mkpath(output_directory)

    bo_ = convert(Vector{Bool}, ta_x_sa_x_nu[1, :])

    me = ke_ar["metric"]

    if me == "signal_to_noise_ratio"

        fu = BioLab.Information.get_signal_to_noise_ratio

    else

        error()

    end

    fe_, sc_ = _compare_and_sort(fu, bo_, fe_x_sa_x_sc, fe_)

    BioLab.Table.write(
        joinpath(output_directory, "gene_x_metric_x_score.tsv"),
        DataFrame("Gene" => fe_, me => sc_),
    )

    se_fe_ = convert(Dict{String, Vector{String}}, BioLab.Dict.read(set_genes_json))

    _filter_set!(
        se_fe_,
        ke_ar["remove_gene_set_genes"],
        fe_,
        ke_ar["minimum_gene_set_size"],
        ke_ar["maximum_gene_set_size"],
    )

    al = _use_algorithm(ke_ar["algorithm"])

    fe = ke_ar["feature_name"]

    sc = ke_ar["score_name"]

    lo = ke_ar["low_text"]

    hi = ke_ar["high_text"]

    ex = ke_ar["exponent"]

    pe = ke_ar["permutation"]

    ra = ke_ar["random_seed"]

    n_pe = ke_ar["number_of_permutations"]

    n_ex = ke_ar["number_of_extreme_gene_sets_to_plot"]

    pl_ = ke_ar["gene_sets_to_plot"]

    if pe == "sample"

        se_en = BioLab.FeatureSetEnrichment.score_set(al, fe_, sc_, se_fe_; ex)

        if 0 < n_pe

            println("Permuting $(pe)s to compute significance")

            seed!(ra)

            se_ra__ = @showprogress [
                BioLab.FeatureSetEnrichment.score_set(
                    al,
                    _compare_and_sort(fu, shuffle!(bo_), fe_x_sa_x_sc, fe_)...,
                    se_fe_;
                    ex,
                ) for _ in 1:n_pe
            ]

        else

            se_ra__ = []

        end

        se_x_st_x_nu = _tabulate_statistic(se_en, se_ra__, output_directory)

        _plot_mountain(
            se_x_st_x_nu,
            fe,
            sc,
            lo,
            hi,
            n_ex,
            pl_,
            al,
            fe_,
            sc_,
            se_fe_,
            ex,
            output_directory,
        )

        se_x_st_x_nu

    elseif pe == "set"

        user_rank(al, fe_, sc_, se_fe_, fe, sc, lo, hi, ex, ra, n_pe, n_ex, pl_, output_directory)

    else

        error("`permutation` is not `sample` or `set`.")

    end

    return nothing

end

"""
The ✨ new ✨ Gene-Set Enrichment Analysis 🧬
"""
@main

end
