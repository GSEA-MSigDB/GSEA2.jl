module GSEA

using Comonicon: @cast, @main

using DataFrames: DataFrame, insertcols!

using ProgressMeter: @showprogress

using Random: seed!, shuffle!

using StatsBase: sample

using BioLab

function _filter_set!(se_fe_, re, it_, mi, ma)

    println("🕺 Before filtering sets")

    BioLab.Dict.print(se_fe_, 0)

    if re

        println("Removing set genes not found in gene-x-sample genes")

        for (se, fe_) in se_fe_

            se_fe_[se] = intersect(fe_, it_)

        end

    end

    println("🎣 Keeping sets: $mi <= size <= $ma")

    for (se, fe_) in se_fe_

        if !(mi <= length(fe_) <= ma)

            pop!(se_fe_, se)

        end

    end

    println("💃 After")

    BioLab.Dict.print(se_fe_, 0)

end

function _make_keyword_argument(ke_ar)

    Dict(
        Symbol(kes) => ke_ar[ke] for (kes, ke) in (("ex", "exponent"), ("n_jo", "number_of_jobs"))
    )

end

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

    se_fe_ = BioLab.Dict.read(set_genes_json)

    _filter_set!(
        se_fe_,
        ke_ar["remove_gene_set_genes"],
        fe_x_sa_x_sc[!, 1],
        ke_ar["minimum_gene_set_size"],
        ke_ar["maximum_gene_set_size"],
    )

    se_x_sa_x_en = BioLab.FeatureSetEnrichment.score_set(
        ke_ar["algorithm"],
        fe_x_sa_x_sc,
        se_fe_;
        _make_keyword_argument(ke_ar)...,
    )

    BioLab.Table.write(
        joinpath(mkpath(output_directory), "set_x_sample_x_enrichment.tsv"),
        se_x_sa_x_en,
    )

end

function _tabulate_statistic(se_en, se_ra_, ou)

    se_ = collect(keys(se_en))

    en_ = collect(values(se_en))

    mkpath(ou)

    if isempty(se_ra_)

        gl_ = gla_ = fill(NaN, length(se_))

    else

        ra__ = [collect(values(se_ra)) for se_ra in se_ra_]

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

    se_x_st_x_nu

end

function _plot_mountain(se_x_st_x_nu, fe, sc, n_ex, pl_, al, fe_, sc_, se_fe_, sy_ar, di)

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

    pop!(sy_ar, :n_jo)

    for se in pl_

        BioLab.FeatureSetEnrichment.score_set(
            al,
            fe_,
            sc_,
            se_fe_[se];
            title_text = se,
            fe = fe,
            sc = sc,
            ht = joinpath(pl, "$(BioLab.Path.clean(se)).html"),
            sy_ar...,
        )

    end

end

function user_rank(al, fe_, sc_, se_fe_, fe, sc, sy_ar, ra, n_pe, n_ex, pl_, ou)

    se_en = BioLab.FeatureSetEnrichment.score_set(al, fe_, sc_, se_fe_; sy_ar...)

    if 0 < n_pe

        println("Permuting sets to compute significance")

        se_si = Dict(se => length(fe_) for (se, fe_) in se_fe_)

        seed!(ra)

        se_ra_ = @showprogress [
            BioLab.FeatureSetEnrichment.score_set(
                al,
                fe_,
                sc_,
                Dict(se => sample(fe_, si, replace = false) for (se, si) in se_si);
                sy_ar...,
            ) for _ in 1:n_pe
        ]

    else

        se_ra_ = []

    end

    se_x_st_x_nu = _tabulate_statistic(se_en, se_ra_, ou)

    _plot_mountain(se_x_st_x_nu, fe, sc, n_ex, pl_, al, fe_, sc_, se_fe_, sy_ar, ou)

    se_x_st_x_nu

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

    fe_, fe_x_me_x_sc =
        BioLab.DataFrame.separate(BioLab.Table.read(gene_x_metric_x_score_tsv))[[2, 4]]

    BioLab.Array.error_duplicate(fe_)

    BioLab.Matrix.error_bad(fe_x_me_x_sc, Real)

    sc_ = fe_x_me_x_sc[:, 1]

    sc_, fe_ = BioLab.Vector.sort_like((sc_, fe_))

    se_fe_ = BioLab.Dict.read(set_genes_json)

    _filter_set!(
        se_fe_,
        ke_ar["remove_gene_set_genes"],
        fe_,
        ke_ar["minimum_gene_set_size"],
        ke_ar["maximum_gene_set_size"],
    )

    user_rank(
        ke_ar["algorithm"],
        fe_,
        sc_,
        se_fe_,
        ke_ar["feature_name"],
        ke_ar["score_name"],
        _make_keyword_argument(ke_ar),
        ke_ar["random_seed"],
        ke_ar["number_of_permutations"],
        ke_ar["number_of_extreme_gene_sets_to_plot"],
        ke_ar["gene_sets_to_plot"],
        output_directory,
    )

end

function _compare_and_sort(bo_, fe_x_sa_x_sc, me, fe_)

    sc_, fes_ = BioLab.Vector.sort_like((BioLab.FeatureXSample.target(bo_, fe_x_sa_x_sc, me), fe_))

    fes_, sc_

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

    ta_, sat_, ta_x_sa_x_nu =
        BioLab.DataFrame.separate(BioLab.Table.read(target_x_sample_x_number_tsv))[[2, 3, 4]]

    BioLab.Array.error_duplicate(ta_)

    BioLab.Matrix.error_bad(ta_x_sa_x_nu, Real)

    fe_, saf_, fe_x_sa_x_sc =
        BioLab.DataFrame.separate(BioLab.Table.read(gene_x_sample_x_score_tsv))[[2, 3, 4]]

    BioLab.Array.error_duplicate(fe_)

    BioLab.Matrix.error_bad(fe_x_sa_x_sc, Real)

    fe_x_sa_x_sc = fe_x_sa_x_sc[:, indexin(sat_, saf_)]

    mkpath(output_directory)

    bo_ = convert(Vector{Bool}, ta_x_sa_x_nu[1, :])

    me = ke_ar["metric"]

    fe_, sc_ = _compare_and_sort(bo_, fe_x_sa_x_sc, me, fe_)

    BioLab.Table.write(
        joinpath(output_directory, "gene_x_metric_x_score.tsv"),
        DataFrame("Gene" => fe_, me => sc_),
    )

    se_fe_ = BioLab.Dict.read(set_genes_json)

    _filter_set!(
        se_fe_,
        ke_ar["remove_gene_set_genes"],
        fe_,
        ke_ar["minimum_gene_set_size"],
        ke_ar["maximum_gene_set_size"],
    )

    al = ke_ar["algorithm"]

    fe = ke_ar["feature_name"]

    sc = ke_ar["score_name"]

    sy_ar = _make_keyword_argument(ke_ar)

    pe = ke_ar["permutation"]

    ra = ke_ar["random_seed"]

    n_pe = ke_ar["number_of_permutations"]

    n_ex = ke_ar["number_of_extreme_gene_sets_to_plot"]

    pl_ = ke_ar["gene_sets_to_plot"]

    if pe == "sample"

        se_en = BioLab.FeatureSetEnrichment.score_set(al, fe_, sc_, se_fe_; sy_ar...)

        if 0 < n_pe

            println("Permuting $(pe)s to compute significance")

            seed!(ra)

            se_ra_ = @showprogress [
                BioLab.FeatureSetEnrichment.score_set(
                    al,
                    _compare_and_sort(shuffle!(bo_), fe_x_sa_x_sc, me, fe_)...,
                    se_fe_;
                    sy_ar...,
                ) for _ in 1:n_pe
            ]

        else

            se_ra_ = []

        end

        se_x_st_x_nu = _tabulate_statistic(se_en, se_ra_, output_directory)

        _plot_mountain(
            se_x_st_x_nu,
            fe,
            sc,
            n_ex,
            pl_,
            al,
            fe_,
            sc_,
            se_fe_,
            sy_ar,
            output_directory,
        )

        se_x_st_x_nu

    elseif pe == "set"

        user_rank(al, fe_, sc_, se_fe_, fe, sc, sy_ar, ra, n_pe, n_ex, pl_, output_directory)

    else

        error("`permutation` is not `sample` or `set`.")

    end

end

"""
The ✨ new ✨ Gene-Set Enrichment Analysis 🧬
"""
@main

end
