module GSEA

using Comonicon: @cast, @main

using ProgressMeter: @showprogress

using Random: seed!, shuffle!

using StatsBase: mean, sample, std

using BioLab

function _read_set(js, it_, mi, ma)

    se_fe1_ = BioLab.Dict.read(js)

    se_ = collect(keys(se_fe1_))

    fe1___ = [convert(Vector{String}, fe1_) for fe1_ in values(se_fe1_)]

    n = length(se_)

    println("🎣 There are $n sets before filtering.")

    ke_ = Vector{Bool}(undef, n)

    for (id, fe1_) in enumerate(fe1___)

        intersect!(fe1_, it_)

        ke_[id] = mi <= length(fe1_) <= ma

    end

    println("🍣 $(sum(ke_)) after.")

    se_[ke_], fe1___[ke_]

end

function _use_algorithm(al)

    if al == "ks"

        BioLab.FeatureSetEnrichment.KS()

    elseif al == "ksa"

        BioLab.FeatureSetEnrichment.KSa()

    elseif al == "kli"

        BioLab.FeatureSetEnrichment.KLi()

    elseif al == "kliop"

        BioLab.FeatureSetEnrichment.KLioP()

    elseif al == "kliom"

        BioLab.FeatureSetEnrichment.KLioM()

    else

        error("`algorithm` is not one of the listed in https://github.com/KwatMDPhD/GSEA.jl.")

    end

end

# TODO: Benchmark if the lower-level functions are getting the correct type.

"""
Run data-rank (single-sample) GSEA.

# Arguments

  - `setting_json`:
  - `feature_x_sample_x_score_tsv`:
  - `set_features_json`:
  - `output_directory`:
"""
@cast function data_rank(
    setting_json,
    feature_x_sample_x_score_tsv,
    set_features_json,
    output_directory,
)

    ke_ar = BioLab.Dict.read(setting_json)

    _fen, fe_::Vector{String}, sa_::Vector{String}, fe_x_sa_x_sc::Matrix{Float64} =
        BioLab.DataFrame.separate(BioLab.Table.read(feature_x_sample_x_score_tsv))

    BioLab.Array.error_duplicate(fe_)

    BioLab.Collection.error_bad(fe_x_sa_x_sc, Real)

    se_, fe1___ =
        _read_set(set_features_json, fe_, ke_ar["minimum_set_size"], ke_ar["maximum_set_size"])

    se_x_sa_x_en = BioLab.FeatureSetEnrichment.enrich(
        _use_algorithm(ke_ar["algorithm"]),
        fe_,
        sa_,
        fe_x_sa_x_sc,
        se_,
        fe1___;
        ex = ke_ar["exponent"],
        n_jo = ke_ar["number_of_jobs"],
    )

    BioLab.Table.write(
        joinpath(mkpath(output_directory), "set_x_sample_x_enrichment.tsv"),
        BioLab.DataFrame.make("Set", se_, sa_, se_x_sa_x_en),
    )

    nothing

end

function _write(
    se_,
    en_,
    se_x_id_x_ra,
    al,
    fe_,
    sc_,
    fe1___,
    ex,
    fe,
    sc,
    lo,
    hi,
    n_pl,
    pl_,
    ou,
    wr,
)

    id_ = sortperm(en_)

    en_ = en_[id_]

    se_ = se_[id_]

    se_x_id_x_ra = se_x_id_x_ra[id_, :]

    n_se = length(se_)

    se_x_st_x_nu = fill(NaN, n_se, 4)

    se_x_st_x_nu[:, 1] = en_

    n_ra = size(se_x_id_x_ra, 2)

    mkpath(ou)

    if !isempty(se_x_id_x_ra)

        if wr

            BioLab.Table.write(
                joinpath(ou, "set_x_index_x_random.tsv"),
                BioLab.DataFrame.make("Set", se_, [string(id) for id in 1:n_ra], se_x_id_x_ra),
            )

        end

        nem_ = Vector{Float64}(undef, n_se)

        pom_ = Vector{Float64}(undef, n_se)

        for id in 1:n_se

            ne_ = Vector{Float64}()

            po_ = Vector{Float64}()

            for ra in se_x_id_x_ra[id, :]

                if ra < 0.0

                    push!(ne_, ra)

                else

                    push!(po_, ra)

                end

            end

            nem_[id] = mean(ne_)

            pom_[id] = mean(po_)

        end

        nei_ = 1:findlast(en < 0.0 for en in en_)

        poi_ = (nei_[end] + 1):n_se

        # TODO: Benchmark against `vcat`.

        enn_ = Vector{Float64}(undef, n_se)

        for id in nei_

            enn_[id] = -en_[id] / nem_[id]

        end

        for id in poi_

            enn_[id] = en_[id] / pom_[id]

        end

        se_x_st_x_nu[:, 2] = enn_

        nen_ = Vector{Float64}()

        pon_ = Vector{Float64}()

        # TODO: Benchmark iteration order.

        for id2 in 1:n_ra

            for id1 in 1:n_se

                ra = se_x_id_x_ra[id1, id2]

                if ra < 0.0

                    push!(nen_, -ra / nem_[id1])

                elseif 0.0 < ra

                    push!(pon_, ra / pom_[id1])

                end

            end

        end

        nep_, nea_ = BioLab.Significance.get_p_value_and_adjust(
            BioLab.Significance.get_p_value_for_less,
            enn_[nei_],
            nen_,
        )

        se_x_st_x_nu[nei_, 3] = nep_

        se_x_st_x_nu[nei_, 4] = nea_

        pop_, poa_ = BioLab.Significance.get_p_value_and_adjust(
            BioLab.Significance.get_p_value_for_more,
            enn_[poi_],
            pon_,
        )

        se_x_st_x_nu[poi_, 3] = pop_

        se_x_st_x_nu[poi_, 4] = poa_

    end

    BioLab.Table.write(
        joinpath(ou, "set_x_statistic_x_number.tsv"),
        BioLab.DataFrame.make(
            "Set",
            se_,
            ["Enrichment", "Normalized Enrichment", "P Value", "Adjusted P Value"],
            se_x_st_x_nu,
        ),
    )

    n_pl = min(n_pl, n_se)

    for id in 1:n_pl

        se = se_[id]

        en = en_[id]

        if en < 0 && !(se in pl_)

            push!(pl_, se)

        end

    end

    for id in n_se:-1:(n_se - n_pl + 1)

        se = se_[id]

        en = en_[id]

        if 0 < en && !(se in pl_)

            push!(pl_, se)

        end

    end

    oup = mkpath(joinpath(ou, "plot"))

    for (se, id) in zip(pl_, indexin(pl_, se_))

        BioLab.FeatureSetEnrichment.enrich(
            al,
            fe_,
            sc_,
            fe1___[id];
            ex,
            title_text = se,
            fe,
            sc,
            lo,
            hi,
            ht = joinpath(oup, "$(BioLab.Path.clean(se)).html"),
        )

    end

    nothing

end

function user_rank(al, fe_, sc_, se_, fe1___, ex, fe, sc, lo, hi, ra, n_pe, n_pl, pl_, ou, wr)

    en_ = BioLab.FeatureSetEnrichment.enrich(al, fe_, sc_, fe1___; ex)

    se_x_id_x_ra = Matrix{Float64}(undef, length(se_), n_pe)

    if 0 < n_pe

        println("Permuting sets to compute significance")

        si_ = [length(fe1_) for fe1_ in fe1___]

        seed!(ra)

        @showprogress for id in 1:n_pe

            se_x_id_x_ra[:, id] = BioLab.FeatureSetEnrichment.enrich(
                al,
                fe_,
                sc_,
                [sample(fe_, si; replace = false) for si in si_];
                ex,
            )

        end

    end

    _write(se_, en_, se_x_id_x_ra, al, fe_, sc_, fe1___, ex, fe, sc, lo, hi, n_pl, pl_, ou, wr)

end

"""
Run user-rank (pre-rank) GSEA.

# Arguments

  - `setting_json`:
  - `feature_x_metric_x_score_tsv`:
  - `set_features_json`:
  - `output_directory`:
"""
@cast function user_rank(
    setting_json,
    feature_x_metric_x_score_tsv,
    set_features_json,
    output_directory,
)

    ke_ar = BioLab.Dict.read(setting_json)

    _fen, fe_, _me_, fe_x_me_x_sc =
        BioLab.DataFrame.separate(BioLab.Table.read(feature_x_metric_x_score_tsv))

    BioLab.Array.error_duplicate(fe_)

    # TODO
    BioLab.Collection.error_bad(fe_x_me_x_sc[:, [1]], Real)

    sc_ = fe_x_me_x_sc[:, 1]

    sc_, fe_ = BioLab.Collection.sort_like((sc_, fe_); ic = false)

    se_, fe1___ =
        _read_set(set_features_json, fe_, ke_ar["minimum_set_size"], ke_ar["maximum_set_size"])

    user_rank(
        _use_algorithm(ke_ar["algorithm"]),
        fe_,
        sc_,
        se_,
        fe1___,
        ke_ar["exponent"],
        ke_ar["feature_name"],
        ke_ar["score_name"],
        ke_ar["low_text"],
        ke_ar["high_text"],
        ke_ar["random_seed"],
        ke_ar["number_of_permutations"],
        ke_ar["number_of_sets_to_plot"],
        ke_ar["more_sets_to_plot"],
        output_directory,
        ke_ar["write_set_x_index_x_random_tsv"],
    )

end

function _get_standard_deviation(nu_, me)

    fr = 0.2

    if me == 0.0

        fr

    else

        if me < 0.0

            me = -me

        end

        max(me * fr, std(nu_; corrected = true))

    end

end

function _get_signal_to_noise_ratio(nu1_, nu2_)

    me1 = mean(nu1_)

    me2 = mean(nu2_)

    (me1 - me2) / (_get_standard_deviation(nu1_, me1) + _get_standard_deviation(nu2_, me2))

end

function _compare_and_sort(fu, bo_, fe_x_sa_x_sc, fe_)

    BioLab.Collection.sort_like(
        (BioLab.FeatureXSample.target(fu, bo_, fe_x_sa_x_sc), fe_);
        ic = false,
    )

end

"""
Run metric-rank (standard) GSEA.

# Arguments

  - `setting_json`:
  - `target_x_sample_x_number_tsv`:
  - `feature_x_sample_x_score_tsv`:
  - `set_features_json`:
  - `output_directory`:
"""
@cast function metric_rank(
    setting_json,
    target_x_sample_x_number_tsv,
    feature_x_sample_x_score_tsv,
    set_features_json,
    output_directory;
    feature2_x_index_x_random = nothing,
)

    ke_ar = BioLab.Dict.read(setting_json)

    _tan, ta_, sat_, ta_x_sa_x_nu =
        BioLab.DataFrame.separate(BioLab.Table.read(target_x_sample_x_number_tsv))

    BioLab.Array.error_duplicate(ta_)

    BioLab.Collection.error_bad(ta_x_sa_x_nu, Real)

    _fen, fe_, saf_, fe_x_sa_x_sc =
        BioLab.DataFrame.separate(BioLab.Table.read(feature_x_sample_x_score_tsv))

    BioLab.Array.error_duplicate(fe_)

    BioLab.Collection.error_bad(fe_x_sa_x_sc, Real)

    fe_x_sa_x_sc = fe_x_sa_x_sc[:, indexin(sat_, saf_)]

    mkpath(output_directory)

    bo_ = [nu == 0 for nu in ta_x_sa_x_nu[1, :]]

    me = ke_ar["metric"]

    if me == "signal_to_noise_ratio"

        fu = _get_signal_to_noise_ratio

    else

        error("`metric` is not one of the listed in https://github.com/KwatMDPhD/GSEA.jl.")

    end

    sc_, fe_ = _compare_and_sort(fu, bo_, fe_x_sa_x_sc, fe_)

    BioLab.Table.write(
        joinpath(output_directory, "feature_x_metric_x_score.tsv"),
        BioLab.DataFrame.make("Feature", fe_, [me], hcat(sc_)),
    )

    se_, fe1___ =
        _read_set(set_features_json, fe_, ke_ar["minimum_set_size"], ke_ar["maximum_set_size"])

    al = _use_algorithm(ke_ar["algorithm"])

    ex = ke_ar["exponent"]

    fe = ke_ar["feature_name"]

    sc = ke_ar["score_name"]

    lo = ke_ar["low_text"]

    hi = ke_ar["high_text"]

    pe = ke_ar["permutation"]

    ra = ke_ar["random_seed"]

    n_pe = ke_ar["number_of_permutations"]

    n_pl = ke_ar["number_of_sets_to_plot"]

    pl_ = ke_ar["more_sets_to_plot"]

    wr = ke_ar["write_set_x_index_x_random_tsv"]

    if pe == "sample"

        en_ = BioLab.FeatureSetEnrichment.enrich(al, fe_, sc_, fe1___; ex)

        se_x_id_x_ra = Matrix{Float64}(undef, length(se_), n_pe)

        if !isnothing(feature2_x_index_x_random)

            _fe2n, fe2_, _id_, fe2_x_id_x_ra = BioLab.DataFrame.separate(feature2_x_index_x_random)

            println("🎰 Using predefined $n_pe $pe permutations to compute significance")

            @showprogress for id in 1:n_pe

                ra_, fer_ = BioLab.Collection.sort_like((fe2_x_id_x_ra[:, id], fe2_); ic = false)

                se_x_id_x_ra[:, id] = BioLab.FeatureSetEnrichment.enrich(al, fer_, ra_, fe1___; ex)

            end

        elseif 0 < n_pe

            println("🎰 Permuting $(pe)s $n_pe times to compute significance")

            seed!(ra)

            @showprogress for id in 1:n_pe

                ra_, fer_ = _compare_and_sort(fu, shuffle!(bo_), fe_x_sa_x_sc, fe_)

                se_x_id_x_ra[:, id] = BioLab.FeatureSetEnrichment.enrich(al, fer_, ra_, fe1___; ex)

            end

        end

        _write(
            se_,
            en_,
            se_x_id_x_ra,
            al,
            fe_,
            sc_,
            fe1___,
            ex,
            fe,
            sc,
            lo,
            hi,
            n_pl,
            pl_,
            output_directory,
            wr,
        )

    elseif pe == "set"

        user_rank(
            al,
            fe_,
            sc_,
            se_,
            fe1___,
            ex,
            fe,
            sc,
            lo,
            hi,
            ra,
            n_pe,
            n_pl,
            pl_,
            output_directory,
            wr,
        )

    else

        error("`permutation` is not one of the listed in https://github.com/KwatMDPhD/GSEA.jl.")

    end

end

"""
🧬 The ✨ new ✨ (Gene-) Set Enrichment Analysis. 📍 Learn more at https://github.com/KwatMDPhD/GSEA.jl.
"""
@main

end
