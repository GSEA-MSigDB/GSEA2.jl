module GSEA

# ----------------------------------------------------------------------------------------------- #

using Comonicon: @cast, @main

using Printf: @sprintf

using ProgressMeter: @showprogress

using Random: seed!, shuffle!

using StatsBase: mean, sample, std

using Omics

# TODO: Generalize.
function _exponentiate(nu, ex)

    ab = abs(nu)

    isone(ex) ? ab : ab^ex

end

# TODO: Generalize.
function _is_in(fe_, f1_)

    map(in(Set(f1_)), fe_)

end

# TODO: Generalize.
function _is_in!(ii_, fe_id, f1_)

    for f1 in f1_

        id = get(fe_id, f1, nothing)

        if !isnothing(id)

            ii_[id] = true

        end

    end

end

function read_cls(cl)

    l1, l2, l3 = readlines(cl)

    np = "Phenotype"

    ph = l2[2:end]

    va_ = split(l3)

    us = lastindex(va_)

    sa_ = Omics.Simulation.label(us, "Sample")

    if l1 == "#numeric"

        # TODO: Benchmark reshape.
        Omics.Table.make(np, ph, sa_, [parse(Float64, va) for _ in 1:1, va in va_])

    else

        l1_ = split(l1)

        u1 = parse(Int, l1_[1])

        if u1 != us

            error("numbers of samples differ: $u1 and $us.")

        end

        u2 = parse(Int, l1_[2])

        ph_ = split(ph)

        up = lastindex(ph_)

        uu = lastindex(unique(va_))

        if !(u2 == up == uu)

            error("numbers of groups differ: $u2, $up, and $uu.")

        end

        ph_id = Dict(ph => id for (id, ph) in enumerate(ph_))

        # TODO: Benchmark reshape.
        Omics.Table.make(np, join(ph_, '_'), sa_, [ph_id[va] for _ in 1:1, va in va_])

    end

end

function read_gct(gc)

    Omics.Table.rea(gc; header = 3, drop = ["Description"])

end

function read_gmt(gm)

    se_me_ = Dict{String, Vector{String}}()

    for li in eachline(gm)

        sp_ = split(li, '\t')

        se = sp_[1]

        if haskey(se_me_, se)

            error("there is more than one $se.")

        end

        se_me_[se] = filter!(!isempty, sp_[3:end])

    end

    se_me_

end

"""
Convert .cls to .tsv.

# Arguments

  - `tsv`:
  - `cls`:
"""
@cast function cls(tsv, cls)

    ta = read_cls(cls)

    na_ = names(ta)

    va = Matrix(ta[!, 2:end])

    Omics.Table.writ(
        tsv,
        Omics.Table.make(na_[1], ta[!, 1], na_[2:end], map!(nu -> nu - 1, va, va)),
    )

end

"""
Convert .gct to .tsv.

# Arguments

  - `tsv`:
  - `gct`:
"""
@cast function gct(tsv, gct)

    ta = read_gct(gct)

    Omics.Table.writ(
        tsv,
        Omics.Table.make("Feature", ta[!, 1], names(ta)[2:end], Matrix(ta[!, 2:end])),
    )

end

"""
Merge .gmts into .json.

# Arguments

  - `json`:
  - `gmt_`:
"""
@cast function gmt(json, gmt_...)

    Omics.Dic.writ(json, reduce(merge!, (read_gmt(gm) for gm in gmt_)))

end

struct KS end

struct KSa end

struct KLioM end

struct KLioP end

struct KLi end

struct KLi1 end

function strin(al)

    string(al)[6:(end - 2)]

end

function _get_normalizer(::Union{KS, KSa}, sc_, ex, ii_)

    s0 = s1 = 0.0

    for id in eachindex(sc_)

        if ii_[id]

            s1 += _exponentiate(sc_[id], ex)

        else

            s0 += 1.0

        end

    end

    -inv(s0), inv(s1)

end

function _get_normalizer(::Union{KLioM, KLioP, KLi}, sc_, ex, ii_)

    sa = s1 = 0.0

    for id in eachindex(sc_)

        am = _exponentiate(sc_[id], ex)

        sa += am

        if ii_[id]

            s1 += am

        end

    end

    inv(sa), inv(s1)

end

function _get_normalizer(na, n1)

    inv(inv(na) - inv(n1))

end

function _get_normalizer(::KLi1, sc_, ex, ii_)

    s1 = 0.0

    for id in eachindex(sc_)

        if ii_[id]

            s1 += _exponentiate(sc_[id], ex)

        end

    end

    inv(s1)

end

function _enrich!(al::KS, sc_, ex, ii_, mo_)

    n0, n1 = _get_normalizer(al, sc_, ex, ii_)

    mo = ba = bm = 0.0

    for id in eachindex(sc_)

        mo += ii_[id] ? _exponentiate(sc_[id], ex) * n1 : n0

        if !isnothing(mo_)

            mo_[id] = mo

        end

        ab = abs(mo)

        if ba < ab

            ba = ab

            bm = mo

        end

    end

    bm

end

function _enrich!(al::KSa, sc_, ex, ii_, mo_)

    n0, n1 = _get_normalizer(al, sc_, ex, ii_)

    mo = ar = 0.0

    for id in eachindex(sc_)

        ar += mo += ii_[id] ? _exponentiate(sc_[id], ex) * n1 : n0

        if !isnothing(mo_)

            mo_[id] = mo

        end

    end

    ar / lastindex(sc_)

end

# TODO: Clip.
const ON = 1.0 + 1e-13

function _enrich!(al::KLioM, sc_, ex, ii_, mo_)

    na, n1 = _get_normalizer(al, sc_, ex, ii_)

    n0 = _get_normalizer(na, n1)

    ra = r0 = r1 = eps()

    la = l0 = l1 = ON

    pa = p0 = p1 = ar = 0.0

    for id in eachindex(sc_)

        am = _exponentiate(sc_[id], ex)

        da = am * na

        if ii_[id]

            d0 = 0.0

            d1 = am * n1

        else

            d0 = am * n0

            d1 = 0.0

        end

        ra += da

        r0 += d0

        r1 += d1

        la -= pa

        l0 -= p0

        l1 -= p1

        ar +=
            mo =
                Omics.Information.get_antisymmetric_kullback_leibler_divergence(
                    r1,
                    r0,
                    ra,
                ) -
                Omics.Information.get_antisymmetric_kullback_leibler_divergence(l1, l0, la)

        if !isnothing(mo_)

            mo_[id] = mo

        end

        pa = da

        p0 = d0

        p1 = d1

    end

    ar / lastindex(sc_)

end

function _enrich!(al::KLioP, sc_, ex, ii_, mo_)

    na, n1 = _get_normalizer(al, sc_, ex, ii_)

    n0 = _get_normalizer(na, n1)

    ra = r0 = r1 = eps()

    la = l0 = l1 = ON

    pa = p0 = p1 = ar = 0.0

    for id in eachindex(sc_)

        am = _exponentiate(sc_[id], ex)

        da = am * na

        if ii_[id]

            d0 = 0.0

            d1 = am * n1

        else

            d0 = am * n0

            d1 = 0.0

        end

        ra += da

        r0 += d0

        r1 += d1

        la -= pa

        l0 -= p0

        l1 -= p1

        ar +=
            mo =
                Omics.Information.get_symmetric_kullback_leibler_divergence(r1, r0, ra) -
                Omics.Information.get_symmetric_kullback_leibler_divergence(l1, l0, la)

        if !isnothing(mo_)

            mo_[id] = mo

        end

        pa = da

        p0 = d0

        p1 = d1

    end

    ar / lastindex(sc_)

end

function _enrich!(al::KLi, sc_, ex, ii_, mo_)

    na, n1 = _get_normalizer(al, sc_, ex, ii_)

    ra = r1 = eps()

    la = l1 = ON

    pa = p1 = ar = 0.0

    for id in eachindex(sc_)

        am = _exponentiate(sc_[id], ex)

        da = am * na

        d1 = ii_[id] ? am * n1 : 0.0

        ra += da

        r1 += d1

        la -= pa

        l1 -= p1

        ar +=
            mo = Omics.Information.get_antisymmetric_kullback_leibler_divergence(
                r1,
                l1,
                ra,
                la,
            )

        if !isnothing(mo_)

            mo_[id] = mo

        end

        pa = da

        p1 = d1

    end

    ar / lastindex(sc_)

end

function _enrich!(al::KLi1, sc_, ex, ii_, mo_)

    uf = lastindex(sc_)

    n1 = _get_normalizer(al, sc_, ex, ii_)

    da = inv(uf)

    ra = r1 = eps()

    la = ON + da

    l1 = ON

    p1 = ar = 0.0

    for id in eachindex(sc_)

        d1 = ii_[id] ? _exponentiate(sc_[id], ex) * n1 : 0.0

        ra += da

        r1 += d1

        la -= da

        l1 -= p1

        ar +=
            mo = Omics.Information.get_antisymmetric_kullback_leibler_divergence(
                r1,
                l1,
                ra,
                la,
            )

        if !isnothing(mo_)

            mo_[id] = mo

        end

        p1 = d1

    end

    ar / uf

end

function _get_extreme(mo_)

    mi, ma = extrema(mo_)

    ai = abs(mi)

    aa = abs(ma)

    if isapprox(ai, aa)

        (mi, ma)

    elseif aa < ai

        (mi,)

    else

        (ma,)

    end

end

function plot(
    ht,
    al,
    fe_,
    sc_,
    me_;
    ex = 1.0,
    nf = "Feature",
    ns = "Score",
    nl = "Low",
    nh = "High",
    la = Dict{String, Any}(),
)

    uf = lastindex(fe_)

    xc_ = collect(1:uf)

    ii_ = _is_in(fe_, me_)

    mo_ = Vector{Float64}(undef, uf)

    en = _enrich!(al, sc_, ex, ii_, mo_)

    tr = Dict("mode" => "lines", "line" => Dict("width" => 0), "fill" => "tozeroy")

    ie_ = map(<(0.0), sc_)

    ip_ = map(>=(0.0), sc_)

    ny = "Δ Enrichment"

    da_ = [
        merge(
            tr,
            Dict(
                "name" => "- Score",
                "y" => sc_[ie_],
                "x" => xc_[ie_],
                "text" => fe_[ie_],
                "fillcolor" => Omics.Color.BL,
            ),
        ),
        merge(
            tr,
            Dict(
                "name" => "+ Score",
                "y" => sc_[ip_],
                "x" => xc_[ip_],
                "text" => fe_[ip_],
                "fillcolor" => Omics.Color.RE,
            ),
        ),
        Dict(
            "yaxis" => "y2",
            "name" => "Set",
            "y" => zeros(sum(ii_)),
            "x" => xc_[ii_],
            "text" => fe_[ii_],
            "mode" => "markers",
            "marker" => Dict(
                "symbol" => "line-ns",
                "size" => 24,
                "line" => Dict(
                    "width" => 2,
                    "color" => Omics.Color.hexify(Omics.Color.SG, 0.72),
                ),
            ),
            "hoverinfo" => "x+text",
        ),
        merge(
            tr,
            Dict(
                "yaxis" => "y3",
                "name" => ny,
                "y" => mo_,
                "x" => xc_,
                "text" => fe_,
                "fillcolor" => "#07fa07",
            ),
        ),
    ]

    if typeof(al) == KS

        ix_ = map(in(_get_extreme(mo_)), mo_)

        push!(
            da_,
            Dict(
                "yaxis" => "y3",
                "name" => "Extrema",
                "y" => mo_[ix_],
                "x" => xc_[ix_],
                "text" => fe_[ix_],
                "mode" => "markers",
                "marker" =>
                    Dict("size" => 32, "color" => Omics.Color.hexify(Omics.Color.HU, 0.72)),
            ),
        )

    end

    an = Dict(
        "y" => 0,
        "font" => Dict("size" => 16),
        "borderpad" => 4.8,
        "borderwidth" => 2.64,
        "bordercolor" => Omics.Color.LI,
        "showarrow" => false,
    )

    ax = uf * 0.008

    if haskey(la, "title") && haskey(la["title"], "text")

        la["title"]["text"] = Omics.Strin.limit(la["title"]["text"], 56)

    end

    Omics.Plot.plot(
        ht,
        da_,
        Omics.Dic.merg(
            Dict(
                "showlegend" => false,
                "yaxis" => Dict("domain" => (0, 0.24), "title" => Dict("text" => ns)),
                "yaxis2" => Dict(
                    "domain" => (0.248, 0.32),
                    "title" => Dict("text" => "Set"),
                    "tickvals" => (),
                ),
                "yaxis3" => Dict("domain" => (0.328, 1), "title" => Dict("text" => ny)),
                "xaxis" => Dict(
                    "title" => Dict("text" => "$nf ($uf)"),
                    "showspikes" => true,
                    "spikemode" => "across",
                    "spikedash" => "solid",
                    "spikethickness" => 0.8,
                    "spikecolor" => "#000000",
                ),
                "annotations" => (
                    Dict(
                        "yref" => "paper",
                        "xref" => "paper",
                        "y" => 1.064,
                        "text" => "Enrichment = <b>$(@sprintf "%.4g" en)</b>",
                        "font" => Dict("size" => 24, "color" => Omics.Color.BR),
                        "showarrow" => false,
                    ),
                    merge(
                        an,
                        Dict(
                            "x" => 1.0 - ax,
                            "xanchor" => "right",
                            "text" => nh,
                            "font" => Dict("color" => Omics.Color.RE),
                        ),
                    ),
                    merge(
                        an,
                        Dict(
                            "x" => uf + ax,
                            "xanchor" => "left",
                            "text" => nl,
                            "font" => Dict("color" => Omics.Color.BL),
                        ),
                    ),
                ),
            ),
            la,
        ),
    )

end

function _select_sort(fe_, sc_)

    id_ = findall(!isnan, sc_)

    fe_ = fe_[id_]

    sc_ = sc_[id_]

    sortperm!(id_, sc_; rev = true)

    fe_[id_], sc_[id_]

end

function enrich(al, fe_, sc_, me___; ex = 1.0, mi = 1, ma = 1000, fr = 0.0)

    fe_, sc_ = _select_sort(fe_, sc_)

    en_ = Vector{Float64}(undef, lastindex(me___))

    ii_ = falses(lastindex(fe_))

    fe_ie = Dict(fe_[ie] => ie for ie in eachindex(fe_))

    for is in eachindex(me___)

        me_ = me___[is]

        _is_in!(ii_, fe_ie, me_)

        ui = sum(ii_)

        en_[is] =
            ui < mi || ma < ui || ui / lastindex(me_) < fr ? NaN :
            _enrich!(al, sc_, ex, ii_, nothing)

        # TODO: Avoid broadcasting.
        ii_[ii_] .= false

    end

    en_

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

    end

end

function _separat(se_me_)

    collect(keys(se_me_)), collect(values(se_me_))

end

function data_rank!(di, al, fe_, sc, ne, se_me_, ns, sa_; st = 0.0, up = 2, ke_ar...)

    if !iszero(st)

        foreach(sc_ -> Omics.XSample.standardize_clamp!(sc_, st), eachcol(sc))

    end

    se_, me___ = _separat(se_me_)

    en = stack((enrich(al, fe_, sc_, me___; ke_ar...) for sc_ in eachcol(sc)))

    ig_ = map(en_ -> all(!isnan, en_), eachrow(en))

    se_ = se_[ig_]

    me___ = me___[ig_]

    en = en[ig_, :]

    pr = joinpath(di, "enrichment")

    Omics.XSample.write_plot(pr, ne, se_, ns, sa_, strin(al), en)

    ig_ = map(!isnan, en)

    for id_ in CartesianIndices(en)[ig_][Omics.Extreme.ge(en[ig_], up)]

        is, ia = Tuple(id_)

        se = se_[is]

        sa = sa_[ia]

        plot(
            joinpath(di, "$(Omics.Numbe.shorten(en[is, ia])).$sa.$se.html"),
            al,
            _select_sort(fe_, sc[:, ia])...,
            me___[is];
            ns = sa,
            la = Dict("title" => Dict("text" => se)),
        )

    end

    se_, en

end

"""
Run data-rank (single-sample) GSEA.

# Arguments

  - `output_directory`:
  - `feature_x_sample_x_score_tsv`:
  - `set_features_json`:

# Options

  - `--standard-deviation`: = 0.0. For normalization by column. 0.0 skips normalization.
  - `--algorithm`: = "ks". "ks" | "ksa" | "kliom" | "kliop" | "kli" | "kli1".
  - `--exponent`: = 1.0.
  - `--minimum-set-size`: = 1.
  - `--maximum-set-size`: = 1000.
  - `--set-fraction`: = 0.0.
  - `--number-of-sets-to-plot`: = 2.
  - `--set-name`: = "Set".
  - `--sample-name`: = "Sample".
"""
@cast function data_rank(
    output_directory,
    feature_x_sample_x_score_tsv,
    set_features_json;
    standard_deviation::Float64 = 0.0,
    algorithm = "ks",
    exponent::Float64 = 1.0,
    minimum_set_size::Int = 1,
    maximum_set_size::Int = 1000,
    set_fraction::Float64 = 0.0,
    number_of_sets_to_plot::Int = 2,
    set_name = "Set",
    sample_name = "Sample",
)

    ta = Omics.Table.rea(feature_x_sample_x_score_tsv)

    data_rank!(
        output_directory,
        _set_algorithm(algorithm),
        ta[!, 1],
        Matrix(ta[!, 2:end]),
        set_name,
        Omics.Dic.rea(set_features_json),
        sample_name,
        names(ta)[2:end];
        st = standard_deviation,
        up = number_of_sets_to_plot,
        ex = exponent,
        mi = minimum_set_size,
        ma = maximum_set_size,
        fr = set_fraction,
    )

end

function _normalize_enrichment(::Union{KS, KSa}, en, mn, mp, ::Any, ::Any)

    en / (en < 0.0 ? -mn : mp)

end

function _normalize_enrichment(::Any, en, mn, mp, sn, sp)

    if en < 0.0

        me = mn

        st = -sn

    else

        me = mp

        st = sp

    end

    # TODO: Check `3 * st` or `(3 * st)`.
    1.0 + (en - me) / 3.0 * st

end

function _normalize_enrichment!(al, en_, ra)

    us = lastindex(en_)

    no_ = Vector{Float64}(undef, us)

    for id in 1:us

        en = en_[id]

        ra_ = ra[id, :]

        rn_, rp_ = Omics.Significance._separate(ra_)

        mn = mean(rn_)

        mp = mean(rp_)

        sn = std(rn_)

        sp = std(rp_)

        no_[id] = _normalize_enrichment(al, en, mn, mp, sn, sp)

        ra[id, :] = map(ra -> _normalize_enrichment(al, ra, mn, mp, sn, sp), ra_)

    end

    no_

end

function _write_plot(di, al, fe_, sc_, ex, se_, me___, en_, ra, up, pl_, nf, ns, nl, nh)

    ig_ = map(!isnan, en_)

    se_ = se_[ig_]

    me___ = me___[ig_]

    en_ = en_[ig_]

    ra = ra[ig_, :]

    no_ = _normalize_enrichment!(al, en_, ra)

    pn_, qn_, pp_, qp_ = Omics.Significance.ge(ra, no_)

    Omics.Table.writ(
        joinpath(di, "result.tsv"),
        Omics.Table.make(
            "Set",
            se_,
            ["Enrichment", "Normalized Enrichment", "P-Value", "Q-Value"],
            stack((en_, no_, vcat(pn_, pp_), vcat(qn_, qp_))),
        ),
    )

    fe_, sc_ = _select_sort(fe_, sc_)

    for is in
        unique!(vcat(Omics.Extreme.ge(en_, up), filter!(!isnothing, indexin(pl_, se_))))

        se = se_[is]

        plot(
            joinpath(di, "$(Omics.Numbe.shorten(en_[is])).$se.html"),
            al,
            fe_,
            sc_,
            me___[is];
            ex,
            nf,
            ns,
            nl,
            nh,
            la = Dict("title" => Dict("text" => se)),
        )

    end

end

function _permute_set(ur, sd, al, fe_, sc_, me___; ke_ar...)

    ra = Matrix{Float64}(undef, lastindex(me___), ur)

    if !iszero(ur)

        um_ = map(lastindex, me___)

        seed!(sd)

        @showprogress for id in 1:ur

            ra[:, id] = enrich(
                al,
                fe_,
                sc_,
                map(um -> sample(fe_, um; replace = false), um_);
                ke_ar...,
            )

        end

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

  - `--algorithm`: = "ks". "ks" | "ksa" | "kliom" | "kliop" | "kli" | "kli1".
  - `--exponent`: = 1.0.
  - `--minimum-set-size`: = 1.
  - `--maximum-set-size`: = 1000.
  - `--set-fraction`: = 0.0.
  - `--number-of-permutations`: = 100.
  - `--random-seed`: = 20150603.
  - `--number-of-sets-to-plot`: = 2.
  - `--more-sets-to-plot`: = "". ;-separated set names.
  - `--feature-name`: = "Gene".
  - `--score-name`: = "My Score".
  - `--low-text`: = "Low".
  - `--high-text`: = "High".
"""
@cast function user_rank(
    output_directory,
    feature_x_metric_x_score_tsv,
    set_features_json;
    algorithm = "ks",
    exponent::Float64 = 1.0,
    minimum_set_size::Int = 1,
    maximum_set_size::Int = 1000,
    set_fraction::Float64 = 0.0,
    number_of_permutations::Int = 100,
    random_seed::Int = 20150603,
    number_of_sets_to_plot::Int = 2,
    more_sets_to_plot = "",
    feature_name = "Gene",
    score_name = "My Score",
    low_text = "Low",
    high_text = "High",
)

    al = _set_algorithm(algorithm)

    ta = Omics.Table.rea(feature_x_metric_x_score_tsv)

    fe_, sc_ = _select_sort(ta[!, 1], ta[!, 2])

    se_, me___ = _separat(Omics.Dic.rea(set_features_json))

    ke_ar = (ex = exponent, mi = minimum_set_size, ma = maximum_set_size, fr = set_fraction)

    _write_plot(
        output_directory,
        al,
        fe_,
        sc_,
        exponent,
        se_,
        me___,
        enrich(al, fe_, sc_, me___; ke_ar...),
        _permute_set(number_of_permutations, random_seed, al, fe_, sc_, me___; ke_ar...),
        number_of_sets_to_plot,
        split(more_sets_to_plot, ';'),
        feature_name,
        score_name,
        low_text,
        high_text,
    )

end

"""
Run metric-rank (standard) GSEA.

# Arguments

  - `output_directory`:
  - `target_x_sample_x_number_tsv`:
  - `feature_x_sample_x_score_tsv`:
  - `set_features_json`:

# Options

  - `--standard-deviation`: = 0.0. For normalization by column. 0.0 skips normalization.
  - `--algorithm`: = "ks". "ks" | "ksa" | "kliom" | "kliop" | "kli" | "kli1".
  - `--exponent`: = 1.0.
  - `--metric`: = "signal-to-noise-ratio". "mean-difference" | "log-ratio" | "signal-to-noise-ratio".
  - `--minimum-set-size`: = 1.
  - `--maximum-set-size`: = 1000.
  - `--set-fraction`: = 0.0.
  - `--permutation`: = "sample". "sample" | "set".
  - `--number-of-permutations`: = 100.
  - `--random-seed`: = 20150603.
  - `--number-of-sets-to-plot`: = 2.
  - `--more-sets-to-plot`: = "". ;-separated set names.
  - `--feature-name`: = "Gene".
  - `--score-name`: = "Signal-to-Noise Ratio".
  - `--low-text`: = "Low".
  - `--high-text`: = "High".
"""
@cast function metric_rank(
    output_directory,
    target_x_sample_x_number_tsv,
    feature_x_sample_x_score_tsv,
    set_features_json;
    standard_deviation::Float64 = 0.0,
    algorithm = "ks",
    exponent::Float64 = 1.0,
    metric = "signal-to-noise-ratio",
    minimum_set_size::Int = 1,
    maximum_set_size::Int = 1000,
    set_fraction::Float64 = 0.0,
    permutation = "sample",
    number_of_permutations::Int = 100,
    random_seed::Int = 20150603,
    number_of_sets_to_plot::Int = 2,
    more_sets_to_plot = "",
    feature_name = "Gene",
    score_name = "Signal-to-Noise Ratio",
    low_text = "Low",
    high_text = "High",
)

    tt = Omics.Table.rea(target_x_sample_x_number_tsv)

    tf = Omics.Table.rea(feature_x_sample_x_score_tsv)

    vt_ = convert(BitVector, collect(tt[1, 2:end]))

    fe_ = tf[!, 1]

    s1 = Matrix(tf[!, indexin(names(tt)[2:end], names(tf))])

    if !iszero(standard_deviation)

        foreach(
            s1_ -> Omics.XSample.standardize_clamp!(s1_, standard_deviation),
            eachcol(s1),
        )

    end

    fu = if metric == "mean-difference"

        Omics.Target.get_mean_difference

    elseif metric == "log-ratio"

        Omics.Target.get_log_ratio

    elseif metric == "signal-to-noise-ratio"

        Omics.Target.get_signal_to_noise_ratio

    end

    s2_ = map(s1_ -> Omics.Target.go(fu, vt_, s1_), eachrow(s1))

    Omics.Table.writ(
        joinpath(output_directory, "metric.tsv"),
        Omics.Table.make("Feature", fe_, [metric], reshape(s2_, :, 1)),
    )

    al = _set_algorithm(algorithm)

    se_, me___ = _separat(Omics.Dic.rea(set_features_json))

    ke_ar = (ex = exponent, mi = minimum_set_size, ma = maximum_set_size, fr = set_fraction)

    if permutation == "set"

        ra =
            _permute_set(number_of_permutations, random_seed, al, fe_, s2_, me___; ke_ar...)

    elseif permutation == "sample"

        ra = Matrix{Float64}(undef, lastindex(se_), number_of_permutations)

        if 0 < number_of_permutations

            seed!(random_seed)

            @showprogress for id in 1:number_of_permutations

                ra[:, id] = enrich(
                    al,
                    fe_,
                    map(s1_ -> Omics.Target.go(fu, shuffle!(vt_), s1_), eachrow(s1)),
                    me___;
                    ke_ar...,
                )

            end

        end

    end

    _write_plot(
        output_directory,
        al,
        fe_,
        s2_,
        exponent,
        se_,
        me___,
        enrich(al, fe_, s2_, me___; ke_ar...),
        ra,
        number_of_sets_to_plot,
        split(more_sets_to_plot, ';'),
        feature_name,
        score_name,
        low_text,
        high_text,
    )

end

"""
"""
@main

end
