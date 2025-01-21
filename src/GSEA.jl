module GSEA

# ----------------------------------------------------------------------------------------------- #

using Comonicon: @cast, @main

using Printf: @sprintf

using ProgressMeter: @showprogress

using Random: seed!, shuffle!

using StatsBase: mean, sample, std

using Omics

# TODO: Generalize.

function _is_in(fe_, f1_)

    map(in(Set(f1_)), fe_)

end

function _is_in!(is_, fe_id, f1_)

    for f1 in f1_

        id = get(fe_id, f1, nothing)

        if !isnothing(id)

            is_[id] = true

        end

    end

end

function _separate(nu_)

    ie_ = findall(<(0), nu_)

    ip_ = findall(>=(0), nu_)

    nu_[ie_], nu_[ip_]

end

function read_cls(cl)

    l1, l2, l3 = readlines(cl)

    ph = l2[2:end]

    va_ = split(l3)

    us = lastindex(va_)

    sa_ = Omics.Simulation.label(us, "Sample")

    if l1 == "#numeric"

        Omics.Table.make("Target", ph, sa_, [parse(Float64, va) for _ in 1:1, va in va_])

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

        Omics.Table.make("Target", join(ph_, '_'), sa_, [ph_id[va] for _ in 1:1, va in va_])

    end

end

function read_gct(gc)

    Omics.Table.rea(gc; header = 3, drop = ["Description"])

end

function read_gmt(gm)

    se_ge_ = Dict{String, Vector{String}}()

    for li in eachline(gm)

        sp_ = split(li, '\t')

        se = sp_[1]

        if haskey(se_ge_, se)

            error("there is more than one $se.")

        end

        se_ge_[se] = filter!(!isempty, sp_[3:lastindex(sp_)])

    end

    se_ge_

end

"""
Convert .cls to .tsv.

# Arguments

  - `tsv`:
  - `cls`:
"""
@cast function cls(tsv, cls)

    ta = read_cls(cls)

    Omics.Table.writ(
        tsv,
        Omics.Table.make(
            "Target",
            ta[:, 1],
            names(ta)[2:end],
            map(nu -> nu - 1.0, Matrix(ta[:, 2:end])),
        ),
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
        Omics.Table.make("Feature", ta[:, 1], names(ta)[2:end], Matrix(ta[:, 2:end])),
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

function select_set(se_me_, fe_, mi, ma, fr)

    me___ = collect(values(se_me_))

    ke_ = BitVector(undef, lastindex(me___))

    for id in eachindex(me___)

        me_ = me___[id]

        ub = lastindex(me_)

        intersect!(me_, fe_)

        ua = lastindex(me_)

        ke_[id] = mi <= ua <= ma && fr <= ua / ub

    end

    collect(keys(se_me_))[ke_], me___[ke_]

end

struct KS end

struct KSa end

struct KLioM end

struct KLioP end

struct KLi end

struct KLi1 end

function _exponentiate(nu, ex)

    ab = abs(nu)

    isone(ex) ? ab : ab^ex

end

function _get_delta(::Union{KS, KSa}, sc_, ex, is_)

    s0 = s1 = 0.0

    for id in eachindex(sc_)

        if is_[id]

            s1 += _exponentiate(sc_[id], ex)

        else

            s0 += 1.0

        end

    end

    -inv(s0), inv(s1)

end

function _get_delta(::Union{KLioM, KLioP, KLi}, sc_, ex, is_)

    sa = s1 = 0.0

    for id in eachindex(sc_)

        ma = _exponentiate(sc_[id], ex)

        sa += ma

        if is_[id]

            s1 += ma

        end

    end

    inv(sa), inv(s1)

end

function _get_delta(na, n1)

    inv(inv(na) - inv(n1))

end

function _get_delta(::KLi1, sc_, ex, is_)

    s1 = 0.0

    for id in eachindex(sc_)

        if is_[id]

            s1 += _exponentiate(sc_[id], ex)

        end

    end

    inv(s1)

end

function _enrich!(al::KS, sc_, ex, is_, mo_)

    n0, n1 = _get_delta(al, sc_, ex, is_)

    mo = be = bs = 0.0

    for id in eachindex(sc_)

        mo += is_[id] ? _exponentiate(sc_[id], ex) * n1 : n0

        if !isnothing(mo_)

            mo_[id] = mo

        end

        cu = abs(mo)

        if be < cu

            be = cu

            bs = mo

        end

    end

    bs

end

function _enrich!(al::KSa, sc_, ex, is_, mo_)

    n0, n1 = _get_delta(al, sc_, ex, is_)

    mo = ar = 0.0

    for id in eachindex(sc_)

        ar += mo += is_[id] ? _exponentiate(sc_[id], ex) * n1 : n0

        if !isnothing(mo_)

            mo_[id] = mo

        end

    end

    ar / lastindex(sc_)

end

const ON = 1.0 + 1e-13

function _enrich!(al::KLioM, sc_, ex, is_, mo_)

    na, n1 = _get_delta(al, sc_, ex, is_)

    n0 = _get_delta(na, n1)

    ra = r0 = r1 = eps()

    pa = p0 = p1 = ar = 0.0

    la = l0 = l1 = ON

    for id in eachindex(sc_)

        ma = _exponentiate(sc_[id], ex)

        da = ma * na

        if is_[id]

            d0 = 0.0

            d1 = ma * n1

        else

            d0 = ma * n0

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

function _enrich!(al::KLioP, sc_, ex, is_, mo_)

    na, n1 = _get_delta(al, sc_, ex, is_)

    n0 = _get_delta(na, n1)

    ra = r0 = r1 = eps()

    pa = p0 = p1 = ar = 0.0

    la = l0 = l1 = ON

    for id in eachindex(sc_)

        ma = _exponentiate(sc_[id], ex)

        da = ma * na

        if is_[id]

            d0 = 0.0

            d1 = ma * n1

        else

            d0 = ma * n0

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

function _enrich!(al::KLi, sc_, ex, is_, mo_)

    na, n1 = _get_delta(al, sc_, ex, is_)

    ra = r1 = eps()

    pa = p1 = ar = 0.0

    la = l1 = ON

    for id in eachindex(sc_)

        ma = _exponentiate(sc_[id], ex)

        da = ma * na

        d1 = is_[id] ? ma * n1 : 0.0

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

function _enrich!(al::KLi1, sc_, ex, is_, mo_)

    uf = lastindex(sc_)

    n1 = _get_delta(al, sc_, ex, is_)

    ra = r1 = eps()

    da = inv(uf)

    p1 = ar = 0.0

    la = ON + da

    l1 = ON

    for id in eachindex(sc_)

        d1 = is_[id] ? _exponentiate(sc_[id], ex) * n1 : 0.0

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

    if ai == aa

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
    nl = "Low",
    nh = "High",
    ns = "Score",
    la = Dict{String, Any}(),
)

    uf = lastindex(fe_)

    xc_ = collect(1:uf)

    is_ = _is_in(fe_, me_)

    mo_ = Vector{Float64}(undef, uf)

    en = _enrich!(al, sc_, ex, is_, mo_)

    tr = Dict("mode" => "lines", "line" => Dict("width" => 0), "fill" => "tozeroy")

    ie_ = findall(<(0.0), sc_)

    ip_ = findall(>=(0.0), sc_)

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
            "y" => zeros(sum(is_)),
            "x" => xc_[is_],
            "text" => fe_[is_],
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

        id_ = findall(in(_get_extreme(mo_)), mo_)

        push!(
            da_,
            Dict(
                "yaxis" => "y3",
                "name" => "Extrema",
                "y" => mo_[id_],
                "x" => xc_[id_],
                "text" => fe_[id_],
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

    ex = uf * 0.008

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
                    "title" => Dict("text" => Omics.Strin.coun(uf, nf)),
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
                            "x" => 1.0 - ex,
                            "xanchor" => "right",
                            "text" => nh,
                            "font" => Dict("color" => Omics.Color.RE),
                        ),
                    ),
                    merge(
                        an,
                        Dict(
                            "x" => uf + ex,
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

function enrich(al, fe_, sc_::AbstractVector, me___; um = 1, ex = 1)

    en_ = Vector{Float64}(undef, lastindex(me___))

    is_ = falses(lastindex(fe_))

    fe_id = Dict(fe_[id] => id for id in eachindex(fe_))

    for id in eachindex(me___)

        _is_in!(is_, fe_id, me___[id])

        en_[id] = sum(is_) < um ? NaN : _enrich!(al, sc_, ex, is_, nothing)

        is_[is_] .= false

    end

    en_

end

function enrich(al, fe_, sc, me___; um = 1, ex = 1)

    us = size(sc, 2)

    en = Matrix{Float64}(undef, lastindex(me___), us)

    is_ = BitVector(undef, lastindex(fe_))

    for id in 1:us

        sc_ = sc[:, id]

        map!(!isnan, is_, sc_)

        so_ = sc_[is_]

        id_ = sortperm(so_; rev = true)

        en[:, id] = enrich(al, fe_[is_][id_], so_[id_], me___; um, ex)

    end

    en

end

function write_plot(ou, fe_, sc, al, se_, me___, ns, sa_, en; ex = 1.0, up = 4)

    pr = joinpath(ou, "enrichment")

    Omics.XSample.write_plot(
        pr,
        "Set",
        se_,
        ns,
        sa_,
        "Enrichment using $(string(al)[6:(end - 2)])",
        en,
    )

    ge_ = map(!isnan, en)

    gs_ = BitVector(undef, lastindex(fe_))

    for id_ in CartesianIndices(en)[ge_][Omics.Extreme.ge(en[ge_], up)]

        ie, ia = Tuple(id_)

        sc_ = sc[:, ia]

        map!(!isnan, gs_, sc_)

        so_ = sc_[gs_]

        id_ = sortperm(so_; rev = true)

        ti = "$(sa_[ia]) Enriching $(se_[ie])"

        plot(
            joinpath(ou, "$(Omics.Strin.shorten(en[ie, ia])).$ti.html"),
            al,
            fe_[gs_][id_],
            so_[id_],
            me___[ie];
            ex,
            la = Dict("title" => Dict("text" => ti)),
        )

    end

end

function _standardize_clamp!(ma, di, st)

    if isone(di)

        ea = eachrow

    elseif di == 2

        ea = eachcol

    end

    foreach(Omics.Normalization.normalize_with_0!, ea(ma))

    clamp!(ma, -st, st)

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

"""
Run data-rank (single-sample) GSEA.

# Arguments

  - `output_directory`:
  - `feature_x_sample_x_score_tsv`:
  - `set_features_json`:

# Options

  - `--normalization-dimension`: = 0. 0 (not normalizing) | 1 | 2.
  - `--normalization-standard-deviation`: = 4.0.
  - `--algorithm`: = "ks". "ks" | "ksa" | "kliom" | "kliop" | "kli" | "kli1".
  - `--exponent`: = 1.0.
  - `--post-skip-minimum-set-size`: = 1.
  - `--minimum-set-size`: = 15.
  - `--maximum-set-size`: = 500.
  - `--minimum-set-fraction`: = 0.0.

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
    algorithm = "ks",
    exponent::Float64 = 1.0,
    post_skip_minimum_set_size::Int = 1,
    minimum_set_size::Int = 15,
    maximum_set_size::Int = 500,
    minimum_set_fraction::Float64 = 0.0,
)

    ta = Omics.Table.rea(feature_x_sample_x_score_tsv)

    fe_ = ta[!, 1]

    sa_ = names(ta)[2:end]

    sc = Matrix(ta[!, 2:end])

    if skip_0

        replace!(sc, 0.0 => NaN)

    end

    if !iszero(normalization_dimension)

        _standardize_clamp!(sc, normalization_dimension, normalization_standard_deviation)

    end

    al = _set_algorithm(algorithm)

    se_, me___ = select_set(
        Omics.Dic.rea(set_features_json),
        fe_,
        minimum_set_size,
        maximum_set_size,
        minimum_set_fraction,
    )

    write_plot(
        output_directory,
        fe_,
        sc,
        al,
        se_,
        me___,
        "Sample",
        sa_,
        enrich(al, fe_, sc, me___; um = post_skip_minimum_set_size, ex = exponent);
        ex = exponent,
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

    1.0 + (en - me) / 3.0 * st

end

function _normalize_enrichment!(al, en_, ra)

    er_ = Vector{Float64}(undef, lastindex(en_))

    for id in eachindex(en_)

        en = en_[id]

        ra_ = ra[id, :]

        ne_, po_ = _separate(ra_)

        mn = mean(ne_)

        mp = mean(po_)

        sn = std(ne_)

        sp = std(po_)

        er_[id] = _normalize_enrichment(al, en, mn, mp, sn, sp)

        ra[id, :] = map(rn -> _normalize_enrichment(al, rn, mn, mp, sn, sp), ra_)

    end

    er_

end

function _write_plot(ou, wr, al, fe_, sc_, ex, se_, me___, en_, ra, up, pl_, nf, nl, nh, ns)

    ue = lastindex(se_)

    re = fill(NaN, ue, 4)

    id_ = sortperm(en_)

    se_ = se_[id_]

    me___ = me___[id_]

    en_ = en_[id_]

    ra = ra[id_, :]

    if wr

        Omics.Table.writ(
            joinpath(ou, "random.tsv"),
            Omics.Table.make("Set", se_, axes(ra, 2), ra),
        )

    end

    er_ = _normalize_enrichment!(al, en_, ra)

    id = findlast(<(0.0), en_)

    ie_ = 1:id

    ip_ = (id + 1):ue

    np_, nq_, pp_, pq_ = Omics.Significance.ge(ra, er_, ie_, ip_)

    re = stack((en_, er_, vcat(np_, pp_), vcat(nq_, pq_)))

    Omics.Table.writ(
        joinpath(ou, "result.tsv"),
        Omics.Table.make(
            "Set",
            se_,
            ["Enrichment", "Normalized Enrichment", "P-Value", "Q-Value"],
            re,
        ),
    )

    for id in unique!(vcat(Omics.Extreme.ge(er_, up), indexin(pl_, se_)))

        if isnothing(id)

            continue

        end

        ti = "$id $(se_[id])"

        plot(
            joinpath(ou, "$ti.html"),
            al,
            fe_,
            sc_,
            me___[id];
            ex,
            nf,
            nl,
            nh,
            ns,
            la = Dict("title" => Dict("text" => ti)),
        )

    end

end

function _permute_set(ur, se, al, fe_, sc_, me___, ex)

    ra = Matrix{Float64}(undef, lastindex(me___), ur)

    if !iszero(ur)

        um_ = map(lastindex, me___)

        seed!(se)

        @showprogress for id in 1:ur

            ra[:, id] =
                enrich(al, fe_, sc_, map(um -> sample(fe_, um; replace = false), um_); ex)

        end

    end

    ra

end

# TODO: Retire
function _use_permutation(permutation, al, fe_, me___, ex)

    ta = Omics.Table.rea(permutation)

    sc = Matrix(ta[indexin(fe_, ta[!, 1]), 2:end])

    ra = Matrix{Float64}(undef, lastindex(me___), size(ra, 2))

    @showprogress for id in axes(sc, 2)

        sc_ = sc[:, id]

        id_ = sortperm(sc_; rev = true)

        ra[:, id] = enrich(al, fe_[id_], sc_[id_], me___; ex)

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
  - `--minimum-set-size`: = 15.
  - `--maximum-set-size`: = 500.
  - `--minimum-set-fraction`: = 0.0.
  - `--permutation`: = "set". "set" | "feature_x_index_x_random.tsv".
  - `--number-of-permutations`: = 100.
  - `--random-seed`: = 20150603.
  - `--number-of-sets-to-plot`: = 4.
  - `--more-sets-to-plot`: = "". Space-separated set names.
  - `--feature-name`: = "Gene".
  - `--low-text`: = "Low".
  - `--high-text`: = "High".
  - `--score-name`: = "User-Defined Score".

# Flags

  - `--write-set-x-index-x-random-tsv`: = false.
"""
@cast function user_rank(
    output_directory,
    feature_x_metric_x_score_tsv,
    set_features_json;
    algorithm = "ks",
    exponent::Float64 = 1.0,
    minimum_set_size::Int = 15,
    maximum_set_size::Int = 500,
    minimum_set_fraction::Float64 = 0.0,
    permutation = "set",
    number_of_permutations::Int = 100,
    random_seed::Int = 20150603,
    write_set_x_index_x_random_tsv::Bool = false,
    number_of_sets_to_plot::Int = 4,
    more_sets_to_plot = "",
    feature_name = "Gene",
    low_text = "Low",
    high_text = "High",
    score_name = "User-Defined Score",
)

    al = _set_algorithm(algorithm)

    ta = Omics.Table.rea(feature_x_metric_x_score_tsv)

    fe_ = ta[!, 1]

    sc_ = ta[!, 2]

    id_ = sortperm(sc_; rev = true)

    fe_ = fe_[id_]

    sc_ = sc_[id_]

    se_, me___ = select_set(
        Omics.Dic.rea(set_features_json),
        fe_,
        minimum_set_size,
        maximum_set_size,
        minimum_set_fraction,
    )

    ra = if permutation == "set"

        _permute_set(number_of_permutations, random_seed, al, fe_, sc_, me___, exponent)

    elseif isfile(permutation)

        _use_permutation(permutation, al, fe_, me___, exponent)

    end

    _write_plot(
        output_directory,
        write_set_x_index_x_random_tsv,
        al,
        fe_,
        sc_,
        exponent,
        se_,
        me___,
        enrich(al, fe_, sc_, me___; ex = exponent),
        ra,
        number_of_sets_to_plot,
        split(more_sets_to_plot),
        feature_name,
        low_text,
        high_text,
        score_name,
    )

end

function _target_sort(fu, is_, sc, fe_)

    re_ = map(sc_ -> Omics.Target.go(fu, is_, sc_), eachrow(sc))

    id_ = sortperm(re_; rev = true)

    fe_[id_], re_[id_]

end

"""
Run metric-rank (standard) GSEA.

# Arguments

  - `output_directory`:
  - `target_x_sample_x_number_tsv`:
  - `feature_x_sample_x_score_tsv`:
  - `set_features_json`:

# Options

  - `--normalization-dimension`: = 0. 0 (not normalizing) | 1 | 2.
  - `--normalization-standard-deviation`: = 4.0.
  - `--algorithm`: = "ks". "ks" | "ksa" | "kliom" | "kliop" | "kli" | "kli1".
  - `--exponent`: = 1.0.
  - `--metric`: = "signal-to-noise-ratio". "mean-difference" | "log-ratio" | "signal-to-noise-ratio".
  - `--minimum-set-size`: = 15.
  - `--maximum-set-size`: = 500.
  - `--minimum-set-fraction`: = 0.0.
  - `--permutation`: = "sample". "sample" | "set" | "feature_x_index_x_random.tsv".
  - `--number-of-permutations`: = 100.
  - `--random-seed`: = 20150603.
  - `--number-of-sets-to-plot`: = 4.
  - `--more-sets-to-plot`: = "". Space-separated set names.
  - `--feature-name`: = "Gene".
  - `--low-text`: = "Low".
  - `--high-text`: = "High".
  - `--score-name`: = "Signal-to-Noise Ratio".

# Flags

  - `--write-set-x-index-x-random-tsv`: = false.
"""
@cast function metric_rank(
    output_directory,
    target_x_sample_x_number_tsv,
    feature_x_sample_x_score_tsv,
    set_features_json;
    normalization_dimension::Int = 0,
    normalization_standard_deviation::Float64 = 4.0,
    algorithm = "ks",
    exponent::Float64 = 1.0,
    metric = "signal-to-noise-ratio",
    minimum_set_size::Int = 15,
    maximum_set_size::Int = 500,
    minimum_set_fraction::Float64 = 0.0,
    permutation = "sample",
    number_of_permutations::Int = 100,
    random_seed::Int = 20150603,
    write_set_x_index_x_random_tsv::Bool = false,
    number_of_sets_to_plot::Int = 4,
    more_sets_to_plot = "",
    feature_name = "Gene",
    low_text = "Low",
    high_text = "High",
    score_name = "Signal-to-Noise Ratio",
)

    tt = Omics.Table.rea(target_x_sample_x_number_tsv)

    tf = Omics.Table.rea(feature_x_sample_x_score_tsv)

    sc = Matrix(tf[!, indexin(names(tt)[2:end], names(tf))])

    if !iszero(normalization_dimension)

        _standardize_clamp!(sc, normalization_dimension, normalization_standard_deviation)

    end

    fu = if metric == "mean-difference"

        Omics.Target.get_mean_difference

    elseif metric == "log-ratio"

        Omics.Target.get_log_ratio

    elseif metric == "signal-to-noise-ratio"

        Omics.Target.get_signal_to_noise_ratio

    end

    is_ = convert(BitVector, collect(tt[1, 2:end]))

    fe_, mt_ = _target_sort(fu, is_, sc, tf[!, 1])

    Omics.Table.writ(
        joinpath(output_directory, "feature_x_metric_x_score.tsv"),
        Omics.Table.make("Feature", fe_, [metric], reshape(mt_, :, 1)),
    )

    se_, me___ = select_set(
        Omics.Dic.rea(set_features_json),
        fe_,
        minimum_set_size,
        maximum_set_size,
        minimum_set_fraction,
    )

    al = _set_algorithm(algorithm)

    if permutation == "sample"

        ra = Matrix{Float64}(undef, lastindex(se_), number_of_permutations)

        if 0 < number_of_permutations

            seed!(random_seed)

            @showprogress for id in 1:number_of_permutations

                ra[:, id] = enrich(
                    al,
                    _target_sort(fu, shuffle!(is_), sc, fe_)...,
                    me___;
                    ex = exponent,
                )

            end

        end

    elseif permutation == "set"

        ra =
            _permute_set(number_of_permutations, random_seed, al, fe_, mt_, me___, exponent)

    elseif isfile(permutation)

        ra = _use_permutation(permutation, al, fe_, me___, exponent)

    end

    _write_plot(
        output_directory,
        write_set_x_index_x_random_tsv,
        al,
        fe_,
        mt_,
        exponent,
        se_,
        me___,
        enrich(al, fe_, mt_, me___; ex = exponent),
        ra,
        number_of_sets_to_plot,
        split(more_sets_to_plot),
        feature_name,
        low_text,
        high_text,
        score_name,
    )

end

"""
Gene set enrichment analysis.
"""
@main

end
