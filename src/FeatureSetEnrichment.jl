module FeatureSetEnrichment

using ProgressMeter: @showprogress

using ..BioLab

struct KS end

struct KSa end

struct KLi end

struct KLioM end

struct KLioP end

function _make_string(al)

    chop(string(al); head = 28, tail = 2)

end

@inline function _index_absolute_exponentiate(sc_, id, ex)

    ab = abs(sc_[id])

    if !isone(ex)

        ab ^= ex

    end

    ab

end

@inline function _sum_01(sc_, ex, is_)

    n = length(sc_)

    su0 = su1 = 0.0

    for id in 1:n

        if is_[id]

            su1 += _index_absolute_exponentiate(sc_, id, ex)

        else

            su0 += 1.0

        end

    end

    n, su0, su1

end

@inline function _sum_all1(sc_, ex, is_)

    n = length(sc_)

    su = su1 = 0.0

    for id in 1:n

        ab = _index_absolute_exponentiate(sc_, id, ex)

        su += ab

        if is_[id]

            su1 += ab

        end

    end

    n, su, su1

end

@inline function _ready_ks(sc_, ex, is_, mo_)

    n, su0, su1 = _sum_01(sc_, ex, is_)

    n, 1 / su0, su1, 0.0, !isnothing(mo_)

end

@inline function _ready_kli(sc_, ex, is_, mo_)

    n, su, su1 = _sum_all1(sc_, ex, is_)

    ep = eps()

    n, su, su1, ep, ep, ep, 1.0, 1.0, 0.0, !isnothing(mo_), 0.0, 0.0

end

function _enrich!(::KS, sc_, ex, is_, mo_)

    n, de, su1, cu, mo = _ready_ks(sc_, ex, is_, mo_)

    et = 0.0

    eta = 0.0

    for id in 1:n

        if is_[id]

            cu += _index_absolute_exponentiate(sc_, id, ex) / su1

        else

            cu -= de

        end

        if mo

            mo_[id] = cu

        end

        cua = abs(cu)

        if eta < cua

            et = cu

            eta = cua

        end

    end

    et

end

function _enrich!(::KSa, sc_, ex, is_, mo_)

    n, de, su1, cu, mo = _ready_ks(sc_, ex, is_, mo_)

    ar = 0.0

    for id in 1:n

        if is_[id]

            cu += _index_absolute_exponentiate(sc_, id, ex) / su1

        else

            cu -= de

        end

        if mo

            mo_[id] = cu

        end

        ar += cu

    end

    ar / n

end

function _enrich!(::KLi, sc_, ex, is_, mo_)

    n, su, su1, ep, ri, ri1, le, le1, ar, mo, ridp, ri1dp = _ready_kli(sc_, ex, is_, mo_)

    for id in 1:n

        ab = _index_absolute_exponentiate(sc_, id, ex)

        rid = ab / su

        if is_[id]

            ri1d = ab / su1

        else

            ri1d = 0.0

        end

        ri += rid

        ri1 += ri1d

        le -= ridp

        le1 -= ri1dp

        if le < ep

            le = ep

        end

        if le1 < ep

            le1 = ep

        end

        en = BioLab.Information.get_antisymmetric_kullback_leibler_divergence(ri1, ri, le1, le)

        ar += en

        if mo

            mo_[id] = en

        end

        ridp = rid

        ri1dp = ri1d

    end

    ar / n

end

function _enrich!(fu::Function, sc_, ex, is_, mo_)

    n, su, su1, ep, ri, ri1, le, le1, ar, mo, ridp, ri1dp = _ready_kli(sc_, ex, is_, mo_)

    su0 = su - su1

    ri0 = ep

    le0 = 1.0

    ri0dp = 0.0

    for id in 1:n

        ab = _index_absolute_exponentiate(sc_, id, ex)

        rid = ab / su

        if is_[id]

            ri1d = ab / su1

            ri0d = 0.0

        else

            ri1d = 0.0

            ri0d = ab / su0

        end

        ri += rid

        ri1 += ri1d

        ri0 += ri0d

        le -= ridp

        le1 -= ri1dp

        le0 -= ri0dp

        if le < ep

            le = ep

        end

        if le1 < ep

            le1 = ep

        end

        if le0 < ep

            le0 = ep

        end

        en = fu(ri1, ri0, ri) - fu(le1, le0, le)

        ar += en

        if mo

            mo_[id] = en

        end

        ridp = rid

        ri1dp = ri1d

        ri0dp = ri0d

    end

    ar / n

end

function _enrich!(::KLioM, sc_, ex, is_, mo_)

    _enrich!(BioLab.Information.get_antisymmetric_kullback_leibler_divergence, sc_, ex, is_, mo_)

end

function _enrich!(::KLioP, sc_, ex, is_, mo_)

    _enrich!(BioLab.Information.get_symmetric_kullback_leibler_divergence, sc_, ex, is_, mo_)

end

function plot(
    ht,
    al,
    fe_,
    sc_,
    fe1_;
    ex = 1,
    title_text = "Set Enrichment",
    naf = "Feature",
    nas = "Score",
    nal = "Low",
    nah = "High",
)

    is_ = in(Set(fe1_)).(fe_)

    n = length(fe_)

    mo_ = Vector{Float64}(undef, n)

    en = _enrich!(al, sc_, ex, is_, mo_)

    x = collect(1:n)

    scatter = Dict(
        "x" => x,
        "text" => fe_,
        "mode" => "lines",
        "line" => Dict("width" => 0.4),
        "fill" => "tozeroy",
    )

    cor = "#ff1992"

    cob = "#1993ff"

    coe1 = "#07fa07"

    coe2 = BioLab.Plot.add_alpha(coe1, 0.32)

    title_text = BioLab.String.limit(title_text, 80)

    yaxis1_domain = (0, 0.24)

    yaxis2_domain = (0.25, 0.31)

    yaxis3_domain = (0.32, 1)

    annotation = Dict("showarrow" => false, "bgcolor" => "#fcfcfc", "borderwidth" => 2.4)

    enf = BioLab.String.format(en)

    annotationhl = BioLab.Dict.merge_recursively(
        annotation,
        Dict(
            "y" => 0,
            "font" => Dict("size" => 16),
            "borderpad" => 4.8,
            "bordercolor" => "#ffd96a",
        ),
    )

    margin = n * 0.008

    BioLab.Plot.plot(
        ht,
        [
            BioLab.Dict.merge_recursively(
                scatter,
                Dict(
                    "name" => "- Score",
                    "y" => ifelse.(sc_ .< 0, sc_, 0),
                    "line" => Dict("color" => cob),
                    "fillcolor" => cob,
                ),
            ),
            BioLab.Dict.merge_recursively(
                scatter,
                Dict(
                    "name" => "+ Score",
                    "y" => ifelse.(0 .< sc_, sc_, 0),
                    "line" => Dict("color" => cor),
                    "fillcolor" => cor,
                ),
            ),
            Dict(
                "yaxis" => "y2",
                "name" => "Set",
                "y" => zeros(sum(is_)),
                "x" => view(x, is_),
                "text" => view(fe_, is_),
                "hoverinfo" => "name+x+text",
                "mode" => "markers",
                "marker" => Dict(
                    "symbol" => "line-ns",
                    "size" => 24,
                    "line" => Dict("width" => 1.28, "color" => "#175e54", "opacity" => 0.8),
                ),
            ),
            BioLab.Dict.merge_recursively(
                scatter,
                Dict(
                    "yaxis" => "y3",
                    "name" => "Δ Enrichment",
                    "y" => mo_,
                    "line" => Dict("width" => 3.2, "color" => coe1),
                    "fillcolor" => coe2,
                ),
            ),
        ],
        Dict(
            "showlegend" => false,
            "title" => Dict(
                "text" => "<b>$title_text</b>",
                "font" => Dict("size" => 32, "family" => "Relaway", "color" => "#2b2028"),
            ),
            "yaxis" => BioLab.Dict.merge_recursively(
                BioLab.Plot.AXIS,
                Dict("domain" => yaxis1_domain, "title" => Dict("text" => "<b>$nas</b>")),
            ),
            "yaxis2" => BioLab.Dict.merge_recursively(
                BioLab.Plot.AXIS,
                Dict(
                    "domain" => yaxis2_domain,
                    "title" => Dict("text" => "<b>Set</b>"),
                    "tickvals" => (),
                ),
            ),
            "yaxis3" => BioLab.Dict.merge_recursively(
                BioLab.Plot.AXIS,
                Dict("domain" => yaxis3_domain, "title" => Dict("text" => "<b>Δ Enrichment</b>")),
            ),
            "xaxis" => reduce(
                BioLab.Dict.merge_recursively,
                (
                    BioLab.Plot.AXIS,
                    BioLab.Plot.SPIKE,
                    Dict("zeroline" => false, "title" => Dict("text" => "<b>$naf (n=$n)</b>")),
                ),
            ),
            "annotations" => (
                BioLab.Dict.merge_recursively(
                    annotation,
                    Dict(
                        "yref" => "paper",
                        "xref" => "paper",
                        "y" => 1.04,
                        "text" => "Enrichment = <b>$enf</b>",
                        "font" => Dict("size" => 20, "color" => "#224634"),
                        "borderpad" => 12.8,
                        "bordercolor" => "#ffd96a",
                    ),
                ),
                BioLab.Dict.merge_recursively(
                    annotationhl,
                    Dict(
                        "x" => 1 - margin,
                        "xanchor" => "right",
                        "text" => nah,
                        "font" => Dict("color" => cor),
                    ),
                ),
                BioLab.Dict.merge_recursively(
                    annotationhl,
                    Dict(
                        "x" => n + margin,
                        "xanchor" => "left",
                        "text" => nal,
                        "font" => Dict("color" => cob),
                    ),
                ),
            ),
        ),
    )

end

function enrich(al, fe_, sc_::AbstractVector, fe1___; mi = 1, ex = 1)

    en_ = Vector{Float64}(undef, length(fe1___))

    fe_id = Dict(fe => id for (id, fe) in enumerate(fe_))

    for (id, fe1_) in enumerate(fe1___)

        is_ = BioLab.Dict.is_in(fe_id, fe1_)

        if sum(is_) < mi

            en = NaN

        else

            en = _enrich!(al, sc_, ex, is_, nothing)

        end

        en_[id] = en

    end

    en_

end

function enrich(al, fe_, fe_x_sa_x_sc::AbstractMatrix, fe1___; mi = 1, ex = 1)

    se_x_sa_x_en = Matrix{Float64}(undef, length(fe1___), size(fe_x_sa_x_sc, 2))

    no_ = BitVector(undef, length(fe_))

    @showprogress for (id, sc_) in enumerate(eachcol(fe_x_sa_x_sc))

        no_ .= .!isnan.(sc_)

        scn_ = sc_[no_]

        so_ = sortperm(scn_; rev = true)

        se_x_sa_x_en[:, id] = enrich(al, view(fe_[no_], so_), view(scn_, so_), fe1___; mi, ex)

    end

    se_x_sa_x_en

end

function plot(di, al, fe_, fe_x_sa_x_sc, fe1___, nac, se_, sa_, se_x_sa_x_en; ex = 1, n_pl = 4)

    BioLab.Error.error_missing(di)

    als = _make_string(al)

    nacc = BioLab.Path.clean(nac)

    BioLab.Plot.plot_heat_map(
        joinpath(di, "set_x_$(nacc)_x_enrichment.html"),
        se_x_sa_x_en,
        se_,
        sa_;
        nar = "Set",
        nac,
        layout = Dict("title" => Dict("text" => "Enrichment using $als")),
    )

    nos_ = BitVector(undef, length(fe_))

    noe = .!isnan.(se_x_sa_x_en)

    for id_ in
        view(view(CartesianIndices(se_x_sa_x_en), noe), sortperm(view(se_x_sa_x_en, noe)))[BioLab.Rank.get_extreme(
        sum(noe),
        n_pl,
    )]

        id1, id2 = Tuple(id_)

        enf = BioLab.String.format(se_x_sa_x_en[id1, id2])

        se = se_[id1]

        sa = sa_[id2]

        sc_ = fe_x_sa_x_sc[:, id2]

        nos_ .= .!isnan.(sc_)

        scn_ = view(sc_, nos_)

        so_ = sortperm(scn_; rev = true)

        plot(
            joinpath(di, "$(sa)_enriching_$se.html"),
            al,
            view(view(fe_, nos_), so_),
            view(scn_, so_),
            fe1___[id1];
            ex,
            title_text = "$sa Enriching $se ($enf)",
        )

    end

    di

end

end
