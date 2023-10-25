module FeatureSetEnrichment

using Printf: @sprintf

using ProgressMeter: @showprogress

using Nucleus

struct KS end

struct KSa end

struct KLi1 end

struct KLi end

struct KLioM end

struct KLioP end

function _make_string(al)

    string(al)[27:(end - 2)]

end

@inline function _absolute_exponentiate(sc, ex)

    ab = abs(sc)

    if !isone(ex)

        ab ^= ex

    end

    ab

end

@inline function _get_0_1_normalizer(sc_, ex, is_)

    n = lastindex(sc_)

    su0 = su1 = 0.0

    for id in 1:n

        if is_[id]

            su1 += _absolute_exponentiate(sc_[id], ex)

        else

            su0 += 1.0

        end

    end

    n, 1 / su0, 1 / su1

end

# TODO: Benchmark.
@inline function _get_1_normalizer(sc_, ex, is_)

    n = lastindex(sc_)

    su1 = 0.0

    for id in 1:n

        if is_[id]

            su1 += _absolute_exponentiate(sc_[id], ex)

        end

    end

    n, 1 / su1

end

@inline function _get_all_1_normalizer(sc_, ex, is_)

    n = lastindex(sc_)

    su = su1 = 0.0

    for id in 1:n

        ab = _absolute_exponentiate(sc_[id], ex)

        su += ab

        if is_[id]

            su1 += ab

        end

    end

    n, 1 / su, 1 / su1

end

@inline function _get_0_normalizer(noa, no1)

    1 / (1 / noa - 1 / no1)

end

@inline function _floor(nu)

    ep = eps()

    if nu < ep

        ep

    else

        nu

    end

end

function _enrich!(::KS, sc_, ex, is_, mo_)

    n, no0, no1 = _get_0_1_normalizer(sc_, ex, is_)

    cu = 0.0

    et = 0.0

    eta = 0.0

    for id in 1:n

        if is_[id]

            cu += _absolute_exponentiate(sc_[id], ex) * no1

        else

            cu -= no0

        end

        if !isnothing(mo_)

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

    n, no0, no1 = _get_0_1_normalizer(sc_, ex, is_)

    cu = 0.0

    ar = 0.0

    for id in 1:n

        if is_[id]

            cu += _absolute_exponentiate(sc_[id], ex) * no1

        else

            cu -= no0

        end

        if !isnothing(mo_)

            mo_[id] = cu

        end

        ar += cu

    end

    ar / n

end

function _enrich!(::KLi1, sc_, ex, is_, mo_)

    n, no1 = _get_1_normalizer(sc_, ex, is_)

    ri = ri1 = eps()

    rid = 1 / n

    le = 1.0 + 1 / n

    le1 = 1.0

    ar = 0.0

    pr1 = 0.0

    for id in 1:n

        if is_[id]

            ri1d = _absolute_exponentiate(sc_[id], ex) * no1

        else

            ri1d = 0.0

        end

        en = Nucleus.Information.get_antisymmetric_kullback_leibler_divergence(
            ri1 += ri1d,
            _floor(le1 -= pr1),
            ri += rid,
            _floor(le -= rid),
        )

        ar += en

        if !isnothing(mo_)

            mo_[id] = en

        end

        pr1 = ri1d

    end

    ar / n

end

function _enrich!(::KLi, sc_, ex, is_, mo_)

    n, noa, no1 = _get_all_1_normalizer(sc_, ex, is_)

    ri = ri1 = eps()

    le = le1 = 1.0

    ar = 0.0

    pr = pr1 = 0.0

    for id in 1:n

        ab = _absolute_exponentiate(sc_[id], ex)

        if is_[id]

            ri1d = ab * no1

        else

            ri1d = 0.0

        end

        rid = ab * noa

        en = Nucleus.Information.get_antisymmetric_kullback_leibler_divergence(
            ri1 += ri1d,
            _floor(le1 -= pr1),
            ri += rid,
            _floor(le -= pr),
        )

        ar += en

        if !isnothing(mo_)

            mo_[id] = en

        end

        pr = rid

        pr1 = ri1d

    end

    ar / n

end

function _enrich!(::KLioM, sc_, ex, is_, mo_)

    n, noa, no1 = _get_all_1_normalizer(sc_, ex, is_)

    no0 = _get_0_normalizer(noa, no1)

    ri = ri1 = ri0 = eps()

    le = le1 = le0 = 1.0

    ar = 0.0

    pr = pr1 = pr0 = 0.0

    for id in 1:n

        ab = _absolute_exponentiate(sc_[id], ex)

        rid = ab * noa

        ri += rid

        le -= pr

        if is_[id]

            ri1d = ab * no1

            ri0d = 0.0

        else

            ri1d = 0.0

            ri0d = ab * no0

        end

        ri1 += ri1d

        le1 -= pr1

        ri0 += ri0d

        le0 -= pr0

        en =
            Nucleus.Information.get_antisymmetric_kullback_leibler_divergence(ri1, ri0, ri) -
            Nucleus.Information.get_antisymmetric_kullback_leibler_divergence(
                _floor(le1),
                _floor(le0),
                _floor(le),
            )

        ar += en

        if !isnothing(mo_)

            mo_[id] = en

        end

        pr = rid

        pr1 = ri1d

        pr0 = ri0d

    end

    ar / n

end

function _enrich!(::KLioP, sc_, ex, is_, mo_)

    n, noa, no1 = _get_all_1_normalizer(sc_, ex, is_)

    no0 = _get_0_normalizer(noa, no1)

    ri = ri1 = ri0 = eps()

    le = le1 = le0 = 1.0

    ar = 0.0

    pr = pr1 = pr0 = 0.0

    for id in 1:n

        ab = _absolute_exponentiate(sc_[id], ex)

        rid = ab * noa

        ri += rid

        le -= pr

        if is_[id]

            ri1d = ab * no1

            ri0d = 0.0

        else

            ri1d = 0.0

            ri0d = ab * no0

        end

        ri1 += ri1d

        le1 -= pr1

        ri0 += ri0d

        le0 -= pr0

        en =
            Nucleus.Information.get_symmetric_kullback_leibler_divergence(ri1, ri0, ri) -
            Nucleus.Information.get_symmetric_kullback_leibler_divergence(
                _floor(le1),
                _floor(le0),
                _floor(le),
            )

        ar += en

        if !isnothing(mo_)

            mo_[id] = en

        end

        pr = rid

        pr1 = ri1d

        pr0 = ri0d

    end

    ar / n

end

function _get_extreme(nu_)

    mi = minimum(nu_)

    ma = maximum(nu_)

    mia = abs(mi)

    maa = abs(ma)

    if isapprox(mia, maa)

        (mi, ma)

    elseif maa < mia

        (mi,)

    else

        (ma,)

    end

end

function _format(nu)

    @sprintf "%.4g" nu

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

    n = lastindex(fe_)

    x = collect(1:n)

    is_ = in(Set(fe1_)).(fe_)

    mo_ = Vector{Float64}(undef, n)

    en = _enrich!(al, sc_, ex, is_, mo_)

    coe1 = "#07fa07"

    coe2 = Nucleus.Color.add_alpha(coe1, 0.32)

    scatter = Dict("x" => x, "text" => fe_, "mode" => "lines", "fill" => "tozeroy")

    if typeof(al) == KS

        fi = Dict("fillcolor" => "#ffffff")

        id_ = findall(in(_get_extreme(mo_)), mo_)

        pe_ = (
            Dict(
                "yaxis" => "y3",
                "y" => mo_[id_],
                "x" => x[id_],
                "mode" => "markers",
                "marker" => Dict(
                    #"symbol" => "circle",
                    "size" => 32,
                    "color" => coe2,
                    "opacity" => 0.72,
                    "line" => Dict("width" => 2, "color" => Nucleus.Color.HEFA),
                ),
            ),
        )

    else

        fi = Dict("fillcolor" => coe2)

        pe_ = ()

    end

    annotation = Dict("showarrow" => false, "bgcolor" => "#fcfcfc", "borderwidth" => 2.4)

    annotationhl = merge(
        annotation,
        Dict(
            "y" => 0,
            "font" => Dict("size" => 16),
            "borderpad" => 4.8,
            "bordercolor" => Nucleus.Color.HEAY,
        ),
    )

    margin = n * 0.008

    Nucleus.Plot.plot(
        ht,
        [
            merge(
                scatter,
                Dict(
                    "name" => "- Score",
                    "y" => ifelse.(sc_ .< 0, sc_, 0),
                    "line" => Dict("width" => 0.4, "color" => Nucleus.Color.HEBL),
                    "fillcolor" => Nucleus.Color.HEBL,
                ),
            ),
            merge(
                scatter,
                Dict(
                    "name" => "+ Score",
                    "y" => ifelse.(0 .< sc_, sc_, 0),
                    "line" => Dict("width" => 0.4, "color" => Nucleus.Color.HERE),
                    "fillcolor" => Nucleus.Color.HERE,
                ),
            ),
            Dict(
                "name" => "Set",
                "yaxis" => "y2",
                "y" => zeros(sum(is_)),
                "x" => view(x, is_),
                "text" => view(fe_, is_),
                "mode" => "markers",
                "marker" => Dict(
                    "symbol" => "line-ns",
                    "size" => 24,
                    "line" => Dict("width" => 1.28, "color" => "#175e54", "opacity" => 0.8),
                ),
            ),
            merge(
                scatter,
                Dict(
                    "name" => "Δ Enrichment",
                    "yaxis" => "y3",
                    "y" => mo_,
                    "line" => Dict("width" => 3.2, "color" => coe1),
                ),
                fi,
            ),
            pe_...,
        ],
        Dict(
            "showlegend" => false,
            "title" => Dict(
                "text" => "<b>$(Nucleus.String.limit(title_text, 80))</b>",
                "font" => Dict("size" => 32, "family" => "Relaway", "color" => "#2b2028"),
            ),
            "yaxis" => Dict(
                "domain" => (0, 0.24),
                "title" => Dict("text" => "<b>$nas</b>"),
                "showgrid" => false,
            ),
            "yaxis2" => Dict(
                "domain" => (0.25, 0.31),
                "title" => Dict("text" => "<b>Set</b>"),
                "tickvals" => (),
                "showgrid" => false,
            ),
            "yaxis3" => Dict(
                "domain" => (0.32, 1),
                "title" => Dict("text" => "<b>Δ Enrichment</b>"),
                "showgrid" => false,
            ),
            "xaxis" => merge(
                Dict(
                    "title" => Dict("text" => "<b>$naf ($n)</b>"),
                    "zeroline" => false,
                    "showgrid" => false,
                ),
                Nucleus.Plot.SPIKE,
            ),
            "annotations" => (
                merge(
                    annotation,
                    Dict(
                        "yref" => "paper",
                        "xref" => "paper",
                        "y" => 1.04,
                        "text" => "Enrichment = <b>$(_format(en))</b>",
                        "font" => Dict("size" => 20, "color" => "#224634"),
                        "borderpad" => 12.8,
                        "bordercolor" => Nucleus.Color.HEAY,
                    ),
                ),
                merge(
                    annotationhl,
                    Dict(
                        "x" => 1 - margin,
                        "xanchor" => "right",
                        "text" => nah,
                        "font" => Dict("color" => Nucleus.Color.HERE),
                    ),
                ),
                merge(
                    annotationhl,
                    Dict(
                        "x" => n + margin,
                        "xanchor" => "left",
                        "text" => nal,
                        "font" => Dict("color" => Nucleus.Color.HEBL),
                    ),
                ),
            ),
        ),
    )

end

function enrich(al, fe_, sc_::AbstractVector, fe1___; mi = 1, ex = 1)

    en_ = Vector{Float64}(undef, lastindex(fe1___))

    fe_id = Dict(fe => id for (id, fe) in enumerate(fe_))

    for (id, fe1_) in enumerate(fe1___)

        is_ = Nucleus.Dict.is_in(fe_id, fe1_)

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

    se_x_sa_x_en = Matrix{Float64}(undef, lastindex(fe1___), size(fe_x_sa_x_sc, 2))

    no_ = BitVector(undef, lastindex(fe_))

    @showprogress for (id, sc_) in enumerate(eachcol(fe_x_sa_x_sc))

        no_ .= .!isnan.(sc_)

        scn_ = view(sc_, no_)

        so_ = sortperm(scn_; rev = true)

        se_x_sa_x_en[:, id] = enrich(al, view(fe_, no_)[so_], scn_[so_], fe1___; mi, ex)

    end

    se_x_sa_x_en

end

function plot(di, al, fe_, fe_x_sa_x_sc, fe1___, nac, se_, sa_, se_x_sa_x_en; ex = 1, n_pl = 4)

    Nucleus.Error.error_missing(di)

    Nucleus.Plot.plot_heat_map(
        joinpath(di, "set_x_$(Nucleus.Path.clean(nac))_x_enrichment.html"),
        se_x_sa_x_en;
        y = se_,
        x = sa_,
        nar = "Set",
        nac,
        layout = Dict("title" => Dict("text" => "Enrichment using $(_make_string(al))")),
    )

    noe = .!isnan.(se_x_sa_x_en)

    no_ = BitVector(undef, lastindex(fe_))

    for id_ in
        view(view(CartesianIndices(se_x_sa_x_en), noe), sortperm(view(se_x_sa_x_en, noe)))[Nucleus.Rank.get_extreme(
        sum(noe),
        n_pl,
    )]

        id1, id2 = Tuple(id_)

        sc_ = view(fe_x_sa_x_sc, :, id2)

        no_ .= .!isnan.(sc_)

        scn_ = view(sc_, no_)

        so_ = sortperm(scn_; rev = true)

        pr = "$(sa_[id2]) Enriching $(se_[id1])"

        plot(
            joinpath(di, "$(Nucleus.Path.clean(pr)).html"),
            al,
            view(fe_, no_)[so_],
            scn_[so_],
            fe1___[id1];
            ex,
            title_text = "$pr ($(_format(se_x_sa_x_en[id1, id2])))",
        )

    end

end

end
