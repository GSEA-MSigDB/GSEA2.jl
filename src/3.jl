function _get_extreme(nu_)

    mi, ma = extrema(nu_)

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

    coe2 = Omics.Color.add_alpha(coe1, 0.32)

    scatter = Dict("x" => x, "text" => fe_, "mode" => "lines", "fill" => "tozeroy")

    ks = typeof(al) == KS

    data = [
        merge(
            scatter,
            Dict(
                "name" => "- Score",
                "y" => (sc -> sc < 0 ? sc : 0).(sc_),
                "line" => Dict("width" => 0.4, "color" => Omics.Color.HEBL),
                "fillcolor" => Omics.Color.HEBL,
            ),
        ),
        merge(
            scatter,
            Dict(
                "name" => "+ Score",
                "y" => (sc -> 0 < sc ? sc : 0).(sc_),
                "line" => Dict("width" => 0.4, "color" => Omics.Color.HERE),
                "fillcolor" => Omics.Color.HERE,
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
                "fillcolor" => ks ? "#ffffff" : coe2,
            ),
        ),
    ]

    if ks

        id_ = findall(in(_get_extreme(mo_)), mo_)

        push!(
            data,
            Dict(
                "yaxis" => "y3",
                "y" => mo_[id_],
                "x" => x[id_],
                "mode" => "markers",
                "marker" => Dict(
                    "size" => 32,
                    "color" => coe2,
                    "opacity" => 0.72,
                    "line" => Dict("width" => 2, "color" => Omics.Color.HEFA),
                ),
            ),
        )

    end

    annotation = Dict("showarrow" => false, "bgcolor" => "#fcfcfc", "borderwidth" => 2.4)

    annotationhl = merge(
        annotation,
        Dict(
            "y" => 0,
            "font" => Dict("size" => 16),
            "borderpad" => 4.8,
            "bordercolor" => Omics.Color.HEAY,
        ),
    )

    margin = n * 0.008

    Omics.Plot.plot(
        ht,
        data,
        Dict(
            "showlegend" => false,
            "title" => Dict(
                "text" => "<b>$(Omics.String.limit(title_text, 80))</b>",
                "font" =>
                    Dict("size" => 32, "family" => "Relaway", "color" => "#2b2028"),
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
            "xaxis" => Dict(
                "title" => Dict("text" => "<b>$naf ($n)</b>"),
                "zeroline" => false,
                "showgrid" => false,
                "showspikes" => true,
                "spikesnap" => "cursor",
                "spikemode" => "across",
                "spikedash" => "solid",
                "spikethickness" => 1,
                "spikecolor" => "#561649",
            ),
            "annotations" => (
                merge(
                    annotation,
                    Dict(
                        "yref" => "paper",
                        "xref" => "paper",
                        "y" => 1.04,
                        "text" => "Enrichment = <b>$(Omics.Number.format4(en))</b>",
                        "font" => Dict("size" => 20, "color" => "#224634"),
                        "borderpad" => 12.8,
                        "bordercolor" => Omics.Color.HEAY,
                    ),
                ),
                merge(
                    annotationhl,
                    Dict(
                        "x" => 1 - margin,
                        "xanchor" => "right",
                        "text" => nah,
                        "font" => Dict("color" => Omics.Color.HERE),
                    ),
                ),
                merge(
                    annotationhl,
                    Dict(
                        "x" => n + margin,
                        "xanchor" => "left",
                        "text" => nal,
                        "font" => Dict("color" => Omics.Color.HEBL),
                    ),
                ),
            ),
        ),
    )

end
