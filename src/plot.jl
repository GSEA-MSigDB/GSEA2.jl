using Printf: @sprintf

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

    is_ = map(in(Set(me_)), fe_)

    mo_ = Vector{Float64}(undef, uf)

    en = _enrich!(al, sc_, ex, is_, mo_)

    he = "#07fa07"

    hx = Omics.Color.hexify(he, 0.32)

    tr = Dict("x" => xc_, "text" => fe_, "mode" => "lines", "fill" => "tozeroy")

    ks = typeof(al) == KS

    da_ = [
        merge(
            tr,
            Dict(
                "name" => "- Score",
                "y" => (sc -> sc < 0 ? sc : 0).(sc_),
                "line" => Dict("width" => 0.4, "color" => Omics.Color.BL),
                "fillcolor" => Omics.Color.BL,
            ),
        ),
        merge(
            tr,
            Dict(
                "name" => "+ Score",
                "y" => (sc -> 0 < sc ? sc : 0).(sc_),
                "line" => Dict("width" => 0.4, "color" => Omics.Color.RE),
                "fillcolor" => Omics.Color.RE,
            ),
        ),
        Dict(
            "name" => "Set",
            "yaxis" => "y2",
            "y" => zeros(sum(is_)),
            "x" => view(xc_, is_),
            "text" => view(fe_, is_),
            "mode" => "markers",
            "marker" => Dict(
                "symbol" => "line-ns",
                "size" => 24,
                "line" => Dict("width" => 1.28, "color" => "#175e54", "opacity" => 0.8),
            ),
        ),
        merge(
            tr,
            Dict(
                "name" => "Δ Enrichment",
                "yaxis" => "y3",
                "y" => mo_,
                "line" => Dict("width" => 3.2, "color" => he),
                "fillcolor" => ks ? "#ffffff" : hx,
            ),
        ),
    ]

    if ks

        id_ = findall(in(extrema(mo_)), mo_)

        push!(
            da_,
            Dict(
                "yaxis" => "y3",
                "y" => mo_[id_],
                "x" => xc_[id_],
                "mode" => "markers",
                "marker" => Dict(
                    "size" => 32,
                    "color" => hx,
                    "opacity" => 0.72,
                    "line" => Dict("width" => 2, "color" => Omics.Color.LI),
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
            "bordercolor" => Omics.Color.YE,
        ),
    )

    ma = uf * 0.008

    Omics.Plot.plot(
        ht,
        da_,
        Omics.Dic.merg(
            Dict(
                "showlegend" => false,
                "title" => Dict(
                    "font" =>
                        Dict("size" => 32, "family" => "Relaway", "color" => "#2b2028"),
                ),
                "yaxis" => Dict(
                    "domain" => (0, 0.24),
                    "title" => Dict("text" => "<b>$ns</b>"),
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
                    "title" => Dict("text" => "<b>$nf ($uf)</b>"),
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
                            "text" => "Enrichment = <b>$(@sprintf "%.4g" en)</b>",
                            "font" => Dict("size" => 20, "color" => "#224634"),
                            "borderpad" => 12.8,
                            "bordercolor" => Omics.Color.YE,
                        ),
                    ),
                    merge(
                        annotationhl,
                        Dict(
                            "x" => 1 - ma,
                            "xanchor" => "right",
                            "text" => nh,
                            "font" => Dict("color" => Omics.Color.RE),
                        ),
                    ),
                    merge(
                        annotationhl,
                        Dict(
                            "x" => uf + ma,
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
