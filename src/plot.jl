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

    tr = Dict("mode" => "lines", "line" => Dict("width" => 0), "fill" => "tozeroy")

    ie_ = findall(<(0.0), sc_)

    ip_ = findall(>=(0.0), sc_)

    da_ = [
        merge(tr, Dict("y" => sc_[ie_], "x" => xc_[ie_], "fillcolor" => Omics.Color.BL)),
        merge(tr, Dict("y" => sc_[ip_], "x" => xc_[ip_], "fillcolor" => Omics.Color.RE)),
        Dict(
            "yaxis" => "y2",
            "y" => zeros(sum(is_)),
            "x" => xc_[is_],
            "mode" => "markers",
            "marker" => Dict(
                "symbol" => "line-ns",
                "size" => 24,
                "line" => Dict(
                    "width" => 1.32,
                    "color" => Omics.Color.hexify(Omics.Color.SG, 0.8),
                ),
            ),
        ),
        merge(tr, Dict("yaxis" => "y3", "y" => mo_, "x" => xc_, "fillcolor" => "#07fa07")),
    ]

    if typeof(al) == KS

        ix_ = findall(in(extrema(mo_)), mo_)

        push!(
            da_,
            Dict(
                "yaxis" => "y3",
                "y" => mo_[ix_],
                "x" => xc_[ix_],
                "mode" => "markers",
                "marker" =>
                    Dict("size" => 32, "color" => Omics.Color.hexify(Omics.Color.YE, 0.96)),
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
                "yaxis" => Dict("domain" => (0, 0.24), "title" => Dict("text" => "$ns")),
                "yaxis2" => Dict(
                    "domain" => (0.25, 0.31),
                    "title" => Dict("text" => "Set"),
                    "tickvals" => (),
                ),
                "yaxis3" =>
                    Dict("domain" => (0.32, 1), "title" => Dict("text" => "Δ Enrichment")),
                "xaxis" => Dict(
                    "title" => Dict("text" => "$nf ($uf)"),
                    "showspikes" => true,
                    "spikesnap" => "cursor",
                    "spikemode" => "across",
                    "spikedash" => "solid",
                    "spikethickness" => 1,
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
