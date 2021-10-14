using PlotlyJS: Layout, SyncPlot, attr, scatter
using Printf: @sprintf

using Kwat.constant: get_golden_ratio
using Kwat.figure: plot

function plot_mountain(
    fe_::Vector{String},
    sc_::Vector{Float64},
    in_::Vector{Float64},
    en_::Vector{Float64},
    en::Float64;
    height::Int64 = 480,
    title_text::String = "Mountain Plot",
    title_font_size::Int64 = 24,
    axis_title_font_size::Int64 = 12,
    names::String = "Score",
    line_width::Float64 = 2.0,
    si::Bool = true,
    pa::String = "",
)::SyncPlot

    width = height * get_golden_ratio()

    yaxis1_domain = [0.0, 0.24]

    yaxis2_domain = [0.24, 0.32]

    yaxis3_domain = [0.32, 1.0]

    annotation = attr(
        xref = "paper",
        yref = "paper",
        yanchor = "middle",
        showarrow = false,
    )

    annotationy = merge(
        annotation,
        attr(xanchor = "right", x = -0.064, font_size = axis_title_font_size),
    )

    annotationx = merge(annotation, attr(xanchor = "center", x = 0.5))

    namee = "Enrichment"

    en = @sprintf "%.3f" en

    n_fe = length(fe_)

    layout = Layout(
        height = height,
        width = width,
        margin = attr(t = trunc(height * 0.16), l = trunc(width * 0.16)),
        legend = attr(
            orientation = "h",
            yanchor = "middle",
            xanchor = "center",
            y = -0.2,
            x = 0.5,
        ),
        yaxis1 = attr(domain = yaxis1_domain, showline = true),
        yaxis2 = attr(
            domain = yaxis2_domain,
            showticklabels = false,
            showgrid = false,
        ),
        yaxis3 = attr(domain = yaxis3_domain, showline = true),
        xaxis = attr(
            zeroline = false,
            showspikes = true,
            spikethickness = 0.8,
            spikedash = "solid",
            spikemode = "across",
        ),
        annotations = [
            merge(
                annotationx,
                attr(
                    y = 1.16,
                    text = string("<b>", title_text, "</b>"),
                    font_size = title_font_size,
                ),
            ),
            merge(
                annotationx,
                attr(
                    y = 1.04,
                    text = string("<b>Statistic = ", en, "</b>"),
                    font = attr(
                        size = title_font_size * 0.64,
                        color = "2a603b",
                    ),
                ),
            ),
            merge(
                annotationy,
                attr(
                    y = get_center(yaxis1_domain...),
                    text = string("<b>", names, "</b>"),
                ),
            ),
            merge(
                annotationy,
                attr(y = get_center(yaxis2_domain...), text = "<b>Set</b>"),
            ),
            merge(
                annotationy,
                attr(
                    y = get_center(yaxis3_domain...),
                    text = string("<b>", namee, "</b>"),
                ),
            ),
            merge(
                annotationx,
                attr(
                    y = -0.088,
                    text = string("<b>Feature Rank (n=", n_fe, ")</b>"),
                ),
            ),
        ],
    )

    x = 1:n_fe

    tracef = scatter(;
        name = names,
        y = sc_,
        x = x,
        text = fe_,
        mode = "lines",
        line = attr(width = 0, color = "20d8ba"),
        fill = "tozeroy",
        hoverinfo = "x+y+text",
    )

    in_ = BitVector(in_)

    traces = scatter(;
        name = "Set",
        yaxis = "y2",
        y = zeros(Int64(sum(in_))),
        x = x[in_],
        text = fe_[in_],
        mode = "markers",
        marker = attr(
            symbol = "line-ns-open",
            size = height * (yaxis2_domain[2] - yaxis2_domain[1]) * 0.64,
            line_width = line_width,
            color = "9017e6",
        ),
        hoverinfo = "x+text",
    )

    trace_ = [tracef, traces]

    if si

        for (name, is_, color) in [
            ["- Enrichment", en_ .< 0.0, "0088ff"],
            ["+ Enrichment", 0.0 .< en_, "ff1968"],
        ]

            push!(
                trace_,
                scatter(;
                    name = name,
                    yaxis = "y3",
                    y = ifelse.(is_, en_, 0.0),
                    x = x,
                    text = fe_,
                    mode = "lines",
                    line = attr(width = 0.0, color = color),
                    fill = "tozeroy",
                    hoverinfo = "x+y+text",
                ),
            )

        end

    else

        push!(
            trace_,
            scatter(;
                name = namee,
                yaxis = "y3",
                y = en_,
                x = x,
                text = fe_,
                mode = "lines",
                line = attr(width = 0.0, color = "#4e40d8"),
                fill = "tozeroy",
                hoverinfo = "x+y+text",
            ),
        )

    end

    return plot(trace_, layout; pa = pa)

end

export plot_mountain
