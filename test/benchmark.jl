using Test: @test

using Printf: @sprintf

using Nucleus

using GSEA

# ---- #

const DIP = joinpath(homedir(), "Desktop", "benchmark")

# ---- #

const DIJ = joinpath(DIP, "json_files")

# ---- #

const DID = joinpath(DIP, "Datasets_and_Phenotypes")

# ---- #

const DIS = joinpath(DIP, "Gene_Sets_Collections")

# ---- #

const DIR = joinpath(DIP, "results")

# ---- #

const AL_ = ("ks", "kli", "kli1", "kliom", "kliop")

# ---- #

const DIB = Nucleus.Path.establish(joinpath(dirname(@__DIR__), "benchmark"))

# ---- #

const TE_ = String[]

# ---- #

const PE_ = Float64[]

# ---- #

function log(te, is_, py, ju)

    if any(is_)

        pe = 100 - 100sum(is_) / size(py, 1)

        pyi = view(py, is_, :)

        jui = view(ju, is_, :)

        @warn "$te $(@sprintf "%.5g" pe)%" pyi jui

    else

        pe = 100.0

        pyi = nothing

        jui = nothing

    end

    push!(TE_, te)

    push!(PE_, pe)

    nothing

end

# ---- #

function compare(jsk, al, st, py, pyi, ju, jui)

    log("$(Nucleus.String.clean(jsk)) $al $st", py[!, pyi] .!= ju[!, jui], py, ju)

end

# ---- #

function compare(jsk, al, st, py, pyi, ju, jui, atol)

    log(
        "$(Nucleus.String.clean(jsk)) $al $st",
        .!isapprox.(
            # TODO
            (fl -> fl isa AbstractString ? parse(Float64, lstrip(fl, '≤')) : fl).(py[!, pyi]),
            ju[!, jui];
            atol,
        ),
        py,
        ju,
    )

end

# ---- #

for (idb, js) in enumerate(Nucleus.Path.read(DIJ))

    @info "$idb $js"

    jsk = view(js, 1:(lastindex(js) - 5))

    ke_va = Nucleus.Dict.read(joinpath(DIJ, js))[jsk]

    di = Nucleus.Path.establish(joinpath(DIB, Nucleus.Path.clean(jsk)))

    dii = Nucleus.Path.establish(joinpath(di, "input"))

    tst = joinpath(dii, "target_x_sample_x_number.tsv")

    tsf = joinpath(dii, "feature_x_sample_x_number.tsv")

    if !isfile(tst) || !isfile(tsf)

        GSEA.convert_cls_gct(tst, tsf, joinpath(DID, ke_va["cls"]), joinpath(DID, ke_va["ds"]))

    end

    jss = joinpath(dii, "set_features.json")

    if !isfile(jss)

        GSEA.convert_gmt(jss, (joinpath(DIS, gm) for gm in ke_va["gene_sets_collections"])...)

    end

    dir = joinpath(DIR, basename(ke_va["results_directory"]))

    minimum_set_size = ke_va["min_gene_set_size"]

    maximum_set_size = ke_va["max_gene_set_size"]

    for (al, pr, no, me) in zip(
        AL_,
        ke_va["results_files_prefix"],
        ke_va["standardize_genes_before_gene_sel"],
        ke_va["gene_selection_metric"],
    )

        @info al no me

        dio = Nucleus.Path.establish(joinpath(di, "output_$al"))

        txm = joinpath(dir, "$(pr)_gene_selection_scores.txt")

        txs = joinpath(dir, "$(pr)_GSEA_results_table.txt")

        tsm = joinpath(dio, "feature_x_metric_x_score.tsv")

        tss = joinpath(dio, "set_x_statistic_x_number.tsv")

        if !isfile(tsm)

            if endswith(me, "diff_means2")

                metric = "mean-difference"

            elseif endswith(me, "signal_to_noise_GSEA4")

                metric = "signal-to-noise-ratio"

            end

            GSEA.metric_rank(
                dio,
                tst,
                tsf,
                jss;
                minimum_set_size,
                maximum_set_size,
                normalization_dimension = Int(no),
                normalization_standard_deviation = 3.0,
                metric,
                algorithm = al,
                number_of_permutations = 0,
                number_of_sets_to_plot = 0,
            )

            rm(tss)

        end

        py = Nucleus.DataFrame.read(txm; select = [1, 2])

        ju = Nucleus.DataFrame.read(tsm)

        sort!(py)

        sort!(ju)

        compare(jsk, al, "Gene", py, 1, ju, 1)

        compare(jsk, al, "Metric", py, 2, ju, 2, 0.0001)

        if !isfile(tss)

            GSEA.user_rank(
                dio,
                txm,
                jss;
                minimum_set_size,
                maximum_set_size,
                algorithm = al,
                permutation = joinpath(dir, "$(pr)_rand_perm_gene_scores.txt"),
                number_of_sets_to_plot = 0,
            )

        end

        py = sort!(Nucleus.DataFrame.read(txs; select = [1, 5, 4, 6, 7]))

        ju = sort!(Nucleus.DataFrame.read(tss))

        compare(jsk, al, "Set", py, 1, ju, 1)

        compare(jsk, al, "Enrichment", py, 3, ju, 2, 0.01)

        compare(jsk, al, "Normalized Enrichment", py, 2, ju, 3, 0.01)

        compare(jsk, al, "P-Value", py, 4, ju, 4, 0.01)

        compare(jsk, al, "Adjusted P-Value", py, 5, ju, 5, 0.01)

    end

end

# ---- #

const ST_ = ("Metric", "Enrichment", "Normalized Enrichment", "P-Value", "Adjusted P-Value")

# ---- #

for al in AL_

    is___ = (findall(endswith("$al $st"), TE_) for st in ST_)

    Nucleus.Plot.plot_bar(
        joinpath(DIB, "summary_$al.html"),
        Tuple(view(PE_, is_) for is_ in is___),
        Tuple(Nucleus.String.split_get.(view(TE_, is_), ' ', 1) for is_ in is___),
        name_ = ST_,
        layout = Dict(
            "title" => Dict("text" => al),
            "yaxis" => Dict("title" => Dict("text" => "% Match")),
            "xaxis" => Dict("title" => Dict("text" => "Benchmark")),
        ),
    )

end
