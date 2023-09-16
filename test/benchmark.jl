using Test: @test

using BioLab

using GSEA

# ---- #

const DIP = joinpath(homedir(), "Desktop", "benchmark")

const DIJ = joinpath(DIP, "json_files")

const DID = joinpath(DIP, "Datasets_and_Phenotypes")

const DIS = joinpath(DIP, "Gene_Sets_Collections")

const DIR = joinpath(DIP, "results")

const AL_ = ("ks", "kli", "kli", "kliom", "kliop")

# ---- #

function parse_float(fl)

    if fl isa AbstractString

        fl = parse(Float64, lstrip(fl, '≤'))

    end

    fl

end

# ---- #

const DIB = joinpath(dirname(@__DIR__), "benchmark")
BioLab.Path.remake_directory(DIB)

# ---- #

tst = tsf = py = ju = nothing

# ---- #

for (id, js) in enumerate(BioLab.Path.read(DIJ))

    if js in (
        "CCLE_STAT3_vs_mRNA.json", # Number of sets differ.
        "CCLE_YAP_vs_mRNA.json", # Number of sets differ.
        "CFC1-overexpressing NGP cells.json", # Target directionalities differ.
        "CRISPR_FOXA1_vs_mRNA.json", # Number of sets differ.
        "CRISPR_NFE2L2_vs_mRNA.json", # Number of sets differ.
        "CRISPR_SOX10_vs_mRNA.json", # Number of sets differ.
        "Cyclin_D1.json", # Number of sets differ.
        "EBV_Arrested.json", # Genes are duplicates.
        "ERbeta.json", # Number of sets differ.
        "European and African American Lung Comparison.json", # Enrichments differ.
        "Gender.json", # Number of sets differ.
        "MD_C1_vs_others.json", # Enrichments differ.
        "MYC mut and wt vs RNA.json", # Enrichments differ.
        "NRF2_Liver_Cancer.json", # Number of sets differ.
        "NRF2_mouse_model.json", # Number of sets differ.
        "PI3K_Inihibtion.json", # Enrichments differ.
        "PIK3CA mut and wt vs RNA.json", # Enrichments differ.
        "Regulation of ZBTB18 in glioblastoma.json", # Target directionalities differ.
        "Stroma_senescence.json", # Genes are duplicates.
    )

        continue

    end

    ke_va = BioLab.Dict.read(joinpath(DIJ, js))[chop(js; tail = 5)]

    @info "$id $js" ke_va

    dib = joinpath(DIB, BioLab.Path.clean(chop(js; tail = 5)))
    BioLab.Path.remake_directory(dib)

    dii = joinpath(dib, "input")
    BioLab.Path.remake_directory(dii)

    tst = joinpath(dii, "target_x_sample_x_number.tsv")

    tsf = joinpath(dii, "feature_x_sample_x_number.tsv")

    if !isfile(tst) || !isfile(tsf)

        @info "convert_cls_gct"

        GSEA.convert_cls_gct(tst, tsf, joinpath(DID, ke_va["cls"]), joinpath(DID, ke_va["ds"]))

    end

    jss = joinpath(dii, "set_features.json")

    if !isfile(jss)

        @info "convert_gmt"

        GSEA.convert_gmt(jss, (joinpath(DIS, gm) for gm in ke_va["gene_sets_collections"])...)

    end

    dir = joinpath(DIR, basename(ke_va["results_directory"]))

    for (al, pr) in zip(AL_, ke_va["results_files_prefix"])

        if al != "ks"

            continue

        end

        @info al

        dio = joinpath(dib, "output_$al")
        BioLab.Path.remake_directory(dio)

        txm = joinpath(dir, "$(pr)_gene_selection_scores.txt")

        txs = joinpath(dir, "$(pr)_GSEA_results_table.txt")

        tsm = joinpath(dio, "feature_x_metric_x_score.tsv")

        tss = joinpath(dio, "set_x_statistic_x_number.tsv")

        if !isfile(tsm)

            @info "metric-rank"

            GSEA.metric_rank(dio, tst, tsf, jss; algorithm = al, number_of_permutations = 0)

            BioLab.Path.remove(tss)

        end

        n_ch = 50

        py = sort!(
            BioLab.DataFrame.read(txm; select = [1, 2]);
            by = [an -> BioLab.String.limit(an, n_ch), an -> an],
        )

        ju =
            sort!(BioLab.DataFrame.read(tsm); by = [an -> BioLab.String.limit(an, n_ch), an -> an])

        @test size(py, 1) === size(ju, 1)

        for id in 1:size(py, 1)

            @test BioLab.String.limit(py[id, 1], n_ch) == BioLab.String.limit(ju[id, 1], n_ch)

            @test isapprox(py[id, 2], ju[id, 2]; atol = 1e-5)

        end

        if !isfile(tss)

            @info "user_rank"

            GSEA.user_rank(
                dio,
                txm,
                jss;
                algorithm = al,
                permutation = joinpath(dir, "$(pr)_rand_perm_gene_scores.txt"),
            )

        end

        py = sort!(BioLab.DataFrame.read(txs; select = [1, 5, 4, 6, 7]))

        ju = sort!(BioLab.DataFrame.read(tss))

        @test size(py, 1) === size(ju, 1)

        for id in 1:size(py, 1)

            @test py[id, 1] == ju[id, 1]

            @test isapprox(py[id, 3], ju[id, 2]; atol = 1e-2)

            @test isapprox(py[id, 2], ju[id, 3]; atol = 1e-2)

            @test isapprox(parse_float(py[id, 4]), ju[id, 4]; atol = 1e-2)

            @test isapprox(parse_float(py[id, 5]), ju[id, 5]; atol = 1e-1)

        end

    end

end

# ---- #

first(py, 10)

first(ju, 10)

# ---- #

is_ = .!isapprox.(py[!, 2], ju[!, 2]; atol = 1e-5);

view(py, is_, :)

view(ju, is_, :)

# ---- #

ta = convert(BitVector, Vector(BioLab.DataFrame.read(tst)[1, 2:end]));

fe = BioLab.DataFrame.read(tsf);

# ---- #

for ge in view(py, is_, 1)

    @info ge

    fe2 = Vector(fe[findfirst(==(ge), fe[!, 1]), 2:end])

    println(GSEA._get_signal_to_noise_ratio(fe2[.!ta], fe2[ta]))

end
