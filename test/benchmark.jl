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
#BioLab.Path.remake_directory(DIB)

# ---- #

function test(st, is_, py, ju)

    if any(is_)

        @error st view(py, is_, :) view(ju, is_, :)

    end

end

function test(st, py, pyi, ju, jui)

    test(st, .!isequal.(py[!, pyi], ju[!, jui]), py, ju)

end

function test(st, py, pyi, ju, jui, atol)

    test(st, .!isapprox.(py[!, pyi], ju[!, jui]; atol), py, ju)

end

# ---- #

for (id, js) in enumerate(BioLab.Path.read(DIJ)[34:34])

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

    @info "$id $js"

    dib = joinpath(DIB, BioLab.Path.clean(chop(js; tail = 5)))
    #BioLab.Path.remake_directory(dib)

    dii = joinpath(dib, "input")
    #BioLab.Path.remake_directory(dii)

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

    for (al, pr) in zip(AL_, ke_va["results_files_prefix"])

        if al != "ks"

            continue

        end

        @info al

        dio = joinpath(dib, "output_$al")
        #BioLab.Path.remake_directory(dio)

        txm = joinpath(dir, "$(pr)_gene_selection_scores.txt")

        txs = joinpath(dir, "$(pr)_GSEA_results_table.txt")

        tsm = joinpath(dio, "feature_x_metric_x_score.tsv")

        tss = joinpath(dio, "set_x_statistic_x_number.tsv")

        if !isfile(tsm)

            GSEA.metric_rank(dio, tst, tsf, jss; algorithm = al, number_of_permutations = 0)

            BioLab.Path.remove(tss)

        end

        n_ch = 50

        py = BioLab.DataFrame.read(txm; select = [1, 2])

        ju = BioLab.DataFrame.read(tsm)

        py[!, 1] = BioLab.String.limit.(py[!, 1], n_ch)

        ju[!, 1] = BioLab.String.limit.(ju[!, 1], n_ch)

        sort!(py)

        sort!(ju)

        @test size(py, 1) === size(ju, 1)

        test("Gene", py, 1, ju, 1)

        test("Signal-to-Noise Ratio", py, 2, ju, 2, 1e-5)

        if !isfile(tss)

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

        test("Set", py, 1, ju, 1)

        test("Enrichment", py, 3, ju, 2, 1e-2)

        test("Normalized Enrichment", py, 2, ju, 3, 1e-2)

        test("P-Value", py, 4, ju, 4, 1e-2)

        test("Adjusted P-Value", py, 5, ju, 5, 1e-2)

    end

end
