using Test: @test

using BioLab

using GSEA

# ---- #

const DII = joinpath(homedir(), "Desktop", "benchmark")

const DIJ = joinpath(DII, "json_files")

const DIT = joinpath(DII, "Datasets_and_Phenotypes")

const DIS = joinpath(DII, "Gene_Sets_Collections")

const DIR = joinpath(DII, "results")

# ---- #

function parse_float(fl)

    if fl isa AbstractString

        fl = parse(Float64, lstrip(fl, '≤'))

    end

    fl

end

# ---- #

function path(di, re)

    if re || !isdir(di)

        BioLab.Path.remake_directory(di)

    end

    di

end

# ---- #

const RE = false

const DIB = path(joinpath(dirname(@__DIR__), "benchmark"), RE)

const AL_ = ("ks", "kli", "kli", "kliom", "kliop")

# ---- #

tst = nothing

tsf = nothing

pda = nothing

jda = nothing

# ---- #

for (id, js) in enumerate(BioLab.Path.read(DIJ))

    @info "$id $js"

    if js in (
        "CCLE_STAT3_vs_mRNA.json", # Number of sets.
        "CCLE_YAP_vs_mRNA.json", # Number of sets.
        "CFC1-overexpressing NGP cells.json", # Direction.
        "CRISPR_FOXA1_vs_mRNA.json", # Number of sets.
        "CRISPR_NFE2L2_vs_mRNA.json", # Number of sets.
        "CRISPR_SOX10_vs_mRNA.json", # Number of sets.
        "Cyclin_D1.json", # Number of sets.
        "EBV_Arrested.json", # Unique genes.
        "ERbeta.json", # Number of sets.
        "European and African American Lung Comparison.json", # Enrichment.
        "Gender.json", # Number of sets.
        "MD_C1_vs_others.json", # Enrichment.
        "MYC mut and wt vs RNA.json", # Enrichment.
        "NRF2_Liver_Cancer.json", # Number of sets.
        "NRF2_mouse_model.json", # Number of sets.
        "PI3K_Inihibtion.json", # Enrichment.
        "PIK3CA mut and wt vs RNA.json", # Enrichment.
        "Regulation of ZBTB18 in glioblastoma.json", # Direction.
        "Stroma_senescence.json", # Unique genes.
    )

        continue

    end

    ke_va = BioLab.Dict.read(joinpath(DIJ, js))[chop(js; tail = 5)]

    dib = path(joinpath(DIB, BioLab.Path.clean(chop(js; tail = 5))), RE)

    dii = path(joinpath(dib, "input"), RE)

    tst = joinpath(dii, "target_x_sample_x_number.tsv")

    tsf = joinpath(dii, "feature_x_sample_x_number.tsv")

    jss = joinpath(dii, "set_features.json")

    if RE || !isfile(tst) || !isfile(tsf)

        @info "Converting `cls` and `gct`"

        GSEA.convert_cls_gct(tst, tsf, joinpath(DIT, ke_va["cls"]), joinpath(DIT, ke_va["ds"]))

    end

    if RE || !isfile(jss)

        @info "Converting `gmt`"

        GSEA.convert_gmt(jss, (joinpath(DIS, gm) for gm in ke_va["gene_sets_collections"])...)

    end

    dip = joinpath(DIR, basename(ke_va["results_directory"]))

    for (al, pr) in zip(AL_, ke_va["results_files_prefix"])

        if al != "ks"

            continue

        end

        @info "Comparing \"$al\""

        dio = path(joinpath(dib, "output_$al"), RE)

        feature_x_metric_x_score_tsv = joinpath(dio, "feature_x_metric_x_score.tsv")

        if RE || !isfile(feature_x_metric_x_score_tsv)

            @info "Computing metrics"

            GSEA.metric_rank(dio, tst, tsf, jss; algorithm = al, number_of_permutations = 0)

        end

        feature_x_metric_x_score = joinpath(dip, "$(pr)_gene_selection_scores.txt")

        pda = BioLab.DataFrame.read(feature_x_metric_x_score; select = [1, 2])

        jda = BioLab.DataFrame.read(feature_x_metric_x_score_tsv)

        @test size(pda, 1) === size(jda, 1)

        pda[!, 1] = [BioLab.String.limit(st, 50) for st in pda[!, 1]]

        jda[!, 1] = [BioLab.String.limit(st, 50) for st in jda[!, 1]]

        pda = sort!(pda)

        jda = sort!(jda)

        for id in 1:size(pda, 1)

            @test pda[id, 1] == jda[id, 1]

            @test isapprox(pda[id, 2], jda[id, 2]; atol = 1e-5)

        end

        set_x_statistic_x_number_tsv = joinpath(dio, "set_x_statistic_x_number.tsv")

        if RE || !isfile(set_x_statistic_x_number_tsv)

            @info "Computing enrichments"

            GSEA.user_rank(
                dio,
                feature_x_metric_x_score,
                jss;
                algorithm = al,
                permutation = joinpath(dip, "$(pr)_rand_perm_gene_scores.txt"),
            )

        end

        pda = sort!(
            BioLab.DataFrame.read(
                joinpath(dip, "$(pr)_GSEA_results_table.txt");
                select = [1, 5, 4, 6, 7],
            ),
        )

        jda = sort!(BioLab.DataFrame.read(set_x_statistic_x_number_tsv))

        @test size(pda, 1) === size(jda, 1)

        for id in 1:size(pda, 1)

            @test pda[id, 1] == jda[id, 1]

            @test isapprox(pda[id, 3], jda[id, 2]; atol = 1e-2)

            #@test isapprox(pda[id, 2], jda[id, 3]; atol = 1e-2)

            #@test isapprox(parse_float(pda[id, 4]), jda[id, 4]; atol = 1e-2)

            #@test isapprox(parse_float(pda[id, 5]), jda[id, 5]; atol = 1e-1)

        end

    end

end

# ---- #
first(pda, 10)
first(jda, 10)
# ---- #
is_ = .!isapprox.(pda[!, 2], jda[!, 2]; atol = 1e-5);
view(pda, is_, :)
view(jda, is_, :)
# ---- #
ta = convert(BitVector, Vector(BioLab.DataFrame.read(tst)[1, 2:end]));
fe = BioLab.DataFrame.read(tsf);
# ---- #
for ge in view(pda, is_, 1)
    @info ge
    fe2 = Vector(fe[findfirst(==(ge), fe[!, 1]), 2:end])
    println(GSEA._get_signal_to_noise_ratio(fe2[.!ta], fe2[ta]))
end
