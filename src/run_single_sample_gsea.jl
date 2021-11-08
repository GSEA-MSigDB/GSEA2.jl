"""
Run single-sample GSEA
"""
@cast function run_single_sample_gsea(
    js::String,
    gm::String,
    ts::String,
    ou::String,
)::DataFrame

    ke_ar = DictExtension.read(js)

    se_fe_ = read_set(gm, ke_ar)

    sc_fe_sa = TableAccess.read(ts)

    en_se_sa =
        score_set(sc_fe_sa, se_fe_; convert_to_keyword_argument(ke_ar)...)

    TableAccess.write(ou, en_se_sa)

    return en_se_sa

end

export run_single_sample_gsea
