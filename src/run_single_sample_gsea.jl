"""
Run single-sample GSEA
"""
@cast function run_single_sample_gsea(
    ke_ar::String,
    se_fe_::String,
    sc_fe_sa::String,
    ou::String,
)::Nothing

    ke_ar = DictExtension.read(ke_ar)

    se_fe_ = read_set(se_fe_, ke_ar)

    sc_fe_sa = table_read(sc_fe_sa)

    en_se_sa =
        score_set(sc_fe_sa, se_fe_; convert_to_keyword_argument(ke_ar)...)

    write(ou, en_se_sa)

    return

end

export run_single_sample_gsea
