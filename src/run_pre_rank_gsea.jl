"""
Run pre-rank GSEA
"""
@cast function run_pre_pank_gsea(
    ke_ar::String,
    se_fe_::String,
    fe_va::String,
    ou::String,
)::Any

    ke_ar = DictExtension.read(ke_ar)

    se_fe_ = read_set(se_fe_, ke_ar)

    fe_sc = table_read(fe_va)

    return ke_ar, se_fe_, fe_sc

end

export run_pre_pank_gsea
