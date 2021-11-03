"""
Run standard GSEA
"""
@cast function run_standard_gsea(
    ke_ar::String,
    se_fe_::String,
    sa_va::String,
    ou::String,
)::Any

    ke_ar = DictExtension.read(ke_ar)

    se_fe_ = read_set(se_fe_, ke_ar)

    sa_va = table_read(sa_va)

    return ke_ar, se_fe_, sa_va

end

export run_standard_gsea
