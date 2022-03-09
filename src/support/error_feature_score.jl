function error_feature_score(fe_, sc_)

    if length(fe_) != length(unique(fe_))

        error("Features are duplicated.")

    end

end

function error_feature_score(fe_x_sa_x_sc)

    error_feature_score(fe_x_sa_x_sc[!, 1], fe_x_sa_x_sc[!, 2])

end
