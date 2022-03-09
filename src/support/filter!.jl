function filter!(se_fe_, re, in_, mi, ma)

    println("Before filtering set-to-genes")

    OnePiece.dict.summarize(se_fe_, n_pr = 0)

    if re

        println("Removing set genes not found in gene-x-sample genes")

        for (se, fe_) in se_fe_

            se_fe_[se] = intersect(fe_, in_)

        end

    end

    for (se, fe_) in se_fe_

        if !(mi <= length(fe_) <= ma)

            pop!(se_fe_, se)

        end

    end

    println("After")

    OnePiece.dict.summarize(se_fe_, n_pr = 0)

end
