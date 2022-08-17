function _include_neighbor(fi)

    di, fi = splitdir(fi)

    for na in readdir(di)

        if (na != fi) && endswith(na, ".jl")

            include(na)

        end

    end

end
