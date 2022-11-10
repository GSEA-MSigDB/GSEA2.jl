macro _include()

    pa = string(__source__.file)

    :(_include($pa))

end

function _include(pa)

    di, fi = splitdir(pa)

    for na in readdir(di)

        if (na != fi) && endswith(na, ".jl")

            include(na)

        end

    end

end
