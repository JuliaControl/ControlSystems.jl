# Get show compatible function with var set
printpolyfun(var) = (io, p, mimetype = MIME"text/plain"()) -> printpolydesc(io, p, var, mimetype)

# Specialized function, copied and adapted from Polynomials.jl
# Ignores var as set in Poly
""" Prints polynomial in descending order, with variable `var`
"""
function printpolydesc(io::IO, p::Poly{T}, var, mimetype = MIME"text/plain"()) where {T}
    printpoly2(io, p, var, mimetype)
end


# Specialized function, copied and adapted from Polynomials.jl
# Ignores var as set in Poly
function printpoly2(io::IO, p::Poly{T}, var, mimetype) where {T}
    first = true
    printed_anything = false
    for i in reverse(eachindex(p))
        printed = showterm2(io, p, i, first, var, mimetype)
        first &= !printed
        printed_anything |= printed
    end
    printed_anything || print(io, zero(T))
    return nothing
end

#Specialized function, copied and adapted from Polynomials.jl
# Ignores var as set in Poly
function showterm2(io::IO, p::Poly{T}, j, first, var, mimetype) where {T}
    pj = p[j]

    pj == zero(T) && return false

    pj = Polynomials.printsign(io, pj, first, mimetype)
    Polynomials.printcoefficient(io, pj, j, mimetype)
    Polynomials.printproductsign(io, pj, j, mimetype)
    Polynomials.printexponent(io, var, j, mimetype)
    true
end
