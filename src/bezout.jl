"""
Calculates Bézout coefficients (see `gcdx`)
"""
function bezout{R}(a::R, b::R)
    rev = a < b
    x, y = rev ? (a,b) : (b,a)

    s0 = one(R)
    s1 = zero(R)
    t0 = zero(R)
    t1 = one(R)

    while x != zero(R)
        q = div(y, x)
        r = y - x * q
        s = s0 - q * s1
        t = t0 - q * t1

        s0 = s1
        s1 = s
        t0 = t1
        t1 = t
        y = x
        x = r
    end

    s, t = !rev ? (s0, t0) : (t0, s0)
    g = y

    if g == a
        s = one(R)
        t = zero(R)
    elseif g == -a
        s = -one(R)
        t = zero(R)
    end

    return s, t, g
end
