#!/usr/bin/env python3
"""Exact scalar and gluing checks for PIC pinching cones in dimensions 7 and 8.

Run: python3 verify_exact.py   (Python 3 and SymPy required)
No local modules, network access, floating-point tests, or positivity assumptions
on the trace variables are used. An unsuccessful check raises an exception.
This verifies algebra, not the geometric proof or imported theorem hypotheses.
"""
from itertools import combinations

import sympy as s

b = s.Symbol("b")
q = s.Rational
checks = 0


def check(name, condition):
    global checks
    if not bool(condition):
        raise AssertionError(name)
    checks += 1
    print("PASS:", name)


def identity(name, lhs, rhs):
    check(name, s.cancel(s.expand(lhs - rhs)) == 0)


def positive_interval(name, expr, end):
    """Exact Sturm root count and sign on (0,end], including denominators."""
    num, den = s.fraction(s.cancel(expr))
    for part in (num, den):
        p = s.Poly(part, b, domain=s.QQ)
        while p.eval(0) == 0:
            p = p.exquo(s.Poly(b, b))
        check(name + ": no roots", p.count_roots(0, end) == 0)
    check(name + ": endpoint sign", s.cancel(expr.subs(b, end)) > 0)


# Scalar checks for the unchanged Chen-type families. The n=7 checks concern
# separate invariance only; there is deliberately no n=7 gluing certificate.
for n, eps, end, second_end in ((8, q(1, 1000), q(1, 18), q(1, 40)),
                                 (7, q(1, 10000), q(131, 2000), q(3, 125))):
    tag = f"n={n} "
    a = (2 + (n-2)*b)**2*b / (2*(2 + (n-3)*b))
    gamma = b / (2 + (n-3)*b)
    rho = s.cancel(b - 2*(n-1)*gamma*(1-2*b)/n**2
          - 2*(n-1)*(1+gamma)*(n**2*b**2 - 2*(n-1)*(a-b)*(1-2*b))
          / (n**2*(1+2*(n-1)*a)))
    A = q(1, (n-1)*(n-4))*(2+8*b) + q(4, n)*(2*b+(n-2)*a)
    P = 2*(1+(n-2)*b)**2 / (1+2*(n-1)*a)
    Q = 2*((1+2*(n-1)*a)**2-(1+(n-2)*b)**2)/(n*(1+2*(n-1)*a))
    g = (1+2*(n-2)*a)**2/(1+2*(n-1)*a)*(2+(n-3)*b)/(1+(n-2)*b)
    h = (1+2*(n-1)*a)**2/((1+2*(n-2)*a)*(1+(n-2)*b)**2)
    R2 = (1-eps)**2 * 27*b*(2+(n-2)*b)/(32*n**2*rho**3)
    omega2 = (1-eps)**2 * q(27, 8)*(2+(n-2)*b)*b*(1+(n-2)*b)**2 \
             / (n**2*rho**3*(2+(n-3)*b)**2)
    identity(tag + "P+nQ", P+n*Q, 2+4*(n-1)*a)
    identity(tag + "omega/R factor", omega2/16,
             R2*((1+(n-2)*b)/(2*(2+(n-3)*b)))**2)
    for name, expr in (("g increasing", s.diff(g, b)),
                       ("h increasing", s.diff(h, b)), ("rho positive", rho),
                       ("rho<b", b-rho), ("rho'>4/9", s.diff(rho, b)-q(4, 9)),
                       ("rho concave", -s.diff(rho, b, 2)),
                       ("R decreasing", -s.diff(R2, b)), ("A positive", A)):
        positive_interval(tag + name, expr, end)
    c4num = n*b**2*(1-2*b)-2*(a-b)*(1-2*b+n*b**2)
    coefficient_poly = 2+(n-8)*b-2*(n+2)*(n-2)*b**2-n*(n-2)**2*b**3
    identity(tag + "remaining coefficient", c4num,
             b**2*coefficient_poly/(2+(n-3)*b))
    positive_interval(tag + "remaining coefficient positive", c4num, end)
    check(tag + "Chen-C comparison gate",
          (1+4*end)*(n-5)**2*(n-2) < (n-1)*(n-4)**2)
    check(tag + "Chen-C bound <1/2", (1+4*end)/(n-1) < q(1, 2))

    if n == 8:
        f, poly_q = 216*b**3+144*b**2+29*b+2, 42*b**2+21*b+2
        identity(tag + "g closed form", g, f**2/((1+6*b)**2*poly_q))
        identity(tag + "h closed form", h, poly_q**2/((2+5*b)*f))
        identity(tag + "rho closed form", rho,
                 3*b*(252*b**2+140*b+19)/(32*poly_q))
        identity(tag + "g' printed formula", s.diff(g, b),
                 f*(108864*b**5+117936*b**4+48672*b**3+9450*b**2+837*b+26)
                 / ((1+6*b)**3*poly_q**2))
        identity(tag + "h' printed formula", s.diff(h, b),
                 4*(1+4*b)*poly_q*(189*b**3+174*b**2+60*b+8)
                 / ((2+5*b)**2*f**2))
        check(tag + "rho endpoint", rho.subs(b, end) == q(31, 712))
        positive_interval(tag + "A>4b", A-4*b, end)
        check(tag + "sqrt48<7", 48 < 7**2)
        check(tag + "12-2sqrt6-6b>0", 12-2*3-6*end > 0)
        positive_interval(tag + "Lambda upper-bound gate",
                          2*(1+4*a)/28+4*a-A, end)
        bounds = (q(5, 4)**2*g, q(5, 4)*g*h)
        margins = (q(25045842231370311401, 34795857494016000000),
                   q(862611084879662129, 259607830464000000))
        for i, (bound, margin) in enumerate(zip(bounds, margins), 1):
            identity(tag + f"exact squared margin {i}",
                     (R2-bound**2).subs(b, end), margin)
            check(tag + f"margin {i} positive", margin > 0)
    else:
        check(tag + "sqrt35 lower enclosure", q(1479, 250)**2 < 35)
        check(tag + "sqrt35 upper enclosure", q(59161, 10000)**2 > 35)
        check(tag + "sqrt5 upper enclosure", q(22361, 10000)**2 > 5)
        positive_interval(tag + "first inequality enclosure",
                          1+q(1479, 500)*A-(1+q(22361, 10000)*b)**2, end)
        U = 1+q(59161, 20000)*A
        Y = 1+12*a
        S2 = U**2*(2+4*b)/(Y*(1+5*b))
        S3 = U*Y*(2+4*b)/(1+5*b)**3
        positive_interval(tag + "S2 increasing", s.diff(S2, b), end)
        positive_interval(tag + "S3 increasing", s.diff(S3, b), end)
        check(tag + "R2 endpoint lower bound", R2.subs(b, end) > q(21567, 1000))
        check(tag + "S2 squared endpoint upper bound", S2.subs(b, end)**2 < q(21529, 1000))
        check(tag + "S3 squared endpoint upper bound", S3.subs(b, end)**2 < q(18623, 1000))

    at = b+q(n-2, 2)*b**2
    zeta = (1+2*(n-1)*at)/(1+2*(n-1)*a.subs(b, end)) \
           * (1+(n-2)*end)/(1+(n-2)*b)*(1+gamma.subs(b, end))
    k = q(n**2-2*n+2, (n-2)**2)
    F = s.cancel(1+(n-2)*(1-zeta)-2*k*zeta**2*b/at)
    for name, expr in (("second c4 positive", n*b**2*(1-2*b)-2*(at-b)*(1-2*b+n*b**2)),
                       ("second Ricci coefficient positive", n**2*b**2-2*(n-1)*(at-b)*(1-2*b)),
                       ("zeta positive", zeta), ("zeta<1", 1-zeta),
                       ("zeta increasing", s.diff(zeta, b)), ("F positive", F),
                       ("F decreasing", -s.diff(F, b))):
        positive_interval(tag + name, expr, second_end)
    if n == 8:
        identity(tag + "F closed form", F,
                 (10780-250761*b-5352228*b**2-28171584*b**3-44881452*b**4)
                 / (7921*(1+3*b)*(1+6*b)**2))
        identity(tag + "F' printed formula", s.diff(F, b),
                 -63*(12823272*b**5+14960484*b**4+6686424*b**3+1453080*b**2+155270*b+6547)
                 / (7921*(1+3*b)**2*(1+6*b)**3))
        check(tag + "F endpoint", F.subs(b, second_end) == q(453196597, 7207159480))
        crude = 1+6*(1-zeta)-2*k*zeta**2
        check(tag + "discarded coupling fails",
              crude.subs(b, second_end) == q(-102165973, 1340866880))

# Endpoint shift and constants, derived from the maps' trace multipliers.
kap = (q(1, 18)-q(1, 40))/(1+6*q(1, 40))
DA = (q(49, 738)-q(43, 1600))/(1+14*q(43, 1600))
Y = DA-kap
R, V = 1+6*kap, 14*DA-6*kap
for name, actual, expected in (("kappa", kap, q(11, 414)),
                               ("Delta_A", DA, q(23333, 812538)),
                               ("Y", Y, q(20054, 9344187)),
                               ("Ricci multiplier", R, q(80, 69)),
                               ("scalar correction", V, q(2266960, 9344187))):
    check(name, actual == expected)
check("Y positive", Y > 0)
check("old scalar gluing estimate fails (squared comparison)",
      (q(2, 69))**2*86 < q(7, 22)**2)

# Free variables: the certificate identity requires no tensor positivity.
t = s.symbols("t12 t3 t4 tW t34 t3W t4W tWW")
t12, t3, t4, tW, t34, t3W, t4W, tWW = t
E = s.Matrix([[0,0,0,1,0,0,0,0], [0,0,0,0,0,0,0,1],
              [0,1,0,1,0,0,0,0], [0,0,0,1,0,0,1,0],
              [1,0,0,2,0,0,0,1], [1,1,0,1,0,1,0,0],
              [1,0,1,1,0,0,1,0], [1,1,1,0,1,0,0,0]])
print("E1,...,E8 coefficient rows in", t, ":", E.tolist())
d, c2, c3 = s.symbols("d c2 c3")
mu = s.Matrix([12*d+14*c2+52*c3-2, 8*c2+44*c3-1,
               2+4*d-6*c2-12*c3, 2*c2+12*c3-4*d,
               1-8*c2-20*c3, 8*c2+16*c3,
               4*d+4*c3-2*c2, 2*c2+4*c3])
U = 2*(t12+t3+t4+4*tW)
Ric33 = 2*t3+t34+4*t3W
scal = 2*t12+4*t3+4*t4+16*tW+2*t34+8*t3W+8*t4W+12*tWW
K2 = t12+2*t3+2*tW+2*d*U+2*c2*Ric33+2*c3*scal
for j, variable in enumerate(t):
    identity("2K certificate coefficient " + str(variable),
             s.expand(K2).coeff(variable), (E.T*mu)[j])

# Positivity in Q(sqrt(86)), certified using rational endpoint enclosures.
z = s.Symbol("sqrt86")
constants = {d: q(2, 69)*z-q(185, 828), c2: q(11, 207),
             c3: q(10027, 9344187)+q(28337, 18688374)*z}
lo, hi = q(927, 100), q(232, 25)
check("radical enclosure squared", lo**2 < 86 < hi**2)
for i, multiplier in enumerate(mu, 1):
    expr = s.expand(multiplier.subs(constants))
    lower = min(expr.subs(z, lo), expr.subs(z, hi))
    check(f"mu{i} rational lower bound {lower}", lower > 0)
    print(f"mu{i} =", expr.subs(z, s.sqrt(86)))

# All averaged frame identities use 28 independent sectional coordinates.
pairs = list(combinations(range(1, 9), 2))
sectionals = dict(zip(pairs, s.symbols("k0:28")))


def sec(i, j):
    return sectionals[tuple(sorted((i, j)))]


def frame(i, j, p, r):
    if len({i, j, p, r}) != 4:
        raise ValueError("not a four-frame")
    return sec(i, p)+sec(i, r)+sec(j, p)+sec(j, r)


def average(values):
    return sum(values)/len(values)


W = range(5, 9)
traces = [sec(1,2), (sec(1,3)+sec(2,3))/2, (sec(1,4)+sec(2,4))/2,
          sum(sec(i,p) for i in (1,2) for p in W)/8, sec(3,4),
          sum(sec(3,p) for p in W)/4, sum(sec(4,p) for p in W)/4,
          sum(sec(p,r) for p,r in combinations(W,2))/6]
averages = [average([frame(1,2,p,r) for p,r in combinations(W,2)])/4,
            (frame(5,6,7,8)+frame(5,7,6,8)+frame(5,8,6,7))/12,
            average([frame(1,2,3,p) for p in W])/2,
            average([frame(i,4,p,r) for i in (1,2) for p,r in combinations(W,2)])/2,
            average([frame(1,p,2,r) for p in W for r in W if p != r]),
            average([frame(i,3,3-i,p) for i in (1,2) for p in W]),
            average([frame(i,4,3-i,p) for i in (1,2) for p in W]),
            average([frame(i,3,3-i,4) for i in (1,2)]
                    + [frame(i,4,3-i,3) for i in (1,2)])]
for i in range(8):
    identity(f"E{i+1} exact sign-averaged frame identity",
             averages[i], (E*s.Matrix(traces))[i])

# Universal Bernstein and splitting identities, not tests at sampled lambda.
lam, beta, gam, x, tau, N0, N1 = s.symbols("lambda beta gamma x tau N0 N1")
Phi = beta+lam**2*gam-2*lam*x
L = Phi+(1-lam**2)*N0+lam**2*N1
quad = Phi/2+(1-lam**2)*(N0-tau/2)+lam**2*N1
identity("tilted-frame splitting", L, (Phi+(1-lam**2)*tau)/2+quad)
q0, q1 = quad.subs(lam, 0), quad.subs(lam, 1)
identity("quadratic Bernstein identity", quad,
         (1-lam)**2*q0+2*lam*(1-lam)*(q0-x/2)+lam**2*q1)
K = (beta-gam)/4-tau/2+N0
identity("middle control minus K", q0-x/2-K, (beta+gam-2*x)/4)
print(f"ALL {checks} EXACT CHECKS PASSED. No geometric proof claim is made.")
