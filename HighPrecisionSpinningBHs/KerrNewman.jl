#KERR SOLUTION
function fKerrN(x::Double64,y::Double64,rh::Double64,χ::Double64,q::Double64=0.0)
    M = 2*rh/sqrt(1-χ^2-q^2)
    Q=q*M

    return (rh^2*(-3 + x)^4*(6*M^2*(-1 + x)^2 - 4*M*rh*(-1 + x)*(5 + (-2 + x)*x) + rh^2*(17 + (-2 + x)*x*(2 + (-2 + x)*x)) + 2*(M^2 - 4*rh^2)*(-1 + x)^2*cos(2*y)))/((8*M^2*(-1 + x)^2 + rh^2*(-3 + x)^2*(1 + x)^2 - 4*M*rh*(-1 + x)*(5 + (-2 + x)*x))^2 - 4*rh^2*(M^2 - 4*rh^2)*(-3 + x)^2*(-1 + x)^2*(1 + x)^2*sin(y)^2)
end
function gKerrN(x::Double64,y::Double64,rh::Double64,χ::Double64,q::Double64=0.0)
    return 1.0
end
function hKerrN(x::Double64,y::Double64,rh::Double64,χ::Double64,q::Double64=0.0)
    M = 2*rh/sqrt(1-χ^2-q^2)
    Q=q*M

    return (6*M^2*(-1 + x)^2 - 4*M*rh*(-1 + x)*(5 + (-2 + x)*x) + rh^2*(17 + (-2 + x)*x*(2 + (-2 + x)*x)) + 2*(M^2 - 4*rh^2)*(-1 + x)^2*cos(2*y))^2/((8*M^2*(-1 + x)^2 + rh^2*(-3 + x)^2*(1 + x)^2 - 4*M*rh*(-1 + x)*(5 + (-2 + x)*x))^2 - 4*rh^2*(M^2 - 4*rh^2)*(-3 + x)^2*(-1 + x)^2*(1 + x)^2*sin(y)^2)
end
function WKerrN(x::Double64,y::Double64,rh::Double64,χ::Double64,q::Double64=0.0)
    M = 2*rh/sqrt(1-χ^2-q^2)
    Q=q*M

    return (64*M*rh^2*sqrt(M^2 - 4*rh^2)*(1 - x)*(-2*M*(-1 + x) + rh*(5 - 2*x + x^2)))/((8*M^2*(-1 + x)^2 + rh^2*(-3 - 2*x + x^2)^2 - 4*M*rh*(-5 + 7*x - 3*x^2 + x^3))^2 - 4*rh^2*(M^2 - 4*rh^2)*(-1 + x)^2*(-3 - 2*x + x^2)^2*sin(y)^2)
end