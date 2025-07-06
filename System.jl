#DO NOT FORGET TO INSTALL ALL THE NECESSARY JULIA LIBRARIES, BOTH THE ONES IN THE PACKAGE HighPrecisionSpinningBHs (CHECK THE FILE SpinningBlackHoles.jl) AND THOSE THAT SHOW UP IN THIS FILE, IF NOT INSTALLED YET

#THE METRIC FUNCTION W USED IN THE CODE DIFFERS FROM THE FUNCTION ω PRESENTED IN THE PAPER, BY A FACTOR OF r^2 FOR CONVENIENCE IN THE NUMERICAL METHOD. ω = W/r^2
using ChangePrecision, DoubleFloats

@changeprecision Double64 begin

#SMALL RESOLUTION TO PROVIDE AN EXAMPLE THAT RUNS QUICKLY
#DEFINE RESOLUTION IN x
global const Nx=60
#DEFINE RESOLUTION IN y
global const Ny=4

include("HighPrecisionSpinningBHs/SpinningBlackHoles.jl")
using .SpinningBlackHoles, NLsolve, DelimitedFiles, Cubature, Printf

LoadSystem()

global const NFields=6
#DEFINE BOUNDARY CONDITIONS AT THE HORIZON
function BC_Horizon!(i::Int64,j::Int64,f::Double64,dfdx::Double64,d2fdx2::Double64,dfdy::Double64,d2fdy2::Double64,d2fdxdy::Double64,g::Double64,dgdx::Double64,d2gdx2::Double64,dgdy::Double64,d2gdy2::Double64,d2gdxdy::Double64,h::Double64,dhdx::Double64,d2hdx2::Double64,dhdy::Double64,d2hdy2::Double64,d2hdxdy::Double64,W::Double64,dWdx::Double64,d2Wdx2::Double64,dWdy::Double64,d2Wdy2::Double64,d2Wdxdy::Double64,T::Double64,dTdx::Double64,d2Tdx2::Double64,dTdy::Double64,d2Tdy2::Double64,d2Tdxdy::Double64,U::Double64,dUdx::Double64,d2Udx2::Double64,dUdy::Double64,d2Udy2::Double64,d2Udxdy::Double64,R::Matrix{Double64},idx::Int64,WBC::Double64,rh::Double64,Mh::Double64,a0::Double64)
    R[idx+0] = dfdx
    R[idx+1] = dgdx
    R[idx+2] = dhdx
    R[idx+3] = -WBC + W
    R[idx+4] = dTdx
    R[idx+5] = dUdx
end

#DEFINE BOUNDARY CONDITIONS AT THE HORIZON - SPIN AS INPUT PARAMETER
function BC_Horizon_Spin!(i::Int64,j::Int64,f::Double64,dfdx::Double64,d2fdx2::Double64,dfdy::Double64,d2fdy2::Double64,d2fdxdy::Double64,g::Double64,dgdx::Double64,d2gdx2::Double64,dgdy::Double64,d2gdy2::Double64,d2gdxdy::Double64,h::Double64,dhdx::Double64,d2hdx2::Double64,dhdy::Double64,d2hdy2::Double64,d2hdxdy::Double64,W::Double64,dWdx::Double64,d2Wdx2::Double64,dWdy::Double64,d2Wdy2::Double64,d2Wdxdy::Double64,T::Double64,dTdx::Double64,d2Tdx2::Double64,dTdy::Double64,d2Tdy2::Double64,d2Tdxdy::Double64,U::Double64,dUdx::Double64,d2Udx2::Double64,dUdy::Double64,d2Udy2::Double64,d2Udxdy::Double64,R::Matrix{Double64},idx::Int64,WBC::Double64,rh::Double64,Mh::Double64,a0::Double64)
    R[idx+0] = dfdx
    R[idx+1] = dgdx
    R[idx+2] = dhdx
    R[idx+3] = W - dWdx
    R[idx+4] = dTdx
    R[idx+5] = dUdx
end

#DEFINE BOUNDARY CONDITIONS AT THE INFINITY
function BC_Infinity!(i::Int64,j::Int64,f::Double64,dfdx::Double64,d2fdx2::Double64,dfdy::Double64,d2fdy2::Double64,d2fdxdy::Double64,g::Double64,dgdx::Double64,d2gdx2::Double64,dgdy::Double64,d2gdy2::Double64,d2gdxdy::Double64,h::Double64,dhdx::Double64,d2hdx2::Double64,dhdy::Double64,d2hdy2::Double64,d2hdxdy::Double64,W::Double64,dWdx::Double64,d2Wdx2::Double64,dWdy::Double64,d2Wdy2::Double64,d2Wdxdy::Double64,T::Double64,dTdx::Double64,d2Tdx2::Double64,dTdy::Double64,d2Tdy2::Double64,d2Tdxdy::Double64,U::Double64,dUdx::Double64,d2Udx2::Double64,dUdy::Double64,d2Udy2::Double64,d2Udxdy::Double64,R::Matrix{Double64},idx::Int64,WBC::Double64,rh::Double64,Mh::Double64,a0::Double64)
    R[idx+0] = -1 + f
    R[idx+1] = -1 + g
    R[idx+2] = -1 + h
    R[idx+3] = W
    R[idx+4] = T
    R[idx+5] = U
end

#DEFINE BOUNDARY CONDITIONS AT THE INFINITY - SPIN AS INPUT PARAMETER
function BC_Infinity_Spin!(i::Int64,j::Int64,f::Double64,dfdx::Double64,d2fdx2::Double64,dfdy::Double64,d2fdy2::Double64,d2fdxdy::Double64,g::Double64,dgdx::Double64,d2gdx2::Double64,dgdy::Double64,d2gdy2::Double64,d2gdxdy::Double64,h::Double64,dhdx::Double64,d2hdx2::Double64,dhdy::Double64,d2hdy2::Double64,d2hdxdy::Double64,W::Double64,dWdx::Double64,d2Wdx2::Double64,dWdy::Double64,d2Wdy2::Double64,d2Wdxdy::Double64,T::Double64,dTdx::Double64,d2Tdx2::Double64,dTdy::Double64,d2Tdy2::Double64,d2Tdxdy::Double64,U::Double64,dUdx::Double64,d2Udx2::Double64,dUdy::Double64,d2Udy2::Double64,d2Udxdy::Double64,R::Matrix{Double64},idx::Int64,WBC::Double64,rh::Double64,Mh::Double64,a0::Double64)
    R[idx+0] = -1 + f
    R[idx+1] = -1 + g
    R[idx+2] = -1 + h
    R[idx+3] = rh*WBC*(2 + dfdx)^2 + dWdx
    R[idx+4] = T
    R[idx+5] = U
end

#DEFINE THE FIELD EQUATIONS
function Field_Eqs!(i::Int64,j::Int64,f::Double64,dfdx::Double64,d2fdx2::Double64,dfdy::Double64,d2fdy2::Double64,d2fdxdy::Double64,g::Double64,dgdx::Double64,d2gdx2::Double64,dgdy::Double64,d2gdy2::Double64,d2gdxdy::Double64,h::Double64,dhdx::Double64,d2hdx2::Double64,dhdy::Double64,d2hdy2::Double64,d2hdxdy::Double64,W::Double64,dWdx::Double64,d2Wdx2::Double64,dWdy::Double64,d2Wdy2::Double64,d2Wdxdy::Double64,T::Double64,dTdx::Double64,d2Tdx2::Double64,dTdy::Double64,d2Tdy2::Double64,d2Tdxdy::Double64,U::Double64,dUdx::Double64,d2Udx2::Double64,dUdy::Double64,d2Udy2::Double64,d2Udxdy::Double64,haxis::Double64,siny::Double64,cosy::Double64,R::Matrix{Double64},idx::Int64,WBC::Double64,rh::Double64,Mh::Double64,a0::Double64,x::Double64,y::Double64)
    R[idx+0] = g*h*(e(x,y,rh,Mh,a0) - T + U)*dfdy + f*(h*(2*T - U)*dgdy + g*(T*dhdy + 2*h*(cosy/siny*T - cosy/siny*U + dTdy)))

    R[idx+1] = (-1 + x)*e(x,y,rh,Mh,a0)*g*h*(-8*f + (-3 - 2*x + x^2)*dfdx) - (1 + x)*(-((3 - 4*x + x^2)*g*h*(T + U)*dfdx) + f*((3 - 4*x + x^2)*h*(T + U)*dgdx + g*(2*(1 + x)*h*(T + U) + (3 - 4*x + x^2)*T*dhdx)))

    R[idx+2] = 128*pi*rh^4*(-3 + x)^5*(1 + x)^2*e(x,y,rh,Mh,a0)*f*g^2*h + 256*rh^2*(-1 + x)^4*(1 + x)*f^2*dgdx + (-3 + x)*(-1 + x)^2*g*(64*rh^2*(1 + x)^2*(dfdy^2 + (-1 + x)^2*dfdx^2) + (-3 + x)^6*(-1 + x)^2*g*siny^2*(4*W^2 + dWdy^2 + 4*(-1 + x)*W*dWdx + (-1 + x)^2*dWdx^2)) + 32*rh^2*(1 + x)*f*(4*pi*rh^2*(-3 + x)^5*(1 + x)*g^2*h*(T + U) - (-1 + x)^2*(-3 - 2*x + x^2)*(dfdy*dgdy + (-1 + x)^2*dfdx*dgdx) - 2*(-1 + x)^2*g*((-3 - 2*x + x^2)*cosy/siny*dfdy + (-3 - 2*x + x^2)*d2fdy2 + (-1 + x)^2*(2*(-1 + x)*dfdx + (-3 - 2*x + x^2)*d2fdx2)))

    R[idx+3] = 4*(-3 - 2*x + x^2)*g*siny*(dfdy*dWdy + (-1 + x)*dfdx*(2*W + (-1 + x)*dWdx)) - f*(3*(-3 - 2*x + x^2)*siny*(dgdy*dWdy + (-1 + x)*dgdx*(2*W + (-1 + x)*dWdx)) + 2*g*(2*(-7 + 6*x + 5*x^2)*siny*W + 3*(-3 - 2*x + x^2)*cosy*dWdy + siny*((-3 - 2*x + x^2)*d2Wdy2 + (-1 + x)*(8*(-2 + x^2)*dWdx + (3 - x - 3*x^2 + x^3)*d2Wdx2))))

    R[idx+4] = -8*pi*rh^2*(-3 + x)^5*(1 + x)*g^3*h*T - (-1 + x)^2*(-3 - 2*x + x^2)*f*(dgdy^2 + (-1 + x)^2*dgdx^2) + 2*(-1 + x)^2*f*g*(2*(-3 - 2*x + x^2)*cosy/siny*dgdy + (-3 - 2*x + x^2)*d2gdy2 + (-1 + x)*((7 - 6*x + 3*x^2)*dgdx + (3 - x - 3*x^2 + x^3)*d2gdx2))

    if j==2
        R[idx+5] = haxis-1
    else
        R[idx+5] = 512*rh^2*(1 + x)*f*g^2*h^2*(pi*rh^2*(-3 + x)^5*(1 + x)*g*h*(T - U) - 2*(-1 + x)^4*dfdx) - (-3 + x)*(-1 + x)^2*g^2*h^2*(-64*rh^2*(1 + x)^2*(dfdy^2 + (-1 + x)^2*dfdx^2) + 3*(-3 + x)^6*(-1 + x)^2*g*siny^2*(4*W^2 + dWdy^2 + 4*(-1 + x)*W*dWdx + (-1 + x)^2*dWdx^2)) + 64*rh^2*(-1 + x)^2*(1 + x)*f^2*(-4*g*h^2*((-3 - 2*x + x^2)*cosy/siny*dgdy + (-5 + 7*x - 3*x^2 + x^3)*dgdx) - (-3 - 2*x + x^2)*h^2*(dgdy^2 + (-1 + x)^2*dgdx^2) - 2*(-3 - 2*x + x^2)*g^2*(dhdy^2 + (-1 + x)^2*dhdx^2 - h*(d2hdy2 + (-1 + x)*(dhdx + (-1 + x)*d2hdx2))))
    end
end

#DEFINE BOUNDARY CONDITIONS AT THE HORIZON
function BC_Horizon_Jac!(i::Int64,j::Int64,f::Double64,dfdx::Double64,d2fdx2::Double64,dfdy::Double64,d2fdy2::Double64,d2fdxdy::Double64,g::Double64,dgdx::Double64,d2gdx2::Double64,dgdy::Double64,d2gdy2::Double64,d2gdxdy::Double64,h::Double64,dhdx::Double64,d2hdx2::Double64,dhdy::Double64,d2hdy2::Double64,d2hdxdy::Double64,W::Double64,dWdx::Double64,d2Wdx2::Double64,dWdy::Double64,d2Wdy2::Double64,d2Wdxdy::Double64,T::Double64,dTdx::Double64,d2Tdx2::Double64,dTdy::Double64,d2Tdy2::Double64,d2Tdxdy::Double64,U::Double64,dUdx::Double64,d2Udx2::Double64,dUdy::Double64,d2Udy2::Double64,d2Udxdy::Double64,J::Matrix{Double64},idx::Int64,WBC::Double64,rh::Double64,Mh::Double64,a0::Double64,v::Int64,funcidx::Int64,type::Int8,ordx::Int64,ordy::Int64,Mx::Array{Double64, 3}=Mx,My::Array{Double64, 4}=My)
    if funcidx==1
        J[v, idx+0] = Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = 0
        J[v, idx+4] = 0
        J[v, idx+5] = 0
    end
    if funcidx==2
        J[v, idx+0] = 0
        J[v, idx+1] = Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+2] = 0
        J[v, idx+3] = 0
        J[v, idx+4] = 0
        J[v, idx+5] = 0
    end
    if funcidx==3
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+3] = 0
        J[v, idx+4] = 0
        J[v, idx+5] = 0
    end
    if funcidx==4
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+4] = 0
        J[v, idx+5] = 0
    end
    if funcidx==5
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = 0
        J[v, idx+4] = Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+5] = 0
    end
    if funcidx==6
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = 0
        J[v, idx+4] = 0
        J[v, idx+5] = Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]
    end
end

function BC_Horizon_Jac_Spin!(i::Int64,j::Int64,f::Double64,dfdx::Double64,d2fdx2::Double64,dfdy::Double64,d2fdy2::Double64,d2fdxdy::Double64,g::Double64,dgdx::Double64,d2gdx2::Double64,dgdy::Double64,d2gdy2::Double64,d2gdxdy::Double64,h::Double64,dhdx::Double64,d2hdx2::Double64,dhdy::Double64,d2hdy2::Double64,d2hdxdy::Double64,W::Double64,dWdx::Double64,d2Wdx2::Double64,dWdy::Double64,d2Wdy2::Double64,d2Wdxdy::Double64,T::Double64,dTdx::Double64,d2Tdx2::Double64,dTdy::Double64,d2Tdy2::Double64,d2Tdxdy::Double64,U::Double64,dUdx::Double64,d2Udx2::Double64,dUdy::Double64,d2Udy2::Double64,d2Udxdy::Double64,J::Matrix{Double64},idx::Int64,WBC::Double64,rh::Double64,Mh::Double64,a0::Double64,v::Int64,funcidx::Int64,type::Int8,ordx::Int64,ordy::Int64,Mx::Array{Double64, 3}=Mx,Cy::Array{Double64, 4}=My)
    if funcidx==1
        J[v, idx+0] = Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = 0
        J[v, idx+4] = 0
        J[v, idx+5] = 0
    end
    if funcidx==2
        J[v, idx+0] = 0
        J[v, idx+1] = Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+2] = 0
        J[v, idx+3] = 0
        J[v, idx+4] = 0
        J[v, idx+5] = 0
    end
    if funcidx==3
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+3] = 0
        J[v, idx+4] = 0
        J[v, idx+5] = 0
    end
    if funcidx==4
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = (Mx[1, 1 + ordx, i] - Mx[2, 1 + ordx, i])*My[1, 1 + ordy, j, type]
        J[v, idx+4] = 0
        J[v, idx+5] = 0
    end
    if funcidx==5
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = 0
        J[v, idx+4] = Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+5] = 0
    end
    if funcidx==6
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = 0
        J[v, idx+4] = 0
        J[v, idx+5] = Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]
    end
end

#DEFINE BOUNDARY CONDITIONS AT THE INFINITY
function BC_Infinity_Jac!(i::Int64,j::Int64,f::Double64,dfdx::Double64,d2fdx2::Double64,dfdy::Double64,d2fdy2::Double64,d2fdxdy::Double64,g::Double64,dgdx::Double64,d2gdx2::Double64,dgdy::Double64,d2gdy2::Double64,d2gdxdy::Double64,h::Double64,dhdx::Double64,d2hdx2::Double64,dhdy::Double64,d2hdy2::Double64,d2hdxdy::Double64,W::Double64,dWdx::Double64,d2Wdx2::Double64,dWdy::Double64,d2Wdy2::Double64,d2Wdxdy::Double64,T::Double64,dTdx::Double64,d2Tdx2::Double64,dTdy::Double64,d2Tdy2::Double64,d2Tdxdy::Double64,U::Double64,dUdx::Double64,d2Udx2::Double64,dUdy::Double64,d2Udy2::Double64,d2Udxdy::Double64,J::Matrix{Double64},idx::Int64,WBC::Double64,rh::Double64,Mh::Double64,a0::Double64,v::Int64,funcidx::Int64,type::Int8,ordx::Int64,ordy::Int64,Mx::Array{Double64, 3}=Mx,My::Array{Double64, 4}=My)
    if funcidx==1
        J[v, idx+0] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = 0
        J[v, idx+4] = 0
        J[v, idx+5] = 0
    end
    if funcidx==2
        J[v, idx+0] = 0
        J[v, idx+1] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+2] = 0
        J[v, idx+3] = 0
        J[v, idx+4] = 0
        J[v, idx+5] = 0
    end
    if funcidx==3
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+3] = 0
        J[v, idx+4] = 0
        J[v, idx+5] = 0
    end
    if funcidx==4
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+4] = 0
        J[v, idx+5] = 0
    end
    if funcidx==5
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = 0
        J[v, idx+4] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+5] = 0
    end
    if funcidx==6
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = 0
        J[v, idx+4] = 0
        J[v, idx+5] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]
    end
end

function BC_Infinity_Jac_Spin!(i::Int64,j::Int64,f::Double64,dfdx::Double64,d2fdx2::Double64,dfdy::Double64,d2fdy2::Double64,d2fdxdy::Double64,g::Double64,dgdx::Double64,d2gdx2::Double64,dgdy::Double64,d2gdy2::Double64,d2gdxdy::Double64,h::Double64,dhdx::Double64,d2hdx2::Double64,dhdy::Double64,d2hdy2::Double64,d2hdxdy::Double64,W::Double64,dWdx::Double64,d2Wdx2::Double64,dWdy::Double64,d2Wdy2::Double64,d2Wdxdy::Double64,T::Double64,dTdx::Double64,d2Tdx2::Double64,dTdy::Double64,d2Tdy2::Double64,d2Tdxdy::Double64,U::Double64,dUdx::Double64,d2Udx2::Double64,dUdy::Double64,d2Udy2::Double64,d2Udxdy::Double64,J::Matrix{Double64},idx::Int64,WBC::Double64,rh::Double64,Mh::Double64,a0::Double64,v::Int64,funcidx::Int64,type::Int8,ordx::Int64,ordy::Int64,Mx::Array{Double64, 3}=Mx,My::Array{Double64, 4}=My)
    if funcidx==1
        J[v, idx+0] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = 2*rh*WBC*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(2 + dfdx)
        J[v, idx+4] = 0
        J[v, idx+5] = 0
    end
    if funcidx==2
        J[v, idx+0] = 0
        J[v, idx+1] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+2] = 0
        J[v, idx+3] = 0
        J[v, idx+4] = 0
        J[v, idx+5] = 0
    end
    if funcidx==3
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+3] = 0
        J[v, idx+4] = 0
        J[v, idx+5] = 0
    end
    if funcidx==4
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+4] = 0
        J[v, idx+5] = 0
    end
    if funcidx==5
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = 0
        J[v, idx+4] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+5] = 0
    end
    if funcidx==6
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = 0
        J[v, idx+4] = 0
        J[v, idx+5] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]
    end
end

#DEFINE THE FIELD EQUATIONS
function Field_Eqs_Jac!(i::Int64,j::Int64,f::Double64,dfdx::Double64,d2fdx2::Double64,dfdy::Double64,d2fdy2::Double64,d2fdxdy::Double64,g::Double64,dgdx::Double64,d2gdx2::Double64,dgdy::Double64,d2gdy2::Double64,d2gdxdy::Double64,h::Double64,dhdx::Double64,d2hdx2::Double64,dhdy::Double64,d2hdy2::Double64,d2hdxdy::Double64,W::Double64,dWdx::Double64,d2Wdx2::Double64,dWdy::Double64,d2Wdy2::Double64,d2Wdxdy::Double64,T::Double64,dTdx::Double64,d2Tdx2::Double64,dTdy::Double64,d2Tdy2::Double64,d2Tdxdy::Double64,U::Double64,dUdx::Double64,d2Udx2::Double64,dUdy::Double64,d2Udy2::Double64,d2Udxdy::Double64,haxis::Double64,siny::Double64,cosy::Double64,J::Matrix{Double64},idx::Int64,WBC::Double64,rh::Double64,Mh::Double64,a0::Double64,v::Int64,funcidx::Int64,type::Int8,ordx::Int64,ordy::Int64,x::Double64,y::Double64,Mx::Array{Double64, 3}=Mx,My::Array{Double64, 4}=My)
    if funcidx==1
        J[v, idx+0] = g*h*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*(e(x,y,rh,Mh,a0) - T + U) + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(h*(2*T - U)*dgdy + g*(T*dhdy + 2*h*(cosy/siny*T - cosy/siny*U + dTdy)))

        J[v, idx+1] = Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*((-1 + x)*(-3 - 2*x + x^2)*e(x,y,rh,Mh,a0)*g*h + (1 + x)*(3 - 4*x + x^2)*g*h*(T + U)) + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-8*(-1 + x)*e(x,y,rh,Mh,a0)*g*h - (1 + x)*((3 - 4*x + x^2)*h*(T + U)*dgdx + g*(2*(1 + x)*h*(T + U) + (3 - 4*x + x^2)*T*dhdx)))

        J[v, idx+2] = -64*rh^2*(-1 + x)^4*(1 + x)*(-3 - 2*x + x^2)*f*g*Mx[3, 1 + ordx, i]*My[1, 1 + ordy, j, type] - 64*rh^2*(-1 + x)^2*(1 + x)*(-3 - 2*x + x^2)*f*g*Mx[1, 1 + ordx, i]*My[3, 1 + ordy, j, type] + Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*(128*rh^2*(-3 + x)*(-1 + x)^2*(1 + x)^2*g*dfdy + 32*rh^2*(1 + x)*f*(-2*(-1 + x)^2*(-3 - 2*x + x^2)*cosy/siny*g - (-1 + x)^2*(-3 - 2*x + x^2)*dgdy)) + Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(128*rh^2*(-3 + x)*(-1 + x)^4*(1 + x)^2*g*dfdx + 32*rh^2*(1 + x)*f*(-4*(-1 + x)^5*g - (-1 + x)^4*(-3 - 2*x + x^2)*dgdx)) + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(128*pi*rh^4*(-3 + x)^5*(1 + x)^2*e(x,y,rh,Mh,a0)*g^2*h + 512*rh^2*(-1 + x)^4*(1 + x)*f*dgdx + 32*rh^2*(1 + x)*(4*pi*rh^2*(-3 + x)^5*(1 + x)*g^2*h*(T + U) - (-1 + x)^2*(-3 - 2*x + x^2)*(dfdy*dgdy + (-1 + x)^2*dfdx*dgdx) - 2*(-1 + x)^2*g*((-3 - 2*x + x^2)*cosy/siny*dfdy + (-3 - 2*x + x^2)*d2fdy2 + (-1 + x)^2*(2*(-1 + x)*dfdx + (-3 - 2*x + x^2)*d2fdx2))))

        J[v, idx+3] = 4*(-3 - 2*x + x^2)*g*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*siny*dWdy + 4*(-1 + x)*(-3 - 2*x + x^2)*g*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*siny*(2*W + (-1 + x)*dWdx) + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-3*(-3 - 2*x + x^2)*siny*(dgdy*dWdy + (-1 + x)*dgdx*(2*W + (-1 + x)*dWdx)) - 2*g*(2*(-7 + 6*x + 5*x^2)*siny*W + 3*(-3 - 2*x + x^2)*cosy*dWdy + siny*((-3 - 2*x + x^2)*d2Wdy2 + (-1 + x)*(8*(-2 + x^2)*dWdx + (3 - x - 3*x^2 + x^3)*d2Wdx2))))

        J[v, idx+4] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-((-1 + x)^2*(-3 - 2*x + x^2)*(dgdy^2 + (-1 + x)^2*dgdx^2)) + 2*(-1 + x)^2*g*(2*(-3 - 2*x + x^2)*cosy/siny*dgdy + (-3 - 2*x + x^2)*d2gdy2 + (-1 + x)*((7 - 6*x + 3*x^2)*dgdx + (3 - x - 3*x^2 + x^3)*d2gdx2)))

        if j==2
            J[v, idx+5] = 0
        else
            J[v, idx+5] = 128*rh^2*(-3 + x)*(-1 + x)^2*(1 + x)^2*g^2*h^2*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*dfdy + Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-1024*rh^2*(-1 + x)^4*(1 + x)*f*g^2*h^2 + 128*rh^2*(-3 + x)*(-1 + x)^4*(1 + x)^2*g^2*h^2*dfdx) + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(512*rh^2*(1 + x)*g^2*h^2*(pi*rh^2*(-3 + x)^5*(1 + x)*g*h*(T - U) - 2*(-1 + x)^4*dfdx) + 128*rh^2*(-1 + x)^2*(1 + x)*f*(-4*g*h^2*((-3 - 2*x + x^2)*cosy/siny*dgdy + (-5 + 7*x - 3*x^2 + x^3)*dgdx) - (-3 - 2*x + x^2)*h^2*(dgdy^2 + (-1 + x)^2*dgdx^2) - 2*(-3 - 2*x + x^2)*g^2*(dhdy^2 + (-1 + x)^2*dhdx^2 - h*(d2hdy2 + (-1 + x)*(dhdx + (-1 + x)*d2hdx2)))))

        end
    end
    if funcidx==2
        J[v, idx+0] = f*h*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*(2*T - U) + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(h*(e(x,y,rh,Mh,a0) - T + U)*dfdy + f*(T*dhdy + 2*h*(cosy/siny*T - cosy/siny*U + dTdy)))

        J[v, idx+1] = -((1 + x)*(3 - 4*x + x^2)*f*h*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(T + U)) + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*((-1 + x)*e(x,y,rh,Mh,a0)*h*(-8*f + (-3 - 2*x + x^2)*dfdx) - (1 + x)*(-((3 - 4*x + x^2)*h*(T + U)*dfdx) + f*(2*(1 + x)*h*(T + U) + (3 - 4*x + x^2)*T*dhdx)))

        J[v, idx+2] = -32*rh^2*(-1 + x)^2*(1 + x)*(-3 - 2*x + x^2)*f*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*dfdy + Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(256*rh^2*(-1 + x)^4*(1 + x)*f^2 - 32*rh^2*(-1 + x)^4*(1 + x)*(-3 - 2*x + x^2)*f*dfdx) + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(256*pi*rh^4*(-3 + x)^5*(1 + x)^2*e(x,y,rh,Mh,a0)*f*g*h + (-3 + x)^7*(-1 + x)^4*g*siny^2*(4*W^2 + dWdy^2 + 4*(-1 + x)*W*dWdx + (-1 + x)^2*dWdx^2) + (-3 + x)*(-1 + x)^2*(64*rh^2*(1 + x)^2*(dfdy^2 + (-1 + x)^2*dfdx^2) + (-3 + x)^6*(-1 + x)^2*g*siny^2*(4*W^2 + dWdy^2 + 4*(-1 + x)*W*dWdx + (-1 + x)^2*dWdx^2)) + 32*rh^2*(1 + x)*f*(8*pi*rh^2*(-3 + x)^5*(1 + x)*g*h*(T + U) - 2*(-1 + x)^2*((-3 - 2*x + x^2)*cosy/siny*dfdy + (-3 - 2*x + x^2)*d2fdy2 + (-1 + x)^2*(2*(-1 + x)*dfdx + (-3 - 2*x + x^2)*d2fdx2))))

        J[v, idx+3] = -3*(-3 - 2*x + x^2)*f*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*siny*dWdy - 3*(-1 + x)*(-3 - 2*x + x^2)*f*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*siny*(2*W + (-1 + x)*dWdx) + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(4*(-3 - 2*x + x^2)*siny*(dfdy*dWdy + (-1 + x)*dfdx*(2*W + (-1 + x)*dWdx)) - 2*f*(2*(-7 + 6*x + 5*x^2)*siny*W + 3*(-3 - 2*x + x^2)*cosy*dWdy + siny*((-3 - 2*x + x^2)*d2Wdy2 + (-1 + x)*(8*(-2 + x^2)*dWdx + (3 - x - 3*x^2 + x^3)*d2Wdx2))))

        J[v, idx+4] = 2*(-1 + x)^3*(3 - x - 3*x^2 + x^3)*f*g*Mx[3, 1 + ordx, i]*My[1, 1 + ordy, j, type] + 2*(-1 + x)^2*(-3 - 2*x + x^2)*f*g*Mx[1, 1 + ordx, i]*My[3, 1 + ordy, j, type] + Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*(4*(-1 + x)^2*(-3 - 2*x + x^2)*cosy/siny*f*g - 2*(-1 + x)^2*(-3 - 2*x + x^2)*f*dgdy) + Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(2*(-1 + x)^3*(7 - 6*x + 3*x^2)*f*g - 2*(-1 + x)^4*(-3 - 2*x + x^2)*f*dgdx) + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-24*pi*rh^2*(-3 + x)^5*(1 + x)*g^2*h*T + 2*(-1 + x)^2*f*(2*(-3 - 2*x + x^2)*cosy/siny*dgdy + (-3 - 2*x + x^2)*d2gdy2 + (-1 + x)*((7 - 6*x + 3*x^2)*dgdx + (3 - x - 3*x^2 + x^3)*d2gdx2)))

        if j==2
            J[v, idx+5] = 0
        else
            J[v, idx+5] = 64*rh^2*(-1 + x)^2*(1 + x)*f^2*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*(-4*(-3 - 2*x + x^2)*cosy/siny*g*h^2 - 2*(-3 - 2*x + x^2)*h^2*dgdy) + 64*rh^2*(-1 + x)^2*(1 + x)*f^2*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-4*(-5 + 7*x - 3*x^2 + x^3)*g*h^2 - 2*(-1 + x)^2*(-3 - 2*x + x^2)*h^2*dgdx) + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(512*pi*rh^4*(-3 + x)^5*(1 + x)^2*f*g^2*h^3*(T - U) + 1024*rh^2*(1 + x)*f*g*h^2*(pi*rh^2*(-3 + x)^5*(1 + x)*g*h*(T - U) - 2*(-1 + x)^4*dfdx) - 3*(-3 + x)^7*(-1 + x)^4*g^2*h^2*siny^2*(4*W^2 + dWdy^2 + 4*(-1 + x)*W*dWdx + (-1 + x)^2*dWdx^2) - 2*(-3 + x)*(-1 + x)^2*g*h^2*(-64*rh^2*(1 + x)^2*(dfdy^2 + (-1 + x)^2*dfdx^2) + 3*(-3 + x)^6*(-1 + x)^2*g*siny^2*(4*W^2 + dWdy^2 + 4*(-1 + x)*W*dWdx + (-1 + x)^2*dWdx^2)) + 64*rh^2*(-1 + x)^2*(1 + x)*f^2*(-4*h^2*((-3 - 2*x + x^2)*cosy/siny*dgdy + (-5 + 7*x - 3*x^2 + x^3)*dgdx) - 4*(-3 - 2*x + x^2)*g*(dhdy^2 + (-1 + x)^2*dhdx^2 - h*(d2hdy2 + (-1 + x)*(dhdx + (-1 + x)*d2hdx2)))))

        end
    end
    if funcidx==3
        J[v, idx+0] = f*g*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*T + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(g*(e(x,y,rh,Mh,a0) - T + U)*dfdy + f*((2*T - U)*dgdy + 2*g*(cosy/siny*T - cosy/siny*U + dTdy)))

        J[v, idx+1] = -((1 + x)*(3 - 4*x + x^2)*f*g*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*T) + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*((-1 + x)*e(x,y,rh,Mh,a0)*g*(-8*f + (-3 - 2*x + x^2)*dfdx) - (1 + x)*(-((3 - 4*x + x^2)*g*(T + U)*dfdx) + f*(2*(1 + x)*g*(T + U) + (3 - 4*x + x^2)*(T + U)*dgdx)))

        J[v, idx+2] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(128*pi*rh^4*(-3 + x)^5*(1 + x)^2*e(x,y,rh,Mh,a0)*f*g^2 + 128*pi*rh^4*(-3 + x)^5*(1 + x)^2*f*g^2*(T + U))

        J[v, idx+3] = 0

        J[v, idx+4] = -8*pi*rh^2*(-3 + x)^5*(1 + x)*g^3*Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*T

        if j==2
            J[v, idx+5] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, 1, 1]
        else
            J[v, idx+5] = 128*rh^2*(-1 + x)^4*(1 + x)*(-3 - 2*x + x^2)*f^2*g^2*h*Mx[3, 1 + ordx, i]*My[1, 1 + ordy, j, type] + 128*rh^2*(-1 + x)^2*(1 + x)*(-3 - 2*x + x^2)*f^2*g^2*h*Mx[1, 1 + ordx, i]*My[3, 1 + ordy, j, type] - 256*rh^2*(-1 + x)^2*(1 + x)*(-3 - 2*x + x^2)*f^2*g^2*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*dhdy - 128*rh^2*(-1 + x)^2*(1 + x)*(-3 - 2*x + x^2)*f^2*g^2*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-((-1 + x)*h) + 2*(-1 + x)^2*dhdx) + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(512*pi*rh^4*(-3 + x)^5*(1 + x)^2*f*g^3*h^2*(T - U) + 1024*rh^2*(1 + x)*f*g^2*h*(pi*rh^2*(-3 + x)^5*(1 + x)*g*h*(T - U) - 2*(-1 + x)^4*dfdx) - 2*(-3 + x)*(-1 + x)^2*g^2*h*(-64*rh^2*(1 + x)^2*(dfdy^2 + (-1 + x)^2*dfdx^2) + 3*(-3 + x)^6*(-1 + x)^2*g*siny^2*(4*W^2 + dWdy^2 + 4*(-1 + x)*W*dWdx + (-1 + x)^2*dWdx^2)) + 64*rh^2*(-1 + x)^2*(1 + x)*f^2*(-8*g*h*((-3 - 2*x + x^2)*cosy/siny*dgdy + (-5 + 7*x - 3*x^2 + x^3)*dgdx) - 2*(-3 - 2*x + x^2)*h*(dgdy^2 + (-1 + x)^2*dgdx^2) + 2*(-3 - 2*x + x^2)*g^2*(d2hdy2 + (-1 + x)*(dhdx + (-1 + x)*d2hdx2))))

        end
    end
    if funcidx==4
        J[v, idx+0] = 0

        J[v, idx+1] = 0

        J[v, idx+2] = 2*(-3 + x)^7*(-1 + x)^4*g^2*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*siny^2*dWdy + (-3 + x)^7*(-1 + x)^4*g^2*Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*siny^2*(8*W + 4*(-1 + x)*dWdx) + (-3 + x)^7*(-1 + x)^4*g^2*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*siny^2*(4*(-1 + x)*W + 2*(-1 + x)^2*dWdx)

        J[v, idx+3] = -2*(-1 + x)*(3 - x - 3*x^2 + x^3)*f*g*Mx[3, 1 + ordx, i]*My[1, 1 + ordy, j, type]*siny - 2*(-3 - 2*x + x^2)*f*g*Mx[1, 1 + ordx, i]*My[3, 1 + ordy, j, type]*siny + Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*(4*(-3 - 2*x + x^2)*g*siny*dfdy - f*(6*(-3 - 2*x + x^2)*cosy*g + 3*(-3 - 2*x + x^2)*siny*dgdy)) + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(8*(-1 + x)*(-3 - 2*x + x^2)*g*siny*dfdx - f*(4*(-7 + 6*x + 5*x^2)*g*siny + 6*(-1 + x)*(-3 - 2*x + x^2)*siny*dgdx)) + Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(4*(-1 + x)^2*(-3 - 2*x + x^2)*g*siny*dfdx - f*(16*(-1 + x)*(-2 + x^2)*g*siny + 3*(-1 + x)^2*(-3 - 2*x + x^2)*siny*dgdx))

        J[v, idx+4] = 0

        if j==2
            J[v, idx+5] = 0
        else
            J[v, idx+5] = -6*(-3 + x)^7*(-1 + x)^4*g^3*h^2*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*siny^2*dWdy - 3*(-3 + x)^7*(-1 + x)^4*g^3*h^2*Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*siny^2*(8*W + 4*(-1 + x)*dWdx) - 3*(-3 + x)^7*(-1 + x)^4*g^3*h^2*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*siny^2*(4*(-1 + x)*W + 2*(-1 + x)^2*dWdx)

        end
    end
    if funcidx==5
        J[v, idx+0] = 2*f*g*h*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type] + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-(g*h*dfdy) + f*(2*h*dgdy + g*(2*cosy/siny*h + dhdy)))

        J[v, idx+1] = -((1 + x)*Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-((3 - 4*x + x^2)*g*h*dfdx) + f*((3 - 4*x + x^2)*h*dgdx + g*(2*(1 + x)*h + (3 - 4*x + x^2)*dhdx))))

        J[v, idx+2] = 128*pi*rh^4*(-3 + x)^5*(1 + x)^2*f*g^2*h*Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]

        J[v, idx+3] = 0

        J[v, idx+4] = -8*pi*rh^2*(-3 + x)^5*(1 + x)*g^3*h*Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]

        if j==2
            J[v, idx+5] = 0
        else
            J[v, idx+5] = 512*pi*rh^4*(-3 + x)^5*(1 + x)^2*f*g^3*h^3*Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]

        end
    end
    if funcidx==6
        J[v, idx+0] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(g*h*dfdy + f*(-2*cosy/siny*g*h - h*dgdy))

        J[v, idx+1] = -((1 + x)*Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-((3 - 4*x + x^2)*g*h*dfdx) + f*(2*(1 + x)*g*h + (3 - 4*x + x^2)*h*dgdx)))

        J[v, idx+2] = 128*pi*rh^4*(-3 + x)^5*(1 + x)^2*f*g^2*h*Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]

        J[v, idx+3] = 0

        J[v, idx+4] = 0

        if j==2
            J[v, idx+5] = 0
        else
            J[v, idx+5] = -512*pi*rh^4*(-3 + x)^5*(1 + x)^2*f*g^3*h^3*Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]

        end
    end
end
#DEFINE THE WHOLE SYSTEM (BCs+FIELD EQS) ON OUR GRID
function residual_jac!(R::Matrix{Double64},J::Matrix{Double64},a::Matrix{Double64},f::Field,g::Field,h::Field,W::Field,T::Field,U::Field,WBC::Double64,rh::Double64,Mh::Double64,a0::Double64,spin::Bool,x::Vector{Double64}=X,y::Vector{Double64}=Y)

    f.a=a[0*Nx+1:1*Nx,:]

    g.a=a[1*Nx+1:2*Nx,:]

    h.a=a[2*Nx+1:3*Nx,:]

    W.a=a[3*Nx+1:4*Nx,:]

    T.a=a[4*Nx+1:5*Nx,:]

    U.a=a[5*Nx+1:6*Nx,:]

    type=[f.type,g.type,h.type,W.type,T.type,U.type]
    ff=Matrix{Double64}(undef,Nx,Ny); dfdx=Matrix{Double64}(undef,Nx,Ny); d2fdx2=Matrix{Double64}(undef,Nx,Ny); dfdy=Matrix{Double64}(undef,Nx,Ny); d2fdy2=Matrix{Double64}(undef,Nx,Ny); d2fdxdy=Matrix{Double64}(undef,Nx,Ny);
    gg=Matrix{Double64}(undef,Nx,Ny); dgdx=Matrix{Double64}(undef,Nx,Ny); d2gdx2=Matrix{Double64}(undef,Nx,Ny); dgdy=Matrix{Double64}(undef,Nx,Ny); d2gdy2=Matrix{Double64}(undef,Nx,Ny); d2gdxdy=Matrix{Double64}(undef,Nx,Ny);
    hh=Matrix{Double64}(undef,Nx,Ny); dhdx=Matrix{Double64}(undef,Nx,Ny); d2hdx2=Matrix{Double64}(undef,Nx,Ny); dhdy=Matrix{Double64}(undef,Nx,Ny); d2hdy2=Matrix{Double64}(undef,Nx,Ny); d2hdxdy=Matrix{Double64}(undef,Nx,Ny);
    WW=Matrix{Double64}(undef,Nx,Ny); dWdx=Matrix{Double64}(undef,Nx,Ny); d2Wdx2=Matrix{Double64}(undef,Nx,Ny); dWdy=Matrix{Double64}(undef,Nx,Ny); d2Wdy2=Matrix{Double64}(undef,Nx,Ny); d2Wdxdy=Matrix{Double64}(undef,Nx,Ny);
    TT=Matrix{Double64}(undef,Nx,Ny); dTdx=Matrix{Double64}(undef,Nx,Ny); d2Tdx2=Matrix{Double64}(undef,Nx,Ny); dTdy=Matrix{Double64}(undef,Nx,Ny); d2Tdy2=Matrix{Double64}(undef,Nx,Ny); d2Tdxdy=Matrix{Double64}(undef,Nx,Ny);
    UU=Matrix{Double64}(undef,Nx,Ny); dUdx=Matrix{Double64}(undef,Nx,Ny); d2Udx2=Matrix{Double64}(undef,Nx,Ny); dUdy=Matrix{Double64}(undef,Nx,Ny); d2Udy2=Matrix{Double64}(undef,Nx,Ny); d2Udxdy=Matrix{Double64}(undef,Nx,Ny);
    haxis=Array{Double64}(undef,Nx);

    for i in 1:Nx
        for j in 1:Ny
            ff[i,j]=f(i,j+1)
            dfdx[i,j]=f(i,j+1,dx=1)
            d2fdx2[i,j]=f(i,j+1,dx=2)
            dfdy[i,j]=f(i,j+1,dy=1)
            d2fdy2[i,j]=f(i,j+1,dy=2)
            d2fdxdy[i,j]=f(i,j+1,dx=1,dy=1)
            gg[i,j]=g(i,j+1)
            dgdx[i,j]=g(i,j+1,dx=1)
            d2gdx2[i,j]=g(i,j+1,dx=2)
            dgdy[i,j]=g(i,j+1,dy=1)
            d2gdy2[i,j]=g(i,j+1,dy=2)
            d2gdxdy[i,j]=g(i,j+1,dx=1,dy=1)
            hh[i,j]=h(i,j+1)
            dhdx[i,j]=h(i,j+1,dx=1)
            d2hdx2[i,j]=h(i,j+1,dx=2)
            dhdy[i,j]=h(i,j+1,dy=1)
            d2hdy2[i,j]=h(i,j+1,dy=2)
            d2hdxdy[i,j]=h(i,j+1,dx=1,dy=1)
            WW[i,j]=W(i,j+1)
            dWdx[i,j]=W(i,j+1,dx=1)
            d2Wdx2[i,j]=W(i,j+1,dx=2)
            dWdy[i,j]=W(i,j+1,dy=1)
            d2Wdy2[i,j]=W(i,j+1,dy=2)
            d2Wdxdy[i,j]=W(i,j+1,dx=1,dy=1)
            TT[i,j]=T(i,j+1)
            dTdx[i,j]=T(i,j+1,dx=1)
            d2Tdx2[i,j]=T(i,j+1,dx=2)
            dTdy[i,j]=T(i,j+1,dy=1)
            d2Tdy2[i,j]=T(i,j+1,dy=2)
            d2Tdxdy[i,j]=T(i,j+1,dx=1,dy=1)
            UU[i,j]=U(i,j+1)
            dUdx[i,j]=U(i,j+1,dx=1)
            d2Udx2[i,j]=U(i,j+1,dx=2)
            dUdy[i,j]=U(i,j+1,dy=1)
            d2Udy2[i,j]=U(i,j+1,dy=2)
            d2Udxdy[i,j]=U(i,j+1,dx=1,dy=1)
        end
        haxis[i]=h(i,1)
    end

    siny=Array{Double64}(undef,Ny);
    cosy=Array{Double64}(undef,Ny);
    for j in 1:Ny
        siny[j]=sin(y[j+1])
        cosy[j]=cos(y[j+1])
    end

    if !(R == nothing)
        #INDEX FOR THE RESIDUAL R
        idx=1
        #LOOP ON ALL X POINTS INCLUDING -1 AND 1
        for i in 1:Nx
        #LOOP ON INTERNAL Y POINTS. NOTE THAT BC ON 0 AND PI/2 ARE AUTOMATICALLY IMPOSED BY OUR CHOICE OF BASIS FUNCTIONS (EVEN COSINES)
            for j in 1:Ny
                if i==1
                    #DEFINE BCS AT THE HORIZON
                    if !spin
                        BC_Horizon!(i,j,ff[i,j],dfdx[i,j],d2fdx2[i,j],dfdy[i,j],d2fdy2[i,j],d2fdxdy[i,j],gg[i,j],dgdx[i,j],d2gdx2[i,j],dgdy[i,j],d2gdy2[i,j],d2gdxdy[i,j],hh[i,j],dhdx[i,j],d2hdx2[i,j],dhdy[i,j],d2hdy2[i,j],d2hdxdy[i,j],WW[i,j],dWdx[i,j],d2Wdx2[i,j],dWdy[i,j],d2Wdy2[i,j],d2Wdxdy[i,j],TT[i,j],dTdx[i,j],d2Tdx2[i,j],dTdy[i,j],d2Tdy2[i,j],d2Tdxdy[i,j],UU[i,j],dUdx[i,j],d2Udx2[i,j],dUdy[i,j],d2Udy2[i,j],d2Udxdy[i,j],R,idx,WBC,rh,Mh,a0)
                    else
                        BC_Horizon_Spin!(i,j,ff[i,j],dfdx[i,j],d2fdx2[i,j],dfdy[i,j],d2fdy2[i,j],d2fdxdy[i,j],gg[i,j],dgdx[i,j],d2gdx2[i,j],dgdy[i,j],d2gdy2[i,j],d2gdxdy[i,j],hh[i,j],dhdx[i,j],d2hdx2[i,j],dhdy[i,j],d2hdy2[i,j],d2hdxdy[i,j],WW[i,j],dWdx[i,j],d2Wdx2[i,j],dWdy[i,j],d2Wdy2[i,j],d2Wdxdy[i,j],TT[i,j],dTdx[i,j],d2Tdx2[i,j],dTdy[i,j],d2Tdy2[i,j],d2Tdxdy[i,j],UU[i,j],dUdx[i,j],d2Udx2[i,j],dUdy[i,j],d2Udy2[i,j],d2Udxdy[i,j],R,idx,WBC,rh,Mh,a0)
                    end
                elseif i==Nx
                    #DEFINE BCS AT INFINITY
                    if !spin
                        BC_Infinity!(i,j,ff[i,j],dfdx[i,j],d2fdx2[i,j],dfdy[i,j],d2fdy2[i,j],d2fdxdy[i,j],gg[i,j],dgdx[i,j],d2gdx2[i,j],dgdy[i,j],d2gdy2[i,j],d2gdxdy[i,j],hh[i,j],dhdx[i,j],d2hdx2[i,j],dhdy[i,j],d2hdy2[i,j],d2hdxdy[i,j],WW[i,j],dWdx[i,j],d2Wdx2[i,j],dWdy[i,j],d2Wdy2[i,j],d2Wdxdy[i,j],TT[i,j],dTdx[i,j],d2Tdx2[i,j],dTdy[i,j],d2Tdy2[i,j],d2Tdxdy[i,j],UU[i,j],dUdx[i,j],d2Udx2[i,j],dUdy[i,j],d2Udy2[i,j],d2Udxdy[i,j],R,idx,WBC,rh,Mh,a0)
                    else
                        BC_Infinity_Spin!(i,j,ff[i,j],dfdx[i,j],d2fdx2[i,j],dfdy[i,j],d2fdy2[i,j],d2fdxdy[i,j],gg[i,j],dgdx[i,j],d2gdx2[i,j],dgdy[i,j],d2gdy2[i,j],d2gdxdy[i,j],hh[i,j],dhdx[i,j],d2hdx2[i,j],dhdy[i,j],d2hdy2[i,j],d2hdxdy[i,j],WW[i,j],dWdx[i,j],d2Wdx2[i,j],dWdy[i,j],d2Wdy2[i,j],d2Wdxdy[i,j],TT[i,j],dTdx[i,j],d2Tdx2[i,j],dTdy[i,j],d2Tdy2[i,j],d2Tdxdy[i,j],UU[i,j],dUdx[i,j],d2Udx2[i,j],dUdy[i,j],d2Udy2[i,j],d2Udxdy[i,j],R,idx,WBC,rh,Mh,a0)
                    end
                else
                    #DEFINE FIELD EQUATIONS EVERYWHERE ELSE ON THE GRID
                    Field_Eqs!(i,j+1,ff[i,j],dfdx[i,j],d2fdx2[i,j],dfdy[i,j],d2fdy2[i,j],d2fdxdy[i,j],gg[i,j],dgdx[i,j],d2gdx2[i,j],dgdy[i,j],d2gdy2[i,j],d2gdxdy[i,j],hh[i,j],dhdx[i,j],d2hdx2[i,j],dhdy[i,j],d2hdy2[i,j],d2hdxdy[i,j],WW[i,j],dWdx[i,j],d2Wdx2[i,j],dWdy[i,j],d2Wdy2[i,j],d2Wdxdy[i,j],TT[i,j],dTdx[i,j],d2Tdx2[i,j],dTdy[i,j],d2Tdy2[i,j],d2Tdxdy[i,j],UU[i,j],dUdx[i,j],d2Udx2[i,j],dUdy[i,j],d2Udy2[i,j],d2Udxdy[i,j],haxis[i],siny[j],cosy[j],R,idx,WBC,rh,Mh,a0,X[i],Y[j])
                end
                idx+=NFields
            end
        end
    end

    if !(J == nothing)
        funcidx=1
        ordx=0
        ordy=0
        for v in 1:Nx*Ny*NFields
            #INDEX FOR THE JACOBIAN J COLUMNS
            idx=1
            #LOOP ON ALL X POINTS INCLUDING -1 AND 1
            for i in 1:Nx
            #LOOP ON INTERNAL Y POINTS. NOTE THAT BC ON 0 AND PI/2 ARE AUTOMATICALLY IMPOSED BY OUR CHOICE OF BASIS FUNCTIONS (EVEN COSINES)
                for j in 1:Ny
                    if i==1
                        #DEFINE BCS AT THE HORIZON
                        if !spin
                            BC_Horizon_Jac!(i,j+1,ff[i,j],dfdx[i,j],d2fdx2[i,j],dfdy[i,j],d2fdy2[i,j],d2fdxdy[i,j],gg[i,j],dgdx[i,j],d2gdx2[i,j],dgdy[i,j],d2gdy2[i,j],d2gdxdy[i,j],hh[i,j],dhdx[i,j],d2hdx2[i,j],dhdy[i,j],d2hdy2[i,j],d2hdxdy[i,j],WW[i,j],dWdx[i,j],d2Wdx2[i,j],dWdy[i,j],d2Wdy2[i,j],d2Wdxdy[i,j],TT[i,j],dTdx[i,j],d2Tdx2[i,j],dTdy[i,j],d2Tdy2[i,j],d2Tdxdy[i,j],UU[i,j],dUdx[i,j],d2Udx2[i,j],dUdy[i,j],d2Udy2[i,j],d2Udxdy[i,j],J,idx,WBC,rh,Mh,a0,v,funcidx,type[funcidx],ordx,ordy)
                        else
                            BC_Horizon_Jac_Spin!(i,j+1,ff[i,j],dfdx[i,j],d2fdx2[i,j],dfdy[i,j],d2fdy2[i,j],d2fdxdy[i,j],gg[i,j],dgdx[i,j],d2gdx2[i,j],dgdy[i,j],d2gdy2[i,j],d2gdxdy[i,j],hh[i,j],dhdx[i,j],d2hdx2[i,j],dhdy[i,j],d2hdy2[i,j],d2hdxdy[i,j],WW[i,j],dWdx[i,j],d2Wdx2[i,j],dWdy[i,j],d2Wdy2[i,j],d2Wdxdy[i,j],TT[i,j],dTdx[i,j],d2Tdx2[i,j],dTdy[i,j],d2Tdy2[i,j],d2Tdxdy[i,j],UU[i,j],dUdx[i,j],d2Udx2[i,j],dUdy[i,j],d2Udy2[i,j],d2Udxdy[i,j],J,idx,WBC,rh,Mh,a0,v,funcidx,type[funcidx],ordx,ordy)
                        end
                    elseif i==Nx
                        #DEFINE BCS AT INFINITY
                        if !spin
                            BC_Infinity_Jac!(i,j+1,ff[i,j],dfdx[i,j],d2fdx2[i,j],dfdy[i,j],d2fdy2[i,j],d2fdxdy[i,j],gg[i,j],dgdx[i,j],d2gdx2[i,j],dgdy[i,j],d2gdy2[i,j],d2gdxdy[i,j],hh[i,j],dhdx[i,j],d2hdx2[i,j],dhdy[i,j],d2hdy2[i,j],d2hdxdy[i,j],WW[i,j],dWdx[i,j],d2Wdx2[i,j],dWdy[i,j],d2Wdy2[i,j],d2Wdxdy[i,j],TT[i,j],dTdx[i,j],d2Tdx2[i,j],dTdy[i,j],d2Tdy2[i,j],d2Tdxdy[i,j],UU[i,j],dUdx[i,j],d2Udx2[i,j],dUdy[i,j],d2Udy2[i,j],d2Udxdy[i,j],J,idx,WBC,rh,Mh,a0,v,funcidx,type[funcidx],ordx,ordy)
                        else
                            BC_Infinity_Jac_Spin!(i,j+1,ff[i,j],dfdx[i,j],d2fdx2[i,j],dfdy[i,j],d2fdy2[i,j],d2fdxdy[i,j],gg[i,j],dgdx[i,j],d2gdx2[i,j],dgdy[i,j],d2gdy2[i,j],d2gdxdy[i,j],hh[i,j],dhdx[i,j],d2hdx2[i,j],dhdy[i,j],d2hdy2[i,j],d2hdxdy[i,j],WW[i,j],dWdx[i,j],d2Wdx2[i,j],dWdy[i,j],d2Wdy2[i,j],d2Wdxdy[i,j],TT[i,j],dTdx[i,j],d2Tdx2[i,j],dTdy[i,j],d2Tdy2[i,j],d2Tdxdy[i,j],UU[i,j],dUdx[i,j],d2Udx2[i,j],dUdy[i,j],d2Udy2[i,j],d2Udxdy[i,j],J,idx,WBC,rh,Mh,a0,v,funcidx,type[funcidx],ordx,ordy)
                        end
                    else
                        #DEFINE FIELD EQUATIONS EVERYWHERE ELSE ON THE GRID
                        Field_Eqs_Jac!(i,j+1,ff[i,j],dfdx[i,j],d2fdx2[i,j],dfdy[i,j],d2fdy2[i,j],d2fdxdy[i,j],gg[i,j],dgdx[i,j],d2gdx2[i,j],dgdy[i,j],d2gdy2[i,j],d2gdxdy[i,j],hh[i,j],dhdx[i,j],d2hdx2[i,j],dhdy[i,j],d2hdy2[i,j],d2hdxdy[i,j],WW[i,j],dWdx[i,j],d2Wdx2[i,j],dWdy[i,j],d2Wdy2[i,j],d2Wdxdy[i,j],TT[i,j],dTdx[i,j],d2Tdx2[i,j],dTdy[i,j],d2Tdy2[i,j],d2Tdxdy[i,j],UU[i,j],dUdx[i,j],d2Udx2[i,j],dUdy[i,j],d2Udy2[i,j],d2Udxdy[i,j],haxis[i],siny[j],cosy[j],J,idx,WBC,rh,Mh,a0,v,funcidx,type[funcidx],ordx,ordy,X[i],Y[j])
                    end
                    idx+=NFields
                end
            end
            ordx+=1
            if mod(v,Nx)==0
                funcidx+=1
                ordx=0
            end
            if mod(v,Nx*NFields)==0
                funcidx=1
                ordx=0
                ordy+=1
            end
        end
        J .= transpose(J)
    end
end

#SOLVE THE NON-LINEAR SYSTEM TO OBTAIN THE SOLUTIONS
function solve_system(f::Field,g::Field,h::Field,W::Field,T::Field,U::Field,tol::Double64,WBC::Double64,rh::Double64,Mh::Double64,a0::Double64,spin::Bool=false,show_trace::Bool=true,iterations::Int64=30,method = :newton)
    return nlsolve(only_fj!((R,J,a)->residual_jac!(R,J,a,f,g,h,W,T,U,WBC,rh,Mh,a0,spin)), vcat(f.a,g.a,h.a,W.a,T.a,U.a), ftol=0.0,xtol=tol, show_trace=show_trace, iterations = iterations, method = method)
end
nothing
end