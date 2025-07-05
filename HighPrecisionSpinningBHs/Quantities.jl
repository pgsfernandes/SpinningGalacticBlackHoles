@changeprecision Double64 begin
#Computes ADM mass of the solution
function GetMass(f::Field,g::Field,h::Field,W::Field,rh::Double64)
    return rh*(2.0+f(Nx,1,dx=1))
end

#Computes angular momentum of the solution. Assumes ds^2 contains (d\phi - W/r^2 dt)^2
function GetJ(f::Field,g::Field,h::Field,W::Field,rh::Double64, Nx::Integer = Nx, Ny::Integer = Ny)
    return -rh*W(Nx,Ny+2,dx=1)
end

#Computes event horizon temperature of the solution
function GetTh(f::Field,g::Field,h::Field,rh::Double64)
    return 1/(16.0*pi*rh)*f(1,1)/sqrt(g(1,1)*h(1,1))
end

#Computes area of the event horizon of the solution
function GetAh(f::Field,g::Field,h::Field,rh::Double64)
    #Th is constant on the horizon
    integral, err = pquadrature(y -> sin(Double64(y))*sqrt(g(Double64("-1.0"),Double64(y))), Double64("0"), Double64(pi), abstol=Double64(1e-14))
    return 2.0*rh/(GetTh(f,g,h,rh)) * integral
end

#THEORETICAL VALUES FOR M, J, T_H, A_H FOR THE KERR BLACK HOLE
function GetωχKerr(WBC::Double64,rh::Double64,spin=true,branch::Integer=1)
    ω=0.0
    χ=0.0
    if spin
        χ=WBC
        if WBC!=0
            ω = (-1+WBC^2+sqrt(1-WBC^2))/(4*WBC)
        end
    else
        ω = WBC/rh
        aux = sqrt(Complex(-3+528*ω^2+768*ω^4))
        if branch == 2
            χ=real((3+(4*ω+(-8*ω*(9+8*ω^2)+3*aux)^(1/3))^2)/(3*(-8*ω*(9+8*ω^2)+3*aux)^(1/3)))
        else
            χ=real(1/6*(16*ω+(im*(im+sqrt(3))*(3+16*ω^2))/(-8*ω*(9+8*ω^2)+3*aux)^(1/3)-(1+im*sqrt(3))*(-8*ω*(9+8*ω^2)+3*aux)^(1/3)))
        end
    end
    return [ω, χ]
end

function GetωχKerrN(WBC::Double64,rh::Double64,q::Double64,spin=true,branch::Integer=1)
    ω=Double64("0.0")
    χ=Double64("0.0")
    if spin
        χ=WBC
        if WBC!=0
            ω = χ*sqrt(1-q^2-χ^2)/(4-2*q^2+4*sqrt(1-q^2-χ^2))
        end
    else
        ω = WBC/rh
        aux = sqrt(Complex(-3+528*ω^2+768*ω^4))
        if branch == 2
            χ=real((3+(4*ω+(-8*ω*(9+8*ω^2)+3*aux)^(1/3))^2)/(3*(-8*ω*(9+8*ω^2)+3*aux)^(1/3)))
        else
            χ=real(1/6*(16*ω+(im*(im+sqrt(3))*(3+16*ω^2))/(-8*ω*(9+8*ω^2)+3*aux)^(1/3)-(1+im*sqrt(3))*(-8*ω*(9+8*ω^2)+3*aux)^(1/3)))
        end
        χ=find_zero(x->ω - x*sqrt(1-q^2-x^2)/(4-2*q^2+4*sqrt(1-q^2-x^2)), χ)
    end
    return [ω, χ]
end

function quantities_kerr(WBC::Double64,rh::Double64,spin=true,branch::Integer=1)
    ω, χ = GetωχKerr(WBC,rh,spin,branch)

    M = 2*rh/sqrt(1-χ^2)
    J = χ*M^2
    Th = 1/(4*pi*M * (1+M/(2*rh)))
    Ah = 8*pi*M*(M+2*rh)
    return [M J χ Th Ah ω/rh]
end

#Returns all the above quantitites along with the dimensionless spin j/M^2
function get_quantities(f::Field,g::Field,h::Field,W::Field,rh::Double64)
    M=GetMass(f,g,h,W,rh)
    j=GetJ(f,g,h,W,rh)
    chi = j/M^2
    Th=GetTh(f,g,h,rh)
    Ah=GetAh(f,g,h,rh)
    Ωh=W(1,1)/rh^2
    return[M j chi Th Ah Ωh]
end

function CircumferencialRadius(x::Double64, f::Field, g::Field, h::Field, W::Field, rh::Double64)
    return sqrt(gphiphi(x,pi/2,f,g,h,W,rh))
end
end