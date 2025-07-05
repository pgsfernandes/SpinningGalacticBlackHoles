function Ergosphere(f::Field, g::Field, h::Field, W::Field, rh::Double64=1.0, file_name::String = "ergosphere.dat", Npoints::Integer=100; q::Double64=0.0, compareMBH=false)
    sol=Array{Double64}(undef,Npoints+2)
    points=sort(append!(collect(pi/(2*Npoints) .* ((0:Npoints-1) .+ 1/2)),[0.0,pi/2]))

    V=Array{Double64}(undef,0,2)

    if !compareMBH
        M=GetMass(f,g,h,W,rh)
        j=GetJ(f,g,h,W,rh)
        χ = j/M^2
        rhK=M/2*sqrt(1-χ^2-q^2)
    else
        M, J, χ, Th, Ah, Ωh = get_quantities(f,g,h,W,rh)
        MBH = 0.5*Th*Ah + 2*J*Ωh
        χ = J/MBH^2
        M=MBH
    end

    function xErgoKerrN(y::Double64,M::Double64,χ::Double64,rh::Double64,q::Double64)
        return rBLKerr2x(M*(1+sqrt(1-q^2-χ^2*cos(y)^2)),rh,χ,q)
    end

    println("Printing ergosphere to file...")

    for i in 1:Npoints+2
        y=points[i]
        
        if i==1
            sol[i] = fzero(x->gtt(x,y,f,g,h,W), -1.0)
        else
            sol[i] = fzero(x->gtt(x,y,f,g,h,W), xErgoKerrN(y,M,χ,rhK,q))
        end

        r=2/(1-sol[i])
        x=r*sin(y)
        z=r*cos(y)

        V=vcat(V,[x z])

    end

    PrintData(file_name,V)
end