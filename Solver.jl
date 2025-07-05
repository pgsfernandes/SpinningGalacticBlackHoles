@changeprecision Double64 begin

function OneSolution(;WBC::Double64, rh::Double64, Mh::Double64, a0::Double64, spin::Bool=true, tol::Double64=1e-11, guess=nothing, branch::Int64=1, ToPrint::Bool=false, ergosphere::Bool=false, light_ring::Bool=false, isco::Bool=false, petrov::Bool=false, sphericity::Bool=false, linvel::Bool=false)
    println()
    println("Solver initated with:")

    if spin
        println("χ=",WBC)
    else
        println("Ωh=",WBC)
    end
    println()

    #DECLARE OUR FIELDS
    f=Field();
    g=Field();
    h=Field();
    W=Field();
    T=Field();
    U=Field();

    fields = [f, g, h, W, T, U]

    if isnothing(guess)
        #IF THERE IS NO GUESS WE USE A COMPARABLE PERTURBED KERR BLACK HOLE AS INITIAL GUESS. WE PERTURB IT WITH THE PERTURBATION CONTROLED BY δ. THIS CAN BE USEFUL WHEN TESTING STUFF 
        ωKN, χKN = GetωχKerrN(WBC,rh,0.0,spin,branch)
        println("χKN = ", χKN)
        interpolate!(f,(x,y)->fKerrN(x,y,rh,χKN,0.0),Nx,Ny)
        interpolate!(g,(x,y)->gKerrN(x,y,rh,χKN,0.0),Nx,Ny)
        interpolate!(h,(x,y)->hKerrN(x,y,rh,χKN,0.0),Nx,Ny)
        interpolate!(W,(x,y)->WKerrN(x,y,rh,χKN,0.0),Nx,Ny)
    else
        #USE guess AS INITIAL GUESS IF THE USER PROVIDES ONE
        Nxguess, Nyguess = size(guess)
        Nxguess /= length(fields)  # Since guess is 6Na x Nb
        Nxguess = floor(Int,Nxguess)

        for (i, var) in enumerate(fields)
            block = guess[(i-1)*Nxguess+1 : i*Nxguess, 1:Nyguess]

            # Prepare resized array of target size
            resized = zeros(Nx, Ny)
            rows = min(Nxguess, Nx)
            cols = min(Nyguess, Ny)
            resized[1:rows, 1:cols] .= block[1:rows, 1:cols]

            var.a .= resized
        end

    end

    #THE LINE BELOW COMPUTES THE PHYSICAL QUANTITIES FOR A COMPARABLE KERR BH
    Mkerr, Jkerr, χkerr, Thkerr, Ahkerr, Ωhkerr = quantities_kerr(WBC,rh,spin,branch)

    #HERE WE SOLVE THE SYSTEM OF PDEs
    sol=solve_system(f,g,h,W,T,U,tol,WBC,rh,Mh,a0,spin,true,30)
    convergence=sol.f_converged || sol.x_converged

    if convergence
        println()
        println("Success! Solution converged!")
        #OBTAIN PHYSICAL QUANTITIES OF OUR SOLUTION
        M, J, χ, Th, Ah, Ωh = get_quantities(f,g,h,W,rh)
        MBH=0.5*Th*Ah + 2*Ωh*J

        χBH=J/MBH^2

        Mhalo = M-MBH

        rh /= MBH
        M /= MBH
        Mhalo /= MBH
        Mhguess = Mh/MBH
        a0 /= MBH
        Th *= MBH
        Ah /= MBH^2
        Ωh *= MBH

        W.a = W.a ./ MBH
        T.a = T.a .* MBH^2
        U.a = U.a .* MBH^2

        println("χBH=$χBH")
        
        if χBH<1.0
            Ahkerr = 8*π*(1+sqrt(1-χBH^2))
            Thkerr = sqrt(1-χBH^2)/(4.0*π*(1.0+sqrt(1.0-χBH^2)))

            println("T/Tkerr-1 = ",Th/Thkerr-1)
            println("A/Akerr-1 = ",Ah/Ahkerr-1)
        end

        println()
        println("Mhalo/MBH = ", Mhalo)
        println("a0/MBH = ", a0)
        println()

        #ERROR ESTIMATES
        println("Mkomar = ",-rh*gtt(1.0,pi/2,f,g,h,W,rh,dx=1))
        println("Madm = ",-rh*grr(1.0,pi/2,f,g,h,W,rh,dx=1))
        δM=grr(1.0,pi/2,f,g,h,W,rh,dx=1)/gtt(1.0,pi/2,f,g,h,W,rh,dx=1)-1.0;
        println("Madm/Mkomar-1 (pi/2) = ",δM)
        println("Madm/Mkomar-1 (0) = ",grr(1.0,0.0,f,g,h,W,rh,dx=1)/gtt(1.0,0.0,f,g,h,W,rh,dx=1)-1.0)

        #println("b/bteo = ",grr(0.4,0.0,f,g,h,W,rh)/b(0.4,rh,Mh*(a0+rh),a0)^4-1.0) #useful in the static case

        Th0=1/(16.0*pi*rh)*f(-1.0,0.0)/sqrt(g(-1.0,0.0)*h(-1.0,0.0))
        Thp2=1/(16.0*pi*rh)*f(-1.0,pi/2.0)/sqrt(g(-1.0,pi/2.0)*h(-1.0,pi/2.0))
        δk=Th0 / Thp2-1.0;

        println("Th(0) / Th(π/2) - 1= ", δk)

        if ergosphere
            Ergosphere(f,g,h,W,rh,"ergosphere.dat",100)
        end

        lrqts=nothing
        if light_ring
            lrqts=LightRing(f,g,h,W,rh,compareMBH=true)
        end

        iscoqts=nothing
        if isco
            iscoqts=ISCO(f,g,h,W,rh,compareMBH=true)
        end

        if ToPrint
            #In units MBH=1
            v=[rh,M,Mhalo,Mhguess,a0,χBH,Th,Ah,Ωh]
            verror=[Nx,Ny,δM,δk]
            dirstring = @sprintf("../Galactic Solutions/j=%.3f, Mh=%.3f, a0=%.3f, Nx=%d, Ny=%d", abs(χBH), Mhalo, a0, Nx, Ny)
            mkpath(dirstring)
            PrintData(string(dirstring,"/","f.dat"),f.a)
            PrintData(string(dirstring,"/","g.dat"),g.a)
            PrintData(string(dirstring,"/","h.dat"),h.a)
            PrintData(string(dirstring,"/","W.dat"),W.a)
            PrintData(string(dirstring,"/","T.dat"),T.a)
            PrintData(string(dirstring,"/","U.dat"),U.a)
            write_data_to_file(string(dirstring,"/","qts.dat"),v)
            write_data_to_file(string(dirstring,"/","error.dat"),verror)
            if !isnothing(iscoqts)
                write_data_to_file(string(dirstring,"/","iscoqts.dat"),iscoqts)
            end
            if !isnothing(lrqts)
                write_data_to_file(string(dirstring,"/","lrqts.dat"),lrqts)
            end
        end
        return convergence, rh, Mhalo, Mhguess, a0, Ωh, f.a, g.a, h.a, W.a, T.a, U.a
    else
        println("Did not converge...")
    end
end

function b(x,rh,c1,c2)
    return ((3 - x)*(1 + (c1*((c2*(-1 + x)*(-24*rh^2 + c2^2*(-3 + x) + c2*rh*(-17 + 11*x)))/(-2*rh + c2*(-1 + x)) + 2*rh*(3*rh*(-1 + x)*log(rh/(c2 + rh)) - 2*(2*c2 + 3*rh)*log((2*rh)/(c2 + 2*rh - c2*x)) + c2*(-1 + x)*log((2*rh^2)/((c2 + rh)*(c2 + 2*rh - c2*x))))))/(2*c2^4*(-3 + x))))/2
end

#ENERGY DENSITY PROFILE
function e(x,y,rh,Mh,a0)

    c1 = Mh*(a0+rh)
    c2 = a0

    return (-1+x)^4 * (1+x)^2 * c1 / ( 16*pi*rh * (2*rh + c2 - x*c2)^3 * b(x,rh,c1,c2)^5 )
end

function read_sol()
    #TO READ GUESS, PLACE THE FOLLOWING DATA FILES IN THE SAME FOLDER AS THIS FILE
    rh, M, Mhalo, Mhguess, a0, χBH, Th, Ah, Ωh=Double64.(readdlm("qts.dat",','));
    guess = Double64.(vcat(readdlm("f.dat"),readdlm("g.dat"),readdlm("h.dat"),readdlm("W.dat"),readdlm("T.dat"),readdlm("U.dat")))
    return rh, M, Mhalo, Mhguess, a0, χBH, Th, Ah, Ωh, guess
end

function SolGuess()
    rh, M, Mhalo, Mhguess, a0, χBH, Th, Ah, Ωh, guess = read_sol()

    WBC = Ωh*rh^2

    @time OneSolution(WBC=WBC,rh=rh,Mh=Mhguess,a0=a0,spin=false,tol=1e-10,branch=1,ToPrint=true,ergosphere=false,light_ring=false,isco=false, guess=guess)
end

function SolNoGuess()
    #USES KERR SOLUTION WITH SAME MASS AND SPIN AS INITIAL GUESS
    a0=100.0
    Mh=10.0
    rh=0.5
    WBC=0.01
    @time OneSolution(WBC=WBC,rh=rh,Mh=Mh,a0=a0,spin=false,tol=1e-10,branch=1,ToPrint=false,ergosphere=false,light_ring=false,isco=false)
end

SolNoGuess()

end
nothing