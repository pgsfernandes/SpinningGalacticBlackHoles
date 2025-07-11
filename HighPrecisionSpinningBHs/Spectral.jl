@changeprecision Double64 begin
#APPROXIMATING A 1D FUNCTION f VIA CHEBYSHEV INTERPOLATION. RETURNS THE SPECTRAL COEFFICIENTS.
function interpolate1D(f::Function,Nx::Integer,Ny::Integer)
    c=zeros(Double64, (Nx,Ny))
    x = Array{Double64, 1}(undef, Nx)
    cheb=Array{Double64, 2}(undef, Nx, Nx)
    for i in 1:Nx
        x[i] = -cos((2*i-1.0)/(2*Nx)*pi)
    end
    for k in 0:Nx-1
        for i in 0:Nx-1
            cheb[k+1,i+1] = chebyshevt(k,x[i+1])
        end
    end
    for j in 0:(Nx-1)
        for x_it in 1:Nx
            c[j+1,1] += 2.0/Nx * f(x[x_it]) * cheb[j+1,x_it]
        end
    end
    c[1,1]/=2.0
    return c
end
function interpolate1D!(F::Field,f::Function,Nx::Integer,Ny::Integer)
    F.a=zeros(Double64, (Nx,Ny))
    x = Array{Double64, 1}(undef, Nx)
    cheb=Array{Double64, 2}(undef, Nx, Nx)
    for i in 1:Nx
        x[i] = -cos((2*i-1.0)/(2*Nx)*pi)
    end
    for k in 0:Nx-1
        for i in 0:Nx-1
            cheb[k+1,i+1] = chebyshevt(k,x[i+1])
        end
    end
    for j in 0:(Nx-1)
        for x_it in 1:Nx
            (F.a)[j+1,1] += 2.0/Nx * f(x[x_it]) * cheb[j+1,x_it]
        end
    end
    (F.a)[1,1]/=2.0
end

#APPROXIMATING A 2D FUNCTION f VIA CHEBYSHEV INTERPOLATION. RETURNS THE SPECTRAL COEFFICIENTS.
function interpolate(f::Function,Nx::Integer,Ny::Integer,type::Integer)
    c=zeros(Double64, (Nx,Ny))
    x = Array{Double64, 1}(undef, Nx)
    y = Array{Double64, 1}(undef, Ny)
    cheb=Array{Double64, 2}(undef, Nx, Nx)
    trig=Array{Double64, 2}(undef, Ny, Ny)
    
    for i in 1:Nx
        x[i] = -cos((2*i-1.0)/(2*Nx)*pi)
    end
    for i in 0:Ny-1
        y[i+1] = pi*(i+0.5)/(2*Ny)
    end

    for k in 0:Nx-1
        for i in 0:Nx-1
            cheb[k+1,i+1] = chebyshevt(k,x[i+1])
        end
    end

    for k in 0:Ny-1
        for i in 0:Ny-1
            trig[k+1,i+1] = dTrig(k,y[i+1],0,type)
        end
    end
    
    for j in 0:(Nx-1)
        for k in 0:(Ny-1)
            for x_it in 1:Nx
                for y_it in 1:Ny
                    if k==0 && (type == 2 || type == 4)
                        c[j+1,k+1] += 8.0/(Nx*Ny) * f(x[x_it],y[y_it]) * cheb[j+1,x_it] * trig[k+1,y_it]
                    else
                        c[j+1,k+1] += 4.0/(Nx*Ny) * f(x[x_it],y[y_it]) * cheb[j+1,x_it] * trig[k+1,y_it]
                    end
                end
            end
        end
    end
    
    for j in 1:Nx
        c[j,1] /= 2.0
    end
    for k in 1:Ny
        c[1,k] /= 2.0
    end
    
    return c
end

function interpolate!(F::Field,f::Function,Nx::Integer,Ny::Integer)
    type=F.type
    x = Array{Double64, 1}(undef, Nx)
    y = Array{Double64, 1}(undef, Ny)
    cheb=Array{Double64, 2}(undef, Nx, Nx)
    trig=Array{Double64, 2}(undef, Ny, Ny)

    aux = type==3 ? 1 : 0
    
    for i in 1:Nx
        x[i] = -cos((2*i-1.0)/(2*Nx)*pi)
    end
    for i in 0:Ny-1
        y[i+1] = pi*(i+0.5)/(2*Ny)
    end

    for k in 0:Nx-1
        for i in 0:Nx-1
            cheb[k+1,i+1] = chebyshevt(k,x[i+1])
        end
    end

    for k in 0:Ny-1
        for i in 0:Ny-1
            trig[k+1,i+1] = dTrig(k+aux,y[i+1],0,type)
        end
    end
    
    coeff_aux = 4.0 / (Nx * Ny)

    for j in 0:(Nx-1)
        for k in 0:(Ny-1)
            coeff = (k == 0 && type!= 1) ? 2.0 * coeff_aux : coeff_aux
            for x_it in 1:Nx
                for y_it in 1:Ny
                    (F.a)[j+1,k+1] += coeff * f(x[x_it],y[y_it]) * cheb[j+1,x_it] * trig[k+1,y_it]
                end
            end
        end
    end
    
    for j in 1:Nx
        (F.a)[j,1] /= 2.0
    end
    for k in 1:Ny
        (F.a)[1,k] /= 2.0
    end
end

#Converts coefficients of 1D approximation such that they can be used in a 2D problem. In our context, takes a static solution and puts it in a form such that it can be used as an approximation for a slowly rotating axisymmetric solution.
function To2D(v::AbstractArray,Ny::Integer)
    return hcat(v,zeros(length(v),Ny-1))
end

#Returns Chebyshev nodes as our grid in the x direction
function xpoints(Nx::Integer)
    # Use Double64 for all calculations and literals
    points = cos.((2 .* (1:Nx-2) .- 1) ./ (2*(Nx-2)) * Double64(pi))
    return sort(append!(points, [Double64(-1.0), Double64(1.0)]))
end

#Equally spaced points in the range [0,pi/2]. Our grid in y direction
function ypoints(Ny::Integer)
    # Use Double64 for all calculations and literals
    interior_points = collect(Double64(pi) / (2 * Ny) .* ((0:Ny-1) .+ Double64(1/2)))
    return sort(append!(interior_points, [Double64(0.0), Double64(pi)/2]))
end

#COMPUTES THE CHEBYSHEV POLYNOMIAL OF ORDER n (OR ITS dx DERIVATIVE) AT POINT x
function dT(n::Integer,x::Double64,dx::Integer=0)
    if n<=0
        if dx==1 || dx==2
            return 0
        else
            return chebyshevt(-n,x)
        end
    else
        if dx==1
            return n*chebyshevu(n-1,x)
        elseif dx==2
            if x==-1
                return (n^4-n^2)/3 * (-1)^n
            elseif x==1
                return (n^4-n^2)/3
            else
                return n*(n*chebyshevt(n,x)-x*chebyshevu(n-1,x))/(-1+x^2)
            end
        else
            return chebyshevt(n,x)
        end
    end 
end

function dTrig(n::Integer, y, dy::Integer=0, type::Integer=1)
    big_n = Double64(n)  # Convert n to Double64
    y = Double64(y)  # Ensure y is Double64

    if dy == 0
        if type==1
            return cos(2*big_n*y)
        elseif type==2
            return cos((2*big_n+1)*y)
        elseif type==3
            return sin(2*big_n*y)
        else
            return sin((2*big_n+1)*y)
        end
    elseif dy==1
        if type==1
            return -2*big_n*sin(2*big_n*y)
        elseif type==2
            return -(2*big_n+1)*sin((2*big_n+1)*y)
        elseif type==3
            return 2*big_n*cos(2*big_n*y)
        else
            return (2*big_n+1)*cos((2*big_n+1)*y)
        end
    else
        if type==1
            return -4*big_n^2*cos(2*big_n*y)
        elseif type==2
            return -(2*big_n+1)^2*cos((2*big_n+1)*y)
        elseif type==3
            return -4*big_n^2*sin(2*big_n*y)
        else
            return -(2*big_n+1)^2*sin((2*big_n+1)*y)
        end
    end
end


#Evaluates all Chebyshev polynomials, first and second derivatives in all points of our x grid. This allows us to not repeat computations when evaluating functions on the grid.
function GetChebyshev(Nx::Integer,x::Vector{Double64})
    C = Double64.(zeros((3,Nx,Nx)))
    for i in 1:3
        for j in 1:Nx
            for k in 1:Nx
                #derivative, poly order, position
                C[i,j,k] = dT(j-1,x[k],i-1)
            end
        end
    end
    return C
end

#Evaluates all even cosines, first and second derivatives in all points of our y grid. This allows us to not repeat computations when evaluating functions on the grid.
function GetTrig(Ny::Integer,y::Vector{Double64})
    C = Double64.(zeros((3,Ny,Ny+2,4)))
    for t in 1:4
        aux = t==3 ? 1 : 0
        for i in 0:2
            for j in 0:Ny-1
                for k in 1:Ny+2
                    C[i+1,j+1,k,t] = dTrig(j+aux,y[k],i,t)
                end
            end
        end
    end
    return C
end
end