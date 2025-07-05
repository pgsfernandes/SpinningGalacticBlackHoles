#Prints data array func to file file_name

function PrintData(file_name::String, func, N::Int=16)
    open(file_name, "w") do io
        #formatted_func = [@sprintf("%.*e", N, x) for x in func]  # Format each value to N decimal places
        formatted_func = [abs(x) < 1e-16 ? @sprintf("%.*e", N, 0.0) : @sprintf("%.*e", N, x) for x in func]
        writedlm(io, formatted_func)
    end
end

function write_data_to_file(filename,data)
    file = open(filename, "w")
        
    # Write the vector to the file
    writedlm(file, data', ", ")
    
    # Close the file
    close(file)
end

function x2r(x::Double64,rh::Double64)
    return 2*rh/(1-x)
end

function r2x(r::Double64,rh::Double64)
    return 1-2*rh/r
end

function rBLKerr2x(rbl::Double64,rh::Double64,χ::Double64=0.0,q::Double64=0.0)
    M=2*rh/sqrt(1-χ^2-q^2)
    riso = 1/2*(-M+rbl+sqrt((M-rbl)^2-4*rh^2))
    return r2x(riso,rh)
end