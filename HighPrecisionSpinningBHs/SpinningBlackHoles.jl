module SpinningBlackHoles

using ClassicalOrthogonalPolynomials, DelimitedFiles, Distributions, Roots, Cubature, ChangePrecision, Printf, DoubleFloats

@changeprecision Double64 include("Field.jl")

@changeprecision Double64 include("Load.jl")

@changeprecision Double64 include("Other.jl")

@changeprecision Double64 include("Spectral.jl")

@changeprecision Double64 include("Quantities.jl")

@changeprecision Double64 include("KerrNewman.jl")

@changeprecision Double64 include("MetricFunctions.jl")

@changeprecision Double64 include("Ergosphere.jl")

@changeprecision Double64 include("LightRing.jl")

@changeprecision Double64 include("ISCO.jl")

export Field, LoadSystem, PrintData, interpolate1D, interpolate1D!, interpolate, interpolate!, GetMass, GetJ, GetTh, GetAh, GetωχKerr, GetωχKerrN, quantities_kerr, get_quantities, fKerrN, gKerrN, hKerrN, WKerrN, Print_Petrov, gtt, gtphi, gphiphi, grr, gthetatheta, Ergosphere, CircumferencialRadius, LightRing, ISCO, write_data_to_file, X, Mx, Y, My, Nx, Ny

end