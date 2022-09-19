using Xfoil, Printf
include("plots_default.jl")

function alterReynolds()
    re = [50000, 100000, 200000, 500000, 1000000]
    n_r = length(re)
    alpha = -10:1:10
    cl = [Vector{Float64}(undef, 5) for _ in 1:length(alpha)]
    cd = [Vector{Float64}(undef, 5) for _ in 1:length(alpha)]
    cm = [Vector{Float64}(undef, 5) for _ in 1:length(alpha)]
    for i = 1:n_r
        coefficients = autoSweep(re[i], "reynolds")
        cl[i] = coefficients[:, 1]
        cd[i] = coefficients[:, 2]
        cm[i] = coefficients[:, 3]
    end

    l = plot(alpha, cl[1:5], xlabel=L"\alpha~\mathrm{(degrees)}", ylabel=L"c_\ell")
    d = plot(alpha, cd[1:5], xlabel=L"\alpha~\mathrm{(degrees)}", ylabel=L"c_d")
    m = plot(alpha, cm[1:5], xlabel=L"\alpha~\mathrm{(degrees)}", ylabel=L"c_m")
    plot(l, d, m)
    savefig("altered-reynolds.pdf")
end

function autoSweep(re = 1e5, mode="coefficients")
    # read airfoil into XFOIL
    open("naca2412.dat", "r") do f
        x = Float64[]
        y = Float64[]
        for line in eachline(f)
            entries = split(chomp(line))
            push!(x, parse(Float64, entries[1]))
            push!(y, parse(Float64, entries[2]))
        end
        Xfoil.set_coordinates(x,y)
    end

    # repanel using XFOIL's `PANE` command
    Xfoil.pane()

    # set operating conditions
    alpha = -10:1:10
    mach = 0.0

    # initialize outputs
    n_a = length(alpha)
    c_l = zeros(n_a)
    c_d = zeros(n_a)
    c_dp = zeros(n_a)
    c_m = zeros(n_a)
    converged = zeros(Bool, n_a)

    coefficients = zeros(n_a, 4)
    for i = 1:n_a
        c_l[i], c_d[i], c_dp[i], c_m[i], converged[i] = Xfoil.solve_alpha(alpha[i], re; mach, iter=100)
    end
    coefficients[:, 1] = c_l
    coefficients[:, 2] = c_d
    coefficients[:, 3] = c_m

    if (mode === "reynolds")
        return coefficients
    end

    # print results
    println("Angle\t\tCl\t\tCd\t\tCm\t\tConverged")
    for i = 1:n_a
      @printf("%8f\t%8f\t%8f\t%8f\t%d\n",alpha[i],c_l[i],c_d[i],c_m[i],converged[i])
    end
end
