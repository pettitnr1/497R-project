using Xfoil
using Printf
include("plots_default.jl")

# function that evaluates airfoils of varying thicknesses and cambers
#  - `request` is a string that determines whether the thickness or camber
#        airfoils will be plotted
function alterThicknessCamber(request)
    # NACA data files for airfoils with different thicknesses and cambers
    alteredCambers = ["naca4412.dat", "naca6412.dat", "naca8412.dat"]
    alteredThicknesses = ["naca2410.dat", "naca2420.dat", "naca2430.dat"]
    n_c = length(alteredCambers)
    n_t = length(alteredThicknesses)
    alpha = -20:1:20
    cld = [Vector{Float64}(undef, n_c) for _ in 1:length(alpha)]
    lcs = [Vector{Float64}(undef, n_c) for _ in 1:length(alpha)]

    # calculate lift and drag coefficients for airfoils with different cambers
    for i = 1:n_c
        cld[i] = getCoefficients(alteredCambers[i], "ld")
        lcs[i] = getCoefficients(alteredCambers[i], "lcs")
    end

    # plot lift over drag and lift for altered camber airfoils
    ac1 = plot(alpha, cld[1:3], xlabel=L"\alpha~\mathrm{(degrees)}", ylabel=L"c_\ell/c_d", title="Altered Camber", titlefontsize=10, label = ["" "" ""])
    ac2 = plot(alpha, lcs[1:3], xlabel=L"\alpha~\mathrm{(degrees)}", ylabel=L"c_\ell", title="Altered Camber", titlefontsize=10, label = ["" "" ""])

    # calculate lift and drag coefficients for airfoils with different thickness
    for i = 1:n_c
        cld[i] = getCoefficients(alteredThicknesses[i], "ld")
        lcs[i] = getCoefficients(alteredThicknesses[i], "lcs")
    end

    # plot lift over drag and lift for altered thickness airfoils
    at1 = plot(alpha, cld[1:3], xlabel=L"\alpha~\mathrm{(degrees)}", ylabel=L"c_\ell/c_d", title="Altered Thickness", titlefontsize=10, label = ["" "" ""])
    at2 = plot(alpha, lcs[1:3], xlabel=L"\alpha~\mathrm{(degrees)}", ylabel=L"c_\ell", title="Altered Thickness", titlefontsize=10, label = ["" "" ""])

    # legend for plots
    st = plot(alpha .* NaN, cld[1:3] .* NaN, label=["10%" "20%" "30%"], showaxis=false)
    sc = plot(alpha .* NaN, cld[1:3] .* NaN, label=["4%" "6%" "8%"], showaxis=false)

    if (request==="thickness")
        plot(at1, at2, st)
    end

    if (request==="camber")
        plot(ac1, ac2, sc)
    end
end

# function that calulates lift, drag, and moment coefficients for given airfoil
#  - `altered` is the NACA airfoil being evaluated
#  - `request` is a string that determines what will be returned
function getCoefficients(altered, request)
    # read airfoil into XFOIL
    open(altered, "r") do f
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
    alpha = -20:1:20
    mach = 0.0
    re = 1e6
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
    coefficients[:, 4] = converged

    if (request === "ld")
        return camber_cld = c_l ./ c_d
    end
    if (request === "lcs")
        return c_l
    end
end

# function that evaluates airfoils with flows that have varying Reynold's numbers
function alterReynolds()
    re = [100000, 500000, 1000000, 2000000, 3000000]
    n_r = length(re)
    alpha = -20:1:20
    cl = [Vector{Float64}(undef, 5) for _ in 1:length(alpha)]
    cd = [Vector{Float64}(undef, 5) for _ in 1:length(alpha)]
    cm = [Vector{Float64}(undef, 5) for _ in 1:length(alpha)]

    # evaluate airfoil under different Reynold's numbers scenarios
    for i = 1:n_r
        coefficients = autoSweep(re[i], "coefficients")
        cl[i] = coefficients[:, 1]
        cd[i] = coefficients[:, 2]
        cm[i] = coefficients[:, 3]
    end

    # plot lift, drag, and moment coefficients for different Reynold's numbers
    l = plot(alpha, cl[1:5], xlabel=L"\alpha~\mathrm{(degrees)}", ylabel=L"c_\ell", label=["" "" "" "" ""])
    d = plot(alpha, cd[1:5], xlabel=L"\alpha~\mathrm{(degrees)}", ylabel=L"c_d", label=["" "" "" "" ""])
    m = plot(alpha, cm[1:5], xlabel=L"\alpha~\mathrm{(degrees)}", ylabel=L"c_m", label=["" "" "" "" ""])
    s = plot(alpha .* NaN, cl[1:5] .* NaN, label=["1e5" "5e5" "1e6" "2e6" "3e6"], showaxis=false)
    plot(l, d, m, s)
end

# function that calculates lift, drag, and moment coefficients for NACA2412 airfoil
#  - `re` is the Reynold's number being used
#  - `mode` is a string that determines what will be returned
function autoSweep(re = 1e5, mode="plot")
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
    alpha = -20:1:20
    mach = 0.0

    # initialize outputs
    n_a = length(alpha)
    c_l = zeros(n_a)
    c_d = zeros(n_a)
    c_dp = zeros(n_a)
    c_m = zeros(n_a)
    converged = zeros(Bool, n_a)

    coefficients = zeros(n_a, 4)

    # perform Xfoil calculations to retrieve coefficients
    for i = 1:n_a
        c_l[i], c_d[i], c_dp[i], c_m[i], converged[i] = Xfoil.solve_alpha(alpha[i], re; mach, iter=100)
    end
    coefficients[:, 1] = c_l
    coefficients[:, 2] = c_d
    coefficients[:, 3] = c_m
    coefficients[:, 4] = converged

    index = 1
    n_vc = length(converged[converged .=== true ])
    validCoefficients = zeros(n_vc, 4)

    # parse through calculated coefficients and only keep valid data
    for row in eachrow(coefficients)
        if (row[4] === 1.0)
            validCoefficients[index, 1:4] = row[1:4]
            index = index + 1
        end
    end

    c_l = validCoefficients[:, 1]
    c_d = validCoefficients[:, 2]
    c_m = validCoefficients[:, 3]
    converged = validCoefficients[:, 4]

    # return validated coefficients
    if (mode === "coefficients")
        return validCoefficients
    end

    # compare experimental to Xfoil data
    if ( mode === "compare")
        airfoilComparison()
    end

    # plot coefficients against angle of attack
    if (mode === "plot")
        cl = plot(alpha, c_l, xlabel=L"\alpha~\mathrm{(degrees)}", ylabel=L"c_\ell", label="")
        cd = plot(alpha, c_d, xlabel=L"\alpha~\mathrm{(degrees)}", ylabel=L"c_d", label="")
        cm = plot(alpha, c_m, xlabel=L"\alpha~\mathrm{(degrees)}", ylabel=L"c_m", label="")
        plot(cl, cd, cm)
    end

    # print results
    # println("Angle\t\tCl\t\tCd\t\tCm\t\tConverged")
    # for i = 1:n_a
    #   @printf("%8f\t%8f\t%8f\t%8f\t%d\n",alpha[i],c_l[i],c_d[i],c_m[i],converged[i])
    # end
end

# function that compares experimental data to data from Xfoil
function airfoilComparison()
    # experimental data
    alpha_exp = [0, 2.5, 7.5, 10, 12.5, 15]
    cl_exp = [0.228, 0.725, 0.88, 1.131, 1.227, 1.093]
    cd_exp = [0.016, 0.0424, 0.053, 0.069, 0.112, 0.378]

    # perform Xfoil calculations
    coefficients = autoSweep(1e5, "coefficients")
    alpha = 0:1:15
    n_a = length(alpha)
    cl_x = coefficients[11:25, 1]
    cd_x = coefficients[11:25, 2]

    # plot experimental data and Xfoil data together for comparison
    plot(alpha_exp, cl_exp, xlabel=L"\alpha~\mathrm{(degrees)}", ylabel=L"c_\ell", label="")
    cl = plot!(alpha, cl_x, label="")

    plot(alpha_exp, cd_exp, xlabel=L"\alpha~\mathrm{(degrees)}", ylabel=L"c_d", label="")
    cd = plot!(alpha, cd_x, label="")

    plot(alpha_exp .* NaN, cl_exp .* NaN, label="experimental", showaxis=false)
    s = plot!(alpha .* NaN, cl_x .* NaN, label="Xfoil")

    plot(cl, cd, s, layout=(2,2))
end
