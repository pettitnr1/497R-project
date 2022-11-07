using VortexLattice
using LinearAlgebra
using Optim

include("plots_default.jl")

"""
    optimize_airframe()

optimizes an airframe with a span of 1.5 m so that it can lift 0.5 kg, is stable, and most efficient

# Arguments
No arguments

# Returns
No returns, but it does print out the ideal chord, span, taper, lift coefficient, lift, and length from wing to tail. It also prints the velocity it takes to produce the necessary lift
"""
function optimize_airframe()

    chords = range(0.1, 1.5, 10)   # chord lengths ranging from 0.1 m to 1.5 m
    lengths = range(0, 2, 10)   # lengths from wing to tail ranging from 0 m to 2 m
    tapers = range(0.1, 1.0, 5)   # taper coefficients ranging from 0.1 to 1

    # initialize variables for use within for loops
    ideal_b = 1.5  # max span is 1.5 m
    ideal_c = [0.1]   # ideal chord length (m)
    CL = [0.0]   # lift coefficient
    ideal_CL = [vortex_lattice([ideal_b, ideal_c[1]], 1.0, 1.0, 1.0, "lift")]  # calculate ideal lift coefficient using ideal span and chord length
    e = 0.0   # inviscid span efficiency
    max_e = 0.0   # max efficiency
    ideal_taper = [1.0]   # ideal taper coefficient
    Sref = [ideal_b*ideal_c[1]]   # calculate reference area using ideal span and chord length (m^2)
    L = 4.905  # what lift needs to be to lift 0.5 kg
    V = [0.0]  # freestream velocity (m/s)
    ideal_V = [100.0]  # ideal freestream velocity (m/s)

    # stability derivatives
    CLa = [0.0]
    Cma = [0.0]
    Clb = [0.0]
    Cnb = [0.0]
    ideal_Cnb = [0.0]
    stability_derivatives = [[]]

    ideal_length = [0.0]   # ideal length from wing to tail (m)

    for c in chords
        for length in lengths
            # calculate lift coefficient for given chord length and distance from wing to tail
            CL[1] = vortex_lattice([ideal_b, c], 1.0, 1.0, length, "lift")
            Sref[1] = ideal_b*c[1]   # calculate reference area

            # throw out lift coefficients that come out as negative
            if CL[1] <= 0
                V[1] = 1000
            else
                V[1] = sqrt(L / (0.5*CL[1]*1.225*Sref[1]))   # calculate lift
            end

            # retrieve stability derivatives
            stability_derivatives = vortex_lattice([ideal_b, c], 1.0, 1.0, length, "stability")
            CLa[1] = stability_derivatives[1][1]
            Cma[1] = stability_derivatives[1][2]
            Clb[1] = stability_derivatives[2][1]
            Cnb[1] = stability_derivatives[2][2]

            # airframe is ideal if it is most stable and velocity needed for necessary lift is low
            if (V[1] < ideal_V[1] && CLa[1] > 0 && Cma[1] < 0 && Clb[1] < 0 && Cnb[1] > ideal_Cnb[1])
                ideal_c[1] = c
                ideal_CL[1] = CL[1]
                ideal_length[1] = length
                ideal_V[1] = V[1]
                ideal_Cnb[1] = Cnb[1]
            end
        end
    end

    for taper in tapers
        mean_c = ideal_c[1]*(1+taper)/2   # calculate mean chord length for given taper

        # retrieve efficiency for airframe with given taper
        e = wing_efficiency(ideal_b, mean_c, taper)

        # taper is ideal if efficiency is closest to, but less than 1
        if (e > max_e && e < 1)
            max_e = e
            ideal_taper[1] = taper
        end
    end

    mean_c = ideal_c[1]*(1+ideal_taper[1])/2   # calculate mean chord length using ideal taper
    Sref[1] = mean_c * ideal_b   # recalculate reference area using mean chord length

    # retrieve lift coefficient for ideal airframe
    ideal_CL[1] = vortex_lattice([ideal_b, mean_c], 1.0, ideal_V[1], ideal_length[1], "lift", ideal_taper[1])
    ideal_V[1] = sqrt(L / (0.5*ideal_CL[1]*1.225*Sref[1]))   # calculate freestream velocity needed for ideal airframe++

    # report if airframe was not able to be made stable
    if (ideal_length[1] == 0)
        println("NOT STABLE")
    end

    # print out all pertinent information for ideal airframe
    println("b = ", ideal_b, " m")
    println("c = ", ideal_c[1], " m")
    println("taper = ", ideal_taper[1])
    println("CL = ", ideal_CL[1])
    println("V = ", ideal_V[1], " m/s")
    println("Lift = ", L, " N")
    println("Length = ", ideal_length[1], " m")
end

"""
    wing_efficiency(b, c, taper)

computes efficiency of wing

# Arguments
- `b`: wing span (m)
- `c`: mean chord length (m)
- `taper`: taper coefficient used to determine tip chord length

# Returns
- `efficiency`: the efficiency of the given wing
"""
function wing_efficiency(b, c, taper)

    coefficients = zeros(2)   # construct array to hold lift and drag coefficients for every aspect ratio
    coefficients = vortex_lattice([b, c], 1.0, 1.0, 1.0, "coefficients", taper)
    cl = coefficients[1]   # get coefficient of lift from previous request
    cd = coefficients[2]   # get coefficient of drag from previous request
    ar = b/c
    efficiency = (cl^2) / (pi * ar * cd)   # compute inviscid span efficiency for given aspect ratio, using coeffcients of lift and drag
    return efficiency
end

"""
    vortex_lattice(airframe=[1.5,0.5], aoa=1.0, v=1.0, length=1.0, request="lift", taper=1.0)

performs all the necessary vortex lattice computations

# Arguments
- `airframe`: array with first value being span and second value being mean chord length, defaults to 1.5 m and 0.5 m respectively
- `aoa`: is the angle of attack of the airframe being evaluated, defaults to 1.0 degree
- `v`: is the velocity of the freestream, defaults to 1.0 m/s
- `length`: the distance from wing to tail, defaults to 1.0 m
- `request`: determines what will be returned by the function, as well as what calculations need to be done
- `taper`: taper coefficient used to determine tip chord length, defaults to 1.0

# Returns
- `CL`: the lift coefficient of the airframe, returned when request is "lift",
- `[long_stability, lat_stability]`: longitudinal and lateral stability derivatives, only includes pertinent derivatives, returned when request is "stability"
- `[CL, CD]`: lift and drag coefficients of airframe, returned when request is "coefficients"
"""
function vortex_lattice(airframe=[1.5,0.5], aoa=1.0, v = 1.0, length=1.0, request = "lift", taper=1.0)
    b = airframe[1]   # span length (m)
    c_ave = airframe[2]   # mean chord length (m)

    # wing
    xle = [0.0, (c_ave-c_ave*taper)/2]
    yle = [0.0, b/2]
    zle = [0.0, 0.0]
    chord = [c_ave, c_ave*taper]
    theta = [0.0, 0.0]
    phi = [0.0, 0.0]
    fc = fill((xc) -> 0, 2) # camberline function for each section
    mirror = true # account for geometric symmetry here because flow is not symmetric

    # horizontal stabilizer
    xle_h = [0.0, 0.0]
    yle_h = [0.0, 0.25]
    zle_h = [0.0, 0.0]
    chord_h = [0.2, 0.2]
    theta_h = [0.0, 0.0]
    phi_h = [0.0, 0.0]
    fc_h = fill((xc) -> 0, 2) # camberline function for each section
    ns_h = 6
    nc_h = 3
    spacing_s_h = Uniform()
    spacing_c_h = Uniform()
    mirror_h = true

    # vertical stabilizer
    xle_v = [0.0, 0.0]
    yle_v = [0.0, 0.0]
    zle_v = [0.0, 0.5]
    chord_v = [0.2, 0.2]
    theta_v = [0.0, 0.0]
    phi_v = [0.0, 0.0]
    fc_v = fill((xc) -> 0, 2) # camberline function for each section
    ns_v = 5
    nc_v = 3
    spacing_s_v = Uniform()
    spacing_c_v = Uniform()
    mirror_v = false

    # discretization parameters
    ns = 12 # number of spanwise panels
    nc = 6 # number of chordwise panels
    spacing_s = Uniform()
    spacing_c = Uniform()

    # length from wing to tail (m)
    l_v = length
    l_h = l_v

    # reference parameters
    Sref = b*c_ave   # calculate reference area
    bref = b
    cref = c_ave # chord
    rref = [chord[1]/4, 0.0, 0.0]   # center of gravity of airframe
    Vinf = v
    ref = Reference(Sref, cref, bref, rref, Vinf)

    # freestream parameters
    alpha = aoa*pi/180
    beta = 1.0*pi/180
    Omega = [0.0; 0.0; 0.0]
    fs = Freestream(Vinf, alpha, beta, Omega)

    # construct surface SOMETHING IS WRONG WITH THIS
    grid, wing = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, fc = fc, spacing_s=spacing_s, spacing_c=spacing_c)

    # generate surface panels for horizontal tail
    hgrid, htail = wing_to_surface_panels(xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h;
        mirror=mirror_h, fc=fc_h, spacing_s=spacing_s_h, spacing_c=spacing_c_h)

    # translate tail for horizontal tail volume ratio
    VortexLattice.translate!(hgrid, [l_h, 0.0, 0.0])
    VortexLattice.translate!(htail, [l_h, 0.0, 0.0])

    # generate surface panels for vertical tail
    vgrid, vtail = wing_to_surface_panels(xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v;
        mirror=mirror_v, fc=fc_v, spacing_s=spacing_s_v, spacing_c=spacing_c_v)

    # translate tail for vertical tail volume ratio
    VortexLattice.translate!(vgrid, [l_v, 0.0, 0.0])
    VortexLattice.translate!(vtail, [l_v, 0.0, 0.0])

    # create vector containing all surfaces
    surfaces = [wing, htail, vtail]
    surface_id = [1, 2, 3]

    # symmetry cannot be used because flow conditions are not symmetric
    symmetric = false

    # perform steady state analysis
    system = steady_analysis(surfaces, ref, fs; symmetric=symmetric, surface_id=surface_id)

    # retrieve near-field forces
    CF, CM = body_forces(system; frame=Wind())

    # perform far-field analysis
    CDiff = far_field_drag(system)

    # retrieve stability derivatives
    dCFs, dCMs = stability_derivatives(system)

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    # stability derivatives
    CDa, CYa, CLa = dCFs.alpha
    Cla, Cma, Cna = dCMs.alpha
    CDb, CYb, CLb = dCFs.beta
    Clb, Cmb, Cnb = dCMs.beta
    CDp, CYp, CLp = dCFs.p
    Clp, Cmp, Cnp = dCMs.p
    CDq, CYq, CLq = dCFs.q
    Clq, Cmq, Cnq = dCMs.q
    CDr, CYr, CLr = dCFs.r
    Clr, Cmr, Cnr = dCMs.r

    long_stability = [CLa, Cma]   # longitudinal stability
    lat_stability = [Clb, Cnb]   # lateral stability

    if request === "lift"
        return CL
    elseif request === "stability"
        return [long_stability, lat_stability]
    elseif request === "coefficients"
        return [CL, CD]
    end
end
