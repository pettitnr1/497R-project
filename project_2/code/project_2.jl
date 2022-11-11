using VortexLattice
using LinearAlgebra

include("plots_default.jl")

"""
    aoa_effect()

computes lift coefficients for airframes with varying angles of attack. Plots
the lift coefficient against the aspect ratios

# Arguments
No arguments

# Returns
No returns, but it does plot the lift coefficient against the angle of attack
"""
function aoa_effect()
    angles = -20:1:20   # angles of attack used range from -20 to 20 degrees
    n_a = length(angles)
    lift = zeros(n_a)   # array to hold lift coeffcients for every angle of attack
    index = 1   # initialize index used for iteration through lift array

    for a in angles
        # retrieve coeffcients of lift and drag for given angle of attack
        cl = vortex_lattice("coefficients", 7.5, 1, a)
        lift[index] = cl[1]   # insert lift coefficient into lift array
        index = index + 1   # set index to be one greater for next iteration
    end

    # plot lift coefficient against angle of attack using angle of attack array and lift array
    plot(angles, lift, xlabel="Angle of Attack(degrees)", ylabel=L"C_\ell", label="")
end

"""
    wing_efficiency()

computes the inviscid span efficiency for airframes with varying aspect ratios.
Plots the efficiency against the aspect ratios

# Arguments
No arguments

# Returns
No returns, but it does plot the efficiency against the aspect ratios
"""
function wing_efficiency()
    ratios = 3:1:15   # initialize aspect ratios to range from 3 to 15
    n_r = length(ratios)
    coefficients = zeros(n_r, 2)   # construct array to hold lift and drag coefficients for every aspect ratio
    efficiency = 0   # initialize efficiency to be 0 so that calculations can be done for it within for loop
    e = Vector{Float64}(undef, n_r)   # construct array for holding efficiencies of airframes with different aspect ratios
    index = 1   # initialize index used for iteration through coefficients array

    for ar in ratios
        # retrieve lift and drag coefficients for given aspect ratio
        coefficients[index, :] .= vortex_lattice("coefficients", ar)
        cl = coefficients[index, 1]   # get coefficient of lift from previous request
        cd = coefficients[index, 2]   # get coefficient of drag from previous request

        efficiency = (cl^2) / (pi * ar * cd)   # compute inviscid span efficiency for given aspect ratio, using coeffcients of lift and drag from above
        e[index] = efficiency   # push calculated efficiency for given aspect ratio to efficiency array
        index = index + 1   # set index to be one greater for next iteration
    end

    # plot efficiency against aspect ratio using efficiency and aspect ratio arrays
    plot(ratios, e, xlabel="Aspect Ratio", ylabel="Efficiency", label="")
end

"""
    stability_derivatives()

computes the stability derivatives for airframes with varying
vertical and horizontal tail volume ratios. Plots those derivatives against the
tail volume ratios

# Arguments
No arguments

# Returns
No returns, but it does plot the stability derivatives against their corresponding tail volume ratios
"""
function stability_derivatives()
    v_ratios = range(0.001, 0.0156, 30)   # vertical tail volume ratios range from 0.001 to 0.0156
    h_ratios = range(0.001, 0.1458, 30)   # horizontal tail volume ratios range from 0.001 to 0.1458
    v_r = length(v_ratios)
    h_r = length(h_ratios)
    v_derivatives = zeros(v_r, 2)   # constuct array for holding stability derivatives affected by changes to vertical tail volume ratios
    h_derivatives = zeros(h_r, 2)   # constuct array for holding stability derivatives affected by changes to horizontal tail volume ratios
    index = 1   # initialize index used for iteration through v_derivatives and h_derivatives arrays

    for r in v_ratios
        # retrieve stability derivatives for given vertical tail volume ratio
        v_derivatives[index, :] .= vortex_lattice("vderivatives", 7.5, r)
        index = index + 1   # set index to be one greater for next iteration
    end

    index = 1   # reset index to equal 1 for next for loop

    for r in h_ratios
        # retrieve stability derivatives for given horizontal tail volume ratio
        h_derivatives[index, :] .= vortex_lattice("hderivatives", 7.5, r)
        index = index + 1   # set index to be one greater for next iteration
    end

    # setup C_lb against vertical tail volume ratios plot
    lb = plot(v_ratios, v_derivatives[:, 1], xlabel="V-tail Volume Ratios", ylabel=L"C_{\ell{b}}", label="")
    # setup C_nb against vertical tail volume ratios plot
    nb = plot(v_ratios, v_derivatives[:, 2], xlabel="V-tail Volume Ratios", ylabel=L"C_{nb}", label="", ylims=(-0.45, 0.1))

    # setup C_La against horizontal tail volume ratios plot
    la = plot(h_ratios, h_derivatives[:, 1], xlabel="H-tail Volume Ratios", ylabel=L"C_{La}", label="", xlims=(0.04375, 0.1458), ylims=(4.806,4.81))
    # setup C_ma against horizontal tail volume ratios plot
    ma = plot(h_ratios, h_derivatives[:, 2], xlabel="H-tail Volume Ratios", ylabel=L"C_{ma}", label="")

    plot(lb, nb, la, ma, layout=(4,1)) # plot all 4 previously defined plots together
end

"""
    vortex_lattice(request, ar=7.5, v=1, aoa=1.0*pi/180)

performs all the necessary vortex lattice computations

# Arguments
- `request`: determines what will be returned by the function, as well as what
    calculations need to be done
- `ar`: is the aspect ratio of the airframe being evaluated, defaults to 7.5
- `v`: is the tail volume ratio of the ariframe being evaluated, defaults to 1
- `aoa`: is the angle of attack of the airframe being evaluated, defaults to
    1.0*pi/180

# Returns
- `coefficients`: an array with the lift and drag coefficient of the given airframe
- `v_derivatives`: the stability derivatives associated with the vertical tail
- `h_derivatives`: the stability derivatives associated with the horizontal tail
"""
function vortex_lattice(request, ar=7.5, v=1, aoa=1.0*pi/180)
    # initialize coefficients array
    coefficients = zeros(1, 2)

    # wing
    xle = [0.0, 0.4]
    yle = [0.0, 7.5]
    zle = [0.0, 0.0]
    chord = [2.2, 1.8]
    theta = [0.0, 0.0]
    phi = [0.0, 0.0]
    fc = fill((xc) -> 0, 2) # camberline function for each section
    if (request === "vderivatives" || request === "hderivatives")
        # account for geometric symmetry here because flow is not symmetric
        mirror = true
    else
        # symmetry is accounted for later
        mirror = false
    end

    # horizontal stabilizer
    xle_h = [0.0, 0.4]
    yle_h = [0.0, 1.25]
    zle_h = [0.0, 0.0]
    chord_h = [0.7, 0.7]
    theta_h = [0.0, 0.0]
    phi_h = [0.0, 0.0]
    fc_h = fill((xc) -> 0, 2) # camberline function for each section
    ns_h = 6
    nc_h = 3
    spacing_s_h = Uniform()
    spacing_c_h = Uniform()
    mirror_h = true

    # vertical stabilizer
    xle_v = [0.0, 0.4]
    yle_v = [0.0, 0.0]
    zle_v = [0.0, 1.0]
    chord_v = [0.7, 0.7]
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

    # reference parameters
    Sref = 30 # area
    bref = sqrt(ar * Sref) # span
    cref = Sref/bref # chord
    rref = [0.50, 0.0, 0.0]
    Vinf = 1.0
    ref = Reference(Sref, cref, bref, rref, Vinf)

    # freestream parameters
    alpha = aoa*pi/180
    beta = 1.0*pi/180
    Omega = [0.0; 0.0; 0.0]
    fs = Freestream(Vinf, alpha, beta, Omega)

    if (request === "vderivatives")
        # vertical tail volume ratio variables
        S_v = zle_v[2]*chord_v[1]   # calculate area of vertical tail
        l_v = (v*Sref*bref)/S_v   # calculate length from wing to tail
        l_h = l_v
    end

    if (request === "hderivatives")
        # horizontal tail volume ratio variables
        S_h = yle_h[2]*chord_h[1]   # calculate area of horizontal tail
        l_h = (v*Sref*(chord[1]+chord[2])/2)/S_h    # calculate length from wing to tail
        l_v = l_h
    end

    # construct surface
    grid, wing = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, fc = fc, spacing_s=spacing_s, spacing_c=spacing_c)

    # now set normal vectors manually
    ncp = avl_normal_vector([xle[2]-xle[1], yle[2]-yle[1], zle[2]-zle[1]], 2.0*pi/180)

    # overwrite normal vector for each wing panel
    for i = 1:length(wing)
        wing[i] = set_normal(wing[i], ncp)
    end

    if (request === "vderivatives" || request === "hderivatives")
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
    else
        # create vector containing all surfaces
        surfaces = [wing]
        surface_id = [1]
    end

    if (request === "vderivatives" || request === "hderivatives")
        # symmetry cannot be used because flow conditions are not symmetric
        symmetric = false
    else
        # we can use symmetry since the geometry and flow conditions are symmetric about the X-Z axis
        symmetric = true
    end

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

    coefficients = [CL, CD]

    v_derivatives = [Clb, Cnb]
    h_derivatives = [CLa, Cma]

    if request === "coefficients"
        # return for aoaEffect() and wingEfficiency()
        return coefficients
    end
    if request === "vderivatives"
        # return for stabilityDerivatives()
        return v_derivatives
    end
    if request === "hderivatives"
        # return for stabilityDerivatives()
        return h_derivatives
    end
end

"""
    avl_normal_vector(ds, theta)

constructs a normal vector the way AVL does

# Arguments
- `ds`: line representing the leading edge
- `theta`: the incidence angle, taken as a rotation (+ by RH rule) about the
    surface's spanwise axis projected onto the Y-Z plane.

# Returns
- `ncp / norm(ncp)`: the normal vector used by AVL
"""
function avl_normal_vector(ds, theta)

    st, ct = sincos(theta)

    # bound vortex vector
    bhat = ds/norm(ds)

    # chordwise strip normal vector
    shat = [0, -ds[3], ds[2]]/sqrt(ds[2]^2+ds[3]^2)

    # camberline vector
    chat = [ct, -st*shat[2], -st*shat[3]]

    # normal vector perpindicular to camberline and bound vortex for entire chordwise strip
    ncp = cross(chat, ds)
    return ncp / norm(ncp) # normal vector used by AVL
end
