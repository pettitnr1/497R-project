using VortexLattice
using LinearAlgebra

include("plots_default.jl")

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
"""
function vortex_lattice(request="", ar=2, v=1, aoa=1.0*pi/180)
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
    mirror = true # account for geometric symmetry here because flow is not symmetric

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

    if (request === "ratios")
        # tail volume ratio variables
        S_v = zle_v[2]*chord_v[1]
        S_h = yle_h[2]*chord_h[1]
        l_v = (v*Sref*bref)/S_v
        l_h = l_v
    else
        l_v = 1
        l_h = l_v
    end
    println(l_v)

    # reference parameters
    bref = 1.5 # span
    if (request === "ratios")
        Sref = (bref^2 * S_h)/S_v # area
    else
        Sref = bref^2 / ar # area
    end
    cref = Sref/bref # chord
    rref = [0.50, 0.0, 0.0]
    Vinf = 1.0
    ref = Reference(Sref, cref, bref, rref, Vinf)

    # freestream parameters
    alpha = aoa*pi/180
    beta = 1.0*pi/180
    Omega = [0.0; 0.0; 0.0]
    fs = Freestream(Vinf, alpha, beta, Omega)

    # construct surface
    grid, wing = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, fc = fc, spacing_s=spacing_s, spacing_c=spacing_c)

    # now set normal vectors manually
    ncp = avl_normal_vector([xle[2]-xle[1], yle[2]-yle[1], zle[2]-zle[1]], 2.0*pi/180)

    # overwrite normal vector for each wing panel
    for i = 1:length(wing)
        wing[i] = set_normal(wing[i], ncp)
    end

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
    println(bref)
    println(cref)
    println(Sref)
    println(CL)
    L = 0.5*CL*1.225*Vinf^2*Sref

    if (request === "ratios")
        return v_derivatives
    else
        return L
    end
end

"""
    avl_normal_vector(ds, theta)

constructs a normal vector the way AVL does

# Arguments
- `ds`: line representing the leading edge
- `theta`: the incidence angle, taken as a rotation (+ by RH rule) about the
    surface's spanwise axis projected onto the Y-Z plane.
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
