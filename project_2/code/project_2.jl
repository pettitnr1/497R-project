using VortexLattice
using LinearAlgebra

include("plots_default.jl")

function aoaEffect()
    angles = -20:1:20
    n_a = length(angles)
    lift = zeros(n_a)
    index = 1
    for a in angles
        cl = vortexLattice("coefficients", 7.5, 1, a)
        lift[index] = cl[1]
        index = index + 1
    end
    plot(angles, lift, xlabel="Angle of Attack(degrees)", ylabel=L"C_\ell", label="")
end

function wingEfficiency()
    ratios = 3:1:15
    n_r = length(ratios)
    coefficients = zeros(n_r, 2)
    efficiency = 0
    e = Vector{Float64}(undef, n_r)
    index = 1
    for ar in ratios
        coefficients[index, :] .= vortexLattice("coefficients", ar)
        cl = coefficients[index, 1]
        cd = coefficients[index, 2]
        efficiency = (cl^2) / (pi * ar * cd)
        e[index] = efficiency
        index = index + 1
    end

    plot(ratios, e, xlabel="Aspect Ratio", ylabel="Efficiency", label="")
end

function stabilityDerivatives()
    v_ratios = range(0.00467, 0.0156, 20)
    h_ratios = range(0.04375, 0.1458, 20)
    v_r = length(v_ratios)
    h_r = length(h_ratios)
    v_derivatives = zeros(v_r, 2)
    h_derivatives = zeros(h_r, 2)
    index = 1
    for r in v_ratios
        v_derivatives[index, :] .= vortexLattice("vderivatives", 7.5, r)
        index = index + 1
    end
    index = 1
    for r in h_ratios
        h_derivatives[index, :] .= vortexLattice("hderivatives", 7.5, r)
        index = index + 1
    end

    lb = plot(v_ratios, v_derivatives[:, 1], xlabel="V-tail Volume Ratios", ylabel=L"C_{\ell{b}}", label="")
    nb = plot(v_ratios, v_derivatives[:, 2], xlabel="V-tail Volume Ratios", ylabel=L"C_{nb}", label="")
    # la = plot(h_ratios, h_derivatives[:, 1], xlabel="H-tail Volume Ratios", ylabel=L"C_{\ell{a}}", label="")
    # ma = plot(h_ratios, h_derivatives[:, 2], xlabel="H-tail Volume Ratios", ylabel=L"C_{ma}", label="")

    plot(lb, nb, layout=(2,1))
end

function vortexLattice(request, ar=7.5, v=1, aoa=1.0*pi/180)
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
        mirror = true
    else
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
        S_v = zle_v[2]*chord_v[1]
        l_v = (v*Sref*bref)/S_v
        l_h = l_v
    end

    if (request === "hderivatives")
        # horizontal tail volume ratio variables
        S_h = yle_h[2]*chord_h[1]
        l_h = (v*Sref*(chord[1]+chord[2])/2)/S_h
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
        VortexLattice.translate!(hgrid, [l_h, 0.0, 0.0])
        VortexLattice.translate!(htail, [l_h, 0.0, 0.0])

        # generate surface panels for vertical tail
        vgrid, vtail = wing_to_surface_panels(xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v;
            mirror=mirror_v, fc=fc_v, spacing_s=spacing_s_v, spacing_c=spacing_c_v)
        VortexLattice.translate!(vgrid, [l_v, 0.0, 0.0])
        VortexLattice.translate!(vtail, [l_v, 0.0, 0.0])

        # create vector containing all surfaces
        surfaces = [wing, htail, vtail]
        surface_id = [1, 2, 3]
    else
        surfaces = [wing]
        surface_id = [1]
    end

    # we can use symmetry since the geometry and flow conditions are symmetric about the X-Z axis
    if (request === "vderivatives" || request === "hderivatives")
        symmetric = false
    else
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
    h_derivatives = [Cla, Cma]

    if request === "coefficients"
        return coefficients
    end
    if request === "vderivatives"
        return v_derivatives
    end
    if request === "hderivatives"
        return h_derivatives
    end

end

# function to construct a normal vector the way AVL does
#  - `ds` is a line representing the leading edge
#  - `theta` is the incidence angle, taken as a rotation (+ by RH rule) about
#        the surface's spanwise axis projected onto the Y-Z plane.
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
