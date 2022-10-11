using VortexLattice
include("plots_default.jl")

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
    ratios = 0.003:0.001:0.008
    n_r = length(ratios)
    derivatives = zeros(n_r, 3)
    index = 1
    for r in ratios
        derivatives[index, :] .= vortexLattice("derivatives", 7.5, r)
        index = index + 1
    end

    d = plot(ratios, derivatives[:, 1], xlabel="Tail Volume Ratios", ylabel=L"\frac{d}{dx}C_d", label="")
    l = plot(ratios, derivatives[:, 3], xlabel="Tail Volume Ratios", ylabel=L"\frac{d}{dx}C_\ell", label="")
    plot(d, l)
end

function vortexLattice(request, ar=7.5, v_v=0.004978)
    coefficients = zeros(1, 2)

    # wing
    xle = [0.0, 0.4]
    yle = [0.0, 7.5]
    zle = [0.0, 0.0]
    chord = [2.2, 1.8]
    theta = [0.0, 0.0]
    phi = [0.0, 0.0]
    fc = fill((xc) -> 0, 2) # camberline function for each section

    # horizontal stabilizer
    xle_h = [0.0, 0.4]
    yle_h = [0.0, 1.25]
    zle_h = [0.0, 0.0]
    chord_h = [0.7, 0.42]
    theta_h = [0.0, 0.0]
    phi_h = [0.0, 0.0]
    fc_h = fill((xc) -> 0, 2) # camberline function for each section
    ns_h = 6
    nc_h = 3
    spacing_s_h = Uniform()
    spacing_c_h = Uniform()
    mirror_h = false

    # vertical stabilizer
    xle_v = [0.0, 0.4]
    yle_v = [0.0, 0.0]
    zle_v = [0.0, 1.0]
    chord_v = [0.7, 0.42]
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
    alpha = 1.0*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    fs = Freestream(Vinf, alpha, beta, Omega)

    # tail volume ratio variables
    S_v = 0.56
    l_v = (v_v*Sref*bref)/S_v

    # construct surface
    grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
        fc = fc, spacing_s=spacing_s, spacing_c=spacing_c)

    # generate surface panels for horizontal tail
    hgrid, htail = wing_to_surface_panels(xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h;
        mirror=mirror_h, fc=fc_h, spacing_s=spacing_s_h, spacing_c=spacing_c_h)
    VortexLattice.translate!(hgrid, [l_v, 0.0, 0.0])
    VortexLattice.translate!(htail, [l_v, 0.0, 0.0])

    # generate surface panels for vertical tail
    vgrid, vtail = wing_to_surface_panels(xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v;
        mirror=mirror_v, fc=fc_v, spacing_s=spacing_s_v, spacing_c=spacing_c_v)
    VortexLattice.translate!(vgrid, [l_v, 0.0, 0.0])
    VortexLattice.translate!(vtail, [l_v, 0.0, 0.0])

    # create vector containing all surfaces
    if request === "coefficients"
        surfaces = [surface]
    end
    if request === "derivatives"
        surfaces = [surface, htail, vtail]
    end

    # we can use symmetry since the geometry and flow conditions are symmetric about the X-Z axis
    symmetric = true

    # perform steady state analysis
    system = steady_analysis(surfaces, ref, fs; symmetric=symmetric)

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

    derivatives = [CDa, CYa, CLa]
    println(derivatives)

    if request === "coefficients"
        return coefficients
    end
    if request === "derivatives"
        return derivatives
    end

end
