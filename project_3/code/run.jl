include("project_3.jl")
include("plots_default.jl")

# retrieve optimized airframe design and pull out specific characteristics
optimized_airframe = optimize_airframe()
b = optimized_airframe[1]
c = optimized_airframe[2]
taper = optimized_airframe[3]
CL = optimized_airframe[4]
v = optimized_airframe[5]
L = optimized_airframe[6]
l = optimized_airframe[7]

# lift coefficient distribution for ideal taper ratio 
cl_dist_plot = vortex_lattice([b, c], 1.0, v, l, "cldist", taper)[3]
savefig(cl_dist_plot, "cl_dist.pdf")

# set values for comparison against ideal airframe
aoa = -10:1:10
chord_compare = 0.6   # mean chord length comparison
length_compare = 1.0    # wing to tail length comparison
taper_compare = 0.5   # taper ratio comparison

# lift coefficient distribution for compared taper ratio
cl_dist_compare = vortex_lattice([b, c], 1.0, v, l, "cldist", taper_compare)[3]
savefig(cl_dist_compare, "cl_dist_compare.pdf")

# retrieve lift and drag coefficients for ideal airframe and comparison airframes
ideal_coeff = aoa_coefficients(b, c, v, l, taper)
chord_coeff = aoa_coefficients(b, chord_compare, v, l, taper)
length_coeff = aoa_coefficients(b, c, v, length_compare, taper)
taper_coeff = aoa_coefficients(b, c, v, l, taper_compare)

# plot lift coefficient and drag coefficient against angle of attack for different designs for comparison
plot(aoa, ideal_coeff[1], xlabel="Angle of Attack (degrees)", ylabel=L"C_\ell", label="ideal")
p1 = plot!(aoa, chord_coeff[1], label="alternate", legend=:topleft)
plot(aoa, ideal_coeff[2], ylabel=L"C_d", label="ideal")
p2 = plot!(aoa, chord_coeff[2], label="alternate", legend=:top)
p12 = plot(p1, p2, layout=(1,2))
savefig(p12, "chord_coeff.pdf")

plot(aoa, ideal_coeff[1], xlabel="Angle of Attack (degrees)", ylabel=L"C_\ell", label="ideal")
p3 = plot!(aoa, length_coeff[1], label="alternate", legend=:topleft)
plot(aoa, ideal_coeff[2], ylabel=L"C_d", label="ideal")
p4 = plot!(aoa, length_coeff[2], label="alternate", legend=:top)
p34 = plot(p3, p4, layout=(1,2))
savefig(p34, "wingtail_coeff.pdf")

plot(aoa, ideal_coeff[1], xlabel="Angle of Attack (degrees)", ylabel=L"C_\ell", label="ideal")
p5 = plot!(aoa, taper_coeff[1], label="alternate", legend=:topleft)
plot(aoa, ideal_coeff[2], ylabel=L"C_d", label="ideal")
p6 = plot!(aoa, taper_coeff[2], label="alternate", legend=:top)
p56 = plot(p5, p6, layout=(1,2))
savefig(p56, "taper_coeff.pdf")

# plot ideal wing
plot([-b/2.0, b/2.0], [l, l], color=:blue, label="")
plot!([-b/2.0, -b/2.0], [l, l-c], color=:blue, label="")
plot!([b/2.0, b/2.0], [l, l-c], color=:blue, label="")
plot!([-b/2.0, b/2.0], [l-c, l-c], color=:blue, label="Wing")

# plot tail
plot!([-0.25, 0.25], [0, 0], color=:red, label="")
plot!([-0.25, -0.25], [0, -0.2], color=:red, label="")
plot!([0.25, 0.25], [0, -0.2], color=:red, label="")
p7 = plot!([-0.25, 0.25], [-0.2, -0.2], color=:red, label="Tail", legend=:bottomright)
savefig(p7, "ideal_design.pdf")

# plot alternate design 1
plot([-b/2.0, b/2.0], [l, l], color=:blue, label="")
plot!([-b/2.0, -b/2.0], [l, l-chord_compare], color=:blue, label="")
plot!([b/2.0, b/2.0], [l, l-chord_compare], color=:blue, label="")
plot!([-b/2.0, b/2.0], [l-chord_compare, l-chord_compare[1]], color=:blue, label="Wing")

# plot tail
plot!([-0.25, 0.25], [0, 0], color=:red, label="")
plot!([-0.25, -0.25], [0, -0.2], color=:red, label="")
plot!([0.25, 0.25], [0, -0.2], color=:red, label="")
p8 = plot!([-0.25, 0.25], [-0.2, -0.2], color=:red, label="Tail", legend=:bottomright)
savefig(p8, "chord_design.pdf")

# plot alternate design 2
plot([-b/2.0, b/2.0], [length_compare, length_compare], color=:blue, label="")
plot!([-b/2.0, -b/2.0], [length_compare, length_compare-c], color=:blue, label="")
plot!([b/2.0, b/2.0], [length_compare, length_compare-c], color=:blue, label="")
plot!([-b/2.0, b/2.0], [length_compare-c, length_compare-c], color=:blue, label="Wing")

# plot tail
plot!([-0.25, 0.25], [0, 0], color=:red, label="")
plot!([-0.25, -0.25], [0, -0.2], color=:red, label="")
plot!([0.25, 0.25], [0, -0.2], color=:red, label="")
p9 = plot!([-0.25, 0.25], [-0.2, -0.2], color=:red, label="Tail", legend=:bottomright)
savefig(p9, "length_design.pdf")

# plot alternate design 3
plot([-b/2.0, 0], [l-(c-c*taper_compare)/2, l], color=:blue, label="")
plot!([0, b/2.0], [l, l-(c-c*taper_compare)/2], color=:blue, label="")
plot!([-b/2.0, -b/2.0], [l-(c-c*taper_compare)/2, l-(c-c*taper_compare)/2-c*taper_compare], color=:blue, label="")
plot!([b/2.0, b/2.0], [l-(c-c*taper_compare)/2, l-(c-c*taper_compare)/2-c*taper_compare], color=:blue, label="")
plot!([-b/2.0, 0], [l-(c-c*taper_compare)/2-c*taper_compare, l-(c-c*taper_compare)-c*taper_compare], color=:blue, label="")
plot!([0, b/2.0], [l-(c-c*taper_compare)-c*taper_compare, l-(c-c*taper_compare)/2-c*taper_compare], color=:blue, label="Wing")

# plot tail
plot!([-0.25, 0.25], [0, 0], color=:red, label="")
plot!([-0.25, -0.25], [0, -0.2], color=:red, label="")
plot!([0.25, 0.25], [0, -0.2], color=:red, label="")
p10 = plot!([-0.25, 0.25], [-0.2, -0.2], color=:red, label="Tail", legend=:bottomright)
savefig(p10, "taper_design.pdf")
