#=
FLOW Lab Undergraduate Onboarding: Introductory Project
Leapfrogging Vortex Rings
=#

using LinearAlgebra
using Plots

"""
    inducedVelocity(vortex, strength, vortexes)

Compute the induced velocity on a vortex
    vortex- position of vortex that induced velocity will be calculated for
    strength- strength of vortexes
    vortexes- positions of other vortexes
"""
function inducedVelocity(vortex, strength, vortexes)
    total = [0, 0, 0]
    for v in vortexes
        distanceVector = vortex - v
        if v[2] >= 0
            single = (LinearAlgebra.cross(strength, distanceVector)) / (2.0 * pi * (LinearAlgebra.norm(distanceVector))^2)
        else
            single = -((LinearAlgebra.cross(strength, distanceVector)) / (2.0 * pi * (LinearAlgebra.norm(distanceVector))^2))
        end
        total = total + single
    end
    return total
end

"""
    leapFroggingVortexes(d, strength, steps, time)

Computes plots for leapfrogging vortex rings and graphs them
    distance- distance between vortexes
    strength- strength of vortexes
    steps- how many times to calculate position
    time- time interval between
"""
function leapFroggingVortexes(distance, strength, steps, time)
    line1 = []
    line2 = []
    line3 = []
    line4 = []
    p1i = [0, -distance/2.0, 0]
    p2i = [0, distance/2.0, 0]
    p3i = [distance, distance/2.0, 0]
    p4i = [distance, -distance/2.0, 0]
    push!(line1, [p1i[1], p1i[2]])
    push!(line2, [p2i[1], p2i[2]])
    push!(line3, [p3i[1], p3i[2]])
    push!(line4, [p4i[1], p4i[2]])
    for step in 1:steps
        p1f = p1i + (time * inducedVelocity(p1i, strength, [p2i, p3i, p4i]))
        p2f = p2i + (time * inducedVelocity(p2i, strength, [p1i, p3i, p4i]))
        p3f = p3i + (time * inducedVelocity(p3i, strength, [p1i, p2i, p4i]))
        p4f = p4i + (time * inducedVelocity(p4i, strength, [p1i, p2i, p3i]))
        push!(line1, [p1f[1], p1f[2]])
        push!(line2, [p2f[1], p2f[2]])
        push!(line3, [p3f[1], p3f[2]])
        push!(line4, [p4f[1], p4f[2]])
        p1i = p1f
        p2i = p2f
        p3i = p3f
        p4i = p4f
    end
    plot(map((point) -> point[1], line1), map((point) -> point[2], line1), color = :blue, label = "Vortex 1", grid = false, xlabel = "Horizontal Position (meters)", ylabel = "Vertical Position (meters)", xlims = (0,12))
    plot!(map((point) -> point[1], line2), map((point) -> point[2], line2), color = :blue, label = "")
    plot!(map((point) -> point[1], line3), map((point) -> point[2], line3), color = :red, label = "Vortex 2")
    plot!(map((point) -> point[1], line4), map((point) -> point[2], line4), color = :red, label = "")
end
