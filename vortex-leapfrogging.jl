"""
FLOW Lab Undergraduate Onboarding: Introductory Project
Vortex Leapfrogging
"""

using DataFrames
using LinearAlgebra
using Plots

function inducedVelocity(vortex = vector, strength = Vector, vortexes = Vector[])
    total = [0, 0, 0]
    for i in vortexes
        distanceVector = vortex - i
        if i[2] >= 0
            single = (cross(strength, distanceVector)) / (2 * pi * (magnitude(distanceVector))^2)
        else
            single = -((cross(strength, distanceVector)) / (2 * pi * (magnitude(distanceVector))^2))
        end
        total = total + single
    end
    return total
end

function magnitude(vector = Vector)
    return sqrt((vector[1])^2 + (vector[2])^2 + (vector[3])^2)
end

function leapFroggingVertexes(d, strength, steps)
    line1x = []
    line1y = []
    line2x = []
    line2y = []
    line3x = []
    line3y = []
    line4x = []
    line4y = []
    p1i = [0, -d/2, 0]
    p2i = [0, d/2, 0]
    p3i = [d, d/2, 0]
    p4i = [d, -d/2, 0]
    push!(line1x, p1i[1])
    push!(line1y, p1i[2])
    push!(line2x, p2i[1])
    push!(line2y, p2i[2])
    push!(line3x, p3i[1])
    push!(line3y, p3i[2])
    push!(line4x, p4i[1])
    push!(line4y, p4i[2])
    i = 1
    t = .01
    while i < steps
        p1f = p1i + (t * inducedVelocity(p1i, strength, [p2i, p3i, p4i]))
        p2f = p2i + (t * inducedVelocity(p2i, strength, [p1i, p3i, p4i]))
        p3f = p3i + (t * inducedVelocity(p3i, strength, [p1i, p2i, p4i]))
        p4f = p4i + (t * inducedVelocity(p4i, strength, [p1i, p2i, p3i]))
        push!(line1x, p1f[1])
        push!(line1y, p1f[2])
        push!(line2x, p2f[1])
        push!(line2y, p2f[2])
        push!(line3x, p3f[1])
        push!(line3y, p3f[2])
        push!(line4x, p4f[1])
        push!(line4y, p4f[2])
        p1i = p1f
        p2i = p2f
        p3i = p3f
        p4i = p4f
        i = i + 1
    end

    plot!(plot!(plot!(plot(line1x, line1y, color = :blue), line2x, line2y, color = :blue), line3x, line3y, color = :red), line4x, line4y, color = :red)
end
