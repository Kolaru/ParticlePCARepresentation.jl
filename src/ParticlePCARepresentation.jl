module ParticlePCARepresentation

using LinearAlgebra
using Makie
using MultivariateStats

import MultivariateStats: fit, projection, predict, reconstruct

export ParticlePCA
export fit, projection, predict, reconstruct
export plot_component3D
export plot_component2D, plot_component2D!
export animate_component2D!

struct ParticlePCA{T}
    model::PCA{T}
    ndims::Int
    nparts::Int
end

function ParticlePCA(data::Array{T, 3}) where T
    ndims, nparts, nobs = size(data)
    ndofs = ndims * nparts
    pca_data = reshape(data, ndofs, nobs)
    model = fit(PCA, pca_data ; pratio=1.0)
    return ParticlePCA(model, ndims, nparts)
end


fit(::Type{<:ParticlePCA}, data) = ParticlePCA(data)

function projection(ppca::ParticlePCA)
    pca_proj = projection(ppca.model)
    return reshape(pca_proj, ppca.ndims, ppca.nparts, :)
end

function predict(ppca::ParticlePCA, particles::Matrix)
    flat = reshape(particles, ppca.ndims * ppca.nparts)
    return predict(ppca.model, flat)
end

function predict(ppca::ParticlePCA, particles::Array{<:Any, 3})
    flat = reshape(particles, ppca.ndims * ppca.nparts, :)
    return predict(ppca.model, flat)
end

function reconstruct(ppca::ParticlePCA, projections::Vector)
    recons = reconstruct(ppca.model, projections)
    return reshape(recons, ppca.ndism, ppca.nparts)
end

function reconstruct(ppca::ParticlePCA, projections::Matrix)
    recons = reconstruct(ppca.model, projections)
    return reshape(recons, ppca.ndism, ppca.nparts, :)
end

function plot_component3D(
        ppca, positions ;
        component = 1,
        flip = false,
        kwargs...)

    comp = projection(ppca)[:, :, component]
    comp = flip ? -comp : comp
    arrow_starts = positions - comp./2

    fig, ax = arrows(
        arrow_starts[1, :], arrow_starts[2, :], arrow_starts[3, :],
        comp[1, :], comp[2, :], comp[3, :] ;
        align = :origin,
        kwargs...
    )

    meshscatter!(
        ax,
        positions[1, :], positions[2, :], positions[3, :] ;
        kwargs...
    )

    return fig, ax
end

function plot_component2D!(
        ax, ppca, positions ;
        component = 1,
        componentscale = 1,
        flip = false,
        xaxis = 1,
        yaxis = 2,
        scalealpha = false,
        color = :black,
        kwargs...)

    comp = projection(ppca)[:, :, component] * componentscale
    comp = flip ? -comp : comp
    arrow_starts = positions - comp./2

    if color isa Vector && scalealpha == true
        norms = norm.(eachcol(comp))
        color = map(zip(norms, color)) do (v, c)
            (c, 0.2 + v/maximum(norms) * 0.8)
        end
    end

    arrows!(
        ax,
        arrow_starts[xaxis, :], arrow_starts[yaxis, :],
        comp[xaxis, :], comp[yaxis, :] ;
        align = :origin,
        color,
        kwargs...
    )

    scatter!(
        ax,
        positions[xaxis, :], positions[yaxis, :] ;
        color,
        kwargs...
    )
end

function plot_component2D(ppca, positions ; kwargs...)
    fig = Figure()
    ax = fig[1, 1] = Axis(fig)
    plot_component2D!(ax, ppca, positions ; kwargs...)
    return fig, ax
end

function animate_component2D!(
        ax, ppca, positions ;
        amplitude = 1,
        component = 1,
        xaxis = 1,
        yaxis = 2,
        rate = 0.5,  # Half a cycle per second
        kwargs...)

    comp = projection(ppca)[:, :, component]
    w = Observable(0.0)
    xx = @lift(positions[xaxis, :] + $w * amplitude * comp[xaxis, :])
    yy = @lift(positions[yaxis, :] + $w * amplitude * comp[yaxis, :])

    scatter!(ax, xx, yy ; kwargs...)

    return function frame(F ; framerate = 24)
        t = F/framerate
        θ = t * rate * 2π 
        w[] = sin(θ)
    end
end

end