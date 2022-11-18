module ParticlePCARepresentation

using Colors
using ColorSchemes
using LinearAlgebra
using Makie
using MultivariateStats

import MultivariateStats: fit, mean, projection, predict, reconstruct

export ParticlePCA
export fit, projection, predict, reconstruct
export plot_blobs, plot_blobs!
export plot_component3D
export plot_component2D, plot_component2D!
export animate_component2D!

struct ParticlePCA{T}
    model::PCA{T}
    ndims::Int
    nparts::Int
end

"""
    ParticlePCA(data::Array{<:Any, 3})

Compute the PCA for data from a particle system.

`data` must have the shape `n_dimensions x n_particles x n_obs`.
"""
function ParticlePCA(data::Array{T, 3}) where T
    ndims, nparts, nobs = size(data)
    ndofs = ndims * nparts
    pca_data = reshape(data, ndofs, nobs)
    model = fit(PCA, pca_data ; pratio = 1.0)
    return ParticlePCA(model, ndims, nparts)
end

"""
    ParticlePCA(Σ::Matrix, μ::Vector ; ndims = 3)

Create a particle PCA representation of data with covariance matrix Σ and
mean μ.
"""
function ParticlePCA(Σ::Matrix, μ::Vector ; ndims = 3)
    nparts = div(length(μ), ndims)
    model = pcacov(Σ, μ ; pratio = 1.0)
    return ParticlePCA(model, ndims, nparts)
end

function ParticlePCA(proj::Matrix, μ::Vector, weights::Vector ; ndims = 3)
    nparts = div(length(μ), ndims)
    model = PCA(μ, proj, weights, sum(weights))
    return ParticlePCA(model, ndims, nparts)
end

fit(::Type{<:ParticlePCA}, data) = ParticlePCA(data)

function covariance(pca::PCA)
    return projection(pca) * Diagonal(principalvars(pca)) * projection(pca)'
end

function covariance(ppca::ParticlePCA)
    Σ = covariance(ppca.model)
    s = reshape(Σ, ppca.ndims, ppca.nparts, ppca.ndims, ppca.nparts)
    return [s[:, i, :, j] for i in 1:ppca.nparts, j in 1:ppca.nparts]
end

mean(ppca::ParticlePCA) = reshape(mean(ppca.model), ppca.ndims, :)

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
        ppca, positions = mean(ppca) ;
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
        ax, ppca, positions = mean(ppca)  ;
        component = 1,
        componentscale = 0.5,
        flip = false,
        xaxis = 1,
        yaxis = 2,
        scalealpha = true,
        color = :black,
        kwargs...)

    ax.aspect = DataAspect()

    rescale = maximum(abs.(positions)) * componentscale

    comp = projection(ppca)[:, :, component]
    comp = comp * rescale / maximum(norm.(eachcol(comp)))
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

function plot_component2D(ppca, positions = mean(ppca) ; kwargs...)
    fig = Figure()
    ax = fig[1, 1] = Axis(fig)
    plot_component2D!(ax, ppca, positions ; kwargs...)
    return fig, ax
end

function plot_blobs!(
        ax, ppca ;
        colors = :black,
        xaxis = 1,
        yaxis = 2)
    
    if !(colors isa Vector)
        colors = Iterators.repeated(colors)
    end

    for (μ, Σ, color) in zip(eachcol(mean(ppca)), diag(covariance(ppca)), colors)
        μ = μ[[xaxis, yaxis]]
        Σ = Σ[[xaxis, yaxis], [xaxis, yaxis]]

        invΣ = inv(Σ)
        gauss(u) = exp( -0.5 * dot(u - μ, invΣ, u - μ) )

        c = parse(Colorant, color)
        cmap = ColorScheme(range(RGBA(RGB(c), 0), c, length=10))
        rad = sqrt(maximum(Σ)) * 2
        xs = μ[1] .+ range(-rad, rad, length=51)
        ys = μ[2] .+ range(-rad, rad, length=51)
        zs = [gauss([x, y]) for x in xs, y in ys]

        contourf!(
            ax, xs, ys, zs ;
            colormap = cmap,
            levels = range(0, 1, length=10))

        scatter!(ax, μ[1], μ[2] ; color=:red)
    end
end

function plot_blobs(ppca ; kwargs...)
    fig = Figure()
    ax = fig[1, 1] = Axis(fig)
    plot_blobs!(ax, ppca ; kwargs...)
    return fig, ax
end

function animate_component2D!(
        ax, ppca, positions = mean(ppca) ;
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