module ParticlePCARepresentation

using Makie
using MultivariateStats

import MultivariateStats: fit, projection

export ParticlePCA
export fit, projection
export plot_component3D
export plot_component2D, plot_component2D!

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
        flip = false,
        xaxis = 1,
        yaxis = 2,
        kwargs...)

    comp = projection(ppca)[:, :, component]
    comp = flip ? -comp : comp
    arrow_starts = positions - comp./2

    scatter!(
        ax,
        positions[xaxis, :], positions[yaxis, :] ;
        kwargs...
    )

    arrows!(
        ax,
        arrow_starts[xaxis, :], arrow_starts[yaxis, :],
        comp[xaxis, :], comp[yaxis, :] ;
        align = :origin,
        kwargs...
    )
end

function plot_component2D(ppca, positions ; kwargs...)
    fig = Figure()
    ax = fig[1, 1] = Axis(fig)
    plot_component2D!(ax, ppca, positions ; kwargs...)
    return fig, ax
end

end