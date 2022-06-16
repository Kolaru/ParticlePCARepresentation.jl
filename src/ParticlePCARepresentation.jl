module ParticlePCARepresentation

using Makie
using MultivariateStats

import MultivariateStats: fit, projection

export ParticlePCA
export fit, projection
export plot_component2D!

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

function plot_component2D!(
        ax, ppca, positions ;
        component = 1,
        xaxis = 1,
        yaxis = 2,
        kwargs...)

    comp = projection(ppca)[:, :, component]

    scatter!(
        ax,
        positions[xaxis, :], positions[yaxis, :],
        kwargs...
    )

    arrows!(
        ax,
        positions[xaxis, :], positions[yaxis, :],
        comp[xaxis, :], comp[yaxis, :] ;
        align = :center,  # TODO Seems to not to what it should
        kwargs...
    )
end

end