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

function covariances(ppca::ParticlePCA)
    Σ = covariance(ppca.model)
    s = reshape(Σ, ppca.ndims, ppca.nparts, ppca.ndims, ppca.nparts)
    return [s[:, i, :, j] for i in 1:ppca.nparts, j in 1:ppca.nparts]
end

function means(ppca::ParticlePCA)
    @chain ppca.model begin
        mean
        reshape(_, ppca.ndims, :)
        eachcol
        collect
    end
end

variances(ppca::ParticlePCA) = diag(covariances(ppca))

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

function blob_model(
        ppca::ParticlePCA ;
        dims = [1, 2],
        atoms = 1:ppca.nparts)

    μs = @chain ppca begin
        means
        _[atoms]
        map(_) do μ
            μ[dims]
        end
    end
    Σs = @chain ppca begin
        variances
        _[atoms]
        map(_) do Σ
            Symmetric(Σ[dims, dims])
        end
    end

    return MixtureModel(MultivariateNormal, collect(zip(μs, Σs)))
end