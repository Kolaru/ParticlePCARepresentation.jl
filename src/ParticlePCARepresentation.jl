module ParticlePCARepresentation

using Chain
using Colors
using ColorSchemes
using Distributions
using LinearAlgebra
using Makie
using MultivariateStats
using StatsBase

import MultivariateStats: fit, mean, projection, predict, reconstruct, principalvars

export ParticlePCA
export fit, projection, predict, reconstruct, principalvars
export sample
export blob_model, plot_blobs, plot_blobs!
export plot_component3D, plot_component3D!
export animate_component3D!, animateable_component3D!
export plot_component2D, plot_component2D!, summarize
export component3D!, trace_component3D!
export animate_component2D!

include("ppca.jl")
include("plot_2D.jl")
include("plot_3D.jl")

end