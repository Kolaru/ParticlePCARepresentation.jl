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
export plot_component2D, plot_component2D!, summarize
export animate_component2D!

include("ppca.jl")
include("plot_2D.jl")
include("plot_3D.jl")

end