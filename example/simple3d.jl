using ParticlePCARepresentation
using GLMakie

data = rand(3, 5, 100)
ppca = ParticlePCA(data)
positions = rand(3, 5)

begin
    fig, ax = plot_component3D(
        ppca, positions ;
        color=[:red, :red, :blue, :blue, :black])

    display(fig)
end