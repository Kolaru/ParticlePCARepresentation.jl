using ParticlePCARepresentation
using GLMakie

data = rand(3, 3, 100)
ppca = ParticlePCA(data)
positions = Float64[
    -1 0 1 ;
    1 0 1 ;
    0 0 0
]

begin
    fig, ax = plot_component2D(
        ppca, positions ;
        color=[:red, :blue, :red])

    ax.xlabel = "x (a.u.)"
    ax.ylabel = "y (a.u.)"

    display(fig)
end
