
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
