
function plot_component2D!(
        ax, ppca ;
        positions = mean(ppca),
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

function plot_component2D(ppca ; kwargs...)
    fig = Figure()
    ax = fig[1, 1] = Axis(fig)
    plot_component2D!(ax, ppca ; kwargs...)
    return fig, ax
end

function plot_blobs!(
        ax, ppca ;
        colors = :black,
        xaxis = 1,
        yaxis = 2,
        cmap = nothing)

    if !(colors isa Vector)
        colors = Iterators.repeated(colors)
    end

    for (μ, Σ, color) in zip(eachcol(mean(ppca)), diag(covariance(ppca)), colors)
        μ = μ[[xaxis, yaxis]]
        Σ = Σ[[xaxis, yaxis], [xaxis, yaxis]]

        invΣ = inv(Σ)
        gauss(u) = exp( -0.5 * dot(u - μ, invΣ, u - μ) )

        if isnothing(cmap)
            c = parse(Colorant, color)
            cmap = ColorScheme(range(RGBA(RGB(c), 0), c, length=10))
        end
        
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


function summarize(gridpos, ppca ;
        components = 1:6,
        positions = mean(ppca),
        labels = nothing,
        kwargs...)
    layout = GridLayout(gridpos, 2, length(components))
    axes = [Axis(layout[i, j]) for i in 3:4, j in eachindex(components)]

    rowgap!(layout, 1, 0)
    hidexdecorations!.(axes[1, :], grid = false, ticks = false)
    hideydecorations!.(axes[:, 2:end], grid = false, ticks = false)

    pad = maximum(abs.(positions))/2
    xlims = (minimum(positions[1, :]) - pad, maximum(positions[1, :]) + pad)
    ylims = (minimum(positions[2, :]) - pad, maximum(positions[2, :]) + pad)
    zlims = (minimum(positions[3, :]) - pad, maximum(positions[3, :]) + pad)

    for (j, component) in enumerate(components)
        label = isnothing(labels) ? "Component $component" : labels[j]
        layout[1, j] = Label(gridpos.layout.parent, label, tellwidth = false)
        evar = round(ppca.model.prinvars[component] / ppca.model.tvar ; sigdigits = 2)
        layout[2, j] = Label(
            gridpos.layout.parent,
            "EVR = $evar",
            tellwidth = false)

        xlims!(axes[1, j], xlims)
        xlims!(axes[2, j], xlims)
        ylims!(axes[1, j], ylims)
        ylims!(axes[2, j], zlims)
        plot_component2D!(axes[1, j], ppca ; component, positions, yaxis = 2, kwargs...)
        plot_component2D!(axes[2, j], ppca ;  component, positions, yaxis = 3, kwargs...)
    end

    return layout, axes
end