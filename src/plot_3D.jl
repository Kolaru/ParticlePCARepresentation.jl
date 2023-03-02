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

function plot_component3D!(
        ax, ppca ;
        component = 1,
        flip = false,
        positions = mean(ppca),
        kwargs...)

    comp = projection(ppca)[:, :, component]
    comp = flip ? -comp : comp

    diam = maximum(comp) - minimum(comp)
    positions ./= diam
    comp ./= diam
    
    arrow_starts = positions - comp./2


    arrows!(
        arrow_starts[1, :], arrow_starts[2, :], arrow_starts[3, :],
        comp[1, :], comp[2, :], comp[3, :] ;
        align = :origin,
        kwargs...
    )

    meshscatter!(
        ax,
        positions[1, :], positions[2, :], positions[3, :] ;
        markersize = diam / 20,
        kwargs...
    )

    return ax
end

function animate_component3D!(
        ax, ppca, positions = mean(ppca) ;
        amplitude = 1,
        component = 1,
        rate = 0.5,  # Half a cycle per second
        kwargs...)

    comp = projection(ppca)[:, :, component]
    w = Observable(0.0)
    xx = @lift(positions[1, :] + $w * amplitude * comp[1, :])
    yy = @lift(positions[2, :] + $w * amplitude * comp[2, :])
    zz = @lift(positions[3, :] + $w * amplitude * comp[3, :])

    scatter!(ax, xx, yy, zz ; kwargs...)

    return function frame(F ; framerate = 24)
        t = F/framerate
        θ = t * rate * 2π 
        w[] = sin(θ)
    end
end
