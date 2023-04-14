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
        fig, ax, ppca,
        savepath, frames ;
        positions = mean(ppca),
        amplitude = 1,
        component = 1,
        rate = 0.5,  # Half a cycle per second
        framerate = 24,
        kwargs...)

    comp = projection(ppca)[:, :, component]
    w = Observable(0.0)
    xx = @lift(positions[1, :] + $w * amplitude * comp[1, :])
    yy = @lift(positions[2, :] + $w * amplitude * comp[2, :])
    zz = @lift(positions[3, :] + $w * amplitude * comp[3, :])

    meshscatter!(ax, xx, yy, zz ; kwargs...)
    
    record(fig, savepath, frames) do F
        t = F/framerate
        θ = t * rate * 2π 
        w[] = sin(θ)
    end
    
    nothing
end

function component3D!(
        ax, ppca, time ;
        positions = reshape(mean(ppca), 3, :),
        amplitude = 0.5*maximum(abs.(positions)),
        component = 1,
        cycle_freq = 0.5,  # Half a cycle per second
        kwargs...,
    )
    comp = projection(ppca)[:, :, component]
    w = @lift sin($time * cycle_freq * 2π)
    xx = @lift(positions[1, :] + $w * amplitude * comp[1, :])
    yy = @lift(positions[2, :] + $w * amplitude * comp[2, :])
    zz = @lift(positions[3, :] + $w * amplitude * comp[3, :])

    return meshscatter!(ax, xx, yy, zz ; kwargs...)
end

animateable_component3D!(args... ; kwargs...) = component3D!(args... ; kwargs...)

function trace_component3D!(
        ax, ppca ;
        nsteps = 20,
        markersize = 1,
        color = :black,
        kwargs...
    )

    if !hasmethod(iterate, Tuple{typeof(color)})
        color = [color]
    end

    for t in range(-0.5, 0.5 ; length = nsteps)
        cc = map(color) do c
            (c, 2t)
        end

        component3D!(
            ax, ppca, t; 
            color = color,
            markersize = (t + 0.5) * markersize,
            kwargs...
        )
    end
end