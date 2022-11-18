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
