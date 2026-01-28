### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 4ca07304-fbe4-11f0-a542-29271494572d
begin
    using Pkg
    Pkg.activate(@__DIR__)
    using NSFD
    using CairoMakie
    CairoMakie.activate!()
end

# ╔═╡ 99fa8e64-3f6b-4675-8265-cf1c077098a9
begin
    x_length = 1.0        # domain size in x-direction
    y_length = 1.0        # domain size in y-direction
    imax = 128            # number of interior cells in x-direction
    jmax = 128            # number of interior cells in y-direction

    geom = Geometry(x_length, y_length, imax, jmax)
end

# ╔═╡ e5e6b3f5-a626-4470-ae13-50e78376a6ff
begin
    t_end = 20.0  # final time
    τ = 0.5       # safety factor for time step size control

    iter_max = 100  # maximal number of pressure iterations in one time step
    ϵₚ = 0.001      # stopping tolerance for pressure iteration
    ω = 1.7         # relaxation parameter for SOR iteration
    γ = 0.9         # upwind differencing factor

    Re = 1000.0  # Reynolds number
    GX = 0.0     # body force in x-direction
    GY = 0.0     # body force in y-direction
    UI = 0.0     # initial velocity in x-direction
    VI = 0.0     # initial velocity in y-direction
    PI = 0.0     # initial pressure
end

# ╔═╡ 434b7f22-83e5-447d-8885-2f408c5ebbfc
begin
    uvp = FlowFields(imax, jmax)

    set!(uvp, UI, VI, PI)

    bc = DomainBC(NoSlip(Domain, North, 1.0),
                  NoSlip(Domain, South),
                  NoSlip(Domain, East),
                  NoSlip(Domain, West))

    fg = FG(imax, jmax)

    # right-hand side for pressure iteration
    rhs = StaggeredField(CC, imax, jmax)
end

# ╔═╡ 43083c3d-cb7c-4526-99ee-15241cb2228a
iter_pressure = PressureIterator(geom, ω, ϵₚ, iter_max)

# ╔═╡ 58b68f61-96cc-4e79-b6a7-1d698da848bf
begin
    t = 0  # set initial time to zero
    while t < t_end
        δt = compute_δt(uvp, geom, Re, τ)
        apply!(bc, uvp, geom)
        fill_fg!(fg, uvp, geom, δt, GX, GY, Re, γ)
        fill_rhs!(rhs, fg, geom, δt)
        it = iter_pressure(uvp, rhs, geom)
        if any(isnan, uvp.p)
            error("NaN in p at t = $t")
        end
        fill_uv_next!(uvp, fg, geom, δt)
        if any(isnan, uvp.u) | any(isnan, uvp.v)
            error("NaN in velocities")
        end
        global t = t + δt
    end
end

# ╔═╡ 33a557f0-df02-45b2-8fb1-fa137291ab9c
begin
    nx_interp = 10
    ny_interp = 10

    x = range(x_length / nx_interp / 2,
              x_length - x_length / nx_interp / 2,
              nx_interp)
    y = range(y_length / ny_interp / 2,
              y_length - y_length / ny_interp / 2,
              ny_interp)

    # meshgrid
    X = repeat(x, 1, length(y))
    Y = repeat(y', length(x), 1)

    U_interp = Matrix{Float64}(undef, nx_interp, ny_interp)
    V_interp = Matrix{Float64}(undef, nx_interp, ny_interp)

    for j in 1:ny_interp
        for i in 1:nx_interp
            U_interp[i, j] = interp(uvp.u, geom.δx, geom.δy, X[i, j], Y[i, j])
            V_interp[i, j] = interp(uvp.v, geom.δx, geom.δy, X[i, j], Y[i, j])
        end
    end

    # magnitude of the interpolated vectors
    mag = sqrt.(U_interp .^ 2 + V_interp .^ 2)

    # maximum magnitude of the vectors
    max_mag = maximum(mag)

    scale = min(step(x), step(y)) / max_mag

    fig = Figure()
    ax = Axis(fig[1, 1]; aspect=DataAspect())
    arrows2d!(ax,
              vec(X), vec(Y),
              vec(scale .* U_interp), vec(scale .* V_interp))
    fig
end

# ╔═╡ Cell order:
# ╠═4ca07304-fbe4-11f0-a542-29271494572d
# ╠═99fa8e64-3f6b-4675-8265-cf1c077098a9
# ╠═e5e6b3f5-a626-4470-ae13-50e78376a6ff
# ╠═434b7f22-83e5-447d-8885-2f408c5ebbfc
# ╠═43083c3d-cb7c-4526-99ee-15241cb2228a
# ╠═58b68f61-96cc-4e79-b6a7-1d698da848bf
# ╠═33a557f0-df02-45b2-8fb1-fa137291ab9c
