using GLMakie

include("nsfd.jl")

x_length = 1.0
y_length = 1.0
imax = 25
jmax = 25
δx = x_length / imax
δy = y_length / jmax

UI = 0.0
VI = 0.0
PI = 0.0

u = StaggeredField(U, imax, jmax)
v = StaggeredField(V, imax, jmax)
p = StaggeredField(P, imax, jmax)

f = StaggeredField(F, imax, jmax)
g = StaggeredField(G, imax, jmax)

rhs = StaggeredField(RHS, imax, jmax)

p_res = StaggeredField(PRes, imax, jmax)

# problem data
Re = 1.0
GX = 0.0
GY = 0.0
γ = 0.9

# time-stepping data
T = 20.0
τ = 0.5

# pressure iteration data
iter_max = 100
ϵₚ = 0.001
ω = 1.7
γ = 0.9

set!(u, UI)
set!(v, VI)
set!(p, PI)

t::Float64 = 0.0

# compute loop
while t < T
    δt = compute_δt(u, v, δx, δy, Re, τ)
    set_bv!(u, v, imax, jmax)
    fill_f!(f, u, v, δx, δy, δt, GX, Re, γ, imax, jmax)
    fill_g!(g, u, v, δx, δy, δt, GY, Re, γ, imax, jmax)
    set_bv!(u, v, f, g, imax, jmax)
    fill_rhs!(rhs, f, g, δx, δy, δt, imax, jmax)
    it = iter_pressure(p, rhs, p_res, δx, δy, ω, ϵₚ, iter_max, imax, jmax)
    fill_u_next!(u, f, p, δx, δt, imax, jmax)
    fill_v_next!(v, g, p, δy, δt, imax, jmax)
    global t = t + δt
end

# interpolate and plot the velocity

nx_interp = 10
ny_interp = 10

x = range(x_length / nx_interp / 2, x_length - x_length / nx_interp / 2, nx_interp)
y = range(y_length / ny_interp / 2, y_length - y_length / ny_interp / 2, ny_interp)

# meshgrid
X = repeat(x, 1, length(y))
Y = repeat(y', length(x), 1)

U_interp = Matrix{Float64}(undef, nx_interp, ny_interp)
V_interp = Matrix{Float64}(undef, nx_interp, ny_interp)

for j in 1:ny_interp
    for i in 1:nx_interp
        U_interp[i, j] = interpolate(u, δx, δy, X[i, j], Y[i, j])
        V_interp[i, j] = interpolate(v, δx, δy, X[i, j], Y[i, j])
    end
end

# magnitude of the interpolated vectors
mag = sqrt.(U_interp .^ 2 + V_interp .^ 2)

# maximum magnitude of the vectors
max_mag = maximum(mag)

# scale by interpolated grid step size and maximum magnitude
# one full grid diagonal is equal to the maximum magnitude
scale = sqrt(step(x)^2 + step(y)^2) / max_mag

f = Figure();
ax = Axis(f[1, 1]; aspect=DataAspect());
arrows2d!(ax, vec(X), vec(Y), vec(scale .* U_interp), vec(scale .* V_interp))
heatmap!(ax, 0:δx:x_length, 0:δy:y_length, p.values[2:(end - 1), 2:(end - 1)]);
f
