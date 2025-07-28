using NSFD, CairoMakie

# grid data
x_length = 1.0
y_length = 1.0
nx = 128
ny = 128
grid_data = GridData(x_length, y_length, nx, ny)

UI = 0.0
VI = 0.0
PI = 0.0

# problem data
Re = 1000.0
GX = 0.0
GY = 0.0
G = StaggeredVector(GX, GY)
problem_data = ProblemData(Re, G)

# time-stepping data
t = 0
T = 20.0
τ = 0.5
Δt = NSFD.compute_Δt(grid_data.Δx, grid_data.Δy, Re, τ)
time_data = TimeSteppingData(t, T, Δt, τ)

# pressure iteration data
iter_max = 100
ϵₚ = 0.001
ω = 1.7
γ = 0.9
p_iter_data = PressureIterData(iter_max, ϵₚ, ω, γ)

initial_state = State(grid_data, problem_data, time_data, p_iter_data)
set!(initial_state.u, (UI, VI))
set!(initial_state.p, PI)

global_bc = GlobalBCond(NNoSlipBCond(1.0), SNoSlipBCond(), ENoSlipBCond(), WNoSlipBCond())
apply!(global_bc, initial_state)

fg = NSFD.FG(initial_state)
fg = NSFD.compute_fg!(fg, initial_state)

rhs = NSFD.RHS(initial_state)
rhs = NSFD.compute_rhs!(rhs, fg, initial_state)

f = Figure(; size=(500, 500))
ax = Axis(f[1, 1])
h, a = plot!(ax, initial_state)

f
