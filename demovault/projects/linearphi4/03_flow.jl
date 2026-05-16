using JLD2, Plots, Interpolations, Random, Statistics

# ── Load Wegner data ──────────────────────────────────────────────────────────
data = load(joinpath(@__DIR__, "data", "wegner.jld2"))
phi_grid_real = data["phi_grid"]   # Vector{Float64}, length N_phi
t_grid        = data["t_grid"]     # Vector{Float64}, length N_t
phi_dot_grid  = data["phi_dot"]    # Matrix{Float64}, size (N_t, N_phi)

# ── Build bilinear interpolant ─────────────────────────────────────────────────
# phi_dot_grid is indexed as [t_index, phi_index]
interp = interpolate((t_grid, phi_grid_real), phi_dot_grid, Gridded(Linear()))
extrap = extrapolate(interp, Flat())

function phi_dot_interp(φ::Real, t::Real)
    t_c = clamp(t, t_grid[1],        t_grid[end])
    φ_c = clamp(φ, phi_grid_real[1], phi_grid_real[end])
    return extrap(t_c, φ_c)
end

# ── RK4 integrator ────────────────────────────────────────────────────────────
function rk4_step(φ::Real, t::Real, dt::Real)
    k1 = phi_dot_interp(φ,             t)
    k2 = phi_dot_interp(φ + dt/2 * k1, t + dt/2)
    k3 = phi_dot_interp(φ + dt/2 * k2, t + dt/2)
    k4 = phi_dot_interp(φ + dt   * k3, t + dt)
    return φ + dt/6 * (k1 + 2k2 + 2k3 + k4)
end

# ── Sample from Gaussian at t=0 ───────────────────────────────────────────────
const N_s     = 10_000
const N_steps = 50
const dt      = 1.0 / N_steps

Random.seed!(42)
phi0 = randn(N_s)

# ── Flow samples t=0 → t=1, recording snapshots ───────────────────────────────
snapshot_ts = [0.0, 0.25, 0.5, 0.75, 1.0]
snapshots   = Dict{Float64, Vector{Float64}}()
snapshots[0.0] = copy(phi0)

phi_now = copy(phi0)
println("Flowing $N_s samples with $N_steps RK4 steps...")
for step in 1:N_steps
    t_now  = (step - 1) * dt
    phi_now .= rk4_step.(phi_now, t_now, dt)
    t_next = step * dt
    for ts in snapshot_ts
        ts == 0.0 && continue
        if abs(t_next - ts) <= dt / 2 + 1e-12
            snapshots[ts] = copy(phi_now)
        end
    end
    step % 10 == 0 && println("  step $step / $N_steps  (t=$(round(t_next, digits=2))) done")
end
println("Done.")

# ── Quick sanity check ────────────────────────────────────────────────────────
final = snapshots[1.0]
println("\nFinal distribution (t=1):")
println("  mean = $(round(mean(final), digits=4))  (expect ≈ 0 by Z₂ symmetry)")
println("  std  = $(round(std(final),  digits=4))  (expect ≈ 2)")

# ── Save ──────────────────────────────────────────────────────────────────────
jldsave(joinpath(@__DIR__, "data", "flow.jld2"); snapshots, phi0)
println("Saved to data/flow.jld2")

# ── Plot: representative trajectories ────────────────────────────────────────
# Track 20 evenly-spaced initial conditions through the flow for visualisation
const N_traj  = 20
const N_track = N_steps + 1
phi_track_0 = range(-3.0, 3.0, length=N_traj) |> collect

trajectories = zeros(Float64, N_traj, N_track)
trajectories[:, 1] = phi_track_0

phi_traj = copy(phi_track_0)
for step in 1:N_steps
    t_now = (step - 1) * dt
    phi_traj .= rk4_step.(phi_traj, t_now, dt)
    trajectories[:, step + 1] = phi_traj
end

t_axis = range(0.0, 1.0, length=N_track)
cmap   = cgrad(:viridis)

p = plot(title="PIK sample trajectories", xlabel="t", ylabel="φ_t",
         legend=false, size=(700, 450))

for n in 1:N_traj
    color_val = (phi_track_0[n] + 3.0) / 6.0
    plot!(p, t_axis, trajectories[n, :], color=cmap[color_val], linewidth=1.2, alpha=0.8)
end

hline!(p, [2.0, -2.0], color=:black, linestyle=:dash, linewidth=0.8, alpha=0.5)
savefig(p, joinpath(@__DIR__, "data", "03_trajectories.png"))
println("Plot saved to data/03_trajectories.png")
