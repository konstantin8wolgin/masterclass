using JLD2, Plots, Interpolations

# ── Load Wegner data ──────────────────────────────────────────────────────────
data = load("data/wegner.jld2")
phi_grid_real = data["phi_grid"]   # Vector{Float64}, length N_phi
t_grid        = data["t_grid"]     # Vector{Float64}, length N_t
phi_dot_grid  = data["phi_dot"]    # Matrix{ComplexF64}, size (N_t, N_phi)

# ── Build interpolant for phi_dot_t(phi) ─────────────────────────────────────
# Separate real and imag parts (Interpolations.jl requires real-valued arrays)
interp_re = interpolate((t_grid, phi_grid_real), real.(phi_dot_grid), Gridded(Linear()))
interp_im = interpolate((t_grid, phi_grid_real), imag.(phi_dot_grid), Gridded(Linear()))

# Extrapolation: clamp to grid boundary
extrap_re = extrapolate(interp_re, Flat())
extrap_im = extrapolate(interp_im, Flat())

function phi_dot_interp(phi_complex, t)
    # Only Re[φ] is used: the velocity field is defined on the real-φ axis, and using Re[φ] is the standard approach for PIKfold visualization.
    φ_r = clamp(real(phi_complex), phi_grid_real[1], phi_grid_real[end])
    t_c = clamp(t, t_grid[1], t_grid[end])
    return extrap_re(t_c, φ_r) + im * extrap_im(t_c, φ_r)
end

# ── RK4 integrator ────────────────────────────────────────────────────────────
function rk4_step(φ, t, dt)
    k1 = phi_dot_interp(φ,        t)
    k2 = phi_dot_interp(φ + dt/2*k1, t + dt/2)
    k3 = phi_dot_interp(φ + dt/2*k2, t + dt/2)
    k4 = phi_dot_interp(φ + dt*k3,   t + dt)
    return φ + dt/6 * (k1 + 2k2 + 2k3 + k4)
end

# ── Flow integration ──────────────────────────────────────────────────────────
const N_init   = 20
const N_steps  = 512
const t_start  = t_grid[1]
const t_end    = t_grid[end]
const dt       = (t_end - t_start) / N_steps

phi0_values = range(-3.0, 3.0, length=N_init)
# trajectories[n, step] = complex phi_t for initial condition n at step
trajectories = zeros(ComplexF64, N_init, N_steps + 1)

println("Integrating flow for $N_init initial conditions, $N_steps RK4 steps...")
for (n, phi0) in enumerate(phi0_values)
    trajectories[n, 1] = phi0
    for step in 1:N_steps
        t_now = t_start + (step-1)*dt
        trajectories[n, step+1] = rk4_step(trajectories[n, step], t_now, dt)
    end
end
println("Done.")

# ── Save ──────────────────────────────────────────────────────────────────────
jldsave("data/flow.jld2";
    trajectories = trajectories,
    phi0_values  = collect(phi0_values),
    t_steps      = collect(range(t_start, t_end, length=N_steps+1)))
println("Saved to data/flow.jld2")

# ── Plot: trajectories in complex phi-plane ───────────────────────────────────
p = plot(title="PIKfold formation: field trajectories in ℂ",
         xlabel="Re[φ]", ylabel="Im[φ]",
         legend=:topright, aspect_ratio=:equal)

cmap = cgrad(:viridis)

for n in 1:N_init
    traj = trajectories[n, :]
    # color segments by t
    for step in 1:N_steps
        color_val = (step-1) / N_steps
        plot!(p,
            [real(traj[step]), real(traj[step+1])],
            [imag(traj[step]), imag(traj[step+1])],
            color=cmap[color_val], linewidth=1.5, alpha=0.8)
    end
end

# Mark start (t=0, real line) and end (t=t_end, PIKfold)
scatter!(p, real.(trajectories[:, 1]),   imag.(trajectories[:, 1]),
         marker=:circle, markersize=4, color=:blue,   label="t=0 (real line)")
scatter!(p, real.(trajectories[:, end]), imag.(trajectories[:, end]),
         marker=:star5,  markersize=6, color=:red,    label="t=$(round(t_end,digits=2)) (PIKfold)")

savefig(p, "data/03_pikfold.png")
println("Plot saved to data/03_pikfold.png")
