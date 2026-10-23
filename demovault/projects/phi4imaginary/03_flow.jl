using JLD2, Plots, Interpolations

# ── Load Wegner data ──────────────────────────────────────────────────────────
data = load("data/wegner.jld2")
phi_grid_real = data["phi_grid"]   # Vector{Float64}, length N_phi
t_grid        = data["t_grid"]     # Vector{Float64}, length N_t
phi_dot_grid  = data["phi_dot"]    # Matrix{ComplexF64}, size (N_t, N_phi)

# ── Build interpolant for phi_dot_t(phi) ─────────────────────────────────────
interp_re = interpolate((t_grid, phi_grid_real), real.(phi_dot_grid), Gridded(Linear()))
interp_im = interpolate((t_grid, phi_grid_real), imag.(phi_dot_grid), Gridded(Linear()))
extrap_re = extrapolate(interp_re, Flat())
extrap_im = extrapolate(interp_im, Flat())

function phi_dot_interp(phi_complex, t)
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
t_steps = collect(range(t_start, t_end, length=N_steps+1))
jldsave("data/flow.jld2";
    trajectories = trajectories,
    phi0_values  = collect(phi0_values),
    t_steps      = t_steps)
println("Saved to data/flow.jld2")

# ── Plot: trajectories in complex phi-plane ───────────────────────────────────
# Wegner flow in ℂ: real initial points (blue) flow into the Picard–Lefschetz
# thimble (PIKfold, red stars) — the steepest-descent contour of the Minkowski action.
p = plot(title="PIKfold formation: Wegner flow in ℂ-plane\n" *
               "real initial conditions → Picard–Lefschetz thimble",
         xlabel="Re[φ]", ylabel="Im[φ]",
         legend=:topright, aspect_ratio=:equal)

# Per-trajectory plot using line_z for proper colourbar in t
for n in 1:N_init
    traj = trajectories[n, :]
    plot!(p, real.(traj), imag.(traj),
          line_z=t_steps, c=:viridis, colorbar_title="t",
          linewidth=1.5, alpha=0.85, label="")
end

# Analytic free-theory thimble at t=t_end: φ → φ·exp(iπ t_end/2)
thimble_phi0 = range(-3.0, 3.0, length=100)
thimble_pts  = thimble_phi0 .* exp(im*π/2*t_end)
plot!(p, real.(thimble_pts), imag.(thimble_pts),
      color=:gray, linestyle=:dash, linewidth=2,
      label="free-theory thimble at t=$(round(t_end,digits=2))")

# Endpoints
scatter!(p, real.(trajectories[:, 1]), imag.(trajectories[:, 1]),
         marker=:circle, markersize=4, color=:blue,
         label="t=0 (real line)")
scatter!(p, real.(trajectories[:, end]), imag.(trajectories[:, end]),
         marker=:star5, markersize=6, color=:red,
         label="t=$(round(t_end,digits=2)) (PIKfold)")

savefig(p, "data/03_pikfold.png")
println("Plot saved to data/03_pikfold.png")
