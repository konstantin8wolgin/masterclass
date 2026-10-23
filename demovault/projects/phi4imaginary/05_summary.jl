using JLD2, Plots, Printf

# ── Load all checkpoint data ──────────────────────────────────────────────────
obs  = load("data/observables.jld2")
flow = load("data/flow.jld2")
weg  = load("data/wegner.jld2")

t_obs      = obs["t_grid"]
phi2       = obs["phi2"]
phi4       = obs["phi4"]
absphi     = obs["absphi"]
avg_sign   = obs["avg_sign"]

trajectories = flow["trajectories"]
phi0_values  = flow["phi0_values"]
t_flow       = flow["t_steps"]

phi_grid_r   = weg["phi_grid"]
t_grid_w     = weg["t_grid"]
phi_dot_grid = weg["phi_dot"]

# ── Printed table ─────────────────────────────────────────────────────────────
println("\n╔══════╦══════════════╦══════════════╦══════════════╦════════════╗")
println("║  t   ║  Re⟨φ²⟩     ║  Im⟨φ²⟩     ║  Re⟨|φ|⟩   ║ |avg sign| ║")
println("╠══════╬══════════════╬══════════════╬══════════════╬════════════╣")
for t_val in [0.0, 0.25, 0.5, 0.75, 1.0]
    idx = argmin(abs.(t_obs .- t_val))
    @printf("║ %.2f ║ %12.6f ║ %12.6f ║ %12.6f ║ %10.6f ║\n",
        t_obs[idx],
        real(phi2[idx]), imag(phi2[idx]),
        real(absphi[idx]),
        abs(avg_sign[idx]))
end
println("╚══════╩══════════════╩══════════════╩══════════════╩════════════╝")

# ── Combined figure (2×3) ─────────────────────────────────────────────────────
M2 = 1.0
plot_t_slices = [0.0, 0.25, 0.5, 0.75, 1.0]

# Panel 1: velocity field Re[φ̇_t] vs φ
p_vel = plot(title="Velocity field Re[φ̇_t(φ)]",
             xlabel="φ", ylabel="Re[φ̇_t]", legend=:topright)
for (tv, col) in zip(plot_t_slices, [:blue,:green,:orange,:red,:purple])
    idx = argmin(abs.(t_grid_w .- tv))
    plot!(p_vel, phi_grid_r, real.(phi_dot_grid[idx, :]),
          label="t=$(round(t_grid_w[idx],digits=2))", color=col)
end

# Panel 2: PIKfold trajectories
p_fold = plot(title="PIKfold formation",
              xlabel="Re[φ]", ylabel="Im[φ]",
              legend=false, aspect_ratio=:equal)
for n in axes(trajectories, 1)
    traj = trajectories[n, :]
    plot!(p_fold, real.(traj), imag.(traj),
          line_z=t_flow, c=:viridis, colorbar=false,
          linewidth=1.2, alpha=0.7, label="")
end
scatter!(p_fold, real.(trajectories[:, 1]), zeros(size(trajectories,1)),
         color=:blue, markersize=3, label="t=0")
scatter!(p_fold, real.(trajectories[:, end]), imag.(trajectories[:, end]),
         color=:red, markersize=4, marker=:star5, label="t=t_end")

# Panel 3: observables ⟨φ²⟩ vs t with free-theory overlay
t_vec = collect(t_obs)
free  = @. exp(im*π/2*t_vec) / M2
p_obs = plot(t_obs, real.(phi2), label="Re⟨φ²⟩", lw=2, color=:blue,
             xlabel="t", ylabel="Observable", title="⟨φ²⟩ vs t (free dotted)")
plot!(p_obs, t_obs, imag.(phi2), label="Im⟨φ²⟩", lw=2, color=:blue, linestyle=:dash)
plot!(p_obs, t_obs, real.(absphi), label="Re⟨|φ|⟩", lw=2, color=:red)
plot!(p_obs, t_vec, real.(free), label="Re free",  lw=1.5, color=:black, linestyle=:dot)
plot!(p_obs, t_vec, imag.(free), label="Im free",  lw=1.5, color=:gray,  linestyle=:dot)

# Panel 4: sign severity
p_sign = plot(t_obs, abs.(avg_sign),
              xlabel="t", ylabel="|avg sign|",
              title="Sign problem severity",
              lw=2, color=:black, ylims=(0, 1.05), legend=false)
hline!(p_sign, [1.0], linestyle=:dash, color=:gray)

# Panel 5: ⟨φ⁴⟩ vs t
p_phi4 = plot(t_obs, real.(phi4), label="Re⟨φ⁴⟩", lw=2, color=:darkorange,
              xlabel="t", ylabel="⟨φ⁴⟩", title="Quartic observable")
plot!(p_phi4, t_obs, imag.(phi4), label="Im⟨φ⁴⟩", lw=2, color=:darkorange, linestyle=:dash)

# Panel 6: |⟨φ²⟩| vs t (rotation magnitude) — should ≈ 1/m² for free theory
p_mag = plot(t_obs, abs.(phi2), label="|⟨φ²⟩|", lw=2, color=:purple,
             xlabel="t", ylabel="|⟨φ²⟩|", title="Rotation magnitude")
hline!(p_mag, [1.0/M2], linestyle=:dash, color=:gray, label="1/m²")

fig = plot(p_vel, p_fold, p_obs, p_sign, p_phi4, p_mag,
           layout=(2,3), size=(1400, 900))
savefig(fig, "data/05_summary.png")
println("\nCombined figure saved to data/05_summary.png")
