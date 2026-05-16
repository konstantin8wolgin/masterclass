using JLD2, Plots, Printf

# ── Load all checkpoint data ──────────────────────────────────────────────────
obs  = load("data/observables.jld2")
flow = load("data/flow.jld2")
weg  = load("data/wegner.jld2")

t_obs      = obs["t_grid"]
phi2       = obs["phi2"]
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

# ── 2×2 combined figure ───────────────────────────────────────────────────────

# Panel 1: velocity field Re[φ̇_t] vs φ
p_vel = plot(title="Velocity field Re[φ̇_t(φ)]",
             xlabel="φ", ylabel="Re[φ̇_t]", legend=:topright)
for (tv, col) in zip([0.0, 0.25, 0.5, 0.75, 0.99], [:blue,:green,:orange,:red,:purple])
    idx = argmin(abs.(t_grid_w .- tv))
    plot!(p_vel, phi_grid_r, real.(phi_dot_grid[idx, :]),
          label="t=$(round(t_grid_w[idx],digits=2))", color=col)
end

# Panel 2: PIKfold trajectories
p_fold = plot(title="PIKfold formation",
              xlabel="Re[φ]", ylabel="Im[φ]",
              legend=false, aspect_ratio=:equal)
cmap = cgrad(:viridis)
for n in axes(trajectories, 1)
    traj = trajectories[n, :]
    for step in 1:length(traj)-1
        c = (step-1) / (length(traj)-1)
        plot!(p_fold,
            [real(traj[step]), real(traj[step+1])],
            [imag(traj[step]), imag(traj[step+1])],
            color=cmap[c], linewidth=1.2, alpha=0.7)
    end
end
scatter!(p_fold, real.(trajectories[:, 1]),   zeros(size(trajectories,1)),
         color=:blue, markersize=3, label="t=0")
scatter!(p_fold, real.(trajectories[:, end]), imag.(trajectories[:, end]),
         color=:red, markersize=4, marker=:star5, label="t=t_end")

# Panel 3: observables vs t
p_obs = plot(t_obs, real.(phi2),   label="Re⟨φ²⟩",   lw=2, color=:blue,
             xlabel="t", ylabel="Observable", title="Observables vs t")
plot!(p_obs, t_obs, real.(absphi), label="Re⟨|φ|⟩",  lw=2, color=:red)

# Panel 4: sign severity
p_sign = plot(t_obs, abs.(avg_sign),
              xlabel="t", ylabel="|avg sign|",
              title="Sign problem severity",
              lw=2, color=:black, ylims=(0, 1.05), legend=false)
hline!(p_sign, [1.0], linestyle=:dash, color=:gray)

fig = plot(p_vel, p_fold, p_obs, p_sign,
           layout=(2,2), size=(1100, 800))
savefig(fig, "data/05_summary.png")
println("\nCombined figure saved to data/05_summary.png")
