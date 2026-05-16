using JLD2, Plots, QuadGK

# ── Load flow data ────────────────────────────────────────────────────────────
data      = load(joinpath(@__DIR__, "data", "flow.jld2"))
snapshots = data["snapshots"]   # Dict{Float64, Vector{Float64}}

# ── Theoretical distribution at t=1 ──────────────────────────────────────────
# S_1(φ) = -φ² + φ⁴/8   (m²(1)=-2, λ(1)=0.5)
S1(φ) = -φ^2 + φ^4 / 8
Z1, _ = quadgk(φ -> exp(-S1(φ)), -10.0, 10.0; rtol=1e-10)
p1_theory(φ) = exp(-S1(φ)) / Z1

# ── Plot 1: t=1 histogram vs theoretical curve ────────────────────────────────
final_samples = snapshots[1.0]

phi_range = range(-5.5, 5.5, length=500)
theory_curve = p1_theory.(phi_range)

p_compare = plot(title="Distribution at t=1: PIK flow vs theory",
                 xlabel="φ", ylabel="p_{t=1}(φ)",
                 legend=:topright, size=(700, 450))

histogram!(p_compare, final_samples;
           bins=60, normalize=:pdf,
           label="PIK samples", fillalpha=0.6, color=:steelblue, linewidth=0)

plot!(p_compare, phi_range, theory_curve;
      label="Theory", color=:orange, linewidth=2.5)

savefig(p_compare, joinpath(@__DIR__, "data", "04_t1_comparison.png"))
println("Saved data/04_t1_comparison.png")

# ── Plot 2: snapshot evolution ────────────────────────────────────────────────
snapshot_ts = [0.0, 0.25, 0.5, 0.75, 1.0]
colors       = [:blue, :cyan, :green, :orange, :red]

p_snap = plot(title="Distribution evolution t=0 → 1",
              xlabel="φ", ylabel="density",
              legend=:topright, size=(700, 450))

for (ts, col) in zip(snapshot_ts, colors)
    histogram!(p_snap, snapshots[ts];
               bins=60, normalize=:pdf,
               label="t=$(ts)", fillalpha=0.35, color=col, linewidth=0)
end

# Add theory at t=1 for reference
plot!(p_snap, phi_range, theory_curve;
      label="Theory (t=1)", color=:black, linewidth=2.0, linestyle=:dash)

savefig(p_snap, joinpath(@__DIR__, "data", "04_snapshots.png"))
println("Saved data/04_snapshots.png")
