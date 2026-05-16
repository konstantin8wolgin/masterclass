using JLD2, Plots, QuadGK, Statistics

# ── Load all data ─────────────────────────────────────────────────────────────
flow_data = load(joinpath(@__DIR__, "data", "flow.jld2"))
snapshots = flow_data["snapshots"]

# ── Theoretical distribution at t=1 ──────────────────────────────────────────
S1(φ) = -φ^2 + φ^4 / 8
Z1, _ = quadgk(φ -> exp(-S1(φ)), -10.0, 10.0; rtol=1e-10)
p1_theory(φ) = exp(-S1(φ)) / Z1

# ── Statistics table ──────────────────────────────────────────────────────────
snapshot_ts = [0.0, 0.25, 0.5, 0.75, 1.0]

function skewness(x)
    m = mean(x); s = std(x)
    mean((xi - m)^3 for xi in x) / s^3
end

println("┌────────┬────────────┬────────────┬────────────┐")
println("│   t    │    mean    │    std     │  skewness  │")
println("├────────┼────────────┼────────────┼────────────┤")
for ts in snapshot_ts
    s = snapshots[ts]
    m = round(mean(s),      digits=4)
    v = round(std(s),       digits=4)
    k = round(skewness(s),  digits=4)
    println("│ $(lpad(ts,6)) │ $(lpad(m,10)) │ $(lpad(v,10)) │ $(lpad(k,10)) │")
end
println("└────────┴────────────┴────────────┴────────────┘")

# ── KS statistic ──────────────────────────────────────────────────────────────
# Compare empirical CDF of final samples vs theoretical CDF
final = sort(snapshots[1.0])
N = length(final)

# Theoretical CDF at each sample point via quadrature
function theory_cdf(φ_val)
    val, _ = quadgk(p1_theory, -10.0, φ_val; rtol=1e-8)
    return val
end

println("\nComputing KS statistic (this may take a moment)...")
ks_vals = Float64[]
step_size = max(1, N ÷ 200)   # subsample for speed
for (i, φ_val) in enumerate(final[1:step_size:end])
    ecdf_val = (i * step_size) / N
    tcdf_val = theory_cdf(φ_val)
    push!(ks_vals, abs(ecdf_val - tcdf_val))
end
ks_stat = maximum(ks_vals)
println("KS statistic = $(round(ks_stat, digits=4))  (target < 0.05)")
println(ks_stat < 0.05 ? "PASS" : "FAIL — distribution mismatch")

# ── Combined figure ───────────────────────────────────────────────────────────
phi_range    = range(-5.5, 5.5, length=500)
theory_curve = p1_theory.(phi_range)
colors       = [:blue, :cyan, :green, :orange, :red]

p_snap = plot(title="Distribution evolution t=0 → 1",
              xlabel="φ", ylabel="density",
              legend=:topright)

for (ts, col) in zip(snapshot_ts, colors)
    histogram!(p_snap, snapshots[ts];
               bins=60, normalize=:pdf,
               label="t=$(ts)", fillalpha=0.35, color=col, linewidth=0)
end
plot!(p_snap, phi_range, theory_curve;
      label="Theory (t=1)", color=:black, linewidth=2.0, linestyle=:dash)

p_compare = plot(title="t=1: PIK vs theory",
                 xlabel="φ", ylabel="p_{t=1}(φ)",
                 legend=:topright)

histogram!(p_compare, snapshots[1.0];
           bins=60, normalize=:pdf,
           label="PIK samples", fillalpha=0.6, color=:steelblue, linewidth=0)
plot!(p_compare, phi_range, theory_curve;
      label="Theory", color=:orange, linewidth=2.5)

combined = plot(p_snap, p_compare, layout=(1,2), size=(1200, 450))
savefig(combined, joinpath(@__DIR__, "data", "05_summary.png"))
println("Saved data/05_summary.png")
