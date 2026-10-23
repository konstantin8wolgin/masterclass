 1. §1 KDE headline gate — relative error of the complex-KDE
     estimator vs analytic on the 601-point [−6, 6] grid below the
     stated MC noise bound (~0.03 · max|e^{-S}|). Both Re and Im.
     2. GPU↔CPU velocity cross-check at the 7 anchor μ-values
     ({-0.675, 0.0, 0.1, 0.2, 1.5, 1.65, 1.75}) with (U=1, β=10):
     max|V_gpu − V_cpu_cached| < 1e-5 after the VEL_MAX_OK=50 clamp.
     Loads data/flow_path1_wegner.jld2.
     3. GL quadrature accuracy at the LHC corner (U=4, β=12, μ=2):
     compare GPU S̄ and V(φ=0, t=1) against an on-the-fly
     QuadGK(...; rtol=1e-8) CPU call. relerr < 1e-6.
     4. Rip-flag consistency at the 7 anchor μ-values: GPU-trajectory
     ripped::Bool matches the existing fermionic_one_site_flow.jl §3
     output exactly; max_im_excursion within 5% relative.
     5. LHC determinism: re-run with the fixed seed produces a
     bit-identical design matrix.
     6. JLD2 roundtrip: immediately after each jldsave in §3-§5, load
     and assert ≈ against the in-memory array (catches silent
     ComplexF64 → ComplexF32 drift).
     7. Visual sanity: open data/10_rip_scatter_3d.png and confirm
     both green and red points are present and not collapsed into a
     single corner. Open data/09_detail_case_1.png (best-success case)
     and confirm Re/Im histogram bars hug the smooth theory curve;
     open data/09_detail_case_4.png (worst rip) and confirm visible
     deviation.