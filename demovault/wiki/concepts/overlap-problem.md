---
title: "Overlap problem"
tags: [overlap-problem, sign-problem, importance-sampling]
sources: ["[[aarts-2026-lattice-sign-problem-review]]"]
---

## Definition

The overlap problem is the inefficiency of Monte Carlo sampling that arises when the simulated distribution (e.g. phase-quenched) has poor overlap with the target distribution (e.g. full QCD at non-zero $\mu$). The dominant configurations of the two theories occupy different regions of configuration space, so reweighting from one to the other requires exponentially many samples.

## Intuition

Importance sampling works by drawing configurations from a tractable distribution $p_0$ and reweighting by $p/p_0$. If $p$ and $p_0$ have negligible overlap — rare regions under $p_0$ dominate $p$ — the variance of the reweighting factor is exponentially large, and the effective sample size is exponentially small. This is a distinct problem from the sign problem (which arises from phase cancellations) but both lead to exponential scaling and are often related.

## Key results

- In QCD: phase-quenched theory at $\mu_B/T \gtrsim 3$ has a Bose-Einstein condensate of pions (due to wrong isospin chemical potential), while full QCD does not. The configurations that dominate $Z_{pq}$ are not relevant to $Z$: exponential overlap problem.
- The SNR for the reweighting factor $e^{i\varphi}$ is $\text{SNR} \sim e^{-(F-F_{pq})/T}$, same exponential decay as the sign problem — both have the same root cause.
- The silver blaze problem ([[silver-blaze-problem]]) is a specific manifestation of the overlap problem at $\mu < \mu_c$.
- Methods that deform the integration manifold (thimbles, PIKfolds, contour deformations) address both the sign and overlap problems simultaneously by moving to a manifold where both $p$ and $p_0$ have the same dominant configurations.

## Connections

- [[sign-problem]] — sign and overlap problems are related and often occur together.
- [[reweighting]] — fails precisely because of the overlap problem.
- [[silver-blaze-problem]] — specific instance of overlap problem at low $\mu$.
- [[lefschetz-thimbles]] — manifold deformation addresses both.
- [[physics-informed-kernels]] — weight-preserving property explicitly avoids the overlap problem.
