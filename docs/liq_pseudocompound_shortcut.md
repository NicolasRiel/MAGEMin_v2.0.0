# Synthesizing redundant solution-phase pseudocompound occurrences

## 1. Problem statement

Today, in `ss_min_LP` (`src/ss_min_function.c`), **every** active pseudocompound
`cp[i]` of a solution phase gets its own full NLopt local minimization each PGE
cycle, even when the same phase (`ph_id`) occurs `nOcc = gv.n_ss_ph[ph_id]`
times simultaneously. At high T, when a lot of melt is stable, `liq` commonly
has `nOcc` in the dozens, and — because all `nOcc` instances share the same
rotated hyperplane (`Gamma`, fixed for the whole PGE cycle) and are refined
from nearby starting points inside `restrict_SS_HyperVolume` boxes — they
routinely converge to the same (or a near-identical) point on the Gibbs
surface. Every one of those solves after the first is redundant work.

This document describes a **composite method** that restructures one PGE
cycle's `ss_min_LP` → `run_LP` sequence into three stages:

**(A) Minimize once per unique phase, not once per occurrence.** For each
solution phase present this cycle, run exactly **one** real NLopt solve
(`rotate_hyperplane` + `restrict_SS_HyperVolume` + `NLopt_opt[ph_id]`),
regardless of how many `cp[i]` occurrences it currently has.

**(B) Re-estimate Gamma by least squares across all the minimized phases.**
Using the one true minimum found per phase in (A) together with each
endmember's chemical potential there, fit a refined `Gamma` by least squares
— sharper than whatever `Gamma` `ss_min_LP` started the cycle with, since it
is informed by every phase's minimum simultaneously rather than derived from
the LP's own dual alone.

**(C) Generate `nOcc` synthetic pseudocompounds per phase, analytically, with
no second minimization.** For each phase with `nOcc > threshold`, construct a
set of points that (i) lie on (a linearization of) the newly-estimated
`Gamma` hyperplane, and (ii) are built by perturbing the *already-minimized*
point from (A) — never by re-running NLopt with the new `Gamma`.

Only then does the cycle proceed into `run_LP` as today (swap pseudocompounds
in/out based on the now-current `Gamma`, update phase fractions, etc.).

## 2. Notation

- `ph_id` — a solution phase (e.g. `liq`).
- `N = nOcc` — number of simultaneous occurrences of `ph_id` this cycle
  (`gv.n_ss_ph[ph_id]`); the shortcut only fires for `N > threshold` (e.g.
  `threshold = 2`, exposed as an option).
- `n_ox` — number of oxide components (dimension of `Gamma` and of each
  endmember's composition vector).
- `n_em(ph_id)` — number of endmembers of phase `ph_id`; `C^{(ph_id)}` is its
  fixed `n_ox × n_em` composition matrix (`Comp[i][j]`, a property of the
  database, not of any particular point), so composition is exactly linear
  in `xeos`: `Comp(x) = C\,x`.
- `x* = x*_{ph_id} ∈ R^{n_xeos}` — the single true local minimum found in
  stage (A) for phase `ph_id` this cycle (one real
  `rotate_hyperplane` + `restrict_SS_HyperVolume` + `NLopt_opt[ph_id]` call).
- `G(x*)`, `∇G(x*)`, `μ(x*)` — the Gibbs energy, its `xeos`-gradient, and the
  vector of endmember chemical potentials at `x*`, all already computed as a
  byproduct of the stage-(A) NLopt solve (no extra evaluation needed).
- `Γ_old` — the `Gamma` `ss_min_LP` started this cycle with (from the
  previous cycle's `run_LP`).
- `Γ_new` — the refined estimate produced by stage (B).
- `x_1, ..., x_N` — the **previous** cycle's converged `xeos` for the `N`
  occurrences of `ph_id` (available as `cp[i].xeos` before this cycle
  overwrites them), `x̄ = (1/N)Σ_k x_k`, `e_k = x_k - x̄` (so `Σ_k e_k = 0`
  exactly — see §5).
- `h` — target step size for the synthetic spread, in `xeos` units
  (proposed ≈ `1e-4`).

## 3. Stage A: minimize each unique phase once

Instead of looping over every `cp[i]` and minimizing each, loop over the
**distinct** `ph_id` values present in `cp[]` this cycle. For each, run the
existing per-occurrence block exactly once, using (for example) the
occurrence with the best current `xeos` guess as the starting point:

```
for k in n_xeos: SS_ref_db[ph_id].iguess[k] = <chosen occurrence>.xeos[k]
SS_ref_db[ph_id] = rotate_hyperplane(gv, SS_ref_db[ph_id])          // uses Γ_old
SS_ref_db[ph_id] = restrict_SS_HyperVolume(gv, SS_ref_db[ph_id], gv.box_size_mode_PGE)
SS_ref_db[ph_id] = (*NLopt_opt[ph_id])(gv, SS_ref_db[ph_id])         // → x*, G(x*), ∇G(x*), μ(x*)
```

This is unchanged machinery — the only change is *how many times* it runs:
once per unique phase (typically single digits), not once per `cp[i]`
occurrence (which for `liq` at high T can be dozens).

## 4. Stage B: least-squares refinement of Gamma

Each minimized phase gives `n_em(ph_id)` endmember chemical potentials
`μ_j(x*)` and their known composition vectors `Comp_j = C^{(ph_id)}_{:,j}
∈ R^{n_ox}`. Stacking every endmember from every minimized phase this cycle
gives an overdetermined linear system in the unknown `Γ ∈ R^{n_ox}`:

$$ \text{Comp}_j \cdot \Gamma \;\approx\; \mu_j, \qquad j = 1, \dots, \sum_{ph} n_{em}(ph) $$

or in matrix form, `M\,\Gamma \approx \mu` where `M` stacks all the `Comp_j^\top`
rows and `μ` stacks all the `μ_j`. This is solved by ordinary least squares:

$$ \Gamma_{\text{new}} \;=\; (M^\top M)^{-1} M^\top \mu. $$

This is well-determined as long as the number of endmembers collected across
all minimized phases is at least `n_ox` and `M` has full column rank (true in
practice once more than a couple of phases are minimized — pure phases
(`pp_min_function`, already minimized earlier in the cycle) can also
contribute their single `Comp/μ` row here for extra conditioning, since they
are exact points on the true hyperplane by construction).

`Γ_new` is therefore informed by *every* phase's minimum jointly, rather than
being the single-phase-at-a-time dual that `rotate_hyperplane` used going
into stage A.

## 5. Stage C: generating `N` points on the `Γ_new` hyperplane around `x*`

This is the key correction from the earlier draft of this document: given
`Γ_new`, we do **not** re-run `rotate_hyperplane` + NLopt per phase. Instead
we build the `N` synthetic points analytically from quantities already in
hand from stage A (`x*`, `G(x*)`, `∇G(x*)`) and from stage B (`Γ_new`).

### 5.1 The rotated (driving-force) energy relative to `Γ_new`

$$ \Delta G(x) \;=\; G(x) \;-\; \Gamma_{\text{new}} \cdot \text{Comp}(x) \;=\; G(x) - \Gamma_{\text{new}}\cdot C\,x. $$

Since `Comp` is exactly linear in `x`, linearizing only `G` around `x*` gives
an expression that is exact in the composition term and first-order accurate
in `G`:

$$ \Delta G(x^{*}+\delta) \;\approx\; \underbrace{\big[G(x^{*}) - \Gamma_{\text{new}}\cdot C\,x^{*}\big]}_{\Delta G_0} \;+\; \underbrace{\big[\nabla G(x^{*}) - C^\top \Gamma_{\text{new}}\big]}_{g_{\text{eff}}} \cdot\ \delta. $$

Both `ΔG_0` (a scalar) and `g_eff` (a vector in `R^{n_xeos}`) are **closed
form** from quantities already computed in stages A and B — no new call to
the objective function is needed to get them.

### 5.2 Landing on the hyperplane: the base point `x**`

Setting `ΔG(x*+δ) = 0` to first order is one linear equation in `δ`. Its
minimum-norm solution is

$$ \delta_0 \;=\; -\frac{\Delta G_0}{\lVert g_{\text{eff}}\rVert^2}\, g_{\text{eff}}, \qquad x^{**} \;=\; x^{*} + \delta_0. $$

`x**` is the new base point: to first order, it sits exactly on the
linearized `Γ_new` hyperplane. `δ_0` is small whenever `Γ_new ≈ Γ_old`
(i.e. Gamma hasn't moved much this cycle, which is the common/converging
case) and is the honest, closed-form correction otherwise.

### 5.3 Spreading `N` points while staying on the hyperplane

Any direction `v` with `g_eff · v = 0` keeps `ΔG ≈ 0` to first order when
stepping away from `x**` — i.e. the hyperplane's own tangent space at `x**`
is exactly the null space of `g_eff`. Reuse the previous occurrence set's
shape for that spread (as in the original construction, §5 of the prior
draft) by projecting its centered deviations `e_k = x_k - x̄` into that null
space:

$$ P \;=\; I \;-\; \frac{g_{\text{eff}}\,g_{\text{eff}}^\top}{\lVert g_{\text{eff}}\rVert^2}, \qquad v_k \;=\; P\,e_k. $$

`P` is a linear projector, so it preserves the zero-sum property exactly:

$$ \sum_{k=1}^N v_k \;=\; P \sum_{k=1}^N e_k \;=\; P\cdot 0 \;=\; 0, $$

even though each individual `v_k` has been rotated to satisfy
`g_eff · v_k = 0` exactly. The final construction:

$$ \boxed{\; x_k \;=\; x^{**} \;+\; s\, P\, e_k \;=\; x^{*} + \delta_0 + s\,P(x_k^{\text{old}} - \bar{x}), \qquad k = 1, \dots, N. \;} $$

with `s = h / max_k‖P e_k‖` (same normalization idea as before, §6 below).
By linearity this guarantees, **exactly**:

$$ \frac{1}{N}\sum_{k=1}^N x_k \;=\; x^{**} \;=\; x^{*} + \delta_0, $$

i.e. every synthesized point satisfies the hyperplane condition to first
order, *and* their average is the known, closed-form point `x**` — both
properties fall out of the same linear construction simultaneously.

**Honest tension to note.** If `δ_0 ≠ 0` (i.e. `Γ_new` actually moved
relative to what `x*` was optimal for), the resulting average `x**` is *not*
exactly the old average `x̄` — it is `x̄`'s natural successor once the phase
picture is updated to the new hyperplane. This is a real, small effect, not
a flaw: it is precisely how much the phase's aggregate composition *should*
move given the newly-refined Gamma. It vanishes whenever Gamma is converged
(`δ_0 → 0`).

## 6. Choosing `s` from the target step size `h`

Exactly as in the original construction, normalize by the largest projected
deviation so the spread magnitude is meaningful regardless of the previous
occurrence set's raw scale:

$$ s \;=\; \frac{h}{\max_k \lVert P e_k \rVert}, \qquad h \approx 10^{-4}. $$

This guarantees `max_k‖s\,P e_k‖ = h` exactly: no synthetic point moves more
than `h` away from `x**` along the hyperplane's tangent directions.

## 7. Filling in G, μ, and composition: evaluate, don't (only) trust the linear model

For each synthetic `x_k`, clamp into `bounds_ref`, then call the real
objective function **once**, directly (no NLopt search), exactly as in the
original draft:

$$ G_k = \texttt{obj\_gh\_}\{ph\_id\}(x_k), \qquad \nabla G_k,\ \mu_k,\ p_k\ \text{all filled in by the same call.} $$

This is important for two reasons, one reused from the original draft and
one new:

- (Reused) It is exact, not an approximation: since `x_k` is near but not at
  `x*`, and `x*` is (to NLopt's tolerance) stationary, `G(x_k)` is
  automatically slightly above `G(x*)` via real local curvature — no
  explicit hyperplane-fitting of `G` itself is needed.
- (New) §5's `ΔG(x_k) ≈ 0` guarantee is only first-order in `G`'s own
  nonlinearity (§5.1 dropped the quadratic term). Evaluating the real
  objective at each `x_k` gives the *true* `ΔG(x_k) = G(x_k) -
  Γ_new·Comp(x_k)`, which will differ from exactly `0` by `O(‖δ‖²)`
  curvature — cheap to check directly (compare against `ΔG_0`) and, if ever
  found to be non-negligible, an early signal that `h` should shrink or that
  `Γ_new` moved too far this cycle for the linearization to be trusted.

Since evaluating the real objective directly costs on the order of a
microsecond (measured for `liq`), the total cost of the composite method for
one phase this cycle is:

$$ \text{1 full NLopt solve (stage A)} \;+\; (N-1)\times(\text{~1 μs direct evaluation, stage C}), $$

versus today's `N` full NLopt solves — plus, across *all* phases, one shared
least-squares solve (stage B) sized `n_ox × n_ox`, negligible next to even a
single NLopt call.

## 8. Edge cases and guards

- **`‖g_eff‖ ≈ 0`.** This means `x*` is already (numerically) stationary for
  `Γ_new` as well as `Γ_old` — the projector `P` is undefined. Fall back to
  `δ_0 = 0`, `P = I` (no projection): the phase's `Γ` didn't really move, so
  the original (pre-hyperplane-correction) construction is already exactly
  what's wanted.
- **Degenerate previous spread.** If `max_k‖P e_k‖` is at or below machine
  epsilon (previous occurrences had already collapsed to numerically
  identical points, or all `e_k` happened to lie along `g_eff` and got
  projected to ~0), fall back to a fixed zero-sum direction set spanning the
  null space of `g_eff` (e.g. an orthonormal basis of `null(g_eff)`, or the
  vertices of a regular `(N-1)`-simplex projected through `P`), scaled to
  `h`.
- **`N = 1`.** No shortcut applies — run the real solve as today.
- **Bounds.** Clamp every `x_k` into `SS_ref_db[ph_id].bounds_ref` after
  construction.
- **Trigger condition.** Gate stage C on `N > threshold` (proposed default
  `2`, exposed as a `gv` option) **and** on this not being the very first
  PGE cycle for this bulk/step (no previous occurrence set `x_1,...,x_N` to
  reuse on iteration 0 — that cycle runs the full `N` solves once, same as
  today, to establish it). Stages A and B apply regardless of `N` (they are
  strictly a win — fewer NLopt calls, better-conditioned `Gamma` — even for
  phases with `N = 1`).

## 9. Where this hooks into the existing code

- **`ss_min_LP`** (`src/ss_min_function.c`): restructure the per-`cp[i]`
  loop into a per-unique-`ph_id` loop for stage A (§3); collect `(Comp_j,
  μ_j)` pairs from each minimized phase (and from already-minimized pure
  phases) for stage B; after solving for `Γ_new` (§4), loop again over
  phases with `N > threshold` to synthesize the `N-1` remaining occurrences
  per §5–§7 (mirroring how `copy_to_Ppc` already calls
  `(*SS_objective[ph_id])(...)` once, unconditionally, to fill in an
  "unrotated" `G` for a given `xeos` — the same direct-call pattern, just
  without a preceding NLopt search).
- **`run_LP`** (`src/PGE_function.c`): unchanged — it continues to consume
  whatever `cp[]` set it's handed (now filled in by a mix of real solves and
  stage-C synthesis) and swap/update as today. The refined `Γ_new` from
  stage B can also seed `run_LP`'s own Gamma update, rather than being
  discarded after stage C.

## 10. Summary

$$ \Delta G_0 = G(x^{*}) - \Gamma_{\text{new}}\cdot C x^{*}, \qquad g_{\text{eff}} = \nabla G(x^{*}) - C^\top \Gamma_{\text{new}} $$

$$ \delta_0 = -\frac{\Delta G_0}{\lVert g_{\text{eff}}\rVert^2} g_{\text{eff}}, \qquad P = I - \frac{g_{\text{eff}} g_{\text{eff}}^\top}{\lVert g_{\text{eff}}\rVert^2} $$

$$ x_k = x^{*} + \delta_0 + \frac{h}{\max_j\lVert P e_j\rVert}\, P\,(x_k^{\text{old}} - \bar x), \qquad k = 1,\dots,N $$

with the exact guarantee (to the linearization in §5.1)

$$ \frac{1}{N}\sum_{k=1}^N x_k = x^{*} + \delta_0 \quad\text{and}\quad \Delta G(x_k) \approx 0 \ \ \forall k. $$

Open parameters to decide: `h` (proposed `1e-4`), the trigger threshold on
`N` (proposed `> 2`), and how pure-phase `(Comp, μ)` rows are weighted (if at
all) relative to solution-phase endmember rows in the stage-B least-squares
fit.
