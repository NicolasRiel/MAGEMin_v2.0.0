# Synthesizing redundant solution-phase pseudocompound occurrences

## 1. Problem statement

Today, in `ss_min_LP` (`src/ss_min_function.c`), **every** active pseudocompound
`cp[i]` of a solution phase gets its own full NLopt local minimization each PGE
cycle, even when the same phase (`ph_id`) occurs `nOcc = gv.n_ss_ph[ph_id]`
times simultaneously. At high T, when a lot of melt is stable, `liq` commonly
has `nOcc` in the dozens, and — because all `nOcc` instances share the same
rotated hyperplane (`gv.gam_tot`, fixed for the whole PGE cycle) and are
refined from nearby starting points inside `restrict_SS_HyperVolume` boxes —
they routinely converge to the same (or a near-identical) point on the Gibbs
surface. Every one of those solves after the first is redundant work.

**Initial deployment scope.** This is implemented for `liq` only, first,
gated on `ph_id == "liq" && N > 2`. The mechanism itself is phase-agnostic
within gh's direct-parametrization phases (§7), but `liq` is the phase where
`nOcc` actually gets large enough to matter, so it is the only phase wired
up initially — extending the gate to other gh phases later is a one-line
change once `liq` is validated.

**What "generate pseudocompounds" means here.** `ss_min_LP` already has a
precedent for exactly this kind of shortcut: when `gv.n_ss_ph[ph_id] == 1`
(`ss_min_function.c:560-591`), after the one real minimization it generates
**two extra candidate pseudocompounds** by blending the new minimum
(`cp[i].xeos_1`) with the pre-minimization guess (`cp[i].xeos_0`) at two
fixed shift fractions, and adds each one to the phase's `Ppc` candidate pool
via `PC_function` + `SS_UPDATE_function` + `copy_to_Ppc` — **not** by writing
new `cp[]` occurrences directly. `copy_to_Ppc` appends into the *pre-allocated*
`SS_ref_db[ph_id].{comp,p,mu,xeos,G,DF,info}_Ppc[...]` arrays (capacity
`n_Ppc`, allocated once at database init), which is exactly the candidate
pool that `run_LP`'s existing swap logic already reads from to decide which
pseudocompounds become active `cp[]` occurrences. This document's "stage C"
(§6) is a **generalization of that exact existing block** — same call
sequence, same target arrays, same allocation-free pattern — from a fixed 2
extra points to `2N` points chosen to lie on a refined Gibbs hyperplane,
evenly distributed around the minimum, instead of a fixed shift array.

The composite method restructures one PGE cycle's `ss_min_LP` → `run_LP`
sequence into three stages:

**(A) Minimize once per unique phase, not once per occurrence.** For each
solution phase present this cycle, run exactly **one** real NLopt solve
(`rotate_hyperplane` + `restrict_SS_HyperVolume` + `NLopt_opt[ph_id]`),
regardless of how many `cp[i]` occurrences it currently has, then store it
via `copy_to_Ppc` exactly as today.

**(B) Re-estimate Gamma by least squares across all the minimized phases.**
Reusing the `Comp`/`mu_Ppc` data that stage A's own `copy_to_Ppc` call just
wrote, fit a refined `Gamma` by least squares — sharper than whatever
`gv.gam_tot` the cycle started with, since it is informed by every phase's
minimum simultaneously.

**(C) Generate `2N` synthetic pseudocompounds per phase, analytically, with
no second minimization**, and add them to the `Ppc` pool exactly the way the
existing `gv.n_ss_ph[ph_id] == 1` block does — `PC_function` +
`SS_UPDATE_function` + `copy_to_Ppc` — except the xeos fed into
`SS_ref_db[ph_id].iguess` are chosen to lie on (a linearization of) the
newly-estimated `Gamma` hyperplane, spread as evenly as possible around the
stage-A minimum (via antipodal `±` pairs, §5), instead of a fixed
shift-blend with the pre-minimization guess.

Only then does the cycle proceed into `run_LP` as today.

## 2. Notation

- `ph_id` — a solution phase (e.g. `liq`), restricted throughout to gh's
  **direct `p = xeos`** phases (see §7 scoping note): `n_xeos = n_em`, and
  composition is exactly linear, `Comp(x) = C\,x`, where `C = SS_ref_db[ph_id].Comp`
  (`n_em × n_ox`, a fixed database table, not recomputed).
- `N = nOcc` — number of simultaneous occurrences of `ph_id` this cycle
  (`gv.n_ss_ph[ph_id]`); the shortcut fires for `ph_id == "liq" && N > 2`
  (§1) — same trigger style as the existing `gv.n_ss_ph[ph_id] == 1` check
  it generalizes, just phase-gated and with a higher threshold.
- `2N` — the number of synthetic pseudocompounds generated per cycle for
  the gated phase (§5-§6): not a replacement set sized to match `N`, but a
  deliberately denser, evenly-distributed set of *candidates* added to the
  `Ppc` pool for `run_LP`'s existing swap logic to choose from.
- `n_ox = gv.len_ox` (≤ `n_ox_all = 16`, the fixed bound already used
  throughout `MAGEMin.h` for oxide-indexed arrays).
- `x* ∈ R^{n_xeos}` — `SS_ref_db[ph_id].xeos` (or `.iguess`, same thing
  post-solve) after the single stage-A NLopt solve for `ph_id`.
- `G(x*)`, `g* = ∇G(x*)` — `SS_ref_db[ph_id].df` and
  `SS_ref_db[ph_id].dfx[]` (the objective value and gradient, already
  computed as a byproduct of the stage-A NLopt call — no extra evaluation).
- `μ*_j`, `j=1..n_em` — the *unrotated* endmember chemical potentials at
  `x*`, i.e. exactly `SS_ref_db[ph_id].mu_Ppc[m][j]` for the pool slot `m`
  that stage A's own `copy_to_Ppc` call just filled (`copy_to_Ppc` already
  calls `non_rot_hyperplane` + a direct objective re-evaluation to get this
  — reused as-is, not recomputed).
- `Γ_old = gv.gam_tot` — the Gamma the cycle started with.
- `Γ_new` — the refined estimate produced by stage B.
- `x_1^{old}, ..., x_N^{old}` — the **previous** cycle's converged `xeos`
  for the `N` occurrences of `ph_id` (`cp[i].xeos` for the `N` `cp[]`
  entries sharing this `ph_id`, read *before* this cycle overwrites them),
  `x̄ = (1/N)Σ_k x_k^{old}`, `e_k = x_k^{old} - x̄` (so `Σ_k e_k = 0` exactly).
- `h` — target step size for the synthetic spread, in `xeos` units
  (proposed ≈ `1e-4`).

## 3. Stage A: minimize each unique phase once

Restructure the per-`cp[i]` loop in `ss_min_LP` into a per-unique-`ph_id`
loop. For each distinct phase present, run the existing block exactly once
(unchanged — same `rotate_hyperplane`, `restrict_SS_HyperVolume`,
`NLopt_opt[ph_id]`, then `PC_function`+`SS_UPDATE_function`+`copy_to_Ppc`
for the real minimum, matching `ss_min_function.c:496-618` as it stands
today), using **the first occurrence encountered in `cp[]` order** for
`liq` as the seed `iguess` (simplest, no pre-scan; confirmed choice for
this iteration). Record the `Ppc` pool slot index `m*` that this call wrote
to (`SS_ref_db[ph_id].id_Ppc - 1` right after the call) so stage B/C can
read `mu_Ppc[m*]`/`xeos_Ppc[m*]` back out — no separate storage needed.

**Deferred generalization (not designed here): solvus phases.** A phase
with a solvus (e.g. two compositionally-distinct stable instances of the
same phase, feldspar unmixing, exsolved pyroxenes) can have `N` occurrences
that are genuinely *not* all near the same minimum — picking one seed and
minimizing once would silently collapse two real physical phases into one.
Handling this properly needs a clustering/splitting step in stage A (detect
distinct occurrence groups by `xeos` separation, minimize once per group)
before the rest of the method applies per group. Out of scope for this
iteration (`liq` has no solvus); noted here so stage A's "one real
minimization" assumption isn't silently misapplied later to a phase that
needs it.

## 4. Stage B: weighted least-squares refinement of Gamma (allocation-free)

Rank is not expected to be a practical concern: `liq` alone already
contributes up to `n_em ≈ n_ox` rows, and every other minimized solution
phase plus every already-minimized pure phase (`pp_min_function`) adds more
— in practice always comfortably overdetermined. What matters more is
**weighting** the rows, since not every endmember's `μ_j` is an equally
trustworthy read on the true hyperplane:

- An endmember with `p_j ≈ 0` in its phase's minimized composition has a
  `μ_j` dominated by the `R\,T\,\log(p_j)` configurational term, which blows
  up/becomes numerically unreliable as `p_j → 0` — it should count for
  little.
- A phase present in only a trace amount is a less robust anchor for the
  hyperplane than one that is genuinely abundant in the current assemblage.

So each endmember row is weighted by its actual molar abundance,
`w_j = p_j \cdot z_{em,j}` (already computed — this is exactly the quantity
`SS_UPDATE_function`/`copy_to_Ppc` use to build `ss_comp` from `p`, reused
here, not recomputed), giving the weighted normal equations

$$ \Gamma_{\text{new}} \;=\; (M^\top W M)^{-1} M^\top W \mu, \qquad W = \text{diag}(w_j). $$

**Never form or store `M`.** Since `M` only enters through `MᵀWM` (`n_ox ×
n_ox`) and `MᵀWμ` (`n_ox`), accumulate both as running weighted rank-1 sums
directly inside the existing per-phase loop, immediately after each
phase's `copy_to_Ppc` call in stage A:

```
for j in 0..n_em(ph_id)-1:
    w = p[j] * z_em[j]                      // molar abundance of endmember j at x*
    for a in 0..n_ox-1:
        Mtmu[a] += w * C[j][a] * mu_Ppc[m*][j]
        for b in 0..n_ox-1:
            MtM[a][b] += w * C[j][a] * C[j][b]
```

`MtM[n_ox_all][n_ox_all]` and `Mtmu[n_ox_all]` are **fixed-size stack
arrays** (`n_ox_all = 16` is already the codebase-wide bound for
oxide-indexed arrays — see `MAGEMin.h:42`), zeroed once at the top of the
cycle. Pure phases already minimized this cycle by `pp_min_function`
contribute their own single `(Comp, μ)` row the same way with `w = 1`
(a pure phase's `μ` has no configurational-entropy term at all, so it needs
no abundance discount).

> This weighting scheme is my best reconstruction of "the weight approach"
> from our earlier discussion of this idea — please confirm it matches what
> you had in mind, or correct it, before it's implemented.

Solving the final small `n_ox × n_ox` system for `Γ_new` is done in place
with a fixed-size Gauss elimination (no malloc) — the same "small, fixed
dimension, in-place, no allocation" style already used for the joint
Newton solves in `GH_cpx_solve_s_from`/`GH_spn_solve_s_from`
(`gh_objective_functions.c`), just `n_ox × n_ox` (≤ 16×16) instead of 2×2/3×3.
Given rank is not expected to be an issue in practice (see above), no
rank-deficiency fallback is included in this design; if it were ever
observed, the simplest safe response would be to keep `Γ_new = Γ_old` for
that cycle rather than trust an ill-conditioned solve.

## 5. Stage C, part 1: the linearized Gamma hyperplane around `x*`

The rotated (driving-force) energy relative to `Γ_new` is

$$ \Delta G(x) \;=\; G(x) \;-\; \Gamma_{\text{new}} \cdot \text{Comp}(x) \;=\; G(x) - \Gamma_{\text{new}}\cdot C\,x. $$

`Comp` is exactly linear in `x` (§2 scoping), so linearizing only `G` around
`x*` gives, using quantities already in hand from §3/§4:

$$ \Delta G(x^{*}+\delta) \;\approx\; \underbrace{\big[G(x^{*}) - \Gamma_{\text{new}}\cdot C\,x^{*}\big]}_{\Delta G_0} \;+\; \underbrace{\big[g^{*} - C^\top \Gamma_{\text{new}}\big]}_{g_{\text{eff}}} \cdot\ \delta. $$

Both are closed-form, `O(n_xeos · n_ox)` work, no new objective evaluation.

**Base point.** The minimum-norm `δ` landing `x*` exactly on the linearized
hyperplane:

$$ \delta_0 \;=\; -\frac{\Delta G_0}{\lVert g_{\text{eff}}\rVert^2}\, g_{\text{eff}}, \qquad x^{**} \;=\; x^{*} + \delta_0. $$

Small whenever `Γ_new ≈ Γ_old` (the converging/common case).

**Spread directions.** Any `v` with `g_eff · v = 0` keeps `ΔG ≈ 0` to first
order from `x**` — the hyperplane's tangent space at `x**` is the null
space of `g_eff`. Reusing the previous occurrence set's shape (§2's `e_k`,
`k = 1..N`) by projecting through that null space,

$$ P \;=\; I \;-\; \frac{g_{\text{eff}}\,g_{\text{eff}}^\top}{\lVert g_{\text{eff}}\rVert^2}, \qquad v_k \;=\; P\,e_k, $$

preserves the zero-sum property exactly (`P` linear ⇒ `Σv_k = P·Σe_k = 0`)
while making each `v_k` satisfy `g_eff · v_k = 0` exactly.

**Doubling to `2N`, evenly, via antipodal pairs.** Rather than emitting the
`N` points `x** + s\,v_k` alone (which inherits whatever asymmetry the
previous occurrence set happened to have), emit **both** `v_k` and `-v_k`
for each `k = 1..N`:

$$ \boxed{\; x_k^{\pm} \;=\; x^{**} \;\pm\; s\,P\,e_k, \qquad k = 1, \dots, N, \;} $$

giving `2N` points total. This is "as evenly distributed as possible around
the minimized point" in the cheapest exact sense available: every direction
sampled is immediately balanced by its own antipode, so the set is
point-symmetric about `x**` by construction — not merely zero-sum on
average, but exactly so pair-by-pair (`x_k^{+} + x_k^{-} = 2x^{**}` for
every `k`), regardless of how skewed or collinear the original `e_k`
directions were. `s = h / max_k‖P e_k‖` (§2's target step, same
normalization idea as the existing shift-array logic).

`g_eff`, `e_k`, `v_k` are fixed-size stack arrays of length `n_xeos` (bounded
by a small new constant, e.g. `GH_LP_MAX_XEOS`, sized to the largest `n_em`
among gh solution phases — currently 13, for `liq`); no allocation.

## 6. Stage C, part 2: adding the points to the Ppc pool (reusing existing code verbatim)

For each of the `2N` points `x_k^{\pm}` from §5, run **exactly** the same
4-line sequence the existing `gv.n_ss_ph[ph_id] == 1` block already uses
(`ss_min_function.c:566-589`), substituting each point for that block's
shift-blend:

```
for (int k = 0; k < cp[i].n_xeos; k++) {
    SS_ref_db[ph_id].iguess[k] = x_k[k];        // clamped to bounds_ref first
}
SS_ref_db[ph_id] = PC_function(        gv, PC_read, SS_ref_db[ph_id], z_b, ph_id);
SS_ref_db[ph_id] = SS_UPDATE_function( gv, SS_ref_db[ph_id], z_b, gv.SS_list[ph_id]);
if (SS_ref_db[ph_id].sf_ok == 1){
    copy_to_Ppc(0, 1, ph_id, gv, SS_objective, SS_ref_db);
}
```

`PC_function` (`TC_database/objective_functions.c:18222`) is the existing
"evaluate the real objective once at `iguess`, no NLopt search" primitive —
this *is* the "call the real objective directly" step from earlier drafts
of this document, just via the already-existing wrapper instead of a new
one. `copy_to_Ppc` writes into the pre-allocated `Ppc` pool exactly as
today. Nothing here allocates memory: `iguess` is the phase's existing
scratch buffer, and the `Ppc` arrays are the existing pre-sized pool.

This also gets the "parallel to the hyperplane but slightly higher" G
behavior for free and gives an honest check on the §5 linearization: since
`PC_function` evaluates the *real* `G(x_k)`, the true `ΔG(x_k) = G(x_k) -
Γ_new·Comp(x_k)` can be compared against the assumed `≈0` — a nonzero
result by more than `O(h²)` would flag that `h` should shrink or that
`Γ_new` moved too far this cycle to trust the linearization.

## 7. Scoping note: gh direct-parametrization phases only

§5's `Comp(x) = C\,x` exactness relies on gh's "direct `p = xeos`"
architecture (single `Σx=1` equality constraint, no reduced-basis site
fraction map) — true for every phase in the gh database. It is **not**
generally true for `tc`'s reduced-basis phases, where `p(xeos)` can be a
nonlinear map (order parameters, site-fraction algebra), making
`Comp(xeos)` nonlinear too. Extending this design to `tc` would require an
extra chain-rule factor (`dComp/dxeos|_{x*} = C · dp/dxeos|_{x*}`) in `g_eff`
— straightforward in principle but out of scope here; this document targets
gh phases (`liq` foremost, but the mechanism is phase-agnostic within gh).

## 8. Edge cases and guards

- **`‖g_eff‖ ≈ 0`.** `x*` is already stationary for `Γ_new` too — set
  `δ_0 = 0`, `P = I` (skip the projection).
- **Degenerate previous spread.** If `max_k‖P e_k‖` is at/below machine
  epsilon, fall back to a fixed zero-sum direction set spanning `null(g_eff)`
  (e.g. an orthonormal basis of it, or a regular `(N-1)`-simplex projected
  through `P`), scaled to `h`.
- **`N ≤ 2` or `ph_id != "liq"`.** No shortcut — run the real per-occurrence
  solve as today (exactly today's behavior, unchanged). Stage C is gated
  strictly on `ph_id == "liq" && N > 2` for this initial deployment (§1).
- **Bounds.** Clamp every `x_k^{\pm}` into `SS_ref_db[ph_id].bounds_ref`
  before `PC_function`, same as the existing block implicitly relies on
  `restrict_SS_HyperVolume` for the real solve.
- **Trigger condition.** Gate stage C's synthesis on `ph_id == "liq" && N >
  2` **and** not the first PGE cycle for this bulk/step (no previous
  `x_k^{old}` set exists yet — that cycle runs `N` real solves once, as
  today). Stages A and B (minimize-once-per-unique-phase, refine Gamma) are
  *not* phase-gated — they apply to every phase regardless of `N`, since
  they are a pure win even for `N = 1` phases (fewer NLopt calls overall,
  better-conditioned Gamma); only stage C's `2N`-point synthesis is
  currently restricted to `liq`.

## 9. Where this hooks into the existing code

- **`ss_min_LP`** (`src/ss_min_function.c`): restructure the per-`cp[i]`
  loop (currently lines 482-649) into a per-unique-`ph_id` loop for stage A
  (§3, reusing lines 496-618 essentially unchanged); accumulate `MtM`/`Mtmu`
  (§4) inline right after each phase's stage-A `copy_to_Ppc` call; after the
  per-phase loop, solve for `Γ_new` once; then, only for the phase where
  `strcmp(gv.SS_list[ph_id], "liq") == 0 && gv.n_ss_ph[ph_id] > 2`,
  synthesize the `2N` points (§5) and add them via the exact existing
  4-line pattern (§6) — this **sits alongside, and for `liq` replaces**, the
  current `if (gv.n_ss_ph[ph_id] == 1) { ... 2 fixed shift-blend points ...
  }` block at lines 560-591 (that block's own logic still applies unchanged
  to every other phase, and to `liq` whenever `N == 1`).
- **`run_LP`** (`src/PGE_function.c`): unchanged — it continues to consume
  whatever the `Ppc` pool holds and swap/update `cp[]` as today. `Γ_new`
  from stage B can also seed `run_LP`'s own Gamma update instead of being
  discarded after stage C.

## 10. Summary

$$ \Delta G_0 = G(x^{*}) - \Gamma_{\text{new}}\cdot C x^{*}, \qquad g_{\text{eff}} = g^{*} - C^\top \Gamma_{\text{new}} $$

$$ \delta_0 = -\frac{\Delta G_0}{\lVert g_{\text{eff}}\rVert^2} g_{\text{eff}}, \qquad P = I - \frac{g_{\text{eff}} g_{\text{eff}}^\top}{\lVert g_{\text{eff}}\rVert^2} $$

$$ x_k^{\pm} = x^{*} + \delta_0 \;\pm\; \frac{h}{\max_j\lVert P e_j\rVert}\, P\,(x_k^{\text{old}} - \bar x), \qquad k = 1,\dots,N $$

(`2N` points total), added to the `Ppc` pool via the existing `PC_function`
+ `SS_UPDATE_function` + `copy_to_Ppc` sequence (§6), gated on `ph_id ==
"liq" && N > 2` (§1, §8), with the exact guarantee (to the §5 linearization)

$$ \frac{1}{2N}\sum_{k=1}^{N}\big(x_k^{+} + x_k^{-}\big) = x^{*} + \delta_0 \quad\text{and}\quad \Delta G(x_k^{\pm}) \approx 0 \ \ \forall k, $$

the average holding not just in aggregate but exactly pair-by-pair
(`x_k^+ + x_k^- = 2(x^*+\delta_0)`), which is what makes the set "as evenly
distributed as possible" around the minimized point rather than merely
zero-sum on average.

Everything is allocation-free: `MtM`/`Mtmu` are fixed `n_ox_all`-sized
stack arrays, `g_eff`/`e_k`/`v_k` are fixed-size stack arrays bounded by a
new `GH_LP_MAX_XEOS` constant, and the `Ppc` pool / `iguess` buffers are the
already-existing pre-allocated per-phase arrays.

**Runtime-tunable vs. compile-time.** `h` (spread step) and the `N`
threshold become new `gv` fields (e.g. `gv.gh_liq_pc_synth_h`, default
`1e-4`; `gv.gh_liq_pc_synth_threshold`, default `2`) — consistent with the
existing convention for tunable numerical knobs (`gv.obj_tol`,
`gv.box_size_mode_LP`, `gv.mSS_df_min_add`, ...). `GH_LP_MAX_XEOS` stays a
compile-time `#define`, not a `gv` field: it sizes fixed stack arrays
(`g_eff`, `e_k`, `v_k`), and the whole point of this design is that those
are allocation-free — a runtime-tunable bound would force either a VLA or a
malloc, defeating that requirement.

Weighted-least-squares row weighting (§4) is included as pure-phase rows
with `w = 1`; the `liq`-only, `N > 2` gate and the first-occurrence seed
choice (§3) are fixed by this iteration of the design, not open.
