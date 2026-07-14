/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : MAGEMin
 **   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 **   Developers   : Nicolas Riel, Boris Kaus
 **   Contributors : Nickolas B. Moccetti, Dominguez, H., Assunção J., Green E., Berlie N., and Rummel L.
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : nriel[at]uni-mainz.de, kaus[at]uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
#ifndef __GH_FLUID_EOS_H_
#define __GH_FLUID_EOS_H_

/**
    Pitzer & Sterner (1994) real-gas equation of state for pure H2O or
    pure CO2, ported from xMELTS-master/sources/fluid.c's fluidPhase()
    (specialized to the pure-component cases x=[1,0] and x=[0,1]).

    References (as cited in the original xMELTS source):
      H2O: Pitzer KS and Sterner SM (1994) Equations of state valid
           continuously from zero to extreme pressures for H2O and CO2.
           J Chem Phys 101: 3111-6
      CO2: Sterner SM and Pitzer KS (1994) An equation of state for
           carbon dioxide valid from zero to extreme pressures.
           Contr Mineral Petrol 117: 362-74

    Returns the standard-state Gibbs energy (J/mol), referenced to the
    elements at 298.15 K/1 bar (Berman 1988 convention - the same
    convention MAGEMin's other endmembers use), at the given T (K) and
    P (bar).
*/
double GH_pitzer_sterner_G(int is_H2O, double T, double Pbar);

/**
    Pure-H2O "h - T*s" recomputation of the Pitzer-Sterner EOS above,
    needed for the two real gibbs.c contexts (pMELTS' standalone "water"
    phase, and - less obviously - pMELTS' own liquid "H2O" basis species,
    via a downstream Birch-Murnaghan-block fallback that recomputes gl
    this way) that use h-T*s instead of GH_pitzer_sterner_G's raw g.
    These disagree by tens of kJ because fluid.c's own reference-state
    shifts for g/h/s are three independently-calibrated constants, not a
    thermodynamically consistent triple. See GH_fluid_eos.c for the full
    derivation (avoids porting fluidPhase()'s dA/dT machinery entirely -
    verified exact against real fluidPhase() output at two (T,P) points).
*/
double GH_pitzer_sterner_H2O_hTs_G(double T, double Pbar);

/**
    Real, mixed H2O-CO2 fluid: the general (composition-dependent) branch
    of fluidPhase() that GH_pitzer_sterner_G() above deliberately dropped.
    x_h2o in [0,1] (x_co2 = 1-x_h2o). Returns total molar Gibbs energy
    G(x,T,P) of the mixture (J/mol, same Berman 1988 reference scale as
    GH_pitzer_sterner_G - reduces to it exactly at x_h2o=1 or x_h2o=0),
    and writes dG/dx_h2o (analytic, at fixed T,P) to *dGdx_h2o.

    The mixing physics ported here (faithful to fluidPhase(), not a
    re-derivation): a single shared fluid density solved from
    composition-linearly-mixed Pitzer-Sterner virial coefficients
    c[i] = x_h2o*c_H2O[i] + x_co2*c_CO2[i], seeded by a van der Waals
    one-fluid Redlich-Kwong initial guess (the only place a genuine H2O-
    CO2 cross term "aMix" enters); the ideal-gas contribution is the
    linear combination x_h2o*Ai_H2O(T,rho) + x_co2*Ai_CO2(T,rho) evaluated
    at that one shared density. Notably (confirmed against fluidPhase()
    itself, not an omission here) there is no separate ideal-mixing
    entropy term (-RT*sum(x*ln x)) anywhere in real MELTS' own fluid
    model - the function's own header comment even flags the c[i] mixing
    rule as "an initial guess at the functional form". Ported as-is.

    dG/dx_h2o is analytic via the envelope theorem: since P is the fixed
    input (rho is solved so the EOS reproduces it), d(A+P/rho)/dx at the
    solution density reduces to the partial dA/dx at FIXED rho (the
    dA/drho * drho/dx cross term vanishes because dA/drho equals P/rho^2
    exactly at the converged density) - mechanically identical in
    structure to fluidPhase()'s own d/dT derivatives (dArdt etc.), just
    substituting c_H2O[i]-c_CO2[i] for dc/dt and dropping the d/drho
    terms, since rho is held fixed rather than T.
*/
double GH_pitzer_sterner_mix_G(double x_h2o, double T, double Pbar, double *dGdx_h2o);

/**
    Haar (1984) H2O equation of state, ported from
    xMELTS-master/sources/water.c's whaar(P, T, ...). Takes an explicit
    pressure since real gibbs.c calls whaar() at TWO different pressures
    depending on which standard state is being built:
      - the liquid's own H2O basis species always uses a fixed 1-bar
        reference (call with P=1.0) - this is what real gibbs.c's liquid
        "H2O" standard state actually calls (NOT GH_pitzer_sterner_G/
        fluidPhase() above, which is a separate model used only for
        pMELTS' H2O and for gh's own "fl" mixed-fluid solution phase) -
        discovered 2026-07-15 while live-verifying gh's liq H2O/CO2 gbase
        against real MELTS output: GH_pitzer_sterner_G had never actually
        been validated against a real reference for this use and was off
        by ~47-50 kJ.
      - the standalone "water" pure phase (GH_gem_function.c's own
        "water" branch) calls whaar() at min(actual P, 10000 bar) - see
        GH_wdh78_G below for the correction needed above that cap. Found
        missing entirely (gh used GH_pitzer_sterner_G there instead,
        off by several kJ) during the 2026-07-15 grid sweep.
    See [[gh-multicalibration-xmelts-rmelts-pmelts]].

    Value-only port (drops all T-derivative and 2nd-rho-derivative
    bookkeeping real whaar() computes only for its Cp/V/dV outputs, which
    gh's standard-state usage never needs) - verified bit-exact (6 decimal
    places) against the real whaar() function at multiple T via a
    standalone harness linked against libMELTSbatch.a.

    Returns G (J/mol) at (T, P), Berman (1988) elements reference
    scale, already including the Berman correction shift real water.c
    applies for MODE_xMELTS/MODE__MELTS/MODE__MELTSandCO2/
    MODE__MELTSandCO2_H2O (the only branch gh needs - xMELTS and rMELTS
    both use it).
*/
double GH_haar_H2O_G(double T, double P);

/**
    High-pressure correction for GH_haar_H2O_G, ported from water.c's own
    wdh78() ("difference in thermodynamic properties of water between P,T
    and 10kbar,T") - needed because whaar() itself is only reliable up to
    10000 bar. Real gibbs.c's standalone "water" branch calls this (added
    on top of whaar(10000,T,...)) whenever the actual pressure exceeds
    10000 bar; not needed anywhere P<=10000 bar (liq's own H2O never
    triggers this, since it always calls whaar at the fixed 1-bar
    reference regardless of the real T,P).

    Value-only port of a simple 15-term (T,P) polynomial - no iterative
    solve involved, unlike whaar() itself.
*/
double GH_wdh78_G(double T, double P);

/**
    Duan (1992) pure-CO2 equation of state at a fixed reference pressure
    of 1 bar, ported from xMELTS-master/sources/fluidPhase.c's
    duanCO2Driver()/idealGasCO2() (the functions real gibbs.c's
    propertiesOfPureCO2() calls), specialized to the pure endpoint
    (x_CO2=1) and value-only (no derivatives - only the fugacity
    coefficient phi(T,p) and ideal-gas h0(T)/s0(T) values are needed for
    G, not Cp/V/their derivatives).

    Verified bit-exact (6 decimal places) against real
    propertiesOfPureCO2() via a standalone harness.

    Returns G (J/mol) at (T, p=1bar), Berman (1988) elements reference
    scale (matches propertiesOfPureCO2's own *g output, which is already
    on this scale via its ideal-gas reference constants).
*/
double GH_duan_CO2_G(double T);

#endif
