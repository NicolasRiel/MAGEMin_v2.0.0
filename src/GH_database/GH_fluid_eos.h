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

#endif
