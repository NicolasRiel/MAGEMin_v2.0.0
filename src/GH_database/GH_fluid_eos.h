/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : MAGEMin
 **   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 **   Developers   : Nicolas Riel, Boris Kaus
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : nriel[at]uni-mainz.de, kaus[at]uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
#ifndef __GH_FLUID_EOS_H_
#define __GH_FLUID_EOS_H_

/**
    Pitzer & Sterner (1994) real-gas equation of state for pure H2O or
    pure CO2, ported from xMELTS-master/sources/fluid.c's fluidPhase()
    (specialized to the pure-component cases x=[1,0] and x=[0,1], since
    that is all "gh" needs the standard state for - the mixture/mixing
    terms of the original coupled H2O-CO2 model are not required and are
    not ported).

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

#endif
