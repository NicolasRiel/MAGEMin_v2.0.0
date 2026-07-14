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
#ifndef __GH_GEM_FUNCTION_H_
#define __GH_GEM_FUNCTION_H_

#include "../gem_function.h"

/** Real gv.EM_database (0=xMELTS/1=rMELTS/2=pMELTS), set once by
    GH_SS_objective_init_function (gh_objective_functions.c) - see
    GH_gem_function.c's own header comment on why this can't just be
    GH_G_EM_function's own (misleadingly-named, wrongly-fed) parameter. */
extern int GH_actual_EM_database;

/** Disambiguates gh's two distinct callers of name=="H2O": real xMELTS
    itself uses two genuinely different gibbs.c branches - "H2O" (the
    liquid's own basis species, fixed 1-bar Haar reference + Robie ideal-
    gas + Oaks-Lange V(T,P)) vs "water" (a real standalone pure phase,
    whaar() at the actual P capped at 10000 bar + wdh78() above that cap,
    no Robie/Oaks-Lange wrapper - see gibbs.c line ~2326). gh only has one
    endmember table entry named "H2O" (shared by the liquid's own
    construction in G_SS_gh_liq_function/G_SS_gh_fluid_function AND the
    standalone gv.PP_list "H2O" pure phase, since gh's PP list was never
    renamed to "water" - renaming was ruled out because toolkit.c/
    dump_function.c key many unrelated, cross-research-group checks off
    the literal PP_list string "H2O"). Default (0) = pure "water" phase
    semantics (used by the standalone PP list AND toolkit.c's system_aH2O
    activity calc, which both conceptually want the real pure-water
    reference state, not the liquid's). Set to 1 only for the narrow
    window around G_SS_gh_liq_function's/G_SS_gh_fluid_function's own
    get_em_data(...,"H2O",...) calls (gh_gss_function.c), which need the
    liquid basis-species formula instead. Found + fixed 2026-07-15 during
    the 10x10 grid sweep - see [[gh-spn-liq-gbase-verification]]. */
extern int GH_H2O_liquid_context;

PP_ref GH_G_EM_function(   int          EM_database,
                            int          len_ox,
                            int         *id,
                            double      *bulk_rock,
                            double      *apo,
                            double       P,
                            double       T,
                            char        *name,
                            char        *state          );

#endif
