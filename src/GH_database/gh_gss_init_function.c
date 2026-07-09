/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : MAGEMin
 **   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 **   Developers   : Nicolas Riel, Boris Kaus
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : nriel[at]uni-mainz.de, kaus[at]uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
#include <stdio.h>
#include <string.h>

#include "../MAGEMin.h"
#include "gh_gss_init_function.h"

/**
    Metadata for the Stage-A Ghiorso/MELTS liquid ("liq"), following the
    same formulation as "sb"'s solution phases (e.g. G_SS_sb11_ol_init_function):
    n_xeos == n_em (the NLopt vector IS the endmember mole-fraction vector
    directly, p[i] = x[i], no reduced/rotated basis), with the Sigma(x)=1
    closure enforced via a single NLopt equality constraint rather than
    algebraic substitution - this matches "sb"'s LP-only solving path,
    which is what "gh" also uses (see SetupDatabase).
    - 13 real basis endmembers (n_em = n_xeos = 13; SO3 and Cl2O-1 were
      dropped - no real reference-state data exists for them anywhere in
      xMELTS - and Ni/Co/P2O5/F were already excluded, see GH_endmembers.c).
    - n_sf = n_em: the "site fractions" here are simply the 13 mole
      fractions themselves (single mixing "site"), used only for the
      generic positivity check in ss_min_function.c's SS_UPDATE_function.
    - symmetric regular solution (Margules W only, no Van Laar), n_w = C(13,2) = 78.
    - n_cat = 0: no site-mixing bookkeeping (C[]/N[]) is populated - like
      "sb", these fields are metadata-only and never read at runtime.
*/
SS_ref G_SS_gh_liq_init_function(SS_ref SS_ref_db, global_variable gv){
    SS_ref_db.is_liq    = 1;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 0;
    SS_ref_db.n_xeos    = 13;
    SS_ref_db.n_em      = 13;
    SS_ref_db.n_sf      = 13;
    SS_ref_db.n_w       = 78;

    return SS_ref_db;
}

/**
    Olivine (Fo-Fa), "sb-trivial" tier: direct p=x, symmetric regular
    solution (single W term, WH1MGFE=WH2MGFE=5075 J in xMELTS'
    sources/olivine.c, i.e. the Mg-Fe binary happens to be symmetric even
    though the full 6-component model is asymmetric/van-Laar elsewhere).
    Real xMELTS olivine also carries 4 internal M1/M2 Fe-Mg-Mn-Co-Ni
    order parameters (order-disorder machinery, "needs order-param
    machinery" tier); this Fo-Fa reduction deliberately drops them and
    uses the disordered-limit treatment instead, exactly like sb's own
    "ol" (2 equivalent M-sites, Sconfig multiplicity 2, no order
    parameter) - matching the "same structure as sb" scope the user chose.
*/
SS_ref G_SS_gh_ol_init_function(SS_ref SS_ref_db, global_variable gv){
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 2;
    SS_ref_db.n_xeos    = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_w       = 1;

    return SS_ref_db;
}

/**
    Feldspar (Ab-An-San), Elkins & Grove (1990) ternary asymmetric
    regular solution (xMELTS sources/feldspar.c) - single framework site,
    no order parameter (Al-Si and albite/sanidine order-disorder
    corrections that xMELTS itself applies in sources/gibbs.c are not
    reproduced, same simplification already made for the SiO2 polymorphs
    and for these same endmembers' own standard states, see
    GH_PP_endmembers.h). Direct p=x, n_xeos=n_em=3.
*/
SS_ref G_SS_gh_fsp_init_function(SS_ref SS_ref_db, global_variable gv){
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 1;
    SS_ref_db.n_xeos    = 3;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_sf      = 3;
    SS_ref_db.n_w       = 7;

    return SS_ref_db;
}

/**
    Biotite (Ann-Phl), Sack & Ghiorso binary regular solution (xMELTS
    sources/biotite.c) - single octahedral site, symmetric (one W term),
    no site-multiplicity factor on the ideal entropy in xMELTS' own
    calibration (ported as-is). Direct p=x, n_xeos=n_em=2.
*/
SS_ref G_SS_gh_bi_init_function(SS_ref SS_ref_db, global_variable gv){
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 3;
    SS_ref_db.n_xeos    = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_w       = 1;

    return SS_ref_db;
}

/**
    Garnet (Gr-Py-Alm), Berman & Koziol (1991) ternary asymmetric
    regular solution (xMELTS sources/garnet.c) - single dodecahedral
    X-site (multiplicity 3), no order parameter. Direct p=x, n_xeos=n_em=3.
*/
SS_ref G_SS_gh_g_init_function(SS_ref SS_ref_db, global_variable gv){
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 3;
    SS_ref_db.n_xeos    = 3;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_sf      = 3;
    SS_ref_db.n_w       = 7;

    return SS_ref_db;
}

void GH_SS_init(            SS_init_type        *SS_init,
                            global_variable      gv              ){

    for (int iss = 0; iss < gv.len_ss; iss++){
        if (strcmp( gv.SS_list[iss], "liq") == 0 ){
            SS_init[iss]  = G_SS_gh_liq_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "ol") == 0 ){
            SS_init[iss]  = G_SS_gh_ol_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "fsp") == 0 ){
            SS_init[iss]  = G_SS_gh_fsp_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "bi") == 0 ){
            SS_init[iss]  = G_SS_gh_bi_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "g") == 0 ){
            SS_init[iss]  = G_SS_gh_g_init_function;
        }
        else{
            printf("\nsolid solution '%s' is not in the 'gh' database, cannot be initiated\n", gv.SS_list[iss]);
        }
    }
}
