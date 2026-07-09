/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : MAGEMin
 **   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 **   Developers   : Nicolas Riel, Boris Kaus
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : nriel[at]uni-mainz.de, kaus[at]uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
#ifndef __GH_init_db_H_
#define __GH_init_db_H_

    #include "../MAGEMin.h"

    /**
        Single "gh" dataset: Stage-A Ghiorso/MELTS silicate liquid (proof of
        concept), plus the common rock-forming pure phases ported from
        xMELTS' own solid-phase database (quartz, cristobalite, tridymite,
        corundum, sillimanite, rutile, sphene, perovskite, calcite,
        aragonite, magnesite, siderite, dolomite, spurrite, muscovite, O2,
        H2O - see GH_PP_endmembers.h) plus andalusite/kyanite (not in
        xMELTS itself, ported from Theriak-Domino's Berman 1988 database,
        cross-validated against xMELTS' sillimanite - see
        GH_PP_endmembers.c); 5 solution phases: "liq", and the 4
        "sb-trivial" mineral solid solutions - olivine ("ol"), feldspar
        ("fsp"), biotite ("bi"), garnet ("g") - whose xMELTS mixing model
        is plain ideal+regular on a single site with no internal order
        parameter, matching sb's n_xeos==n_em pattern.
    **/
    typedef struct gh_datasets {
        int     ds_version;
        int     n_ox;
        int     n_pp;
        int     n_ss;
        char    ox[13][20];
        char    PP[19][20];
        char    SS[5][20];

        int     verifyPC[5];
        int     n_SS_PC[5];
        double  SS_PC_stp[5];

        double  PC_df_add;
        double  solver_switch_T;
        double  min_melt_T;

        double  inner_PGE_ite;
        double  max_n_phase;
        double  max_g_phase;
        double  max_fac;

        double  merge_value;
        double  re_in_n;

        double  obj_tol;

    } gh_dataset;

    global_variable global_variable_GH_init(   global_variable      gv,
                                                bulk_info           *z_b    );

    global_variable get_bulk_gh( global_variable gv );

#endif
