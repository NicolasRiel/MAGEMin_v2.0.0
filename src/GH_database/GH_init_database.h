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
        GH_PP_endmembers.c); 10 solution phases: "liq", the 4 "sb-trivial"
        mineral solid solutions - olivine ("ol"), feldspar ("fsp"), biotite
        ("bi"), garnet ("g") - whose xMELTS mixing model is plain
        ideal+regular on a single site with no internal order parameter,
        matching sb's n_xeos==n_em pattern, plus hornblende ("hb",
        pargasite-ferropargasite-magnesiohastingsite) and leucite ("lc",
        leucite-analcime-na-leucite), both NS=0 (no internal order
        parameter) in real xMELTS too, so they follow the exact same
        direct-p=x pattern despite their richer (multi-site or coupled)
        entropy-of-mixing formulas, plus melilite ("mel", NS=1, the first
        phase needing an embedded order-parameter Newton/bisection solve),
        cummingtonite ("cum", cummingtonite-grunerite, NS=2, a genuine 2D
        embedded order-parameter solve for M4/M1+M3/M2 site Fe-Mg
        ordering), and spinel ("spn", chromite-hercynite-magnetite-spinel-
        ulvospinel, NR=4/NS=3, gh's first 3D embedded order-parameter
        solve, using nonlinear Gauss-Seidel over robust 1D bisections
        rather than a joint Newton step - see gh_objective_functions.c;
        clinopyroxene/orthopyroxene ("cpx"/"opx", Di-Cen-Hed-CaTs(Al)-
        CaTs(Fe3+)-Ess-Jd, NR=6); and the real mixed H2O-CO2 fluid
        ("fluid", Pitzer & Sterner 1994, from xMELTS' sources/fluid.c) -
        gh's first phase that is a genuine real-gas EOS mixture rather
        than a discrete-endmember Margules solution, see
        GH_pitzer_sterner_mix_G's header comment (GH_fluid_eos.h).
    **/
    typedef struct gh_datasets {
        int     ds_version;
        int     n_ox;
        int     n_pp;
        int     n_ss;
        char    ox[13][20];
        char    PP[22][20];
        char    SS[16][20];

        int     verifyPC[16];
        int     n_SS_PC[16];
        double  SS_PC_stp[16];

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
