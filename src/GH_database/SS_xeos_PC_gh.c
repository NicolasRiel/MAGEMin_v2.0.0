/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : MAGEMin
 **   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 **   Developers   : Nicolas Riel, Boris Kaus
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : nriel[at]uni-mainz.de, kaus[at]uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
#include <string.h>

#include "../MAGEMin.h"
#include "SS_xeos_PC_gh.h"

/**
    Coarse pseudocompound composition grid for "liq" (Stage-A proof of
    concept): direct 13-dimensional mole-fraction vectors (p[i], summing
    to 1), following "sb"'s formulation (n_xeos == n_em, xeos IS p). The
    pure-SiO2 corner, two points along each of the other 12 composition
    axes (x_i = 0.3, 0.6, SiO2 making up the remainder), and one interior
    "generic andesitic" guess. This is deliberately small/coarse - unlike
    SB_database's production-quality barycentric grids - since this
    database's primary verification target is the local-minimization
    objective function itself (see project plan), not levelling accuracy.
    p order: SiO2,TiO2,Al2O3,Fe2O3,MgCr2O4,Fe2SiO4,MnSi0.5O2,Mg2SiO4,
             CaSiO3,Na2SiO3,KAlSiO4,CO2,H2O
*/
struct ss_pc gh_liq_pc_xeos[26] = {
    {{ 1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000 }},
    {{ 0.7000, 0.3000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000 }},
    {{ 0.4000, 0.6000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000 }},
    {{ 0.7000, 0.0000, 0.3000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000 }},
    {{ 0.4000, 0.0000, 0.6000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000 }},
    {{ 0.7000, 0.0000, 0.0000, 0.3000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000 }},
    {{ 0.4000, 0.0000, 0.0000, 0.6000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000 }},
    {{ 0.7000, 0.0000, 0.0000, 0.0000, 0.3000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000 }},
    {{ 0.4000, 0.0000, 0.0000, 0.0000, 0.6000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000 }},
    {{ 0.7000, 0.0000, 0.0000, 0.0000, 0.0000, 0.3000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000 }},
    {{ 0.4000, 0.0000, 0.0000, 0.0000, 0.0000, 0.6000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000 }},
    {{ 0.7000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.3000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000 }},
    {{ 0.4000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.6000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000 }},
    {{ 0.7000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.3000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000 }},
    {{ 0.4000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.6000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000 }},
    {{ 0.7000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.3000, 0.0000, 0.0000, 0.0000, 0.0000 }},
    {{ 0.4000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.6000, 0.0000, 0.0000, 0.0000, 0.0000 }},
    {{ 0.7000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.3000, 0.0000, 0.0000, 0.0000 }},
    {{ 0.4000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.6000, 0.0000, 0.0000, 0.0000 }},
    {{ 0.7000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.3000, 0.0000, 0.0000 }},
    {{ 0.4000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.6000, 0.0000, 0.0000 }},
    {{ 0.7000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.3000, 0.0000 }},
    {{ 0.4000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.6000, 0.0000 }},
    {{ 0.7000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.3000 }},
    {{ 0.4000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.6000 }},
    {{ 0.2100, 0.0200, 0.1500, 0.0000, 0.0000, 0.1000, 0.0100, 0.1000, 0.1500, 0.0800, 0.0300, 0.0000, 0.1500 }},
};

/** Olivine (Fo-Fa): p order fo,fa. */
struct ss_pc gh_ol_pc_xeos[9] = {
    {{ 0.9000, 0.1000 }},
    {{ 0.8000, 0.2000 }},
    {{ 0.7000, 0.3000 }},
    {{ 0.6000, 0.4000 }},
    {{ 0.5000, 0.5000 }},
    {{ 0.4000, 0.6000 }},
    {{ 0.3000, 0.7000 }},
    {{ 0.2000, 0.8000 }},
    {{ 0.1000, 0.9000 }},
};

/** Biotite (Ann-Phl): p order ann,phl. */
struct ss_pc gh_bi_pc_xeos[9] = {
    {{ 0.9000, 0.1000 }},
    {{ 0.8000, 0.2000 }},
    {{ 0.7000, 0.3000 }},
    {{ 0.6000, 0.4000 }},
    {{ 0.5000, 0.5000 }},
    {{ 0.4000, 0.6000 }},
    {{ 0.3000, 0.7000 }},
    {{ 0.2000, 0.8000 }},
    {{ 0.1000, 0.9000 }},
};

/** Feldspar (Ab-An-San): p order ab,an,san - corners, edge midpoints, interior. */
struct ss_pc gh_fsp_pc_xeos[15] = {
    {{ 0.9000, 0.0500, 0.0500 }},
    {{ 0.0500, 0.9000, 0.0500 }},
    {{ 0.0500, 0.0500, 0.9000 }},
    {{ 0.5000, 0.5000, 0.0000 }},
    {{ 0.5000, 0.0000, 0.5000 }},
    {{ 0.0000, 0.5000, 0.5000 }},
    {{ 0.7000, 0.2000, 0.1000 }},
    {{ 0.7000, 0.1000, 0.2000 }},
    {{ 0.2000, 0.7000, 0.1000 }},
    {{ 0.1000, 0.7000, 0.2000 }},
    {{ 0.2000, 0.1000, 0.7000 }},
    {{ 0.1000, 0.2000, 0.7000 }},
    {{ 0.3333, 0.3333, 0.3334 }},
    {{ 0.6000, 0.3000, 0.1000 }},
    {{ 0.6000, 0.1000, 0.3000 }},
};

/** Garnet (Gr-Py-Alm): p order gr,py,alm - corners, edge midpoints, interior. */
struct ss_pc gh_g_pc_xeos[15] = {
    {{ 0.9000, 0.0500, 0.0500 }},
    {{ 0.0500, 0.9000, 0.0500 }},
    {{ 0.0500, 0.0500, 0.9000 }},
    {{ 0.5000, 0.5000, 0.0000 }},
    {{ 0.5000, 0.0000, 0.5000 }},
    {{ 0.0000, 0.5000, 0.5000 }},
    {{ 0.7000, 0.2000, 0.1000 }},
    {{ 0.7000, 0.1000, 0.2000 }},
    {{ 0.2000, 0.7000, 0.1000 }},
    {{ 0.1000, 0.7000, 0.2000 }},
    {{ 0.2000, 0.1000, 0.7000 }},
    {{ 0.1000, 0.2000, 0.7000 }},
    {{ 0.3333, 0.3333, 0.3334 }},
    {{ 0.6000, 0.3000, 0.1000 }},
    {{ 0.6000, 0.1000, 0.3000 }},
};

/**
    Small value used to shift PC-grid corners/edges off exact 0.0/1.0
    (matches gv.bnd_val, the same eps used for gh's bounds_ref[eps,1-eps]),
    so that the unconditional Sconfig += p[i]*log(p[i]+d_em[i]) form in
    gh_objective_functions.c never sees log(0) for a non-boiled-out
    endmember whose grid coordinate happens to be an exact corner value.
*/
#define GH_PC_EPS 1e-7

/**
    Push any exact 0.0 entry in a PC-grid row to GH_PC_EPS, then renormalize
    the whole row back to Sigma=1 (rows are hand-written to sum to exactly
    1.0, so this only nudges values by ~n*GH_PC_EPS - negligible for this
    deliberately coarse levelling grid).
*/
static void GH_shift_row(double *row, int n){
    double sum = 0.0;
    for (int i = 0; i < n; i++){
        if (row[i] == 0.0){ row[i] = GH_PC_EPS; }
        sum += row[i];
    }
    for (int i = 0; i < n; i++){ row[i] /= sum; }
}

/** Apply GH_shift_row to every row of every gh PC grid, exactly once per process. */
static int gh_pc_shifted = 0;
static void GH_shift_pc_grids_once(void){
    if (gh_pc_shifted){ return; }
    for (int k = 0; k < 26; k++){ GH_shift_row(gh_liq_pc_xeos[k].xeos_pc, 13); }
    for (int k = 0; k < 9;  k++){ GH_shift_row(gh_ol_pc_xeos[k].xeos_pc,   2); }
    for (int k = 0; k < 9;  k++){ GH_shift_row(gh_bi_pc_xeos[k].xeos_pc,   2); }
    for (int k = 0; k < 15; k++){ GH_shift_row(gh_fsp_pc_xeos[k].xeos_pc,  3); }
    for (int k = 0; k < 15; k++){ GH_shift_row(gh_g_pc_xeos[k].xeos_pc,    3); }
    gh_pc_shifted = 1;
}

/**
    Per-bulk working copies of the (boundary-shifted) static grids: rebuilt
    every call from the current phase's z_em, zeroing the column of any
    boiled-out endmember and renormalizing the row back to Sigma=1, so a PC
    point that concentrated weight on a since-absent endmember doesn't get
    scored with an under-summed composition during levelling.
*/
static struct ss_pc gh_liq_pc_xeos_work[26];
static struct ss_pc gh_ol_pc_xeos_work[9];
static struct ss_pc gh_bi_pc_xeos_work[9];
static struct ss_pc gh_fsp_pc_xeos_work[15];
static struct ss_pc gh_g_pc_xeos_work[15];

static void GH_build_pc_work(struct ss_pc *src, struct ss_pc *work, int n_row, int n_col, double *z_em){
    for (int k = 0; k < n_row; k++){
        double sum = 0.0;
        for (int i = 0; i < n_col; i++){
            double v = (z_em[i] == 0.0) ? 0.0 : src[k].xeos_pc[i];
            work[k].xeos_pc[i] = v;
            sum += v;
        }
        if (sum > 0.0){
            for (int i = 0; i < n_col; i++){ work[k].xeos_pc[i] /= sum; }
        }
    }
}

void GH_pc_init_function(  PC_ref  *SS_pc_xeos,
                            int      iss,
                            char    *name,
                            double  *z_em           ){
    GH_shift_pc_grids_once();

    if (strcmp(name, "liq") == 0){
        GH_build_pc_work(gh_liq_pc_xeos, gh_liq_pc_xeos_work, 26, 13, z_em);
        SS_pc_xeos[iss].ss_pc_xeos = gh_liq_pc_xeos_work;
    }
    else if (strcmp(name, "ol") == 0){
        GH_build_pc_work(gh_ol_pc_xeos, gh_ol_pc_xeos_work, 9, 2, z_em);
        SS_pc_xeos[iss].ss_pc_xeos = gh_ol_pc_xeos_work;
    }
    else if (strcmp(name, "fsp") == 0){
        GH_build_pc_work(gh_fsp_pc_xeos, gh_fsp_pc_xeos_work, 15, 3, z_em);
        SS_pc_xeos[iss].ss_pc_xeos = gh_fsp_pc_xeos_work;
    }
    else if (strcmp(name, "bi") == 0){
        GH_build_pc_work(gh_bi_pc_xeos, gh_bi_pc_xeos_work, 9, 2, z_em);
        SS_pc_xeos[iss].ss_pc_xeos = gh_bi_pc_xeos_work;
    }
    else if (strcmp(name, "g") == 0){
        GH_build_pc_work(gh_g_pc_xeos, gh_g_pc_xeos_work, 15, 3, z_em);
        SS_pc_xeos[iss].ss_pc_xeos = gh_g_pc_xeos_work;
    }
}
