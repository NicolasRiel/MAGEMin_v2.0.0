/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : MAGEMin
 **   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 **   Developers   : Nicolas Riel, Boris Kaus
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : nriel[at]uni-mainz.de, kaus[at]uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
/**
    Objective function for the Stage-A Ghiorso/MELTS liquid ("liq"): a
    15-endmember symmetric regular solution (ideal mole-fraction mixing +
    Margules excess energy), following the exact same formulation as
    SB_database's solution models (e.g. obj_sb11_ol): n_xeos == n_em,
    p[i] = x[i] directly (no reduced/rotated basis), with the Sigma(p)=1
    closure enforced by an NLopt equality constraint (see
    GH_NLopt_opt_function.c) rather than algebraic elimination. This
    matches "gh"'s LP-only solving path, which reuses the same shared
    LP/PGE machinery as "sb".

    The regular-solution excess-energy loop (mu_Gex[i] via the eye[][]
    identity-matrix trick) is the identical pattern used by every
    SB_database solution model, generalized here from 2 endmembers
    (forsterite-fayalite) to 15.

    Ideal-mixing entropy uses TC_database's own d_em mechanism (e.g.
    TC_database/objective_functions.c: "R*T*creal(clog(sf[4]+d_em[6]))"):
    d_em[i] is 0 for a normal endmember and 1 for one boiled out because
    the bulk-rock is missing an oxide it needs (set in
    gh_gss_function.c's reference-state functions, alongside pinning that
    xeos to bounds=[0,0] and z_em[i]=0 - see GH_boiled_out there), so
    log(p[i]+d_em[i])=log(1)=0 there regardless of p[i]. Unlike an earlier
    version of this file, the log() below is unconditional (no guard on
    p[i]+d_em[i]>0), matching tc's own convention: this is safe because
    the pseudocompound grids that seed p[i] during levelling (SS_xeos_PC_gh.c)
    are shifted off exact 0.0/1.0 and, per bulk, renormalized to Sigma=1
    after any boiled-out endmember's column is zeroed (see
    GH_pc_init_function there), so a non-boiled-out endmember (d_em=0)
    never reaches p[i]=0 through that path, while NLopt's own box bounds
    (bounds_ref=[eps,1-eps]) keep it that way during refinement too. Every
    other obj_gh_* function below follows this same pattern.
*/
#include <math.h>
#include <string.h>

#include "../MAGEMin.h"
#include "gh_objective_functions.h"

double obj_gh_liq(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    int n_em    = d->n_em;
    double R    = d->R;
    double T    = d->T;
    double *p   = d->p;
    double *gb  = d->gb_lvl;
    double *mu_Gex = d->mu_Gex;

    for (int i = 0; i < n_em; i++){
        p[i] = x[i];
    }

    /* symmetric regular-solution excess energy (identical pattern to obj_sb11_ol,
       generalized from 2 to n_em endmembers)                                     */
    for (int i = 0; i < n_em; i++){
        double Gex = 0.0;
        int it = 0;
        for (int j = 0; j < n_em; j++){
            double tmp = d->eye[i][j] - p[j];
            for (int k = j+1; k < n_em; k++){
                Gex -= tmp*(d->eye[i][k]-p[k])*(d->W[it]);
                it += 1;
            }
        }
        mu_Gex[i] = Gex/1000.0;   /* J -> kJ */
        d->sf[i]  = p[i];         /* mixing fraction, used only for the generic sf>=0 check */
    }

    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    /* d_em[i]=1 for a boiled-out endmember makes log(p[i]+d_em[i])=log(1)=0
       regardless of p[i]; a non-boiled-out endmember never reaches p[i]=0
       (see this file's header comment) so the log below is unconditional. */
    double Sconfig = 0.0;
    double dSi[13];
    for (int i = 0; i < n_em; i++){
        Sconfig += p[i]*log(p[i] + d->d_em[i]);
        dSi[i]   = R*T*(log(p[i] + d->d_em[i])+1.0);
    }
    Sconfig *= R*T;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += (mu_Gex[i] + gb[i])*p[i];
    }
    d->df_raw += Sconfig;
    d->df = d->df_raw * d->factor;

    if (grad){
        for (int i = 0; i < n_em; i++){
            grad[i] = (mu_Gex[i] + gb[i] + dSi[i])*d->factor
                      - (d->df_raw*d->factor*(d->ape[i]/d->sum_apep));
        }
    }

    return d->df;
}

/**
    Olivine (Fo-Fa): a genuine reduction of xMELTS' full 6-component,
    order-parameter-bearing olivine model (sources/olivine.c: NR=5, NS=4,
    Taylor-expansion in r[]/s[] "Darken equation" coefficients), not an
    ad hoc simplification. Substituting r0=r2=r3=-1 (no Mn/Co), r4=0 (no
    Ca) and evaluating the reduced 2-variable (r1=2*xFa-1, s1=M1/M2 Fe-Mg
    ordering) Taylor polynomial gives:
        H(r1,s1) = C0 + C2*r1^2 + C5*s1^2   (C1=C3=C4=0 identically)
    with C0 = -C2 (verified numerically: H=0 at both pure endmembers) and
    C5 >= 0 (0 at Pr=1bar, growing slowly with P). Since both H and the
    ideal site-mixing entropy S(s1) are even functions of s1 minimized/
    maximized at s1=0 for every composition, the equilibrium order
    parameter is s1*=0 identically across the whole Fo-Fa join (real,
    not assumed: physically consistent with olivine showing negligible
    M1/M2 Fe-Mg ordering, unlike e.g. orthopyroxene) - so no NLopt order-
    parameter dimension is needed here, only the disordered 2-site
    entropy already used, PLUS the corrected effective regular-solution
    energy W(P) = 4*C0(P) = 20300 + 0.015*(Pbar-1) J (NOT the raw
    WH1MGFE/WH2MGFE=5075 J ordering-reaction parameter previously used
    by mistake - that constant is the M1/M2 EXCHANGE energy entering the
    Taylor coefficients, not the bulk disordered-limit regular-solution
    W; verified by direct numerical evaluation of the full xMELTS
    coefficient set reduced to the Fo-Fa limit).
*/
double obj_gh_ol(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    int n_em    = d->n_em;
    double R    = d->R;
    double T    = d->T;
    double *p   = d->p;
    double *gb  = d->gb_lvl;
    double *mu_Gex = d->mu_Gex;

    /* W(P) already computed in G_SS_gh_ol_function (gh_gss_function.c) */
    double W = d->W[0];

    for (int i = 0; i < n_em; i++){
        p[i] = x[i];
    }

    /* mu_i = W*p_j^2 (j != i), the same eye[][]-trick chemical potential
       used everywhere else, now with the corrected W */
    mu_Gex[0] = W*p[1]*p[1]/1000.0;
    mu_Gex[1] = W*p[0]*p[0]/1000.0;
    d->sf[0] = p[0]; d->sf[1] = p[1];

    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    /* unconditional log() - see obj_gh_liq's header comment */
    double Sconfig = 2.0*R*T*(p[0]*log(p[0] + d->d_em[0]) + p[1]*log(p[1] + d->d_em[1]));

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += (mu_Gex[i] + gb[i])*p[i];
    }
    d->df_raw += Sconfig;
    d->df = d->df_raw * d->factor;

    if (grad){
        double dS0 = 2.0*R*T*(log(p[0] + d->d_em[0])+1.0);
        double dS1 = 2.0*R*T*(log(p[1] + d->d_em[1])+1.0);
        grad[0] = (dS0 + mu_Gex[0] + gb[0])*d->factor - (d->df_raw*d->factor*(d->ape[0]/d->sum_apep));
        grad[1] = (dS1 + mu_Gex[1] + gb[1])*d->factor - (d->df_raw*d->factor*(d->ape[1]/d->sum_apep));
    }
    return d->df;
}

/**
    Biotite (Ann-Phl): identical pattern, single symmetric W term, no
    site-multiplicity factor on the ideal entropy (matches xMELTS'
    sources/biotite.c calibration as-is - see G_SS_gh_bi_init_function).
*/
double obj_gh_bi(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    int n_em    = d->n_em;
    double R    = d->R;
    double T    = d->T;
    double *p   = d->p;
    double *gb  = d->gb_lvl;
    double *mu_Gex = d->mu_Gex;

    for (int i = 0; i < n_em; i++){
        p[i] = x[i];
    }

    for (int i = 0; i < n_em; i++){
        double Gex = 0.0;
        int it = 0;
        for (int j = 0; j < n_em; j++){
            double tmp = d->eye[i][j] - p[j];
            for (int k = j+1; k < n_em; k++){
                Gex -= tmp*(d->eye[i][k]-p[k])*(d->W[it]);
                it += 1;
            }
        }
        mu_Gex[i] = Gex/1000.0;
        d->sf[i]  = p[i];
    }

    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    /* unconditional log() - see obj_gh_liq's header comment */
    double Sconfig = R*T*(p[0]*log(p[0] + d->d_em[0]) + p[1]*log(p[1] + d->d_em[1]));

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += (mu_Gex[i] + gb[i])*p[i];
    }
    d->df_raw += Sconfig;
    d->df = d->df_raw * d->factor;

    if (grad){
        double dS0 = R*T*(log(p[0] + d->d_em[0])+1.0);
        double dS1 = R*T*(log(p[1] + d->d_em[1])+1.0);
        grad[0] = (dS0 + mu_Gex[0] + gb[0])*d->factor - (d->df_raw*d->factor*(d->ape[0]/d->sum_apep));
        grad[1] = (dS1 + mu_Gex[1] + gb[1])*d->factor - (d->df_raw*d->factor*(d->ape[1]/d->sum_apep));
    }
    return d->df;
}

/**
    Feldspar (Ab-An-San): Elkins & Grove (1990) asymmetric ternary
    Margules, ported from xMELTS' sources/feldspar.c H/S/V macros
    (real published W parameters). Unlike the symmetric phases above,
    Wij != Wji here, so the simple eye[][] trick (which only supports one
    W per pair) cannot represent this - instead the raw partial
    derivatives of H, S_excess and V (each expanded as an explicit
    polynomial in the 3 independent p[i], no variable eliminated) are
    combined into the standard partial-molar transform
    mu_i = Gex + sum_k (delta_ik - p_k) * dGex/dp_k, which is the general
    n-component form of the same identity the eye[][] trick specializes
    to for simple symmetric regular solutions (verified: for Gex=W*x0*x1
    it reduces to exactly W*x1^2, matching obj_gh_ol/obj_sb11_ol). This
    keeps the same "direct p=x + Sigma(p)=1 equality constraint" NLopt
    formulation as every other gh/sb phase - only the excess-energy
    formula itself is richer. xMELTS' own Al-Si/albite-monalbite
    order-disorder correction (sources/gibbs.c) is not reproduced, same
    simplification as these endmembers' own standard states.
*/
double obj_gh_fsp(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    double R = d->R;
    double T = d->T;
    double Pbar = d->P*1000.0;   /* kbar -> bar, matching W calibration units */
    double *p   = d->p;
    double *gb  = d->gb_lvl;
    double *mu_Gex = d->mu_Gex;

    for (int i = 0; i < 3; i++){ p[i] = x[i]; }
    double xab = p[0], xan = p[1], xor_ = p[2];

    const double whabor=18810.0, wsabor=10.3, wvabor=0.4602;
    const double whorab=27320.0, wsorab=10.3, wvorab=0.3264;
    const double whaban=7924.0,  whanab=0.0;
    const double whanor=38974.0, wvanor=-0.1037;
    const double whoran=40317.0;
    const double whabanor=12545.0, wvabanor=-1.095;

    double KH = 0.5*(whaban+whanab+whabor+whorab+whanor+whoran) + whabanor;
    double KS = 0.5*(wsabor+wsorab);
    double KV = 0.5*(wvabor+wvorab+wvanor) + wvabanor;

    double H  = whaban*xab*xan*xan + whanab*xab*xab*xan + whabor*xab*xor_*xor_
              + whorab*xab*xab*xor_ + whanor*xan*xor_*xor_ + whoran*xan*xan*xor_
              + xab*xan*xor_*KH;
    double Sx = wsabor*xab*xor_*xor_ + wsorab*xab*xab*xor_ + xab*xan*xor_*KS;
    double V  = wvabor*xab*xor_*xor_ + wvorab*xab*xab*xor_ + wvanor*xan*xor_*xor_
              + xab*xan*xor_*KV;

    double dH[3], dS[3], dV[3];
    dH[0] = whaban*xan*xan + 2.0*whanab*xab*xan + whabor*xor_*xor_ + 2.0*whorab*xab*xor_ + xan*xor_*KH;
    dH[1] = 2.0*whaban*xab*xan + whanab*xab*xab + whanor*xor_*xor_ + 2.0*whoran*xan*xor_ + xab*xor_*KH;
    dH[2] = 2.0*whabor*xab*xor_ + whorab*xab*xab + 2.0*whanor*xan*xor_ + whoran*xan*xan + xab*xan*KH;

    dS[0] = wsabor*xor_*xor_ + 2.0*wsorab*xab*xor_ + xan*xor_*KS;
    dS[1] = xab*xor_*KS;
    dS[2] = 2.0*wsabor*xab*xor_ + wsorab*xab*xab + xab*xan*KS;

    dV[0] = 2.0*wvorab*xab*xor_ + wvabor*xor_*xor_ + xan*xor_*KV;
    dV[1] = wvanor*xor_*xor_ + xab*xor_*KV;
    dV[2] = 2.0*wvabor*xab*xor_ + wvorab*xab*xab + 2.0*wvanor*xan*xor_ + xab*xan*KV;

    double Gex = H - T*Sx + (Pbar-1.0)*V;
    double dGex[3];
    for (int k = 0; k < 3; k++){ dGex[k] = dH[k] - T*dS[k] + (Pbar-1.0)*dV[k]; }

    for (int i = 0; i < 3; i++){
        double mu = Gex;
        for (int k = 0; k < 3; k++){
            double delta_ik = (i==k) ? 1.0 : 0.0;
            mu += (delta_ik - p[k])*dGex[k];
        }
        mu_Gex[i] = mu/1000.0;   /* J -> kJ */
        d->sf[i]  = p[i];
    }

    d->sum_apep = 0.0;
    for (int i = 0; i < 3; i++){ d->sum_apep += d->ape[i]*p[i]; }
    d->factor = d->fbc/d->sum_apep;

    /* unconditional log() - see obj_gh_liq's header comment */
    double Sconfig = 0.0;
    double dSi[3];
    for (int i = 0; i < 3; i++){
        Sconfig += p[i]*log(p[i] + d->d_em[i]);
        dSi[i]   = R*T*(log(p[i] + d->d_em[i])+1.0);
    }
    Sconfig *= R*T;

    d->df_raw = 0.0;
    for (int i = 0; i < 3; i++){ d->df_raw += (mu_Gex[i] + gb[i])*p[i]; }
    d->df_raw += Sconfig;
    d->df = d->df_raw * d->factor;

    if (grad){
        for (int i = 0; i < 3; i++){
            grad[i] = (mu_Gex[i] + gb[i] + dSi[i])*d->factor - (d->df_raw*d->factor*(d->ape[i]/d->sum_apep));
        }
    }
    return d->df;
}

/**
    Garnet (Gr-Py-Alm): Berman & Koziol (1991) asymmetric ternary
    Margules, ported from xMELTS' sources/garnet.c H/S/V macros (real
    published W parameters, "1==Grossular, 2==Pyrope, 3==Almandine").
    Same mu-transform approach as obj_gh_fsp above, since Wij != Wji here
    too. Unlike feldspar, garnet's own convention uses a plain P (not
    P-1) in its volume term (see sources/garnet.c's G macro) - kept as-is.
*/
double obj_gh_g(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    double R = d->R;
    double T = d->T;
    double Pbar = d->P*1000.0;   /* kbar -> bar, matching W calibration units */
    double *p   = d->p;
    double *gb  = d->gb_lvl;
    double *mu_Gex = d->mu_Gex;

    for (int i = 0; i < 3; i++){ p[i] = x[i]; }
    double gr = p[0], py = p[1], al = p[2];

    const double WH112=21560.0, WS112=18.79, WV112=0.10;
    const double WH122=69200.0, WS122=18.79, WV122=0.10;
    const double WH113=20320.0, WS113=5.08,  WV113=0.17;
    const double WH133=2620.0,  WS133=5.08,  WV133=0.09;
    const double WH223=230.0,   WS223=0.0,   WV223=0.01;
    const double WH233=3720.0,  WS233=0.0,   WV233=0.06;
    const double WH123=0.0,     WS123=0.0,   WV123=0.0;

    double KH = 0.5*(WH112+WH122+WH113+WH133+WH223+WH233) + WH123;
    double KS = 0.5*(WS112+WS122+WS113+WS133+WS223+WS233) + WS123;
    double KV = 0.5*(WV112+WV122+WV113+WV133+WV223+WV233) + WV123;

    double H  = WH112*gr*gr*py + WH122*gr*py*py + WH113*gr*gr*al + WH133*gr*al*al
              + WH223*py*py*al + WH233*py*al*al + gr*py*al*KH;
    double Sx = WS112*gr*gr*py + WS122*gr*py*py + WS113*gr*gr*al + WS133*gr*al*al
              + WS223*py*py*al + WS233*py*al*al + gr*py*al*KS;
    double V  = WV112*gr*gr*py + WV122*gr*py*py + WV113*gr*gr*al + WV133*gr*al*al
              + WV223*py*py*al + WV233*py*al*al + gr*py*al*KV;

    double dH[3], dS[3], dV[3];
    dH[0] = 2.0*WH112*gr*py + WH122*py*py + 2.0*WH113*gr*al + WH133*al*al + py*al*KH;
    dH[1] = WH112*gr*gr + 2.0*WH122*gr*py + 2.0*WH223*py*al + WH233*al*al + gr*al*KH;
    dH[2] = WH113*gr*gr + 2.0*WH133*gr*al + WH223*py*py + 2.0*WH233*py*al + gr*py*KH;

    dS[0] = 2.0*WS112*gr*py + WS122*py*py + 2.0*WS113*gr*al + WS133*al*al + py*al*KS;
    dS[1] = WS112*gr*gr + 2.0*WS122*gr*py + 2.0*WS223*py*al + WS233*al*al + gr*al*KS;
    dS[2] = WS113*gr*gr + 2.0*WS133*gr*al + WS223*py*py + 2.0*WS233*py*al + gr*py*KS;

    dV[0] = 2.0*WV112*gr*py + WV122*py*py + 2.0*WV113*gr*al + WV133*al*al + py*al*KV;
    dV[1] = WV112*gr*gr + 2.0*WV122*gr*py + 2.0*WV223*py*al + WV233*al*al + gr*al*KV;
    dV[2] = WV113*gr*gr + 2.0*WV133*gr*al + WV223*py*py + 2.0*WV233*py*al + gr*py*KV;

    double Gex = H - T*Sx + Pbar*V;
    double dGex[3];
    for (int k = 0; k < 3; k++){ dGex[k] = dH[k] - T*dS[k] + Pbar*dV[k]; }

    for (int i = 0; i < 3; i++){
        double mu = Gex;
        for (int k = 0; k < 3; k++){
            double delta_ik = (i==k) ? 1.0 : 0.0;
            mu += (delta_ik - p[k])*dGex[k];
        }
        mu_Gex[i] = mu/1000.0;   /* J -> kJ */
        d->sf[i]  = p[i];
    }

    d->sum_apep = 0.0;
    for (int i = 0; i < 3; i++){ d->sum_apep += d->ape[i]*p[i]; }
    d->factor = d->fbc/d->sum_apep;

    /* unconditional log() - see obj_gh_liq's header comment */
    double Sconfig = 0.0;
    double dSi[3];
    for (int i = 0; i < 3; i++){
        Sconfig += p[i]*log(p[i] + d->d_em[i]);
        dSi[i]   = 3.0*R*T*(log(p[i] + d->d_em[i])+1.0);
    }
    Sconfig *= 3.0*R*T;

    d->df_raw = 0.0;
    for (int i = 0; i < 3; i++){ d->df_raw += (mu_Gex[i] + gb[i])*p[i]; }
    d->df_raw += Sconfig;
    d->df = d->df_raw * d->factor;

    if (grad){
        for (int i = 0; i < 3; i++){
            grad[i] = (mu_Gex[i] + gb[i] + dSi[i])*d->factor - (d->df_raw*d->factor*(d->ape[i]/d->sum_apep));
        }
    }
    return d->df;
}

void GH_SS_objective_init_function(    obj_type            *SS_objective,
                                        global_variable      gv                  ){
    for (int iss = 0; iss < gv.len_ss; iss++){
        if (strcmp( gv.SS_list[iss], "liq") == 0 ){
            SS_objective[iss] = obj_gh_liq;
        }
        else if (strcmp( gv.SS_list[iss], "ol") == 0 ){
            SS_objective[iss] = obj_gh_ol;
        }
        else if (strcmp( gv.SS_list[iss], "fsp") == 0 ){
            SS_objective[iss] = obj_gh_fsp;
        }
        else if (strcmp( gv.SS_list[iss], "bi") == 0 ){
            SS_objective[iss] = obj_gh_bi;
        }
        else if (strcmp( gv.SS_list[iss], "g") == 0 ){
            SS_objective[iss] = obj_gh_g;
        }
    }
}

void GH_PC_init(                       PC_type             *PC_read,
                                        global_variable      gv                  ){
    for (int iss = 0; iss < gv.len_ss; iss++){
        if (strcmp( gv.SS_list[iss], "liq") == 0 ){
            PC_read[iss] = obj_gh_liq;
        }
        else if (strcmp( gv.SS_list[iss], "ol") == 0 ){
            PC_read[iss] = obj_gh_ol;
        }
        else if (strcmp( gv.SS_list[iss], "fsp") == 0 ){
            PC_read[iss] = obj_gh_fsp;
        }
        else if (strcmp( gv.SS_list[iss], "bi") == 0 ){
            PC_read[iss] = obj_gh_bi;
        }
        else if (strcmp( gv.SS_list[iss], "g") == 0 ){
            PC_read[iss] = obj_gh_g;
        }
    }
}

/**
    p-to-xeos map for the "pure endmember as LP-swap candidate" mechanism
    (swap_pure_endmembers, simplex_levelling.c): unlike tc's per-phase
    algebra (e.g. p2x_ig_fsp drops a reference endmember to build a reduced
    basis), every gh phase is direct p=x with no basis reduction, so the
    map is just identity + clamp to bounds. One shared function suffices
    for all 5 phases (liq/ol/fsp/bi/g).
*/
void p2x_gh_generic(void *SS_ref_db, double eps){
    SS_ref *d = (SS_ref *) SS_ref_db;
    for (int i = 0; i < d->n_xeos; i++){
        d->iguess[i] = d->p[i];
        if (d->iguess[i] < d->bounds[i][0]){ d->iguess[i] = d->bounds[i][0]; }
        if (d->iguess[i] > d->bounds[i][1]){ d->iguess[i] = d->bounds[i][1]; }
    }
}

void GH_P2X_init(                      P2X_type            *P2X_read,
                                        global_variable      gv                  ){
    for (int iss = 0; iss < gv.len_ss; iss++){
        if (strcmp( gv.SS_list[iss], "liq") == 0 ||
            strcmp( gv.SS_list[iss], "ol")  == 0 ||
            strcmp( gv.SS_list[iss], "fsp") == 0 ||
            strcmp( gv.SS_list[iss], "bi")  == 0 ||
            strcmp( gv.SS_list[iss], "g")   == 0){
            P2X_read[iss] = p2x_gh_generic;
        }
    }
}

/** Evaluate a candidate pseudocompound's G and reconstruct its oxide composition
    (mirrors SB_PC_function's structure). */
SS_ref GH_PC_function(     global_variable      gv,
                            PC_type             *PC_read,
                            SS_ref               SS_ref_db,
                            bulk_info            z_b,
                            int                  ph_id                  ){

    double G0 = (*PC_read[ph_id])( SS_ref_db.n_xeos,
                                    SS_ref_db.iguess,
                                    SS_ref_db.dfx,
                                    &SS_ref_db                     );
    SS_ref_db.df = G0;

    for (int j = 0; j < gv.len_ox; j++){
        SS_ref_db.ss_comp[j] = 0.0;
    }
    for (int i = 0; i < SS_ref_db.n_em; i++){
        for (int j = 0; j < gv.len_ox; j++){
            SS_ref_db.ss_comp[j] += SS_ref_db.Comp[i][j]*SS_ref_db.p[i]*SS_ref_db.z_em[i];
        }
    }

    SS_ref_db.sf_ok = 1;

    return SS_ref_db;
}
