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
/**
    Closed-form standard-state Gibbs energy for Ghiorso/MELTS liquid basis
    species (research group "gh"). Unlike the THERMOCALC ("tc") and
    Stixrude ("sb") endmember EOS, which need an iterative volume/EOS
    solve, the MELTS liquid volume model used here (a linear-in-(T,P)
    "Kress" polynomial, see GH_endmembers.h) integrates in closed form, so
    no Newton iteration is needed.

    G(T,P) is built the way MELTS itself builds it: from a *solid*
    reference state at Tr=298.15K, Pr=1bar (H, S, 4-term Berman Cp),
    integrated up to the species' fusion temperature Tfus, plus a fusion
    correction (deltaH_fus = Tfus*Sfus), plus the liquid's own constant Cp
    integrated from Tfus to T, plus a pressure correction from the Kress
    liquid-volume polynomial integrated from Pr to P.

    All internal arithmetic is done in MELTS' native units (J, bar, K,
    matching xMELTS' R=8.3143 J/K and Tr=298.15K/Pr=1bar reference state);
    the result is converted to MAGEMin's convention (kJ, kbar) only in the
    final assignment to PP_ref_db.gbase, mirroring how SB_G_EM_function
    converts its own bar/J-based EOS result via "/kbar2bar" at the very end.
*/

#include <math.h>
#include <string.h>

#include "../MAGEMin.h"
#include "GH_endmembers.h"
#include "GH_PP_endmembers.h"
#include "GH_fluid_eos.h"
#include "GH_gem_function.h"

#define GH_Tr    298.15   /* reference temperature (K), matches xMELTS TR      */
#define GH_Pr    1.0      /* reference pressure (bar), matches xMELTS PR      */
#define GH_kbar2bar 1000.0

/**
    Integrate the Berman (1988) solid Cp polynomial (k0,k1,k2,k3) from Tr to
    T, adding the order-disorder lambda-transition contribution (Tt, deltaH,
    l1, l2) when the phase has one (Tt != 0) - ported from xMELTS'
    sources/gibbs.c (the generic CP_BERMAN branch of intCpsolid, shared by
    every Berman-calibrated solid that isn't given its own dedicated
    pressure-shifted transition curve, e.g. nepheline/kalsilite there and
    quartz/tridymite/cristobalite here, see GH_PP_endmembers.h).
*/
static void GH_berman_HS(  double T, double H0, double S0, const double cp[8],
                            double *H_out, double *S_out               ){
    double k0 = cp[0], k1 = cp[1], k2 = cp[2], k3 = cp[3];
    double Tt = cp[4], deltaH = cp[5], l1 = cp[6], l2 = cp[7];

    double H = H0 + k0*(T-GH_Tr) + 2.0*k1*(sqrt(T)-sqrt(GH_Tr))
             - k2*(1.0/T - 1.0/GH_Tr) - 0.5*k3*(1.0/(T*T) - 1.0/(GH_Tr*GH_Tr));
    double S = S0 + k0*log(T/GH_Tr)
             - 2.0*k1*(1.0/sqrt(T) - 1.0/sqrt(GH_Tr))
             - 0.5*k2*(1.0/(T*T) - 1.0/(GH_Tr*GH_Tr))
             - (1.0/3.0)*k3*(1.0/(T*T*T) - 1.0/(GH_Tr*GH_Tr*GH_Tr));

    if (Tt != 0.0){
        if (T > Tt){
            H += deltaH + 0.5*l1*l1*(Tt*Tt-GH_Tr*GH_Tr)
               + (2.0/3.0)*l1*l2*(Tt*Tt*Tt-GH_Tr*GH_Tr*GH_Tr)
               + 0.25*l2*l2*(Tt*Tt*Tt*Tt-GH_Tr*GH_Tr*GH_Tr*GH_Tr);
            S += deltaH/Tt + l1*l1*(Tt-GH_Tr) + l1*l2*(Tt*Tt-GH_Tr*GH_Tr)
               + (1.0/3.0)*l2*l2*(Tt*Tt*Tt-GH_Tr*GH_Tr*GH_Tr);
        }
        else{
            H += 0.5*l1*l1*(T*T-GH_Tr*GH_Tr)
               + (2.0/3.0)*l1*l2*(T*T*T-GH_Tr*GH_Tr*GH_Tr)
               + 0.25*l2*l2*(T*T*T*T-GH_Tr*GH_Tr*GH_Tr*GH_Tr);
            S += l1*l1*(T-GH_Tr) + l1*l2*(T*T-GH_Tr*GH_Tr)
               + (1.0/3.0)*l2*l2*(T*T*T-GH_Tr*GH_Tr*GH_Tr);
        }
    }
    *H_out = H;
    *S_out = S;
}

/**
    Berman (1988) solid volume-EOS pressure integral (dG = V dP from Pr to
    P), ported from xMELTS' sources/gibbs.c (intEOSsolid, EOS_BERMAN case):
    V(T,P) = V0*(1 + v1*(P-Pr) + v2*(P-Pr)^2 + v3*(T-Tr) + v4*(T-Tr)^2).
    Closed form, no iterative solve needed (unlike sb's Birch-Murnaghan-
    Debye or tc's Modified Tait volume models).
*/
static double GH_berman_EOS_dG(double T, double P, double V0, const double eos[4]){
    double v1 = eos[0], v2 = eos[1], v3 = eos[2], v4 = eos[3];
    return V0*( (v1/2.0-v2)*(P*P-GH_Pr*GH_Pr) + v2*(P*P*P-GH_Pr*GH_Pr*GH_Pr)/3.0
              + (1.0-v1+v2+v3*(T-GH_Tr)+v4*(T-GH_Tr)*(T-GH_Tr))*(P-GH_Pr) );
}

/**
    Ideal-gas standard state for O2, ported from xMELTS' sources/gibbs.c
    (the unconditional "ideal gas treatment" fallback in gibbs(), used
    whenever none of that file's SPECIAL_O2/CONSTANT_P_O2/ZERO_O2 build
    macros are set - the vanilla case). H0=0, S0=205.15 J/K are O2's own
    elemental reference values (O2 gas *is* the reference state for
    oxygen), not the Tr,Pr=298.15K/1bar solid convention used elsewhere.
*/
static double GH_O2_G(double T, double P){
    double R = 8.3143;
    double H = 23.10248*(T-GH_Tr) + 2.0*804.8876*(sqrt(T)-sqrt(GH_Tr))
             - 1762835.0*(1.0/T - 1.0/GH_Tr) - 18172.9196*log(T/GH_Tr)
             + 0.5*0.002676*(T*T-GH_Tr*GH_Tr);
    double S = 205.15 + 23.10248*log(T/GH_Tr)
             - 2.0*804.8876*(1.0/sqrt(T) - 1.0/sqrt(GH_Tr))
             - 0.5*1762835.0*(1.0/(T*T) - 1.0/(GH_Tr*GH_Tr))
             + 18172.9196*(1.0/T - 1.0/GH_Tr) + 0.002676*(T-GH_Tr)
             - R*log(P/GH_Pr);
    return H - T*S;
}

/**
    Real MELTS treatment for the SiO2 polymorphs quartz, tridymite and
    cristobalite (ported from xMELTS' sources/gibbs.c special cases for
    each name, VALUE only - no lambdaV volume-derivative terms, since
    those never feed back into gs/hs in the source either, and MAGEMin
    gets volume/moduli from finite differences on G, not analytics).
    Below the pressure-shifted transition temperature cp_t = Tt + dTdP*
    (P-1), integrates the alpha-phase Cp + lambda transition (delt =
    Tt-cp_t shifted, reducing exactly to the plain lambda formula used
    elsewhere in gh when delt=0, e.g. tridymite at any P since its
    dTdP=0); above cp_t, switches to the beta phase's own H0/S0 (Berman
    1988) reintegrated with the SAME Cp coefficients from Tr, and its own
    V0/eos_berman for the pressure term.
*/
static double GH_SiO2_polymorph_G(int beta_id, double T, double P,
                                   const PP_db_gh *alpha, const PP_db_gh_beta *beta){
    double k0 = alpha->cp_berman[0], k1 = alpha->cp_berman[1];
    double k2 = alpha->cp_berman[2], k3 = alpha->cp_berman[3];
    double Tt = alpha->cp_berman[4], l1 = alpha->cp_berman[6], l2 = alpha->cp_berman[7];
    double cp_t = Tt + beta->dTdP*(P-1.0);

    double H, S, V0; const double *eos;
    if (T > cp_t){
        H = beta->H + k0*(T-GH_Tr) + 2.0*k1*(sqrt(T)-sqrt(GH_Tr))
          - k2*(1.0/T-1.0/GH_Tr) - 0.5*k3*(1.0/(T*T)-1.0/(GH_Tr*GH_Tr));
        S = beta->S + k0*log(T/GH_Tr)
          - 2.0*k1*(1.0/sqrt(T)-1.0/sqrt(GH_Tr))
          - 0.5*k2*(1.0/(T*T)-1.0/(GH_Tr*GH_Tr))
          - (1.0/3.0)*k3*(1.0/(T*T*T)-1.0/(GH_Tr*GH_Tr*GH_Tr));
        V0 = beta->V; eos = beta->eos_berman;
    } else {
        H = alpha->H + k0*(T-GH_Tr) + 2.0*k1*(sqrt(T)-sqrt(GH_Tr))
          - k2*(1.0/T-1.0/GH_Tr) - 0.5*k3*(1.0/(T*T)-1.0/(GH_Tr*GH_Tr));
        S = alpha->S + k0*log(T/GH_Tr)
          - 2.0*k1*(1.0/sqrt(T)-1.0/sqrt(GH_Tr))
          - 0.5*k2*(1.0/(T*T)-1.0/(GH_Tr*GH_Tr))
          - (1.0/3.0)*k3*(1.0/(T*T*T)-1.0/(GH_Tr*GH_Tr*GH_Tr));

        double delt = Tt - cp_t;
        double x1 = l1*l1*delt + 2.0*l1*l2*delt*delt + l2*l2*delt*delt*delt;
        double x2 = l1*l1 + 4.0*l1*l2*delt + 3.0*l2*l2*delt*delt;
        double x3 = 2.0*l1*l2 + 3.0*l2*l2*delt;
        double x4 = l2*l2;
        double Tr_shift = GH_Tr - delt;
        H += x1*(T-Tr_shift) + x2*(T*T-Tr_shift*Tr_shift)/2.0
           + x3*(T*T*T-Tr_shift*Tr_shift*Tr_shift)/3.0
           + x4*(T*T*T*T-Tr_shift*Tr_shift*Tr_shift*Tr_shift)/4.0;
        S += x1*(log(T)-log(Tr_shift)) + x2*(T-Tr_shift)
           + x3*(T*T-Tr_shift*Tr_shift)/2.0 + x4*(T*T*T-Tr_shift*Tr_shift*Tr_shift)/3.0;

        V0 = alpha->V; eos = alpha->eos_berman;
    }
    return (H - T*S) + GH_berman_EOS_dG(T, P, V0, eos);
}

/**
    Sanidine (K-feldspar) Al-Si disordering correction, ported from
    xMELTS' sources/gibbs.c "sanidine" special case (value only - the
    extra Cp/pressure-cross terms in the source only refine hs/ss/cps for
    T<1436K and never feed back into gs beyond what's captured here via
    the td=min(1436,T) cap). Zero below Tr (T<=Tr).
*/
static double GH_sanidine_ordering_G(double T, double P){
    double td = (T < 1436.0) ? T : 1436.0;
    double d0=282.98, d1=-4.83e3, d2=36.21e5, d3=-15.733e-2, d4=34.770e-6, d5=41.063e4;

    if (T <= GH_Tr) return 0.0;

    double dhdis = d0*(td-GH_Tr) + 2.0*d1*(sqrt(td)-sqrt(GH_Tr)) - d2*(1.0/td-1.0/GH_Tr)
                 + d3*(td*td-GH_Tr*GH_Tr)/2.0 + d4*(td*td*td-GH_Tr*GH_Tr*GH_Tr)/3.0;
    double dsdis = d0*(log(td)-log(GH_Tr)) - 2.0*d1*(1.0/sqrt(td)-1.0/sqrt(GH_Tr))
                 - d2*(1.0/(td*td)-1.0/(GH_Tr*GH_Tr))/2.0 + d3*(td-GH_Tr)
                 + d4*(td*td-GH_Tr*GH_Tr)/2.0;
    double dvdis = dhdis/d5;

    return dhdis - T*dsdis + dvdis*(P-1.0);
}

/**
    Albite Al-Si tricritical ordering correction (Salje 1985), ported
    from xMELTS' sources/albite.c: a 2-order-parameter (q, qod) Landau
    model, Newton-solved here directly (2x2 analytic Hessian, same
    iteration albite.c's own order() performs) since this is a single
    pure-phase standard-state evaluation, not a joint solution-phase
    optimization - no NLopt/order-parameter machinery needed. Returns 0
    when the disordered solution (q=qod=0) is the equilibrium state
    (matches albite.c's "skip" test), i.e. above ~tc=1251K.
*/
static double GH_albite_ordering_G(double T, double P){
    const double a0=5.479, b=6854.0, aod0=41.620, bod=-9301.0, cod=43600.0;
    const double d0=-2.171, d1=-3.043, d2=-0.001569, d3=0.000002109;
    const double tc=1251.0, tod=824.1;
    const double vscale = 335282.925;

    double s0 = 0.6, s1 = 0.9;
    for (int iter = 0; iter < 200; iter++){
        double pf = 1.0 + (P-1.0)/vscale;
        double dgds0 = (-a0*tc*s0 + b*s0*s0*s0 + (d0-d2*T*T-2.0*d3*T*T*T)*s1)*pf
                     + T*(a0*s0 + (d1+2.0*d2*T+3.0*d3*T*T)*s1);
        double dgds1 = (-aod0*tod*s1 + bod*s1*s1*s1 + cod*s1*s1*s1*s1*s1
                        + (d0-d2*T*T-2.0*d3*T*T*T)*s0)*pf
                     + T*(aod0*s1 + (d1+2.0*d2*T+3.0*d3*T*T)*s0);

        double j00 = (-a0*tc + 3.0*b*s0*s0)*pf + T*a0;
        double j01 = (d0-d2*T*T-2.0*d3*T*T*T)*pf + T*(d1+2.0*d2*T+3.0*d3*T*T);
        double j11 = (-aod0*tod + 3.0*bod*s1*s1 + 5.0*cod*s1*s1*s1*s1)*pf + T*aod0;

        double det = j00*j11 - j01*j01;
        if (fabs(det) < 1e-30) break;
        double inv00 =  j11/det, inv01 = -j01/det, inv11 =  j00/det;

        double ds0 = inv00*dgds0 + inv01*dgds1;
        double ds1 = inv01*dgds0 + inv11*dgds1;
        s0 -= ds0; s1 -= ds1;

        if (fabs(ds0) < 1e-10 && fabs(ds1) < 1e-10) break;
    }

    if (fabs(s0) < 1.4901161193847656e-8 && fabs(s1) < 1.4901161193847656e-8){
        return 0.0;   /* disordered (high albite/monalbite): no correction */
    }

    double pf = 1.0 + (P-1.0)/vscale;
    double S = -(0.5*a0*s0*s0 + 0.5*aod0*s1*s1 + (d1+2.0*d2*T+3.0*d3*T*T)*s0*s1);
    double H = -0.5*a0*tc*s0*s0 + 0.25*b*s0*s0*s0*s0
             - 0.5*aod0*tod*s1*s1 + 0.25*bod*s1*s1*s1*s1
             + (cod/6.0)*s1*s1*s1*s1*s1*s1
             + (d0-d2*T*T-2.0*d3*T*T*T)*s0*s1;
    double V = H/vscale;
    return H - T*S + (P-1.0)*V;
}

static PP_ref GH_pack_PP_ref(char *name, int len_ox, const double *Comp,
                              double *bulk_rock, double *apo, double gbase_J){
    PP_ref PP_ref_db;
    strcpy(PP_ref_db.Name, name);
    for (int i = 0; i < len_ox; i++){
        PP_ref_db.Comp[i] = Comp[i];
    }
    double fbc = 0.0;
    for (int i = 0; i < len_ox; i++){
        fbc += bulk_rock[i]*apo[i];
    }
    double ape = 0.0;
    for (int i = 0; i < len_ox; i++){
        ape += Comp[i]*apo[i];
    }
    PP_ref_db.gbase             = gbase_J/GH_kbar2bar;   /* J -> kJ */
    PP_ref_db.factor             = fbc/ape;
    PP_ref_db.phase_shearModulus = 0.0;
    PP_ref_db.phase_bulkModulus  = 0.0;
    PP_ref_db.phase_expansivity  = 0.0;
    PP_ref_db.phase_cp           = 0.0;
    return PP_ref_db;
}

PP_ref GH_G_EM_function(   int          EM_database,
                            int          len_ox,
                            int         *id,
                            double      *bulk_rock,
                            double      *apo,
                            double       Pkbar,
                            double       T,
                            char        *name,
                            char        *state          ){

    double P  = Pkbar * GH_kbar2bar;   /* kbar -> bar */

    /* Common rock-forming pure phases (research group "gh"'s own thermodynamic
       basis, ported from xMELTS' solid-phase database - see
       GH_PP_endmembers.h), not part of the 13-species liquid endmember table
       looked up via find_EM_id() below: must be special-cased before that
       lookup, since they aren't in gh's endmember hashtable.                  */
    /* SiO2 polymorphs: real MELTS dual alpha/beta reference-state
       treatment (pressure-shifted transition), not the generic Berman-
       solid path below - see GH_SiO2_polymorph_G(). */
    int beta_id = GH_find_SiO2_beta_id(name);
    if (beta_id >= 0 && (strcmp(name,"q")==0 || strcmp(name,"crst")==0 || strcmp(name,"trd")==0)){
        int pp_id_sio2 = GH_find_PP_id(name);
        PP_db_gh alpha = Access_GH_PP_DB(pp_id_sio2);
        PP_db_gh_beta beta = Access_GH_SiO2_beta_DB(beta_id);
        double gbase_J = GH_SiO2_polymorph_G(beta_id, T, P, &alpha, &beta);
        return GH_pack_PP_ref(name, len_ox, alpha.Comp, bulk_rock, apo, gbase_J);
    }

    int pp_id = GH_find_PP_id(name);
    if (pp_id >= 0){
        PP_db_gh PP_return = Access_GH_PP_DB(pp_id);
        double H_T, S_T;
        GH_berman_HS(T, PP_return.H, PP_return.S, PP_return.cp_berman, &H_T, &S_T);
        double gbase_J = (H_T - T*S_T) + GH_berman_EOS_dG(T, P, PP_return.V, PP_return.eos_berman);
        /* real MELTS Al-Si order-disorder corrections for these two
           feldspar endmember standard states (sources/gibbs.c/albite.c) */
        if (strcmp(name, "ab") == 0){
            gbase_J += GH_albite_ordering_G(T, P);
        }
        else if (strcmp(name, "san") == 0){
            gbase_J += GH_sanidine_ordering_G(T, P);
        }
        return GH_pack_PP_ref(name, len_ox, PP_return.Comp, bulk_rock, apo, gbase_J);
    }
    if (strcmp(name, "O2") == 0){
        double O2_Comp[16] = {0,0,0,0,0,0,0,0,2.0,0,0,0,0,0,0,0};   /* O2 = 2 "O" (excess-oxygen axis) */
        double gbase_J = GH_O2_G(T, P);
        return GH_pack_PP_ref(name, len_ox, O2_Comp, bulk_rock, apo, gbase_J);
    }

    /* endmember lookup (implicit declaration of find_EM_id, resolved at
       link time against src/hash_init.h's definition - same pattern used
       by SB_G_EM_function/TC_G_EM_function, neither of which includes
       hash_init.h directly either)                                        */
    int p_id            = find_EM_id(name);
    EM_db_gh EM_return   = Access_GH_EM_DB(p_id);

    /* H2O and CO2 are real supercritical/gas-like fluids at magmatic
       conditions, not "hypothetical liquids" built from a crystalline
       reference + fusion the way the oxide components are - MELTS itself
       gets their standard state from a dedicated real-gas EOS (Pitzer &
       Sterner 1994), not the Berman-Cp/Kress-volume construction below.
       gbase_J from GH_pitzer_sterner_G() is already on the same
       elements-referenced (Berman 1988) scale as every other endmember. */
    if (strcmp(name, "H2O") == 0 || strcmp(name, "CO2") == 0){
        double gbase_J = GH_pitzer_sterner_G(strcmp(name,"H2O")==0, T, P);
        return GH_pack_PP_ref(name, len_ox, EM_return.Comp, bulk_rock, apo, gbase_J);
    }

    double Tfus = EM_return.Tfus;
    double Sfus = EM_return.Sfus;
    double Cpl  = EM_return.Cpl;
    double Vl    = EM_return.Vl;
    double dvdt   = EM_return.kress[0];
    double dvdp   = EM_return.kress[1];
    double d2vdtp = EM_return.kress[2];
    double d2vdp2 = EM_return.kress[3];

    /* integrate the SOLID's Berman Cp (+ lambda transition, if any) from Tr to Tfus */
    double H_sol_Tfus, S_sol_Tfus;
    GH_berman_HS(Tfus, EM_return.H, EM_return.S, EM_return.cp_berman, &H_sol_Tfus, &S_sol_Tfus);

    /* fusion correction: deltaG_fus(Tfus) = 0 => deltaH_fus = Tfus*deltaS_fus */
    double H_liq_Tfus = H_sol_Tfus + Tfus*Sfus;
    double S_liq_Tfus = S_sol_Tfus + Sfus;

    /* integrate the LIQUID's constant Cp from Tfus to T */
    double H_liq_T = H_liq_Tfus + Cpl*(T - Tfus);
    double S_liq_T = S_liq_Tfus + Cpl*log(T/Tfus);

    double G_Pr = H_liq_T - T*S_liq_T;

    /* pressure correction: integral of the Kress liquid-volume polynomial
       V(T,P) = Vl + dvdt*(T-Tr) + dvdp*(P-Pr) + d2vdtp*(T-Tr)*(P-Pr) + 0.5*d2vdp2*(P-Pr)^2
       from Pr to P, at fixed T                                              */
    double dT = T - GH_Tr;
    double dP = P - GH_Pr;
    double dG_P =  Vl*dP
                 + dvdt*dT*dP
                 + 0.5*dvdp*dP*dP
                 + 0.5*d2vdtp*dT*dP*dP
                 + (d2vdp2/6.0)*dP*dP*dP;

    double gbase_J = G_Pr + dG_P;

    PP_ref PP_ref_db = GH_pack_PP_ref(name, len_ox, EM_return.Comp, bulk_rock, apo, gbase_J);
    PP_ref_db.phase_cp = Cpl/GH_kbar2bar;

    return PP_ref_db;
}
