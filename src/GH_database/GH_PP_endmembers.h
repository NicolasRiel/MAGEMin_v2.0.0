/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : MAGEMin
 **   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 **   Developers   : Nicolas Riel, Boris Kaus
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : nriel[at]uni-mainz.de, kaus[at]uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
#ifndef __PP_ENDMEMBER_DATABASE_GH_H_
#define __PP_ENDMEMBER_DATABASE_GH_H_

    /**
        Store standard-state thermodynamic data for a common rock-forming
        pure (solid) phase, in the "gh" (Ghiorso/MELTS) research group.

        Data ported verbatim from xMELTS-master/includes/sol_struct_data.h
        (the canonical MELTS solid-phase database also used to build the
        liquid model's fusion reference states, see GH_endmembers.c): H, S,
        V at Tr=298.15K/Pr=1bar and the 8-term Berman (1988) Cp polynomial
        (k0,k1,k2,k3, plus an optional order-disorder lambda transition
        Tt,deltaH,l1,l2) and 4-term Berman volume-EOS (v1,v2,v3,v4). See
        GH_gem_function.c for the closed-form G(T,P) integration.

        quartz/tridymite/cristobalite: the plain generic lambda formula
        above is NOT what real MELTS uses for these three - it special-
        cases each with its own pressure-shifted transition temperature
        and a distinct high-T (beta) reference state (own H0/S0/V0/EOS),
        ported faithfully in GH_gem_function.c's GH_SiO2_polymorph_G()
        (see GH_SiO2_beta[] below for the beta-phase data and
        GH_gem_function.c for the dTdP pressure-shift slopes). The extra
        "lambdaV" volume-derivative cross terms real MELTS also computes
        are NOT reproduced, since MAGEMin gets pure-phase volume/moduli
        from finite differences on G (toolkit.c), not analytic derivatives
        - only the Gibbs energy value itself needs to be exact.

        Also holds the pure-endmember standard states of the mineral solid
        solutions ("sb-trivial" tier: olivine, feldspar, biotite, garnet -
        i.e. those whose xMELTS mixing model is plain ideal+regular on a
        single site, matching sb's n_xeos==n_em pattern, with no internal
        order parameter): forsterite/fayalite, albite/anorthite/sanidine,
        annite/phlogopite, almandine/grossular/pyrope. These are not part
        of gh's PP_list (they are not standalone pure phases in xMELTS,
        only solution-phase endmembers) but are looked up the same way,
        via GH_find_PP_id(), so the solid-solution objective functions can
        fetch their gbase(T,P) exactly like gh_liq's endmembers already do
        (see gh_gss_function.c's get_em_data calls). albite's monalbite
        reference state gets xMELTS' real Salje (1985) tricritical Al-Si
        ordering correction (2-parameter Landau model, Newton-solved once
        per gbase call - GH_albite_ordering_G() in GH_gem_function.c) and
        sanidine gets its real disordering correction (GH_sanidine_ordering_G());
        both ported from sources/gibbs.c/albite.c, not approximated.
    **/
    typedef struct PP_db_gh_ {
        char   Name[16];
        double Comp[16];        /** oxide composition, same axis order as GH_endmembers.h */
        double H;                /** reference enthalpy at Tr,Pr (J)                       */
        double S;                /** reference entropy at Tr,Pr (J/K)                      */
        double V;                /** reference volume at Tr,Pr (J/bar)                     */
        double cp_berman[8];     /** Berman Cp coefficients: k0,k1,k2,k3,Tt,deltah,l1,l2    */
        double eos_berman[4];    /** Berman volume EOS coefficients: v1,v2,v3,v4            */
    } PP_db_gh;

    /**
        Beta (high-T) polymorph reference states for quartz/tridymite/
        cristobalite - own H0,S0 (reintegrated with the SAME k0-k3 Cp as
        the alpha phase, see arr_pp_db_gh) and own V0/eos_berman, used by
        GH_SiO2_polymorph_G() when T exceeds the pressure-shifted
        transition temperature. dTdP is the transition curve's slope
        (K/bar): quartz 0.0237, tridymite 0.0 (no shift), cristobalite
        0.0480 (Berman 1988).
    **/
    typedef struct PP_db_gh_beta_ {
        char   Name[16];
        double H;
        double S;
        double V;
        double eos_berman[4];
        double dTdP;
    } PP_db_gh_beta;

    PP_db_gh_beta Access_GH_SiO2_beta_DB(int id);
    int           GH_find_SiO2_beta_id(char *name);

    #define GH_N_PP 27

    PP_db_gh Access_GH_PP_DB(int id);
    int      GH_find_PP_id(char *name);

#endif
