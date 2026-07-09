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

#include "GH_PP_endmembers.h"

/**
    Comp[] axis order (identical to GH_endmembers.h / GH_init_database.c):
    0 SiO2, 1 Al2O3, 2 CaO, 3 MgO, 4 FeO, 5 K2O, 6 Na2O, 7 TiO2, 8 O,
    9 MnO, 10 Cr2O3, 11 H2O, 12 CO2

    O2 and H2O pure phases are handled separately in GH_gem_function.c
    (a dedicated ideal-gas formula for O2, ported from xMELTS'
    sources/gibbs.c; the same Pitzer & Sterner (1994) EOS already used for
    the liquid's H2O basis species) and so are not part of this table.
**/
static const PP_db_gh arr_pp_db_gh[GH_N_PP] = {
    /* quartz - SiO2 */
    { "q", { 1.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -910700.0, 41.460, 2.269,
        { 80.01, -2.403E2, -35.467E5, 49.157E7, 848.0, 0.0, -9.187E-2, 24.607E-5 },
        { -2.434E-6, 10.137E-12, 23.895E-6, 0.0 } },
    /* cristobalite - SiO2 */
    { "crst", { 1.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -907753.0, 43.394, 2.587,
        { 83.51, -3.747E2, -24.554E5, 28.007E7, 535.0, 0.0, -14.216E-2, 44.142E-5 },
        { -2.515E-6, 0.0, 20.824E-6, 0.0 } },
    /* tridymite - SiO2 */
    { "trd", { 1.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -907750.0, 43.770, 2.675,
        { 75.37, 0.0, -59.581E5, 95.825E7, 383.0, 130.0, 42.670E-2, -144.575E-5 },
        { -2.508E-6, 0.0, 19.339E-6, 0.0 } },
    /* corundum - Al2O3 */
    { "cor", { 0,1.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -1675700.0, 50.820, 2.558,
        { 155.02, -8.284E2, -38.614E5, 40.908E7, 0.0, 0.0, 0.0, 0.0 },
        { -0.385E-6, 0.375E-12, 21.342E-6, 47.180E-10 } },
    /* sillimanite - Al2SiO5 = Al2O3 + SiO2 */
    { "sill", { 1.0,1.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -2586091.0, 95.930, 4.983,
        { 256.73, -18.872E2, -29.774E5, 25.096E7, 0.0, 0.0, 0.0, 0.0 },
        { -0.753E-6, 0.0, 13.431E-6, 0.0 } },
    /* andalusite - Al2SiO5 = Al2O3 + SiO2. xMELTS has no data for this
       (not modeled by MELTS itself); H/S/V/Cp ported from Theriak-Domino's
       JUN92d.bs Berman (1988) database (../theriak-domino/src/JUN92d.bs),
       cross-validated exactly against xMELTS' own sillimanite entry above
       (same H/S/V/Cp to 5+ sig figs - confirms same underlying Berman
       calibration). Theriak's volume-EOS line uses a different internal
       parameterization than xMELTS' (v1,v2,v3,v4) polynomial that could
       not be safely reverse-mapped (checked against corundum: the
       per-term scale factors were inconsistent) - reuses sillimanite's
       own real EOS terms as a physically-reasonable stand-in, since all
       three Al2SiO5 polymorphs have similar compressibility and this only
       affects a small pressure correction to G. */
    { "and", { 1.0,1.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -2589972.17, 91.4337, 5.147,
        { 236.47818, -1102.941, -7526810.0, 936442368.0, 0.0, 0.0, 0.0, 0.0 },
        { -0.753E-6, 0.0, 13.431E-6, 0.0 } },
    /* kyanite - Al2SiO5 = Al2O3 + SiO2. Same source/caveat as andalusite
       above (Theriak-Domino JUN92d.bs, Berman 1988; EOS borrowed from
       sillimanite). */
    { "ky", { 1.0,1.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -2594220.46, 82.4300, 4.412,
        { 262.68478, -2001.407, -1999740.0, -63181880.0, 0.0, 0.0, 0.0, 0.0 },
        { -0.753E-6, 0.0, 13.431E-6, 0.0 } },
    /* rutile - TiO2 */
    { "ru", { 0,0,0,0,0,0,0,1.0,0,0,0,0,0,0,0,0 },
        -944750.0, 50.460, 1.882,
        { 77.84, 0.0, -33.678E5, 40.294E7, 0.0, 0.0, 0.0, 0.0 },
        { -0.454E-6, 0.584E-12, 25.716E-6, 15.409E-10 } },
    /* sphene (titanite) - CaTiSiO5 = CaO + TiO2 + SiO2 */
    { "sph", { 1.0,0,1.0,0,0,0,0,1.0,0,0,0,0,0,0,0,0 },
        -2596652.0, 129.290, 5.565,
        { 234.62, -10.403E2, -51.183E5, 59.146E7, 0.0, 0.0, 0.0, 0.0 },
        { -0.590E-6, 0.0, 25.200E-6, 0.0 } },

    /* --- olivine endmembers (Fo-Fa, "sb-trivial" tier) --- */
    /* forsterite - Mg2SiO4 */
    { "fo", { 1.0,0,0,2.0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -2174420.0, 94.010, 4.366,
        { 238.64, -20.013E2, 0.0, -11.624E7, 0.0, 0.0, 0.0, 0.0 },
        { -0.791E-6, 1.351E-12, 29.464E-6, 88.633E-10 } },
    /* fayalite - Fe2SiO4 */
    { "fa", { 1.0,0,0,0,2.0,0,0,0,0,0,0,0,0,0,0,0 },
        -1479360.0, 150.930, 4.630,
        { 248.93, -19.239E2, 0.0, -13.910E7, 0.0, 0.0, 0.0, 0.0 },
        { -0.730E-6, 0.0, 26.546E-6, 79.482E-10 } },

    /* --- feldspar endmembers ("sb-trivial" tier) --- */
    /* albite - NaAlSi3O8 = 0.5 Na2O + 0.5 Al2O3 + 3 SiO2 (monalbite reference state) */
    { "ab", { 3.0,0.5,0,0,0,0,0.5,0,0,0,0,0,0,0,0,0 },
        -3921618.0, 224.412, 10.083,
        { 393.64, -24.155E2, -78.928E5, 107.064E7, 0.0, 0.0, 0.0, 0.0 },
        { -1.945E-6, 4.861E-12, 26.307E-6, 32.407E-10 } },
    /* anorthite - CaAl2Si2O8 = CaO + Al2O3 + 2 SiO2 */
    { "an", { 2.0,1.0,1.0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -4228730.0+3.7*4184.0, 200.186+3.7*4184.0/2200.0, 10.075,
        { 439.37, -37.341E2, 0.0, -31.702E7, 0.0, 0.0, 0.0, 0.0 },
        { -1.272E-6, 3.176E-12, 10.918E-6, 41.985E-10 } },
    /* sanidine - KAlSi3O8 = 0.5 K2O + 0.5 Al2O3 + 3 SiO2 (fully-disordered reference state) */
    { "san", { 3.0,0.5,0,0,0,0.5,0,0,0,0,0,0,0,0,0,0 },
        -3970791.0, 214.145, 10.869,
        { 381.37, -19.410E2, -120.373E5, 183.643E7, 0.0, 0.0, 0.0, 0.0 },
        { -1.805E-6, 5.112E-12, 15.145E-6, 54.850E-10 } },

    /* --- biotite endmembers ("sb-trivial" tier) --- */
    /* annite - KFe3AlSi3O10(OH)2 = 0.5 K2O + 0.5 Al2O3 + 3 FeO + 3 SiO2 + H2O */
    { "ann", { 3.0,0.5,0,0,3.0,0.5,0,0,0,0,0,1.0,0,0,0,0 },
        -5142800.0, 420.0, 15.408,
        { 727.208, -47.75040E2, -138.319E5, 211.906E7, 0.0, 0.0, 0.0, 0.0 },
        { -1.6969784E-6, 0.0, 34.4473262E-6, 0.0 } },
    /* phlogopite - KMg3AlSi3O10(OH)2 = 0.5 K2O + 0.5 Al2O3 + 3 MgO + 3 SiO2 + H2O */
    { "phl", { 3.0,0.5,0,3.0,0,0.5,0,0,0,0,0,1.0,0,0,0,0 },
        -6210391.0, 334.346, 14.977,
        { 610.37988, -20.83781E2, -215.33008E5, 284.1040896E7, 0.0, 0.0, 0.0, 0.0 },
        { -1.6969784E-6, 0.0, 34.4473262E-6, 0.0 } },

    /* --- garnet endmembers ("sb-trivial" tier) --- */
    /* almandine - Fe3Al2Si3O12 = 3 FeO + Al2O3 + 3 SiO2 */
    { "alm", { 3.0,1.0,0,0,3.0,0,0,0,0,0,0,0,0,0,0,0 },
        -5267216.0, 340.007, 11.511,
        { 573.96, -14.831E2, -292.920E5, 502.208E7, 0.0, 0.0, 0.0, 0.0 },
        { -0.558E-6, 0.321E-12, 18.613E-6, 74.539E-10 } },
    /* grossular - Ca3Al2Si3O12 = 3 CaO + Al2O3 + 3 SiO2 */
    { "gr", { 3.0,1.0,3.0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -6632859.0, 255.150, 12.538,
        { 573.43, -20.394E2, -188.872E5, 231.931E7, 0.0, 0.0, 0.0, 0.0 },
        { -0.654E-6, 1.635E-12, 18.994E-6, 79.756E-10 } },
    /* pyrope - Mg3Al2Si3O12 = 3 MgO + Al2O3 + 3 SiO2 */
    { "py", { 3.0,1.0,0,3.0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -6286548.0, 266.359, 11.316,
        { 640.72, -45.421E2, -47.019E5, 0.0E7, 0.0, 0.0, 0.0, 0.0 },
        { -0.576E-6, 0.442E-12, 22.519E-6, 37.044E-10 } },

    /* --- additional common MELTS pure phases --- */
    /* perovskite - CaTiO3 = CaO + TiO2 */
    { "perov", { 0,0,1.0,0,0,0,0,1.0,0,0,0,0,0,0,0,0 },
        -1660630.0, 93.64, 3.3626,
        { 150.49, -6.213E2, 0.0, -43.010E7, 1530.0, 550.0*4.184, 0.0, 0.0 },
        { 0.0, 0.0, 0.0, 0.0 } },
    /* calcite - CaCO3 = CaO + CO2 */
    { "cc", { 0,0,1.0,0,0,0,0,0,0,0,0,0,1.0,0,0,0 },
        -1206819.0, 91.725, 3.690,
        { 178.19, -16.577E2, -4.827E5, 16.660E7, 0.0, 0.0, 0.0, 0.0 },
        { -1.4E-6, 0.0, 8.907E-6, 227.402E-10 } },
    /* aragonite - CaCO3 = CaO + CO2 (polymorph of calcite) */
    { "arag", { 0,0,1.0,0,0,0,0,0,0,0,0,0,1.0,0,0,0 },
        -1206819.0+1100.0, 88.0, 3.415,
        { 166.62, -14.994E2, 0.0, 5.449E7, 0.0, 0.0, 0.0, 0.0 },
        { -1.4E-6, 0.0, 8.907E-6, 227.402E-10 } },
    /* magnesite - MgCO3 = MgO + CO2 */
    { "mgs", { 0,0,0,1.0,0,0,0,0,0,0,0,0,1.0,0,0,0 },
        -1113636.0, 65.210, 2.803,
        { 162.30, -11.093E2, -48.826E5, 87.466E7, 0.0, 0.0, 0.0, 0.0 },
        { -0.890E-6, 2.212E-12, 18.436E-6, 415.968E-10 } },
    /* siderite - FeCO3 = FeO + CO2 */
    { "sid", { 0,0,0,0,1.0,0,0,0,0,0,0,0,1.0,0,0,0 },
        -755900.0, 95.5, 2.938,
        { 177.36, -16.694E2, -3.551E5, 15.078E7, 0.0, 0.0, 0.0, 0.0 },
        { -0.890E-6, 2.212E-12, 18.436E-6, 415.968E-10 } },
    /* dolomite - CaMg(CO3)2 = CaO + MgO + 2 CO2 */
    { "dol", { 0,0,1.0,1.0,0,0,0,0,0,0,0,0,2.0,0,0,0 },
        -2324500.0+1100.0, 155.2, 6.434,
        { 368.02, -37.508E2, 0.0, 18.079E7, 0.0, 0.0, 0.0, 0.0 },
        { -1.4E-6, 0.0, 8.907E-6, 227.402E-10 } },
    /* spurrite - Ca5Si2O8CO3 = 2 SiO2 + 5 CaO + CO2 */
    { "spu", { 2.0,0,5.0,0,0,0,0,0,0,0,0,0,1.0,0,0,0 },
        -5840200.0, 331.0, 14.712,
        { 597.163, -36.929E2, -50.5712E5, 43.382E7, 0.0, 0.0, 0.0, 0.0 },
        { -0.890E-6, 2.212E-12, 18.436E-6, 415.968E-10 } },
    /* muscovite - KAl2(AlSi3O10)(OH)2 = 3 SiO2 + 1.5 Al2O3 + 0.5 K2O + H2O */
    { "mu", { 3.0,1.5,0,0,0,0.5,0,0,0,0,0,1.0,0,0,0,0 },
        -5976740.0, 293.157, 14.087,
        { 651.49, -38.732E2, -185.232E5, 274.247E7, 0.0, 0.0, 0.0, 0.0 },
        { -1.717E-6, 4.295E-12, 33.527E-6, 0.0 } },
};

/**
    Beta (high-T) polymorph reference states for the 3 SiO2 phases whose
    alpha-form data is in arr_pp_db_gh[] above (Berman 1988), used by
    GH_SiO2_polymorph_G() (GH_gem_function.c) when T exceeds the
    pressure-shifted transition temperature Tt+dTdP*(P-1).
**/
static const PP_db_gh_beta arr_pp_db_gh_beta[3] = {
    /* beta-quartz */
    { "q",    -908627.0, 44.207, 2.370, { -1.238E-6,  10.137E-12, 0.0,     0.0 }, 0.0237 },
    /* beta-cristobalite */
    { "crst", -906377.0, 46.029, 2.730, { -1.100E-6,   5.535E-12, 3.189E-6, 0.0 }, 0.0480 },
    /* beta-tridymite (no pressure shift on Tt: dTdP=0) */
    { "trd",  -907045.0, 45.524, 2.737, { -0.740E-6,   3.735E-12, 4.829E-6, 0.0 }, 0.0 },
};

PP_db_gh Access_GH_PP_DB(int id){
    return arr_pp_db_gh[id];
}

int GH_find_PP_id(char *name){
    for (int i = 0; i < GH_N_PP; i++){
        if (strcmp(arr_pp_db_gh[i].Name, name) == 0){
            return i;
        }
    }
    return -1;
}

PP_db_gh_beta Access_GH_SiO2_beta_DB(int id){
    return arr_pp_db_gh_beta[id];
}

int GH_find_SiO2_beta_id(char *name){
    for (int i = 0; i < 3; i++){
        if (strcmp(arr_pp_db_gh_beta[i].Name, name) == 0){
            return i;
        }
    }
    return -1;
}
