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
    Pitzer & Sterner (1994) real-gas EOS for pure H2O and pure CO2, ported
    from xMELTS-master/sources/fluid.c's fluidPhase(). That original
    function handles a general H2O-CO2 MIXTURE (composition x[2]); this
    port keeps only the two pure-component special cases (x=[1,0] and
    x=[0,1]) that "gh" actually needs for the liquid's H2O/CO2 standard
    states, dropping the cross/mixing terms (aMix, the general branch of
    idealGas(), and the higher-order (2nd/3rd) T,rho-derivatives that
    fluidPhase() only needed for its Cp/V outputs, not G).

    References (as cited in the original xMELTS source):
      H2O: Pitzer KS and Sterner SM (1994) Equations of state valid
           continuously from zero to extreme pressures for H2O and CO2.
           J Chem Phys 101: 3111-6
      CO2: Sterner SM and Pitzer KS (1994) An equation of state for
           carbon dioxide valid from zero to extreme pressures.
           Contr Mineral Petrol 117: 362-74

    The physics/equations themselves (virial-type Pitzer-Sterner
    coefficients, density solved by Newton iteration from a Redlich-Kwong
    initial guess, then G = A + P/rho with a final shift onto the
    Berman (1988) elemental reference-state scale) are a faithful,
    unabridged port - not a re-derivation.
*/
#include <math.h>
#include <float.h>

#include "GH_fluid_eos.h"

#define nCOEFF 11

/* Pitzer & Sterner (1994) coefficients. Columns per row: [unused, T^-4, T^-2, T^-1, const, T, T^2] */
static const double h2o_c[nCOEFF][7] = {
    { 0.0,            0.0,             0.0,            0.0,            0.0,            0.0,           0.0 },
    { 0.0,            0.0,             0.0,   0.24657688e+6,  0.51359951e+2,           0.0,           0.0 },
    { 0.0,            0.0,             0.0,   0.58638965e+0, -0.28646939e-2, 0.31375577e-4,           0.0 },
    { 0.0,            0.0,             0.0,  -0.62783840e+1,  0.14791599e-1, 0.35779579e-3, 0.15432925e-7 },
    { 0.0,            0.0,             0.0,            0.0, -0.42719875e+0, -0.16325155e-4,           0.0 },
    { 0.0,            0.0,             0.0,   0.56654978e+4, -0.16580167e+2, 0.76560762e-1,           0.0 },
    { 0.0,            0.0,             0.0,            0.0,  0.10917883e+0,           0.0,           0.0 },
    { 0.0,   0.38878656e+13,  -0.13494878e+9,  0.30916564e+6,  0.75591105e+1,           0.0,           0.0 },
    { 0.0,            0.0,             0.0,  -0.65537898e+5,  0.18810675e+3,           0.0,           0.0 },
    { 0.0,  -0.14182435e+14,   0.18165390e+9, -0.19769068e+6, -0.23530318e+2,           0.0,           0.0 },
    { 0.0,            0.0,             0.0,   0.92093375e+5,  0.12246777e+3,           0.0,           0.0 },
};

static const double co2_c[nCOEFF][7] = {
    { 0.0,            0.0,             0.0,            0.0,            0.0,            0.0,           0.0 },
    { 0.0,            0.0,             0.0,   0.18261340e+7,  0.79224365e+2,           0.0,           0.0 },
    { 0.0,            0.0,             0.0,            0.0,  0.66560660e-4, 0.57152798e-5, 0.30222363e-9 },
    { 0.0,            0.0,             0.0,            0.0,  0.59957845e-2, 0.71669631e-4, 0.62416103e-8 },
    { 0.0,            0.0,             0.0,  -0.13270279e+1, -0.15210731e+0, 0.53654244e-3, -0.71115142e-7 },
    { 0.0,            0.0,             0.0,   0.12456776e+0,  0.49045367e+1, 0.98220560e-2, 0.55962121e-5 },
    { 0.0,            0.0,             0.0,            0.0,  0.75522299e+0,           0.0,           0.0 },
    { 0.0,  -0.39344644e+12,   0.90918237e+8,  0.42776716e+6, -0.22347856e+2,           0.0,           0.0 },
    { 0.0,            0.0,             0.0,   0.40282608e+3,  0.11971627e+3,           0.0,           0.0 },
    { 0.0,            0.0,    0.22995650e+8, -0.78971817e+5, -0.63376456e+2,           0.0,           0.0 },
    { 0.0,            0.0,             0.0,   0.95029765e+5,  0.18038071e+2,           0.0,           0.0 },
};

/* Ideal-gas contribution coefficients (Cooper 1982, for H2O; Berman 1988 /
   Stull & Prophet, for CO2), same values fluidPhase()'s idealGas() uses. */
static const double H2O_b[9]    = { 0.0, 0.134865, -5.005143, 4.006320, 0.012436, 0.973150, 1.279500, 0.969560, 0.248730 };
static const double H2O_beta[9] = { 0.0, 0.0,       0.0,      0.0,      833.0,    2289.0,   5009.0,   5982.0,   17800.0 };

/**
    Redlich-Kwong initial density guess (Edminster 1968 cubic solution),
    ported verbatim from fluid.c's redlichKwong(). Returns the gas-branch
    root (index 0), which is what fluidPhase() always uses as its initial
    guess before Newton-refining against the full Pitzer-Sterner EOS.
*/
static double gh_redlichKwong_Z(double t, double p, double b, double a2b){
    double n, m, arg;

    if (a2b <= 0.0) a2b = 0.001;

    n = b*p*(a2b-b*p-1.0)/3.0 - a2b*b*p*b*p - 2.0/27.0;
    m = b*p*(a2b-b*p-1.0) - 1.0/3.0;

    arg = n*n/4.0 + m*m*m/27.0;
    if (arg > 0.0){
        double term1 = - n/2.0 + sqrt(arg);
        double term2 = - n/2.0 - sqrt(arg);
        double z = (term1 >= 0.0) ?  pow( term1, 1.0/3.0) : -pow(-term1, 1.0/3.0);
        z        += (term2 >= 0.0) ?  pow( term2, 1.0/3.0) : -pow(-term2, 1.0/3.0);
        z        += 1.0/3.0;
        return z;
    }
    else if (arg == 0.0){
        return 1.0;
    }
    else{
        double PI = acos(-1.0);
        double cosPhi = (n > 0) ? -sqrt(- (n*n/4.0) / (m*m*m/27.0)) : sqrt(- (n*n/4.0) / (m*m*m/27.0));
        double phi = acos(cosPhi);
        double r1  = 2.0*sqrt(-m/3.0)*cos(phi/3.0) + 1.0/3.0;
        double r2  = 2.0*sqrt(-m/3.0)*cos(phi/3.0 + 2.0*PI/3.0) + 1.0/3.0;
        double r3  = 2.0*sqrt(-m/3.0)*cos(phi/3.0 + 4.0*PI/3.0) + 1.0/3.0;
        double zgas = r1;
        if (r2 > zgas) zgas = r2;
        if (r3 > zgas) zgas = r3;
        return zgas;
    }
}

/**
    Pitzer & Sterner (1994) pure-component Gibbs energy, ported from
    fluidPhase()'s density-solve + residual-term + ideal-gas-term +
    "G = A + P/rho" logic, specialized to a single pure component
    (dropping the H2O-CO2 mixing terms entirely).
    T in K, Pbar in bar. Returns G in J/mol (Berman 1988 elemental
    reference scale, matching the rest of MAGEMin's endmembers).
*/
double GH_pitzer_sterner_G(int is_H2O, double T, double Pbar){

    const double R  = 83.14241; /* cm^3-bar/mol-K */
    const double P0 = 1.01325;  /* reference pressure (bar) */

    double a, b;
    if (is_H2O){
        a = (111.3057 + 50.70033*exp(-0.982646e-2*(T-273.15)))*10.0e5;
        b = 14.6;
    }
    else{
        a = (73.03 - 0.0714*(T-273.15) + 2.157e-5*(T-273.15)*(T-273.15))*10.0e5;
        b = 29.7;
    }

    double c[nCOEFF];
    for (int i = 1; i < nCOEFF; i++){
        const double *row = is_H2O ? h2o_c[i] : co2_c[i];
        c[i] = row[1]/(T*T*T*T) + row[2]/(T*T) + row[3]/T + row[4] + row[5]*T + row[6]*T*T;
    }

    double Zgas = gh_redlichKwong_Z(T, Pbar, b/(82.05*T), a/(b*82.05*pow(T,1.5)));
    double rh   = Pbar/(Zgas*R*T);

    double dp = DBL_MAX;
    for (int count = 1; count <= 200 && fabs(dp) > 10.0*DBL_EPSILON; count++){
        double temp1  = c[3] + 2.0*c[4]*rh + 3.0*c[5]*rh*rh + 4.0*c[6]*rh*rh*rh;
        double dtemp1 = 2.0*c[4] + 6.0*c[5]*rh + 12.0*c[6]*rh*rh;
        double temp2  = c[2] + c[3]*rh + c[4]*rh*rh + c[5]*rh*rh*rh + c[6]*rh*rh*rh*rh;
        double dtemp2 = c[3] + 2.0*c[4]*rh + 3.0*c[5]*rh*rh + 4.0*c[6]*rh*rh*rh;

        double pr  = rh + c[1]*rh*rh - rh*rh*temp1/(temp2*temp2)
                   + c[7]*rh*rh*exp(-c[8]*rh) + c[9]*rh*rh*exp(-c[10]*rh);
        pr *= R*T;
        pr -= Pbar;

        double dpr = 1.0 + 2.0*c[1]*rh - 2.0*rh*temp1/(temp2*temp2)
                   - rh*rh*(temp2*temp2*dtemp1 - temp1*2.0*temp2*dtemp2)/(temp2*temp2*temp2*temp2)
                   + 2.0*c[7]*rh*exp(-c[8]*rh) - c[7]*c[8]*rh*rh*exp(-c[8]*rh)
                   + 2.0*c[9]*rh*exp(-c[10]*rh) - c[9]*c[10]*rh*rh*exp(-c[10]*rh);
        dpr *= R*T;

        dp  = -pr/dpr;
        rh += dp;
        if (rh < 0.0) rh = 10.0*DBL_EPSILON;
    }

    double term1      = c[2] + c[3]*rh + c[4]*rh*rh + c[5]*rh*rh*rh + c[6]*rh*rh*rh*rh;
    double dterm1drh  = c[3] + 2.0*c[4]*rh + 3.0*c[5]*rh*rh + 4.0*c[6]*rh*rh*rh;

    double term2a     = c[7]/c[8];
    double dterm2bdrh = -c[8]*exp(-c[8]*rh);
    double term2      = term2a*(exp(-c[8]*rh) - 1.0);
    double dterm2drh  = term2a*dterm2bdrh;

    double term3a     = c[9]/c[10];
    double dterm3bdrh = -c[10]*exp(-c[10]*rh);
    double term3      = term3a*(exp(-c[10]*rh) - 1.0);
    double dterm3drh  = term3a*dterm3bdrh;

    double Ar     = c[1]*rh + 1.0/term1 - 1.0/c[2] - term2 - term3;
    double dArdrh = c[1] - dterm1drh/(term1*term1) - dterm2drh - dterm3drh;

    double Ai, dAidrh;
    if (is_H2O){
        double f = -1.0 + log(rh*R*T/P0) + H2O_b[1] + H2O_b[2]/T + H2O_b[3]*(1.0-log(T));
        for (int i = 4; i <= 8; i++){
            f += H2O_b[i]*log(1.0 - exp(-H2O_beta[i]/T));
        }
        Ai     = R*T*f;
        dAidrh = R*T/rh;
    }
    else{
        const double k0 = 93.0*10.0, k1 = -13.409e2*10.0, k2 = 1.238e5*10.0, k4 = -0.002876*10.0, k5 = 6336.2*10.0;
        const double Gref = -94265.0*4.184*10.0, Sref = 51.072*4.184*10.0, Tref = 298.15;

        Ai = -k2/(2.0*T) - 2.0*k1*T/sqrt(Tref) - k2*T/(2.0*Tref*Tref) - k5*T/Tref
           - k0*T*log(T) + k5*log(T) + R*T*log(R*T*rh) + 4.0*k1*sqrt(T) - k4*T*T/2.0
           + k0*T*log(Tref) - R*T*log(P0) - R*T - Sref*T + k0*T + k4*Tref*T + k2/Tref
           - k5*log(Tref) + Tref*Sref + Gref - Tref*k0 - 2.0*k1*sqrt(Tref) - k4*Tref*Tref/2.0 + k5;
        dAidrh = R*T/rh;
    }

    double A     = R*T*Ar + Ai;
    double dAdrh = R*T*dArdrh + dAidrh;
    double G     = (A + rh*dAdrh)/10.0; /* bar*cm3/mol -> J/mol */

    /* shift the EOS's own 298.15K/1bar value onto the Berman (1988)
       elemental reference scale, exactly as fluidPhase() does          */
    if (is_H2O){
        const double refH2O_g = -228538.00;
        G += refH2O_g - (-46493.8016496949);
    }
    else{
        const double refCO2_g = -394341.00;
        G += refCO2_g - (-394450.0);
    }

    return G;
}

/**
    Pure-H2O Pitzer-Sterner "h - T*s" recomputation, needed for the two
    real gibbs.c contexts that use this instead of fluidPhase()'s own raw
    *g output (pMELTS' standalone "water" phase, gibbs.c ~line 2336's
    "gH2O = hH2O - t*sH2O", commented there as "used gH2O for calibration";
    and, less obviously, pMELTS' own LIQUID "H2O" basis species too - gh's own H2O
    liquid-table row has all-zero Kress dV/dP terms, matching real
    xMELTS' own table, which means the *downstream*, unconditional
    Birch-Murnaghan block (gibbs.c ~line 1493) always takes its "no
    pressure correction" else-branch for H2O specifically, and THAT
    branch unconditionally recomputes gl = hl - t*sl from whatever hl/sl
    the H2O-specific branch left behind - so the liquid's own gl ALSO
    ends up as h-T*s, not the raw g the H2O branch itself computed. Both
    findings made 2026-07-16 by direct comparison against a real gibbs()
    dispatch (not assumed from the source alone) - see
    [[gh-spn-liq-gbase-verification]].

    Derivation (avoids porting fluidPhase()'s full dA/dT machinery):
    internally, *g=A+p/rh, *h=A+p/rh-t*dAdt, *s=-dAdt satisfy h-t*s=g
    EXACTLY before fluidPhase()'s own reference-state shift - so the
    only reason h-T*s differs from g in the FINAL (shifted) output is
    that fluid.c applies three INDEPENDENT additive shifts to g, h, s
    (its own calibrated constants, not a thermodynamically consistent
    triple: shift_g != shift_h - T*shift_s). Since GH_pitzer_sterner_G
    already returns g_raw+shift_g, this recovers h_shifted-T*s_shifted
    as g_raw + shift_h - T*shift_s = GH_pitzer_sterner_G(...) +
    (shift_h-shift_g) - T*shift_s, using fluid.c's own refH2O
    g/h/s constants (-228538.00/-241816.00/188.72) and its own
    298.15K/1bar EOS-only baseline values (-46493.8016496949/
    9430.96281231262/187.572582951064). Verified exact (<0.001 kJ) at
    two independent (T,P) against real fluidPhase() output directly.
*/
double GH_pitzer_sterner_H2O_hTs_G(double T, double Pbar){
    const double refH2O_g = -228538.00,   baseline_g = -46493.8016496949;
    const double refH2O_h = -241816.00,   baseline_h =   9430.96281231262;
    const double refH2O_s =     188.72,   baseline_s =    187.572582951064;

    double shift_g = refH2O_g - baseline_g;
    double shift_h = refH2O_h - baseline_h;
    double shift_s = refH2O_s - baseline_s;

    double g_current = GH_pitzer_sterner_G(1, T, Pbar);
    return g_current + (shift_h - shift_g) - T*shift_s;
}

/* Cooper (1982) H2O ideal-gas term, factored out of GH_pitzer_sterner_G's
   is_H2O branch above (same formula, unchanged) so GH_pitzer_sterner_mix_G
   below can evaluate it at the mixture's shared density. */
static void GH_ideal_gas_H2O(double T, double rh, double *Ai, double *dAidrh){
    const double R  = 83.14241;
    const double P0 = 1.01325;
    double f = -1.0 + log(rh*R*T/P0) + H2O_b[1] + H2O_b[2]/T + H2O_b[3]*(1.0-log(T));
    for (int i = 4; i <= 8; i++){
        f += H2O_b[i]*log(1.0 - exp(-H2O_beta[i]/T));
    }
    *Ai     = R*T*f;
    *dAidrh = R*T/rh;
}

/* Stull & Prophet CO2 ideal-gas term, factored out of GH_pitzer_sterner_G's
   !is_H2O branch above (same formula, unchanged). */
static void GH_ideal_gas_CO2(double T, double rh, double *Ai, double *dAidrh){
    const double R  = 83.14241;
    const double P0 = 1.01325;
    const double k0 = 93.0*10.0, k1 = -13.409e2*10.0, k2 = 1.238e5*10.0, k4 = -0.002876*10.0, k5 = 6336.2*10.0;
    const double Gref = -94265.0*4.184*10.0, Sref = 51.072*4.184*10.0, Tref = 298.15;

    *Ai = -k2/(2.0*T) - 2.0*k1*T/sqrt(Tref) - k2*T/(2.0*Tref*Tref) - k5*T/Tref
        - k0*T*log(T) + k5*log(T) + R*T*log(R*T*rh) + 4.0*k1*sqrt(T) - k4*T*T/2.0
        + k0*T*log(Tref) - R*T*log(P0) - R*T - Sref*T + k0*T + k4*Tref*T + k2/Tref
        - k5*log(Tref) + Tref*Sref + Gref - Tref*k0 - 2.0*k1*sqrt(Tref) - k4*Tref*Tref/2.0 + k5;
    *dAidrh = R*T/rh;
}

double GH_pitzer_sterner_mix_G(double x_h2o, double T, double Pbar, double *dGdx_h2o){

    const double R  = 83.14241; /* cm^3-bar/mol-K */
    const double x0 = x_h2o, x1 = 1.0 - x_h2o;

    /* van der Waals one-fluid Redlich-Kwong mixing, used only to seed the
       initial density guess - the only place a genuine H2O-CO2 cross
       term enters (ported verbatim from fluidPhase()).                  */
    double aCO2  = (73.03 - 0.0714*(T-273.15) + 2.157e-5*(T-273.15)*(T-273.15))*10.0e5;
    double bCO2  = 29.7;
    double aH2O  = (111.3057 + 50.70033*exp(-0.982646e-2*(T-273.15)))*10.0e5;
    double bH2O  = 14.6;
    double a0CO2 = 46.0e6;
    double a0H2O = exp(4.881243 + 0.1823047e-2*(T-273.15) - 0.1712269e-4*(T-273.15)*(T-273.15)
                       + 6.479419e-8*(T-273.15)*(T-273.15)*(T-273.15))*10.0e5;
    double k     = exp(-11.071 + 5953.0/T - 2.746e6/(T*T) + 4.646e8/(T*T*T));
    double aMix  = sqrt(a0H2O*a0CO2) + k*82.05*82.05*pow(T, 2.5);

    double aRK = x0*x0*aH2O + x1*x1*aCO2 + 2.0*x0*x1*aMix;
    double bRK = x0*bH2O + x1*bCO2;

    double Zgas = gh_redlichKwong_Z(T, Pbar, bRK/(82.05*T), aRK/(bRK*82.05*pow(T,1.5)));
    double rh   = Pbar/(Zgas*R*T);

    /* mixture virial coefficients: mole-fraction-linear combination of the
       pure H2O/CO2 coefficients ("an initial guess at the functional
       form" per fluidPhase()'s own comment - ported as-is). c[] mixes;
       dc[] = c_H2O-c_CO2 is the x0-derivative of that linear mix, needed
       below for dG/dx_h2o.                                              */
    double c[nCOEFF], dc[nCOEFF];
    for (int i = 1; i < nCOEFF; i++){
        const double *rowH = h2o_c[i];
        const double *rowC = co2_c[i];
        double cH2O = rowH[1]/(T*T*T*T) + rowH[2]/(T*T) + rowH[3]/T + rowH[4] + rowH[5]*T + rowH[6]*T*T;
        double cCO2 = rowC[1]/(T*T*T*T) + rowC[2]/(T*T) + rowC[3]/T + rowC[4] + rowC[5]*T + rowC[6]*T*T;
        c[i]  = x0*cH2O + x1*cCO2;
        dc[i] = cH2O - cCO2;
    }

    double dp = DBL_MAX;
    for (int count = 1; count <= 200 && fabs(dp) > 10.0*DBL_EPSILON; count++){
        double temp1  = c[3] + 2.0*c[4]*rh + 3.0*c[5]*rh*rh + 4.0*c[6]*rh*rh*rh;
        double dtemp1 = 2.0*c[4] + 6.0*c[5]*rh + 12.0*c[6]*rh*rh;
        double temp2  = c[2] + c[3]*rh + c[4]*rh*rh + c[5]*rh*rh*rh + c[6]*rh*rh*rh*rh;
        double dtemp2 = c[3] + 2.0*c[4]*rh + 3.0*c[5]*rh*rh + 4.0*c[6]*rh*rh*rh;

        double pr  = rh + c[1]*rh*rh - rh*rh*temp1/(temp2*temp2)
                   + c[7]*rh*rh*exp(-c[8]*rh) + c[9]*rh*rh*exp(-c[10]*rh);
        pr *= R*T;
        pr -= Pbar;

        double dpr = 1.0 + 2.0*c[1]*rh - 2.0*rh*temp1/(temp2*temp2)
                   - rh*rh*(temp2*temp2*dtemp1 - temp1*2.0*temp2*dtemp2)/(temp2*temp2*temp2*temp2)
                   + 2.0*c[7]*rh*exp(-c[8]*rh) - c[7]*c[8]*rh*rh*exp(-c[8]*rh)
                   + 2.0*c[9]*rh*exp(-c[10]*rh) - c[9]*c[10]*rh*rh*exp(-c[10]*rh);
        dpr *= R*T;

        dp  = -pr/dpr;
        rh += dp;
        if (rh < 0.0) rh = 10.0*DBL_EPSILON;
    }

    double term1      = c[2] + c[3]*rh + c[4]*rh*rh + c[5]*rh*rh*rh + c[6]*rh*rh*rh*rh;
    double dterm1drh  = c[3] + 2.0*c[4]*rh + 3.0*c[5]*rh*rh + 4.0*c[6]*rh*rh*rh;
    double term1_x    = dc[2] + dc[3]*rh + dc[4]*rh*rh + dc[5]*rh*rh*rh + dc[6]*rh*rh*rh*rh;

    double term2a     = c[7]/c[8];
    double term2a_x   = dc[7]/c[8] - c[7]*dc[8]/(c[8]*c[8]);
    double term2b     = exp(-c[8]*rh) - 1.0;
    double dterm2bdrh = -c[8]*exp(-c[8]*rh);
    double term2b_x   = -dc[8]*rh*exp(-c[8]*rh);
    double term2      = term2a*term2b;
    double dterm2drh  = term2a*dterm2bdrh;
    double term2_x    = term2a_x*term2b + term2a*term2b_x;

    double term3a     = c[9]/c[10];
    double term3a_x   = dc[9]/c[10] - c[9]*dc[10]/(c[10]*c[10]);
    double term3b     = exp(-c[10]*rh) - 1.0;
    double dterm3bdrh = -c[10]*exp(-c[10]*rh);
    double term3b_x   = -dc[10]*rh*exp(-c[10]*rh);
    double term3      = term3a*term3b;
    double dterm3drh  = term3a*dterm3bdrh;
    double term3_x    = term3a_x*term3b + term3a*term3b_x;

    double Ar     = c[1]*rh + 1.0/term1 - 1.0/c[2] - term2 - term3;
    double dArdrh = c[1] - dterm1drh/(term1*term1) - dterm2drh - dterm3drh;
    double Ar_x   = dc[1]*rh - term1_x/(term1*term1) + dc[2]/(c[2]*c[2]) - term2_x - term3_x;

    double Ai_h2o, dAidrh_h2o, Ai_co2, dAidrh_co2;
    GH_ideal_gas_H2O(T, rh, &Ai_h2o, &dAidrh_h2o);
    GH_ideal_gas_CO2(T, rh, &Ai_co2, &dAidrh_co2);
    double Ai      = x0*Ai_h2o + x1*Ai_co2;
    double dAidrh  = x0*dAidrh_h2o + x1*dAidrh_co2;
    double Ai_x    = Ai_h2o - Ai_co2;

    double A      = R*T*Ar   + Ai;
    double dAdrh  = R*T*dArdrh + dAidrh;
    double A_x    = R*T*Ar_x + Ai_x;   /* dA/dx_h2o at fixed rho (envelope theorem) */

    double G      = (A + rh*dAdrh)/10.0; /* bar*cm3/mol -> J/mol */
    double dGdx   = A_x/10.0;

    /* reference-state shift, weighted by composition (reduces exactly to
       GH_pitzer_sterner_G's own shift at x_h2o=1 or x_h2o=0)             */
    const double refH2O_g = -228538.00;
    const double refCO2_g = -394341.00;
    G    += x0*(refH2O_g - (-46493.8016496949)) + x1*(refCO2_g - (-394450.0));
    dGdx += (refH2O_g - (-46493.8016496949)) - (refCO2_g - (-394450.0));

    *dGdx_h2o = dGdx;
    return G;
}

/* ============================================================================
   Haar (1984) H2O EOS at fixed p=1 bar, ported from water.c's whaar().
   Value + 1st-rho-derivative only (drops T-derivative and 2nd-rho-
   derivative bookkeeping only needed for real whaar()'s unused Cp/V
   outputs). See GH_fluid_eos.h for verification/provenance notes.
   ============================================================================ */
static void GH_kubik(double b, double c, double d, double *x1, double *x2, double *x2i, double *x3){
    double q, p, r, phi3, ff;
    const double pi = 3.14159263538979;
    *x2 = 0.0; *x2i = 0.0; *x3 = 0.0;
    if (c == 0.0 && d == 0.0) { *x1 = -b; return; }
    q = (2.0*b*b*b/27.0 - b*c/3.0 + d)/2.0;
    p = (3.0*c - b*b)/9.0;
    ff = fabs(p);
    r = sqrt(ff);
    ff = r*q;
    if (ff < 0.0) r = -r;
    ff = q/(r*r*r);
    if (p > 0.0) {
        phi3 = log(ff + sqrt(ff*ff+1.0))/3.0;
        *x1 = -r*(exp(phi3) - exp(-phi3)) - b/3.0;
        *x2i = 1.0;
    } else if (q*q + p*p*p > 0.0) {
        phi3 = log(ff + sqrt(ff*ff-1.0))/3.0;
        *x1 = -r*(exp(phi3) + exp(-phi3)) - b/3.0;
        *x2i = 1.0;
    } else {
        phi3 = atan(sqrt(1.0-ff*ff)/ff)/3.0;
        *x1 = -2.0*r*cos(phi3) - b/3.0;
        *x2 = 2.0*r*cos(pi/3.0-phi3) - b/3.0;
        *x2i = 0.0;
        *x3 = 2.0*r*cos(pi/3.0+phi3) - b/3.0;
    }
}

static double GH_psat2(double t){
    double w, wsq, v, ff;
    int i;
    static const double a[9] = { 0.0,
        -7.8889166, 2.5514255, -6.716169, 33.239495,
        -105.38479, 174.35319, -148.39348, 48.631602
    };
    if (t <= 314.0) return exp(6.3573118-8858.843/t + 607.56335/pow(t, 0.6));
    v = t/647.25;
    w = fabs(1.0-v);
    wsq = sqrt(w);
    ff = 0.0;
    for(i=1;i<=8;i++) { ff = ff + a[i]*w; w = w*wsq; }
    return 220.93*exp(ff/v);
}

double GH_haar_H2O_G(double t, double p){
    static const double gi[41] = { 0.0,
        -.53062968529023e4,  .22744901424408e5,  .78779333020687e4,
        -.69830527374994e3,  .17863832875422e6, -.39514731563338e6,
         .33803884280753e6, -.13855050202703e6, -.25637436613260e7,
         .48212575981415e7, -.34183016969660e7,  .12223156417448e7,
         .11797433655832e8, -.21734810110373e8,  .10829952168620e8,
        -.25441998064049e7, -.31377774947767e8,  .52911910757704e8,
        -.13802577177877e8, -.25109914369001e7,  .46561826115608e8,
        -.72752773275387e8,  .41774246148294e7,  .14016358244614e8,
        -.31555231392127e8,  .47929666384584e8,  .40912664781209e7,
        -.13626369388386e8,  .69625220862664e7, -.10834900096447e8,
        -.22722827401688e7,  .38365486000660e7,  .68833257944332e5,
         .21757245522644e6, -.26627944829770e5, -.70730418082074e6,
        -.225e1,            -1.68e1,             .055e1,
        -93.0e1
    };
    static const int ki[41] = { 0,
        1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5,
        6, 6, 6, 6, 7, 7, 7, 7, 9, 9, 9, 9, 3, 3, 1, 5, 2, 2, 2, 4
    };
    static const int li[41] = { 0,
        1, 2, 4, 6, 1, 2, 4, 6, 1, 2, 4, 6, 1, 2, 4, 6, 1, 2, 4, 6,
        1, 2, 4, 6, 1, 2, 4, 6, 1, 2, 4, 6, 0, 3, 3, 3, 0, 2, 0, 0
    };
    static const double rhoi[41] = { 0.0,
        0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0, 0.319, 0.310, 0.310, 1.55
    };
    static const double ttti[41] = { 0.0,
        0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0, 640.0, 640.0, 641.6, 270.0
    };
    static const double alpi[41] = { 0.0,
        0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0, 34.0, 40.0, 30.0, 1050.0
    };
    static const double beti[41] = { 0.0,
        0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0, 2.0e4, 2.0e4, 4.0e4, 25.0
    };
    static const double bi[6]  = { 0.7478629, -0.3540782, 0.0, 0.007159876, 0.0, -0.003528426 };
    static const double bbi[6] = { 1.1278334, -0.5944001, -5.010996, 0.0, 0.63684256, 0.0 };
    static const double ci[19] = { 0.0,
        .19730271018e2,      .209662681977e2,   -.483429455355,
        .605743189245e1,   22.56023885,        -9.87532442,
       -.43135538513e1,      .458155781,        -.47754901883e-1,
        .41238460633e-2,    -.27929052852e-3,    .14481695261e-4,
       -.56473658748e-6,     .16200446e-7,      -.3303822796e-9,
        .451916067368e-11,  -.370734122708e-13,  .137546068238e-15
    };

    const double r     = 4.6152;
    const double gref  = -54955.23970014;
    const double t0    = 647.073;
    const double rr    = 8.31441;
    const double alpha = 11.0;
    const double beta  = 133.0/3.0;
    const double gammaC= 7.0/2.0;
    const double P0    = 1.01325;
    /* p is now a caller-supplied parameter (was hardcoded to 1.0) - real
       gibbs.c calls whaar() at TWO different pressures depending on which
       standard state is being built: the liquid's own H2O basis species
       always uses the fixed 1-bar reference (matches the original port),
       while the standalone "water" pure phase (GH_gem_function.c's own
       "water" branch) calls whaar() at the ACTUAL pressure (capped at
       10000 bar) - see GH_wdh78_G below for the correction above that cap. */

    double taui[7];
    int i, count;

    taui[0] = 1.0; taui[1] = t/t0; for(i=2;i<=6;i++) taui[i] = taui[i-1]*taui[1];

    double b  = bi[1]*log(taui[1]) + bi[0] + bi[3]/taui[3] + bi[5]/taui[5];
    double bb = bbi[0] + bbi[1]/taui[1] + bbi[2]/taui[2] + bbi[4]/taui[4];

    double ps = 220.55;
    if (t <= 647.25) ps = GH_psat2(t);

    double ark = 1.279186e8 - 2.241415e4 * t;
    double brk = 1.428062e1 + 6.092237e-4 * t;
    double oft = ark/(p*sqrt(t));
    double buk = -10.0*rr*t/p;
    double cuk = oft - brk*brk + brk*buk;
    double duk = - brk*oft;
    double x1, x2, x2i, x3;
    GH_kubik(buk, cuk, duk, &x1, &x2, &x2i, &x3);

    double vol;
    if (x2i != 0.0) vol = x1;
    else            vol = (p < ps) ? fmax(x1, fmax(x2, x3)) : fmin(x1, fmin(x2, x3));

    double rhn = (vol <= 0.0) ? 1.9 : (1.0/vol)*18.0152;

    double dp = DBL_MAX, dr = DBL_MAX, rh = rhn;
    for (count=1; count<=100 && (dp > 10.0*DBL_EPSILON || dr > 10.0*DBL_EPSILON); count++){
        double y, ermi[10], prv, dpr;
        rh = rhn;
        if (rh <= 0.0) rh = 1.e-8;
        if (rh >  1.9) rh = 1.9;
        y = rh*b/4.0;
        ermi[0] = 1.0; ermi[1] = 1.0-exp(-rh);
        for (i=2;i<=9;i++) ermi[i] = ermi[i-1]*ermi[1];

        prv = 0.0; dpr = 0.0;
        for (i=1; i<=36; i++) {
            prv = prv + gi[i]/taui[li[i]]*ermi[ki[i]-1];
            dpr = dpr + (2.0+rh*(ki[i]*exp(-rh)-1.0)/ermi[1])*gi[i]/taui[li[i]]*ermi[ki[i]-1];
        }
        for(i=37; i<=40; i++) {
            double del, tau, abc, q10, qm;
            del = rh/rhoi[i] - 1.0;
            tau = t/ttti[i]  - 1.0;
            abc = -alpi[i] * pow(del, (double) ki[i]) - beti[i] * tau*tau;
            q10 = (abc > -100.00) ? gi[i] * pow(del, (double) li[i]) * exp(abc) : 0.0;
            qm = li[i]/del - ki[i]*alpi[i]*pow(del, (double) (ki[i]-1));
            prv = prv + q10*qm*rh*rh/rhoi[i];
            dpr = dpr + (q10*qm*rh*rh/rhoi[i]) * (2.0/rh+qm/rhoi[i])
                      - rh*rh/(rhoi[i]*rhoi[i])*q10*
                        (li[i]/del/del + ki[i]*(ki[i]-1)*alpi[i]*pow(del, (double) (ki[i]-2)));
        }

        prv = rh*(rh*exp(-rh)*prv + r*t*((1.0 + alpha*y + beta*y*y)/pow(1.0-y,3.0)
                                + 4.0*y*(bb/b - gammaC)));
        dpr = rh*exp(-rh)*dpr
            + r*t*( (1.0 + 2.0*alpha*y + 3.0*beta*y*y)/pow(1.0-y,3.0)
                   + 3.0*y*(1.0 + alpha*y + beta*y*y)/pow(1.0-y,4.0)
                   + 2.0*4.0*y*(bb/b - gammaC));

        if (dpr <= 0.0) rhn *= (p <= ps) ? 0.95 : 1.05;
        else {
            double x;
            if (dpr < 0.01) dpr = 0.01;
            x = (p - prv)/dpr;
            if (fabs(x) > 0.1) x = 0.1*x/fabs(x);
            rhn = rh + x;
        }
        dp = fabs(1.0 - prv/p);
        dr = fabs(1.0 - rhn/rh);
    }
    rh = rhn;

    /* base function Z, value + rho-derivative only */
    double y = rh*b/4.0;
    double dydrh = b/4.0;

    double Z = - log(1.0-y) - (beta-1.0)/(1.0-y)
             + (alpha+beta+1.0)/(2.0*(1.0-y)*(1.0-y)) + 4*y*(bb/b - gammaC)
             - (alpha-beta+3.0)/2.0 + log(rh*r*t/P0);
    double dZdrh = 1.0/rh + 4.0*(bb/b-gammaC)*dydrh + dydrh/(1.0-y)
             + (alpha+beta+1.0)*dydrh/pow(1.0-y,3.0)
             - (beta-1.0)*dydrh/pow(1.0-y,2.0);

    double Ab = r*t*Z;
    double dAbdrh = r*t*dZdrh;

    /* residual function Ar, value + rho-derivative only */
    double ermi[10];
    ermi[0] = 1.0; ermi[1] = 1.0-exp(-rh);
    for (i=2; i<=9; i++) ermi[i] = ermi[i-1]*ermi[1];

    double Ar = 0.0, dArdrh = 0.0;
    for (i=1; i<=36; i++){
        Ar     += gi[i]/ki[i]/taui[li[i]]*ermi[ki[i]];
        dArdrh += gi[i]/taui[li[i]]*ermi[ki[i]-1]*exp(-rh);
    }
    for(i=37; i<=40; i++){
        double del = rh/rhoi[i] - 1.0;
        double tau = t/ttti[i] - 1.0;
        double Q = -alpi[i]*pow(del, (double) ki[i]) - beti[i]*tau*tau;
        double dQdrh = (ki[i] == 0) ? 0.0 : -alpi[i]*ki[i]*pow(del, (double) (ki[i]-1))/rhoi[i];
        double expQ = (Q > -100.0) ? exp(Q) : 0.0;
        Ar += gi[i]*pow(del, (double) li[i])*expQ;
        dArdrh += (li[i] == 0) ? gi[i]*expQ*dQdrh
                : gi[i]*li[i]*pow(del, (double) (li[i]-1))*expQ/rhoi[i]
                  + gi[i]*pow(del, (double) li[i])*expQ*dQdrh;
    }

    /* ideal gas Ai, T-only value (no rho-dependence) */
    double tr = t/1.0e2;
    double Zi = 1.0 + (ci[1]/tr + ci[2])*log(tr);
    for (i=3; i<=18; i++) Zi += ci[i]*pow(tr, (double) (i-6));
    double Ai = -r*t*Zi;

    double A     = Ab + Ar + Ai;
    double dAdrh = dAbdrh + dArdrh;

    double pcalc = rh*rh*dAdrh;
    double gH2O  = A + pcalc/rh;

    gH2O *= 1.80152; /* -> J/mol */

    /* Berman shift, MODE_xMELTS/MODE__MELTS/MODE__MELTSandCO2/MODE__MELTSandCO2_H2O branch (the only one gh needs) */
    gH2O += -285829.96 - (298.15*69.9146) - gref;

    return gH2O;
}

/* ============================================================================
   wdh78 high-pressure correction for water.c's whaar() EOS, ported from
   water.c's own wdh78() ("Returns difference in thermodynamic properties
   of water between p,t and 10kb,t"). Real gibbs.c's standalone "water"
   pure-phase branch calls whaar() at min(P,10000) then, if the ACTUAL
   pressure exceeds 10000 bar, adds this delta to extend the EOS beyond
   whaar()'s own reliable range - found missing entirely from gh (which
   only had the liquid's fixed-1-bar H2O, never the standalone "water"
   phase) during the 2026-07-15 grid sweep. Value-only port (drops the
   h/s/cp/v/dvdt/etc: outputs real wdh78() also computes, none needed for
   a standard-state G).
   ============================================================================ */
static double GH_wdh78_poly(double t_celsius, double p_bar){
    static const double a[5][5] = {
        { -5.6130073e+04,  3.8101798e-01, -2.1167697e-06,  2.0266445e-11, -8.3225572e-17 },
        { -1.5285559e+01,  1.3752390e-04, -1.5586868e-09,  6.6329577e-15,  0.0           },
        { -2.6092451e-02,  3.5988857e-08, -2.7916588e-14,  0.0,            0.0           },
        {  1.7140501e-05, -1.6860893e-11,  0.0,            0.0,            0.0           },
        { -6.0126987e-09,  0.0,            0.0,            0.0,            0.0           }
    };
    double g = 0.0;
    for (int j = 0; j < 5; j++){
        for (int l = 0; l < 5-j; l++){
            g += a[j][l]*pow(t_celsius, (double) j)*pow(p_bar, (double) l);
        }
    }
    return g;
}

double GH_wdh78_G(double t, double p){
    double t_celsius = t - 273.15;
    double g10kb = GH_wdh78_poly(t_celsius, 10000.0);
    double g_p   = GH_wdh78_poly(t_celsius, p);
    return (g_p - g10kb) * 4.184; /* cal -> J, matching real wdh78()'s own *4.184 */
}

/* ============================================================================
   Duan (1992) pure-CO2 EOS at fixed p=1 bar, ported from fluidPhase.c's
   duanCO2Driver()+idealGasCO2(), specialized to the pure endpoint and
   value-only. See GH_fluid_eos.h for verification/provenance notes.
   ============================================================================ */
static const double GH_CO2Tc = 304.1282;
static const double GH_CO2Pc = 73.773;
#define GH_CO2Vc (8.314467*GH_CO2Tc/GH_CO2Pc)

static const double GH_CO2La1  =  1.14400435E-01;
static const double GH_CO2La2  = -9.38526684E-01;
static const double GH_CO2La3  =  7.21857006E-01;
static const double GH_CO2La4  =  8.81072902E-03;
static const double GH_CO2La5  =  6.36473911E-02;
static const double GH_CO2La6  = -7.70822213E-02;
static const double GH_CO2La7  =  9.01506064E-04;
static const double GH_CO2La8  = -6.81834166E-03;
static const double GH_CO2La9  =  7.32364258E-03;
static const double GH_CO2La10 = -1.10288237E-04;
static const double GH_CO2La11 =  1.26524193E-03;
static const double GH_CO2La12 = -1.49730823E-03;
static const double GH_CO2La   =  7.81940730E-03;
static const double GH_CO2Lb   = -4.22918013E+00;
static const double GH_CO2Lc   =  1.58500000E-01;

/* Berman(1988)-style ideal-gas Cp/H/S polynomial coefficients for CO2,
   ported from fluidPhase.c's idealCoeff[][CO2] column (index 1)          */
static const double GH_idealCoeff_CO2[13] = {
    -1.8188731,       12.903022,       -9.6634864,       4.2251879,
    -1.0421640,       0.12683515,      -0.0049939675,    2.4950242,
    -0.82723750,      0.15372481,      -0.015861243,     0.00086017150,
    -0.000019222165
};

double GH_duan_CO2_G(double t){
    const double R_duan = 8.314467;
    const double p = 1.0; /* fixed reference pressure - matches propertiesOfPureCO2(t, 1.0, ...) call site in real gibbs.c */
    double CO2Tr = t/GH_CO2Tc;

    double bEnd = GH_CO2La1 + GH_CO2La2/CO2Tr/CO2Tr + GH_CO2La3/CO2Tr/CO2Tr/CO2Tr;
    double cEnd = GH_CO2La4 + GH_CO2La5/CO2Tr/CO2Tr + GH_CO2La6/CO2Tr/CO2Tr/CO2Tr;
    double dEnd = GH_CO2La7 + GH_CO2La8/CO2Tr/CO2Tr + GH_CO2La9/CO2Tr/CO2Tr/CO2Tr;
    double eEnd = GH_CO2La10 + GH_CO2La11/CO2Tr/CO2Tr + GH_CO2La12/CO2Tr/CO2Tr/CO2Tr;
    double fEnd = GH_CO2La/CO2Tr/CO2Tr/CO2Tr;

    double bv = bEnd*GH_CO2Vc;
    double cv = cEnd*GH_CO2Vc*GH_CO2Vc;
    double dv = dEnd*GH_CO2Vc*GH_CO2Vc*GH_CO2Vc*GH_CO2Vc;
    double ev = eEnd*GH_CO2Vc*GH_CO2Vc*GH_CO2Vc*GH_CO2Vc*GH_CO2Vc;
    double fv = fEnd*GH_CO2Vc*GH_CO2Vc;
    double beta = GH_CO2Lb;
    double gammav = GH_CO2Lc*GH_CO2Vc*GH_CO2Vc;

    double bvPrime = 2.0*bv;
    double cvPrime = 3.0*cv;
    double dvPrime = 5.0*dv;
    double evPrime = 6.0*ev;
    double fvPrime = 2.0*fv;
    double betaPrime = beta;
    double gammavPrime = 3.0*gammav;

    double v, z;
    int iter = 0;
    double delv = 1.0, vPrevious = 1.0, delvPrevious = 1.0;
    v = R_duan*t/p;
    while (iter < 200) {
        z = 1.0 + bv/v + cv/v/v + dv/v/v/v/v + ev/v/v/v/v/v + (fv/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v);
        delv = z*R_duan*t/p - v;
        if ( ((iter > 1) && (delv*delvPrevious < 0.0)) || (fabs(delv) < v*100.0*DBL_EPSILON) ) break;
        vPrevious = v;
        delvPrevious = delv;
        v = (z*R_duan*t/p + v)/2.0;
        iter++;
    }
    if (fabs(delv) > v*100.0*DBL_EPSILON && iter < 200) {
        double dx;
        double rtb = (delv < 0.0) ? (dx = vPrevious-v,v) : (dx = v-vPrevious,vPrevious);
        iter = 0;
        while (iter < 200) {
            v = rtb + (dx *= 0.5);
            z = 1.0 + bv/v + cv/v/v + dv/v/v/v/v + ev/v/v/v/v/v + (fv/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v);
            delv = z*R_duan*t/p - v;
            if (delv <= 0.0) rtb = v;
            if ( (fabs(dx) < 100.0*DBL_EPSILON) || (delv == 0.0) ) break;
            iter++;
        }
    }

    double lnPhiCO2 = 0.0;
    lnPhiCO2 += -log(z);
    lnPhiCO2 += bvPrime/v;
    lnPhiCO2 += cvPrime/2.0/v/v;
    lnPhiCO2 += dvPrime/4.0/v/v/v/v;
    lnPhiCO2 += evPrime/5.0/v/v/v/v/v;
    lnPhiCO2 += ((fvPrime*beta + betaPrime*fv)/2.0/gammav)*(1.0-exp(-gammav/v/v));
    lnPhiCO2 += ((fvPrime*gammav+gammavPrime*fv-fv*beta*(gammavPrime-gammav))/2.0/gammav/gammav)
               *(1.0 - (gammav/v/v + 1.0)*exp(-gammav/v/v));
    lnPhiCO2 += ((gammavPrime-gammav)*fv/2.0/gammav/gammav)*(-2.0 + (gammav*gammav/v/v/v/v + 2.0*gammav/v/v + 2.0)*exp(-gammav/v/v));

    double phi = exp(lnPhiCO2);

    /* idealGasCO2: h0(t), s0(t) integrated Cp polynomial */
    double s0=0.0, h0=0.0;
    int i;
    for (i=0; i<7; i++) h0 += GH_idealCoeff_CO2[i]*pow(t/1000.0, (double) (i+1))/((double) (i+1));
    h0 += GH_idealCoeff_CO2[7]*log(t/1000.0);
    for (i=8; i<13; i++) h0 += GH_idealCoeff_CO2[i]/pow(t/1000.0, (double) (i-7))/((double) (7-i));
    s0 = GH_idealCoeff_CO2[0]*log(t/1000.0);
    for (i=1; i<7; i++) s0 += GH_idealCoeff_CO2[i]*pow(t/1000.0, (double) i)/((double) i);
    for (i=7; i<13; i++) s0 += GH_idealCoeff_CO2[i]/pow(t/1000.0, (double) (i-6))/((double) (6-i));
    h0 *= 8.31451*1000.0;
    s0 *= 8.31451;
    h0 += -385358.2260;
    s0 +=  210.0304;

    const double R_outer = 8.3143;
    double g = h0 - t*s0 + R_outer*t*log(phi*p);
    return g;
}
