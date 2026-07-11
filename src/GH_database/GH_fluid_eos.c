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
