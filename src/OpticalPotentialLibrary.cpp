// OpticalPotentialLibrary.cpp — named optical-model potentials.
//
// Ported verbatim from Cleopatra's potentials.h. Each function reproduces the
// original arithmetic exactly (same operation order, same literals), so values
// are bit-identical to the original globals-based code; the only change is that
// results are returned in an OMPset instead of assigned to file-scope globals.

#include "OpticalPotentialLibrary.h"
#include <cmath>

using std::string;
using std::pow;
using std::exp;
using std::sqrt;

// ====================================================================
// Citations (verbatim from potentials.h::potentialRef)
// ====================================================================
string potentialRef(const string& name) {
    //======== Deuteron
    if (name == "A") return "An and Cai (2006) E < 183 | 12 < A < 238 | http://dx.doi.org/10.1103/PhysRevC.73.054605";
    if (name == "H") return "Han, Shi, Shen (2006) E < 200 | 12 < A < 209 | http://dx.doi.org/10.1103/PhysRevC.74.044615";
    if (name == "B") return "Bojowald et al.(1988) 50 < E < 80 | 27 < A < 208 | http://dx.doi.org/10.1103/PhysRevC.38.1153";
    if (name == "D") return "Daehnick, Childs, Vrcelj (1980) 11.8 < E < 80 | 27 < A < 238 (REL) | http://dx.doi.org/10.1103/PhysRevC.21.2253";
    if (name == "C") return "Daehnick, Childs, Vrcelj (1980) 11.8 < E < 80 | 27 < A < 238 (NON-REL) | http://dx.doi.org/10.1103/PhysRevC.21.2253";
    if (name == "L") return "Lohr and Haeberli (1974) 9 < E < 13 | 40 < A | http://dx.doi.org/10.1016/0375-9474(74)90627-7";
    if (name == "Q") return "Perey and Perey (1963) 12 < E < 25  | 40 < A | http://dx.doi.org/10.1016/0370-1573(91)90039-O";
    if (name == "Z") return "Zhang, Pang, Lou (2016) 5 < E < 170  | A < 18 | https://doi.org/10.1103/PhysRevC.94.014619";

    //======= Proton
    if (name == "K") return "Koning and Delaroche (2009) 0.001 < E < 200 | 24 < A < 209 | Iso. Dep. | http://dx.doi.org/10.1016/S0375-9474(02)01321-0";
    if (name == "V") return "Varner et al., (CH89) (1991) 16 < E < 65 | 4 < A < 209 | http://dx.doi.org/10.1016/0370-1573(91)90039-O";
    if (name == "M") return "Menet et al. (1971) 30 <  E < 60 | 40 < A | http://dx/doi.org/10.1016/0092-640X(76)90007-3";
    if (name == "G") return "Becchetti and Greenlees (1969) E < 50 | 40 < A | http://dx.doi.org/10.1103/PhysRev.182.1190";
    if (name == "P") return "Perey (1963) E < 20 | 30 < A < 100 | http://dx/doi.org/10.1016/0092-640X(76)90007-3";

    //====== A = 3
    if (name == "x") return "XU, GUO, HAN, SHEN (2011) E < 250 | 20 < A < 209 | http://dx.doi.org/10.1007/s11433-011-4488-5";
    if (name == "l") return "Liang, Li, Cai (2009) E < 270 | All masses | http://dx.doi.org/10.1088/0954-3899/36/8/085104";
    if (name == "p") return "Pang et al., (2009) All E | All masses | Isospin dep. | http://dx.doi.org/10.1103/PhysRevC.79.024615";
    if (name == "c") return "Li, Liang, Cai, (2007) E < 40 | All masses | 48 < A < 232 | Tritons   | http://dx.doi.org/10.1016/j.nuclphysa.2007.03.004";
    if (name == "t") return "Trost et al., (1987) 10 < E < 220 | 10 < A < 208 | http://dx.doi.org/10.1016/0375-9474(87)90551-3";
    if (name == "h") return "Hyakutake et al., (1980) 90 < E < 120 | About 58 < A < 92 | http://dx.doi.org/10.1016/0375-9474(80)90013-5";
    if (name == "b") return "Becchetti and Greenlees, (1971) E < 40 | 40 < A | Iso. Dep.";

    //======= alpha
    if (name == "s") return "Su and Han, (2015) E < 398 | 20 < A < 209 | http://dx.doi/org/10.1142/S0218301315500925";
    if (name == "a") return "Avrigeanu et al., (2009) E ??? | A ??? | http://dx.doi/org/10.1016/j.adt.2009.02.001";
    if (name == "f") return "(FIXED) Bassani and Picard, (1969) 24 < E < 31 | A = 90 | https://doi.org/10.1016/0375-9474(69)90601-0";

    //====== custom
    if (name == "Y") return "Bardayan Parameters PRC 78 052801(R) (2008)";
    if (name == "X") return "Bardayan Parameters PRC 78 052801(R) (2008)";

    return "";
}

// ====================================================================
// custom
// ====================================================================
static OMPset CustomXPotential(int /*A*/, int /*Z*/, double /*E*/) {
    OMPset p;
    p.v = 54.19; p.r0 = 1.25; p.a = 0.65;
    p.vi = 0.0;  p.ri0 = 1.25; p.ai = 0.65;
    p.vsi = 13.5; p.rsi0 = 1.25; p.asi = 0.47;
    p.vso = 7.5;  p.rso0 = 1.25; p.aso = 0.47;
    p.vsoi = 0.0; p.rsoi0 = 1.25; p.asoi = 0.47;
    p.rc0 = 1.25;
    p.ok = true;
    return p;
}

static OMPset CustomYPotential(int /*A*/, int /*Z*/, double /*E*/) {
    OMPset p;
    p.v = 85.31; p.r0 = 1.15; p.a = 0.81;
    p.vi = 0.0;  p.ri0 = 1.15; p.ai = 0.81;
    p.vsi = 16.0; p.rsi0 = 1.34; p.asi = 0.68;
    p.vso = 0.0;  p.rso0 = 1.0;  p.aso = 1.0;
    p.vsoi = 0.0; p.rsoi0 = 1.0; p.asoi = 1.0;
    p.rc0 = 1.15;
    p.ok = true;
    return p;
}

// ====================================================================
// deuteron
// ====================================================================
static OMPset AnCaiPotential(int A, int Z, double E) {
    OMPset p;
    double A3 = pow(A, 1./3.);

    p.v  = 91.85 - 0.249*E + 0.000116*pow(E,2) + 0.642 * Z / A3;
    p.r0 = 1.152 - 0.00776 / A3;
    p.a  = 0.719 + 0.0126 * A3;

    p.vi  = 1.104 + 0.0622 * E;
    p.ri0 = 1.305 + 0.0997 / A3;
    p.ai  = 0.855 - 0.1 * A3;

    p.vsi  = 10.83 - 0.0306 * E;
    p.rsi0 = 1.334 + 0.152 / A3;
    p.asi  = 0.531 + 0.062 * A3;

    p.vso  = 3.557;
    p.rso0 = 0.972;
    p.aso  = 1.011;

    p.vsoi  = 0.0;
    p.rsoi0 = 0.0;
    p.asoi  = 0.0;

    p.rc0 = 1.303;
    p.ok = true;
    return p;
}

static OMPset HSSPotential(int A, int Z, double E) {
    OMPset p;
    int N = A-Z;
    double A3 = pow(A, 1./3.);

    double VIcond = -4.916 + (0.0555*E) + 0.0000442 * pow(E,2) +  35. * (N-Z)/A;

    p.v  = 82.18 - 0.148 * E - 0.000886 * pow(E,2) - 34.811*(N-Z)/A + 1.058*Z/A3;
    p.r0 = 1.174;
    p.a  = 0.809;

    p.vi  = VIcond > 0 ? VIcond : 0 ;
    p.ri0 = VIcond > 0 ? 1.563 : 0;
    p.ai  = VIcond > 0 ? 0.7 + 0.045 * A3 : 0;

    p.vsi  = 20.968 - 0.0794 * E - 43.398 * (N-Z)/A;
    p.rsi0 = 1.328;
    p.asi  = 0.465 + 0.045*A3;

    p.vso  = 3.703;
    p.rso0 = 1.234;
    p.aso  = 0.813;

    p.vsoi  = -0.206;
    p.rsoi0 = 1.234;
    p.asoi  = 0.813;

    p.rc0 = 1.698;
    p.ok = true;
    return p;
}

static OMPset BojowaldPotential(int A, int Z, double E) {
    OMPset p;
    double A3 = pow(A, 1./3.);

    double VIcond = 0.132 * (E - 45.);

    p.v  = 81.33 - 0.24 * E + 1.43 * Z / A3;
    p.r0 = 1.18;
    p.a  = 0.636 + 0.035 * A3;

    p.vi  = VIcond > 0 ? VIcond : 0 ;
    p.ri0 = 1.27;
    p.ai  = 0.768 + 0.021 * A3;

    p.vsi  = 7.8 + 1.04 * A3 - 0.712 * p.vi;
    p.rsi0 = 1.27;
    p.asi  = 0.768 + 0.021 * A3;

    p.vso  = 6.0;
    p.rso0 = 0.78 + 0.038 * A3;
    p.aso  = 0.78 + 0.038 * A3;

    p.vsoi  = 0;
    p.rsoi0 = 0;
    p.asoi  = 0;

    p.rc0 = 1.3;
    p.ok = true;
    return p;
}

static OMPset DaehnickPotential(int A, int Z, double E) {
    OMPset p;
    int N = A-Z;
    double A3 = pow(A, 1./3.);

    double beta = -1.* pow(E/100.,2);

    double MU1 = exp(-1*pow((8. - N)/2.,2));
    double MU2 = exp(-1*pow((20. - N)/2.,2));
    double MU3 = exp(-1*pow((28. - N)/2.,2));
    double MU4 = exp(-1*pow((50. - N)/2.,2));
    double MU5 = exp(-1*pow((82. - N)/2.,2));
    double MU6 = exp(-1*pow((126. - N)/2.,2));

    double MU = MU1 + MU2 + MU3 + MU4 + MU5 + MU6;

    p.v  = 88. - 0.283 * E + 0.88 * Z / A3;
    p.r0 = 1.17;
    p.a  = 0.717 + 0.0012 * E;

    p.vi  = (12. + 0.031 * E )*(1 - exp(beta));
    p.ri0 = 1.376 - 0.01 * sqrt(E);
    p.ai  = 0.52 + 0.07 * A3  - 0.04 * MU;

    p.vsi  = (12. + 0.031 * E) * exp(beta);
    p.rsi0 = p.ri0;
    p.asi  = p.ai;

    p.vso  = 7.2 - 0.032 * E;
    p.rso0 = 1.07;
    p.aso  = 0.66;

    p.vsoi  = 0;
    p.rsoi0 = 0;
    p.asoi  = 0;

    p.rc0 = 1.3;
    p.ok = true;
    return p;
}

static OMPset LohrPotential(int A, int Z, double E) {
    OMPset p;
    double A3 = pow(A, 1./3.);

    p.v  = 91.13 + 2.2 * Z / A3;
    p.r0 = 1.05;
    p.a  = 0.86;

    p.vi  = 0;
    p.ri0 = 0;
    p.ai  = 0;

    p.vsi  = 218./pow(A,2./3.);
    p.rsi0 = 1.43;
    p.asi  = 0.5 + 0.013 * pow(A,2./3.);

    p.vso  = 7;
    p.rso0 = 0.75;
    p.aso  = 0.5;

    p.vsoi  = 0;
    p.rsoi0 = 0;
    p.asoi  = 0;

    p.rc0 = 1.3;
    p.ok = true;
    return p;
}

static OMPset PereyPereyPotential(int A, int Z, double E) {
    OMPset p;
    double A3 = pow(A, 1./3.);

    p.v  = 81 - 0.22*E + 2 * Z / A3;
    p.r0 = 1.15;
    p.a  = 0.81;

    p.vi  = 0;
    p.ri0 = 0;
    p.ai  = 0;

    p.vsi  = 14.4 + 0.24 * E;
    p.rsi0 = 1.34;
    p.asi  = 0.68;

    p.vso  = 0.;
    p.rso0 = 0.;
    p.aso  = 0.;

    p.vsoi  = 0;
    p.rsoi0 = 0;
    p.asoi  = 0;

    p.rc0 = 1.15;
    p.ok = true;
    return p;
}

static OMPset ZhangPangLouPotential(int A, int Z, double E, double Zproj) {
    OMPset p;
    p.vso = 0;
    p.rso0 = 0;
    p.aso = 0;
    p.vsoi = 0;
    p.rsoi0 = 0;
    p.rc0 = 1.3;

    double A3 = pow(A, 1./3.);

    double RC = p.rc0 * A3;
    double Ec = 6 * Zproj * Z * 1.44 / 5 / RC;

    if (A == 6 && Z == 3) {
        p.r0  = 1.62 - 0.0122 * (E - Ec) / A3;
        p.ri0 = 2.83 - 0.0911 * (E - Ec) / A3;
        p.rsi0 = p.ri0;

        p.a   = 0.876;
        p.ai  = 0.27;
        p.asi = p.ai;

        p.v   = 47.9 + 2.37 * (E - Ec);
        p.vi  = 0;
        p.vsi = 11.3 + 3.44 * (E - Ec);

    } else if (A == 7 && Z == 3) {
        p.r0  = 1.45  + 0.097 * (E - Ec) / A3;
        p.ri0 = 2.12  + 0.022 * (E - Ec) / A3;
        p.rsi0 = p.ri0;

        p.a   = 0.844;
        p.ai  = 0.261;
        p.asi = p.ai;

        p.v   = 26.1 + 1.19 * (E - Ec);
        p.vi  = 0;
        p.vsi = 215.0 - 16.1 * (E - Ec);

    } else {
        p.r0   = 1.11  - 0.167 / A3 + 0.00117 * (E - Ec) / A3;
        p.ri0  = 0.561 + 3.07  / A3 - 0.00449 * (E - Ec) / A3;
        p.rsi0 = p.ri0;

        p.a   = 0.776;
        p.ai  = 0.744;
        p.asi = p.ai;

        p.v   = 98.9 - 0.279 * (E - Ec);
        p.vi  = 11.5 / ( 1 + exp((18.1 - (E - Ec))/5.97));
        p.vsi = 7.56 / ( 1 + exp(((E - Ec) - 14.3)/4.55));
    }

    p.ok = true;
    return p;
}

// ====================================================================
// proton
// ====================================================================
static OMPset KoningPotential(int A, int Z, double E, double Zproj) {
    OMPset p;
    int N   = A-Z;
    double A3 = pow(A, 1./3.);

    double vp1 = 59.3 + 21.*(N-Z)/A - 0.024*A;
    double vn1 = 59.3 - 21.*(N-Z)/A - 0.024*A;

    double vp2 = 0.007067 + 0.00000423*A;
    double vn2 = 0.007228 - 0.00000148*A;

    double vp3 = 0.00001729 + 0.00000001136 * A;
    double vn3 = 0.00001994 - 0.00000002 * A;

    double vp4 = 7e-9; // = vn4
    double vn4 = vp4;

    double wp1 = 14.667 + 0.009629*A;
    double wn1 = 12.195 + 0.0167*A;

    double wp2 = 73.55 + 0.0795*A; // = wn2
    double wn2 = wp2;

    double dp1 = 16 + 16.*(N-Z)/A;
    double dn1 = 16 - 16.*(N-Z)/A;

    double dp2 = 0.018 + 0.003802/(1 + exp((A-156.)/8)); // = dn2
    double dn2 = dp2;

    double dp3 = 11.5 ; // = dn3
    double dn3 = dp3;

    double vso1 = 5.922 + 0.003 * A;
    double vso2 = 0.004;

    double wso1 = -3.1;
    double wso2 = 160;

    double epf = -8.4075 + 0.01378 *A;
    double enf = -11.2814 + 0.02646 *A;

    double rc = 1.198 + 0.697/pow(A3,2) + 12.995/pow(A3,5);
    double vc = 1.73/rc * Z / A3;

    p.v  = vp1*(1 - vp2*(E-epf) + vp3*pow(E-epf,2) - vp4*pow(E-epf,3)) + vc * vp1 * (vp2 - 2*vp3*(E-epf) + 3*vp4*pow(E-epf,2));
    //neutron
    if (Zproj == 0) p.v  = vn1*(1 - vn2*(E-enf) + vn3*pow(E-enf,2) - vn4*pow(E-enf,3));

    p.r0 = 1.3039 - 0.4054 / A3;
    p.a  = 0.6778 - 0.000148 * A;

    p.vi  = wp1 * pow(E-epf,2)/(pow(E-epf,2) + pow(wp2,2));
    if (Zproj == 0) p.vi  = wn1 * pow(E-enf,2)/(pow(E-enf,2) + pow(wn2,2));

    p.ri0 = 1.3039 - 0.4054 / A3;
    p.ai  = 0.6778 - 0.000148 * A;

    p.vsi  = dp1 * pow(E-epf,2)/(pow(E-epf,2)+pow(dp3,2)) * exp(-dp2*(E-epf));
    if (Zproj == 0)   p.vsi  = dn1 * pow(E-enf,2)/(pow(E-enf,2)+pow(dn3,2)) * exp(-dn2*(E-enf));

    p.rsi0 = 1.3424 - 0.01585 * A3;
    p.asi  = 0.5187 + 0.0005205 * A;
    if (Zproj == 0) p.asi = 0.5446 - 0.0001656 * A;

    p.vso  = vso1 * exp(-vso2 * (E-epf));
    if (Zproj == 0) p.vso = vso1 * exp(-vso2 * (E-enf));

    p.rso0 = 1.1854 - 0.647/A3;
    p.aso  = 0.59;

    p.vsoi  = wso1 * pow(E-epf,2)/(pow(E-epf,2)+pow(wso2,2));
    if (Zproj == 0)   p.vsoi  = wso1 * pow(E-enf,2)/(pow(E-enf,2)+pow(wso2,2));

    p.rsoi0 = 1.1854 - 0.647/A3;
    p.asoi  = 0.59;

    p.rc0 = rc;
    p.ok = true;
    return p;
}

static OMPset VarnerPotential(int A, int Z, double E) {
    OMPset p;
    int N   = A-Z;
    double A3 = pow(A, 1./3.);

    double Rc = 1.24 * A3 + 0.12;
    double EC = 1.73 * Z / Rc;

    p.v  = 52.9 + 13.1 * (N-Z)/A - 0.299 * (E - EC);
    p.r0 = 1.25 - 0.225/A3;
    p.a  = 0.690;

    p.vi  = 7.8/(1 + exp((35. - E + EC)/16.) );
    p.ri0 = 1.33 - 0.42/A3;
    p.ai  = 0.69;

    p.vsi  = (10 + 18.*(N-Z)/A)/(1+exp((E-EC - 36.)/37.));
    p.rsi0 = 1.33 - 0.42/A3;
    p.asi  = 0.69;

    p.vso  = 5.9;
    p.rso0 = 1.34 - 1.2/A3;
    p.aso  = 0.63;

    p.vsoi  = 0;
    p.rsoi0 = 0;
    p.asoi  = 0;

    p.rc0 = Rc/A3;
    p.ok = true;
    return p;
}

static OMPset MenetPotential(int A, int Z, double E) {
    OMPset p;
    int N   = A-Z;
    double A3 = pow(A, 1./3.);

    p.v  = 49.9 - 0.22 * E + 26.4 * (N-Z) / A + 0.4 * Z / A3;
    p.r0 = 1.16;
    p.a  = 0.75;

    p.vi  = 1.2 + 0.09 * E;
    p.ri0 = 1.37;
    p.ai  = 0.74 - 0.008 * E + 1.*(N-Z)/A;

    p.vsi  = 4.2 - 0.05 * E + 15.5 * (N-Z)/A;
    p.rsi0 = 1.37;
    p.asi  = 0.74 - 0.008 * E + 1.*(N-Z)/A;

    p.vso  = 6.04;
    p.rso0 = 1.064;
    p.aso  = 0.78;

    p.vsoi  = 0;
    p.rsoi0 = 0;
    p.asoi  = 0;

    p.rc0 = 1.25;
    p.ok = true;
    return p;
}

static OMPset BecchettiPotential(int A, int Z, double E) {
    OMPset p;
    int N   = A-Z;
    double A3 = pow(A, 1./3.);

    p.v  = 54 - 0.32 * E + 24. * (N-Z) / A + 0.4 * Z / A3;
    p.r0 = 1.17;
    p.a  = 0.75;

    p.vi  = 0.22 * E - 2.7 < 0 ? 0 : 0.22 * E - 2.7;
    p.ri0 = 1.32;
    p.ai  = 0.51 + 0.7 * (N-Z)/A;

    p.vsi = 11.8 - 0.258 * E + 12. * (N-Z)/A;
    if (p.vsi < 0) {
        p.vsi = 0;
    }
    p.rsi0 = 1.320;
    p.asi  = 0.51 + 0.7 * (N-Z)/A;

    p.vso  = 6.2;
    p.rso0 = 1.1;
    p.aso  = 0.75;

    p.vsoi  = 0;
    p.rsoi0 = 0;
    p.asoi  = 0;

    p.rc0 = 1.3;
    p.ok = true;
    return p;
}

static OMPset PereyPotential(int A, int Z, double E) {
    OMPset p;
    int N   = A-Z;
    double A3 = pow(A, 1./3.);

    p.v  = 53.3 - 0.55 * E + 27. * (N-Z) / A + 0.4 * Z / A3;
    p.r0 = 1.25;
    p.a  = 0.65;

    p.vi  = 0.;
    p.ri0 = 0.;
    p.ai  = 0.;

    p.vsi  = 13.5;
    p.rsi0 =  1.25;
    p.asi  = 0.47;

    p.vso  = 7.5;
    p.rso0 = 1.25;
    p.aso  = 0.47;

    p.vsoi  = 0;
    p.rsoi0 = 0;
    p.asoi  = 0;

    p.rc0 = 1.25;
    p.ok = true;
    return p;
}

// ====================================================================
// A = 3
// ====================================================================
static OMPset XuPotential(int A, int Z, double E) {
    OMPset p;
    int N   = A-Z;
    double A3 = pow(A, 1./3.);

    double vsiTest = 33.26647 - 0.16975 * E - 12.0 * (N-Z)/A;
    double viTest = -2 + 0.10645*E - 0.00016156 * pow(E,2);

    p.v  = 136.34988 - 0.20315 * E - 0.00030147 * pow(E,2) - 24.0 * (N-Z) / A + 0.4 * Z / A3;
    p.r0 = 1.14963;
    p.a  = 0.78836;

    p.vi = 0.0; if (viTest > 0) p.vi = viTest;
    p.ri0 = 1.61807;
    p.ai  = 0.66485;

    p.vsi = 0.0; if (vsiTest > 0) p.vsi = vsiTest;
    p.rsi0 =  1.20655;
    p.asi  = 0.73593;

    p.vso  = 3.0;
    p.rso0 = 1.26864;
    p.aso  = 0.89999;

    p.vsoi  = 0;
    p.rsoi0 = 0;
    p.asoi  = 0;

    p.rc0 = 1.25;
    p.ok = true;
    return p;
}

static OMPset LiangPotential(int A, int Z, double E) {
    OMPset p;
    int N   = A-Z;
    double A3 = pow(A, 1./3.);

    p.v  = 118.36 - 0.2071 * E + 0.000063961 * pow(E,2) + 26.001 * (N-Z) / A + 0.5668 * Z / A3;
    p.r0 = 1.1657 + 0.0401 / A3;
    p.a  = 0.6641 + 0.0305 * A3;

    p.vi  = -6.8871 + 0.3115 * E - 0.00068096 * pow(E,2);
    p.ri0 = 1.4022 + 0.0418 / A3;
    p.ai  = 0.7732 + 0.0219 * A3;

    p.vsi  = 20.119 - 0.1626 * E - 5.4067 * (N-Z) / A + 1.2087 * A3;
    p.rsi0 = 1.1802 + 0.0587 / A3;
    p.asi  = 0.6292 + 0.0657 * A3;

    p.vso  = 2.0491 + 0.0099804 * A3;
    p.rso0 = 0.7211 + 0.0586 / A3;
    p.aso  = 0.7643 + 0.0535 * A3;

    p.vsoi  = -1.1591;
    p.rsoi0 = 0.7211 + 0.0586 * A3;
    p.asoi  = 0.7643 + 0.0535 * A3;

    p.rc0 = 1.289;
    p.ok = true;
    return p;
}

static OMPset PangPotential(int A, int Z, double E, int Zproj) {
    OMPset p;
    int N   = A-Z;
    double A3 = pow(A, 1./3.);

    double rc = 1.24 * A3 + 0.12;
    double EC = 1.728 * Z * Zproj / rc;
    double ETA = 1.0 * (N-Z) / A;
    double vsiAsym = 35.0 - 34.2 * ETA;
    if (Zproj == 2.0) vsiAsym = 35 + 34.2* ETA;

    p.v  = 118.3 - 0.13 * (E - EC);
    p.r0 = 1.3  - 0.48 / A3;
    p.a  = 0.82;

    p.vi  = 38.5/(1.0 + exp( (156.1 - E + EC)/52.4 ));
    p.ri0 = 1.31 - 0.13 / A3;
    p.ai  = 0.84;

    p.vsi  = vsiAsym / ( 1.0 + exp( (E - EC - 30.8) / 106.4 ));
    p.rsi0 = 1.31 - 0.13 / A3;
    p.asi  = 0.84;

    p.vso  = 0;
    if (E < 85) p.vso = 1.7 - 0.02 * E;
    p.rso0 = 0.64 + 1.18 / A3;
    p.aso  = 0.13;

    p.vsoi  = 0.;
    p.rsoi0 = 0.;
    p.asoi  = 0.;

    p.rc0 = rc/ A3;
    p.ok = true;
    return p;
}

static OMPset LiLiangCaiPotential(int A, int Z, double E) {
    OMPset p;
    int N   = A-Z;
    double A3 = pow(A, 1./3.);

    p.v  = 137.6 - 0.1456 * E + 0.0436 * pow(E,2) + 4.3751 * (N-Z) / A + 1.0474 * Z / A3;
    p.r0 = 1.1201 - 0.1504 / A3;
    p.a  = 0.6833 + 0.0191 * A3;

    p.vi  = 7.383 + 0.5025 * E - 0.0097 * pow(E,2);
    p.ri0 = 1.3202 - 0.1776 / A3;
    p.ai  = 1.119  + 0.01913 * A3;

    p.vsi  = 37.06 - 0.6451 * E - 47.19 * (N-Z) / A;
    p.rsi0 = 1.251 - 0.4622 / A3;
    p.asi  = 0.8114 + 0.01159 * A3;

    p.vso  = 1.9029;
    p.rso0 = 0.46991 + 0.1294 / A3;
    p.aso  = 0.3545 - 0.0522 * A3;

    p.vsoi  = 0.0;
    p.rsoi0 = 0.0;
    p.asoi  = 0.0;

    p.rc0 = 1.422;
    p.ok = true;
    return p;
}

static OMPset TrostPotential(int A, int Z, double E) {
    OMPset p;
    double A3 = pow(A, 1./3.);

    double JR = 272.33 * (1.0 + 0.002029 * A) * ( 1.0 - 0.001453 * E ) * (1.0 + 3.4931 * pow(A3,-2)) * ( 1.0+(-0.825165+exp(0.92059-0.079154*A))*exp(-0.065066*E) );
    double AP = 3 ;

    p.r0 = 1.150;
    double RR = A3 * p.r0;
    p.a  = 0.64*(1+(0.0004*A))*(1+0.25*(1-exp(-0.2*A))*(1-exp(-0.06*E)));
    double pi = 3.141592653589793;
    p.v  = 3.0/(4 * pi) * JR * AP * A * pow(RR, -3.0) / (1 + pow((pi * p.a)/RR,2) );

    p.vi  = 0.0;
    p.ri0 = 0.0;
    p.ai  = 0.0;

    p.vsi  = 24.5*(1+(1-(0.011*A))*(-0.0018)*E)*(1-exp(-(1.0+(0.003*E))*0.1*A))*(1-exp(-0.1*E));
    p.rsi0 = 1.26*(1-(0.00055*E))*(1+exp(-0.9163-(0.005*A))*exp(-0.09*E));
    p.asi  = 0.8;

    p.vso  = 0.0;
    p.rso0 = 0.0;
    p.aso  = 0.0;

    p.vsoi  = 0.0;
    p.rsoi0 = 0.0;
    p.asoi  = 0.0;

    p.rc0 = 1.4;
    p.ok = true;
    return p;
}

static OMPset HyakutakePotential(int A, int Z, double E) {
    OMPset p;
    int N   = A-Z;
    double A3 = pow(A, 1./3.);

    p.v  = 111.4 - 0.173 * E + 14.9 * (N-Z) / A + 1.1 * Z / A3;
    p.r0 = 1.21;
    p.a  = 0.76;

    p.vi  = 0.0;
    p.ri0 = 0.0;
    p.ai  = 0.0;

    p.vsi  = 24.8 - 0.028 * E;
    p.rsi0 = 1.17;
    p.asi  = 0.754 + 0.78 * (N-Z) / A;

    p.vso  = 0.0;
    p.rso0 = 0.0;
    p.aso  = 0.0;

    p.vsoi  = 0.0;
    p.rsoi0 = 0.0;
    p.asoi  = 0.0;

    p.rc0 = 1.300;
    p.ok = true;
    return p;
}

static OMPset BecchettiA3Potential(int A, int Z, double E, int Zproj) {
    OMPset p;
    int N   = A-Z;
    double A3 = pow(A, 1./3.);

    p.v  = 165 - 0.17 * E - 6.4 * (N-Z) / A;
    if (Zproj == 2) p.v = 151.9 - 0.17 * E + 50. * (N-Z) / A;

    p.r0 = 1.20;
    p.a  = 0.72;

    p.vi  = 46.0 - 0.33 * E - 110. * (N-Z) / A;
    if (Zproj == 2) p.vi = 41.7 - 0.33 * E + 44.0 * (N-Z) / A;
    p.ri0 = 1.4;
    p.ai  = 0.84;
    if (Zproj == 2) p.ai = 0.88;

    p.vsi  = 0.0;
    p.rsi0 = 0.0;
    p.asi  = 0.0;

    p.vso  = 2.5;
    p.rso0 = 1.2;
    p.aso  = 0.72;

    p.vsoi  = 0.0;
    p.rsoi0 = 0.0;
    p.asoi  = 0.0;

    p.rc0 = 1.300;
    p.ok = true;
    return p;
}

// ====================================================================
// alpha
// ====================================================================
static OMPset SuAndHanPotential(int A, int Z, double E) {
    OMPset p;
    int N   = A-Z;
    double A3 = pow(A, 1./3.);

    double vsiCOND = 27.5816 - 0.0797 * E + 48.0*(N-Z)/A;
    double  viCOND = -4.0174 + 0.1409 * E ;

    p.v  = 175.0881 - 0.6236 * E + 0.0006*E*E + 30.*(N-Z)/A - 0.236 * Z/A3;
    p.r0 = 1.3421;
    p.a  = 0.6578;

    p.vi  = viCOND; if (viCOND < 0) p.vi = 0.0;
    p.ri0 = 1.4259;
    p.ai  = 0.6578;

    p.vsi  = vsiCOND; if (vsiCOND < 0) p.vsi = 0.0;
    p.rsi0 = 1.2928;
    p.asi  = 0.6359;

    p.vso  = 0.0;
    p.rso0 = 1.2686;
    p.aso  = 0.85;

    p.vsoi  = 0.0;
    p.rsoi0 = 0.0;
    p.asoi  = 0.0;

    p.rc0 = 1.350;
    p.ok = true;
    return p;
}

static OMPset AvrigeanuPotential(int A, int Z, double E) {
    OMPset p;
    double A3 = pow(A, 1./3.);

    double e3 = 23.6 + 0.181 * Z/A3;
    double e2 = (2.59 + 10.4/A)*Z/(2.66+1.36*A3);
    double e1 = -3.03 + 0.762 * A3 + 0.24 + e2;
    (void)e1;

    p.v  = 116.5 + 0.337 * Z /A3 + 0.453*E;
    if (E < e3) p.v = 168 + 0.733 * Z / A3 - 2.64 *E;
    p.r0 = 0.00;
    p.a  = 0.00;

    p.vi  = 2.73 - 2.88 * A3 + 1.11*E;
    p.ri0 = 0.00;
    p.ai  = 0.00;

    p.vsi  = 0.0;
    p.rsi0 = 0.0;
    p.asi  = 0.0;

    p.vso  = 0.0;
    p.rso0 = 0.0;
    p.aso  = 0.0;

    p.vsoi  = 0.0;
    p.rsoi0 = 0.0;
    p.asoi  = 0.0;

    p.rc0 = 1.300;
    p.ok = true;
    return p;
}

static OMPset BassaniPicardPotential(int /*A*/, int /*Z*/, double /*E*/) {
    OMPset p;
    p.v  = 207;
    p.r0 = 1.3;
    p.a  = 0.65;

    p.vi  = 28;
    p.ri0 = 1.3;
    p.ai  = 0.52;

    p.vsi  = 0.0;
    p.rsi0 = 0.0;
    p.asi  = 0.0;

    p.vso  = 0.0;
    p.rso0 = 0.0;
    p.aso  = 0.0;

    p.vsoi  = 0.0;
    p.rsoi0 = 0.0;
    p.asoi  = 0.0;

    p.rc0 = 1.400;
    p.ok = true;
    return p;
}

// ====================================================================
// dispatcher
// ====================================================================
OMPset callPotential(const string& potName, int A, int Z, double E, int Zproj) {
    if (potName == "A") return AnCaiPotential(A, Z, E);
    if (potName == "H") return HSSPotential(A, Z, E);
    if (potName == "B") return BojowaldPotential(A, Z, E);
    if (potName == "D") return DaehnickPotential(A, Z, E);
    if (potName == "L") return LohrPotential(A, Z, E);
    if (potName == "Q") return PereyPereyPotential(A, Z, E);
    if (potName == "Z") return ZhangPangLouPotential(A, Z, E, Zproj);

    if (potName == "K") return KoningPotential(A, Z, E, Zproj);
    if (potName == "V") return VarnerPotential(A, Z, E);
    if (potName == "M") return MenetPotential(A, Z, E);
    if (potName == "G") return BecchettiPotential(A, Z, E);
    if (potName == "P") return PereyPotential(A, Z, E);

    if (potName == "x") return XuPotential(A, Z, E);
    if (potName == "l") return LiangPotential(A, Z, E);
    if (potName == "p") return PangPotential(A, Z, E, Zproj);
    if (potName == "c") return LiLiangCaiPotential(A, Z, E);
    if (potName == "t") return TrostPotential(A, Z, E);
    if (potName == "h") return HyakutakePotential(A, Z, E);
    if (potName == "b") return BecchettiA3Potential(A, Z, E, Zproj);

    if (potName == "s") return SuAndHanPotential(A, Z, E);
    if (potName == "a") return AvrigeanuPotential(A, Z, E);
    if (potName == "f") return BassaniPicardPotential(A, Z, E);

    if (potName == "X") return CustomXPotential(A, Z, E);
    if (potName == "Y") return CustomYPotential(A, Z, E);

    return OMPset{};  // ok = false
}
