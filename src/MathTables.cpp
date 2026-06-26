// MathTables.cpp — static-local-init storage for the sqrt-factorial
// and log-factorial tables. The storage is hidden inside this TU;
// callers reach the tables through the factorialTable() /
// logFactorialTable() accessor functions declared in MathTables.h.
//
// Was previously a pair of header-extern globals (`FactorialTable
// factorialTable` and `LogFactorialTable logFactorialTable`) populated
// by FacdumInit / LogdumInit static-init structs in MathAngular.cpp.

#include "MathTables.h"

const FactorialTable& factorialTable() {
    static const FactorialTable t = []{
        FactorialTable f = {};
        f.maxFactorial = 96;
        f.factTable[ 1] =  0.1000000000000000e+01;
        f.factTable[ 2] =  0.1000000000000000e+01;
        f.factTable[ 3] =  0.1414213562373095e+01;
        f.factTable[ 4] =  0.2449489742783178e+01;
        f.factTable[ 5] =  0.4898979485566356e+01;
        f.factTable[ 6] =  0.1095445115010332e+02;
        f.factTable[ 7] =  0.2683281572999748e+02;
        f.factTable[ 8] =  0.7099295739719539e+02;
        f.factTable[ 9] =  0.2007984063681781e+03;
        f.factTable[10] =  0.6023952191045344e+03;
        f.factTable[11] =  0.1904940943966505e+04;
        f.factTable[12] =  0.6317974358922328e+04;
        f.factTable[13] =  0.2188610518114176e+05;
        f.factTable[14] =  0.7891147445080469e+05;
        f.factTable[15] =  0.2952597012800765e+06;
        f.factTable[16] =  0.1143535905863913e+07;
        f.factTable[17] =  0.4574143623455652e+07;
        f.factTable[18] =  0.1885967730625315e+08;
        f.factTable[19] =  0.8001483428544984e+08;
        f.factTable[20] =  0.3487765766344294e+09;
        f.factTable[21] =  0.1559776268628498e+10;
        f.factTable[22] =  0.7147792818185865e+10;
        f.factTable[23] =  0.3352612008237171e+11;
        f.factTable[24] =  0.1607856235454059e+12;
        f.factTable[25] =  0.7876854713229383e+12;
        f.factTable[26] =  0.3938427356614691e+13;
        f.factTable[27] =  0.2008211794424596e+14;
        f.factTable[28] =  0.1043497458090740e+15;
        f.factTable[29] =  0.5521669535672285e+15;
        f.factTable[30] =  0.2973510046012911e+16;
        f.factTable[31] =  0.1628658527169496e+17;
        f.factTable[32] =  0.9067986906793549e+17;
        f.factTable[33] =  0.5129628026803635e+18;
        f.factTable[34] =  0.2946746955341073e+19;
        f.factTable[35] =  0.1718233974287565e+20;
        f.factTable[36] =  0.1016520927791757e+21;
        f.factTable[37] =  0.6099125566750542e+21;
        f.factTable[38] =  0.3709953246501409e+22;
        f.factTable[39] =  0.2286968774309350e+23;
        f.factTable[40] =  0.1428211541796153e+24;
        f.factTable[41] =  0.9032802905233224e+24;
        f.factTable[42] =  0.5783815921445271e+25;
        f.factTable[43] =  0.3748341123420972e+26;
        f.factTable[44] =  0.2457951648494613e+27;
        f.factTable[45] =  0.1630420674178431e+28;
        f.factTable[46] =  0.1093719437815202e+29;
        f.factTable[47] =  0.7417966136220958e+29;
        f.factTable[48] =  0.5085501366740237e+30;
        f.factTable[49] =  0.3523338699662023e+31;
        f.factTable[50] =  0.2466337089763416e+32;
        f.factTable[51] =  0.1743963680863606e+33;
        f.factTable[52] =  0.1245439180886559e+34;
        f.factTable[53] =  0.8980989654316716e+34;
        f.factTable[54] =  0.6538259159791714e+35;
        f.factTable[55] =  0.4804619624270389e+36;
        f.factTable[56] =  0.3563201278858420e+37;
        f.factTable[57] =  0.2666455677120592e+38;
        f.factTable[58] =  0.2013129889124823e+39;
        f.factTable[59] =  0.1533154046820762e+40;
        f.factTable[60] =  0.1177637968756484e+41;
        f.factTable[61] =  0.9121944481710788e+41;
        f.factTable[62] =  0.7124466393192018e+42;
        f.factTable[63] =  0.5609810447812647e+43;
        f.factTable[64] =  0.4452649004137245e+44;
        f.factTable[65] =  0.3562119203309796e+45;
        f.factTable[66] =  0.2871872314724746e+46;
        f.factTable[67] =  0.2333120097803461e+47;
        f.factTable[68] =  0.1909741105966688e+48;
        f.factTable[69] =  0.1574812859496909e+49;
        f.factTable[70] =  0.1308137807832727e+50;
        f.factTable[71] =  0.1094466613011557e+51;
        f.factTable[72] =  0.9222139602976428e+51;
        f.factTable[73] =  0.7825244940376377e+52;
        f.factTable[74] =  0.6685892207860282e+53;
        f.factTable[75] =  0.5751421947239992e+54;
        f.factTable[76] =  0.4980877514193197e+55;
        f.factTable[77] =  0.4342228346904444e+56;
        f.factTable[78] =  0.3810289910601106e+57;
        f.factTable[79] =  0.3365156932181068e+58;
        f.factTable[80] =  0.2991016905800262e+59;
        f.factTable[81] =  0.2675246849288189e+60;
        f.factTable[82] =  0.2407722164359370e+61;
        f.factTable[83] =  0.2180285150390389e+62;
        f.factTable[84] =  0.1986334304622628e+63;
        f.factTable[85] =  0.1820505461284133e+64;
        f.factTable[86] =  0.1678423103505356e+65;
        f.factTable[87] =  0.1556505553593457e+66;
        f.factTable[88] =  0.1451811729660402e+67;
        f.factTable[89] =  0.1361920123419132e+68;
        f.factTable[90] =  0.1284832874770429e+69;
        f.factTable[91] =  0.1218899489080934e+70;
        f.factTable[92] =  0.1162756005221389e+71;
        f.factTable[93] =  0.1115276380752381e+72;
        f.factTable[94] =  0.1075533591796017e+73;
        f.factTable[95] =  0.1042768505784838e+74;
        f.factTable[96] =  0.1016365017512855e+75;
        f.factTable[97] =  0.9958302741285533e+75;
        return f;
    }();
    return t;
}

LogFactorialTable& logFactorialTable() {
    // Seeded with maxLf=0 and lf[1]=log(0!)=0.0; setLog() extends lf lazily.
    static LogFactorialTable t = []{
        LogFactorialTable l = {};
        l.maxLf = 0;
        l.lf[1] = 0.0;
        return l;
    }();
    return t;
}
