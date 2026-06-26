//
// Linkule for Argonne v18 deuteron wave function.
//
// May be used to calculate both wavefunction and potential,
// or just the effective potential. Uses tabulated S and D state
// wave functions and potentials with 3-point interpolation.
//
// 3/16/03 - new linkule based on Reid somewhat
//

#include "ptolemy_types.h"
#include <cstdio>
#include <cmath>
#include "Reaction.h"
#include "LinkulePlugin.h"

// AV18LinkulePlugin — Argonne v18 deuteron wavefunction/effective-potential
struct AV18LinkulePlugin : LinkulePlugin {
    void run(char8 alias, int* /*linkuleInts*/, int potType, int requestCode,
             int& callStatus, int L, double& J, double rStart, double stepSize,
             int nPts, double* array1, double* array2, Reaction& reaction) override;
};

void AV18LinkulePlugin::run(char8 alias, int* /*linkuleInts*/, int potType, int requestCode,
          int& callStatus, int L, double& J, double rStart,
          double stepSize, int nPts, double* array1, double* array2,
          Reaction& reaction)
{
    // implicit real*8 (a-h, o-z)
    int l = L;          // run() param name -> former av18 param name
    double& jp = J;     // bound-state j out-param

    static const char sd[3] = { ' ', 'S', 'D' };  // 1-based: sd[1]='S', sd[2]='D'

    // Grid parameters
    constexpr int numGrid = 241;
    constexpr double gridStep = 0.05;
    constexpr double gridMax = (numGrid - 1) * gridStep;

    // Tabulated wave functions: phis[j][caseIndex], 0-based j index, 1-based caseIndex
    // phis(0:num_grid-1, num_case)
    static real4 phis[3][numGrid] = {
        { 0 },  // [0] unused
        // S-state (caseIndex=1)
        {
            .08152f, .083498f, .089482f, .09958f, .1139f, .13244f, .15504f,
            .18127f, .21044f, .24155f, .2734f, .30473f, .33428f, .36098f, .38402f,
            .40291f, .41745f, .42774f, .43406f, .43685f, .43658f, .43378f,
            .42891f, .42243f, .4147f, .40605f, .39675f, .38701f, .37699f, .36684f,
            .35665f, .34652f, .3365f, .32664f, .31698f, .30754f, .29833f, .28938f,
            .28069f, .27226f, .2641f, .2562f, .24856f, .24118f, .23405f, .22716f,
            .22051f, .21409f, .2079f, .20192f, .19615f, .19058f, .18521f, .18002f,
            .17501f, .17018f, .16551f, .161f, .15665f, .15244f, .14837f, .14444f,
            .14064f, .13697f, .13341f, .12997f, .12664f, .12342f, .1203f, .11728f,
            .11436f, .11152f, .10878f, .10611f, .10353f, .10103f, .098599f,
            .096243f, .093957f, .091737f, .089583f, .087491f, .085459f, .083485f,
            .081567f, .079703f, .077891f, .07613f, .074417f, .072751f, .071131f,
            .069554f, .06802f, .066527f, .065073f, .063658f, .06228f, .060938f,
            .059631f, .058357f, .057116f, .055907f, .054728f, .053579f, .052459f,
            .051367f, .050301f, .049262f, .048249f, .047259f, .046295f, .045353f,
            .044434f, .043537f, .042661f, .041806f, .040971f, .040155f, .039359f,
            .038581f, .037821f, .037078f, .036353f, .035644f, .03495f, .034273f,
            .033611f, .032964f, .032331f, .031712f, .031107f, .030515f, .029936f,
            .029369f, .028815f, .028273f, .027743f, .027224f, .026716f, .026219f,
            .025733f, .025257f, .024791f, .024334f, .023887f, .02345f, .023022f,
            .022602f, .022191f, .021789f, .021395f, .021008f, .02063f, .02026f,
            .019897f, .019541f, .019192f, .018851f, .018516f, .018188f, .017866f,
            .017551f, .017242f, .016939f, .016642f, .016351f, .016065f, .015785f,
            .015511f, .015241f, .014977f, .014718f, .014465f, .014215f, .013971f,
            .013731f, .013496f, .013265f, .013039f, .012817f, .012599f, .012385f,
            .012175f, .01197f, .011767f, .011569f, .011375f, .011183f, .010996f,
            .010812f, .010631f, .010454f, .01028f, .010109f, .0099407f,
            .0097759f, .0096141f, .0094552f, .0092991f, .0091459f, .0089954f,
            .0088476f, .0087024f, .0085599f, .0084199f, .0082823f, .0081472f,
            .0080145f, .0078842f, .0077561f, .0076303f, .0075068f, .0073853f,
            .007266f, .0071488f, .0070336f, .0069205f, .0068093f, .0067f,
            .0065926f, .0064871f, .0063834f, .0062815f, .0061813f, .0060829f,
            .0059861f, .005891f, .0057976f, .0057057f, .0056154f, .0055266f,
            .0054393f, .0053535f, .0052692f, .0051862f, .0051047f, .0050246f,
            .0049457f, .0048683f, .0047921f, .0047172f
        },
        // D-state (caseIndex=2)
        {
            .0f, .0052663f, .019706f, .042113f, .072017f, .10918f, .15328f,
            .20356f, .25871f, .31689f, .37582f, .43311f, .48642f, .53383f,
            .57393f, .6059f, .62952f, .64505f, .65313f, .6546f, .65043f, .64161f,
            .62907f, .61368f, .59616f, .57715f, .55717f, .53665f, .51593f,
            .49526f, .47487f, .45488f, .43543f, .41658f, .39838f, .38088f,
            .36408f, .34798f, .33259f, .31789f, .30386f, .29049f, .27774f,
            .26561f, .25405f, .24305f, .23257f, .22261f, .21312f, .20409f, .1955f,
            .18732f, .17953f, .17211f, .16504f, .15831f, .15189f, .14578f,
            .13995f, .13439f, .12909f, .12403f, .1192f, .11459f, .11019f, .10599f,
            .10197f, .098132f, .094463f, .090953f, .087597f, .084385f, .081312f,
            .07837f, .075552f, .072854f, .070268f, .06779f, .065415f, .063137f,
            .060953f, .058857f, .056846f, .054916f, .053063f, .051283f, .049574f,
            .047931f, .046353f, .044835f, .043376f, .041973f, .040623f, .039324f,
            .038074f, .036871f, .035712f, .034596f, .033521f, .032485f, .031488f,
            .030526f, .029598f, .028704f, .027841f, .027009f, .026206f, .025431f,
            .024683f, .023961f, .023264f, .02259f, .021939f, .021311f, .020703f,
            .020116f, .019548f, .018999f, .018468f, .017955f, .017458f, .016977f,
            .016512f, .016062f, .015626f, .015204f, .014795f, .014399f, .014015f,
            .013644f, .013284f, .012935f, .012596f, .012268f, .01195f, .011642f,
            .011343f, .011052f, .010771f, .010498f, .010233f, .0099752f,
            .0097254f, .0094828f, .0092473f, .0090185f, .0087964f, .0085805f,
            .0083709f, .0081671f, .0079691f, .0077767f, .0075896f, .0074077f,
            .0072309f, .007059f, .0068918f, .0067291f, .0065709f, .006417f,
            .0062672f, .0061215f, .0059797f, .0058416f, .0057072f, .0055764f,
            .005449f, .005325f, .0052042f, .0050866f, .004972f, .0048603f,
            .0047516f, .0046456f, .0045424f, .0044417f, .0043436f, .004248f,
            .0041548f, .0040639f, .0039754f, .003889f, .0038047f, .0037225f,
            .0036424f, .0035642f, .0034879f, .0034135f, .0033409f, .00327f,
            .0032009f, .0031334f, .0030675f, .0030032f, .0029404f, .0028791f,
            .0028192f, .0027608f, .0027037f, .002648f, .0025935f, .0025404f,
            .0024884f, .0024377f, .0023881f, .0023396f, .0022923f, .002246f,
            .0022008f, .0021566f, .0021134f, .0020712f, .0020299f, .0019896f,
            .0019501f, .0019116f, .0018738f, .001837f, .0018009f, .0017656f,
            .0017311f, .0016973f, .0016643f, .001632f, .0016004f, .0015695f,
            .0015392f, .0015096f, .0014806f, .0014523f, .0014245f, .0013974f,
            .0013708f, .0013448f, .0013193f, .0012943f, .0012699f, .001246f,
            .0012226f, .0011997f, .0011772f
        }
    };

    // Tabulated potentials: vs[caseIndex][j], 0-based j, 1-based caseIndex
    static real4 vs[3][numGrid] = {
        { 0 },  // [0] unused
        // S-state potential (caseIndex=1)
        {
            2408.f, 2368.4f, 2247.6f, 2055.6f, 1812.6f, 1541.7f, 1263.9f, 995.83f,
            749.81f, 533.59f, 351.2f, 203.56f, 89.222f, 5.0502f, -53.16f,
            -90.118f, -110.55f, -118.84f, -118.74f, -113.3f, -104.84f, -95.015f,
            -84.918f, -75.227f, -66.31f, -58.323f, -51.294f, -45.172f, -39.872f,
            -35.292f, -31.332f, -27.901f, -24.917f, -22.312f, -20.028f, -18.018f,
            -16.24f, -14.664f, -13.262f, -12.011f, -10.892f, -9.8895f, -8.99f,
            -8.1816f, -7.454f, -6.7984f, -6.2072f, -5.6734f, -5.191f, -4.7547f,
            -4.3597f, -4.0017f, -3.677f, -3.3823f, -3.1144f, -2.8707f, -2.6489f,
            -2.4467f, -2.2622f, -2.0937f, -1.9397f, -1.7987f, -1.6695f, -1.5511f,
            -1.4423f, -1.3424f, -1.2504f, -1.1657f, -1.0877f, -1.0156f, -.94907f,
            -.88754f, -.8306f, -.77785f, -.72895f, -.68357f, -.64142f, -.60224f,
            -.56579f, -.53185f, -.50022f, -.47073f, -.44322f, -.41752f, -.3935f,
            -.37104f, -.35002f, -.33034f, -.31191f, -.29462f, -.27841f, -.26319f,
            -.24889f, -.23547f, -.22284f, -.21097f, -.19979f, -.18927f, -.17936f,
            -.17003f, -.16122f, -.15292f, -.14509f, -.1377f, -.13071f, -.12412f,
            -.11789f, -.11199f, -.10642f, -.10115f, -.096165f, -.091444f,
            -.086974f, -.082739f, -.078728f, -.074926f, -.071322f, -.067904f,
            -.064662f, -.061586f, -.058667f, -.055896f, -.053266f, -.050767f,
            -.048394f, -.046139f, -.043995f, -.041958f, -.040021f, -.038179f,
            -.036428f, -.034761f, -.033175f, -.031666f, -.030229f, -.028861f,
            -.027559f, -.026318f, -.025137f, -.024011f, -.022939f, -.021916f,
            -.020942f, -.020013f, -.019127f, -.018283f, -.017478f, -.016709f,
            -.015976f, -.015277f, -.01461f, -.013973f, -.013365f, -.012785f,
            -.012231f, -.011702f, -.011197f, -.010714f, -.010253f, -.0098131f,
            -.0093925f, -.0089906f, -.0086066f, -.0082396f, -.0078888f,
            -.0075535f, -.007233f, -.0069265f, -.0066335f, -.0063533f,
            -.0060853f, -.005829f, -.0055839f, -.0053494f, -.0051251f,
            -.0049104f, -.0047051f, -.0045086f, -.0043205f, -.0041405f,
            -.0039683f, -.0038034f, -.0036455f, -.0034944f, -.0033497f,
            -.0032112f, -.0030786f, -.0029516f, -.0028299f, -.0027134f,
            -.0026018f, -.002495f, -.0023926f, -.0022945f, -.0022005f,
            -.0021105f, -.0020242f, -.0019416f, -.0018624f, -.0017865f,
            -.0017137f, -.001644f, -.0015772f, -.0015131f, -.0014517f,
            -.0013928f, -.0013364f, -.0012823f, -.0012304f, -.0011807f,
            -.001133f, -.0010873f, -.0010434f, -.0010014f, -9.6106e-4f,
            -9.2238e-4f, -8.8528e-4f, -8.497e-4f, -8.1557e-4f, -7.8283e-4f,
            -7.5142e-4f, -7.2129e-4f, -6.9238e-4f, -6.6465e-4f, -6.3804e-4f,
            -6.1251e-4f, -5.8802e-4f, -5.6451e-4f, -5.4196e-4f, -5.2032e-4f,
            -4.9955e-4f, -4.7961e-4f, -4.6048e-4f, -4.4212e-4f, -4.245e-4f,
            -4.0759e-4f, -3.9136e-4f, -3.7577e-4f, -3.6082e-4f, -3.4646e-4f,
            -3.3268e-4f
        },
        // D-state potential (caseIndex=2)
        {
            1763.5f, -7088.1f, -2610.3f, -1239.7f, -678.63f, -451.01f, -384.97f,
            -398.84f, -447.45f, -504.36f, -554.5f, -590.14f, -608.53f, -610.1f,
            -597.15f, -572.89f, -540.68f, -503.63f, -464.35f, -424.85f, -386.57f,
            -350.43f, -316.96f, -286.41f, -258.8f, -234.02f, -211.89f, -192.19f,
            -174.67f, -159.09f, -145.25f, -132.92f, -121.94f, -112.12f, -103.34f,
            -95.457f, -88.365f, -81.968f, -76.182f, -70.934f, -66.161f, -61.809f,
            -57.831f, -54.185f, -50.837f, -47.755f, -44.911f, -42.283f, -39.85f,
            -37.592f, -35.494f, -33.542f, -31.723f, -30.024f, -28.437f, -26.953f,
            -25.562f, -24.257f, -23.032f, -21.882f, -20.799f, -19.781f, -18.821f,
            -17.915f, -17.061f, -16.254f, -15.492f, -14.771f, -14.088f, -13.442f,
            -12.83f, -12.249f, -11.699f, -11.176f, -10.68f, -10.209f, -9.7612f,
            -9.3353f, -8.9301f, -8.5446f, -8.1775f, -7.8279f, -7.4947f, -7.1772f,
            -6.8745f, -6.5858f, -6.3103f, -6.0474f, -5.7963f, -5.5566f, -5.3277f,
            -5.1089f, -4.8998f, -4.6999f, -4.5088f, -4.3259f, -4.1511f, -3.9837f,
            -3.8236f, -3.6702f, -3.5234f, -3.3829f, -3.2482f, -3.1192f, -2.9957f,
            -2.8772f, -2.7637f, -2.6549f, -2.5506f, -2.4506f, -2.3546f, -2.2626f,
            -2.1743f, -2.0896f, -2.0084f, -1.9304f, -1.8556f, -1.7837f, -1.7148f,
            -1.6486f, -1.585f, -1.524f, -1.4654f, -1.4091f, -1.355f, -1.3031f,
            -1.2532f, -1.2053f, -1.1592f, -1.115f, -1.0725f, -1.0316f, -.99234f,
            -.95459f, -.91831f, -.88344f, -.84991f, -.81768f, -.7867f, -.75691f,
            -.72826f, -.70072f, -.67423f, -.64876f, -.62426f, -.6007f, -.57805f,
            -.55625f, -.53529f, -.51513f, -.49573f, -.47708f, -.45913f, -.44186f,
            -.42525f, -.40926f, -.39389f, -.37909f, -.36485f, -.35116f, -.33798f,
            -.32529f, -.31309f, -.30134f, -.29004f, -.27916f, -.26869f, -.25862f,
            -.24892f, -.23959f, -.23061f, -.22196f, -.21364f, -.20564f, -.19793f,
            -.19051f, -.18337f, -.17649f, -.16987f, -.1635f, -.15737f, -.15147f,
            -.14579f, -.14032f, -.13505f, -.12998f, -.1251f, -.1204f, -.11588f,
            -.11153f, -.10733f, -.1033f, -.099412f, -.095671f, -.09207f,
            -.088603f, -.085264f, -.08205f, -.078956f, -.075976f, -.073108f,
            -.070346f, -.067687f, -.065126f, -.062661f, -.060288f, -.058003f,
            -.055803f, -.053684f, -.051645f, -.049681f, -.04779f, -.045969f,
            -.044216f, -.042528f, -.040903f, -.039339f, -.037832f, -.036382f,
            -.034985f, -.03364f, -.032346f, -.031099f, -.029899f, -.028743f,
            -.027631f, -.026559f, -.025528f, -.024535f, -.023579f, -.022658f,
            -.021772f, -.020919f, -.020097f, -.019307f, -.018545f, -.017812f,
            -.017106f, -.016427f, -.015773f, -.015143f
        }
    };

    static real4 spams[3]    = { 0.f, .97069f, -.24034f };     // 1-based
    static real4 eBounds[3]  = { 0.f, -2.22457f, -2.22457f };  // 1-based
    static real4 phiTails[3][3] = {
        { 0 },
        { 0.f, .21114f, .2316f },      // phiTails(:,1)
        { 0.f, .021156f, .23113f }      // phiTails(:,2)
    };
    static real4 vTails[3][3] = {
        { 0 },
        { 0.f, -18.019f, .72772f },     // vTails(:,1)
        { 0.f, -656.41f, .71085f }       // vTails(:,2)
    };
    static real4 vCouls[3] = { 0.f, 0.f, 0.f };  // 1-based

    // SAVE (all statics)

    // Skip to setup, printing, or calculation code.

    callStatus = 0;

    // Get the case number which points to the tables above
    int caseIndex = l / 2 + 1;
    double undef = reaction.internalState.undefValue;

    if (requestCode < 2) {
    // =========================================================================
    // Setup (requestCode = 1): check for error in l
    // =========================================================================
    if (l != 0 && l != 2) {
        printf("\n****** error in linkule %.8s\n"
               " ******  L =%5d   L must = 0 or 2\n", alias.data, l);
        callStatus = -1;
        return;
    }

    // Set jp, nodes, v, r, a, e, j, jsp, jst, spam
    // All arrays are 1-based flat copies: arr[N] = Fortran arr(N)
    jp = l + 1;
    reaction.angMom.nNodes = 0;
    // The potential depth is always set to 1.0
    reaction.opticalPotentialParams.V = 1.0;          // flt[62] = V
    if (reaction.integrationGrid.R == undef) reaction.integrationGrid.R = 1.0;    // flt[43] = R
    if (reaction.opticalPotentialParams.A == undef) reaction.opticalPotentialParams.A = 0.4;    // flt[1] = A
    if (reaction.energies.E == undef) reaction.energies.E = eBounds[caseIndex]; // flt[12] = E
    reaction.angMom.J = 2.0;                // jblock[1] = J
    reaction.angMom.spinProj = 1.0;              // jblock[8] = spinProj
    reaction.angMom.spinTarget = 1.0;              // jblock[9] = spinTarget
    // The default spec. amplitudes
    if (reaction.spec.spam == undef) reaction.spec.spam = spams[caseIndex];
    // "spamp" or "spamt", as appropriate, is always set to "spam".
    {
        int channelIndex = reaction.flags.hasNextBlock;
        if (channelIndex == 1) reaction.spec.specAmpProj = reaction.spec.spam;
        if (channelIndex == 2) reaction.spec.specAmpTgt = reaction.spec.spam;
    }
    return;

    } else if (requestCode == 2) {
    // =========================================================================
    // Printing (requestCode = 2)
    // =========================================================================
    if (potType == 6) {
        // Printout for wave function calculation.
        printf("%10sArgonne v18 deuteron wave function,%5s%c-state probability =%7.4f\n",
               "", "", sd[caseIndex], spams[caseIndex] * spams[caseIndex]);
    } else {
        // Printout for potential calculation
        printf("%10sEffective potential derived from Argonne v18 deuteron wavefunction,\n"
               "%5smultiplier =%9.6f\n", "", "", reaction.opticalPotentialParams.V);
    }

    // Printout for both wf and potential calculations.
    printf("%10sspectroscopic amplitude (\"spam\") =%8.5f\n", "", reaction.spec.spam);
    return;

    }
    // =========================================================================
    // Calculation (requestCode = 3)
    // =========================================================================
    {
        bool wfSwitch = (potType == 6);

        double rt   = rStart;
        double rEnd = rStart + (nPts - 1) * stepSize;
        int endIndex = (int)((std::min(rEnd, gridMax) - rStart) / stepSize);
        double hi   = 1.0 / gridStep;

        // Loop through radii, up to the highest tabulated value.
        for (int i = 1; i <= endIndex; i++) {
            // Three-point interpolation
            double frac = rt * hi;
            int j = (int)(frac + 0.5);
            frac = frac - j;
            double fracSquared = frac * frac;
            double x0, xPlus, xMinus;

            if (j == 0) {
                j = 1;
                x0  = 2.0 * frac - fracSquared;
                xPlus = -0.5 * (frac - fracSquared);
                xMinus = 1.0 - 1.5 * frac + 0.5 * fracSquared;
            } else {
                x0  = 1.0 - fracSquared;
                xPlus = 0.5 * (fracSquared + frac);
                xMinus = 0.5 * (fracSquared - frac);
            }

            if (wfSwitch) {
                double phi = x0 * phis[caseIndex][j] + xPlus * phis[caseIndex][j + 1]
                           + xMinus * phis[caseIndex][j - 1];
                array1[i] = rt * phi;
            }

            double vPot = x0 * vs[caseIndex][j] + xPlus * vs[caseIndex][j + 1]
                         + xMinus * vs[caseIndex][j - 1];
            array2[i] = -vPot;

            rt = rt + stepSize;
        }

        // Add the tails
        for (int i = endIndex + 1; i <= nPts; i++) {
            if (wfSwitch) {
                double x  = rt * phiTails[caseIndex][2];
                double xi = 1.0 / x;
                double f  = std::exp(-x) * xi;
                if (l == 1) f = (1.0 + xi) * f;
                if (l == 2) f = (1.0 + 3.0 * xi + 3.0 * xi * xi) * f;
                array1[i] = rt * phiTails[caseIndex][1] * f;
            }

            double x     = rt * vTails[caseIndex][2];
            double vPot  = vTails[caseIndex][1] * std::exp(-x) / x;
            vPot = vPot + vCouls[caseIndex] / rt;
            array2[i] = -vPot;

            rt = rt + stepSize;
        }
    }

    return;
}

std::unique_ptr<LinkulePlugin> makeAV18Plugin() {
    return std::make_unique<AV18LinkulePlugin>();
}
