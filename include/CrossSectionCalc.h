#pragma once
// cross_section_calc.h — CrossSectionCalc class
// Holds the cross-section calculation methods (angleSet/xSection/linterp/...).
// Phase-4/5 array migration done elsewhere: SMATR/SMATI live on
// inelastic.smatrArr/smatiArr; H on gridData.hsA12Arr; HINTEGRL/HABSINT
// on dwbaGrid.hint/habs; XLAMBDA/NLAMBDA on dwbaGrid.xlam/nlam. This class
// is methods-only — no owned storage.
//
// Thread safety: not thread-safe (ptolemy is single-threaded).

#include <vector>

class Reaction;

class CrossSectionCalc {
public:
    explicit CrossSectionCalc(Reaction& reaction) : reaction_(reaction) {}

    Reaction& reaction_;

    // -----------------------------------------------------------------------
    // -----------------------------------------------------------------------
    void angleSet();
    void xSection();
    void linterp();
    void ampCalc(double angle, int isElastic, int nSpline, int lMn, int lMx, int lSkp,
                int lxMax, int identicalParticles, double eta, double kWave, double sigZero,
                int returnFLowHigh, int lHigh,
                int* jtocs, float* betas, double* aLowFc,
                double* F, double* fLow, double* fHigh, double* fError,
                double* fCoul, double* plm, double* contR, double* contI,
                double* fEpsLow, int& flopCount);
    void phasePrint(bool printSwitch);
    void analyzingPower(double angleMin, double angleMax, double angleStep, int jA, int jB,
                int jResidual, int jBigB, int nSpline, int lxParity,
                int debugSwitch, double eLab, const char* channelName,
                const std::vector<double>& F_in, const int* tocsBase);
    void betCalc(int isElastic, double kWave, int spinProj, int nSpline, int lMn, int lMx, int lSkp,
                int lxMin, int lxMax, int lDeltaMax, int tempsCount, int statsCode, int printSwitch,
                double& sigTotal, double& sigReaction,
                int* jtocs, double* S, float* sMag, float* sPhase, double* sigIn, double* sigOut,
                float* betas, double* totLx, double* totMx, double* aLowFc, double* temps);
    void elasticDcs(double kWave, double eta, double angMinArg, double angMaxArg, double angStepArg,
               int lMn, int lMx, int lSkp, int statsCode, int spinProj,
               const double* smatBase, const int* tocsBase, int nSpline, const double* sigBase, double eLab, double tau,
               int printSwitch, std::vector<double>& torutOut,
               int keepFAmplitude, std::vector<double>& F_out, double& sigReaction,
               std::vector<double>& ruthOut, std::vector<double>& crossOut);
    void muElCoupling(int kA, int jA, int jB,
                int jResidual, int jBigB, int qaCount,
                int nSpline, int lxParity, int* tocsPointer, double* coefficient, double* coEp,
                double* coEt, double* coEx, double* termK, int debugSwitch);
};
