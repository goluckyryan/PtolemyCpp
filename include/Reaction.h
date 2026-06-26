#pragma once
// reaction.h — Reaction class
// Pure particle kinematics: masses, charges, energies, spins, reaction label.
// This is the physics input — what particles are involved and at what energy.

#include "Kinematics.h"
#include "DWBAGrid.h"
#include "LinkuleData.h"
#include "BoundState.h"
#include "DistortedWave.h"
#include "InternalState.h"
#include "InelasticData.h"
#include "Timing.h"
#include "GridData.h"
#include "NamedStorage.h"
#include "ReactionParams.h"
#include "Spectroscopy.h"
#include "Energies.h"
#include "IntegrationGrid.h"
#include "OpticalPotentialParams.h"
#include "Masses.h"
#include "Charges.h"
#include "Flags.h"
#include "AngularMomentum.h"
#include <array>

class Reaction {
public:
    // Reaction identification
    char reactStr[46] = {};   // reaction string (NUL-padded)
    char header[66] = {};  // user header line (NUL-padded)

    // Kinematics
    Kinematics kin = {};

    // Orthogonal control flags (struct in Flags.h).
    Flags flags = {};

    // Angular-momentum / parity quantum numbers (struct in AngularMomentum.h).
    AngularMomentum angMom = {};

    // Reaction parameters
    ReactionParams rxn = {};

    // Spectroscopic / spin factors
    Spectroscopy spec = {};

    // Energies / Q-value / excitation
    Energies energies = {};

    // Integration-grid parameters
    IntegrationGrid integrationGrid = {};

    // Linkule user-parameter block (struct in LinkuleData.h).
    LinkuleParams linkuleParams = {};

    // Optical potential (Woods-Saxon)
    OpticalPotentialParams opticalPotentialParams = {};

    // Reaction masses (struct in Masses.h).
    Masses masses = {};

    // Nuclear charges (Z) — struct in Charges.h.
    Charges charges = {};

    // DWBA integration grid (std::vector storage).
    DWBAGrid dwbaGrid;

    // Linkule (linked-potential) address block. Linkules tie groups of
    // optical-potential parameters to a single scaling parameter for the fitter.
    LinkuleData linkuleData = {};

    // Bound-state class. Owns the form-factor work array (VPHI) and per-vertex
    // bound-state wavefunctions/potentials.
    BoundState boundState = {};

    // Distorted-wave class. Owns the two scattering channels (incoming/outgoing)
    // and the ScatteringSolver work arrays.
    DistortedWave distortedWave = {};

    // Internal calculation state. Sentinels (undefValue/notDefSentinel), per-calc
    // channel indices, energy bounds, completion flags, wasSet array, etc.
    InternalState internalState = {};

    // Inelastic / collective-DWBA scattering control. Holds the L_x sums,
    // S-matrix scratch (smatr/smati/sMag/sPhase), Coulomb integrals (CL1FF..),
    // and the betas / BETANRAT / ATERM vectors.
    InelasticData inelastic = {};

    // CPU timing / perf counters. Per-calculation timers; output in
    // input_reader's run summary.
    Timing timing = {};

    // DWBA integration grid: integration grids and pointer caches for the
    // inelastic/transfer DWBA integral. Separate from `dwbaGrid` above — gridData
    // holds the index/pointer layout, dwbaGrid holds the std::vector storage.
    GridData gridData = {};

    // Named-storage pool. Holds user-DEFINEd arrays (DEFINE input cards) and
    // linkule SHAPE potential-shape buffers, keyed by 8-char-style names.
    NamedStorage named = {};

    // Reaction setup — channel/potential and prologue helpers.
    void setChannel(int channelIndex, int& returnCode);
    void probePrint(int& returnCode);

    // Initialize all defaults.
    void applyDefaults();

    // Clear channel-specific parameters between channel calculations.
    void clearChannel(int channelIndex);

    // Parse a channel specification ("PROJ + TARG" or "PROJ + TARG -> ...")
    // from the input buffer. Returns true on success, false on parse error.
    bool setupChannel();

    // Parse the general "40CA(O16, 12C(2+ 3.26))44TI" reaction string
    // from the input buffer. Returns true on success, false on parse error.
    bool parseReactionString();

    // Apply a canned PARAMETERSET (e.g. CA60A, ALPHA1, EL1, INELOCA1) from
    // the input buffer to the reaction's grid/optical-potential fields.
    // Returns true on success, false if the data-name doesn't match.
    bool applyParameterSet();

    // Process a LINKULE input keyword: parse the linkule name from the input
    // buffer, look it up in the built-in linkule table, and store the address
    // in linkuleData.linkuleAddr for later dispatch by LINKUL.
    void loadLinkule(char8 linkuleKey);

    // Resolve and validate optical-potential parameters for one channel.
    // Reads from the just-parsed optical-potential block and fills in default
    // radii/diffusenesses / Coulomb radius using r0Mass. Returns true if every
    // required parameter is present, false if any of R/A/rC/ZP/ZT is missing.
    bool setupOpticalPotential();

    // Compute the full optical potential grid (real / imag / spin-orbit /
    // Coulomb / surface absorption) for channel channelIndexIn, populating the
    // rlvPointer/imvPointer/centPointer arrays on distortedWave.channel[channelIndex].
    void makePotential(int channelIndexIn, int& returnCode);

    // Set up the wavefunction-side Coulomb work arrays (F/G/NF1S/NG1S/sigma)
    // and S-matrix index/storage on the active waveChannel. Called once per
    // channel before the partial-wave scattering loop. Returns true on success,
    // false on RCWFN failure.
    bool setupWavefunctionPotential();

    // Three-pass build of the inelastic angular-momentum index tables (INDXS /
    // TOCS) on inelastic.indxsArr / tocsArr. pass=1 returns size into
    // counter; pass=2 fills INDXS and computes n_spl; pass=3 fills TOCS.
    // For CC, callers add the per-channel offset to the *Base pointers.
    void setupInelasticAngMomTable(int& counter, int* indxsBase, int* tocsBase, int pass);
};
