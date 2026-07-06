# PtolemyCpp (d,t)/(p,d) Cross-Section Bug — Root Cause Analysis

**Date:** 2026-06-27  
**Analyser:** Master HELIOS subagent  
**Benchmark case:** `90Zr(d,t)89Zr(9/2+, Ex=3.14 MeV)`, Elab=37.6 MeV, l=4, n=0

---

## (a) First Divergent Numerical Quantity

The **first** observable difference between PtolemyCpp and f2c/Cleopatra is in the printed
optical-model potential for the **OUTGOING channel**:

| Code | Channel | Surface Absorption (vSi) | Radius (rSi) | Diffuseness (aSi) |
|------|---------|--------------------------|--------------|-------------------|
| f2c / Cleopatra | OUTGOING (t + ⁸⁹Zr) | **absent** (0.000) | — | — |
| PtolemyCpp      | OUTGOING (t + ⁸⁹Zr) | **9.6790 MeV** (WRONG) | 6.1078 fm | 0.8090 fm |

The input `OUTGOING` block explicitly specifies `vsi=0.000`.  
The `INCOMING` block specifies `vsi=9.679`.

PtolemyCpp applies the **INCOMING** surface imaginary to the **OUTGOING** channel,
adding ~12 MeV of spurious imaginary surface absorption to the triton+⁸⁹Zr potential.
This increases outgoing-channel absorption → **reduces** the transfer cross section.

**Numerical consequence:**

| Code | Total σ (mb/sr × … mb) | Ratio |
|------|------------------------|-------|
| f2c / Cleopatra | 0.052750 | 1.000 (reference) |
| PtolemyCpp      | 0.032715 | 0.6206 (factor ~1.61 too small) |

The error is uniform across all angles (verified by per-angle ratios 1.4–1.7), confirming
it is a normalization/potential error, not a computational instability.

**Confirming cross-check:** Running both codes with the same vSi in both channels
(`vsi=9.679` explicitly in the OUTGOING block) gives identical results (0.032715)
for both codes — proving the sole source of discrepancy is the erroneous vSi in the
outgoing channel of PtolemyCpp.

**Why only (d,t) and (p,d) — not (d,p):** The `(d,p)` outgoing potentials (proton OMP)
routinely have `vsi ≠ 0` (e.g., Koning-Delaroche has `vsi=5.659`). Since that value is
nonzero, `setOMparams` applies it correctly and overwrites the leaked INCOMING vSi.
The bug only manifests when the OUTGOING block has `vsi=0` (explicitly zero).

---

## (b) Root Cause

### Background: PtolemyCpp vs. f2c input processing

**In the original Ptolemy Fortran (f2c oracle):**  
The `CONTRL` state machine re-reads the input text sequentially. After processing the
`INCOMING` channel, `CLRCHN(3)` is called. For the DWBA case with `ICHAN=3` and
`IDONE=3`, the `goto L400` branch is taken — which does **not** clear `VSI` to 0.
However, the code then re-reads the `OUTGOING` block from the input text, and
`SRREAD` unconditionally assigns `VSI = 0.0` from the literal `vsi=0.000` in the text.
Thus `VSI` is correctly zero for the outgoing channel.

**In PtolemyCpp:**  
All input is pre-parsed once into a `ParsedInput` struct. The `OMParams::vSi` field for
the OUTGOING block is set to `0.0` (parsed from `vsi=0.000`). Then `setOMparams()` is
called to apply these params to `reaction.opticalPotentialParams`. The critical code:

```cpp
// src/InputParser.cpp, line 716:
if (om.vSi != 0.0) reaction.opticalPotentialParams.vSi = om.vSi;
```

Since `om.vSi == 0.0`, this condition is **false** and the assignment is **skipped**.
`reaction.opticalPotentialParams.vSi` retains its value from the INCOMING channel: **9.679**.

`clearChannel(3)` (called between INCOMING and OUTGOING) does not clear `vSi` either
— it goes through the `ICHAN>=3 && IDONE/4!=3` branch, matching the f2c behavior.
So after `clearChannel(3)`, `vSi = 9.679` persists. The f2c code recovers because it
re-reads the literal input and assigns 0 directly; PtolemyCpp has no such recovery.

### Exact divergent location in PtolemyCpp

**File:** `src/InputParser.cpp`  
**Function:** `setOMparams()`  
**Line:** 716  
**Code:**
```cpp
if (om.vSi != 0.0) reaction.opticalPotentialParams.vSi = om.vSi;
```
**Problem:** Treats `om.vSi == 0.0` as "not specified", but `0.0` is a valid
explicit user-provided value.

The `OMParams` struct (see `include/InputParser.h`, comment at line ~99):
```cpp
double rC0 = 0.0;   // 0 = not set (setOMparams skips when 0)
```
This "0 = not set" convention is safe for geometric parameters (zero radius/diffuseness
is unphysical) but **incorrect for absorptive strength terms** like `vSi`, `vI`, `vSo`,
`vSoi`, which are legitimately zero.

### Corresponding correct logic in f2c oracle

**File:** `Ptolemy-f2c/Cpp/src/contrl_translated.cpp` — CONTRL state machine  
**File:** `Ptolemy-f2c/Cpp/src/clrchn_translated.cpp` (line ~168): `VSI = 0;` (when the
condition `ICHAN>=3 && IDONE/4!=3` is NOT met — i.e., for non-scattering channels)  
**The key:** In f2c, `SRREAD` assigns `VSI` from input unconditionally. There is no
"skip if zero" guard. When `vsi=0` appears in the OUTGOING block, `VSI` is set to `0`.

---

## (c) Precise Code Fix

### Option A — Minimal fix (vSi only, targeted to known bug)

Add an explicit-set flag for `vSi` in `OMParams`:

**`include/InputParser.h`** — add `vSiSet` to `OMParams`:
```diff
     struct OMParams {
         double V    = 0.0;   double r0   = 0.0;   double a    = 0.0;
         double vI   = 0.0;   double rI0  = 0.0;   double aI   = 0.0;
         double vSi  = 0.0;   double rSi0 = 0.0;   double aSi  = 0.0;
         double vSo  = 0.0;   double rSo0 = 0.0;   double aSo  = 0.0;
         double vSoi = 0.0;   double rSoi0= 0.0;   double aSoi = 0.0;
         double rC0  = 0.0;   // 0 = not set (setOMparams skips when 0)
         bool   set  = false; // true if the block was present in input
+        bool   vSiSet = false;  // true if VSI was explicitly specified (even if 0)
     };
```

**`src/InputParser.cpp`** — set `vSiSet` when `VSI` is parsed:
```diff
-        if (getDouble(tok, "VSI",  dv)) { om.vSi  = dv; continue; }
+        if (getDouble(tok, "VSI",  dv)) { om.vSi  = dv; om.vSiSet = true; continue; }
```
and for the bare-key form:
```diff
-        if (tryBare("VSI",  om.vSi))  continue;
+        if (tryBare("VSI",  om.vSi))  { om.vSiSet = true; continue; }
```

**`src/InputParser.cpp`** — apply vSi when explicitly set or nonzero:
```diff
-    if (om.vSi  != 0.0) reaction.opticalPotentialParams.vSi  = om.vSi;
+    if (om.vSiSet || om.vSi != 0.0) reaction.opticalPotentialParams.vSi  = om.vSi;
```

### Option B — Complete fix (all absorptive/strength terms)

The same latent bug exists for `vI`, `vSo`, and `vSoi` — any of these could be
explicitly zeroed in a later channel and the zero would be dropped. The complete fix
follows the same pattern for all four:

```diff
 struct OMParams {
     double V    = 0.0;   double r0   = 0.0;   double a    = 0.0;
     double vI   = 0.0;   double rI0  = 0.0;   double aI   = 0.0;
     double vSi  = 0.0;   double rSi0 = 0.0;   double aSi  = 0.0;
     double vSo  = 0.0;   double rSo0 = 0.0;   double aSo  = 0.0;
     double vSoi = 0.0;   double rSoi0= 0.0;   double aSoi = 0.0;
     double rC0  = 0.0;
     bool   set  = false;
+    // Explicit-set flags for strength params that can be legitimately zero:
+    bool   vISet   = false;
+    bool   vSiSet  = false;
+    bool   vSoSet  = false;
+    bool   vSoiSet = false;
 };
```

In `parseKeyvals()`, set the flag alongside the value:
```diff
-    if (getDouble(tok, "VI",   dv)) { om.vI   = dv; continue; }
+    if (getDouble(tok, "VI",   dv)) { om.vI   = dv; om.vISet   = true; continue; }
-    if (getDouble(tok, "VSI",  dv)) { om.vSi  = dv; continue; }
+    if (getDouble(tok, "VSI",  dv)) { om.vSi  = dv; om.vSiSet  = true; continue; }
-    if (getDouble(tok, "VSO",  dv)) { om.vSo  = dv; continue; }
+    if (getDouble(tok, "VSO",  dv)) { om.vSo  = dv; om.vSoSet  = true; continue; }
-    if (getDouble(tok, "VSOI", dv)) { om.vSoi = dv; continue; }
+    if (getDouble(tok, "VSOI", dv)) { om.vSoi = dv; om.vSoiSet = true; continue; }
```
(and same for the `tryBare` counterparts)

In `setOMparams()`:
```diff
-    if (om.vI   != 0.0) reaction.opticalPotentialParams.vI   = om.vI;
+    if (om.vISet   || om.vI   != 0.0) reaction.opticalPotentialParams.vI   = om.vI;
-    if (om.vSi  != 0.0) reaction.opticalPotentialParams.vSi  = om.vSi;
+    if (om.vSiSet  || om.vSi  != 0.0) reaction.opticalPotentialParams.vSi  = om.vSi;
-    if (om.vSo  != 0.0) reaction.opticalPotentialParams.vSo  = om.vSo;
+    if (om.vSoSet  || om.vSo  != 0.0) reaction.opticalPotentialParams.vSo  = om.vSo;
-    if (om.vSoi != 0.0) reaction.opticalPotentialParams.vSoi = om.vSoi;
+    if (om.vSoiSet || om.vSoi != 0.0) reaction.opticalPotentialParams.vSoi = om.vSoi;
```

**Recommendation:** Implement Option B for completeness. The cost is 4 additional
bool fields in `OMParams` and 8 additional flag-set assignments in `parseKeyvals()`.

---

## (d) Confidence Level

**VERY HIGH (>99%).**

- The bug is confirmed by direct output comparison: PtolemyCpp shows  
  `SURFACE ABSORPTION  9.6790` in the OUTGOING channel; f2c does not.
- The physical direction is correct: extra imaginary surface absorption in the  
  outgoing channel reduces flux → cross section is smaller in PtolemyCpp. ✓
- Cross-check: when the OUTGOING block explicitly has `vsi=9.679` (matching  
  INCOMING), both codes give identical results (0.032715). ✓
- The code path in `setOMparams()` with the "skip-if-zero" guard at line 716 is  
  unambiguous and mechanically explains the failure mode.
- The Fortran oracle behavior (direct assignment from re-read input, no skip-if-zero)  
  is confirmed in `contrl_translated.cpp` and `clrchn_translated.cpp`.

---

## Summary Table

| Item | Value |
|------|-------|
| **First divergent quantity** | `vSi` in OUTGOING optical potential: PtolemyCpp=9.679, correct=0.000 |
| **Root cause file:line** | `src/InputParser.cpp:716`, function `setOMparams()` |
| **Bug mechanism** | `if (om.vSi != 0.0)` skips explicit zero values; outgoing inherits incoming vSi |
| **Fix** | Add `vSiSet` boolean; apply if `vSiSet || om.vSi != 0` |
| **Also affected** | `vI`, `vSo`, `vSoi` (same latent bug; fix all four in Option B) |
| **Reactions affected** | Any reaction where OUTGOING has `vsi=0` and INCOMING has `vsi≠0` |
| **Reactions unaffected** | (d,p): outgoing OMP has `vsi≠0` → no leak; (d,3He): same |
| **Confidence** | Very high (>99%) |

---

## FOLLOW-UP (2026-06-27): Second bug found in 3-way benchmark

After the setOMparams fix, a 3-way benchmark (Cpp vs f2c vs Cleopatra, 800 cases)
gave **792 PASS / 4 FAIL**. f2c-vs-Cleopatra disagreements = 0 (references agree).

The 4 remaining FAILs are ALL high-excitation 18O pickup reactions:
- 18O(d,t)17O ex=4.75 (cases 133, 537, 626)
- 18O(p,d)17O ex=4.75 (case 251)

### Root cause
When the excitation energy is ABOVE the neutron separation energy (Sn=4.143 MeV
for 18O), gen_input.py emits an explicit **fictitious binding energy** in the
TARGET block:

    TARGET
    nodes=0  l=1  jp=3/2  E=-.2  $ above Sn, fictitious binding
    ...

- **f2c / Cleopatra (correct):** honor the explicit `E=-0.2 MeV`
  (bound-state energy), KAPPA=0.09545 — a weakly-bound near-threshold form factor.
- **PtolemyCpp (wrong):** IGNORES the explicit `E=-.2` and runs its energy/well-depth
  search instead, landing on E=-12.794 MeV, KAPPA=0.76344 — a deeply-bound form
  factor. Total cross section 0.470 vs correct 0.115 (~4x too high).

Same FAMILY as the setOMparams bug: an explicitly-specified input field
(here the TARGET bound-state `E=`) is being treated as "unset" and overridden
by a default/search path.

### Where to look
InputParser TARGET-block parsing + BoundState well-depth/energy search.
Whatever decides "search for binding energy" vs "use the user-supplied E"
is not seeing the explicit E= from the TARGET block. Likely an analogous
`if (E != 0.0)` / "was-set" flag gap to the setOMparams fix.

### Reproduce
    18O(d,t)17O  E=68.7  ex=4.75  l=1  n=0  j=1.5
    -> f2c=0.11516  cleo=0.11515  cpp=0.47003  (WRONG)
Ground state (ex=0, bound) works: cpp=f2c=0.97143.

Confidence: high (root cause directly visible in bound-state printout: E and KAPPA differ).
