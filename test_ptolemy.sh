#!/bin/bash
# PtolemyCpp regression test suite vs Maple (source of truth).
#
# Stronger check (2026-06-04): primary regression criterion is now
# bit-identical *angular DCS curves*, not just the integrated total.
# This catches single-channel directionality bugs and per-angle drift
# that the previous total-only check would have missed.
#
# How the curve diff works:
#   1. Run PtolemyCpp + Maple on the same input
#   2. From each output, grep all lines matching the standard Ptolemy
#      DCS table format `^\s+<angle>\s+<dsigma/dOmega>\s+...`
#      (this captures every angle row across every section, for transfer
#      / inelastic / elastic / multi-LX outputs)
#   3. diff line-by-line; ANY difference == FAIL
#   4. Also extract the integrated total (legacy diagnostic, shown alongside)
#
# All 35 non-CC tests must produce bit-identical curves vs Maple.
# Maple is verified 35/35 bit-identical to 32-bit Cleopatra.
#
# Coupled-channels tests (cc_*.in) are skipped — PtolemyCpp has no CC support yet.

set -u
cd "$(dirname "$0")"

# Inputs are vendored in ./test_inputs (override with TESTDIR=...).
# The Maple oracle is for the bit-identical cross-check; advanced users can
# point Ptolemy_f2c=... at their own build. If it's missing we still run PtolemyCpp and
# just skip the curve diff (see HAVE_MAPLE below).
Ptolemy_f2c=${Ptolemy_f2c:-../Ptolemy-f2c/ptolemy}
PTOLEMY=./ptolemy
TESTDIR=${TESTDIR:-./test_inputs}
VERBOSE=${1:-""}

# XFAIL lists empty — every test passes bit-identically as of 2026-06-04.
XFAIL_CURVE=()

is_in_list() {
    local needle="$1"; shift
    for x in "$@"; do [ "$x" = "$needle" ] && return 0; done
    return 1
}

echo "Building PtolemyCpp..."
make -j$(nproc) 2>&1 | tail -1
echo ""

HAVE_MAPLE=1
if [ ! -x "$Ptolemy_f2c" ]; then
    HAVE_MAPLE=0
    echo "WARNING: Ptolemy_f2c binary not found at $Ptolemy_f2c."
    echo "         Running PtolemyCpp only — the bit-identical curve diff is SKIPPED."
    echo "         Set Ptolemy_f2c=/path/to/ptolemy to enable the cross-check."
    echo ""
fi

extract_total() {
    local out="$1"
    local v
    v=$(echo "$out" | grep "^0TOTAL:" | tail -1 | awk '{print $2}')
    [ -n "$v" ] && { echo "$v"; return; }
    v=$(echo "$out" | grep "TOTAL REACTION CROSS SECTION" | grep -o '[0-9][0-9]*\.[0-9]*')
    [ -n "$v" ] && { echo "$v"; return; }
    echo ""
}

# Pull every DCS-table data row from a Ptolemy output stream.
# Matches: leading whitespace, decimal angle, then at least one numeric column.
# Filters out the few cosmetic banner lines that can confuse a strict diff.
extract_curve() {
    grep -E "^\s+[0-9]+\.[0-9]+\s+[-+0-9.E ]" "$1" \
        | grep -v 'HEADER.*DESTROYED' \
        | grep -v 'OBJECT.*LOCATION'
}

PASS=0; XFAIL=0; REGRESS=0; ERR=0; SKIP=0; TOTAL=0; NODIFF=0
TOTAL_S=0; TOTAL_M=0

printf "%-45s %-15s %-15s %-12s %-7s %s\n" "Test" "PtolemyCpp" "Maple" "Time(s/m)" "ΔLines" "Status"
echo "-----------------------------------------------------------------------------------------------------------"
for f in $TESTDIR/*.in; do
    name=$(basename "$f" .in)
    case "$name" in cc_*) SKIP=$((SKIP+1)); continue;; esac
    TOTAL=$((TOTAL+1))

    rm -f fort.*
    t_s=$( { time timeout 60 $PTOLEMY < "$f" > /tmp/_s.out 2>/dev/null; } 2>&1 | grep real | sed 's/.*0m//;s/s//')
    s_rc=$?; rm -f fort.*
    if [ "$HAVE_MAPLE" = "1" ]; then
        t_m=$( { time timeout 60 $Ptolemy_f2c < "$f" > /tmp/_m.out 2>/dev/null; } 2>&1 | grep real | sed 's/.*0m//;s/s//')
        m_rc=$?; rm -f fort.*
    else
        t_m="-"; m_rc=0; : > /tmp/_m.out
    fi

    s_total=$(extract_total "$(cat /tmp/_s.out)")
    m_total=$(extract_total "$(cat /tmp/_m.out)")

    TOTAL_S=$(awk -v a=$TOTAL_S -v b=$t_s 'BEGIN{printf "%.3f", a+b}')
    TOTAL_M=$(awk -v a=$TOTAL_M -v b=$t_m 'BEGIN{printf "%.3f", a+b}')

    actual="UNKNOWN"; ndiff="-"
    if [ "$s_rc" = "124" ]; then
        actual="TIMEOUT"
    elif [ "$s_rc" = "134" ] || [ "$s_rc" = "139" ]; then
        actual="ABORT(rc=$s_rc)"
    elif [ ! -s /tmp/_s.out ]; then
        actual="NO_OUTPUT(rc=$s_rc)"
    elif [ "$HAVE_MAPLE" = "0" ]; then
        # No oracle available — PtolemyCpp ran and produced output; can't diff.
        actual="OK_NODIFF"
    elif [ ! -s /tmp/_m.out ]; then
        actual="MAPLE_ERR"
    else
        # The real test: bit-identical curve
        extract_curve /tmp/_s.out > /tmp/_s.curve
        extract_curve /tmp/_m.out > /tmp/_m.curve
        ndiff=$(diff /tmp/_s.curve /tmp/_m.curve | wc -l)
        if [ "$ndiff" = "0" ]; then
            actual="PASS"
        else
            actual="CURVE_DIFF"
        fi
    fi

    # Classify vs XFAIL list
    status="?"
    if [ "$actual" = "OK_NODIFF" ]; then
        status="OK(no-maple)"
        NODIFF=$((NODIFF+1))
    elif [ "$actual" = "PASS" ]; then
        if is_in_list "$name" "${XFAIL_CURVE[@]}"; then
            status="UNEXPECTED_PASS"
        else
            status="PASS"
        fi
        PASS=$((PASS+1))
    elif [ "$actual" = "CURVE_DIFF" ]; then
        if is_in_list "$name" "${XFAIL_CURVE[@]}"; then
            status="XFAIL($actual)"
            XFAIL=$((XFAIL+1))
        else
            status="REGRESSION($actual)"
            REGRESS=$((REGRESS+1))
        fi
    elif [[ "$actual" == ABORT* ]] || [ "$actual" = "TIMEOUT" ] || [[ "$actual" == NO_OUTPUT* ]]; then
        if is_in_list "$name" "${XFAIL_CURVE[@]}"; then
            status="XFAIL($actual)"
            XFAIL=$((XFAIL+1))
        else
            status="REGRESSION($actual)"
            REGRESS=$((REGRESS+1))
        fi
    else
        status="ERR($actual)"
        ERR=$((ERR+1))
    fi

    printf "%-45s %-15s %-15s %-12s %-7s %s\n" \
        "$name" "${s_total:--}" "${m_total:--}" "${t_s}s/${t_m}s" "$ndiff" "$status"

    # When verbose, print first 5 diff lines for any non-PASS
    if [ -n "$VERBOSE" ] && [ "$actual" != "PASS" ] && [ -s /tmp/_s.curve ] && [ -s /tmp/_m.curve ]; then
        echo "    --- first 5 diff lines (PtolemyCpp | Maple) ---"
        diff /tmp/_s.curve /tmp/_m.curve | head -10 | sed 's/^/      /'
    fi
done
echo "-----------------------------------------------------------------------------------------------------------"
if [ "$HAVE_MAPLE" = "1" ]; then
    printf "  Results: %d PASS, %d XFAIL, %d REGRESSION, %d ERR / %d (skipped %d CC)\n" \
        "$PASS" "$XFAIL" "$REGRESS" "$ERR" "$TOTAL" "$SKIP"
    printf "  Timing : PtolemyCpp %ss  Maple %ss\n" "$TOTAL_S" "$TOTAL_M"
    printf "  Check  : bit-identical angular DCS curves vs Cleopatra-verified Maple\n"
else
    printf "  Results: %d ran OK (no diff), %d REGRESSION, %d ERR / %d (skipped %d CC)\n" \
        "$NODIFF" "$REGRESS" "$ERR" "$TOTAL" "$SKIP"
    printf "  Timing : PtolemyCpp %ss  (Maple oracle absent)\n" "$TOTAL_S"
    printf "  Check  : PtolemyCpp ran without crashing; bit-identical diff SKIPPED (no Maple)\n"
fi
echo ""

if [ "$REGRESS" -gt 0 ] || [ "$ERR" -gt 0 ]; then
    echo "  ❌ REGRESSION DETECTED — revert before continuing."
    exit 1
fi
if [ "$HAVE_MAPLE" = "1" ]; then
    echo "  ✅ All currently-passing tests still pass curves bit-identical. $XFAIL XFAILs pending."
else
    echo "  ✅ PtolemyCpp ran on all $NODIFF inputs without crashing (no bit-identical check — Maple absent)."
fi
