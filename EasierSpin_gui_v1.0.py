"""
EasySpin Input Builder v6.5 (single-file build)
================================================
A desktop GUI (customtkinter) that helps assemble Sys/Exp/Opt parameter
structures for EasySpin (the MATLAB EPR simulation toolbox) and generates
ready-to-run .m scripts, with dataset loading (Bruker .DSC + text data),
a Residual Builder, presets tied to compound IDs, hardened tooltips,
zoom hotkeys on the preview plot, a d-block electron-count/spin-state
helper, a full multi-component esfit (interactive fitting) block with
per-nucleus Vary A/weight and an auto-generated Supporting Information
summary + caption, a Line Broadening section, automatic LIVE
same-element isotope recognition, and a periodic-table isotope picker.

Supports: garlic (liquid), pepper (powder/frozen).

IMPORTANT ACCURACY NOTE
------------------------
Field names/units/behavior were cross-checked against EasySpin's own
official documentation pages (esfit.html and isotopetable.html,
retrieved 2026-08-17). Isotope spin/abundance values are transcribed
verbatim from EasySpin's own internal nuclear isotope database - see
chem_helpers.ISOTOPE_TABLE. The Supporting Information summary reads
fit_result.argsfit, whose exact shape for multi-component fits is
inferred from esfit's docs but not independently confirmed - wrapped in
try/catch so a mismatch doesn't break the CSV/MAT saves, just skips the
summary. Always sanity-check generated scripts against your installed
EasySpin version before relying on them for real work.

Run with:  python easyspin_gui_v6_5.py
Requires:  customtkinter, matplotlib   (pip install customtkinter matplotlib)

v6.5 changes:
  - NEW: Zero-field splitting on the Pepper tab ("Zero-Field Splitting
    (S >= 1)" section under g). Enter D (cm^-1, GHz or MHz) and E/D
    (0 to 1/3); the script writes Sys.D = [D E] in MHz - the documented
    2-element form (spinsystem.html) - with an explicit numeric unit
    conversion (1 cm^-1 = 29979.2458 MHz). Generation is blocked if S < 1
    or E/D is outside 0..1/3.
  - NEW: optional D-strain (Sys.DStrain = [FWHM_D FWHM_E], MHz, per
    broadenings.html), entered in the same unit as D.
  - NEW: Fit section gets "Vary ZFS D", "Vary ZFS E/D" (and "Vary
    D-strain") rows when ZFS is on. esfit can only vary Sys.D's real
    elements, so these are merged into SysVary.D = [dD dE] (MHz), with
    dE = d(E/D) x |D_start|. Like g, D is only varied on the first
    component unless that component has "Independent g" on.
  - NEW: Part 2's SI summary/caption reports fitted D (cm^-1) and E/D.
  - Presets save/load ZFS settings; older presets load with ZFS off.
  - NOT added to the Garlic tab: garlic.html restricts garlic to S=1/2
    doublets, where ZFS doesn't exist.
  - FIXED (caught while testing v6.5): since v6.3, the nucleus row's I=0
    warning label was created inside _toggle_override_fields() (an
    indentation slip) instead of __init__, so typing any isotope raised
    a background AttributeError and the "hide A for I=0 isotopes"
    behavior silently didn't run unless "Own starting g/lwpp" had been
    toggled first. Moved back into __init__.

v6.4 changes:
  - FIXED (real bug, reported directly): the periodic-table isotope
    picker's grab_set() (used to make it modal) was never released
    before destroy() - neither when picking an isotope nor when closing
    via the window's own X button. On some platforms that left the
    modal grab dangling, so the popup could appear to "not close" (or
    the whole app stop responding to clicks) even though destroy() had
    technically run. Both paths now release the grab first.
  - FIXED (real bug, reported directly): v6.3's "Own starting g/lwpp"
    checkbox only changed a component's STARTING g - it did nothing to
    let that g actually vary, since that was gated by a completely
    separate checkbox (v6.2's "Independent g", in a different section).
    Result: give a component its own starting g, generate, and the
    esfit table still showed only ONE g parameter. There's little
    reason to override a component's starting g and then force it
    right back onto the shared value during fitting anyway, so an
    active g override now automatically implies independent-g too - no
    need to find and check a second box in a different section for the
    common case. (The separate "Independent g" checkbox still exists
    for the narrower case of wanting a component to vary independently
    from the SAME starting point as everything else.)

v6.3 changes:
  - NEW: each nucleus row (below the isotope/n/A/weight fields) now has
    an "Own starting g/lwpp for this component" checkbox. Unchecked by
    default, so every component still just inherits the shared Sys.g/
    Sys.lwpp exactly as before (via "Sys2 = Sys;"). Checking it reveals
    two fields - g and lwpp (mT; type one number for a plain Gaussian,
    or two space-separated numbers for [Gaussian Lorentzian]) - and the
    generated script writes "Sys2.g = ...;"/"Sys2.lwpp = ...;" right
    after the "Sys2 = Sys;" copy line, overriding just that component's
    starting point.
  - This is the complementary piece to v6.2's "Independent g" (which
    lets an already-shared g float independently ONCE FITTING STARTS).
    v6.3 is about the STARTING value itself, before any fitting - for
    testing whether a shoulder/satellite might be a distinct species
    with its own g and/or linewidth from the outset, not just letting a
    shared starting point drift during optimization.
  - Persists through Save/Load preset (older presets without this key
    load with no override, same as leaving the checkbox unchecked -
    fully backward compatible).

v6.2 changes:
  - NEW: each nucleus row in the Fit section now has an "Independent g:
    <label>" checkbox (unchecked by default). The v6.0 fix made g free
    ONLY on the first component, since isotope satellites of the same
    species must share one true g and letting each vary independently
    was a real bug. But sometimes the point of a fit is exactly to
    TEST whether a component is a genuinely different species with
    its own g - not a same-species isotope satellite - and there was
    no way to do that without hand-editing the generated script every
    time. Checking this box for a component lets it keep g free (using
    the same delta as the main "Vary g-value" field) instead of having
    it stripped. This is a deliberate, visible opt-in - the default
    stays the safe v6.0 behavior - and the generated script includes
    an inline comment noting which component(s) have it on and why, so
    it's never a silent surprise later.

v6.1 changes:
  - FIXED (caught directly, twice, while fitting 47Ti hyperfine to a
    residual spectrum): the per-nucleus "Vary A" delta defaulted to an
    unusably tight +/-0.1 mT. If the true A is off from the starting
    guess by more than a fraction of an mT - likely, since there's
    rarely a converged starting value yet - esfit can never reach it.
    Default widened to +/-3 mT (still editable per row as always).
  - NEW: when fitting is enabled and a nucleus with real hyperfine
    (I>0) has "Vary A" left unchecked, the generated script now
    includes an inline comment warning that A stays fixed and esfit
    cannot determine the actual hyperfine splitting. This doesn't
    block generation - fixing A to a known literature value is a
    legitimate choice - it just stops the omission from being silent,
    since the generated script otherwise looks complete (Sys.A is
    still there, just never varied) and it's easy to miss that the
    checkbox was left unchecked.

v6.0 changes:
  - FIXED (real bug, caught while adding 47Ti/49Ti hyperfine satellite
    components to a Ti(III) garlic fit): the Fit section's global "Vary
    g-value(s)" (and, on the Pepper tab, "Vary g-strain") checkbox was
    being applied to EVERY generated Sys component, not just the first.
    Since all isotope components start as copies of the same base Sys,
    this seemed harmless, but esfit has no way to tie separate Sys/
    Sys2/Sys3 structs' parameters together during optimization - each
    SysVary.g entry is a fully independent free parameter once fitting
    starts. With multiple isotope components sharing one true
    electronic g-value (isotopes don't change g), that let low-weight
    satellite components' g drift away from the main component's g
    almost for free, since their small weight means their own fit
    quality barely affects the total residual - a likely cause of
    satellite lines converging far from their real position in
    multi-component fits. g/gStrain are now only left free on the
    FIRST component in both the Garlic and Pepper tabs; every other
    component keeps g/gStrain fixed at its copied starting value.
  - NOTE: APP_VERSION had been left at "5.5" in the title bar since
    that release despite v5.6-v5.9 shipping real functional changes
    (crop-control propagation, the Exp.Range/CenterSweep crash fix,
    SysVary.gStrain array-shape expansion, and the debounce/Copy-button
    fix) - this was never bumped alongside those fixes. Now corrected
    to display "6.0" to match this file.

v5.5 changes:
  - NEW: the Preview plot now has a "Crop field range (mT)" control - two
    numeric boxes (min/max) plus "Crop" and "Full view" buttons, right
    under the existing Z/F/H zoom hotkeys. Typing exact bounds is faster
    and more precise than dragging a zoom rectangle every time you want
    to inspect the same region, and the chosen crop is remembered and
    re-applied automatically the next time the plot redraws (e.g.
    reloading a dataset), instead of resetting to full view each time.

v5.4 changes:
  - FIXED (real bug, reported directly): selecting "axial" g symmetry in
    pepper still produced 3 independently-varying g parameters in esfit's
    table (g(1,1), g(1,2), g(1,3)) instead of 2. Root cause: axial mode
    expanded Sys.g to a 3-element "[gperp gperp gpar]" (duplicating
    gperp) rather than EasySpin's own native 2-element axial form
    "[gperp gpar]" (per garlic.html/pepper.html: Sys.g accepts 1, 2, or 3
    elements). Since SysVary has no way to tie two array elements
    together, that duplicated third slot let esfit vary g(1,1) and
    g(1,2) completely independently during fitting, silently breaking
    the axial symmetry the whole point of choosing "axial" was meant to
    enforce. Fixed by using the native 2-element form for both Sys.g and
    the matching SysVary.g - axial now genuinely has only 2 free g
    parameters, with no redundant slot for esfit to split away from.
  - NOTE (not yet fixed, flagging for visibility): Sys.gStrain is still
    always written as a hardcoded 3 values regardless of g mode. Per
    EasySpin's convention, gStrain should match Sys.g's own element count
    (1 for isotropic, 2 for axial, 3 for rhombic) - so enabling g-strain
    in isotropic or axial mode will likely hit a similar shape-mismatch
    error to the one just fixed. Isotropic-g pepper systems are
    uncommon, and this hasn't been reported yet, so it's flagged rather
    than rushed - ask if you hit it or want it fixed proactively.

v5.3 changes:
  - FIXED (real esfit error: "Field .A is missing in input 1."): v5.2
    made an I=0 nucleus's component omit Sys.Nucs/Sys.A entirely (correct
    - there's nothing to give it hyperfine for), but the Fit section's
    per-nucleus "Vary A" checkbox still applied to that same row, so
    checking it told esfit to vary a Sys.A field that plain doesn't
    exist on that component. The checkbox is now replaced with a plain
    "I=0, no A to vary" label for any such row - there's nothing left to
    check - and generate() has a matching safety net that refuses to
    emit a Vary spec for A on a zero-spin nucleus regardless of any
    saved checkbox state from before this fix.

v5.2 changes:
  - CHANGED: nucleus rows are single-line again (isotope/n/A(s)/weight%/
    "->100%"/remove all in one row) instead of v5.1's two-row layout -
    the two-row split was a workaround for the "->100%" button being
    unreachable, but felt like a confusing "bend". Fixed properly this
    time by shrinking widths and switching the row to a pack-based layout
    (needed anyway for the A-field show/hide below), which fits everything
    on one line without the button falling off-screen.
  - CHANGED: for an isotope with I=0 (no nuclear spin), the A field(s) are
    no longer just grayed out at a locked 0 - they're removed from the row
    entirely while that isotope is entered, so there's no box left for a
    value to be mistakenly typed into. Enter 48Ti/47Ti/49Ti as three rows
    and only 47Ti/49Ti (I>0) show an A field; 48Ti's row simply doesn't
    offer one. This also carries through to the generated MATLAB script:
    a component whose isotope is I=0 now omits Sys.Nucs/Sys.A entirely
    (previously it wrote a Sys.A of 0, which worked but wasn't necessary).
  - CHANGED: layout - the Preview plot is bigger (300px -> 440px tall) and
    the two script panels are shorter (280px -> 170px tall each), per
    direct request, so the plot is the visual focus with the scripts as
    supporting detail underneath.
  - REMOVED: the "Generate Script" button is gone. The script now updates
    live as you type/click, so a manual trigger button no longer serves a
    purpose. The live-status line under where the button was now shows
    validation feedback directly (green = valid, red = error message).
  - CHANGED: live regeneration is now event-driven instead of a blind
    700ms polling timer - it only re-runs shortly (~220ms debounce) after
    an actual change, rather than constantly regenerating whether or not
    anything changed. Covers text entries (KeyRelease), checkboxes (their
    BooleanVar), and buttons/dropdowns (wrapping their existing command),
    walked once at startup - plus per-row bindings on every nucleus row
    (isotope/n/A/weight) and explicit triggers after dataset load and
    preset load, so rows added later via "+Add nucleus" and programmatic
    field updates are covered too, not just the widgets present at
    startup.

v5.1 changes:
  (numbering note: v4.10 should have been called v5.0 given the scope of
  change since v4.x - this release is numbered v5.1 to correct that, no
  v5.0 was separately published)
  - FIXED: on a nucleus row, the "->100%" weight-balance button (and even
    the remove button) could end up past the visible/scrollable area on
    narrower windows, making them unreachable. The isotope/n/A fields now
    stay on their own row with tighter widths, and weight%/->100%/remove
    moved to a second row, left-aligned from column 0 - so those controls
    are always reachable regardless of window width.
  - NEW: nucleus rows with a spinless isotope (I=0, e.g. 48Ti) now FORCE
    their A field(s) to 0 and lock them (grayed out) automatically, rather
    than just warning. This means a mixed list like 48Ti + 47Ti + 49Ti
    works exactly as intended: A only needs to be filled in on the I>0
    rows (47Ti, 49Ti) - 48Ti's A is handled for you and can't be left at
    a stale nonzero value.
  - NEW: the script now regenerates automatically every ~700ms instead of
    only on clicking "Generate Script" - no more clicking the button after
    every small edit. Auto-regeneration reports problems in a small status
    line under the button (green when valid, red with the message when
    not) instead of a popup; the button itself still works for a manual
    refresh and is now the only path that can show a blocking dialog.
  - CHANGED: the Preview plot and both generated-script panels (Part 1,
    Part 2) are laid out together without an outer scrollbar - the plot
    sits on top, and Part 1/Part 2 sit side-by-side underneath, so both
    are visible at once without scrolling the whole page (matching the
    layout requested directly). The settings form on the left is still
    independently scrollable, since that list is long by nature.

v4.10 changes:
  - NEW: nucleus rows now warn live when the entered isotope has I=0 (no
    nuclear spin) - e.g. typing '48Ti' shows "has I=0 (no nuclear spin) -
    its A value(s) will have NO effect on the simulation" right under the
    row. Generating a script is also now BLOCKED if any nucleus with I=0
    has a nonzero A entered, with the same explanation plus guidance:
    zero it out (or leave it - EasySpin ignores it either way) if it's
    meant as the no-hyperfine background, or add a separate row for one
    of that element's I>0 isotopes if real hyperfine splitting is wanted.
    (Came up directly: a generated '48Ti' component had A values attached
    that were silently inert, since 48Ti's nuclear spin is I=0 - only
    47Ti/49Ti actually carry hyperfine for titanium.)

v4.9 changes:
  - CONFIRMED AGAINST OFFICIAL DOCS (garlic.html, pepper.html - the user
    provided these directly): pepper's "Anisotropic linewidth (3
    principal values)" option was based on a WRONG assumption from the
    start, not just its v4.7 combination with the Gaussian/Lorentzian
    split. garlic.html documents Sys.lwpp/lw (a field shared by garlic
    AND pepper) as accepting ONLY 1 element (Gaussian) or 2 elements
    ([Gaussian Lorentzian]) - there is no 3-element anisotropic form,
    ever. This matches the user's real EasySpin error ("System.lwpp has
    wrong size") even with the v4.8 "fix" that tried a plain 1x3 vector.
  - pepper.html's own worked example confirms the REAL mechanism for
    anisotropic/orientation-dependent broadening is Sys.HStrain (already
    present in this tool as a separate checkbox) - not Sys.lwpp.
  - Fix: removed the "Anisotropic linewidth" checkbox from the Pepper
    tab entirely. Pepper's Line Broadening section is now identical in
    shape to Garlic's: a base lwpp plus an optional Gaussian+Lorentzian
    split, with HStrain/gStrain as the separate (and correct) anisotropic
    mechanisms.

v4.8 changes:
  - CONFIRMED BUG (via real esfit/pepper error "System.lwpp has wrong
    size" on a generated script): the combined anisotropic + Gaussian/
    Lorentzian 3x2 Sys.lwpp form ("[Gxx Lxx; Gyy Lyy; Gzz Lzz]") that
    v4.1-v4.7 flagged as "not independently doc-verified" is confirmed
    INVALID against a real installed EasySpin (6.0.12). pepper's
    Sys.lwpp only accepts a scalar, a 1x2 [Gaussian Lorentzian] vector,
    or a 1x3 [lwxx lwyy lwzz] anisotropic vector - never a combined
    3x2 matrix.
  - Fix: "Anisotropic linewidth" and "Split into Gaussian + Lorentzian
    components" are now mutually exclusive checkboxes on the Pepper tab
    (checking one unchecks the other), so the invalid combined form can
    no longer be generated. If you need both anisotropic AND mixed-
    lineshape broadening, use Sys.HStrain/gStrain (Gaussian, additive)
    alongside a single scalar or 1x2 lwpp instead.

v4.7 changes:
  - Part 2 of the esfit workflow now builds a Supporting Information
    summary from the FITTED parameters (via fit_result.argsfit, not the
    initial guesses baked into the script) - number of components, each
    component's isotope/weight/S/g/lwpp(Gaussian+Lorentzian
    breakdown)/A, mwFreq, and Temperature if provided. This is written
    as '%'-prefixed comment lines at the top of the SAME CSV the trace
    is saved to (any reader that skips non-numeric header lines, and
    this tool's own Residual Builder, reads past them automatically).
  - Also prints a ready-to-edit SI figure caption in the format "Figure
    SXX. Solution X-band EPR spectrum of <COMPOUND> recorded at ... K
    (black) and the corresponding simulation (red) performed using the
    EasySpin <function> function (...)." to the MATLAB console.
  - Added an optional Temperature (K) field to the Garlic tab (Pepper
    already had one) - garlic likely doesn't need it for the simulation
    itself given AutoScale absorbs overall intensity effects, but it's
    standard to report for a solution EPR measurement and now flows
    into both Exp.Temperature and the SI summary/caption.
  - Fixed a real MATLAB syntax bug caught during testing: the SI caption
    sprintf had a string literal split across two generated lines
    without a closing quote or '...' continuation - would have caused
    an "unterminated string" error in MATLAB. Verified fixed via a
    bracket/quote balance check across all four generated script
    variants (garlic/pepper x Part 1/Part 2).

v4.6 changes:
  - Added "Vary weight" as a per-nucleus fitting option

v4.5 changes:
  - Fixed a real bug in the Residual Builder's Exp value handling

v4.4 changes - Residual Builder added, for isolating a minor component

v4.3 changes:
  - Fixed a spacing bug; pepper's lwpp now supports anisotropic +
    Gaussian/Lorentzian split, combinable

v4.2 changes:
  - Hardened tooltips against getting stuck open

v4.1 changes:
  - Line Broadening section added to both tabs

v4.0 changes - MAJOR: multi-component esfit fitting now fully supported

v3.3 changes:
  - Added a "->100%" button per Nuclei row to balance weights to 100%

v3.2 changes:
  - Same-element weight auto-fill is now LIVE

v3.1 changes:
  - Fixed weight auto-fill being blocked by full-row validation; fixed
    the periodic table picker opening behind the main window

v3.0 changes:
  - Removed the Chili (slow-motion) tab entirely

v2.9 changes:
  - Replaced the dedicated Isotope Mixture panel with automatic
    same-element recognition

v2.8 changes:
  - Every nucleus row is its own separate weighted Sys component

v2.7 changes:
  - Isotope Mixture feature added to the Garlic tab

v2.6 changes:
  - Isotope Mixture abundance percentages editable/overridable

v2.5 changes:
  - Explicit Target/weighting control added to Fit (esfit)

v2.4 changes:
  - Periodic-table isotope picker (click any isotope field)

v2.3 changes:
  - "Natural-abundance isotope mixture" section added to pepper/chili

v2.2 changes:
  - Hyperfine (A) fitting added to esfit

v2.1 changes:
  - Split esfit output into two scripts (Part 1 / Part 2)

v2.0 changes:
  - Post-fit export opens a native "Save As" dialog (uiputfile)

v1.9 changes:
  - Fixed tooltips floating on top of OTHER applications

v1.8 changes:
  - Post-fit export writes a proper CSV with labeled columns

v1.7 changes:
  - Post-export section uses esfit's REAL export naming convention

v1.6 changes:
  - Removed the pre-fit preview plot; added post-fit export section

v1.5 changes:
  - Fixed "Unrecognized field name 'Arg'"; fixed a mixed-separator path bug

v1.4 changes:
  - esfit block points to the experimental data FILE at MATLAB runtime

v1.3 changes:
  - Added an optional "Fit (esfit)" section on all three tabs

v1.2 changes:
  - Fixed a 200px row-height layout bug; tooltips added; electron-count
    suggestions now live; larger preview plot
"""

import os
import re
import json
import time
import datetime

import customtkinter as ctk
from tkinter import filedialog, messagebox
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg



##############################################################################
# chem_helpers.py content (d-block electron counting, official EasySpin isotope table, periodic table data, same-element auto-recognition)
##############################################################################

import re

# Periodic group number for d-block elements (3-12), used for d-electron
# counting. Only first-, second-, and third-row transition metals included.
D_BLOCK_GROUPS = {
    # First row
    "Sc": 3, "Ti": 4, "V": 5, "Cr": 6, "Mn": 7, "Fe": 8, "Co": 9,
    "Ni": 10, "Cu": 11, "Zn": 12,
    # Second row
    "Y": 3, "Zr": 4, "Nb": 5, "Mo": 6, "Tc": 7, "Ru": 8, "Rh": 9,
    "Pd": 10, "Ag": 11, "Cd": 12,
    # Third row
    "La": 3, "Hf": 4, "Ta": 5, "W": 6, "Re": 7, "Os": 8, "Ir": 9,
    "Pt": 10, "Au": 11, "Hg": 12,
}

# Ground-state-term f-electron counts for common lanthanides/actinides.
# Value = number of f electrons in the ion (element - oxidation-derived).
# This is used ONLY to compute a Hund's-rule S hint; see module docstring.
F_BLOCK_BASE_ELECTRONS = {
    # Lanthanides (neutral atom 4f count; typical +3 ion is base-2 relative
    # to neutral in simple counting, but exceptions are common -
    # this table is intentionally coarse and flagged as a hint only)
    "Ce": 1, "Pr": 2, "Nd": 3, "Pm": 4, "Sm": 5, "Eu": 6, "Gd": 7,
    "Tb": 8, "Dy": 9, "Ho": 10, "Er": 11, "Tm": 12, "Yb": 13, "Lu": 14,
    # Actinides (5f count for neutral atom, same caveat)
    "Th": 0, "Pa": 2, "U": 3, "Np": 4, "Pu": 5, "Am": 6, "Cm": 7,
    "Bk": 8, "Cf": 9, "Es": 10, "Fm": 11, "Md": 12, "No": 13, "Lr": 14,
}

F_BLOCK_ELEMENTS = set(F_BLOCK_BASE_ELECTRONS.keys())


def d_electron_count(element, oxidation_state):
    """
    Return the d-electron count for a transition-metal ion, or None if the
    element isn't in the d-block table.
    """
    group = D_BLOCK_GROUPS.get(element)
    if group is None:
        return None
    d_count = group - oxidation_state
    return d_count


# Mapping: d_count -> list of (S, label) for octahedral fields.
# "unambiguous" entries have exactly one item; d4-d7 and d8 have two.
_OCTAHEDRAL_TABLE = {
    0: [(0, "d0 - diamagnetic, no unpaired electrons")],
    1: [(0.5, "d1 - S=1/2, only one possibility")],
    2: [(1, "d2 - S=1, only one possibility")],
    3: [(1.5, "d3 - S=3/2, only one possibility")],
    4: [(2, "d4 high-spin, S=2"), (1, "d4 low-spin, S=1")],
    5: [(2.5, "d5 high-spin, S=5/2"), (0.5, "d5 low-spin, S=1/2")],
    6: [(2, "d6 high-spin, S=2"), (0, "d6 low-spin, S=0 (diamagnetic)")],
    7: [(1.5, "d7 high-spin, S=3/2"), (0.5, "d7 low-spin, S=1/2")],
    8: [(1, "d8 octahedral, S=1")],
    9: [(0.5, "d9 - S=1/2, only one possibility")],
    10: [(0, "d10 - diamagnetic, no unpaired electrons")],
}

_TETRAHEDRAL_TABLE = {
    0: [(0, "d0 - diamagnetic")],
    1: [(0.5, "d1 - S=1/2")],
    2: [(1, "d2 - S=1")],
    3: [(1.5, "d3 - S=3/2 (tetrahedral d3 is uncommon but this is the "
              "only spin value)")],
    4: [(2, "d4 tetrahedral high-spin (weak-field default), S=2"),
        (1, "d4 tetrahedral low-spin, S=1 (rare - needs an unusually "
            "strong-field ligand set)")],
    5: [(2.5, "d5 tetrahedral high-spin (weak-field default), S=5/2"),
        (0.5, "d5 tetrahedral low-spin, S=1/2 (rare - needs an unusually "
              "strong-field ligand set)")],
    6: [(2, "d6 tetrahedral high-spin (weak-field default), S=2"),
        (0, "d6 tetrahedral low-spin, S=0 (rare - needs an unusually "
            "strong-field ligand set)")],
    7: [(1.5, "d7 tetrahedral high-spin (weak-field default), S=3/2"),
        (0.5, "d7 tetrahedral low-spin, S=1/2 (rare - needs an unusually "
              "strong-field ligand set)")],
    8: [(1, "d8 tetrahedral - S=1")],
    9: [(0.5, "d9 - S=1/2")],
    10: [(0, "d10 - diamagnetic")],
}

_SQUARE_PLANAR_TABLE = {
    0: [(0, "d0 - diamagnetic")],
    1: [(0.5, "d1 - S=1/2")],
    2: [(1, "d2 square planar - S=1 (unusual geometry for d2)")],
    3: [(1.5, "d3 square planar - S=3/2 (unusual geometry for d3)")],
    4: [(2, "d4 square planar high-spin, S=2"),
        (1, "d4 square planar low-spin, S=1")],
    5: [(2.5, "d5 square planar high-spin, S=5/2"),
        (0.5, "d5 square planar low-spin, S=1/2")],
    6: [(2, "d6 square planar high-spin, S=2"),
        (0, "d6 square planar low-spin, S=0")],
    7: [(1.5, "d7 square planar high-spin, S=3/2"),
        (0.5, "d7 square planar low-spin, S=1/2")],
    8: [(0, "d8 square planar - low-spin, S=0 (diamagnetic, common for "
            "e.g. Ni(II), Pd(II), Pt(II))")],
    9: [(0.5, "d9 - S=1/2")],
    10: [(0, "d10 - diamagnetic")],
}

# Five-coordinate geometries: generally intermediate-field, closer to
# octahedral-like splitting than tetrahedral, but real splitting patterns
# vary more with the specific ligand set than for Oh/Td/D4h - HS is
# somewhat more common in practice, but both are offered since either can
# occur depending on the ligand field.
_TRIGONAL_BIPYRAMIDAL_TABLE = _OCTAHEDRAL_TABLE
_SQUARE_PYRAMIDAL_TABLE = _OCTAHEDRAL_TABLE

_GEOMETRY_TABLES = {
    "octahedral": _OCTAHEDRAL_TABLE,
    "tetrahedral": _TETRAHEDRAL_TABLE,
    "square_planar": _SQUARE_PLANAR_TABLE,
    "trigonal_bipyramidal": _TRIGONAL_BIPYRAMIDAL_TABLE,
    "square_pyramidal": _SQUARE_PYRAMIDAL_TABLE,
}

GEOMETRY_OPTIONS = [
    "any", "octahedral", "tetrahedral", "square_planar",
    "trigonal_bipyramidal", "square_pyramidal",
]

GEOMETRY_CAVEATS = {
    "trigonal_bipyramidal": (
        "5-coordinate geometry: crystal-field splitting patterns are more "
        "ligand-set-dependent than for octahedral/tetrahedral. Both "
        "high-spin and low-spin options are listed, but verify against "
        "literature for your specific ligand set."
    ),
    "square_pyramidal": (
        "5-coordinate geometry: crystal-field splitting patterns are more "
        "ligand-set-dependent than for octahedral/tetrahedral. Both "
        "high-spin and low-spin options are listed, but verify against "
        "literature for your specific ligand set."
    ),
}


def spin_state_suggestions(element, oxidation_state, geometry="any"):
    """
    Return a list of (S_value, description) tuples suggesting possible
    electron-spin quantum numbers for a d-block ion.

    geometry: one of GEOMETRY_OPTIONS. "any" returns the union of
    suggestions across all known geometries, each labeled with which
    geometry/geometries it applies to.

    Returns None if the element isn't in the d-block table (e.g. main
    group or f-block - use f_block_hint() for f-block elements instead).
    """
    d_count = d_electron_count(element, oxidation_state)
    if d_count is None or d_count < 0 or d_count > 10:
        return None

    if geometry != "any":
        table = _GEOMETRY_TABLES.get(geometry, _OCTAHEDRAL_TABLE)
        return table.get(d_count, [])

    # "any": union across geometries, merged by S value with geometry
    # names appended to the label so the user can see which geometries
    # give which answer.
    merged = {}  # S -> set of (geometry, label)
    for geom_name, table in _GEOMETRY_TABLES.items():
        for S, label in table.get(d_count, []):
            merged.setdefault(S, []).append((geom_name, label))

    out = []
    for S in sorted(merged.keys()):
        geoms = merged[S]
        geom_names = ", ".join(g for g, _ in geoms)
        # use the first label's descriptive text, prefixed with which
        # geometries it applies to
        out.append((S, f"S={S:g} - applies to: {geom_names}"))
    return out


def f_block_hint(element):
    """
    Return a clearly-flagged Hund's-rule-based hint string for a
    lanthanide/actinide element, or None if not in the table.

    THIS IS NOT A CONFIDENT SUGGESTION. Strong spin-orbit coupling in
    f-elements typically means the EPR-active effective spin differs from
    the free-ion Hund's-rule S (e.g. many actinide/lanthanide EPR signals
    reflect a Kramers doublet, effective S=1/2, rather than the free-ion
    ground term). Always verify against literature or the spectrum itself.
    """
    n_f = F_BLOCK_BASE_ELECTRONS.get(element)
    if n_f is None:
        return None

    # Hund's rule: maximize S first. Max S for n electrons in 7 f-orbitals:
    if n_f <= 7:
        S = n_f / 2
    else:
        S = (14 - n_f) / 2

    return (
        f"UNVERIFIED HINT ONLY: naive Hund's-rule ground state for "
        f"{element} with {n_f} f-electrons gives free-ion S={S:g}. "
        f"Actual EPR-active effective spin for f-block ions is frequently "
        f"different (spin-orbit coupling, crystal-field splitting, "
        f"Kramers doublets). Verify against literature or the spectrum "
        f"before using this value."
    )


def is_f_block(element):
    return element in F_BLOCK_ELEMENTS


def is_d_block(element):
    return element in D_BLOCK_GROUPS


# --------------------------------------------------------------------------- #
# Periodic table picker data: standard IUPAC grid position (row, col) for
# every element - conventional chemistry layout, not a fact that needs
# independent verification the way isotope abundances do.
# --------------------------------------------------------------------------- #
PERIODIC_TABLE_LAYOUT = {
    "H": (1, 1), "He": (1, 18),
    "Li": (2, 1), "Be": (2, 2), "B": (2, 13), "C": (2, 14), "N": (2, 15),
    "O": (2, 16), "F": (2, 17), "Ne": (2, 18),
    "Na": (3, 1), "Mg": (3, 2), "Al": (3, 13), "Si": (3, 14), "P": (3, 15),
    "S": (3, 16), "Cl": (3, 17), "Ar": (3, 18),
    "K": (4, 1), "Ca": (4, 2), "Sc": (4, 3), "Ti": (4, 4), "V": (4, 5),
    "Cr": (4, 6), "Mn": (4, 7), "Fe": (4, 8), "Co": (4, 9), "Ni": (4, 10),
    "Cu": (4, 11), "Zn": (4, 12), "Ga": (4, 13), "Ge": (4, 14), "As": (4, 15),
    "Se": (4, 16), "Br": (4, 17), "Kr": (4, 18),
    "Rb": (5, 1), "Sr": (5, 2), "Y": (5, 3), "Zr": (5, 4), "Nb": (5, 5),
    "Mo": (5, 6), "Tc": (5, 7), "Ru": (5, 8), "Rh": (5, 9), "Pd": (5, 10),
    "Ag": (5, 11), "Cd": (5, 12), "In": (5, 13), "Sn": (5, 14), "Sb": (5, 15),
    "Te": (5, 16), "I": (5, 17), "Xe": (5, 18),
    "Cs": (6, 1), "Ba": (6, 2),
    "Hf": (6, 4), "Ta": (6, 5), "W": (6, 6), "Re": (6, 7), "Os": (6, 8),
    "Ir": (6, 9), "Pt": (6, 10), "Au": (6, 11), "Hg": (6, 12), "Tl": (6, 13),
    "Pb": (6, 14), "Bi": (6, 15), "Po": (6, 16), "At": (6, 17), "Rn": (6, 18),
    "Fr": (7, 1), "Ra": (7, 2),
    "Rf": (7, 4), "Db": (7, 5), "Sg": (7, 6), "Bh": (7, 7), "Hs": (7, 8),
    "Mt": (7, 9), "Ds": (7, 10), "Rg": (7, 11), "Cn": (7, 12), "Nh": (7, 13),
    "Fl": (7, 14), "Mc": (7, 15), "Lv": (7, 16), "Ts": (7, 17), "Og": (7, 18),
    # Lanthanides (row 9)
    "La": (9, 3), "Ce": (9, 4), "Pr": (9, 5), "Nd": (9, 6), "Pm": (9, 7),
    "Sm": (9, 8), "Eu": (9, 9), "Gd": (9, 10), "Tb": (9, 11), "Dy": (9, 12),
    "Ho": (9, 13), "Er": (9, 14), "Tm": (9, 15), "Yb": (9, 16), "Lu": (9, 17),
    # Actinides (row 10)
    "Ac": (10, 3), "Th": (10, 4), "Pa": (10, 5), "U": (10, 6), "Np": (10, 7),
    "Pu": (10, 8), "Am": (10, 9), "Cm": (10, 10), "Bk": (10, 11), "Cf": (10, 12),
    "Es": (10, 13), "Fm": (10, 14), "Md": (10, 15), "No": (10, 16), "Lr": (10, 17),
}

ISOTOPE_TABLE = {
    # (mass_number, spin_I_string, abundance_percent_or_None, note)
    # Source: EasySpin's OWN internal nuclear isotope database (verbatim from
    # easyspin/private/isotopedata.txt, per the official documentation page
    # https://easyspin.org/easyspin/documentation/isotopetable.html), not a
    # separately-compiled reference - chosen specifically so labels/abundances
    # match what EasySpin itself uses internally. abundance is None for
    # elements with no specific isotope in EasySpin's database at all.
    "H": [(1, '1/2', 99.9885, ''), (2, '1', 0.0115, ''), (3, '1/2', 0.0, 'synthetic, not naturally occurring')],
    "He": [(3, '1/2', 0.000137, ''), (4, '0', 99.999863, 'I=0, no hyperfine')],
    "Li": [(6, '1', 7.59, ''), (7, '3/2', 92.41, '')],
    "Be": [(9, '3/2', 100.0, '')],
    "B": [(10, '3', 19.9, ''), (11, '3/2', 80.1, '')],
    "C": [(12, '0', 98.93, 'I=0, no hyperfine'), (13, '1/2', 1.07, ''), (14, '0', 0.0, 'I=0, no hyperfine')],
    "N": [(14, '1', 99.632, ''), (15, '1/2', 0.368, '')],
    "O": [(16, '0', 99.757, 'I=0, no hyperfine'), (17, '5/2', 0.038, ''), (18, '0', 0.205, 'I=0, no hyperfine')],
    "F": [(19, '1/2', 100.0, '')],
    "Ne": [(20, '0', 90.48, 'I=0, no hyperfine'), (21, '3/2', 0.27, ''), (22, '0', 9.25, 'I=0, no hyperfine')],
    "Na": [(22, '3', 0.0, 'synthetic, not naturally occurring'), (23, '3/2', 100.0, '')],
    "Mg": [(24, '0', 78.99, 'I=0, no hyperfine'), (25, '5/2', 10.0, ''), (26, '0', 11.01, 'I=0, no hyperfine')],
    "Al": [(27, '5/2', 100.0, '')],
    "Si": [(28, '0', 92.2297, 'I=0, no hyperfine'), (29, '1/2', 4.6832, ''), (30, '0', 3.0872, 'I=0, no hyperfine')],
    "P": [(31, '1/2', 100.0, '')],
    "S": [(32, '0', 94.93, 'I=0, no hyperfine'), (33, '3/2', 0.76, ''), (34, '0', 4.29, 'I=0, no hyperfine'), (36, '0', 0.02, 'I=0, no hyperfine')],
    "Cl": [(35, '3/2', 75.78, ''), (36, '2', 0.0, 'synthetic, not naturally occurring'), (37, '3/2', 24.22, '')],
    "Ar": [(36, '0', 0.3365, 'I=0, no hyperfine'), (38, '0', 0.0632, 'I=0, no hyperfine'), (39, '7/2', 0.0, 'synthetic, not naturally occurring'), (40, '0', 99.6003, 'I=0, no hyperfine')],
    "K": [(39, '3/2', 93.2581, ''), (40, '4', 0.0117, ''), (41, '3/2', 6.7302, '')],
    "Ca": [(40, '0', 96.941, 'I=0, no hyperfine'), (41, '7/2', 0.0, 'synthetic, not naturally occurring'), (42, '0', 0.647, 'I=0, no hyperfine'), (43, '7/2', 0.135, ''), (44, '0', 2.086, 'I=0, no hyperfine'), (46, '0', 0.004, 'I=0, no hyperfine'), (48, '0', 0.187, 'I=0, no hyperfine')],
    "Sc": [(45, '7/2', 100.0, '')],
    "Ti": [(46, '0', 8.25, 'I=0, no hyperfine'), (47, '5/2', 7.44, ''), (48, '0', 73.72, 'I=0, no hyperfine'), (49, '7/2', 5.41, ''), (50, '0', 5.18, 'I=0, no hyperfine')],
    "V": [(50, '6', 0.25, ''), (51, '7/2', 99.75, '')],
    "Cr": [(50, '0', 4.345, 'I=0, no hyperfine'), (52, '0', 83.789, 'I=0, no hyperfine'), (53, '3/2', 9.501, ''), (54, '0', 2.365, 'I=0, no hyperfine')],
    "Mn": [(53, '7/2', 0.0, 'synthetic, not naturally occurring'), (55, '5/2', 100.0, '')],
    "Fe": [(54, '0', 5.845, 'I=0, no hyperfine'), (56, '0', 91.754, 'I=0, no hyperfine'), (57, '1/2', 2.119, ''), (58, '0', 0.282, 'I=0, no hyperfine')],
    "Co": [(59, '7/2', 100.0, ''), (60, '5', 0.0, 'synthetic, not naturally occurring')],
    "Ni": [(58, '0', 68.0769, 'I=0, no hyperfine'), (60, '0', 26.2231, 'I=0, no hyperfine'), (61, '3/2', 1.1399, ''), (62, '0', 3.6345, 'I=0, no hyperfine'), (64, '0', 0.9256, 'I=0, no hyperfine')],
    "Cu": [(63, '3/2', 69.17, ''), (65, '3/2', 30.83, '')],
    "Zn": [(64, '0', 48.63, 'I=0, no hyperfine'), (66, '0', 27.9, 'I=0, no hyperfine'), (67, '5/2', 4.1, ''), (68, '0', 18.75, 'I=0, no hyperfine'), (70, '0', 0.62, 'I=0, no hyperfine')],
    "Ga": [(69, '3/2', 60.108, ''), (71, '3/2', 39.892, '')],
    "Ge": [(70, '0', 20.84, 'I=0, no hyperfine'), (72, '0', 27.54, 'I=0, no hyperfine'), (73, '9/2', 7.73, ''), (74, '0', 36.28, 'I=0, no hyperfine'), (76, '0', 7.61, 'I=0, no hyperfine')],
    "As": [(75, '3/2', 100.0, '')],
    "Se": [(74, '0', 0.89, 'I=0, no hyperfine'), (76, '0', 9.37, 'I=0, no hyperfine'), (77, '1/2', 7.63, ''), (78, '0', 23.77, 'I=0, no hyperfine'), (79, '7/2', 0.0, 'synthetic, not naturally occurring'), (80, '0', 49.61, 'I=0, no hyperfine'), (82, '0', 8.73, 'I=0, no hyperfine')],
    "Br": [(79, '3/2', 50.69, ''), (81, '3/2', 49.31, '')],
    "Kr": [(78, '0', 0.35, 'I=0, no hyperfine'), (80, '0', 2.28, 'I=0, no hyperfine'), (82, '0', 11.58, 'I=0, no hyperfine'), (83, '9/2', 11.49, ''), (84, '0', 57.0, 'I=0, no hyperfine'), (85, '9/2', 0.0, 'synthetic, not naturally occurring'), (86, '0', 17.3, 'I=0, no hyperfine')],
    "Rb": [(85, '5/2', 72.17, ''), (87, '3/2', 27.83, '')],
    "Sr": [(84, '0', 0.56, 'I=0, no hyperfine'), (86, '0', 9.86, 'I=0, no hyperfine'), (87, '9/2', 7.0, ''), (88, '0', 82.58, 'I=0, no hyperfine')],
    "Y": [(89, '1/2', 100.0, '')],
    "Zr": [(90, '0', 51.45, 'I=0, no hyperfine'), (91, '5/2', 11.22, ''), (92, '0', 17.15, 'I=0, no hyperfine'), (94, '0', 17.38, 'I=0, no hyperfine'), (96, '0', 2.8, 'I=0, no hyperfine')],
    "Nb": [(93, '9/2', 100.0, '')],
    "Mo": [(92, '0', 14.84, 'I=0, no hyperfine'), (94, '0', 9.25, 'I=0, no hyperfine'), (95, '5/2', 15.92, ''), (96, '0', 16.68, 'I=0, no hyperfine'), (97, '5/2', 9.55, ''), (98, '0', 24.13, 'I=0, no hyperfine'), (100, '0', 9.63, 'I=0, no hyperfine')],
    "Tc": [(99, '9/2', 0.0, 'synthetic, not naturally occurring')],
    "Ru": [(96, '0', 5.54, 'I=0, no hyperfine'), (98, '0', 1.87, 'I=0, no hyperfine'), (99, '5/2', 12.76, ''), (100, '0', 12.6, 'I=0, no hyperfine'), (101, '5/2', 17.06, ''), (102, '0', 31.55, 'I=0, no hyperfine'), (104, '0', 18.62, 'I=0, no hyperfine')],
    "Rh": [(103, '1/2', 100.0, '')],
    "Pd": [(102, '0', 1.02, 'I=0, no hyperfine'), (104, '0', 11.14, 'I=0, no hyperfine'), (105, '5/2', 22.33, ''), (106, '0', 27.33, 'I=0, no hyperfine'), (108, '0', 26.46, 'I=0, no hyperfine'), (110, '0', 11.72, 'I=0, no hyperfine')],
    "Ag": [(107, '1/2', 51.839, ''), (109, '1/2', 48.161, '')],
    "Cd": [(106, '0', 1.25, 'I=0, no hyperfine'), (108, '0', 0.89, 'I=0, no hyperfine'), (110, '0', 12.49, 'I=0, no hyperfine'), (111, '1/2', 12.8, ''), (112, '0', 24.13, 'I=0, no hyperfine'), (113, '1/2', 12.22, ''), (114, '0', 28.73, 'I=0, no hyperfine'), (116, '0', 7.49, 'I=0, no hyperfine')],
    "In": [(113, '9/2', 4.29, ''), (115, '9/2', 95.71, '')],
    "Sn": [(112, '0', 0.97, 'I=0, no hyperfine'), (114, '0', 0.66, 'I=0, no hyperfine'), (115, '1/2', 0.34, ''), (116, '0', 14.54, 'I=0, no hyperfine'), (117, '1/2', 7.68, ''), (118, '0', 24.22, 'I=0, no hyperfine'), (119, '1/2', 8.59, ''), (120, '0', 32.58, 'I=0, no hyperfine'), (122, '0', 4.63, 'I=0, no hyperfine'), (124, '0', 5.79, 'I=0, no hyperfine')],
    "Sb": [(121, '5/2', 57.21, ''), (123, '7/2', 42.79, ''), (125, '7/2', 0.0, 'synthetic, not naturally occurring')],
    "Te": [(120, '0', 0.09, 'I=0, no hyperfine'), (122, '0', 2.55, 'I=0, no hyperfine'), (123, '1/2', 0.89, ''), (124, '0', 4.74, 'I=0, no hyperfine'), (125, '1/2', 7.07, ''), (126, '0', 18.84, 'I=0, no hyperfine'), (128, '0', 31.74, 'I=0, no hyperfine'), (130, '0', 34.08, 'I=0, no hyperfine')],
    "I": [(127, '5/2', 100.0, ''), (129, '7/2', 0.0, 'synthetic, not naturally occurring')],
    "Xe": [(124, '0', 0.09, 'I=0, no hyperfine'), (126, '0', 0.09, 'I=0, no hyperfine'), (128, '0', 1.92, 'I=0, no hyperfine'), (129, '1/2', 26.44, ''), (130, '0', 4.08, 'I=0, no hyperfine'), (131, '3/2', 21.18, ''), (132, '0', 26.89, 'I=0, no hyperfine'), (134, '0', 10.44, 'I=0, no hyperfine'), (136, '0', 8.87, 'I=0, no hyperfine')],
    "Cs": [(133, '7/2', 100.0, ''), (134, '4', 0.0, 'synthetic, not naturally occurring'), (135, '7/2', 0.0, 'synthetic, not naturally occurring'), (137, '7/2', 0.0, 'synthetic, not naturally occurring')],
    "Ba": [(130, '0', 0.106, 'I=0, no hyperfine'), (132, '0', 0.101, 'I=0, no hyperfine'), (133, '1/2', 0.0, 'synthetic, not naturally occurring'), (134, '0', 2.417, 'I=0, no hyperfine'), (135, '3/2', 6.592, ''), (136, '0', 7.854, 'I=0, no hyperfine'), (137, '3/2', 11.232, ''), (138, '0', 71.698, 'I=0, no hyperfine')],
    "La": [(137, '7/2', 0.0, 'synthetic, not naturally occurring'), (138, '5', 0.09, ''), (139, '7/2', 99.91, '')],
    "Ce": [(136, '0', 0.185, 'I=0, no hyperfine'), (138, '0', 0.251, 'I=0, no hyperfine'), (140, '0', 88.45, 'I=0, no hyperfine'), (142, '0', 11.114, 'I=0, no hyperfine')],
    "Pr": [(141, '5/2', 100.0, '')],
    "Nd": [(142, '0', 27.2, 'I=0, no hyperfine'), (143, '7/2', 12.2, ''), (144, '0', 23.8, 'I=0, no hyperfine'), (145, '7/2', 8.3, ''), (146, '0', 17.2, 'I=0, no hyperfine'), (148, '0', 5.7, 'I=0, no hyperfine'), (150, '0', 5.6, 'I=0, no hyperfine')],
    "Pm": [(147, '7/2', 0.0, 'synthetic, not naturally occurring')],
    "Sm": [(144, '0', 3.07, 'I=0, no hyperfine'), (147, '7/2', 14.99, ''), (148, '0', 11.24, 'I=0, no hyperfine'), (149, '7/2', 13.82, ''), (150, '0', 7.38, 'I=0, no hyperfine'), (151, '5/2', 0.0, 'synthetic, not naturally occurring'), (152, '0', 26.75, 'I=0, no hyperfine'), (154, '0', 22.75, 'I=0, no hyperfine')],
    "Eu": [(151, '5/2', 47.81, ''), (152, '3', 0.0, 'synthetic, not naturally occurring'), (153, '5/2', 52.19, ''), (154, '3', 0.0, 'synthetic, not naturally occurring'), (155, '5/2', 0.0, 'synthetic, not naturally occurring')],
    "Gd": [(152, '0', 0.2, 'I=0, no hyperfine'), (154, '0', 2.18, 'I=0, no hyperfine'), (155, '3/2', 14.8, ''), (156, '0', 20.47, 'I=0, no hyperfine'), (157, '3/2', 15.65, ''), (158, '0', 24.84, 'I=0, no hyperfine'), (160, '0', 21.86, 'I=0, no hyperfine')],
    "Tb": [(157, '3/2', 0.0, 'synthetic, not naturally occurring'), (159, '3/2', 100.0, ''), (160, '3', 0.0, 'synthetic, not naturally occurring')],
    "Dy": [(156, '0', 0.06, 'I=0, no hyperfine'), (158, '0', 0.1, 'I=0, no hyperfine'), (160, '0', 2.34, 'I=0, no hyperfine'), (161, '5/2', 18.91, ''), (162, '0', 25.51, 'I=0, no hyperfine'), (163, '5/2', 24.9, ''), (164, '0', 28.18, 'I=0, no hyperfine')],
    "Ho": [(165, '7/2', 100.0, '')],
    "Er": [(162, '0', 0.14, 'I=0, no hyperfine'), (164, '0', 1.61, 'I=0, no hyperfine'), (166, '0', 33.61, 'I=0, no hyperfine'), (167, '7/2', 22.93, ''), (168, '0', 26.78, 'I=0, no hyperfine'), (170, '0', 14.93, 'I=0, no hyperfine')],
    "Tm": [(169, '1/2', 100.0, ''), (171, '1/2', 0.0, 'synthetic, not naturally occurring')],
    "Yb": [(168, '0', 0.13, 'I=0, no hyperfine'), (170, '0', 3.04, 'I=0, no hyperfine'), (171, '1/2', 14.28, ''), (172, '0', 21.83, 'I=0, no hyperfine'), (173, '5/2', 16.13, ''), (174, '0', 31.83, 'I=0, no hyperfine'), (176, '0', 12.76, 'I=0, no hyperfine')],
    "Lu": [(173, '7/2', 0.0, 'synthetic, not naturally occurring'), (174, '1', 0.0, 'synthetic, not naturally occurring'), (175, '7/2', 97.41, ''), (176, '7', 2.59, '')],
    "Hf": [(174, '0', 0.16, 'I=0, no hyperfine'), (176, '0', 5.26, 'I=0, no hyperfine'), (177, '7/2', 18.6, ''), (178, '0', 27.28, 'I=0, no hyperfine'), (179, '9/2', 13.62, ''), (180, '0', 35.08, 'I=0, no hyperfine')],
    "Ta": [(180, '9', 0.012, ''), (181, '7/2', 99.988, '')],
    "W": [(180, '0', 0.12, 'I=0, no hyperfine'), (182, '0', 26.5, 'I=0, no hyperfine'), (183, '1/2', 14.31, ''), (184, '0', 30.64, 'I=0, no hyperfine'), (186, '0', 28.43, 'I=0, no hyperfine')],
    "Re": [(185, '5/2', 37.4, ''), (187, '5/2', 62.6, '')],
    "Os": [(184, '0', 0.02, 'I=0, no hyperfine'), (186, '0', 1.59, 'I=0, no hyperfine'), (187, '1/2', 1.96, ''), (188, '0', 13.24, 'I=0, no hyperfine'), (189, '3/2', 16.15, ''), (190, '0', 26.26, 'I=0, no hyperfine'), (192, '0', 40.78, 'I=0, no hyperfine')],
    "Ir": [(191, '3/2', 37.3, ''), (193, '3/2', 62.7, '')],
    "Pt": [(190, '0', 0.014, 'I=0, no hyperfine'), (192, '0', 0.784, 'I=0, no hyperfine'), (194, '0', 32.967, 'I=0, no hyperfine'), (195, '1/2', 33.832, ''), (196, '0', 25.242, 'I=0, no hyperfine'), (198, '0', 7.163, 'I=0, no hyperfine')],
    "Au": [(197, '3/2', 100.0, '')],
    "Hg": [(196, '0', 0.15, 'I=0, no hyperfine'), (198, '0', 9.97, 'I=0, no hyperfine'), (199, '1/2', 16.87, ''), (200, '0', 23.1, 'I=0, no hyperfine'), (201, '3/2', 13.18, ''), (202, '0', 29.86, 'I=0, no hyperfine'), (204, '0', 6.87, 'I=0, no hyperfine')],
    "Tl": [(203, '1/2', 29.524, ''), (204, '2', 0.0, 'synthetic, not naturally occurring'), (205, '1/2', 70.476, '')],
    "Pb": [(204, '0', 1.4, 'I=0, no hyperfine'), (206, '0', 24.1, 'I=0, no hyperfine'), (207, '1/2', 22.1, ''), (208, '0', 52.4, 'I=0, no hyperfine')],
    "Bi": [(207, '9/2', 0.0, 'synthetic, not naturally occurring'), (209, '9/2', 100.0, '')],
    "Po": [(209, '1/2', 0.0, 'synthetic, not naturally occurring')],
    "At": [(0, 'UNVERIFIED - could not confirm', None, 'no specific isotope in EasySpin\'s own database')],
    "Rn": [(0, 'UNVERIFIED - could not confirm', None, 'no specific isotope in EasySpin\'s own database')],
    "Fr": [(0, 'UNVERIFIED - could not confirm', None, 'no specific isotope in EasySpin\'s own database')],
    "Ra": [(0, 'UNVERIFIED - could not confirm', None, 'no specific isotope in EasySpin\'s own database')],
    "Ac": [(227, '3/2', 0.0, 'synthetic, not naturally occurring')],
    "Th": [(229, '5/2', 0.0, 'synthetic, not naturally occurring'), (232, '0', 100.0, 'I=0, no hyperfine')],
    "Pa": [(231, '3/2', 100.0, 'radioactive, naturally occurring trace')],
    "U": [(234, '0', 0.0055, 'I=0, no hyperfine'), (235, '7/2', 0.72, 'radioactive, naturally occurring trace'), (238, '0', 99.2745, 'I=0, no hyperfine')],
    "Np": [(237, '5/2', 0.0, 'synthetic, not naturally occurring')],
    "Pu": [(239, '1/2', 0.0, 'synthetic, not naturally occurring')],
    "Am": [(243, '5/2', 0.0, 'synthetic, not naturally occurring')],
    "Cm": [(0, 'UNVERIFIED - could not confirm', None, 'no specific isotope in EasySpin\'s own database')],
    "Bk": [(0, 'UNVERIFIED - could not confirm', None, 'no specific isotope in EasySpin\'s own database')],
    "Cf": [(0, 'UNVERIFIED - could not confirm', None, 'no specific isotope in EasySpin\'s own database')],
    "Es": [(0, 'UNVERIFIED - could not confirm', None, 'no specific isotope in EasySpin\'s own database')],
    "Fm": [(0, 'UNVERIFIED - could not confirm', None, 'no specific isotope in EasySpin\'s own database')],
    "Md": [(0, 'UNVERIFIED - could not confirm', None, 'no specific isotope in EasySpin\'s own database')],
    "No": [(0, 'UNVERIFIED - could not confirm', None, 'no specific isotope in EasySpin\'s own database')],
    "Lr": [(0, 'UNVERIFIED - could not confirm', None, 'no specific isotope in EasySpin\'s own database')],
    "Rf": [(0, 'UNVERIFIED - could not confirm', None, 'no specific isotope in EasySpin\'s own database')],
    "Db": [(0, 'UNVERIFIED - could not confirm', None, 'no specific isotope in EasySpin\'s own database')],
    "Sg": [(0, 'UNVERIFIED - could not confirm', None, 'no specific isotope in EasySpin\'s own database')],
    "Bh": [(0, 'UNVERIFIED - could not confirm', None, 'no specific isotope in EasySpin\'s own database')],
    "Hs": [(0, 'UNVERIFIED - could not confirm', None, 'no specific isotope in EasySpin\'s own database')],
    "Mt": [(0, 'UNVERIFIED - could not confirm', None, 'no specific isotope in EasySpin\'s own database')],
    "Ds": [(0, 'UNVERIFIED - could not confirm', None, 'no specific isotope in EasySpin\'s own database')],
    "Rg": [(0, 'UNVERIFIED - could not confirm', None, 'no specific isotope in EasySpin\'s own database')],
    "Cn": [(0, 'UNVERIFIED - could not confirm', None, 'no specific isotope in EasySpin\'s own database')],
    "Nh": [(0, 'UNVERIFIED - could not confirm', None, 'no specific isotope in EasySpin\'s own database')],
    "Fl": [(0, 'UNVERIFIED - could not confirm', None, 'no specific isotope in EasySpin\'s own database')],
    "Mc": [(0, 'UNVERIFIED - could not confirm', None, 'no specific isotope in EasySpin\'s own database')],
    "Lv": [(0, 'UNVERIFIED - could not confirm', None, 'no specific isotope in EasySpin\'s own database')],
    "Ts": [(0, 'UNVERIFIED - could not confirm', None, 'no specific isotope in EasySpin\'s own database')],
    "Og": [(0, 'UNVERIFIED - could not confirm', None, 'no specific isotope in EasySpin\'s own database')],
}


def get_isotope_table(element):
    """Return the isotope list for `element`, or None if not in the table."""
    return ISOTOPE_TABLE.get(element)


# --------------------------------------------------------------------------- #
# Same-element auto-recognition for the Nuclei list.
#
# Rather than a separate "Isotope Mixture" panel, this lets the plain
# Nuclei list itself recognize when two or more rows are different
# isotopes of the SAME element (e.g. 47Ti and 49Ti both entered as
# separate rows) and automatically:
#   1. Fill in each row's weight % from its natural abundance (only for
#      rows still at the untouched default "100", so a deliberately
#      customized weight is never silently overwritten)
#   2. Compute the implied "background" weight (the I=0 majority
#      isotopes of that element, which contribute no hyperfine at all)
#      as an extra component to add to the simulation
#
# This only activates for a given element once 2+ rows of it are present
# AND every one of those rows is a recognized isotope with a known
# numeric abundance - if any entered isotope for that element is
# unrecognized, this is skipped entirely for that element (no partial/
# guessed background) rather than risk a wrong number silently.
# --------------------------------------------------------------------------- #
_ISOTOPE_LABEL_RE = re.compile(r'^(\d+)([A-Za-z]+)$')


def parse_isotope_label(label):
    """Split '47Ti' into (mass_number=47, element='Ti'). Returns (None,
    None) if the label doesn't match the expected <mass><symbol> form."""
    m = _ISOTOPE_LABEL_RE.match(label.strip())
    if not m:
        return None, None
    return int(m.group(1)), m.group(2)


def lookup_isotope_spin(element, mass_number):
    """Return the spin-I string (e.g. '0', '1/2', '5/2') for one specific
    isotope, or None if that element/mass combination isn't in the table."""
    isotopes = ISOTOPE_TABLE.get(element)
    if not isotopes:
        return None
    for mass, spin, abund, note in isotopes:
        if mass == mass_number:
            return spin
    return None


def is_zero_spin_isotope(iso_label):
    """True if `iso_label` (e.g. '48Ti') parses to a known isotope with
    I=0 - used at script-generation time to omit Sys.Nucs/Sys.A entirely
    for such a nucleus, since hyperfine coupling has no effect on a
    spinless nucleus (see NucleusRow's matching UI behavior, which hides
    the A field(s) for the same isotopes)."""
    mass, element = parse_isotope_label(iso_label)
    if mass is None:
        return False
    return lookup_isotope_spin(element, mass) == '0'


def lookup_isotope_abundance(element, mass_number):
    """Return the abundance percent (float) for one specific isotope, or
    None if that element/mass combination isn't in the table or has no
    known numeric abundance (e.g. a synthetic isotope)."""
    isotopes = ISOTOPE_TABLE.get(element)
    if not isotopes:
        return None
    for mass, spin, abund, note in isotopes:
        if mass == mass_number:
            return abund
    return None


def detect_same_element_groups(nuclei):
    """
    nuclei: list of (isotope_label, n, a_vals, weight_str) as returned by
    NucleusList.collect().

    Returns (updated_nuclei, backgrounds):
      updated_nuclei - same shape as `nuclei`, but any row still at the
        default weight "100" within a detected group gets auto-filled
        with its natural abundance percentage.
      backgrounds - {element: (background_fraction_0_to_1,
        [missing_active_isotope_labels])} for every element with 2+
        recognized rows and a computable (non-negative) background.
        missing_active_isotope_labels lists that element's OTHER
        naturally-active isotopes (I != 0, known abundance) that weren't
        entered as rows - included so the caller can flag this as an
        approximation if it's non-empty.
    """
    groups = {}
    for i, (iso, n, a, weight) in enumerate(nuclei):
        mass, element = parse_isotope_label(iso)
        if element:
            groups.setdefault(element, []).append(i)

    updated = list(nuclei)
    backgrounds = {}

    for element, indices in groups.items():
        if len(indices) < 2:
            continue
        isotopes_known = ISOTOPE_TABLE.get(element)
        if not isotopes_known:
            continue

        entered_masses = []
        all_known = True
        for i in indices:
            iso, n, a, weight = updated[i]
            mass, _ = parse_isotope_label(iso)
            entered_masses.append(mass)
            if lookup_isotope_abundance(element, mass) is None:
                all_known = False
        if not all_known:
            continue

        for i in indices:
            iso, n, a, weight = updated[i]
            mass, _ = parse_isotope_label(iso)
            if weight.strip() == "100":
                abund = lookup_isotope_abundance(element, mass)
                updated[i] = (iso, n, a, f"{abund:.6g}")

        total_frac = sum(float(updated[i][3]) / 100.0 for i in indices)
        background = 1.0 - total_frac

        active_isotopes = [
            f"{mass}{element}" for mass, spin, abund, note in isotopes_known
            if spin != "0" and abund is not None and mass not in entered_masses
        ]

        if background > 1e-6:
            backgrounds[element] = (background, active_isotopes)

    return updated, backgrounds


##############################################################################
# dsc_parser.py content (Bruker .DSC / data file parsing, fit-export CSV parsing, residual data-file writing)
##############################################################################

import os
import re


def find_data_companion(dsc_path):
    """
    Given a specific .DSC file path, find its companion data file (same
    basename, no extension, in the same folder). Returns the data path, or
    None if not found.
    """
    folder = os.path.dirname(dsc_path)
    stem = os.path.splitext(os.path.basename(dsc_path))[0]
    data_path = os.path.join(folder, stem)
    if os.path.isfile(data_path):
        return data_path
    return None


def find_dataset_pair(folder):
    """
    Look in `folder` for a .DSC file and its companion data file (same
    basename, no extension). Returns (dsc_path, data_path) or (None, None)
    if not found.

    NOTE: this scans for the FIRST .DSC file in the folder - if a folder
    has multiple datasets, prefer letting the user pick the exact .DSC
    file (see DatasetPanel.browse in widgets_common.py, which uses
    askopenfilename + find_data_companion instead of this function).
    """
    dsc_path = None
    for fname in os.listdir(folder):
        if fname.lower().endswith(".dsc"):
            dsc_path = os.path.join(folder, fname)
            break
    if dsc_path is None:
        return None, None

    return dsc_path, find_data_companion(dsc_path)


def parse_dsc(dsc_path):
    """
    Parse a .DSC file for the fields relevant to Exp struct auto-fill.

    Returns a dict with any of the following keys that could be found:
        mwFreq_GHz, center_G, sweep_G, harmonic, modamp_G, npoints, title
    Missing values are simply absent from the dict - caller should fall
    back to defaults / leave fields blank for the user to fill in.
    """
    result = {}

    with open(dsc_path, "r", errors="ignore") as f:
        text = f.read()

    # Title
    m = re.search(r"^TITL\s+'(.*)'", text, re.MULTILINE)
    if m:
        result["title"] = m.group(1)

    # Number of points
    m = re.search(r"^XPTS\s+(\d+)", text, re.MULTILINE)
    if m:
        result["npoints"] = int(m.group(1))

    # Microwave frequency - MWFQ is in Hz in the standard parameter layer
    m = re.search(r"^MWFQ\s+([\d.eE+-]+)", text, re.MULTILINE)
    if m:
        freq_hz = float(m.group(1))
        result["mwFreq_GHz"] = freq_hz / 1e9

    # Alternative / higher-precision frequency from device-specific layer
    m = re.search(r"FrequencyMon\s+([\d.eE+-]+)\s*GHz", text)
    if m:
        result["mwFreq_GHz"] = float(m.group(1))

    # Center field and sweep width, standard parameter layer (in Gauss,
    # though A1CT is sometimes recorded in different units - the
    # device-specific CenterField/SweepWidth lines are more reliable and
    # explicitly carry a unit, so prefer those if present).
    m = re.search(r"CenterField\s+([\d.eE+-]+)\s*G", text)
    if m:
        result["center_G"] = float(m.group(1))

    m = re.search(r"SweepWidth\s+([\d.eE+-]+)\s*G", text)
    if m:
        result["sweep_G"] = float(m.group(1))

    # Harmonic
    m = re.search(r"^Harmonic\s+(\d+)", text, re.MULTILINE)
    if m:
        result["harmonic"] = int(m.group(1))

    # Modulation amplitude
    m = re.search(r"ModAmp\s+([\d.eE+-]+)\s*G", text)
    if m:
        result["modamp_G"] = float(m.group(1))

    return result


def load_experimental_trace(data_path):
    """
    Load the plain-text data file (index, Field [G], Intensity columns).

    Returns (field_mT, intensity) as two lists of floats. Raises ValueError
    if the file doesn't look like the expected format.
    """
    field_G = []
    intensity = []

    with open(data_path, "r", errors="ignore") as f:
        lines = f.readlines()

    data_started = False
    for line in lines:
        parts = line.split()
        if len(parts) != 3:
            continue
        try:
            idx = int(parts[0])
            field = float(parts[1])
            inten = float(parts[2])
        except ValueError:
            continue  # header row or malformed line
        field_G.append(field)
        intensity.append(inten)
        data_started = True

    if not data_started:
        raise ValueError(
            f"Could not parse any numeric data rows from {data_path}. "
            "Expected columns: index, Field [G], Intensity."
        )

    field_mT = [g / 10.0 for g in field_G]
    return field_mT, intensity


def compound_id_from_filename(path):
    """Derive a suggested compound ID from a dataset file/folder name."""
    stem = os.path.splitext(os.path.basename(path))[0]
    return stem


def parse_fit_export_csv(csv_path):
    """
    Parse a 'Field_mT,Experimental,Fit' CSV - the exact format this
    tool's own esfit Part 2 export produces (see build_esfit_block's
    post_lines). Returns (field_mT, experimental, fit) as three lists of
    floats. Tolerant of a header row; raises ValueError if no numeric
    data rows are found.
    """
    field_mT, experimental, fit = [], [], []
    with open(csv_path, "r", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = [p.strip() for p in line.split(",")]
            if len(parts) < 3:
                continue
            try:
                fld = float(parts[0])
                exp = float(parts[1])
                ft = float(parts[2])
            except ValueError:
                continue  # header row
            field_mT.append(fld)
            experimental.append(exp)
            fit.append(ft)
    if not field_mT:
        raise ValueError(
            f"Could not parse any numeric data rows from {csv_path}. "
            "Expected columns: Field_mT, Experimental, Fit (comma-separated)."
        )
    return field_mT, experimental, fit


def write_residual_data_file(path, field_mT, residual):
    """
    Write (field_mT, residual) in the SAME plain-text 3-column format
    (index, Field [G], Intensity) that load_experimental_trace() expects,
    so a residual can be fed back in as a new "experimental" dataset for
    fitting a minor component separately.
    """
    with open(path, "w") as f:
        for i, (fld_mT, val) in enumerate(zip(field_mT, residual), start=1):
            field_G = fld_mT * 10.0
            f.write(f"{i}\t{field_G:.6f}\t{val:.10g}\n")


##############################################################################
# widgets_common.py content (shared widgets, hardened tooltips, plot preview, presets, multi-component esfit helpers with per-nucleus Vary A/weight and SI summary generation, periodic table picker, line broadening section, residual builder)
##############################################################################

import os
import json
import time
import datetime

import customtkinter as ctk
from tkinter import filedialog, messagebox
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg



# --------------------------------------------------------------------------- #
# Preview plot with Z/F/H hotkeys (zoom-rectangle / full view / height-fit)
# --------------------------------------------------------------------------- #
class PlotPreview(ctk.CTkFrame):
    """
    Experimental-trace preview plot with keyboard shortcuts:
      Z - click-and-drag a rectangle to zoom into that region
      F - reset to full view (both axes)
      H - autoscale the y-axis (height) to whatever's visible on x, without
          changing the current x-zoom

    Click on the plot once to give it keyboard focus before using hotkeys.
    """

    def __init__(self, master, height=300, width=900, on_crop_change=None, **kwargs):
        super().__init__(master, **kwargs)
        self.field_mT = None
        self.intensity = None
        self._full_xlim = None
        self._full_ylim = None
        self._zoom_mode = False
        self._zoom_start_px = None
        self._zoom_rect_id = None
        # v5.6: optional callback fired whenever the crop range actually
        # changes (Crop button / Enter in the boxes / Full view). Lets the
        # owning tab sync its own Exp.Range fields (rangelo/rangehi for
        # pepper, center/sweep for garlic) to what's on-screen, instead of
        # the plot crop being a purely cosmetic view-zoom that Exp.Range
        # never finds out about.
        self._on_crop_change = on_crop_change

        # The actual on-screen size is driven by the matplotlib Figure's
        # size in inches * dpi, NOT by a "height=" kwarg passed to this
        # frame - without an explicitly larger Figure, the frame just
        # shrinks back down to fit a small canvas regardless of any
        # frame-level height setting.
        dpi = 100
        self.figure = Figure(figsize=(width / dpi, height / dpi), dpi=dpi)
        self.ax = self.figure.add_subplot(111)
        self._show_placeholder()

        self.canvas = FigureCanvasTkAgg(self.figure, master=self)
        self.tk_widget = self.canvas.get_tk_widget()
        self.tk_widget.pack(fill="both", expand=True)
        self.canvas.draw()

        hint = ctk.CTkLabel(
            self, text="Click plot, then: Z = zoom-drag, F = full view, H = fit height",
            font=ctk.CTkFont(size=10), text_color="gray50",
        )
        hint.pack(anchor="w", pady=(2, 0))

        # v5.5: numeric field-range crop, added directly (no mouse-drag
        # needed each time) - re-dragging a zoom rectangle every time the
        # plot is worth inspecting was the specific complaint. Typing
        # exact bounds is faster and more precise, and the chosen crop is
        # remembered (self._crop_xlim) and re-applied automatically the
        # next time plot_experimental() redraws (e.g. loading a new
        # dataset), so it doesn't need re-entering constantly either.
        crop_row = ctk.CTkFrame(self, fg_color="transparent")
        crop_row.pack(anchor="w", pady=(4, 0), fill="x")
        ctk.CTkLabel(crop_row, text="Crop field range (mT):", font=ctk.CTkFont(size=11),
                     text_color="gray70").pack(side="left", padx=(0, 6))
        self.crop_lo_entry = ctk.CTkEntry(crop_row, width=70, placeholder_text="min")
        self.crop_lo_entry.pack(side="left", padx=(0, 4))
        self.crop_hi_entry = ctk.CTkEntry(crop_row, width=70, placeholder_text="max")
        self.crop_hi_entry.pack(side="left", padx=(0, 8))
        ctk.CTkButton(crop_row, text="Crop", width=50,
                      command=self._apply_crop_from_entries).pack(side="left", padx=(0, 6))
        ctk.CTkButton(crop_row, text="Full view", width=70,
                      command=self._clear_crop).pack(side="left")
        self.crop_lo_entry.bind("<Return>", lambda e: self._apply_crop_from_entries())
        self.crop_hi_entry.bind("<Return>", lambda e: self._apply_crop_from_entries())
        self._crop_xlim = None

        self.tk_widget.bind("<Button-1>", self._on_press)
        self.tk_widget.bind("<B1-Motion>", self._on_drag)
        self.tk_widget.bind("<ButtonRelease-1>", self._on_release)
        self.tk_widget.bind("<Key>", self._on_key)
        self.tk_widget.bind("<Enter>", lambda e: self.tk_widget.focus_set())

    def _show_placeholder(self):
        self.ax.clear()
        self.ax.set_xlabel("Field (mT)")
        self.ax.set_ylabel("Intensity")
        self.ax.text(0.5, 0.5, "No dataset loaded", ha="center", va="center",
                     transform=self.ax.transAxes, color="gray")

    def plot_experimental(self, field_mT, intensity, title="Experimental trace"):
        self.field_mT = field_mT
        self.intensity = intensity
        self.ax.clear()
        self.ax.plot(field_mT, intensity, linewidth=0.8, color="#4FC3F7")
        self.ax.set_xlabel("Field (mT)")
        self.ax.set_ylabel("Intensity")
        self.ax.set_title(title, fontsize=10)
        self.figure.tight_layout()
        self._full_xlim = self.ax.get_xlim()
        self._full_ylim = self.ax.get_ylim()
        if self._crop_xlim is not None:
            # Re-apply whatever crop range was last typed in, so reloading
            # a dataset (or any other redraw) doesn't silently discard it.
            self.ax.set_xlim(self._crop_xlim)
            self._fit_height()
        self.canvas.draw()

    def _apply_crop_from_entries(self):
        lo_str = self.crop_lo_entry.get().strip()
        hi_str = self.crop_hi_entry.get().strip()
        if not lo_str or not hi_str:
            return
        try:
            lo, hi = float(lo_str), float(hi_str)
        except ValueError:
            return
        if lo >= hi:
            return
        self._crop_xlim = (lo, hi)
        self.ax.set_xlim(lo, hi)
        self._fit_height()
        self.canvas.draw()
        if self._on_crop_change is not None:
            self._on_crop_change(lo, hi)

    def _clear_crop(self):
        self._crop_xlim = None
        self.crop_lo_entry.delete(0, "end")
        self.crop_hi_entry.delete(0, "end")
        self._reset_full_view()
        if self._on_crop_change is not None:
            # None, None tells the owning tab "crop was cleared, don't
            # touch Exp.Range" - full-view is just the plot's own zoom-out,
            # it shouldn't silently blow away whatever the user has typed
            # into Exp.Range/center+sweep.
            self._on_crop_change(None, None)

    # ------------------------------------------------------------------ #
    # Hotkeys
    # ------------------------------------------------------------------ #
    def _on_key(self, event):
        key = (event.keysym or "").lower()
        if key == "z":
            self._zoom_mode = True
        elif key == "f":
            self._reset_full_view()
        elif key == "h":
            self._fit_height()

    def _reset_full_view(self):
        if self._full_xlim is None:
            return
        self.ax.set_xlim(self._full_xlim)
        self.ax.set_ylim(self._full_ylim)
        self.canvas.draw()

    def _fit_height(self):
        if self.field_mT is None:
            return
        xlim = self.ax.get_xlim()
        visible_y = [
            y for x, y in zip(self.field_mT, self.intensity)
            if xlim[0] <= x <= xlim[1]
        ]
        if not visible_y:
            return
        ymin, ymax = min(visible_y), max(visible_y)
        pad = (ymax - ymin) * 0.05 or 0.01
        self.ax.set_ylim(ymin - pad, ymax + pad)
        self.canvas.draw()

    # ------------------------------------------------------------------ #
    # Zoom-drag (rubber-band rectangle)
    # ------------------------------------------------------------------ #
    def _on_press(self, event):
        self.tk_widget.focus_set()
        if self._zoom_mode:
            self._zoom_start_px = (event.x, event.y)

    def _on_drag(self, event):
        if not self._zoom_mode or self._zoom_start_px is None:
            return
        if self._zoom_rect_id:
            self.tk_widget.delete(self._zoom_rect_id)
        x0, y0 = self._zoom_start_px
        self._zoom_rect_id = self.tk_widget.create_rectangle(
            x0, y0, event.x, event.y, outline="#4FC3F7"
        )

    def _on_release(self, event):
        if not self._zoom_mode or self._zoom_start_px is None:
            return
        x0_px, y0_px = self._zoom_start_px
        x1_px, y1_px = event.x, event.y

        if self._zoom_rect_id:
            self.tk_widget.delete(self._zoom_rect_id)
            self._zoom_rect_id = None
        self._zoom_start_px = None
        self._zoom_mode = False

        if abs(x1_px - x0_px) < 4 or abs(y1_px - y0_px) < 4:
            return  # too small a drag - treat as a click, not a zoom

        inv = self.ax.transData.inverted()
        # Tk canvas y is measured top-down; matplotlib pixel coords are
        # bottom-up, so flip y before converting to data coordinates.
        height = self.tk_widget.winfo_height()
        (xd0, yd0) = inv.transform((x0_px, height - y0_px))
        (xd1, yd1) = inv.transform((x1_px, height - y1_px))

        self.ax.set_xlim(min(xd0, xd1), max(xd0, xd1))
        self.ax.set_ylim(min(yd0, yd1), max(yd0, yd1))
        self.canvas.draw()


# --------------------------------------------------------------------------- #
# Tooltips (hover info icon)
# --------------------------------------------------------------------------- #
class Tooltip:
    """Small floating popup shown on hover over `widget`.

    Hiding is driven by a single application-wide polling loop rather than
    per-widget <Leave> events - Tkinter's Leave event doesn't always fire
    reliably on composite widgets (like CTkLabel, which wraps internal
    children). Note: customtkinter's bind_all() is deliberately disabled
    (raises AttributeError), so the global check uses a recurring
    `.after()` timer on the root window instead of a <Motion> binding.

    Deliberately does NOT use "-topmost" - that flag pins the popup above
    EVERY window system-wide, including completely different applications
    (e.g. MATLAB), and since hiding only checks raw mouse coordinates (not
    which app has OS focus), switching to another app without moving the
    mouse could leave the tooltip floating on top of it indefinitely. The
    poll also checks whether this app still has focus at all, and hides
    immediately if not.
    """

    _active = []          # all Tooltip instances currently showing
    _watcher_installed = False
    _root = None

    def __init__(self, widget, text):
        self.widget = widget
        self.text = text
        self.tipwindow = None
        self._shown_at = None
        widget.bind("<Enter>", self.show)
        widget.bind("<Leave>", self.hide)
        Tooltip._ensure_watcher(widget)

    @classmethod
    def _ensure_watcher(cls, any_widget):
        if cls._watcher_installed:
            return
        cls._watcher_installed = True
        cls._root = any_widget.winfo_toplevel()
        # Safety net #1: hide every active tooltip on ANY click anywhere in
        # the app. Position-tracking can miss edge cases (composite-widget
        # event quirks, a widget resized/repositioned after the tooltip
        # opened, etc.) - a click is a near-universal "I've moved on" signal
        # that doesn't depend on getting hover geometry exactly right.
        try:
            cls._root.bind_all("<Button-1>", cls._hide_all, add="+")
        except Exception:
            pass  # if unavailable, the other two safety nets still apply
        cls._tick()

    @classmethod
    def _hide_all(cls, event=None):
        for tt in list(cls._active):
            tt.hide()

    @classmethod
    def _tick(cls):
        # Safety net #2: the reschedule at the end is wrapped in
        # try/finally so a single unexpected exception anywhere above
        # can't silently kill the polling loop forever, which would leave
        # any tooltip open at that moment stuck for good (no more checks
        # would ever run again).
        try:
            app_focused = True
            try:
                focused_widget = cls._root.focus_displayof()
                app_focused = focused_widget is not None
            except Exception:
                app_focused = True  # fail open - fall back to position check

            for tt in list(cls._active):
                if not app_focused:
                    tt.hide()
                else:
                    tt._check_and_maybe_hide()
        finally:
            if cls._root:
                cls._root.after(100, cls._tick)

    def show(self, event=None):
        if self.tipwindow or not self.text:
            return
        x = self.widget.winfo_rootx() + 16
        y = self.widget.winfo_rooty() + self.widget.winfo_height() + 4
        self.tipwindow = tw = ctk.CTkToplevel(self.widget)
        tw.wm_overrideredirect(True)
        # Deliberately no "-topmost" - see class docstring.
        tw.wm_geometry(f"+{x}+{y}")
        label = ctk.CTkLabel(
            tw, text=self.text, fg_color="#2B2B2B", text_color="#F0F0F0",
            corner_radius=6, wraplength=300, justify="left",
            font=ctk.CTkFont(size=11), padx=10, pady=8,
        )
        label.pack()
        if self not in Tooltip._active:
            Tooltip._active.append(self)
        self._shown_at = time.monotonic()

    def _check_and_maybe_hide(self):
        if not self.tipwindow:
            return
        # Safety net #3: a hard maximum lifetime regardless of position
        # tracking. A tooltip has no legitimate reason to stay open this
        # long - if hover/click tracking somehow missed it, this is the
        # last line of defense.
        if self._shown_at is not None and time.monotonic() - self._shown_at > 8:
            self.hide()
            return
        try:
            px = self.widget.winfo_pointerx()
            py = self.widget.winfo_pointery()
            wx0 = self.widget.winfo_rootx()
            wy0 = self.widget.winfo_rooty()
            wx1 = wx0 + self.widget.winfo_width()
            wy1 = wy0 + self.widget.winfo_height()
            if not (wx0 <= px <= wx1 and wy0 <= py <= wy1):
                self.hide()
        except Exception:
            self.hide()

    def hide(self, event=None):
        if self in Tooltip._active:
            Tooltip._active.remove(self)
        if self.tipwindow:
            try:
                self.tipwindow.destroy()
            except Exception:
                pass
            self.tipwindow = None
        self._shown_at = None


def info_icon(parent, tooltip_text):
    """A small 'i' badge that shows `tooltip_text` on hover. Pack it yourself."""
    icon = ctk.CTkLabel(
        parent, text="\u24d8", width=16, text_color="#6FA8DC",
        font=ctk.CTkFont(size=13), cursor="hand2",
    )
    Tooltip(icon, tooltip_text)
    return icon


# --------------------------------------------------------------------------- #
# Small layout helpers
# --------------------------------------------------------------------------- #
def labeled_entry(parent, label, default="", width=140, placeholder="", tooltip=None):
    frame = ctk.CTkFrame(parent, fg_color="transparent")
    frame.pack(anchor="w", fill="x", pady=3)
    label_frame = ctk.CTkFrame(frame, fg_color="transparent", width=190, height=24)
    label_frame.pack(side="left")
    label_frame.pack_propagate(False)
    ctk.CTkLabel(label_frame, text=label, anchor="w").pack(side="left")
    if tooltip:
        info_icon(label_frame, tooltip).pack(side="left", padx=(4, 0))
    entry = ctk.CTkEntry(frame, width=width, placeholder_text=placeholder)
    if default:
        entry.insert(0, default)
    entry.pack(side="left", padx=6)
    return entry, frame


def section_title(parent, text):
    ctk.CTkLabel(parent, text=text, font=ctk.CTkFont(size=15, weight="bold")).pack(
        anchor="w", pady=(14, 4)
    )


def info_note(parent, text):
    ctk.CTkLabel(
        parent, text=text, font=ctk.CTkFont(size=11), text_color="gray60",
        wraplength=380, justify="left", anchor="w",
    ).pack(anchor="w", pady=(0, 6))


# --------------------------------------------------------------------------- #
# Periodic table picker: click an element, see its isotopes (mass number,
# nuclear spin I, natural abundance), click one to fill in an isotope
# entry field. Data compiled from general chemistry knowledge for
# commonly EPR-relevant elements - see chem_helpers.ISOTOPE_TABLE for the
# accuracy caveat. Elements outside that table show a clear "no data"
# message rather than a guess.
# --------------------------------------------------------------------------- #
class PeriodicTablePicker(ctk.CTkToplevel):
    def __init__(self, master, on_select):
        super().__init__(master)
        self.title("Periodic Table - pick an isotope")
        self.geometry("980x560")
        self.on_select = on_select

        self.grid_frame = ctk.CTkFrame(self, fg_color="transparent")
        self.grid_frame.pack(fill="both", expand=True, padx=10, pady=10)

        self.isotope_frame = ctk.CTkFrame(self, fg_color="transparent")
        # isotope_frame is only packed when showing an element's isotopes

        note = ctk.CTkLabel(
            self, text="Isotope data from a verified IUPAC/CIAAW/NIST reference "
                       "table (spin I and natural abundance) - see chem_helpers.py "
                       "for sources. Still worth a spot-check for critical work.",
            font=ctk.CTkFont(size=11), text_color="gray60",
            wraplength=940, justify="left",
        )
        note.pack(anchor="w", padx=14, pady=(0, 8))

        self._build_grid()

        # Without this, some window managers open a new Toplevel BEHIND its
        # parent instead of in front - transient() ties it to the parent
        # window, lift() raises it, focus_force() gives it keyboard focus,
        # and grab_set() makes it modal so it can't get lost behind the
        # main window again while it's open.
        self.transient(master.winfo_toplevel())
        self.lift()
        self.focus_force()
        self.grab_set()

        # v6.4: FIXED (real bug, reported directly): selecting an isotope
        # would call self.destroy() without first releasing the grab_set()
        # above. On some platforms that leaves the modal grab dangling -
        # the popup can look like it "won't close" (or the app stops
        # responding to clicks) even though destroy() technically ran,
        # since the OS/window-manager-level input grab was never handed
        # back. Also wired the window's own close (X) button through the
        # same release-then-destroy path, since closing that way had the
        # identical problem.
        self.protocol("WM_DELETE_WINDOW", self._close)

    def _close(self):
        self.grab_release()
        self.destroy()

    def _build_grid(self):
        for symbol, (row, col) in PERIODIC_TABLE_LAYOUT.items():
            has_data = get_isotope_table(symbol) is not None
            color = "#2E7D32" if has_data else "#3A3A3A"
            hover = "#245F27" if has_data else "#333333"
            btn = ctk.CTkButton(
                self.grid_frame, text=symbol, width=48, height=36,
                fg_color=color, hover_color=hover,
                font=ctk.CTkFont(size=12),
                command=lambda s=symbol: self._show_isotopes(s),
            )
            btn.grid(row=row, column=col, padx=1, pady=1)

    def _show_isotopes(self, element):
        self.grid_frame.pack_forget()
        for w in self.isotope_frame.winfo_children():
            w.destroy()

        back_btn = ctk.CTkButton(self.isotope_frame, text="< Back to periodic table",
                                  width=180, command=self._show_grid)
        back_btn.pack(anchor="w", pady=(0, 10))

        isotopes = get_isotope_table(element)
        if isotopes is None:
            ctk.CTkLabel(
                self.isotope_frame,
                text=f"No isotope data available for {element} in this tool yet.\n"
                     f"Check a reference table (e.g. NIST isotope data) for its "
                     f"nuclear spin and abundance values.",
                text_color="gray60", justify="left",
            ).pack(anchor="w", pady=10)
        else:
            ctk.CTkLabel(self.isotope_frame, text=f"{element} isotopes",
                         font=ctk.CTkFont(size=16, weight="bold")).pack(anchor="w", pady=(0, 8))
            header = ctk.CTkFrame(self.isotope_frame, fg_color="transparent")
            header.pack(anchor="w", fill="x")
            for text, w in [("Isotope", 100), ("Spin I", 80), ("Abundance", 100), ("Note", 340)]:
                ctk.CTkLabel(header, text=text, width=w, anchor="w",
                             font=ctk.CTkFont(size=12, weight="bold")).pack(side="left", padx=4)

            for mass, spin, abundance, note in isotopes:
                if mass == 0:
                    # Placeholder entry (EasySpin's own database has no
                    # specific isotope at all for this element) - show an
                    # explanation instead of a selectable (fake) isotope.
                    ctk.CTkLabel(
                        self.isotope_frame,
                        text=f"EasySpin's own database has no specific isotope "
                             f"for {element} - it can't be used in a Sys.Nucs "
                             f"field as-is.",
                        text_color="gray60", justify="left",
                    ).pack(anchor="w", pady=10)
                    continue
                label = f"{mass}{element}"
                row = ctk.CTkFrame(self.isotope_frame, fg_color="transparent")
                row.pack(anchor="w", fill="x", pady=2)
                ctk.CTkButton(row, text=label, width=100,
                              command=lambda l=label: self._select(l)).pack(side="left", padx=4)
                ctk.CTkLabel(row, text=f"I = {spin}", width=80, anchor="w").pack(side="left", padx=4)
                abund_text = f"{abundance:.3g}%" if abundance is not None else "n/a"
                ctk.CTkLabel(row, text=abund_text, width=100, anchor="w").pack(side="left", padx=4)
                ctk.CTkLabel(row, text=note, width=340, anchor="w",
                             text_color="gray60", font=ctk.CTkFont(size=11)).pack(side="left", padx=4)

        self.isotope_frame.pack(fill="both", expand=True, padx=10, pady=10)

    def _show_grid(self):
        self.isotope_frame.pack_forget()
        self.grid_frame.pack(fill="both", expand=True, padx=10, pady=10)

    def _select(self, isotope_label):
        if self.on_select:
            self.on_select(isotope_label)
        self._close()


# --------------------------------------------------------------------------- #
# One row describing a nucleus: isotope, count, hyperfine value(s) in mT
# --------------------------------------------------------------------------- #
class NucleusRow(ctk.CTkFrame):
    """
    mode="iso"    -> single scalar A value (mT) [garlic]
    mode="tensor" -> three principal values Axx,Ayy,Azz (mT) [pepper/chili]

    Each row also carries a weight % field: every nucleus row is treated
    as its OWN separate weighted Sys component (Sys = {Sys1, Sys2, ...}
    per EasySpin's documented mixture convention), not combined onto one
    shared molecule. Defaults to 100% for the common single-row case.
    "Vary" controls used to live here but were moved to the Fit (esfit)
    section, since which parameters to fit is a fitting-workflow concern,
    not a spin-system-definition one.
    """

    def __init__(self, master, mode="iso", on_remove=None, on_change=None,
                 balance_callback=None, **kwargs):
        super().__init__(master, fg_color="transparent", **kwargs)
        self.mode = mode
        self.on_remove = on_remove
        self.on_change = on_change
        self.balance_callback = balance_callback

        # v5.2: switched from grid to pack for the control row, specifically
        # so hiding/showing the A field(s) for I=0 vs I!=0 isotopes (see
        # _update_a_visibility) can use pack's before=/after= to reinsert a
        # widget at the right spot - grid would leave a dead column gap
        # instead of the rest of the row sliding over to fill the space.
        # Also single-line by construction now: everything below is one
        # pack(side="left") chain in one frame, so the row can never
        # "bend" into a second line the way the old two-row grid layout
        # did (that layout existed only to keep the "->100%" button
        # reachable - fixed here instead by shrinking widths).
        fields_row = ctk.CTkFrame(self, fg_color="transparent")
        fields_row.grid(row=0, column=0, sticky="w")

        self.isotope_entry = ctk.CTkEntry(fields_row, placeholder_text="e.g. 14N", width=58)
        self.isotope_entry.pack(side="left", padx=(0, 3))
        self.isotope_entry.bind("<Button-1>", self._open_periodic_table)
        self.isotope_entry.bind("<KeyRelease>", self._fire_on_change)

        self.n_entry = ctk.CTkEntry(fields_row, placeholder_text="n", width=30)
        self.n_entry.insert(0, "1")
        self.n_entry.pack(side="left", padx=3)
        self.n_entry.bind("<KeyRelease>", self._fire_on_change)

        if mode == "iso":
            self.a_entries = [ctk.CTkEntry(fields_row, placeholder_text="A (mT)", width=60)]
            self.a_entries[0].insert(0, "1")
        else:
            self.a_entries = [
                ctk.CTkEntry(fields_row, placeholder_text=lbl, width=42)
                for lbl in ("Axx", "Ayy", "Azz")
            ]
            for e in self.a_entries:
                e.insert(0, "1")
        # A entries start packed (shown); _check_spin_warning below hides
        # them immediately if the (empty) isotope field resolves to I=0,
        # which it won't yet since it's blank - shown by default is right.
        for e in self.a_entries:
            e.pack(side="left", padx=3)
            e.bind("<KeyRelease>", self._fire_on_change)

        self.weight_entry = ctk.CTkEntry(fields_row, width=42, placeholder_text="100")
        self.weight_entry.insert(0, "100")
        self.weight_entry.pack(side="left", padx=(6, 1))
        self.weight_entry.bind("<KeyRelease>", self._fire_on_change)
        ctk.CTkLabel(fields_row, text="%", text_color="gray60").pack(side="left", padx=(0, 6))

        balance_btn = ctk.CTkButton(
            fields_row, text="\u2192100%", width=46, fg_color="#3A5A78", hover_color="#2E4A63",
            command=self._request_balance,
        )
        balance_btn.pack(side="left", padx=(0, 6))
        Tooltip(balance_btn, "Set this row's % to whatever makes all rows sum to "
                              "100% (100 minus the total of every other row). "
                              "Handy for folding an unlisted no-hyperfine isotope "
                              "into an existing row with the same signature, "
                              "instead of leaving an invisible extra component.")

        self.remove_btn = ctk.CTkButton(
            fields_row, text="\u2715", width=24, fg_color="#8B2E2E", hover_color="#6E2222",
            command=self._remove,
        )
        self.remove_btn.pack(side="left")

        # v6.5: FIXED - this block had been accidentally indented into
        # _toggle_override_fields() since v6.3, so self.spin_warning only
        # existed after toggling the override checkbox; until then every
        # isotope edit raised AttributeError and I=0 A-field hiding never ran.
        # v4.9/v5.2: when the entered isotope has NO nuclear spin (I=0),
        # its A field(s) would have zero effect on the simulation - since
        # that can't be fixed by just warning (caught out an earlier
        # '48Ti' + A-value script that looked reasonable but did nothing),
        # v5.2 goes further than a lock/gray-out and removes the A box(es)
        # from view entirely for such a row, and the generated script
        # omits Sys.A for it too (see generate()) - there's no field left
        # for a value to be typed into, so it can't be gotten wrong.
        self.spin_warning = ctk.CTkLabel(
            self, text="", text_color="#D9A441", font=ctk.CTkFont(size=11))
        self.spin_warning.grid(row=1, column=0, sticky="w", padx=(2, 0))

        # v6.3: NEW - optional per-component g/lwpp override. By default a
        # row's component just inherits the shared Sys.g/Sys.lwpp (via
        # "Sys2 = Sys;" in generate()) - fine for isotope satellites of
        # the same species, which really do share one true g. But testing
        # whether a shoulder/satellite is a genuinely DIFFERENT species
        # needs its own starting g and/or lwpp, not just a shared value
        # that's later allowed to drift during fitting (that's what the
        # Fit section's "Independent g" checkbox, v6.2, does - this is
        # the complementary piece: a different STARTING point, before any
        # fitting happens at all). Unchecked by default - no behavior
        # change unless explicitly opted into.
        override_row = ctk.CTkFrame(self, fg_color="transparent")
        override_row.grid(row=2, column=0, sticky="w", padx=(2, 0), pady=(2, 0))
        self.override_var = ctk.BooleanVar(value=False)
        override_check = ctk.CTkCheckBox(
            override_row, text="Own starting g/lwpp for this component",
            variable=self.override_var, width=220, command=self._toggle_override_fields)
        override_check.pack(side="left")
        self.g_override_entry = ctk.CTkEntry(override_row, width=70, placeholder_text="g (blank=shared)")
        self.lwpp_override_entry = ctk.CTkEntry(override_row, width=110,
                                                  placeholder_text="lwpp mT (blank=shared)")
        self.g_override_entry.bind("<KeyRelease>", self._fire_on_change)
        self.lwpp_override_entry.bind("<KeyRelease>", self._fire_on_change)
        # Fields start hidden (unchecked default) - packed in
        # _toggle_override_fields when the checkbox is turned on.

    def _toggle_override_fields(self):
        if self.override_var.get():
            self.g_override_entry.pack(side="left", padx=(8, 3))
            self.lwpp_override_entry.pack(side="left", padx=(3, 0))
        else:
            self.g_override_entry.pack_forget()
            self.lwpp_override_entry.pack_forget()
        self._fire_on_change()


    def _check_spin_warning(self):
        mass, element = parse_isotope_label(self.isotope_entry.get())
        if mass is None:
            self.spin_warning.configure(text="")
            self._set_a_visible(True)
            return
        spin = lookup_isotope_spin(element, mass)
        if spin == '0':
            # v5.2: A has no effect for a spinless nucleus. Rather than
            # locking a visible field at 0 (still a box a user could stare
            # at and second-guess), the A field(s) are removed from the
            # row entirely for as long as the isotope resolves to I=0 -
            # there's no box left to put a value into, so a mixed list
            # like 48Ti/47Ti/49Ti simply never offers A on the 48Ti row,
            # only on the I>0 rows (47Ti, 49Ti). The underlying entries
            # still hold "0" internally (used by generate(), which also
            # skips writing Sys.A at all for an I=0 component - see
            # generate()'s nucs_block handling).
            for e in self.a_entries:
                if e.get().strip() != "0":
                    e.delete(0, "end")
                    e.insert(0, "0")
            self._set_a_visible(False)
            self.spin_warning.configure(
                text=f"\u26a0 {mass}{element} has I=0 (no nuclear spin) - no A needed.")
        else:
            self._set_a_visible(True)
            self.spin_warning.configure(text="")

    def _set_a_visible(self, visible):
        if visible == getattr(self, "_a_visible", True) and hasattr(self, "_a_visible"):
            return
        self._a_visible = visible
        if visible:
            for e in self.a_entries:
                if not e.winfo_ismapped():
                    e.pack(side="left", padx=3, before=self.weight_entry)
        else:
            for e in self.a_entries:
                e.pack_forget()

    def _request_balance(self):
        if self.balance_callback:
            self.balance_callback(self)

    def _fire_on_change(self, event=None):
        self._check_spin_warning()
        if self.on_change:
            self.on_change()

    def _remove(self):
        if self.on_remove:
            self.on_remove(self)
        self.destroy()

    def _open_periodic_table(self, event=None):
        def fill_isotope(label):
            self.isotope_entry.delete(0, "end")
            self.isotope_entry.insert(0, label)
            self._fire_on_change()
        PeriodicTablePicker(self, on_select=fill_isotope)
        return "break"  # prevent the click from also placing a text cursor

    def get_override(self):
        """Return (g_override_str_or_None, lwpp_override_str_or_None).
        Both None if the "own starting g/lwpp" checkbox is unchecked, or
        if a field is checked but left blank (falls back to inherited
        shared value either way - blank isn't treated as "override to
        empty"). Kept separate from get_values() deliberately, so the
        many existing (iso, n, a, weight) unpacking call sites throughout
        generate() don't all need to change shape for this add-on."""
        if not self.override_var.get():
            return None, None
        g_val = self.g_override_entry.get().strip() or None
        lwpp_val = self.lwpp_override_entry.get().strip() or None
        return g_val, lwpp_val

    def set_override(self, g_override=None, lwpp_override=None):
        """For preset load. Leaves the checkbox unchecked/fields blank if
        both are None/empty (backward-compatible with presets saved
        before this existed)."""
        has_override = bool(g_override) or bool(lwpp_override)
        self.override_var.set(has_override)
        self._toggle_override_fields()
        if g_override:
            self.g_override_entry.delete(0, "end")
            self.g_override_entry.insert(0, str(g_override))
        if lwpp_override:
            self.lwpp_override_entry.delete(0, "end")
            self.lwpp_override_entry.insert(0, str(lwpp_override))

    def get_values(self):
        """Return (isotope, n, [A values in mT as strings, possibly empty],
        weight_pct_str) or None if the isotope field itself is blank.
        Does NOT validate that A is filled in - that's checked later in
        generate(), deliberately AFTER same-element weight auto-detection
        runs, so the % field still auto-fills even if A isn't filled in
        yet (auto-detection only needs the isotope label, not A)."""
        isotope = self.isotope_entry.get().strip()
        if not isotope:
            return None
        n = self.n_entry.get().strip() or "1"
        a_vals = [e.get().strip() for e in self.a_entries]
        weight = self.weight_entry.get().strip() or "100"
        return isotope, n, a_vals, weight

    def set_values(self, isotope, n, a_vals, weight="100"):
        self.isotope_entry.delete(0, "end")
        self.isotope_entry.insert(0, isotope)
        self.n_entry.delete(0, "end")
        self.n_entry.insert(0, str(n))
        for e, v in zip(self.a_entries, a_vals):
            e.delete(0, "end")
            e.insert(0, str(v))
        self.weight_entry.delete(0, "end")
        self.weight_entry.insert(0, str(weight))
        self._check_spin_warning()


class NucleusList(ctk.CTkFrame):
    def __init__(self, master, mode="iso", on_change=None, **kwargs):
        super().__init__(master, fg_color="transparent", **kwargs)
        self.mode = mode
        self.on_change = on_change
        self.rows = []

        self.rows_frame = ctk.CTkFrame(self, fg_color="transparent")
        self.rows_frame.pack(fill="x", expand=True)

        self.add_btn = ctk.CTkButton(self, text="+ Add nucleus", width=140, command=self.add_row)
        self.add_btn.pack(anchor="w", pady=(6, 0))

        self.add_row()

    def add_row(self):
        row = NucleusRow(self.rows_frame, mode=self.mode, on_remove=self._remove_row,
                          on_change=self.on_change, balance_callback=self._balance_row)
        row.pack(anchor="w", pady=2)
        self.rows.append(row)
        return row

    def _balance_row(self, target_row):
        """Set target_row's weight % to whatever makes all rows sum to
        100% (100 minus the total of every other row's current weight)."""
        total_others = 0.0
        for row in self.rows:
            if row is target_row:
                continue
            try:
                total_others += float(row.weight_entry.get().strip() or "0")
            except ValueError:
                pass
        remainder = 100.0 - total_others
        target_row.weight_entry.delete(0, "end")
        target_row.weight_entry.insert(0, f"{remainder:.6g}")
        if self.on_change:
            self.on_change()

    def _remove_row(self, row):
        if row in self.rows:
            self.rows.remove(row)
        if self.on_change:
            self.on_change()

    def clear(self):
        for row in list(self.rows):
            row.destroy()
        self.rows = []

    def collect(self):
        """Return list of (isotope, n, [A values in mT], weight_pct_str)
        for filled rows."""
        out = []
        for row in self.rows:
            vals = row.get_values()
            if vals is not None:
                out.append(vals)
        return out

    def collect_with_overrides(self):
        """Same order/filtering as collect(), but also returns each row's
        (g_override, lwpp_override) from get_override() (None, None if
        that row's checkbox is unchecked). Returns a list of
        (values_tuple, (g_override, lwpp_override)) pairs, kept separate
        from collect() so existing callers of collect() are unaffected."""
        out = []
        for row in self.rows:
            vals = row.get_values()
            if vals is not None:
                out.append((vals, row.get_override()))
        return out

    def load(self, nuclei_data):
        """nuclei_data: list of dicts {isotope, n, a_mT: [...], weight}
        (the "weight" key is optional, for backward compatibility with
        presets saved before per-nucleus weight existed - defaults to
        100; "g_override"/"lwpp_override" are likewise optional, for
        presets saved before v6.3 - default to no override, same as
        leaving the checkbox unchecked)."""
        self.clear()
        if not nuclei_data:
            self.add_row()
            return
        for entry in nuclei_data:
            row = self.add_row()
            row.set_values(entry["isotope"], entry["n"], entry["a_mT"],
                            entry.get("weight", "100"))
            row.set_override(entry.get("g_override"), entry.get("lwpp_override"))


# --------------------------------------------------------------------------- #
# Dataset loader panel: browse folder, parse .DSC, plot experimental trace
# --------------------------------------------------------------------------- #
class ResidualBuilderDialog(ctk.CTkToplevel):
    """
    Loads a 'Field_mT,Experimental,Fit' CSV - the exact format this
    tool's own esfit Part 2 export produces - computes
    Residual = Experimental - Fit, previews all three traces, and lets
    you save the residual as a plain data file (index, Field [G],
    Intensity - the same format this tool's own dataset loader expects)
    so it can be fed back in as a new "experimental" dataset for fitting
    a minor component separately (e.g. isolating a minority isotope's
    hyperfine satellites after subtracting the majority isotope's fit).
    """

    def __init__(self, master, on_use_residual=None):
        super().__init__(master)
        self.title("Build Residual (Experimental - Fit)")
        self.geometry("920x680")
        self.on_use_residual = on_use_residual
        self.field_mT = None
        self.experimental = None
        self.fit = None
        self.residual = None

        top = ctk.CTkFrame(self, fg_color="transparent")
        top.pack(fill="x", padx=14, pady=(14, 6))
        ctk.CTkButton(top, text="Browse for Fit-export CSV...",
                      command=self._browse).pack(side="left")
        self.file_label = ctk.CTkLabel(top, text="No file loaded", text_color="gray60")
        self.file_label.pack(side="left", padx=10)

        info_note(self, "Expects the 'Field_mT,Experimental,Fit' CSV format this "
                         "tool's own esfit Part 2 export produces (any CSV with "
                         "those three comma-separated columns, in that order, works "
                         "too).")

        self.plot_frame = ctk.CTkFrame(self, fg_color="transparent")
        self.plot_frame.pack(fill="both", expand=True, padx=14, pady=6)
        self._show_placeholder()

        btn_row = ctk.CTkFrame(self, fg_color="transparent")
        btn_row.pack(fill="x", padx=14, pady=(6, 14))
        self.save_btn = ctk.CTkButton(btn_row, text="Save residual as data file...",
                                       command=self._save, state="disabled")
        self.save_btn.pack(side="left", padx=(0, 8))
        self.use_btn = ctk.CTkButton(
            btn_row, text="Use as active dataset for fitting", command=self._use,
            state="disabled", fg_color="#2E7D32", hover_color="#245F27")
        self.use_btn.pack(side="left")

        self.transient(master.winfo_toplevel())
        self.lift()
        self.focus_force()

    def _show_placeholder(self):
        for w in self.plot_frame.winfo_children():
            w.destroy()
        ctk.CTkLabel(self.plot_frame, text="Browse for a Fit-export CSV to preview "
                                            "Experimental / Fit / Residual here.",
                     text_color="gray60").pack(expand=True)

    def _browse(self):
        path = filedialog.askopenfilename(
            title="Select Fit-export CSV",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
        if not path:
            return
        try:
            field_mT, experimental, fit = parse_fit_export_csv(path)
        except Exception as e:
            messagebox.showerror("Parse error", str(e))
            return
        self.field_mT = field_mT
        self.experimental = experimental
        self.fit = fit
        self.residual = [e - f for e, f in zip(experimental, fit)]
        self.file_label.configure(text=f"{os.path.basename(path)} ({len(field_mT)} points)")
        self._plot()
        self.save_btn.configure(state="normal")
        self.use_btn.configure(state="normal")

    def _plot(self):
        for w in self.plot_frame.winfo_children():
            w.destroy()
        fig = Figure(figsize=(8.5, 5), dpi=100)
        ax = fig.add_subplot(111)
        ax.plot(self.field_mT, self.experimental, label="Experimental",
                color="black", linewidth=0.8)
        ax.plot(self.field_mT, self.fit, label="Fit (major component)",
                color="tab:blue", linewidth=0.8)
        ax.plot(self.field_mT, self.residual, label="Residual (minor component)",
                color="tab:red", linewidth=0.8)
        ax.set_xlabel("Field (mT)")
        ax.set_ylabel("Intensity (a.u.)")
        ax.legend(fontsize=9)
        fig.tight_layout()
        canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)

    def _save(self):
        path = filedialog.asksaveasfilename(
            title="Save residual data file (plain text, no extension needed)")
        if not path:
            return
        write_residual_data_file(path, self.field_mT, self.residual)
        messagebox.showinfo(
            "Saved",
            f"Residual saved to:\n{path}\n\n"
            "This is the same plain-text format (index, Field [G], Intensity) "
            "this tool's dataset loader expects. The experimental parameters "
            "(mwFreq, field range, nPoints) are unchanged from the original "
            "measurement - either use 'Use as active dataset' here, or load it "
            "manually by browsing for the ORIGINAL .DSC file and then replacing "
            "the loaded data with this file's contents."
        )

    def _use(self):
        if self.on_use_residual:
            self.on_use_residual(self.field_mT, self.residual)
        self.destroy()


class DatasetPanel(ctk.CTkFrame):
    """
    A shared panel for loading an experimental dataset folder. Parses the
    .DSC file for Exp auto-fill values, loads the experimental trace for
    preview, and derives a compound ID (overridable) used for preset
    filenames.

    on_loaded(exp_values: dict, compound_id: str, folder: str) is called
    after a successful load so the parent tab can auto-fill its Exp fields.
    """

    def __init__(self, master, on_loaded=None, **kwargs):
        super().__init__(master, **kwargs)
        self.on_loaded = on_loaded
        self.dataset_folder = None
        self.data_path = None
        self.field_mT = None
        self.intensity = None
        self.exp_values = {}

        section_title(self, "Dataset")

        browse_row = ctk.CTkFrame(self, fg_color="transparent")
        browse_row.pack(anchor="w", fill="x", pady=3)
        ctk.CTkButton(browse_row, text="Browse for .DSC file...",
                       command=self.browse).pack(side="left")
        self.folder_label = ctk.CTkLabel(browse_row, text="(no dataset loaded)",
                                          text_color="gray60", wraplength=280,
                                          justify="left", anchor="w")
        self.folder_label.pack(side="left", padx=10)

        residual_row = ctk.CTkFrame(self, fg_color="transparent")
        residual_row.pack(anchor="w", fill="x", pady=3)
        ctk.CTkButton(residual_row, text="Build Residual from Fit...",
                       command=self._open_residual_builder).pack(side="left")
        info_icon(residual_row, "Load a 'Field_mT, Experimental, Fit' CSV (the format "
                                 "this tool's own esfit Part 2 export produces) to "
                                 "compute Residual = Experimental - Fit, then use that "
                                 "residual directly as the active dataset here for "
                                 "fitting a minor component separately."
                  ).pack(side="left", padx=(6, 0))

        id_row = ctk.CTkFrame(self, fg_color="transparent")
        id_row.pack(anchor="w", fill="x", pady=3)
        ctk.CTkLabel(id_row, text="Compound ID:", width=190, anchor="w").pack(side="left")
        self.compound_id_entry = ctk.CTkEntry(id_row, width=180,
                                               placeholder_text="e.g. MS-1-33")
        self.compound_id_entry.pack(side="left", padx=6)

        self.status_label = ctk.CTkLabel(self, text="", text_color="gray60",
                                          wraplength=400, justify="left", anchor="w")
        self.status_label.pack(anchor="w", pady=(4, 4))

    def browse(self):
        dsc_path = filedialog.askopenfilename(
            title="Select a .DSC file (pick the specific dataset, not just a folder)",
            filetypes=[("EPR parameter file", "*.DSC *.dsc"), ("All files", "*.*")],
        )
        if not dsc_path:
            return
        self.load_dsc_file(dsc_path)

    def load_dsc_file(self, dsc_path):
        """Load a specific .DSC file the user picked (not a folder scan)."""
        if not os.path.isfile(dsc_path):
            messagebox.showerror("File not found", f"Could not find:\n{dsc_path}")
            return

        try:
            exp_values = parse_dsc(dsc_path)
        except Exception as e:
            messagebox.showerror("Error parsing .DSC", str(e))
            return

        data_path = find_data_companion(dsc_path)
        folder = os.path.dirname(dsc_path)

        self.dataset_folder = folder
        self.data_path = data_path
        self.exp_values = exp_values
        self.folder_label.configure(text=os.path.basename(dsc_path))

        compound_id = exp_values.get("title") or (
            compound_id_from_filename(data_path) if data_path
            else compound_id_from_filename(dsc_path)
        )
        self.compound_id_entry.delete(0, "end")
        self.compound_id_entry.insert(0, compound_id)

        if data_path:
            try:
                self.field_mT, self.intensity = load_experimental_trace(data_path)
                self.status_label.configure(
                    text=f"Loaded {len(self.field_mT)} points from experimental trace. "
                         f"Exp fields auto-filled below (override as needed)."
                )
            except Exception as e:
                self.field_mT, self.intensity = None, None
                self.status_label.configure(
                    text=f"Found data file but could not parse trace: {e}"
                )
        else:
            self.field_mT, self.intensity = None, None
            self.status_label.configure(
                text="No companion data file found alongside .DSC; "
                     "Exp fields auto-filled from .DSC only."
            )

        if self.on_loaded:
            self.on_loaded(exp_values, compound_id, folder)

    def get_compound_id(self):
        return self.compound_id_entry.get().strip()

    def has_experimental_trace(self):
        return self.field_mT is not None and self.intensity is not None

    def _open_residual_builder(self):
        ResidualBuilderDialog(self, on_use_residual=self.use_residual)

    def use_residual(self, field_mT, residual):
        """
        Save `residual` to a new plain-text data file and switch this
        panel's active dataset to it, re-applying the ORIGINAL Exp values
        (mwFreq, field range, nPoints) to the tab's actual entry fields -
        the residual is the same physical measurement with a component's
        fit subtracted, not a new experiment, but the residual CSV itself
        carries no Exp metadata at all, so those fields must come from
        whatever .DSC was loaded earlier in THIS session. If none was
        (self.exp_values is empty), refuses with a clear message rather
        than silently leaving the Exp fields at their generic defaults -
        that mismatch (e.g. nPoints=1024 default vs. a 15000-point
        residual) is exactly what produces esfit's "must match exactly"
        error. Returns True if a file was actually saved.
        """
        if not self.exp_values:
            messagebox.showerror(
                "No original dataset loaded",
                "The residual file itself doesn't carry experimental "
                "parameters (mwFreq, field range, nPoints) - only the field/"
                "intensity trace does. Please load the ORIGINAL .DSC file "
                "first (via 'Browse for .DSC file...' above) in this "
                "session, so those values are known, before using a "
                "residual as the active dataset."
            )
            return False

        default_name = f"{self.get_compound_id() or 'dataset'}_residual"
        path = filedialog.asksaveasfilename(
            title="Save residual as the new active dataset (no extension needed)",
            initialfile=default_name,
            initialdir=self.dataset_folder or None,
        )
        if not path:
            return False

        write_residual_data_file(path, field_mT, residual)

        self.data_path = path
        self.field_mT = list(field_mT)
        self.intensity = list(residual)
        self.dataset_folder = os.path.dirname(path)

        current_id = self.get_compound_id()
        new_id = current_id + "_residual" if current_id and not current_id.endswith(
            "_residual") else (current_id or "dataset_residual")
        self.compound_id_entry.delete(0, "end")
        self.compound_id_entry.insert(0, new_id)

        self.folder_label.configure(text=f"{os.path.basename(path)} (residual)")
        self.status_label.configure(
            text=f"Now using the residual ({len(self.field_mT)} points) as the "
                 f"active dataset - Exp values re-applied from the original "
                 f".DSC. Fit this for the minor component below."
        )

        # Re-push the ORIGINAL Exp values into the tab's actual entry
        # fields - guards against them having been edited/cleared between
        # loading the .DSC and building the residual, and is what actually
        # fixes the mismatch (not just updating self.exp_values, which
        # nothing downstream reads directly - generate() reads the entry
        # widgets themselves).
        if self.on_loaded:
            self.on_loaded(self.exp_values, new_id, self.dataset_folder)

        return True


# --------------------------------------------------------------------------- #
# Preset save/load (JSON, stored alongside the dataset folder)
# --------------------------------------------------------------------------- #
def preset_path(folder, compound_id, function_name):
    safe_id = compound_id.replace("/", "_").replace("\\", "_").strip() or "untitled"
    return os.path.join(folder, f"{safe_id}_{function_name}.json")


def save_preset(folder, compound_id, function_name, data):
    if not folder:
        raise ValueError("No dataset folder loaded - load a dataset before saving a preset.")
    path = preset_path(folder, compound_id, function_name)
    data = dict(data)
    data["_saved"] = datetime.datetime.now().isoformat(timespec="seconds")
    data["_function"] = function_name
    data["_compound_id"] = compound_id
    with open(path, "w") as f:
        json.dump(data, f, indent=2)
    return path


def load_preset_dialog(folder, function_name):
    initial_dir = folder if folder else os.path.expanduser("~")
    path = filedialog.askopenfilename(
        title=f"Load {function_name} preset",
        initialdir=initial_dir,
        filetypes=[("EasySpin preset", "*.json")],
    )
    if not path:
        return None
    with open(path, "r") as f:
        data = json.load(f)
    return data


# --------------------------------------------------------------------------- #
# esfit (interactive fitting) section - shared by garlic/pepper/chili tabs
#
# Built from a researched easyspin_esfit_reference.md (EasySpin 6.0.12,
# checked 2026-08-14). Key confirmed facts this relies on:
#   - Call form: esfit(data, fcn, {Sys0,Exp0}, {Vary}, FitOpt)
#   - No output argument requested -> opens the interactive GUI (this is
#     the actual documented trigger, not a separate flag)
#   - Vary fields mirror Sys0 field names; each value is a symmetric +/-
#     absolute delta (same units as the Sys0 field), NOT a percentage
#   - data must be a plain intensity vector; Exp.nPoints/Range must match
#     the real acquisition exactly (no resampling is documented)
#   - chili should vary logtcorr/logDiff, not tcorr/Diff, per the docs
#   - FitOpt.Method combines an algorithm keyword AND a target keyword in
#     one space-separated string (e.g. 'simplex fcn'); esfit's own
#     internal default target has been observed as NOT 'fcn' (its GUI
#     showed 'integral' when no target was specified), so a target is
#     always written explicitly rather than left to chance
# --------------------------------------------------------------------------- #
FITOPT_METHODS = ["simplex", "levmar", "montecarlo", "genetic", "grid"]
FITOPT_TARGETS = ["fcn", "int", "dint", "diff", "fft"]
FITOPT_TARGET_LABELS = {
    "fcn": "fcn (fit the raw spectrum as-is - default)",
    "int": "int (fit the integral of the spectrum)",
    "dint": "dint (fit the double integral)",
    "diff": "diff (fit the derivative)",
    "fft": "fft (fit the Fourier transform)",
}


class FitSection(ctk.CTkFrame):
    """
    Reusable "Fit (esfit)" UI block: an enable checkbox, one
    checkbox+delta-entry row per fittable Sys field, a FitOpt.Method
    algorithm selector, and a Target selector (what esfit actually
    compares - raw spectrum by default, or its integral/derivative/FFT).
    Call collect_vary() to get the {field: delta} dict of currently-
    checked fields.
    """

    def __init__(self, parent, fields, **kwargs):
        """
        fields: list of (sys_field_name, display_label, default_delta,
        default_unit_note) tuples, e.g.
        [("g", "g-value", "0.02", ""), ("lwpp", "Linewidth", "0.05", "mT")]
        """
        super().__init__(parent, fg_color="transparent", **kwargs)

        self.enabled_var = ctk.BooleanVar(value=True)
        top_row = ctk.CTkFrame(self, fg_color="transparent")
        top_row.pack(anchor="w", fill="x", pady=(0, 4))
        ctk.CTkCheckBox(top_row, text="Add interactive fit (esfit) to the script",
                         variable=self.enabled_var).pack(side="left")
        info_icon(top_row, "Adds an esfit(...) call using your loaded experimental "
                            "data. Called with no output argument, which is what "
                            "opens EasySpin's interactive fitting GUI (sliders, "
                            "algorithm choice, live comparison to your data). "
                            "Requires a dataset to be loaded above."
                  ).pack(side="left", padx=(6, 0))

        self.field_rows = {}
        for name, label, default_delta, unit in fields:
            row = ctk.CTkFrame(self, fg_color="transparent")
            row.pack(anchor="w", fill="x", pady=2)
            var = ctk.BooleanVar(value=True)
            ctk.CTkCheckBox(row, text=f"Vary {label}", variable=var, width=170).pack(side="left")
            entry = ctk.CTkEntry(row, width=70, placeholder_text=default_delta)
            entry.insert(0, default_delta)
            entry.pack(side="left", padx=(6, 4))
            unit_text = f"({unit}, +/-)" if unit else "(+/-)"
            ctk.CTkLabel(row, text=unit_text, text_color="gray60",
                         font=ctk.CTkFont(size=11)).pack(side="left")
            self.field_rows[name] = (var, entry)

        # Dynamic extra vary rows for optional Sys fields that only exist
        # when a corresponding UI toggle is on (e.g. HStrain/gStrain from
        # the Line Broadening section) - unlike field_rows above, this set
        # can change at runtime, so it's rebuilt via sync_extra_fields().
        self.extra_vary_frame = ctk.CTkFrame(self, fg_color="transparent")
        self.extra_vary_frame.pack(anchor="w", fill="x")
        self._extra_vary_widgets = {}  # name -> (var, entry)
        self._extra_vary_state = {}    # name -> {"checked":bool,"delta":str}

        method_row = ctk.CTkFrame(self, fg_color="transparent")
        method_row.pack(anchor="w", fill="x", pady=(6, 0))
        ctk.CTkLabel(method_row, text="Algorithm (FitOpt.Method):", width=190,
                     anchor="w").pack(side="left")
        self.method_var = ctk.StringVar(value="simplex")
        ctk.CTkOptionMenu(method_row, variable=self.method_var,
                           values=FITOPT_METHODS).pack(side="left")
        info_icon(method_row, "simplex (Nelder-Mead) is the EasySpin default and a "
                               "solid general choice. levmar can converge faster near "
                               "a good starting guess. montecarlo/genetic/grid are "
                               "slower but better for exploring far from the starting "
                               "point."
                  ).pack(side="left", padx=(6, 0))

        target_row = ctk.CTkFrame(self, fg_color="transparent")
        target_row.pack(anchor="w", fill="x", pady=(4, 0))
        ctk.CTkLabel(target_row, text="Target/weighting:", width=190,
                     anchor="w").pack(side="left")
        self.target_var = ctk.StringVar(value=FITOPT_TARGET_LABELS["fcn"])
        ctk.CTkOptionMenu(target_row, variable=self.target_var,
                           values=list(FITOPT_TARGET_LABELS.values()),
                           width=280).pack(side="left")
        info_icon(target_row, "What esfit actually compares against your data. "
                               "Left at its default (fcn), this fits the raw "
                               "spectrum directly - the normal choice. Without this "
                               "set explicitly, esfit's own internal default may "
                               "NOT be 'fcn' (its GUI has been observed defaulting "
                               "to 'integral' when no target is specified), so this "
                               "tool always writes one explicitly rather than "
                               "leaving it to chance."
                  ).pack(side="left", padx=(6, 0))

        section_title(self, "Vary hyperfine A (per nucleus)")
        info_note(self, "Only available when there's exactly one nucleus/weighted "
                         "component total - esfit's behavior with a multi-component "
                         "Sys hasn't been verified, so fitting is disabled rather "
                         "than guessed at whenever more than one component exists.")
        self.nucleus_vary_frame = ctk.CTkFrame(self, fg_color="transparent")
        self.nucleus_vary_frame.pack(anchor="w", fill="x")
        self._nucleus_vary_state = {}  # index -> {"checked": bool, "delta": str}
        self._nucleus_vary_widgets = []  # [(var, entry), ...] in current order
        self._nucleus_labels = []

    def sync_nucleus_rows(self, labels):
        """Rebuild the per-nucleus Vary-A / Vary-weight checkboxes to match
        `labels` (the current isotope labels from the tab's NucleusList, in
        order). No-op if `labels` matches what's already displayed - this
        is called on every generate(), and rebuilding unconditionally would
        destroy the very checkbox/delta state the user just set, before
        it's ever read via get_nucleus_vary_list()."""
        if labels == self._nucleus_labels and len(self._nucleus_vary_widgets) == len(labels):
            return

        self._nucleus_labels = list(labels)
        for w in self.nucleus_vary_frame.winfo_children():
            w.destroy()
        self._nucleus_vary_widgets = []

        for i, label in enumerate(labels):
            prior = self._nucleus_vary_state.get(i, {
                "checked": True, "delta": "3",
                "weight_checked": False, "weight_delta": "3",
                "indep_g_checked": False,
            })
            row = ctk.CTkFrame(self.nucleus_vary_frame, fg_color="transparent")
            row.pack(anchor="w", pady=2)

            # v5.3: a nucleus with I=0 (no nuclear spin) has no Sys.A field
            # at all in the generated script (see generate()'s zero-spin
            # handling) - offering a "Vary A" checkbox for it let esfit be
            # told to vary a field that plain doesn't exist on that
            # component, producing "Field .A is missing in input N." No
            # checkbox/entry is created for such a row; var/entry stay
            # None and get_nucleus_vary_list() below treats that as
            # "never vary A" regardless of any stale saved state.
            if is_zero_spin_isotope(label):
                ctk.CTkLabel(row, text=f"{label}: I=0, no A to vary",
                             text_color="gray60", width=200, anchor="w").pack(side="left")
                var, entry = None, None
            else:
                var = ctk.BooleanVar(value=prior["checked"])
                ctk.CTkCheckBox(row, text=f"Vary A: {label}", variable=var, width=150).pack(side="left")
                entry = ctk.CTkEntry(row, width=60, placeholder_text="3")
                entry.insert(0, prior["delta"])
                entry.pack(side="left", padx=(4, 4))
                ctk.CTkLabel(row, text="(mT, +/-)", text_color="gray60",
                             font=ctk.CTkFont(size=11)).pack(side="left", padx=(0, 14))

            weight_var = ctk.BooleanVar(value=prior["weight_checked"])
            ctk.CTkCheckBox(row, text=f"Vary weight: {label}", variable=weight_var,
                             width=160).pack(side="left")
            weight_entry = ctk.CTkEntry(row, width=60, placeholder_text="3")
            weight_entry.insert(0, prior["weight_delta"])
            weight_entry.pack(side="left", padx=(4, 4))
            ctk.CTkLabel(row, text="(%, +/-)", text_color="gray60",
                         font=ctk.CTkFont(size=11)).pack(side="left", padx=(0, 14))

            # v6.2: NEW - opt-in override of the v6.0 shared-g restriction.
            # By default g is only free on the FIRST component (isotope
            # satellites of the same species must share one true g - see
            # the v6.0 note in generate()). But sometimes the point of the
            # fit IS to test whether a component is a genuinely different
            # species with its OWN g (e.g. checking if a shoulder is a
            # distinct minor species rather than a Ti isotope satellite).
            # Checking this lets that one component keep the SAME delta
            # typed into the main "Vary g-value" field above, instead of
            # having it stripped. Off by default - this should be a
            # deliberate choice, not a silent default, since letting every
            # component's g float independently is exactly the v6.0 bug.
            indep_g_var = ctk.BooleanVar(value=prior["indep_g_checked"])
            ctk.CTkCheckBox(row, text=f"Independent g: {label}", variable=indep_g_var,
                             width=170).pack(side="left")

            self._nucleus_vary_widgets.append((var, entry, weight_var, weight_entry, indep_g_var))

        if not labels:
            ctk.CTkLabel(self.nucleus_vary_frame, text="(no nuclei defined yet)",
                         text_color="gray60").pack(anchor="w")

    def get_nucleus_vary_list(self):
        """Return [(A_checked, A_delta, weight_checked, weight_delta_pct,
        indep_g_checked), ...] in the same order as the last
        sync_nucleus_rows() call, and persist the current state so it
        survives the next rebuild. weight_delta_pct is a percentage
        string (e.g. "3" for +/-3%) - convert to a fraction (divide by
        100) before writing it as a Sys.weight-scale delta. A row whose
        isotope is I=0 has var/entry set to None (see sync_nucleus_rows)
        - always reported as unchecked here, since that component has no
        Sys.A field for esfit to vary in the first place. indep_g_checked
        (v6.2) is False by default - see the matching note in
        sync_nucleus_rows for what it does."""
        out = []
        for i, (var, entry, weight_var, weight_entry, indep_g_var) in enumerate(self._nucleus_vary_widgets):
            checked = var.get() if var is not None else False
            delta = (entry.get().strip() or "3") if entry is not None else "3"
            weight_checked = weight_var.get()
            weight_delta = weight_entry.get().strip() or "3"
            indep_g_checked = indep_g_var.get()
            self._nucleus_vary_state[i] = {
                "checked": checked, "delta": delta,
                "weight_checked": weight_checked, "weight_delta": weight_delta,
                "indep_g_checked": indep_g_checked,
            }
            out.append((checked, delta, weight_checked, weight_delta, indep_g_checked))
        return out


    def is_enabled(self):
        return self.enabled_var.get()

    def collect_vary(self):
        """Return {sys_field_name: delta_string} for every checked field
        with a non-empty delta value (includes the dynamic extra fields
        from sync_extra_fields(), e.g. HStrain/gStrain when active)."""
        out = {}
        for name, (var, entry) in self.field_rows.items():
            if var.get():
                val = entry.get().strip()
                if val:
                    out[name] = val
        out.update(self.get_extra_vary())
        return out

    def sync_extra_fields(self, fields):
        """Rebuild the dynamic extra-field vary rows to match `fields` -
        a list of (sys_field_name, display_label, default_delta, unit)
        for whichever optional Sys fields are CURRENTLY active (e.g. only
        HStrain if that broadening checkbox is on, both if gStrain is
        also on, neither otherwise). No-op if the active field set hasn't
        changed, so it doesn't wipe user edits on every generate() call -
        same reasoning as sync_nucleus_rows(). State (checked/delta) is
        preserved by field NAME across rebuilds, safe here since these
        are a small fixed set of named fields, not user-typed labels."""
        target_names = [f[0] for f in fields]
        if target_names == list(self._extra_vary_widgets.keys()):
            return

        for name, (var, entry) in self._extra_vary_widgets.items():
            self._extra_vary_state[name] = {
                "checked": var.get(), "delta": entry.get().strip() or "0.1"}

        for w in self.extra_vary_frame.winfo_children():
            w.destroy()
        self._extra_vary_widgets = {}

        for name, label, default_delta, unit in fields:
            prior = self._extra_vary_state.get(
                name, {"checked": True, "delta": default_delta})
            row = ctk.CTkFrame(self.extra_vary_frame, fg_color="transparent")
            row.pack(anchor="w", pady=2)
            var = ctk.BooleanVar(value=prior["checked"])
            ctk.CTkCheckBox(row, text=f"Vary {label}", variable=var, width=170).pack(side="left")
            entry = ctk.CTkEntry(row, width=70, placeholder_text=default_delta)
            entry.insert(0, prior["delta"])
            entry.pack(side="left", padx=(6, 4))
            unit_text = f"({unit}, +/-)" if unit else "(+/-)"
            ctk.CTkLabel(row, text=unit_text, text_color="gray60",
                         font=ctk.CTkFont(size=11)).pack(side="left")
            self._extra_vary_widgets[name] = (var, entry)

    def get_extra_vary(self):
        """Return {sys_field_name: delta_string} for every checked
        dynamic extra field (independent of collect_vary(), so callers
        that need JUST the extras - e.g. to decide whether to add a
        Sys.HStrain line at all - can call this directly)."""
        out = {}
        for name, (var, entry) in self._extra_vary_widgets.items():
            if var.get():
                val = entry.get().strip()
                if val:
                    out[name] = val
        return out

    def get_method(self):
        return self.method_var.get()

    def get_target(self):
        label = self.target_var.get()
        for keyword, lbl in FITOPT_TARGET_LABELS.items():
            if lbl == label:
                return keyword
        return "fcn"


def format_matlab_vector(values, per_line=8, indent="    "):
    """Format a list of floats as a MATLAB row-vector literal spanning
    multiple lines with '...' continuation, e.g. for embedding a large
    experimental data array directly in a generated script."""
    lines = []
    for i in range(0, len(values), per_line):
        chunk = values[i:i + per_line]
        chunk_str = ", ".join(f"{v:.10g}" for v in chunk)
        is_last = (i + per_line) >= len(values)
        suffix = "" if is_last else " ..."
        lines.append(f"{indent}{chunk_str}{suffix}")
    return "[\n" + "\n".join(lines) + "\n]"


def build_esfit_block(fcn_name, sim_sys_expr, vary_specs, method, data_path,
                       field_lo, field_hi, npoints, chili_mode=False,
                       has_opt=False, target="fcn"):
    """
    Build the MATLAB code block (as a list of lines) for:
      1. loading the experimental data directly from its file at MATLAB
         runtime (mirrors this GUI's own parser: keep only rows with
         exactly 3 numeric columns, discarding headers/blank lines - so
         it doesn't depend on a fixed header-line count)
      2. the esfit(...) call itself (no output arg -> opens the
         interactive GUI), using Sys/Exp (and Opt, if has_opt) already
         defined earlier in the script.

    Returns (main_lines, post_export_lines) - two SEPARATE scripts meant
    to be run one after the other in the same MATLAB session: the first
    loads data and opens the interactive fit; the second (run only after
    you've fit interactively and clicked "Export" in the GUI) pulls out
    the optimized trace and saves it. They're split because MATLAB
    variables from the first script (Sys, Exp, B, field_axis, data) must
    still be in the workspace when the second one runs, but the GUI
    itself doesn't block execution, so they can't be one linear script.

    sim_sys_expr: the exact Sys expression already used for the plain
    simulation call - either "Sys" (single component) or
    "{Sys, Sys2, Sys3}" (multi-component). Reused as-is for esfit's own
    Sys0 argument, since per esfit's documentation the same {Sys1,Sys2}
    cell-array form is valid there too: esfit(data,@pepper,
    {{Sys1,Sys2},Exp0},{{varySys1,varySys2}}).

    vary_specs: list of per-component {field: delta} dicts, in the SAME
    order as the components inside sim_sys_expr (one dict per component,
    even if sim_sys_expr is a single plain "Sys" - then this list has
    exactly one entry). An empty dict for a component becomes a literal
    [] in the generated Vary cell array, per esfit's documented
    convention for "nothing varies in this component".

    data_path: absolute path to the plain-text data file (index, Field
    [G], Intensity columns), embedded as a MATLAB string literal. Loaded
    at runtime rather than embedding the values themselves, so the
    script stays small and stays in sync if the data file is re-exported.

    has_opt: True if the script already defines an Opt struct (pepper's
    GridSize) that should be passed inside the {Sys,Exp,Opt} cell, per
    the documented esfit(data,@pepper,{Sys0,Exp0,Opt},{varySys}) form.

    chili_mode: if True, silently maps tcorr/Diff Vary entries to
    logtcorr/logDiff per the documented chili-specific fitting guidance.
    """
    # MATLAB single-quoted strings don't interpret backslashes, so a
    # literal Windows path can be inserted as-is with no escaping. Still
    # normalize it first, since os.path.join() can otherwise produce a
    # mixed-separator path (e.g. forward slashes from a Tk file dialog
    # joined with a backslash from os.path.join on Windows).
    data_path = os.path.normpath(data_path)
    matlab_path = data_path.replace("'", "''")  # escape literal quotes only

    lines = []
    lines.append("% --- Experimental data (for esfit) ---")
    lines.append("% Loaded directly from your dataset's data file at runtime (kept")
    lines.append("% out of the script itself so it stays small and stays in sync if")
    lines.append("% you re-export the data). Parsing mirrors the GUI's own reader:")
    lines.append("% any row that isn't exactly 3 numeric columns is skipped, so this")
    lines.append("% doesn't depend on a fixed header-line count. If your lab's export")
    lines.append("% format differs (different column count/order), adjust below.")
    lines.append(f"data_file = '{matlab_path}';")
    lines.append("fid = fopen(data_file, 'r');")
    lines.append("raw_lines = textscan(fid, '%s', 'Delimiter', '\\n');")
    lines.append("fclose(fid);")
    lines.append("raw_lines = raw_lines{1};")
    lines.append("")
    lines.append("field_G = [];")
    lines.append("data = [];")
    lines.append("for i = 1:numel(raw_lines)")
    lines.append("    parts = strsplit(strtrim(raw_lines{i}));")
    lines.append("    if numel(parts) == 3")
    lines.append("        vals = str2double(parts);")
    lines.append("        if ~any(isnan(vals))")
    lines.append("            field_G(end+1) = vals(2);")
    lines.append("            data(end+1) = vals(3);")
    lines.append("        end")
    lines.append("    end")
    lines.append("end")
    lines.append("field_axis = field_G / 10;  % Gauss -> mT")
    lines.append("data = data(:);  % ensure column vector")
    lines.append("")
    lines.append("% --- Restrict to the fitted field window ---")
    lines.append("% The block above loads EVERY row from the raw data file, spanning")
    lines.append("% whatever field width was originally recorded - NOT just the window")
    lines.append("% you're fitting over. Without this step, esfit would compare a")
    lines.append("% simulation computed over that (possibly cropped) window against")
    lines.append("% experimental data spanning the ORIGINAL full width: same array")
    lines.append("% length, different physical field width per index - which silently")
    lines.append("% stretches/compresses the apparent lineshape instead of raising an")
    lines.append("% error. The bounds below are embedded as literal numbers (matching")
    lines.append("% whatever Exp.Range/Exp.CenterSweep was set to above) rather than")
    lines.append("% read back off the Exp struct, since garlic uses Exp.CenterSweep")
    lines.append("% while pepper/chili use Exp.Range - reading 'Exp.Range' here")
    lines.append("% unconditionally broke garlic scripts with 'Unrecognized field")
    lines.append("% name \"Range\"'. Slicing both field_axis and data here, and")
    lines.append("% deriving Exp.nPoints from what's actually left, guarantees the")
    lines.append("% data and simulation axes always agree - regardless of whatever")
    lines.append("% Exp.nPoints was typed into the GUI.")
    lines.append(f"fit_field_lo = {field_lo:.6g};  % mT")
    lines.append(f"fit_field_hi = {field_hi:.6g};  % mT")
    lines.append("in_range = field_axis >= fit_field_lo & field_axis <= fit_field_hi;")
    lines.append("field_axis = field_axis(in_range);")
    lines.append("data = data(in_range);")
    lines.append("Exp.nPoints = numel(data);")
    lines.append("")
    lines.append("% --- Interactive fit (esfit) ---")
    lines.append("% Vary fields mirror Sys field names; each value is a symmetric +/-")
    lines.append("% delta around the Sys value above (NOT a percentage). Fields not")
    lines.append("% listed are held fixed at their Sys value.")

    n_components = len(vary_specs)
    vary_component_names = []
    any_vary_at_all = False
    for i, vdict in enumerate(vary_specs):
        if not vdict:
            vary_component_names.append("[]")
            continue
        any_vary_at_all = True
        varname = "SysVary" if i == 0 else f"SysVary{i + 1}"
        lines.append(f"{varname} = struct();")
        for field, delta in vdict.items():
            mapped_field = field
            if chili_mode and field in ("tcorr", "Diff"):
                mapped_field = "log" + field
                lines.append(f"% NOTE: chili should be fit via log{field}, not {field}, per")
                lines.append(f"% EasySpin's fitting guide - using {varname}.{mapped_field} instead.")
            lines.append(f"{varname}.{mapped_field} = {delta};")
        vary_component_names.append(varname)

    lines.append("")
    if not any_vary_at_all:
        lines.append("% (no fields checked as 'vary' anywhere - add at least one before fitting)")
    if n_components > 1:
        lines.append("%")
        lines.append("% NOTE: when the esfit GUI opens, it may pre-select only the FIRST")
        lines.append("% component's parameters as active (checked) - the ranges above are")
        lines.append("% correctly defined for every component regardless, but if you want ALL")
        lines.append("% of them varied together, click 'all' in the GUI's Selection toolbar")
        lines.append("% before starting the fit.")

    lines.append(f"FitOpt.Method = '{method} {target}';")
    lines.append("")

    p0_cell = f"{{{sim_sys_expr}, Exp, Opt}}" if has_opt else f"{{{sim_sys_expr}, Exp}}"
    if n_components == 1:
        vary_expr = "{" + vary_component_names[0] + "}"
    else:
        # Per esfit's docs: {{varySys1,varySys2}} - one extra outer cell
        # wrapping the per-component list, matching how Sys0 itself
        # nests as {{Sys1,Sys2},Exp0} for multi-component mixtures.
        vary_expr = "{{" + ", ".join(vary_component_names) + "}}"

    lines.append(f"esfit(data, @{fcn_name}, {p0_cell}, {vary_expr}, FitOpt);")
    lines.append("% No output argument requested above -> opens the interactive")
    lines.append("% fitting GUI (sliders, algorithm/target choice, live comparison).")
    lines.append("% This is the documented trigger: esfit opens the GUI iff no")
    lines.append("% output is requested, regardless of a trailing semicolon.")
    lines.append("%")
    lines.append("% NOTE: if you see 'Unrecognized field name \"Arg\"' or another")
    lines.append("% argument-parsing error here, check that at least one Vary field")
    lines.append("% above is non-empty - an empty SysVary struct can trigger this.")
    lines.append("% To run headlessly instead and capture results directly:")
    lines.append(f"% fit = esfit(data, @{fcn_name}, {p0_cell}, {vary_expr}, FitOpt);")

    post_lines = []
    post_lines.append("%% --- Part 2: run this AFTER exporting from the esfit GUI ---")
    post_lines.append("% This is a separate script from the one that opened the esfit GUI.")
    post_lines.append("% Run it in the SAME MATLAB session (so 'Sys', 'Exp', 'B', 'field_axis'")
    post_lines.append("% and 'data' are still in your workspace from running Part 1) - after")
    post_lines.append("% you've fit interactively and clicked 'Export' in the GUI window.")
    post_lines.append("%")
    post_lines.append("% MATLAB's esfit GUI auto-names each exported result 'fitN', where N")
    post_lines.append("% is the 'ID' column of the GUI's 'Fit parameter sets' table (e.g.")
    post_lines.append("% exporting the row with ID 3 prints \"Fit parameter set 3 assigned")
    post_lines.append("% to variable 'fit3'\" and creates a variable called fit3). EDIT the")
    post_lines.append("% line below to match whichever ID you actually exported:")
    post_lines.append("fit_result = fit1;   % <-- change 'fit1' to match your own export (e.g. fit2, fit3, ...)")
    post_lines.append("")
    post_lines.append("fitted_spc = fit_result.fit;    % simulated spectrum, scaled to match data")
    post_lines.append("% (fit_result.fitraw is the unscaled version straight from the")
    post_lines.append("% simulation function, if you want that instead of the scaled 'fit')")
    post_lines.append("")
    post_lines.append("figure")
    post_lines.append("plot(field_axis, data, 'k'); hold on")
    post_lines.append("plot(B, fitted_spc, 'r'); hold off")
    post_lines.append("xlabel('Magnetic field (mT)')")
    post_lines.append("ylabel('Intensity (a.u.)')")
    post_lines.append("legend('Experimental data', 'Optimized fit')")
    post_lines.append(f"title('{fcn_name}: experimental data vs. optimized fit')")
    post_lines.append("")
    post_lines.append("% --- Build a Supporting Information summary from the FITTED parameters ---")
    post_lines.append("% Reads fit_result.argsfit, which per esfit's own documentation contains")
    post_lines.append("% \"EasySpin-style structures containing the fitted parameters\", mirroring")
    post_lines.append("% the shape of the original p0 argument passed to esfit. For a")
    post_lines.append("% multi-component fit that means argsfit{1} should be a cell array of")
    post_lines.append("% fitted Sys structs and argsfit{2} the fitted Exp struct - this exact")
    post_lines.append("% shape is inferred from the docs, not independently confirmed against a")
    post_lines.append("% real multi-component esfit output, so it's wrapped in try/catch: if it")
    post_lines.append("% doesn't match, the CSV/MAT saves below still work, you'll just see a")
    post_lines.append("% warning instead of the summary/caption.")
    post_lines.append("try")
    post_lines.append("    argsfit = fit_result.argsfit;")
    post_lines.append("    SysFit = argsfit{1};")
    post_lines.append("    if ~iscell(SysFit)")
    post_lines.append("        SysFit = {SysFit};")
    post_lines.append("    end")
    post_lines.append("    ExpFit = argsfit{2};")
    post_lines.append("")
    post_lines.append("    si_lines = {};")
    post_lines.append("    si_lines{end+1} = '% --- Simulation summary (for Supporting Information) ---';")
    post_lines.append(f"    si_lines{{end+1}} = sprintf('%% Function: {fcn_name}');")
    post_lines.append("    si_lines{end+1} = sprintf('%% Number of components: %d', numel(SysFit));")
    post_lines.append("    si_lines{end+1} = sprintf('%% Exp.mwFreq = %.6f GHz', ExpFit.mwFreq);")
    post_lines.append("    if isfield(ExpFit, 'Temperature')")
    post_lines.append("        si_lines{end+1} = sprintf('%% Exp.Temperature = %.4g K', ExpFit.Temperature);")
    post_lines.append("    end")
    post_lines.append("")
    post_lines.append("    comp_descs = {};")
    post_lines.append("    for k = 1:numel(SysFit)")
    post_lines.append("        Sk = SysFit{k};")
    post_lines.append("        wt = 1;")
    post_lines.append("        if isfield(Sk, 'weight'); wt = Sk.weight; end")
    post_lines.append("        nucs = 'none';")
    post_lines.append("        if isfield(Sk, 'Nucs'); nucs = Sk.Nucs; end")
    post_lines.append("        Astr = 'n/a';")
    post_lines.append("        if isfield(Sk, 'A')")
    post_lines.append("            A_mT = unitconvert(Sk.A, 'MHz->mT');")
    post_lines.append("            Astr = mat2str(round(A_mT, 4));")
    post_lines.append("        end")
    post_lines.append("        Sstr = '1/2';")
    post_lines.append("        if isfield(Sk, 'S'); Sstr = mat2str(Sk.S); end")
    post_lines.append("        gstr = mat2str(round(Sk.g, 4));")
    post_lines.append("        % Describe lwpp in mT, distinguishing pure/mixed Gaussian-Lorentzian")
    post_lines.append("        if numel(Sk.lwpp) == 1")
    post_lines.append("            lw_desc = sprintf('%.4g mT (Gaussian)', Sk.lwpp);")
    post_lines.append("        elseif numel(Sk.lwpp) == 2 && Sk.lwpp(1) == 0")
    post_lines.append("            lw_desc = sprintf('%.4g mT (Lorentzian)', Sk.lwpp(2));")
    post_lines.append("        elseif numel(Sk.lwpp) == 2 && Sk.lwpp(2) == 0")
    post_lines.append("            lw_desc = sprintf('%.4g mT (Gaussian)', Sk.lwpp(1));")
    post_lines.append("        elseif numel(Sk.lwpp) == 2")
    post_lines.append("            lw_desc = sprintf('%.4g mT Gaussian + %.4g mT Lorentzian', Sk.lwpp(1), Sk.lwpp(2));")
    post_lines.append("        else")
    post_lines.append("            lw_desc = [mat2str(round(Sk.lwpp, 4)) ' mT'];")
    post_lines.append("        end")
    # v6.5: zero-field splitting, reported as D (cm^-1) and E/D
    post_lines.append("        zfs_str = '';")
    post_lines.append("        if isfield(Sk, 'D') && numel(Sk.D) == 2 && Sk.D(1) ~= 0")
    post_lines.append("            zfs_str = sprintf(', D = %.4g cm^-1, E/D = %.3g', ...")
    post_lines.append("                Sk.D(1)/29979.2458, Sk.D(2)/Sk.D(1));")
    post_lines.append("            if isfield(Sk, 'DStrain')")
    post_lines.append("                zfs_str = [zfs_str sprintf(' (D-strain FWHM = %s cm^-1)', ...")
    post_lines.append("                    mat2str(round(Sk.DStrain/29979.2458, 4)))];")
    post_lines.append("            end")
    post_lines.append("        end")
    post_lines.append("        si_lines{end+1} = sprintf(['%% Component %d: Nucs=%s, weight=%.4g, ' ...")
    post_lines.append("            'S=%s, g=%s%s, lwpp=%s, A=%s mT'], k, nucs, wt, Sstr, gstr, zfs_str, lw_desc, Astr);")
    post_lines.append("        comp_descs{end+1} = sprintf('%s: S = %s, g = %s%s, linewidth = %s, A = %s mT', ...")
    post_lines.append("            nucs, Sstr, gstr, zfs_str, lw_desc, Astr);")
    post_lines.append("    end")
    post_lines.append("")
    post_lines.append("    temp_str = 'the recorded temperature';")
    post_lines.append("    if isfield(ExpFit, 'Temperature')")
    post_lines.append("        temp_str = sprintf('%.4g K', ExpFit.Temperature);")
    post_lines.append("    end")
    post_lines.append("    si_caption = sprintf(['Figure SXX. Solution X-band EPR spectrum of ' ...")
    post_lines.append("        '<COMPOUND> recorded at %s (black) and the corresponding simulation ' ...")
    post_lines.append(f"        '(red) performed using the EasySpin {fcn_name} function (%d ' ...")
    post_lines.append("        'component(s): %s).'], temp_str, numel(SysFit), strjoin(comp_descs, '; '));")
    post_lines.append("    disp(' ');")
    post_lines.append("    disp('--- Suggested Supporting Information caption (edit figure number,');")
    post_lines.append("    disp('    compound name, and temperature/prose as needed) ---');")
    post_lines.append("    disp(si_caption);")
    post_lines.append("    disp(' ');")
    post_lines.append("catch ME")
    post_lines.append("    si_lines = {'% Could not auto-extract fit_result.argsfit - inspect it manually.'};")
    post_lines.append("    warning('SI summary generation failed: %s', ME.message);")
    post_lines.append("end")
    post_lines.append("")
    post_lines.append("% Save the optimized trace as a CSV, with the Supporting Information")
    post_lines.append("% summary above written as '%'-prefixed comment lines before the data")
    post_lines.append("% table (any CSV/text reader that skips non-numeric header lines, and")
    post_lines.append("% this tool's own Residual Builder, will read past them automatically):")
    post_lines.append("[csv_name, csv_dir] = uiputfile('*.csv', 'Save optimized trace as CSV', ...")
    post_lines.append(f"    '{fcn_name}_fit_result.csv');")
    post_lines.append("if isequal(csv_name, 0)")
    post_lines.append("    disp('CSV save cancelled.');")
    post_lines.append("else")
    post_lines.append("    csv_path = fullfile(csv_dir, csv_name);")
    post_lines.append("    fid_out = fopen(csv_path, 'w');")
    post_lines.append("    for i = 1:numel(si_lines)")
    post_lines.append("        fprintf(fid_out, '%s\\n', si_lines{i});")
    post_lines.append("    end")
    post_lines.append("    fprintf(fid_out, 'Field_mT,Experimental,Fit\\n');")
    post_lines.append("    for i = 1:numel(field_axis)")
    post_lines.append("        fprintf(fid_out, '%.6f,%.10g,%.10g\\n', field_axis(i), data(i), fitted_spc(i));")
    post_lines.append("    end")
    post_lines.append("    fclose(fid_out);")
    post_lines.append("end")
    post_lines.append("")
    post_lines.append("% Save the full fit structure (parameters, uncertainties, etc.) -")
    post_lines.append("% also opens a Save As dialog:")
    post_lines.append("[mat_name, mat_dir] = uiputfile('*.mat', 'Save fit structure as .mat', ...")
    post_lines.append(f"    '{fcn_name}_fit_result.mat');")
    post_lines.append("if isequal(mat_name, 0)")
    post_lines.append("    disp('MAT save cancelled.');")
    post_lines.append("else")
    post_lines.append("    save(fullfile(mat_dir, mat_name), 'fit_result');")
    post_lines.append("end")

    return lines, post_lines


# --------------------------------------------------------------------------- #
# Line broadening section.
#
# mode="tensor" (pepper) - lwpp is the always-applied base broadening;
# HStrain and gStrain are OPTIONAL anisotropic-strain mechanisms (a
# spread in tensor values across the sample), each toggled on/off with
# a checkbox that reveals a 3-value tensor entry. These genuinely need
# an anisotropic/powder model to mean anything, AND (confirmed against
# pepper's own doc example) are EasySpin's actual documented mechanism
# for anisotropic/orientation-dependent broadening.
#
# lwpp/lw itself (mode="iso" garlic AND mode="tensor" pepper - this is
# a SHARED spin-system field per garlic.html's docs) only ever accepts
# 1 element (Gaussian) or 2 elements ([Gaussian Lorentzian]) - splitting
# into Gaussian+Lorentzian components is the only extra option either
# model actually supports for lwpp itself. A previous version's
# "Anisotropic linewidth (3 principal values)" option for pepper's lwpp
# was based on a wrong assumption and has been removed (v4.9) - it was
# never a real EasySpin capability; HStrain/gStrain above are what
# actually deliver anisotropic broadening.
#
# A-strain (needs per-nucleus UI, since it modifies a specific nucleus's
# hyperfine spread) and D-strain (needs Sys.D/E fields, not yet exposed
# anywhere in this tool) are intentionally NOT implemented here.
# --------------------------------------------------------------------------- #
class BroadeningSection(ctk.CTkFrame):
    def __init__(self, master, mode="tensor", **kwargs):
        super().__init__(master, fg_color="transparent", **kwargs)
        self.mode = mode

        base_row = ctk.CTkFrame(self, fg_color="transparent")
        base_row.pack(anchor="w", fill="x", pady=2)
        ctk.CTkLabel(base_row, text="Linewidth, lwpp (mT):", width=190,
                     anchor="w").pack(side="left")
        self.lwpp_entry = ctk.CTkEntry(base_row, width=90, placeholder_text="1.0")
        self.lwpp_entry.insert(0, "1.0")
        self.lwpp_entry.pack(side="left", padx=4)
        info_icon(base_row, "Base broadening, always applied - peak-to-peak "
                             "linewidth in mT."
                  ).pack(side="left", padx=(6, 0))

        if mode == "iso":
            gl_row = ctk.CTkFrame(self, fg_color="transparent")
            gl_row.pack(anchor="w", fill="x", pady=(6, 2))
            self.gl_var = ctk.BooleanVar(value=False)
            ctk.CTkCheckBox(gl_row, text="Split into Gaussian + Lorentzian components",
                             variable=self.gl_var, command=self._refresh_gl).pack(side="left")
            info_icon(gl_row, "garlic's Sys.lwpp accepts a 2-element [Gaussian "
                               "Lorentzian] vector instead of one scalar, for a "
                               "mixed lineshape. HStrain/gStrain (anisotropic tensor "
                               "spread) don't apply to garlic's isotropic fast-motion "
                               "model - this is the broadening option it actually "
                               "supports."
                      ).pack(side="left", padx=(6, 0))
            self.gl_frame = ctk.CTkFrame(self, fg_color="transparent", height=1)
            self.gl_frame.pack(anchor="w", fill="x")
            self.gauss_entry = None
            self.lorentz_entry = None
        else:
            # v4.9: the "Anisotropic linewidth (3 principal values)"
            # checkbox that used to live here has been REMOVED - confirmed
            # against official EasySpin docs (garlic.html, which defines
            # the shared Sys.lwpp/lw field also used by pepper) that
            # Sys.lwpp/lw only ever accepts 1 or 2 elements ([Gaussian] or
            # [Gaussian Lorentzian]), never a 3-element anisotropic form.
            # pepper's own doc example achieves anisotropic broadening via
            # Sys.HStrain instead (see the HStrain checkbox below) - that
            # was already the intended real mechanism for this, so nothing
            # is lost, only a wrong extra option is gone.
            gl_row = ctk.CTkFrame(self, fg_color="transparent")
            gl_row.pack(anchor="w", fill="x", pady=(6, 2))
            self.gl_var = ctk.BooleanVar(value=False)
            ctk.CTkCheckBox(gl_row, text="Split into Gaussian + Lorentzian components",
                             variable=self.gl_var, command=self._refresh_lw_extra).pack(side="left")
            info_icon(gl_row, "Mixed lineshape instead of a pure Gaussian - "
                               "Sys.lwpp = [Gaussian Lorentzian], per EasySpin's "
                               "documented Sys.lwpp/lw field (confirmed 1- or 2-"
                               "element only, shared by garlic and pepper). For "
                               "orientation-dependent (anisotropic) broadening, use "
                               "HStrain/gStrain below instead - that's EasySpin's "
                               "actual mechanism for it, per pepper's own doc example."
                      ).pack(side="left", padx=(6, 0))

            self.lw_extra_frame = ctk.CTkFrame(self, fg_color="transparent", height=1)
            self.lw_extra_frame.pack(anchor="w", fill="x")
            self.lw_extra_entries = []  # flat list; single (Gaussian, Lorentzian) pair when split is on

            h_row = ctk.CTkFrame(self, fg_color="transparent")
            h_row.pack(anchor="w", fill="x", pady=(8, 2))
            self.hstrain_var = ctk.BooleanVar(value=False)
            ctk.CTkCheckBox(h_row, text="Add field (H) strain - Sys.HStrain",
                             variable=self.hstrain_var, command=self._refresh_hstrain).pack(side="left")
            info_icon(h_row, "Gaussian broadening from unresolved environmental "
                              "heterogeneity, directly in field units (mT). Three "
                              "principal values, same axes as g."
                      ).pack(side="left", padx=(6, 0))
            self.hstrain_frame = ctk.CTkFrame(self, fg_color="transparent", height=1)
            self.hstrain_frame.pack(anchor="w", fill="x")
            self.hstrain_entries = []

            g_row = ctk.CTkFrame(self, fg_color="transparent")
            g_row.pack(anchor="w", fill="x", pady=2)
            self.gstrain_var = ctk.BooleanVar(value=False)
            ctk.CTkCheckBox(g_row, text="Add g-strain - Sys.gStrain",
                             variable=self.gstrain_var, command=self._refresh_gstrain).pack(side="left")
            info_icon(g_row, "Gaussian spread in the principal g-values themselves "
                              "(different molecules sitting in slightly different g "
                              "environments). Three principal values matching g's "
                              "axes."
                      ).pack(side="left", padx=(6, 0))
            self.gstrain_frame = ctk.CTkFrame(self, fg_color="transparent", height=1)
            self.gstrain_frame.pack(anchor="w", fill="x")
            self.gstrain_entries = []

            info_note(self, "A-strain isn't available here yet - it needs per-"
                             "nucleus UI (it modifies a specific nucleus's "
                             "hyperfine spread). D-strain (v6.5) lives in the "
                             "Zero-Field Splitting section above, since it only "
                             "means anything when Sys.D is on.")

        if mode == "iso":
            self._refresh_gl()
        else:
            self._refresh_lw_extra()
            self._refresh_hstrain()
            self._refresh_gstrain()

    def _refresh_lw_extra(self):
        # Rebuilds the Gaussian/Lorentzian split fields for pepper mode.
        # v4.9: simplified back down to a single (Gaussian, Lorentzian)
        # pair now that the invalid "anisotropic linewidth" option is gone
        # - independent from HStrain/gStrain's own refreshes, so toggling
        # one never wipes the other.
        for w in self.lw_extra_frame.winfo_children():
            w.destroy()
        self.lw_extra_entries = []

        if not self.gl_var.get():
            return

        base = self.lwpp_entry.get().strip() or "1.0"
        row = ctk.CTkFrame(self.lw_extra_frame, fg_color="transparent")
        row.pack(anchor="w", pady=2)
        ctk.CTkLabel(row, text="Gaussian (mT):", width=110).pack(side="left")
        g_entry = ctk.CTkEntry(row, width=70)
        g_entry.insert(0, base)
        g_entry.pack(side="left", padx=4)
        ctk.CTkLabel(row, text="Lorentzian (mT):", width=110).pack(
            side="left", padx=(10, 0))
        l_entry = ctk.CTkEntry(row, width=70)
        l_entry.insert(0, "0")
        l_entry.pack(side="left", padx=4)
        self.lw_extra_entries.append((g_entry, l_entry))

    def _refresh_gl(self):
        for w in self.gl_frame.winfo_children():
            w.destroy()
        if self.gl_var.get():
            row = ctk.CTkFrame(self.gl_frame, fg_color="transparent")
            row.pack(anchor="w", pady=2)
            ctk.CTkLabel(row, text="Gaussian (mT):", width=110).pack(side="left")
            self.gauss_entry = ctk.CTkEntry(row, width=70)
            self.gauss_entry.insert(0, self.lwpp_entry.get().strip() or "1.0")
            self.gauss_entry.pack(side="left", padx=4)
            ctk.CTkLabel(row, text="Lorentzian (mT):", width=110).pack(
                side="left", padx=(10, 0))
            self.lorentz_entry = ctk.CTkEntry(row, width=70)
            self.lorentz_entry.insert(0, "0")
            self.lorentz_entry.pack(side="left", padx=4)
        else:
            self.gauss_entry = None
            self.lorentz_entry = None

    def _refresh_hstrain(self):
        # Only rebuilds hstrain_frame - toggling gStrain never touches
        # this, so HStrain values entered here survive a gStrain toggle.
        for w in self.hstrain_frame.winfo_children():
            w.destroy()
        if self.hstrain_var.get():
            row = ctk.CTkFrame(self.hstrain_frame, fg_color="transparent")
            row.pack(anchor="w", pady=2)
            self.hstrain_entries = [ctk.CTkEntry(row, width=70) for _ in range(3)]
            for lbl, e in zip(("Hx", "Hy", "Hz"), self.hstrain_entries):
                ctk.CTkLabel(row, text=f"{lbl} (mT):", width=65).pack(side="left")
                e.insert(0, "0.5")
                e.pack(side="left", padx=4)
        else:
            self.hstrain_entries = []

    def _refresh_gstrain(self):
        # Only rebuilds gstrain_frame - toggling HStrain never touches
        # this, so gStrain values entered here survive an HStrain toggle.
        for w in self.gstrain_frame.winfo_children():
            w.destroy()
        if self.gstrain_var.get():
            row = ctk.CTkFrame(self.gstrain_frame, fg_color="transparent")
            row.pack(anchor="w", pady=2)
            self.gstrain_entries = [ctk.CTkEntry(row, width=70) for _ in range(3)]
            for lbl, e in zip(("gx", "gy", "gz"), self.gstrain_entries):
                ctk.CTkLabel(row, text=f"{lbl}:", width=45).pack(side="left")
                e.insert(0, "0.01")
                e.pack(side="left", padx=4)
        else:
            self.gstrain_entries = []

    def get_lwpp(self):
        return self.lwpp_entry.get().strip()

    def get_lwpp_full(self):
        """Sys.lwpp expression:
          garlic (mode="iso") and pepper (mode="tensor") are now IDENTICAL
          here: plain scalar, or "[Gaussian Lorentzian]" if the G+L split
          is checked. Confirmed against garlic.html's documented Sys.lwpp/
          lw field (shared by garlic and pepper): only 1- or 2-element
          forms are valid - no anisotropic 3-element form exists. (v4.9:
          the old pepper-only "anisotropic linewidth" branch was removed
          entirely - see the v4.9 changelog.)"""
        base = self.get_lwpp()
        if self.mode == "iso":
            if self.gl_var.get() and self.lorentz_entry is not None:
                gauss = self.gauss_entry.get().strip()
                lorentz = self.lorentz_entry.get().strip()
                return f"[{gauss} {lorentz}]"
            return base

        if self.gl_var.get() and self.lw_extra_entries:
            g_entry, l_entry = self.lw_extra_entries[0]
            return f"[{g_entry.get().strip()} {l_entry.get().strip()}]"
        return base

    def is_gl_split_enabled(self):
        return self.gl_var.get()

    def get_lw_extra_raw_values(self):
        """Flat list of every raw string value currently in the dynamic
        lwpp-extension fields (for validation - checking none are blank),
        regardless of which combination is active."""
        out = []
        for item in self.lw_extra_entries:
            if isinstance(item, tuple):
                out.extend(e.get().strip() for e in item)
            else:
                out.append(item.get().strip())
        return out

    def is_hstrain_enabled(self):
        return self.mode == "tensor" and self.hstrain_var.get()

    def get_hstrain_values(self):
        return [e.get().strip() for e in self.hstrain_entries]

    def is_gstrain_enabled(self):
        return self.mode == "tensor" and self.gstrain_var.get()

    def get_gstrain_values(self):
        return [e.get().strip() for e in self.gstrain_entries]


# --------------------------------------------------------------------------- #
# v6.5: Zero-field splitting (pepper only).
#
# Verified against EasySpin's own docs (spinsystem.html, broadenings.html,
# garlic.html - retrieved 2026-10-05):
#   - Sys.D is in MHz ("1 cm-1 = 29979 MHz") and accepts a scalar D, a
#     2-element [D E], 3 principal values, or a full 3x3 matrix. This tool
#     always writes the 2-element [D E] form, so esfit sees exactly two
#     ZFS parameters (D and E) - same reasoning as the v5.4 axial-g fix:
#     no redundant array slot for esfit to split apart.
#   - Sys.DStrain = FWHM_D or [FWHM_D FWHM_E], Gaussian widths in MHz.
#   - garlic is restricted to S=1/2 doublets, so ZFS is NOT offered there.
#
# The GUI takes D in a chosen unit (cm^-1 / GHz / MHz) and E/D (unitless,
# conventionally 0 <= E/D <= 1/3) and the generated script converts to MHz
# explicitly with a numeric factor (29979.2458 MHz per cm^-1, from c),
# rather than relying on a unitconvert() unit string.
#
# Fitting: esfit can only vary the actual Sys.D array elements, i.e. D
# and E (in MHz) - not the E/D ratio itself. The E/D delta typed into the
# Fit section is converted to an E delta as  dE = d(E/D) * |D_start|.
# The SI summary in Part 2 converts the fitted [D E] back to D (cm^-1)
# and E/D.
# --------------------------------------------------------------------------- #
CM1_TO_MHZ = 29979.2458  # 1 cm^-1 in MHz (c = 299792458 m/s)


def _parse_spin_value(s):
    """'1', '1.5', '3/2', '[1 1/2]' (first spin) -> float, or None."""
    inner = s.strip().strip("[]").strip()
    if not inner:
        return None
    tok = inner.split()[0].rstrip(",;")
    try:
        if "/" in tok:
            num, den = tok.split("/", 1)
            return float(num) / float(den)
        return float(tok)
    except (ValueError, ZeroDivisionError):
        return None
ZFS_UNIT_TO_MHZ = {"cm^-1": CM1_TO_MHZ, "GHz": 1000.0, "MHz": 1.0}


class ZFSSection(ctk.CTkFrame):
    """Sys.D (as D + E/D) and optional Sys.DStrain for the Pepper tab.
    All widgets are created ONCE at construction and only shown/hidden
    via pack/pack_forget, so the tab's one-time _wire_live_regenerate()
    walk hooks every one of them."""

    def __init__(self, master, **kwargs):
        super().__init__(master, fg_color="transparent", **kwargs)

        top = ctk.CTkFrame(self, fg_color="transparent")
        top.pack(anchor="w", fill="x", pady=2)
        self.enabled_var = ctk.BooleanVar(value=False)
        ctk.CTkCheckBox(top, text="Include zero-field splitting (Sys.D)",
                         variable=self.enabled_var, command=self._refresh).pack(side="left")
        info_icon(top, "Zero-field splitting for S >= 1 systems. Enter the axial "
                        "parameter D and the rhombicity E/D (0 = axial, 1/3 = "
                        "fully rhombic). Written as Sys.D = [D E] in MHz, the "
                        "documented 2-element form. Note: the SIGN of D only "
                        "affects a powder spectrum through Boltzmann populations "
                        "- set Exp.Temperature if you want the sign to matter."
                  ).pack(side="left", padx=(6, 0))

        self.body = ctk.CTkFrame(self, fg_color="transparent")
        # (packed/unpacked by _refresh)

        d_row = ctk.CTkFrame(self.body, fg_color="transparent")
        d_row.pack(anchor="w", fill="x", pady=2)
        ctk.CTkLabel(d_row, text="D:", width=190, anchor="w").pack(side="left")
        self.d_entry = ctk.CTkEntry(d_row, width=90, placeholder_text="e.g. 0.1")
        self.d_entry.insert(0, "0.1")
        self.d_entry.pack(side="left", padx=4)
        self.unit_var = ctk.StringVar(value="cm^-1")
        ctk.CTkOptionMenu(d_row, variable=self.unit_var, width=90,
                           values=list(ZFS_UNIT_TO_MHZ.keys()),
                           command=lambda _v: self._refresh()).pack(side="left", padx=(4, 0))

        e_row = ctk.CTkFrame(self.body, fg_color="transparent")
        e_row.pack(anchor="w", fill="x", pady=2)
        ctk.CTkLabel(e_row, text="E/D (0 to 1/3):", width=190, anchor="w").pack(side="left")
        self.eod_entry = ctk.CTkEntry(e_row, width=90, placeholder_text="e.g. 0.05")
        self.eod_entry.insert(0, "0")
        self.eod_entry.pack(side="left", padx=4)
        info_icon(e_row, "Rhombicity. E is computed as (E/D) x D, so E carries the "
                          "same sign as D. Convention keeps 0 <= E/D <= 1/3; values "
                          "outside that range describe the same tensor with relabeled "
                          "axes, so they're rejected to keep fits unambiguous."
                  ).pack(side="left", padx=(6, 0))

        ds_row = ctk.CTkFrame(self.body, fg_color="transparent")
        ds_row.pack(anchor="w", fill="x", pady=(6, 2))
        self.dstrain_var = ctk.BooleanVar(value=False)
        ctk.CTkCheckBox(ds_row, text="Add D-strain - Sys.DStrain",
                         variable=self.dstrain_var, command=self._refresh).pack(side="left")
        info_icon(ds_row, "Gaussian distributions (FWHM) of D and E across the "
                           "sample - per broadenings.html, Sys.DStrain = [FWHM_D "
                           "FWHM_E] in MHz. Entered here in the same unit as D "
                           "above and converted for you. Often the dominant "
                           "broadening for high-spin frozen solutions."
                   ).pack(side="left", padx=(6, 0))
        self.ds_fields = ctk.CTkFrame(self.body, fg_color="transparent")
        ctk.CTkLabel(self.ds_fields, text="FWHM_D:", width=70).pack(side="left")
        self.ds_d_entry = ctk.CTkEntry(self.ds_fields, width=70)
        self.ds_d_entry.insert(0, "0.01")
        self.ds_d_entry.pack(side="left", padx=4)
        ctk.CTkLabel(self.ds_fields, text="FWHM_E:", width=70).pack(side="left", padx=(10, 0))
        self.ds_e_entry = ctk.CTkEntry(self.ds_fields, width=70)
        self.ds_e_entry.insert(0, "0")
        self.ds_e_entry.pack(side="left", padx=4)
        self.ds_unit_label = ctk.CTkLabel(self.ds_fields, text="(cm^-1)", text_color="gray60")
        self.ds_unit_label.pack(side="left", padx=(4, 0))

        self.note = ctk.CTkLabel(
            self.body, font=ctk.CTkFont(size=11), text_color="gray60",
            wraplength=380, justify="left", anchor="w",
            text="Requires S >= 1. For S = 1/2 there is no zero-field splitting "
                 "and the script won't generate with this on.")
        self.note.pack(anchor="w", pady=(2, 6))

        self._refresh()

    def _refresh(self):
        if self.enabled_var.get():
            self.body.pack(anchor="w", fill="x")
        else:
            self.body.pack_forget()
        if self.dstrain_var.get():
            self.ds_fields.pack(anchor="w", pady=2, before=self.note)
        else:
            self.ds_fields.pack_forget()
        self.ds_unit_label.configure(text=f"({self.unit_var.get()})")

    # --- accessors ---------------------------------------------------- #
    def is_enabled(self):
        return self.enabled_var.get()

    def is_dstrain_enabled(self):
        return self.enabled_var.get() and self.dstrain_var.get()

    def unit_factor(self):
        return ZFS_UNIT_TO_MHZ.get(self.unit_var.get(), CM1_TO_MHZ)

    def get_values(self):
        """Validated (D, E/D, unit) as floats/str. Raises ValueError."""
        d_raw = self.d_entry.get().strip()
        eod_raw = self.eod_entry.get().strip()
        if not d_raw or not eod_raw:
            raise ValueError("Zero-field splitting is on but D or E/D is empty.")
        try:
            d = float(d_raw)
            eod = float(eod_raw)
        except ValueError:
            raise ValueError("Zero-field splitting: D and E/D must be numbers.")
        if not (0.0 <= eod <= 1.0 / 3.0 + 1e-9):
            raise ValueError(f"E/D = {eod:g} is outside the conventional 0 to 1/3 "
                              "range. (A value outside it describes the same tensor "
                              "with relabeled axes - pick the equivalent D, E/D "
                              "inside the range.)")
        return d, eod, self.unit_var.get()

    def get_dstrain_values(self):
        vals = [self.ds_d_entry.get().strip(), self.ds_e_entry.get().strip()]
        if any(v == "" for v in vals):
            raise ValueError("D-strain is enabled but FWHM_D/FWHM_E aren't both filled in.")
        try:
            return [float(v) for v in vals]
        except ValueError:
            raise ValueError("D-strain FWHM values must be numbers.")

    def collect_state(self):
        return {
            "enabled": self.enabled_var.get(),
            "D": self.d_entry.get().strip(),
            "EoverD": self.eod_entry.get().strip(),
            "unit": self.unit_var.get(),
            "dstrain_enabled": self.dstrain_var.get(),
            "dstrain": [self.ds_d_entry.get().strip(), self.ds_e_entry.get().strip()],
        }

    def apply_state(self, data):
        data = data or {}
        self.enabled_var.set(bool(data.get("enabled", False)))
        self.d_entry.delete(0, "end"); self.d_entry.insert(0, data.get("D", "0.1"))
        self.eod_entry.delete(0, "end"); self.eod_entry.insert(0, data.get("EoverD", "0"))
        unit = data.get("unit", "cm^-1")
        self.unit_var.set(unit if unit in ZFS_UNIT_TO_MHZ else "cm^-1")
        self.dstrain_var.set(bool(data.get("dstrain_enabled", False)))
        ds = data.get("dstrain", ["0.01", "0"])
        self.ds_d_entry.delete(0, "end"); self.ds_d_entry.insert(0, ds[0] if ds else "0.01")
        self.ds_e_entry.delete(0, "end"); self.ds_e_entry.insert(0, ds[1] if len(ds) > 1 else "0")
        self._refresh()


##############################################################################
# electron_helper.py content (spin-state suggestion dialog)
##############################################################################

import customtkinter as ctk


class ElectronCountDialog(ctk.CTkToplevel):
    def __init__(self, master, on_select=None):
        super().__init__(master)
        self.title("Electron count / spin-state helper")
        self.geometry("560x620")
        self.minsize(480, 480)
        self.on_select = on_select

        ctk.CTkLabel(self, text="Element (e.g. Fe, Ti, U):").pack(anchor="w", padx=16, pady=(14, 2))
        self.element_entry = ctk.CTkEntry(self, width=100)
        self.element_entry.pack(anchor="w", padx=16)
        self.element_entry.bind("<KeyRelease>", lambda e: self.suggest())

        ctk.CTkLabel(self, text="Oxidation state (e.g. 2 for +2):").pack(anchor="w", padx=16, pady=(10, 2))
        self.oxstate_entry = ctk.CTkEntry(self, width=100)
        self.oxstate_entry.pack(anchor="w", padx=16)
        self.oxstate_entry.bind("<KeyRelease>", lambda e: self.suggest())

        ctk.CTkLabel(self, text="Geometry (d-block only):").pack(anchor="w", padx=16, pady=(10, 2))
        self.geometry_var = ctk.StringVar(value="any")
        ctk.CTkOptionMenu(self, variable=self.geometry_var, values=GEOMETRY_OPTIONS,
                           command=lambda _: self.suggest()).pack(anchor="w", padx=16)
        ctk.CTkLabel(
            self, text="'any' shows every possibility across all geometries below - "
                       "narrow it down once you know your compound's coordination geometry.",
            font=ctk.CTkFont(size=11), text_color="gray60", wraplength=500, justify="left",
        ).pack(anchor="w", padx=16, pady=(2, 10))

        self.results_frame = ctk.CTkScrollableFrame(self, width=500, height=300)
        self.results_frame.pack(fill="both", expand=True, padx=16, pady=(0, 10))

        self.note_label = ctk.CTkLabel(
            self, text="", text_color="gray60", wraplength=520, justify="left"
        )
        self.note_label.pack(anchor="w", padx=16, pady=(0, 14))

    def _clear_results(self):
        for w in self.results_frame.winfo_children():
            w.destroy()

    def suggest(self):
        self._clear_results()
        self.note_label.configure(text="")

        element = self.element_entry.get().strip()
        ox_str = self.oxstate_entry.get().strip()
        geometry = self.geometry_var.get()

        if not element or not ox_str:
            return  # quietly wait for both fields while the user is typing

        try:
            oxidation_state = int(ox_str)
        except ValueError:
            return  # quietly wait - e.g. user has only typed "-" so far

        if is_d_block(element):
            suggestions = spin_state_suggestions(element, oxidation_state, geometry)
            if not suggestions:
                self.note_label.configure(
                    text=f"No suggestions available for {element} in oxidation state "
                         f"{oxidation_state} (resulting d-count out of range)."
                )
                return
            for S, label in suggestions:
                self._add_suggestion_button(S, label)

            caveat = GEOMETRY_CAVEATS.get(geometry)
            if caveat:
                self.note_label.configure(text=caveat)

        elif is_f_block(element):
            hint = f_block_hint(element)
            self.note_label.configure(text=hint)
            # Still offer the naive Hund's-rule S as a selectable (but
            # clearly flagged) candidate, since it may be a reasonable
            # starting point even though it needs verification.
            n_f = F_BLOCK_BASE_ELECTRONS.get(element)
            if n_f is not None:
                S = n_f / 2 if n_f <= 7 else (14 - n_f) / 2
                self._add_suggestion_button(
                    S, f"Hund's-rule hint only - S={S:g} (VERIFY before using)"
                )

        else:
            self.note_label.configure(
                text=f"'{element}' not recognized as a d-block or f-block element in "
                     f"this tool's table. Enter S manually."
            )

    def _add_suggestion_button(self, S, label):
        row = ctk.CTkFrame(self.results_frame, fg_color="transparent")
        row.pack(anchor="w", fill="x", pady=4)
        ctk.CTkButton(row, text=f"Use S={S:g}", width=100,
                      command=lambda: self._select(S)).pack(side="left", padx=(0, 10))
        ctk.CTkLabel(row, text=label, wraplength=360, justify="left", anchor="w").pack(side="left")

    def _select(self, S):
        if self.on_select:
            self.on_select(S)
        self.destroy()


##############################################################################
# tab_garlic.py content
##############################################################################

import datetime

import customtkinter as ctk
from tkinter import filedialog, messagebox



class GarlicTab(ctk.CTkFrame):
    def __init__(self, master):
        super().__init__(master, fg_color="transparent")

        left = ctk.CTkScrollableFrame(self, width=440)
        left.pack(side="left", fill="both", padx=(0, 10), pady=10)

        right = ctk.CTkFrame(self, fg_color="transparent")  # v5.1: no outer scroll
        right.pack(side="left", fill="both", expand=True, pady=10)

        # --- Dataset panel ---
        self.dataset_panel = DatasetPanel(left, on_loaded=self._on_dataset_loaded)
        self.dataset_panel.pack(anchor="w", fill="x", pady=(0, 4))

        preset_row = ctk.CTkFrame(left, fg_color="transparent")
        preset_row.pack(anchor="w", fill="x", pady=(0, 10))
        ctk.CTkButton(preset_row, text="Save preset", command=self.save_preset).pack(
            side="left", padx=(0, 8))
        ctk.CTkButton(preset_row, text="Load preset", command=self.load_preset).pack(side="left")

        # --- Spin system ---
        section_title(left, "Spin System (Sys)")
        info_note(left, "garlic is restricted to S=1/2 doublet radicals by EasySpin "
                         "itself - Sys.S is not exposed here.")
        self.g_entry, _ = labeled_entry(
            left, "Isotropic g-value:", default="2.0023",
            tooltip="The electron g-factor. Free electron = 2.0023. Organic "
                    "radicals typically 2.00-2.01. Verify against literature "
                    "for your specific radical/species."
        )
        section_title(left, "Line Broadening")
        self.broadening = BroadeningSection(left, mode="iso")
        self.broadening.pack(anchor="w", fill="x")
        self.lwpp_entry = self.broadening.lwpp_entry  # alias for existing code below
        self.lwpp_entry.delete(0, "end")
        self.lwpp_entry.insert(0, "0.3")  # garlic-typical default (sharper lines than pepper)
        self.tcorr_entry, _ = labeled_entry(
            left, "Rotational corr. time, tcorr (s):", default="",
            placeholder="optional, e.g. 1e-11",
            tooltip="Rotational correlation time in seconds - only needed if "
                    "you're seeing motional line-broadening. Small molecules in "
                    "low-viscosity solvents: ~1e-11 to 1e-10 s. Above roughly "
                    "1e-9 s, garlic's fast-motion approximation may become "
                    "inaccurate - consider chili instead."
        )

        nuc_header = ctk.CTkFrame(left, fg_color="transparent")
        nuc_header.pack(anchor="w", fill="x", pady=(10, 2))
        ctk.CTkLabel(nuc_header, text="Nuclei (hyperfine A entered in mT, isotropic per nucleus):",
                     anchor="w").pack(side="left")
        info_icon(nuc_header, "Isotope label (e.g. 14N, 1H, 55Mn), how many equivalent "
                               "nuclei of that isotope (n), and the isotropic hyperfine "
                               "coupling A in mT. Typical 1H: 0.05-0.2 mT. Typical 14N "
                               "(e.g. nitroxide radicals): 0.5-1.6 mT. Always check "
                               "against literature for your specific radical."
                  ).pack(side="left", padx=(6, 0))
        self.nuc_list = NucleusList(left, mode="iso", on_change=self._on_nuclei_changed)
        self.nuc_list.pack(anchor="w", fill="x")
        info_note(left, "If you enter 2+ isotopes of the SAME element as separate rows "
                         "(e.g. 47Ti and 49Ti), their weights auto-fill from natural "
                         "abundance and the remaining I=0 background is added "
                         "automatically - no separate setup needed.")

        # --- Experimental parameters ---
        section_title(left, "Experiment (Exp)")
        self.mwfreq_entry, _ = labeled_entry(
            left, "Microwave freq (GHz):", default="9.5",
            tooltip="X-band EPR (most common benchtop instruments): ~9.2-9.9 "
                    "GHz. Q-band: ~34 GHz. S-band: ~3-4 GHz. Use the value from "
                    "your .DSC file (MWFQ) for real data."
        )
        self.center_entry, _ = labeled_entry(
            left, "Center field (mT):", default="339",
            tooltip="Center of the field sweep, in mT. For X-band (9.5 GHz) "
                    "with g~2, resonance falls near 339 mT."
        )
        self.sweep_entry, _ = labeled_entry(
            left, "Sweep width (mT):", default="10",
            tooltip="Total field range swept, in mT. Narrow-line radicals: "
                    "5-20 mT is common. Wider for more complex hyperfine "
                    "patterns."
        )
        self.npoints_entry, _ = labeled_entry(
            left, "Number of points:", default="1024",
            tooltip="Number of data points across the sweep. 1024-4096 is "
                    "typical; higher gives finer resolution but simulates "
                    "slower."
        )
        self.harmonic_entry, _ = labeled_entry(
            left, "Harmonic:", default="1",
            tooltip="Detection harmonic. 1 = first-derivative (standard cw "
                    "EPR display, almost always what you want). 0 = absorption "
                    "(rare). 2 = second-derivative (rare, resolution "
                    "enhancement)."
        )
        self.modamp_entry, _ = labeled_entry(
            left, "Modulation amplitude (mT):", default="0.1",
            tooltip="Field modulation amplitude in mT. Should generally stay "
                    "well below the linewidth (roughly a third or less) to "
                    "avoid over-modulation broadening/distortion of the "
                    "lineshape."
        )
        self.temp_entry, _ = labeled_entry(
            left, "Temperature (K, optional):", default="",
            placeholder="e.g. 298",
            tooltip="Sample temperature. For a simple S=1/2 system fit with "
                    "AutoScale='lsq', this likely won't meaningfully change "
                    "the simulated lineshape (any overall intensity effect "
                    "gets absorbed into the autoscaling) - included mainly "
                    "because it's standard to report for a solution EPR "
                    "measurement, and gets added to the script's Supporting "
                    "Information summary."
        )

        # v5.2: no more "Generate Script" button - the script now updates
        # live as soon as you change anything (see _wire_live_regenerate),
        # so a manual trigger isn't needed. This status line shows live
        # validation instead: green when the current inputs produce a
        # valid script, red with the message when they don't (auto ticks
        # never pop a blocking dialog).
        self.gen_status = ctk.CTkLabel(
            left, text="", text_color="#5DCAA5", font=ctk.CTkFont(size=11),
            wraplength=400, justify="left", anchor="w")
        self.gen_status.pack(anchor="w", fill="x", pady=(12, 8))


        # --- Fit (esfit) ---
        section_title(left, "Fit to experimental data (esfit)")
        self.fit_section = FitSection(left, fields=[
            ("g", "g-value", "1", ""),
            ("lwpp", "Linewidth", "1", "mT"),
        ])
        self.fit_section.pack(anchor="w", fill="x")

        # --- Output: plot preview + script (v5.1: side-by-side script
        # panels, non-scrolling right pane, so the plot and both script
        # parts are visible together without having to scroll down) ---
        section_title(right, "Preview")
        self.plot_preview = PlotPreview(right, height=440, on_crop_change=self._on_plot_crop)
        self.plot_preview.pack(fill="x", padx=10, pady=(0, 6))

        code_row = ctk.CTkFrame(right, fg_color="transparent")
        code_row.pack(fill="both", expand=True, padx=10, pady=(0, 10))

        part1_col = ctk.CTkFrame(code_row, fg_color="transparent")
        part1_col.pack(side="left", fill="both", expand=True, padx=(0, 5))
        part2_col = ctk.CTkFrame(code_row, fg_color="transparent")
        part2_col.pack(side="left", fill="both", expand=True, padx=(5, 0))

        section_title(part1_col, "Part 1 (run first)")
        self.output_box = ctk.CTkTextbox(part1_col, font=("Courier", 11), height=170)
        self.output_box.pack(fill="both", expand=True, pady=(0, 6))
        btn_row = ctk.CTkFrame(part1_col, fg_color="transparent")
        btn_row.pack(fill="x")
        ctk.CTkButton(btn_row, text="Copy to Clipboard", command=self.copy_to_clipboard).pack(
            side="left", padx=(0, 8))
        ctk.CTkButton(btn_row, text="Save As .m ...", command=self.save_as).pack(side="left")

        section_title(part2_col, "Part 2 (after esfit export)")
        self.output_box2 = ctk.CTkTextbox(part2_col, font=("Courier", 11), height=170)
        self.output_box2.pack(fill="both", expand=True, pady=(0, 6))
        btn_row2 = ctk.CTkFrame(part2_col, fg_color="transparent")
        btn_row2.pack(fill="x")
        ctk.CTkButton(btn_row2, text="Copy to Clipboard", command=self.copy_to_clipboard2).pack(
            side="left", padx=(0, 8))
        ctk.CTkButton(btn_row2, text="Save As .m ...", command=self.save_as2).pack(side="left")

        self._wire_live_regenerate(self)

    # ------------------------------------------------------------------ #
    # Plot preview
    # ------------------------------------------------------------------ #
    def _plot_experimental(self):
        panel = self.dataset_panel
        if not panel.has_experimental_trace():
            return
        self.plot_preview.plot_experimental(panel.field_mT, panel.intensity)

    # ------------------------------------------------------------------ #
    # Dataset auto-fill
    # ------------------------------------------------------------------ #
    def _on_dataset_loaded(self, exp_values, compound_id, folder):
        if "mwFreq_GHz" in exp_values:
            self.mwfreq_entry.delete(0, "end")
            self.mwfreq_entry.insert(0, f"{exp_values['mwFreq_GHz']:.6f}")
        if "center_G" in exp_values:
            self.center_entry.delete(0, "end")
            self.center_entry.insert(0, f"{exp_values['center_G'] / 10:.4f}")
        if "sweep_G" in exp_values:
            self.sweep_entry.delete(0, "end")
            self.sweep_entry.insert(0, f"{exp_values['sweep_G'] / 10:.4f}")
        if "npoints" in exp_values:
            self.npoints_entry.delete(0, "end")
            self.npoints_entry.insert(0, str(exp_values["npoints"]))
        if "harmonic" in exp_values:
            self.harmonic_entry.delete(0, "end")
            self.harmonic_entry.insert(0, str(exp_values["harmonic"]))
        if "modamp_G" in exp_values:
            self.modamp_entry.delete(0, "end")
            self.modamp_entry.insert(0, f"{exp_values['modamp_G'] / 10:.4f}")

        self._plot_experimental()
        self._debounced_generate()

    # ------------------------------------------------------------------ #
    # Preset save/load
    # ------------------------------------------------------------------ #
    def _collect_state(self):
        return {
            "g": self.g_entry.get().strip(),
            "lwpp": self.lwpp_entry.get().strip(),
            "tcorr": self.tcorr_entry.get().strip(),
            "nuclei": [
                {"isotope": iso, "n": n, "a_mT": a, "weight": weight,
                 "g_override": g_ov, "lwpp_override": lwpp_ov}
                for (iso, n, a, weight), (g_ov, lwpp_ov) in self.nuc_list.collect_with_overrides()
            ],
            "mwFreq": self.mwfreq_entry.get().strip(),
            "center": self.center_entry.get().strip(),
            "sweep": self.sweep_entry.get().strip(),
            "npoints": self.npoints_entry.get().strip(),
            "harmonic": self.harmonic_entry.get().strip(),
            "modamp": self.modamp_entry.get().strip(),
            "temperature": self.temp_entry.get().strip(),
        }

    def _apply_state(self, data):
        self.g_entry.delete(0, "end"); self.g_entry.insert(0, data.get("g", ""))
        self.lwpp_entry.delete(0, "end"); self.lwpp_entry.insert(0, data.get("lwpp", ""))
        self.tcorr_entry.delete(0, "end"); self.tcorr_entry.insert(0, data.get("tcorr", ""))
        self.nuc_list.load(data.get("nuclei", []))
        self.mwfreq_entry.delete(0, "end"); self.mwfreq_entry.insert(0, data.get("mwFreq", ""))
        self.center_entry.delete(0, "end"); self.center_entry.insert(0, data.get("center", ""))
        self.sweep_entry.delete(0, "end"); self.sweep_entry.insert(0, data.get("sweep", ""))
        self.npoints_entry.delete(0, "end"); self.npoints_entry.insert(0, data.get("npoints", ""))
        self.harmonic_entry.delete(0, "end"); self.harmonic_entry.insert(0, data.get("harmonic", ""))
        self.modamp_entry.delete(0, "end"); self.modamp_entry.insert(0, data.get("modamp", ""))
        self.temp_entry.delete(0, "end"); self.temp_entry.insert(0, data.get("temperature", ""))
        self._debounced_generate()

    def save_preset(self):
        folder = self.dataset_panel.dataset_folder
        compound_id = self.dataset_panel.get_compound_id()
        if not folder:
            messagebox.showerror("No dataset loaded", "Load a dataset folder first.")
            return
        if not compound_id:
            messagebox.showerror("No compound ID", "Enter a compound ID before saving.")
            return
        try:
            path = save_preset(folder, compound_id, "garlic", self._collect_state())
            messagebox.showinfo("Preset saved", f"Saved to:\n{path}")
        except Exception as e:
            messagebox.showerror("Error saving preset", str(e))

    def load_preset(self):
        try:
            data = load_preset_dialog(self.dataset_panel.dataset_folder, "garlic")
        except Exception as e:
            messagebox.showerror("Error loading preset", str(e))
            return
        if data is None:
            return
        self._apply_state(data)

    # ------------------------------------------------------------------ #
    # Script generation
    # ------------------------------------------------------------------ #
    def _on_nuclei_changed(self):
        # v5.2: on_change from the nucleus list used to only re-run weight
        # auto-fill; now it also triggers a live regeneration, so editing
        # an isotope label, n, A, or weight anywhere in the nucleus list
        # (including rows added later via "+Add nucleus", which the
        # one-time _wire_live_regenerate walk never sees) still updates
        # the script promptly.
        self._auto_fill_weights()
        self._debounced_generate()

    def _auto_fill_weights(self):
        """Collect nuclei rows, auto-detect same-element groups, write any
        auto-filled weight back to the visible % field immediately (called
        live on every isotope edit, not just when generating), and return
        (nuclei, backgrounds, overrides) for reuse in generate(). overrides
        (v6.3) is a list of (g_override, lwpp_override) in the SAME order
        as nuclei - detect_same_element_groups() preserves list order/
        length exactly (only rewrites the weight field in place), so
        row_value_pairs's order stays valid for zipping against it."""
        row_value_pairs = []
        for row in self.nuc_list.rows:
            vals = row.get_values()
            if vals is not None:
                row_value_pairs.append((row, vals))
        nuclei = [v for row, v in row_value_pairs]

        nuclei, backgrounds = detect_same_element_groups(nuclei)
        overrides = [row.get_override() for row, _ in row_value_pairs]

        for (row, orig_vals), (iso, n, a, weight) in zip(row_value_pairs, nuclei):
            if weight != orig_vals[3]:
                row.weight_entry.delete(0, "end")
                row.weight_entry.insert(0, weight)

        return nuclei, backgrounds, overrides

    def _on_plot_crop(self, lo, hi):
        # v5.6: the plot's Crop control used to be view-only - it changed
        # what you SAW but never touched Exp.Range, so the generated
        # script (and esfit's actual compute cost) didn't shrink at all.
        # garlic's Exp.Range is expressed as center+sweep rather than
        # lo/hi directly, so convert. lo=hi=None means "crop was cleared"
        # (Full view) - deliberately NOT touching center/sweep in that
        # case, since zooming back out to look at the whole trace
        # shouldn't silently blow away a range you'd already dialed in.
        if lo is None or hi is None:
            return
        center = (lo + hi) / 2.0
        sweep = hi - lo
        self.center_entry.delete(0, "end"); self.center_entry.insert(0, f"{center:.4f}")
        self.sweep_entry.delete(0, "end"); self.sweep_entry.insert(0, f"{sweep:.4f}")
        # Same reasoning as PepperTab: shrink Exp.nPoints to match the real
        # data actually inside the crop window, so the plain simulation
        # call stays fast too (esfit's own block self-corrects nPoints
        # regardless, see build_esfit_block).
        fm = self.dataset_panel.field_mT
        if fm is not None:
            n_in_range = sum(1 for x in fm if lo <= x <= hi)
            if n_in_range > 1:
                self.npoints_entry.delete(0, "end")
                self.npoints_entry.insert(0, str(n_in_range))
        self._debounced_generate()

    def _set_gen_status(self, text, error=False):
        self.gen_status.configure(text=text, text_color="#E24B4A" if error else "#5DCAA5")

    def _debounced_generate(self, *_args):
        # v5.2: replaces the old blind 700ms polling timer (which kept
        # regenerating even with nothing changed) with a short debounce
        # after an ACTUAL change - see _wire_live_regenerate, which hooks
        # every relevant widget to call this. 220ms is short enough to
        # feel immediate but avoids regenerating on every single keystroke
        # while still typing a number.
        if getattr(self, "_gen_after_id", None) is not None:
            self.after_cancel(self._gen_after_id)
        self._gen_after_id = self.after(220, self._do_auto_generate)

    def _do_auto_generate(self):
        self._gen_after_id = None
        self.generate(auto=True)

    def _wire_live_regenerate(self, root):
        """Walk every descendant of `root` once (after the whole tab is
        built) and hook it so any change triggers _debounced_generate:
        CTkEntry via <KeyRelease>, CTkCheckBox via its BooleanVar, and any
        widget with its own `command=` (buttons, option menus, segmented
        buttons - covers +Add/remove/->100% nucleus buttons, dropdowns,
        etc.) by wrapping the existing command so both still run."""
        def hook(widget):
            if isinstance(widget, ctk.CTkEntry):
                widget.bind("<KeyRelease>", self._debounced_generate, add="+")
            elif isinstance(widget, ctk.CTkCheckBox):
                var = getattr(widget, "_variable", None)
                if var is not None:
                    var.trace_add("write", self._debounced_generate)
            else:
                try:
                    existing = widget.cget("command")
                except Exception:
                    existing = None
                if existing not in (None, ""):
                    def wrapped(*a, _orig=existing, **kw):
                        _orig(*a, **kw)
                        self._debounced_generate()
                    try:
                        widget.configure(command=wrapped)
                    except Exception:
                        pass
            for child in widget.winfo_children():
                hook(child)
        hook(root)


    def generate(self, auto=False):
        try:
            g = self.g_entry.get().strip()
            lwpp = self.broadening.get_lwpp_full()
            tcorr = self.tcorr_entry.get().strip()

            # Auto-recognize same-element isotope groups (e.g. 47Ti + 49Ti
            # entered as separate rows): auto-fills natural-abundance
            # weights for rows still at the default "100" (also runs live
            # on every isotope edit, not just here), and computes the
            # implied I=0 background as an extra weightless component.
            nuclei, backgrounds, nuclei_overrides = self._auto_fill_weights()

            mwfreq = self.mwfreq_entry.get().strip()
            center = self.center_entry.get().strip()
            sweep = self.sweep_entry.get().strip()
            npoints = self.npoints_entry.get().strip()
            harmonic = self.harmonic_entry.get().strip()
            modamp = self.modamp_entry.get().strip()
            temp = self.temp_entry.get().strip()

            if not all([g, lwpp, mwfreq, center, sweep, npoints, harmonic]):
                raise ValueError("Please fill in all required fields (g, lwpp, "
                                  "mwFreq, center field, sweep width, nPoints, harmonic).")

            if self.broadening.is_gl_split_enabled():
                gauss = self.broadening.gauss_entry.get().strip()
                lorentz = self.broadening.lorentz_entry.get().strip()
                if not gauss or not lorentz:
                    raise ValueError("Gaussian/Lorentzian split is enabled but one of "
                                      "the two values isn't filled in.")

            # NOW validate that every nucleus row actually has its hyperfine
            # value filled in - deliberately checked AFTER the weight
            # auto-fill above, so % fields still update even if you haven't
            # filled in A yet.
            for iso, n, a_vals, weight in nuclei:
                if any(v == "" for v in a_vals):
                    raise ValueError(f"Missing hyperfine value for nucleus '{iso}' - "
                                      "fill in its A (mT) field before generating "
                                      "(the % weight above has already been updated).")

            # v4.9: block generating a script where A is set on a nucleus
            # with NO nuclear spin (I=0) - that A value would silently
            # have zero effect on the simulation (there's no nuclear spin
            # for the electron to couple to). Caught a real case: '48Ti'
            # (I=0) had an A value entered that did nothing; the actual
            # hyperfine-active isotopes are 47Ti/49Ti and need their OWN
            # separate weighted rows instead.
            for iso, n, a_vals, weight in nuclei:
                mass, element = parse_isotope_label(iso)
                if mass is None:
                    continue
                spin = lookup_isotope_spin(element, mass)
                if spin == '0' and any(float(v) != 0 for v in a_vals if v != ""):
                    raise ValueError(
                        f"'{iso}' has I=0 (no nuclear spin) - its A value(s) would "
                        "have NO effect on the simulation; hyperfine coupling only "
                        "applies to isotopes with nonzero nuclear spin. If you meant "
                        "this row as the no-hyperfine background, set its A to 0 (or "
                        "just leave it - EasySpin ignores A for I=0 regardless). If "
                        "you're after real hyperfine splitting for this element, add "
                        "a separate row for one of its I>0 isotopes instead.")

            # Every nuclei row is its own separate weighted Sys component
            # (Sys = {Sys1, Sys2, ...}), never combined onto one shared
            # molecule - count how many total components this needs.
            num_components = len(nuclei) + len(backgrounds)

            self.fit_section.sync_nucleus_rows([iso for iso, n, a, weight in nuclei])

            def nucs_block(var_name, entries):
                """Build Sys.Nucs/n/A lines for `entries` = [(iso, n, A_mT), ...],
                writing into the MATLAB variable `var_name` (e.g. 'Sys' or 'Sys2').
                Garlic supports Sys.n natively, unlike pepper/chili."""
                isotopes = ",".join(iso for iso, n, a in entries)
                counts = ", ".join(n for iso, n, a in entries)
                a_vals = ", ".join(a for iso, n, a in entries)
                return [
                    f"{var_name}.Nucs = '{isotopes}';",
                    f"{var_name}.n = [{counts}];",
                    f"{var_name}.A = unitconvert([{a_vals}], 'mT->MHz');  "
                    "% entered in mT, converted via EasySpin's own unitconvert()",
                ]

            lines = []
            lines.append("%% Auto-generated by EasySpin Input Builder (garlic)")
            lines.append(f"%% Generated: {datetime.datetime.now().isoformat(timespec='seconds')}")
            lines.append("%% NOTE: cross-checked against EasySpin 6.0.12 documentation as of")
            lines.append("%% 2026-08-14, but please verify against your installed version.")
            lines.append("")
            lines.append("clear Sys Exp Opt")
            lines.append("")
            lines.append("% --- Spin system ---")
            lines.append("% garlic is restricted to S=1/2 doublet radicals; Sys.S is not set.")
            lines.append(f"Sys.g = {g};")

            if num_components <= 1:
                # Simple, single-component case - plain Sys, no weight needed.
                if nuclei:
                    iso, n, a, weight = nuclei[0]
                    if is_zero_spin_isotope(iso):
                        lines.append(f"% '{iso}' has I=0 (no nuclear spin) - Sys.Nucs/Sys.A "
                                      "omitted (hyperfine has no effect for a spinless nucleus).")
                    else:
                        lines.extend(nucs_block("Sys", [(iso, n, a[0])]))
                lines.append(f"Sys.lwpp = {lwpp};  % mT, peak-to-peak linewidth")
                if tcorr:
                    lines.append(f"Sys.tcorr = {tcorr};  % s, rotational correlation time")
                sim_sys_expr = "Sys"
            else:
                lines.append("%")
                lines.append("% Each nucleus row is its own separate weighted Sys component -")
                lines.append("% Sys = {Sys1, Sys2, ...} per EasySpin's documented mixture")
                lines.append("% convention (per esfit's own docs, weight is RELATIVE, not")
                lines.append("% absolute - [0.7 0.3] is equivalent to [14 6] - so these don't")
                lines.append("% need to sum to 1 or 100).")
                lines.append(f"Sys.lwpp = {lwpp};  % mT, peak-to-peak linewidth")
                if tcorr:
                    lines.append(f"Sys.tcorr = {tcorr};  % s, rotational correlation time")

                component_vars = []
                first = True
                for i, (iso, n, a, weight) in enumerate(nuclei):
                    var = "Sys" if first else f"Sys{len(component_vars) + 1}"
                    if not first:
                        lines.append("")
                        lines.append(f"{var} = Sys;")
                        # v6.3: NEW - optional per-component starting g/lwpp
                        # override (the "Own starting g/lwpp" checkbox on
                        # this nucleus row), applied AFTER the "{var} = Sys;"
                        # copy so it overwrites the inherited shared value.
                        # Only meaningful for non-first components (the
                        # first component IS Sys itself - there's no copy
                        # step to override after; use the main g/lwpp
                        # fields directly for that one instead).
                        g_ov, lwpp_ov = (nuclei_overrides[i] if i < len(nuclei_overrides)
                                         else (None, None))
                        if g_ov:
                            lines.append(f"{var}.g = {g_ov};  % own starting g "
                                          "(overrides the shared value above)")
                        if lwpp_ov:
                            lwpp_parts = lwpp_ov.split()
                            if len(lwpp_parts) == 1:
                                lines.append(f"{var}.lwpp = {lwpp_parts[0]};  % mT - own "
                                              "starting lwpp (overrides the shared value above)")
                            else:
                                lines.append(f"{var}.lwpp = [{' '.join(lwpp_parts[:2])}];  "
                                              "% mT [Gaussian Lorentzian] - own starting lwpp "
                                              "(overrides the shared value above)")
                    if is_zero_spin_isotope(iso):
                        lines.append(f"% '{iso}' has I=0 (no nuclear spin) - Sys.Nucs/Sys.A "
                                      "omitted (hyperfine has no effect for a spinless nucleus).")
                    else:
                        lines.extend(nucs_block(var, [(iso, n, a[0])]))
                    weight_frac = float(weight) / 100.0
                    lines.append(f"{var}.weight = {weight_frac:.6g};  % {weight}% of the signal")
                    component_vars.append(var)
                    first = False

                for element, (background, missing) in backgrounds.items():
                    bg_var = "Sys" if first else f"Sys{len(component_vars) + 1}"
                    if not first:
                        lines.append("")
                        lines.append(f"{bg_var} = Sys;")
                    lines.append(f"{bg_var}.weight = {background:.6g};  "
                                  f"% {background*100:.4g}% - I=0 {element} isotopes "
                                  "(auto-detected), no hyperfine")
                    if missing:
                        lines.append(f"% NOTE: {element}'s other active isotope(s) "
                                      f"{', '.join(missing)} weren't entered as rows and "
                                      "are being folded into this background - add them "
                                      "as separate rows if they matter for your system.")
                    component_vars.append(bg_var)
                    first = False

                sim_sys_expr = "{" + ", ".join(component_vars) + "}"

            lines.append("")
            lines.append("% --- Experimental parameters ---")
            lines.append(f"Exp.mwFreq = {mwfreq};  % GHz")
            lines.append(f"Exp.CenterSweep = [{center} {sweep}];  % [center, sweep width] in mT")
            lines.append(f"Exp.nPoints = {npoints};")
            lines.append(f"Exp.Harmonic = {harmonic};")
            if modamp:
                lines.append(f"Exp.ModAmp = {modamp};  % mT")
            if temp:
                lines.append(f"Exp.Temperature = {temp};  % K")
            lines.append("")
            lines.append("% --- Run simulation ---")
            lines.append(f"[B, spc] = garlic({sim_sys_expr}, Exp);")
            lines.append("")
            lines.append("figure")
            lines.append("plot(B, spc)")
            lines.append("xlabel('Magnetic field (mT)')")
            lines.append("ylabel('Intensity (a.u.)')")
            lines.append("title('Garlic simulation')")

            post_lines = []
            if self.fit_section.is_enabled():
                if not self.dataset_panel.has_experimental_trace():
                    raise ValueError("Fit (esfit) is enabled but no experimental "
                                      "dataset with a data trace is loaded - load a "
                                      "dataset above first, or uncheck the fit option.")
                global_vary = self.fit_section.collect_vary()  # {"g":delta, "lwpp":delta}
                if "lwpp" in global_vary and self.broadening.is_gl_split_enabled():
                    # Sys.lwpp is a 2-vector [Gaussian Lorentzian] when the
                    # split is active - match that shape for the vary delta
                    # too, so esfit's Sys0/Vary shapes align.
                    delta = global_vary["lwpp"]
                    global_vary["lwpp"] = f"[{delta} {delta}]"
                nucleus_vary = self.fit_section.get_nucleus_vary_list()

                # Build one vary spec per component, in the SAME order as
                # sim_sys_expr's components, plus each nucleus row's own
                # "Vary A"/"Vary weight" if checked. Weight's delta is
                # entered as a percentage but Sys.weight itself is a
                # fraction, so convert before writing it out. The
                # auto-detected background component(s) get no A/weight
                # vary (no A field, and no associated UI row), just
                # whatever global vary survives the g/gStrain restriction
                # below. esfit's own docs confirm this multi-component
                # Vary shape: {{varySys1,varySys2}} paired with
                # Sys0={{Sys1,Sys2},Exp0}.
                vary_specs = []
                keep_g_flags = []
                for i, (iso, n, a, weight) in enumerate(nuclei):
                    spec = dict(global_vary)
                    indep_g = False
                    if i < len(nucleus_vary):
                        checked, delta, w_checked, w_delta, indep_g = nucleus_vary[i]
                        if checked and not is_zero_spin_isotope(iso):
                            spec["A"] = f"unitconvert([{delta}], 'mT->MHz')"
                        if w_checked:
                            spec["weight"] = f"{float(w_delta) / 100.0:.6g}"
                    vary_specs.append(spec)
                    # v6.4: FIXED (real bug, reported directly) - giving a
                    # component its own starting g (v6.3's "Own starting
                    # g/lwpp" checkbox) did NOT also let that g vary - it
                    # only changed the STARTING point, while the separate
                    # v6.2 "Independent g" checkbox (in a different
                    # section entirely) was the only thing that let it
                    # actually fit independently. That split produced
                    # exactly this bug: set a distinct starting g, still
                    # only see ONE g in the esfit table. There's little
                    # reason to override a component's starting g and then
                    # force it right back onto the shared value anyway, so
                    # an active g override now ALSO implies independent-g,
                    # without needing the separate checkbox checked too.
                    g_override_here = (nuclei_overrides[i][0]
                                        if i < len(nuclei_overrides) else None)
                    keep_g_flags.append(indep_g or bool(g_override_here))
                for _element in backgrounds:
                    vary_specs.append(dict(global_vary))
                    keep_g_flags.append(False)
                if not vary_specs:
                    # No nuclei/backgrounds at all - the lone "Sys" still
                    # needs exactly one (possibly empty) vary spec.
                    vary_specs = [dict(global_vary)]
                    keep_g_flags = [False]

                # v6.0: FIXED (real bug, caught while adding 47Ti/49Ti
                # satellite components to a Ti(III) garlic fit): "g" (and
                # "gStrain") was being copied into EVERY component's vary
                # spec via dict(global_vary) above, since all components
                # start life as copies of the same base Sys. But esfit has
                # no way to tie separate Sys/Sys2/Sys3 structs' parameters
                # together during optimization - each SysVary.g entry is a
                # fully independent free parameter. With multiple isotope
                # components sharing one electronic g-value (isotopes
                # don't change g), that let low-weight satellite
                # components' g drift away from the main component's g
                # with almost no fit-quality penalty (their weight in the
                # total spectrum is small), which is a likely cause of
                # satellite lines converging far from their true position
                # once more than one component is present. g (and
                # gStrain, same reasoning) is now only left free on the
                # FIRST component; every other component keeps it fixed
                # at its copied starting value instead of getting its own
                # independent vary range.
                #
                # v6.2: UPDATED - a component can opt out of this
                # restriction via its "Independent g" checkbox (see
                # sync_nucleus_rows) when the fit is deliberately testing
                # whether that component is a distinct species rather
                # than an isotope satellite of the same one. keep_g_flags
                # tracks that per component; index 0 always keeps g
                # regardless (it's the reference every other g gets
                # compared against).
                if len(vary_specs) > 1:
                    for idx, spec in enumerate(vary_specs[1:], start=1):
                        if not keep_g_flags[idx]:
                            spec.pop("g", None)
                            spec.pop("gStrain", None)

                # v6.1: NEW - warn (don't block) when a nucleus with real
                # hyperfine (I>0) is present but its "Vary A" checkbox is
                # unchecked while fitting is enabled. A is the ONLY
                # parameter that produces hyperfine splitting - leaving it
                # fixed makes it structurally impossible for esfit to find
                # the satellite's actual position, and this is easy to
                # miss since the generated script still looks complete
                # (Sys.A is there, just never varied). Not a hard error,
                # since deliberately fixing A to a known literature value
                # is legitimate - just flagged inline so it's not silent.
                unvaried_hyperfine = [
                    iso for i, (iso, n, a, weight) in enumerate(nuclei)
                    if not is_zero_spin_isotope(iso)
                    and not (i < len(nucleus_vary) and nucleus_vary[i][0])
                ]
                if unvaried_hyperfine:
                    lines.append("")
                    lines.append("% NOTE: the following nucleus/nuclei have hyperfine (I>0) but "
                                 "'Vary A'")
                    lines.append(f"% is unchecked - {', '.join(unvaried_hyperfine)}: A stays fixed "
                                 "at its starting")
                    lines.append("% value below and esfit CANNOT determine the actual hyperfine "
                                 "splitting.")
                    lines.append("% Check 'Vary A' in the Fit section above if you want this fit.")

                indep_g_labels = [iso for i, (iso, n, a, weight) in enumerate(nuclei)
                                   if i > 0 and i < len(keep_g_flags) and keep_g_flags[i]]
                if indep_g_labels:
                    lines.append("")
                    lines.append(f"% NOTE: 'Independent g' is ON for {', '.join(indep_g_labels)} - "
                                 "its g is fit")
                    lines.append("% separately from the first component's g, overriding this "
                                 "tool's usual")
                    lines.append("% shared-g assumption. Only meaningful if you're deliberately "
                                 "testing whether")
                    lines.append("% that component is a distinct species, not a normal isotope "
                                 "satellite fit.")

                if not any(vary_specs):
                    raise ValueError("Fit (esfit) is enabled but no 'Vary' checkbox is "
                                      "checked anywhere - esfit needs at least one "
                                      "parameter to fit. Check 'Vary g-value', 'Vary "
                                      "Linewidth', or a nucleus's 'Vary A' above (or "
                                      "uncheck the fit option entirely).")
                method = self.fit_section.get_method()
                target = self.fit_section.get_target()
                center_f = float(center)
                sweep_f = float(sweep)
                field_lo = center_f - sweep_f / 2
                field_hi = center_f + sweep_f / 2
                lines.append("")
                fit_lines, post_lines = build_esfit_block(
                    fcn_name="garlic",
                    sim_sys_expr=sim_sys_expr,
                    vary_specs=vary_specs,
                    method=method,
                    target=target,
                    data_path=self.dataset_panel.data_path,
                    field_lo=field_lo,
                    field_hi=field_hi,
                    npoints=npoints,
                )
                lines.extend(fit_lines)

            self.output_box.delete("1.0", "end")
            self.output_box.insert("1.0", "\n".join(lines))
            self.output_box2.delete("1.0", "end")
            self.output_box2.insert("1.0", "\n".join(post_lines))

        except ValueError as e:
            if auto:
                self._set_gen_status(str(e), error=True)
            else:
                messagebox.showerror("Input error", str(e))
            return
        else:
            self._set_gen_status("", error=False)

    def copy_to_clipboard(self):
        text = self.output_box.get("1.0", "end").strip()
        if not text:
            messagebox.showinfo("Nothing to copy", "Generate a script first.")
            return
        self.clipboard_clear()
        self.clipboard_append(text)

    def save_as(self):
        text = self.output_box.get("1.0", "end").strip()
        if not text:
            messagebox.showinfo("Nothing to save", "Generate a script first.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".m",
                                             filetypes=[("MATLAB script", "*.m")],
                                             initialfile="garlic_sim.m")
        if path:
            with open(path, "w") as f:
                f.write(text)

    def copy_to_clipboard2(self):
        text = self.output_box2.get("1.0", "end").strip()
        if not text:
            messagebox.showinfo("Nothing to copy", "Generate a script first (with fit enabled).")
            return
        self.clipboard_clear()
        self.clipboard_append(text)

    def save_as2(self):
        text = self.output_box2.get("1.0", "end").strip()
        if not text:
            messagebox.showinfo("Nothing to save", "Generate a script first (with fit enabled).")
            return
        path = filedialog.asksaveasfilename(defaultextension=".m",
                                             filetypes=[("MATLAB script", "*.m")],
                                             initialfile="garlic_fit_export.m")
        if path:
            with open(path, "w") as f:
                f.write(text)


##############################################################################
# tab_pepper.py content
##############################################################################

import datetime

import customtkinter as ctk
from tkinter import filedialog, messagebox



class PepperTab(ctk.CTkFrame):
    def __init__(self, master):
        super().__init__(master, fg_color="transparent")

        left = ctk.CTkScrollableFrame(self, width=460)
        left.pack(side="left", fill="both", padx=(0, 10), pady=10)

        right = ctk.CTkFrame(self, fg_color="transparent")  # v5.1: no outer scroll
        right.pack(side="left", fill="both", expand=True, pady=10)

        # --- Dataset panel ---
        self.dataset_panel = DatasetPanel(left, on_loaded=self._on_dataset_loaded)
        self.dataset_panel.pack(anchor="w", fill="x", pady=(0, 4))

        preset_row = ctk.CTkFrame(left, fg_color="transparent")
        preset_row.pack(anchor="w", fill="x", pady=(0, 10))
        ctk.CTkButton(preset_row, text="Save preset", command=self.save_preset).pack(
            side="left", padx=(0, 8))
        ctk.CTkButton(preset_row, text="Load preset", command=self.load_preset).pack(side="left")

        # --- Spin system ---
        section_title(left, "Spin System (Sys)")

        s_row = ctk.CTkFrame(left, fg_color="transparent")
        s_row.pack(anchor="w", fill="x", pady=3)
        ctk.CTkLabel(s_row, text="Electron spin, S:", width=190, anchor="w").pack(side="left")
        self.s_entry = ctk.CTkEntry(s_row, width=80, placeholder_text="e.g. 0.5")
        self.s_entry.insert(0, "0.5")
        self.s_entry.pack(side="left", padx=6)
        ctk.CTkButton(s_row, text="Suggest S...", width=100,
                      command=self.open_electron_helper).pack(side="left", padx=(6, 0))
        info_icon(s_row, "Total electron spin quantum number. S=1/2 for most "
                          "organic radicals and many d1/d9 metal ions. Higher "
                          "S (1, 3/2, 2, 5/2...) for multi-unpaired-electron "
                          "metal centers - use 'Suggest S...' if unsure."
                  ).pack(side="left", padx=(6, 0))

        self.g_mode = ctk.StringVar(value="isotropic")
        mode_frame = ctk.CTkFrame(left, fg_color="transparent")
        mode_frame.pack(anchor="w", pady=(6, 6))
        ctk.CTkLabel(mode_frame, text="g symmetry:", width=190, anchor="w").pack(side="left")
        ctk.CTkOptionMenu(mode_frame, variable=self.g_mode,
                           values=["isotropic", "axial", "rhombic"],
                           command=lambda _: self._refresh_g_fields()).pack(side="left")

        self.g_fields_frame = ctk.CTkFrame(left, fg_color="transparent")
        self.g_fields_frame.pack(anchor="w", fill="x")
        self.g_entries = []
        self._refresh_g_fields()

        # v6.5: NEW - zero-field splitting (D, E/D) + optional D-strain.
        section_title(left, "Zero-Field Splitting (S >= 1)")
        self.zfs = ZFSSection(left)
        self.zfs.pack(anchor="w", fill="x")

        section_title(left, "Line Broadening")
        self.broadening = BroadeningSection(left, mode="tensor")
        self.broadening.pack(anchor="w", fill="x")
        self.lwpp_entry = self.broadening.lwpp_entry  # alias for existing code below

        nuc_header = ctk.CTkFrame(left, fg_color="transparent")
        nuc_header.pack(anchor="w", fill="x", pady=(10, 2))
        ctk.CTkLabel(nuc_header, text="Nuclei (hyperfine tensor entered in mT):",
                     anchor="w").pack(side="left")
        info_icon(nuc_header, "Isotope label (e.g. 14N, 1H, 55Mn), how many equivalent "
                               "nuclei (n - auto-expanded, see note below), and the "
                               "principal hyperfine tensor values Axx, Ayy, Azz in mT. "
                               "For an axially symmetric tensor, set Axx = Ayy."
                  ).pack(side="left", padx=(6, 0))
        info_note(left, "pepper does not support Sys.n - if count > 1, this tool will "
                         "auto-expand into repeated rows in the generated script.")
        self.nuc_list = NucleusList(left, mode="tensor", on_change=self._on_nuclei_changed)
        self.nuc_list.pack(anchor="w", fill="x")
        info_note(left, "If you enter 2+ isotopes of the SAME element as separate rows "
                         "(e.g. 47Ti and 49Ti), their weights auto-fill from natural "
                         "abundance and the remaining I=0 background is added "
                         "automatically - no separate setup needed.")

        # --- Experimental parameters ---
        section_title(left, "Experiment (Exp)")
        self.mwfreq_entry, _ = labeled_entry(
            left, "Microwave freq (GHz):", default="9.5",
            tooltip="X-band EPR (most common benchtop instruments): ~9.2-9.9 "
                    "GHz. Q-band: ~34 GHz. Use the value from your .DSC file "
                    "(MWFQ) for real data."
        )
        self.rangelo_entry, _ = labeled_entry(
            left, "Field range low (mT):", default="300",
            tooltip="Lower bound of the field sweep, in mT."
        )
        self.rangehi_entry, _ = labeled_entry(
            left, "Field range high (mT):", default="380",
            tooltip="Upper bound of the field sweep, in mT. Powder patterns "
                    "with large g- or hyperfine-anisotropy may need a wide "
                    "range to capture all turning points."
        )
        self.npoints_entry, _ = labeled_entry(
            left, "Number of points:", default="1024",
            tooltip="Number of data points across the sweep. 1024-4096 is "
                    "typical."
        )
        self.harmonic_entry, _ = labeled_entry(
            left, "Harmonic:", default="1",
            tooltip="Detection harmonic. 1 = first-derivative (standard cw "
                    "EPR display). 0 = absorption (rare). 2 = second-"
                    "derivative (rare)."
        )
        self.temp_entry, _ = labeled_entry(
            left, "Temperature (K, optional):", default="",
            tooltip="Sample temperature in K. Affects relative Boltzmann "
                    "populations across spin sublevels - mainly matters when "
                    "S > 1/2 or at very low temperatures."
        )

        section_title(left, "Options (Opt)")
        self.gridsize_entry, _ = labeled_entry(
            left, "Powder grid size, GridSize:", default="91",
            placeholder="orientation grid resolution",
            tooltip="Number of orientations sampled over the powder average. "
                    "Higher = smoother but slower. 91 is a commonly used "
                    "starting value - increase if you see jagged/discretized "
                    "artifacts in the simulated spectrum."
        )
        info_note(left, "No numeric default is documented for Opt.GridSize in EasySpin "
                         "6.0.12 - 91 is a commonly used starting value from the docs' own "
                         "examples, not an official default.")

        # v5.2: no more "Generate Script" button - the script now updates
        # live as soon as you change anything (see _wire_live_regenerate),
        # so a manual trigger isn't needed. This status line shows live
        # validation instead: green when the current inputs produce a
        # valid script, red with the message when they don't (auto ticks
        # never pop a blocking dialog).
        self.gen_status = ctk.CTkLabel(
            left, text="", text_color="#5DCAA5", font=ctk.CTkFont(size=11),
            wraplength=400, justify="left", anchor="w")
        self.gen_status.pack(anchor="w", fill="x", pady=(12, 8))


        # --- Fit (esfit) ---
        section_title(left, "Fit to experimental data (esfit)")
        self.fit_section = FitSection(left, fields=[
            ("g", "g-value(s)", "1", ""),
            ("lwpp", "Linewidth", "1", "mT"),
        ])
        self.fit_section.pack(anchor="w", fill="x")

        # --- Output: plot preview + script (v5.1: side-by-side script
        # panels, non-scrolling right pane, so the plot and both script
        # parts are visible together without having to scroll down) ---
        section_title(right, "Preview")
        self.plot_preview = PlotPreview(right, height=440, on_crop_change=self._on_plot_crop)
        self.plot_preview.pack(fill="x", padx=10, pady=(0, 6))

        code_row = ctk.CTkFrame(right, fg_color="transparent")
        code_row.pack(fill="both", expand=True, padx=10, pady=(0, 10))

        part1_col = ctk.CTkFrame(code_row, fg_color="transparent")
        part1_col.pack(side="left", fill="both", expand=True, padx=(0, 5))
        part2_col = ctk.CTkFrame(code_row, fg_color="transparent")
        part2_col.pack(side="left", fill="both", expand=True, padx=(5, 0))

        section_title(part1_col, "Part 1 (run first)")
        self.output_box = ctk.CTkTextbox(part1_col, font=("Courier", 11), height=170)
        self.output_box.pack(fill="both", expand=True, pady=(0, 6))
        btn_row = ctk.CTkFrame(part1_col, fg_color="transparent")
        btn_row.pack(fill="x")
        ctk.CTkButton(btn_row, text="Copy to Clipboard", command=self.copy_to_clipboard).pack(
            side="left", padx=(0, 8))
        ctk.CTkButton(btn_row, text="Save As .m ...", command=self.save_as).pack(side="left")

        section_title(part2_col, "Part 2 (after esfit export)")
        self.output_box2 = ctk.CTkTextbox(part2_col, font=("Courier", 11), height=170)
        self.output_box2.pack(fill="both", expand=True, pady=(0, 6))
        btn_row2 = ctk.CTkFrame(part2_col, fg_color="transparent")
        btn_row2.pack(fill="x")
        ctk.CTkButton(btn_row2, text="Copy to Clipboard", command=self.copy_to_clipboard2).pack(
            side="left", padx=(0, 8))
        ctk.CTkButton(btn_row2, text="Save As .m ...", command=self.save_as2).pack(side="left")

        self._wire_live_regenerate(self)

    # ------------------------------------------------------------------ #
    # Electron-count helper
    # ------------------------------------------------------------------ #
    def open_electron_helper(self):
        ElectronCountDialog(self, on_select=self._apply_suggested_s)

    def _apply_suggested_s(self, S):
        self.s_entry.delete(0, "end")
        self.s_entry.insert(0, f"{S:g}")

    # ------------------------------------------------------------------ #
    # g-field mode
    # ------------------------------------------------------------------ #
    def _refresh_g_fields(self):
        for w in self.g_fields_frame.winfo_children():
            w.destroy()
        self.g_entries = []

        mode = self.g_mode.get()
        row = ctk.CTkFrame(self.g_fields_frame, fg_color="transparent")
        row.pack(anchor="w", pady=3)
        ctk.CTkLabel(row, text="g values:", width=190, anchor="w").pack(side="left")

        if mode == "isotropic":
            labels = ["g"]
        elif mode == "axial":
            labels = ["g_perp", "g_par"]
        else:
            labels = ["gx", "gy", "gz"]

        for lbl in labels:
            e = ctk.CTkEntry(row, width=80, placeholder_text=lbl)
            e.insert(0, "2.00")
            e.pack(side="left", padx=4)
            self.g_entries.append(e)

    # ------------------------------------------------------------------ #
    # Plot preview
    # ------------------------------------------------------------------ #
    def _plot_experimental(self):
        panel = self.dataset_panel
        if not panel.has_experimental_trace():
            return
        self.plot_preview.plot_experimental(panel.field_mT, panel.intensity)

    # ------------------------------------------------------------------ #
    # Dataset auto-fill
    # ------------------------------------------------------------------ #
    def _on_dataset_loaded(self, exp_values, compound_id, folder):
        if "mwFreq_GHz" in exp_values:
            self.mwfreq_entry.delete(0, "end")
            self.mwfreq_entry.insert(0, f"{exp_values['mwFreq_GHz']:.6f}")
        if "center_G" in exp_values and "sweep_G" in exp_values:
            center_mT = exp_values["center_G"] / 10
            sweep_mT = exp_values["sweep_G"] / 10
            lo = center_mT - sweep_mT / 2
            hi = center_mT + sweep_mT / 2
            self.rangelo_entry.delete(0, "end"); self.rangelo_entry.insert(0, f"{lo:.4f}")
            self.rangehi_entry.delete(0, "end"); self.rangehi_entry.insert(0, f"{hi:.4f}")
        if "npoints" in exp_values:
            self.npoints_entry.delete(0, "end")
            self.npoints_entry.insert(0, str(exp_values["npoints"]))
        if "harmonic" in exp_values:
            self.harmonic_entry.delete(0, "end")
            self.harmonic_entry.insert(0, str(exp_values["harmonic"]))

        self._plot_experimental()
        self._debounced_generate()

    # ------------------------------------------------------------------ #
    # Preset save/load
    # ------------------------------------------------------------------ #
    def _collect_state(self):
        return {
            "S": self.s_entry.get().strip(),
            "g_mode": self.g_mode.get(),
            "g_values": [e.get().strip() for e in self.g_entries],
            "lwpp": self.lwpp_entry.get().strip(),
            "nuclei": [
                {"isotope": iso, "n": n, "a_mT": a, "weight": weight,
                 "g_override": g_ov, "lwpp_override": lwpp_ov}
                for (iso, n, a, weight), (g_ov, lwpp_ov) in self.nuc_list.collect_with_overrides()
            ],
            "mwFreq": self.mwfreq_entry.get().strip(),
            "rangelo": self.rangelo_entry.get().strip(),
            "rangehi": self.rangehi_entry.get().strip(),
            "npoints": self.npoints_entry.get().strip(),
            "harmonic": self.harmonic_entry.get().strip(),
            "temperature": self.temp_entry.get().strip(),
            "gridsize": self.gridsize_entry.get().strip(),
            "zfs": self.zfs.collect_state(),  # v6.5
        }

    def _apply_state(self, data):
        self.s_entry.delete(0, "end"); self.s_entry.insert(0, data.get("S", "0.5"))
        # v6.5: older presets have no "zfs" key -> ZFS off (backward compatible)
        self.zfs.apply_state(data.get("zfs"))
        self.g_mode.set(data.get("g_mode", "isotropic"))
        self._refresh_g_fields()
        for e, v in zip(self.g_entries, data.get("g_values", [])):
            e.delete(0, "end"); e.insert(0, v)
        self.lwpp_entry.delete(0, "end"); self.lwpp_entry.insert(0, data.get("lwpp", ""))
        self.nuc_list.load(data.get("nuclei", []))
        self.mwfreq_entry.delete(0, "end"); self.mwfreq_entry.insert(0, data.get("mwFreq", ""))
        self.rangelo_entry.delete(0, "end"); self.rangelo_entry.insert(0, data.get("rangelo", ""))
        self.rangehi_entry.delete(0, "end"); self.rangehi_entry.insert(0, data.get("rangehi", ""))
        self.npoints_entry.delete(0, "end"); self.npoints_entry.insert(0, data.get("npoints", ""))
        self.harmonic_entry.delete(0, "end"); self.harmonic_entry.insert(0, data.get("harmonic", ""))
        self.temp_entry.delete(0, "end"); self.temp_entry.insert(0, data.get("temperature", ""))
        self.gridsize_entry.delete(0, "end"); self.gridsize_entry.insert(0, data.get("gridsize", ""))
        self._debounced_generate()

    def save_preset(self):
        folder = self.dataset_panel.dataset_folder
        compound_id = self.dataset_panel.get_compound_id()
        if not folder:
            messagebox.showerror("No dataset loaded", "Load a dataset folder first.")
            return
        if not compound_id:
            messagebox.showerror("No compound ID", "Enter a compound ID before saving.")
            return
        try:
            path = save_preset(folder, compound_id, "pepper", self._collect_state())
            messagebox.showinfo("Preset saved", f"Saved to:\n{path}")
        except Exception as e:
            messagebox.showerror("Error saving preset", str(e))

    def load_preset(self):
        try:
            data = load_preset_dialog(self.dataset_panel.dataset_folder, "pepper")
        except Exception as e:
            messagebox.showerror("Error loading preset", str(e))
            return
        if data is None:
            return
        self._apply_state(data)

    # ------------------------------------------------------------------ #
    # Script generation
    # ------------------------------------------------------------------ #
    def _on_nuclei_changed(self):
        # v5.2: on_change from the nucleus list used to only re-run weight
        # auto-fill; now it also triggers a live regeneration, so editing
        # an isotope label, n, A, or weight anywhere in the nucleus list
        # (including rows added later via "+Add nucleus", which the
        # one-time _wire_live_regenerate walk never sees) still updates
        # the script promptly.
        self._auto_fill_weights()
        self._debounced_generate()

    def _auto_fill_weights(self):
        """Collect nuclei rows, auto-detect same-element groups, write any
        auto-filled weight back to the visible % field immediately (called
        live on every isotope edit, not just when generating), and return
        (nuclei, backgrounds, overrides) for reuse in generate(). overrides
        (v6.3) is a list of (g_override, lwpp_override) in the SAME order
        as nuclei - detect_same_element_groups() preserves list order/
        length exactly (only rewrites the weight field in place), so
        row_value_pairs's order stays valid for zipping against it."""
        row_value_pairs = []
        for row in self.nuc_list.rows:
            vals = row.get_values()
            if vals is not None:
                row_value_pairs.append((row, vals))
        nuclei = [v for row, v in row_value_pairs]

        nuclei, backgrounds = detect_same_element_groups(nuclei)
        overrides = [row.get_override() for row, _ in row_value_pairs]

        for (row, orig_vals), (iso, n, a, weight) in zip(row_value_pairs, nuclei):
            if weight != orig_vals[3]:
                row.weight_entry.delete(0, "end")
                row.weight_entry.insert(0, weight)

        return nuclei, backgrounds, overrides

    def _on_plot_crop(self, lo, hi):
        # v5.6: same fix as GarlicTab, but pepper stores Exp.Range as
        # explicit lo/hi entries already, so no center/sweep conversion
        # is needed - just push the values straight through and
        # regenerate. lo=hi=None (Full view) intentionally leaves
        # rangelo/rangehi alone.
        if lo is None or hi is None:
            return
        self.rangelo_entry.delete(0, "end"); self.rangelo_entry.insert(0, f"{lo:.4f}")
        self.rangehi_entry.delete(0, "end"); self.rangehi_entry.insert(0, f"{hi:.4f}")
        # Also shrink Exp.nPoints to match how many REAL loaded data points
        # fall inside the new window - otherwise the plain (non-esfit)
        # simulation call earlier in the script keeps using the old,
        # full-width point count over a now-narrower range: not the
        # data-mismatch bug (esfit's own block now self-corrects that, see
        # build_esfit_block), just wasted compute on an over-resolved
        # simulation that defeats the point of cropping for speed.
        fm = self.dataset_panel.field_mT
        if fm is not None:
            n_in_range = sum(1 for x in fm if lo <= x <= hi)
            if n_in_range > 1:
                self.npoints_entry.delete(0, "end")
                self.npoints_entry.insert(0, str(n_in_range))
        self._debounced_generate()

    def _set_gen_status(self, text, error=False):
        self.gen_status.configure(text=text, text_color="#E24B4A" if error else "#5DCAA5")

    def _debounced_generate(self, *_args):
        # v5.2: replaces the old blind 700ms polling timer (which kept
        # regenerating even with nothing changed) with a short debounce
        # after an ACTUAL change - see _wire_live_regenerate, which hooks
        # every relevant widget to call this. 220ms is short enough to
        # feel immediate but avoids regenerating on every single keystroke
        # while still typing a number.
        if getattr(self, "_gen_after_id", None) is not None:
            self.after_cancel(self._gen_after_id)
        self._gen_after_id = self.after(220, self._do_auto_generate)

    def _do_auto_generate(self):
        self._gen_after_id = None
        self.generate(auto=True)

    def _wire_live_regenerate(self, root):
        """Walk every descendant of `root` once (after the whole tab is
        built) and hook it so any change triggers _debounced_generate:
        CTkEntry via <KeyRelease>, CTkCheckBox via its BooleanVar, and any
        widget with its own `command=` (buttons, option menus, segmented
        buttons - covers +Add/remove/->100% nucleus buttons, dropdowns,
        etc.) by wrapping the existing command so both still run."""
        def hook(widget):
            if isinstance(widget, ctk.CTkEntry):
                widget.bind("<KeyRelease>", self._debounced_generate, add="+")
            elif isinstance(widget, ctk.CTkCheckBox):
                var = getattr(widget, "_variable", None)
                if var is not None:
                    var.trace_add("write", self._debounced_generate)
            else:
                try:
                    existing = widget.cget("command")
                except Exception:
                    existing = None
                if existing not in (None, ""):
                    def wrapped(*a, _orig=existing, **kw):
                        _orig(*a, **kw)
                        self._debounced_generate()
                    try:
                        widget.configure(command=wrapped)
                    except Exception:
                        pass
            for child in widget.winfo_children():
                hook(child)
        hook(root)


    def generate(self, auto=False):
        try:
            S = self.s_entry.get().strip()
            g_vals = [e.get().strip() for e in self.g_entries]
            if any(v == "" for v in g_vals):
                raise ValueError("Please fill in all g-value fields.")
            if not S:
                raise ValueError("Please fill in Sys.S (electron spin).")

            mode = self.g_mode.get()
            if mode == "isotropic":
                g_str = g_vals[0]
            elif mode == "axial":
                # v5.4: FIXED - this used to expand to a 3-element
                # "[gp gp ga]" (duplicating g_perp), which forced esfit to
                # treat g(1,1) and g(1,2) as two INDEPENDENT parameters
                # during fitting (SysVary has no way to tie array elements
                # together) - exactly the reported bug: selecting "axial"
                # still produced 3 separately-varying g parameters. Per
                # EasySpin's own docs (garlic.html/pepper.html: "1, 2, or
                # 3 elements... 2 = axial"), Sys.g natively accepts a
                # 2-element [gperp gpar] for axial symmetry - no
                # duplication needed, and now there's no redundant third
                # slot for esfit to split away from the other two.
                gp, ga = g_vals
                g_str = f"[{gp} {ga}]"
            else:
                g_str = f"[{g_vals[0]} {g_vals[1]} {g_vals[2]}]"

            lwpp = self.broadening.get_lwpp_full()

            # Auto-recognize same-element isotope groups (e.g. 47Ti + 49Ti
            # entered as separate rows): auto-fills natural-abundance
            # weights for rows still at the default "100" (also runs live
            # on every isotope edit, not just here), and computes the
            # implied I=0 background as an extra weightless component.
            nuclei, backgrounds, nuclei_overrides = self._auto_fill_weights()

            mwfreq = self.mwfreq_entry.get().strip()
            rlo = self.rangelo_entry.get().strip()
            rhi = self.rangehi_entry.get().strip()
            npoints = self.npoints_entry.get().strip()
            harmonic = self.harmonic_entry.get().strip()
            temp = self.temp_entry.get().strip()
            gridsize = self.gridsize_entry.get().strip()

            if not all([lwpp, mwfreq, rlo, rhi, npoints, harmonic]):
                raise ValueError("Please fill in all required fields (lwpp, mwFreq, "
                                  "field range, nPoints, harmonic).")

            # NOW validate that every nucleus row has its full hyperfine
            # tensor filled in - deliberately checked AFTER the weight
            # auto-fill above, so % fields still update even if Axx/Ayy/Azz
            # aren't filled in yet.
            for iso, n, a_vals, weight in nuclei:
                if any(v == "" for v in a_vals):
                    raise ValueError(f"Missing hyperfine value for nucleus '{iso}' - "
                                      "fill in its Axx/Ayy/Azz (mT) fields before "
                                      "generating (the % weight above has already "
                                      "been updated).")

            # v4.9: same I=0 check as garlic's path above - A(xx/yy/zz) on
            # a nucleus with no nuclear spin would have zero effect.
            for iso, n, a_vals, weight in nuclei:
                mass, element = parse_isotope_label(iso)
                if mass is None:
                    continue
                spin = lookup_isotope_spin(element, mass)
                if spin == '0' and any(float(v) != 0 for v in a_vals if v != ""):
                    raise ValueError(
                        f"'{iso}' has I=0 (no nuclear spin) - its Axx/Ayy/Azz "
                        "value(s) would have NO effect on the simulation; "
                        "hyperfine coupling only applies to isotopes with nonzero "
                        "nuclear spin. If you meant this row as the no-hyperfine "
                        "background, set its A values to 0 (or just leave them - "
                        "EasySpin ignores A for I=0 regardless). If you're after "
                        "real hyperfine splitting for this element, add a separate "
                        "row for one of its I>0 isotopes instead.")

            # Every nuclei row is its own separate weighted Sys component
            # (Sys = {Sys1, Sys2, ...}), never combined onto one shared
            # molecule - count how many total components this needs.
            num_components = len(nuclei) + len(backgrounds)

            self.fit_section.sync_nucleus_rows([iso for iso, n, a, weight in nuclei])

            # Line broadening: Gaussian+Lorentzian split, HStrain, gStrain
            # are all optional and global to every component (they're all
            # copies of the same base system) - validate their values if
            # enabled, and sync the Fit section's dynamic vary rows for
            # them.
            if self.broadening.is_gl_split_enabled():
                lw_extra_vals = self.broadening.get_lw_extra_raw_values()
                if any(v == "" for v in lw_extra_vals):
                    raise ValueError("Gaussian+Lorentzian split "
                                      "is enabled but one of its fields isn't filled in.")

            extra_broadening_fields = []
            if self.broadening.is_hstrain_enabled():
                hvals = self.broadening.get_hstrain_values()
                if any(v == "" for v in hvals):
                    raise ValueError("H-strain is enabled but Hx/Hy/Hz aren't all filled in.")
                extra_broadening_fields.append(("HStrain", "H-strain", "0.1", "mT"))
            if self.broadening.is_gstrain_enabled():
                gsvals = self.broadening.get_gstrain_values()
                if any(v == "" for v in gsvals):
                    raise ValueError("g-strain is enabled but gx/gy/gz aren't all filled in.")
                extra_broadening_fields.append(("gStrain", "g-strain", "0.001", ""))

            # v6.5: zero-field splitting validation + its Fit vary rows.
            zfs_on = self.zfs.is_enabled()
            zfs_D = zfs_EoD = zfs_factor = None
            zfs_unit = ""
            zfs_dstrain = None
            if zfs_on:
                S_val = _parse_spin_value(S)
                if S_val is None:
                    raise ValueError(f"Can't interpret Sys.S = '{S}' as a number - needed "
                                      "to check that zero-field splitting applies.")
                if S_val < 1:
                    raise ValueError(f"Zero-field splitting is on but S = {S} - ZFS only "
                                      "exists for S >= 1. Set S (e.g. 1, 3/2, 2, 5/2) or "
                                      "uncheck 'Include zero-field splitting'.")
                zfs_D, zfs_EoD, zfs_unit = self.zfs.get_values()
                zfs_factor = self.zfs.unit_factor()
                if self.zfs.is_dstrain_enabled():
                    zfs_dstrain = self.zfs.get_dstrain_values()
                extra_broadening_fields.append(("D", "ZFS D", "0.02", zfs_unit))
                extra_broadening_fields.append(("EoverD", "ZFS E/D", "0.03", "unitless"))
                if zfs_dstrain is not None:
                    extra_broadening_fields.append(("DStrain", "D-strain", "0.01", zfs_unit))
            self.fit_section.sync_extra_fields(extra_broadening_fields)



            def expand_row(iso, n_str, a_vals):
                """Pepper doesn't support Sys.n - a count of 3 becomes three
                repeated Nucs entries within THIS component."""
                try:
                    count = max(1, int(n_str))
                except ValueError:
                    count = 1
                return [(iso, a_vals) for _ in range(count)]

            def nucs_block(var_name, entries):
                """Build Sys.Nucs/A lines (no Sys.n - pepper doesn't support it)
                for `entries` = [(iso, [Axx,Ayy,Azz]), ...], writing into the
                MATLAB variable `var_name` (e.g. 'Sys' or 'Sys2')."""
                out = []
                isotopes = ",".join(iso for iso, a in entries)
                out.append(f"{var_name}.Nucs = '{isotopes}';")
                out.append(f"{var_name}.A = unitconvert([")
                for iso, a in entries:
                    out.append(f"    {a[0]}  {a[1]}  {a[2]}   % {iso}, values in mT")
                out.append("], 'mT->MHz');")
                return out

            lines = []
            lines.append("%% Auto-generated by EasySpin Input Builder (pepper)")
            lines.append(f"%% Generated: {datetime.datetime.now().isoformat(timespec='seconds')}")
            lines.append("%% NOTE: cross-checked against EasySpin 6.0.12 documentation as of")
            lines.append("%% 2026-08-14, but please verify against your installed version.")
            lines.append("")
            lines.append("clear Sys Exp Opt")
            lines.append("")
            lines.append("% --- Spin system ---")
            lines.append(f"Sys.S = {S};")
            lines.append(f"Sys.g = {g_str};")
            if zfs_on:
                # v6.5: explicit numeric conversion to MHz (Sys.D's unit per
                # spinsystem.html) - written out so D and E/D can be read/
                # edited directly at the top of the script.
                lines.append("% --- Zero-field splitting (Sys.D = [D E] in MHz) ---")
                lines.append(f"ZFS_D = {zfs_D:.10g};       % D in {zfs_unit}")
                lines.append(f"ZFS_EoverD = {zfs_EoD:.10g};  % E/D (0 = axial, 1/3 = fully rhombic)")
                lines.append(f"ZFS_to_MHz = {zfs_factor:.10g};  % {zfs_unit} -> MHz"
                             + ("  (1 cm^-1 = 29979.2458 MHz)" if zfs_unit == "cm^-1" else ""))
                lines.append("Sys.D = [ZFS_D, ZFS_EoverD*ZFS_D] * ZFS_to_MHz;  % [D E], MHz")
                if zfs_dstrain is not None:
                    lines.append(f"Sys.DStrain = [{zfs_dstrain[0]:.10g} {zfs_dstrain[1]:.10g}] "
                                 "* ZFS_to_MHz;  % [FWHM_D FWHM_E], MHz")
                if not temp:
                    lines.append("% NOTE: no Exp.Temperature set - all sublevels equally "
                                 "populated, so the SIGN of D")
                    lines.append("% has no effect on the spectrum. Set a temperature if the "
                                 "sign matters.")

            def broadening_lines():
                """Base lwpp plus optional HStrain/gStrain, applied once on
                the base 'Sys' - every other component inherits these
                through the 'SysN = Sys;' copy pattern, same as lwpp
                already does."""
                out = [f"Sys.lwpp = {lwpp};  % mT"]
                if self.broadening.is_hstrain_enabled():
                    hx, hy, hz = self.broadening.get_hstrain_values()
                    out.append(f"Sys.HStrain = [{hx} {hy} {hz}];  "
                                "% mT, Gaussian FWHM field strain")
                if self.broadening.is_gstrain_enabled():
                    gx, gy, gz = self.broadening.get_gstrain_values()
                    out.append(f"Sys.gStrain = [{gx} {gy} {gz}];  % g-value strain")
                return out

            if num_components <= 1:
                if nuclei:
                    iso, n, a_vals, weight = nuclei[0]
                    if is_zero_spin_isotope(iso):
                        lines.append(f"% '{iso}' has I=0 (no nuclear spin) - Sys.Nucs/Sys.A "
                                      "omitted (hyperfine has no effect for a spinless nucleus).")
                    else:
                        lines.append("% pepper does not support Sys.n - equivalent nuclei are")
                        lines.append("% expanded into repeated rows (one row per nucleus).")
                        lines.extend(nucs_block("Sys", expand_row(iso, n, a_vals)))
                lines.extend(broadening_lines())
                sim_sys_expr = "Sys"
            else:
                lines.append("% pepper does not support Sys.n - equivalent nuclei are")
                lines.append("% expanded into repeated rows within each component.")
                lines.append("%")
                lines.append("% Each nucleus row is its own separate weighted Sys component -")
                lines.append("% Sys = {Sys1, Sys2, ...} per EasySpin's documented mixture")
                lines.append("% convention (per esfit's own docs, weight is RELATIVE, not")
                lines.append("% absolute - [0.7 0.3] is equivalent to [14 6] - so these don't")
                lines.append("% need to sum to 1 or 100).")
                lines.extend(broadening_lines())

                component_vars = []
                first = True
                for i, (iso, n, a_vals, weight) in enumerate(nuclei):
                    var = "Sys" if first else f"Sys{len(component_vars) + 1}"
                    if not first:
                        lines.append("")
                        lines.append(f"{var} = Sys;")
                        # v6.3: NEW - same per-component starting g/lwpp
                        # override as the Garlic tab (see the matching
                        # v6.3 note there). g here may be a vector
                        # (axial/rhombic) - typed as-is, matching whatever
                        # form the main g field above uses.
                        g_ov, lwpp_ov = (nuclei_overrides[i] if i < len(nuclei_overrides)
                                         else (None, None))
                        if g_ov:
                            lines.append(f"{var}.g = {g_ov};  % own starting g "
                                          "(overrides the shared value above)")
                        if lwpp_ov:
                            lwpp_parts = lwpp_ov.split()
                            if len(lwpp_parts) == 1:
                                lines.append(f"{var}.lwpp = {lwpp_parts[0]};  % mT - own "
                                              "starting lwpp (overrides the shared value above)")
                            else:
                                lines.append(f"{var}.lwpp = [{' '.join(lwpp_parts[:2])}];  "
                                              "% mT [Gaussian Lorentzian] - own starting lwpp "
                                              "(overrides the shared value above)")
                    if is_zero_spin_isotope(iso):
                        lines.append(f"% '{iso}' has I=0 (no nuclear spin) - Sys.Nucs/Sys.A "
                                      "omitted (hyperfine has no effect for a spinless nucleus).")
                    else:
                        lines.extend(nucs_block(var, expand_row(iso, n, a_vals)))
                    weight_frac = float(weight) / 100.0
                    lines.append(f"{var}.weight = {weight_frac:.6g};  % {weight}% of the signal")
                    component_vars.append(var)
                    first = False

                for element, (background, missing) in backgrounds.items():
                    bg_var = "Sys" if first else f"Sys{len(component_vars) + 1}"
                    if not first:
                        lines.append("")
                        lines.append(f"{bg_var} = Sys;")
                    lines.append(f"{bg_var}.weight = {background:.6g};  "
                                  f"% {background*100:.4g}% - I=0 {element} isotopes "
                                  "(auto-detected), no hyperfine")
                    if missing:
                        lines.append(f"% NOTE: {element}'s other active isotope(s) "
                                      f"{', '.join(missing)} weren't entered as rows and "
                                      "are being folded into this background - add them "
                                      "as separate rows if they matter for your system.")
                    component_vars.append(bg_var)
                    first = False

                sim_sys_expr = "{" + ", ".join(component_vars) + "}"

            lines.append("")
            lines.append("% --- Experimental parameters ---")
            lines.append(f"Exp.mwFreq = {mwfreq};  % GHz")
            lines.append(f"Exp.Range = [{rlo} {rhi}];  % mT")
            lines.append(f"Exp.nPoints = {npoints};")
            lines.append(f"Exp.Harmonic = {harmonic};")
            if temp:
                lines.append(f"Exp.Temperature = {temp};  % K")
            lines.append("")
            lines.append("% --- Options ---")
            if gridsize:
                lines.append(f"Opt.GridSize = {gridsize};  % powder orientation grid resolution")
            lines.append("")
            lines.append("% --- Run simulation ---")
            lines.append(f"[B, spc] = pepper({sim_sys_expr}, Exp, Opt);")
            lines.append("")
            lines.append("figure")
            lines.append("plot(B, spc)")
            lines.append("xlabel('Magnetic field (mT)')")
            lines.append("ylabel('Intensity (a.u.)')")
            lines.append("title('Pepper simulation')")

            post_lines = []
            if self.fit_section.is_enabled():
                if not self.dataset_panel.has_experimental_trace():
                    raise ValueError("Fit (esfit) is enabled but no experimental "
                                      "dataset with a data trace is loaded - load a "
                                      "dataset above first, or uncheck the fit option.")
                global_vary = self.fit_section.collect_vary()  # {"g":delta, "lwpp":delta}
                if "g" in global_vary and mode != "isotropic":
                    # v5.4: FIXED to match the corrected Sys.g shape above -
                    # axial is 2 elements [gperp gpar], not 3. Previously
                    # this always forced a 3-element SysVary.g even for
                    # axial, which is what let g(1,1) and g(1,2) drift
                    # apart independently during fitting despite the user
                    # selecting "axial" - there's no axial tie to break
                    # anymore since there's no duplicate third slot.
                    delta = global_vary["g"]
                    if mode == "axial":
                        global_vary["g"] = f"[{delta} {delta}]"
                    else:
                        global_vary["g"] = f"[{delta} {delta} {delta}]"
                if "lwpp" in global_vary:
                    split = self.broadening.is_gl_split_enabled()
                    delta = global_vary["lwpp"]
                    if split:
                        global_vary["lwpp"] = f"[{delta} {delta}]"
                if "gStrain" in global_vary:
                    # Sys.gStrain is always written as a 3-element [gx gy
                    # gz] array (see the v5.4 note above - not yet made to
                    # track g's own symmetry-dependent element count), but
                    # the Fit section only exposes ONE delta box for it -
                    # same situation "g" was in before the v5.4 fix.
                    # Without this expansion, SysVary.gStrain came out as a
                    # bare scalar against a 3-element Sys.gStrain, so esfit
                    # only exposed a single combined gStrain row instead of
                    # independent gx/gy/gz parameters.
                    delta = global_vary["gStrain"]
                    global_vary["gStrain"] = f"[{delta} {delta} {delta}]"
                # v6.5: ZFS - esfit varies Sys.D's actual elements [D E]
                # (MHz), so the separate D and E/D deltas are merged into
                # one 2-element SysVary.D. A zero element means "don't vary
                # that one". dE is derived from d(E/D) as d(E/D)*|D_start|.
                d_delta_raw = global_vary.pop("D", None)
                eod_delta_raw = global_vary.pop("EoverD", None)
                dstrain_delta_raw = global_vary.pop("DStrain", None)
                if zfs_on and (d_delta_raw or eod_delta_raw):
                    try:
                        dD = float(d_delta_raw) if d_delta_raw else 0.0
                        dEoD = float(eod_delta_raw) if eod_delta_raw else 0.0
                    except ValueError:
                        raise ValueError("ZFS D and E/D vary deltas must be numbers.")
                    if dEoD and zfs_D == 0:
                        raise ValueError("Can't vary E/D with D = 0 - E would be 0 regardless. "
                                          "Give D a nonzero starting value.")
                    dD_MHz = abs(dD) * zfs_factor
                    dE_MHz = abs(dEoD) * abs(zfs_D) * zfs_factor
                    if dD_MHz or dE_MHz:
                        global_vary["D"] = (f"[{dD_MHz:.6g} {dE_MHz:.6g}];  "
                                            f"% MHz: dD = {dD:g} {zfs_unit}, "
                                            f"dE = {dEoD:g} x |D|")
                if zfs_on and zfs_dstrain is not None and dstrain_delta_raw:
                    try:
                        dS = abs(float(dstrain_delta_raw)) * zfs_factor
                    except ValueError:
                        raise ValueError("D-strain vary delta must be a number.")
                    # Don't vary an element that starts at 0 - a +/- range
                    # around 0 would let esfit try a negative FWHM.
                    dS_D = dS if zfs_dstrain[0] != 0 else 0.0
                    dS_E = dS if zfs_dstrain[1] != 0 else 0.0
                    if dS_D or dS_E:
                        global_vary["DStrain"] = (f"[{dS_D:.6g} {dS_E:.6g}];  % MHz "
                                                  "(0 = element held fixed)")
                nucleus_vary = self.fit_section.get_nucleus_vary_list()

                # Build one vary spec per component, in the SAME order as
                # sim_sys_expr's components, plus each nucleus row's own
                # "Vary A" if checked (one shared checkbox per nucleus,
                # applied to all three tensor components uniformly - not
                # per Axx/Ayy/Azz individually, since that control moved
                # out of the Nuclei list into this Fit section). "Vary
                # weight"'s delta is entered as a percentage but
                # Sys.weight itself is a fraction, so convert before
                # writing it out. The auto-detected background
                # component(s) get no A/weight vary (no A field, no
                # associated UI row), just whatever global vary survives
                # the g/gStrain restriction below. esfit's own docs
                # confirm this multi-component Vary shape:
                # {{varySys1,varySys2}} paired with Sys0={{Sys1,Sys2},Exp0}.
                vary_specs = []
                keep_g_flags = []
                for i, (iso, n, a_vals, weight) in enumerate(nuclei):
                    spec = dict(global_vary)
                    indep_g = False
                    if i < len(nucleus_vary):
                        checked, delta, w_checked, w_delta, indep_g = nucleus_vary[i]
                        if checked and not is_zero_spin_isotope(iso):
                            spec["A"] = f"unitconvert([{delta} {delta} {delta}], 'mT->MHz')"
                        if w_checked:
                            spec["weight"] = f"{float(w_delta) / 100.0:.6g}"
                    vary_specs.append(spec)
                    # v6.4: FIXED - same bug and fix as the Garlic tab (see
                    # the matching v6.4 note there for the full
                    # explanation): an active g override now also implies
                    # independent-g, so the two checkboxes don't need to
                    # be found and checked separately for this to work.
                    g_override_here = (nuclei_overrides[i][0]
                                        if i < len(nuclei_overrides) else None)
                    keep_g_flags.append(indep_g or bool(g_override_here))
                for _element in backgrounds:
                    vary_specs.append(dict(global_vary))
                    keep_g_flags.append(False)
                if not vary_specs:
                    vary_specs = [dict(global_vary)]
                    keep_g_flags = [False]

                # v6.0: FIXED - same bug and fix as the Garlic tab (see the
                # matching v6.0 note there for the full explanation): "g"
                # and "gStrain" were being copied into every component's
                # vary spec, letting low-weight isotope satellite
                # components' g/gStrain drift independently during
                # fitting since esfit can't tie separate Sys structs
                # together. Now only the FIRST component keeps them free.
                #
                # v6.2: UPDATED - same "Independent g" opt-out as the
                # Garlic tab (see the matching v6.2 note there).
                if len(vary_specs) > 1:
                    for idx, spec in enumerate(vary_specs[1:], start=1):
                        if not keep_g_flags[idx]:
                            spec.pop("g", None)
                            spec.pop("gStrain", None)
                            # v6.5: ZFS is a property of the electronic
                            # species just like g - isotope satellites
                            # share it, so only the first component (or
                            # an 'Independent g' distinct species) varies it.
                            spec.pop("D", None)
                            spec.pop("DStrain", None)

                # v6.1: NEW - same warning as the Garlic tab (see the
                # matching v6.1 note there for the full explanation).
                unvaried_hyperfine = [
                    iso for i, (iso, n, a_vals, weight) in enumerate(nuclei)
                    if not is_zero_spin_isotope(iso)
                    and not (i < len(nucleus_vary) and nucleus_vary[i][0])
                ]
                if unvaried_hyperfine:
                    lines.append("")
                    lines.append("% NOTE: the following nucleus/nuclei have hyperfine (I>0) but "
                                 "'Vary A'")
                    lines.append(f"% is unchecked - {', '.join(unvaried_hyperfine)}: A stays fixed "
                                 "at its starting")
                    lines.append("% value below and esfit CANNOT determine the actual hyperfine "
                                 "splitting.")
                    lines.append("% Check 'Vary A' in the Fit section above if you want this fit.")

                indep_g_labels = [iso for i, (iso, n, a_vals, weight) in enumerate(nuclei)
                                   if i > 0 and i < len(keep_g_flags) and keep_g_flags[i]]
                if indep_g_labels:
                    lines.append("")
                    lines.append(f"% NOTE: 'Independent g' is ON for {', '.join(indep_g_labels)} - "
                                 "its g is fit")
                    lines.append("% separately from the first component's g, overriding this "
                                 "tool's usual")
                    lines.append("% shared-g assumption. Only meaningful if you're deliberately "
                                 "testing whether")
                    lines.append("% that component is a distinct species, not a normal isotope "
                                 "satellite fit.")

                if not any(vary_specs):
                    raise ValueError("Fit (esfit) is enabled but no 'Vary' checkbox is "
                                      "checked anywhere - esfit needs at least one "
                                      "parameter to fit. Check 'Vary g-value(s)', 'Vary "
                                      "Linewidth', or a nucleus's 'Vary A' above (or "
                                      "uncheck the fit option entirely).")
                method = self.fit_section.get_method()
                target = self.fit_section.get_target()
                lines.append("")
                fit_lines, post_lines = build_esfit_block(
                    fcn_name="pepper",
                    sim_sys_expr=sim_sys_expr,
                    vary_specs=vary_specs,
                    method=method,
                    target=target,
                    data_path=self.dataset_panel.data_path,
                    field_lo=float(rlo),
                    field_hi=float(rhi),
                    npoints=npoints,
                    has_opt=True,
                )
                lines.extend(fit_lines)

            self.output_box.delete("1.0", "end")
            self.output_box.insert("1.0", "\n".join(lines))
            self.output_box2.delete("1.0", "end")
            self.output_box2.insert("1.0", "\n".join(post_lines))

        except ValueError as e:
            if auto:
                self._set_gen_status(str(e), error=True)
            else:
                messagebox.showerror("Input error", str(e))
            return
        else:
            self._set_gen_status("", error=False)

    def copy_to_clipboard(self):
        text = self.output_box.get("1.0", "end").strip()
        if not text:
            messagebox.showinfo("Nothing to copy", "Generate a script first.")
            return
        self.clipboard_clear()
        self.clipboard_append(text)

    def save_as(self):
        text = self.output_box.get("1.0", "end").strip()
        if not text:
            messagebox.showinfo("Nothing to save", "Generate a script first.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".m",
                                             filetypes=[("MATLAB script", "*.m")],
                                             initialfile="pepper_sim.m")
        if path:
            with open(path, "w") as f:
                f.write(text)

    def copy_to_clipboard2(self):
        text = self.output_box2.get("1.0", "end").strip()
        if not text:
            messagebox.showinfo("Nothing to copy", "Generate a script first (with fit enabled).")
            return
        self.clipboard_clear()
        self.clipboard_append(text)

    def save_as2(self):
        text = self.output_box2.get("1.0", "end").strip()
        if not text:
            messagebox.showinfo("Nothing to save", "Generate a script first (with fit enabled).")
            return
        path = filedialog.asksaveasfilename(defaultextension=".m",
                                             filetypes=[("MATLAB script", "*.m")],
                                             initialfile="pepper_fit_export.m")
        if path:
            with open(path, "w") as f:
                f.write(text)


##############################################################################
# main.py content (entry point)
##############################################################################

import customtkinter as ctk


ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("blue")

APP_VERSION = "6.5"


class EasySpinBuilderApp(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title(f"EasySpin Input Builder v{APP_VERSION}")
        self.geometry("1250x820")

        header = ctk.CTkLabel(
            self, text=f"EasySpin Input Builder v{APP_VERSION}",
            font=ctk.CTkFont(size=20, weight="bold"),
        )
        header.pack(pady=(14, 0))

        subheader = ctk.CTkLabel(
            self,
            text="Fill in spin-system and experimental parameters, generate a MATLAB script for EasySpin.",
            font=ctk.CTkFont(size=12), text_color="gray70",
        )
        subheader.pack(pady=(0, 10))

        self.tabview = ctk.CTkTabview(self)
        self.tabview.pack(fill="both", expand=True, padx=14, pady=(0, 14))

        self.tabview.add("Garlic (Liquid)")
        self.tabview.add("Pepper (Powder/Frozen)")

        GarlicTab(self.tabview.tab("Garlic (Liquid)")).pack(fill="both", expand=True)
        PepperTab(self.tabview.tab("Pepper (Powder/Frozen)")).pack(fill="both", expand=True)


if __name__ == "__main__":
    app = EasySpinBuilderApp()
    app.mainloop()
