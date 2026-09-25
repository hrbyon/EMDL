"""cof_xrd — powder-XRD analysis for covalent organic frameworks (COFs).

Peak search, hexagonal/tetragonal indexing, Pawley refinement, CIF-based pattern
calculation and comparison of a simulation (or structure model) with experiment.
It implements the same procedures as ``xrd-workbench.html``.

Quick start (command line)::

    python cof_xrd.py data.csv                      # auto-detect experiment / simulation columns
    python cof_xrd.py data.csv --cif model.cif --scan --plot
    python cof_xrd.py --selftest

Library::

    import cof_xrd as cx
    series = cx.load_xrd("data.csv")               # list of Series(name, x, y, is_sim, has_neg)
    exp = next(s for s in series if not s.is_sim)
    fit = cx.pawley(exp.x, exp.y, lattice="hex", a0=30.0, tt_range=(2.3, 12))
    print(fit.summary())

Dependencies: numpy, scipy (matplotlib only for ``--plot``).

Notes on COF data
-----------------
* All reflections sit at low angle, so the lattice parameter and the zero shift are
  strongly correlated. Use ``zero_scan`` (or ``--scan``) and, ideally, an internal
  standard.
* Below ~5° the Lorentz–polarization factor changes appreciably across a broad peak
  and skews it to lower angle; ``lp_profile=True`` (default) accounts for this.
* Laboratory CSVs are often background-subtracted and smoothed. The residuals are then
  serially correlated (Durbin–Watson ≪ 2); e.s.d.s are inflated by √(2/DW).
"""
from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from dataclasses import dataclass, field
from fractions import Fraction

import numpy as np
from scipy.optimize import least_squares, lsq_linear
from scipy.signal import find_peaks as _find_peaks

LAM_CU = 1.5406
D2R = math.pi / 180
_trapz = getattr(np, "trapezoid", None) or np.trapz   # numpy ≥ 2 renamed trapz

# Cromer–Mann coefficients (International Tables Vol. C): a1 b1 a2 b2 a3 b3 a4 b4 c
CROMER_MANN = {
    "H": (0.489918, 20.6593, 0.262003, 7.74039, 0.196767, 49.5519, 0.049879, 2.20159, 0.001305),
    "Li": (1.1282, 3.9546, 0.7508, 1.0524, 0.6175, 85.3905, 0.4653, 168.261, 0.0377),
    "B": (2.0545, 23.2185, 1.3326, 1.021, 1.0979, 60.3498, 0.7068, 0.1403, -0.1932),
    "C": (2.31, 20.8439, 1.02, 10.2075, 1.5886, 0.5687, 0.865, 51.6512, 0.2156),
    "N": (12.2126, 0.0057, 3.1322, 9.8933, 2.0125, 28.9975, 1.1663, 0.5826, -11.529),
    "O": (3.0485, 13.2771, 2.2868, 5.7011, 1.5463, 0.3239, 0.867, 32.9089, 0.2508),
    "F": (3.5392, 10.2825, 2.6412, 4.2944, 1.517, 0.2615, 1.0243, 26.1476, 0.2776),
    "Na": (4.7626, 3.285, 3.1736, 8.8422, 1.2674, 0.3136, 1.1128, 129.424, 0.676),
    "Mg": (5.4204, 2.8275, 2.1735, 79.2611, 1.2269, 0.3808, 2.3073, 7.1937, 0.8584),
    "Al": (6.4202, 3.0387, 1.9002, 0.7426, 1.5936, 31.5472, 1.9646, 85.0886, 1.1151),
    "Si": (6.2915, 2.4386, 3.0353, 32.3337, 1.9891, 0.6785, 1.541, 81.6937, 1.1407),
    "P": (6.4345, 1.9067, 4.1791, 27.157, 1.78, 0.526, 1.4908, 68.1645, 1.1149),
    "S": (6.9053, 1.4679, 5.2034, 22.2151, 1.4379, 0.2536, 1.5863, 56.172, 0.8669),
    "Cl": (11.4604, 0.0104, 7.1962, 1.1662, 6.2556, 18.5194, 1.6455, 47.7784, -9.5574),
    "K": (8.2186, 12.7949, 7.4398, 0.7748, 1.0519, 213.187, 0.8659, 41.6841, 1.4228),
    "Ca": (8.6266, 10.4421, 7.3873, 0.6599, 1.5899, 85.7484, 1.0211, 178.437, 1.3751),
    "Ti": (9.7595, 7.8508, 7.3558, 0.5, 1.6991, 35.6338, 1.9021, 116.105, 1.2807),
    "Fe": (11.7695, 4.7611, 7.3573, 0.3072, 3.5222, 15.3535, 2.3045, 76.8805, 1.0369),
    "Co": (12.2841, 4.2791, 7.3409, 0.2784, 4.0034, 13.5359, 2.3488, 71.1692, 1.0118),
    "Ni": (12.8376, 3.8785, 7.292, 0.2565, 4.4438, 12.1763, 2.38, 66.3421, 1.0341),
    "Cu": (13.338, 3.5828, 7.1676, 0.247, 5.6158, 11.3966, 1.6735, 64.8126, 1.191),
    "Zn": (14.0743, 3.2655, 7.0318, 0.2333, 5.1652, 10.3163, 2.41, 58.7097, 1.3041),
    "Br": (17.1789, 2.1723, 5.2358, 16.5796, 5.6377, 0.2609, 3.9851, 41.4328, 2.9557),
    "I": (20.1472, 4.347, 18.9949, 0.3814, 7.5138, 27.766, 2.2735, 66.8776, 4.0712),
}
ATOMIC_Z = {"H": 1, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Na": 11, "Mg": 12, "Al": 13,
            "Si": 14, "P": 15, "S": 16, "Cl": 17, "K": 19, "Ca": 20, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25,
            "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Br": 35, "Pd": 46, "Ag": 47, "I": 53, "Pt": 78, "Au": 79}


def f0(element, s):
    """Atomic scattering factor at s = sinθ/λ (Å⁻¹). Unknown elements are scaled from carbon by Z."""
    c = CROMER_MANN.get(element)
    if c is None:
        return ATOMIC_Z.get(element, 6) / 6.0 * f0("C", s)
    s2 = np.asarray(s, float) ** 2
    return c[0]*np.exp(-c[1]*s2) + c[2]*np.exp(-c[3]*s2) + c[4]*np.exp(-c[5]*s2) + c[6]*np.exp(-c[7]*s2) + c[8]


def lp_factor(tt):
    """Lorentz–polarization factor for powder diffraction, (1 + cos²2θ)/(sin²θ cosθ)."""
    th = np.asarray(tt, float) * D2R / 2
    return (1 + np.cos(2*th)**2) / (np.sin(th)**2 * np.cos(th))


def two_theta(d, lam=LAM_CU):
    return 2 * np.degrees(np.arcsin(lam / (2 * np.asarray(d, float))))


def d_spacing(tt, lam=LAM_CU):
    return lam / (2 * np.sin(np.asarray(tt, float) * D2R / 2))


def scherrer(tt, fwhm, lam=LAM_CU, instrumental=0.05, K=0.9):
    """Crystallite size in nm from FWHM in degrees (instrumental width removed in quadrature)."""
    beta = math.sqrt(max(fwhm**2 - instrumental**2, 1e-8)) * D2R
    return K * lam / (beta * math.cos(tt * D2R / 2)) / 10


def inv_metric(a, b, c, alpha, beta, gamma):
    ca, cb, cg = (math.cos(v * D2R) for v in (alpha, beta, gamma))
    G = np.array([[a*a, a*b*cg, a*c*cb], [a*b*cg, b*b, b*c*ca], [a*c*cb, b*c*ca, c*c]])
    return np.linalg.inv(G)


# ───────────────────────── data loading ─────────────────────────
@dataclass
class Series:
    name: str
    x: np.ndarray
    y: np.ndarray
    is_sim: bool = False
    has_neg: bool = False


def _split(line):
    for sep in ("\t", ",", ";"):
        if sep in line:
            return line.split(sep)
    return line.split()


def _isnum(t):
    try:
        float(str(t).strip())
        return str(t).strip() != ""
    except ValueError:
        return False


def parse_xrd_text(text, name="data"):
    """Parse a CSV/TXT/XY file into a list of Series.

    Handles '(2θ, I)(2θ, I)…' pair layouts whose columns have different lengths (e.g. a
    0.05°-step simulation next to a 0.01°-step measurement) and '2θ + several I' layouts.
    A series is flagged as a simulation when ≥5 % of its points are exactly zero or its
    header contains sim/calc/model.
    """
    lines = [ln for ln in text.splitlines() if ln.strip() and not re.match(r"^\s*[#!]", ln)]
    cells = [_split(ln) for ln in lines]
    ncol = max(len(r) for r in cells)
    header = []
    for r in cells:
        if sum(_isnum(t) for t in r) >= 2 and not any(t.strip() and not _isnum(t) for t in r):
            break
        header.append(r)

    def colname(c):
        for hr in header:
            t = hr[c].strip() if c < len(hr) else ""
            if t and not _isnum(t):
                return t
        return None

    def col(c):
        return np.array([float(r[c]) if c < len(r) and _isnum(r[c]) else np.nan for r in cells])

    def monotonic(v):
        v = v[np.isfinite(v)]
        return len(v) >= 10 and np.mean(np.diff(v) > 0) > 0.98

    paired = ncol >= 4 and ncol % 2 == 0 and all(monotonic(col(c)) for c in range(0, ncol, 2))
    out = []
    pairs = [(c, c + 1) for c in range(0, ncol, 2)] if paired else [(0, c) for c in range(1, ncol)]
    for k, (cx_, cy_) in enumerate(pairs):
        X, Y = col(cx_), col(cy_)
        m = np.isfinite(X) & np.isfinite(Y)
        if m.sum() < 10:
            continue
        nm = colname(cy_) if not paired else (colname(cx_) or colname(cy_))
        if nm and re.search(r"theta|2θ|angle|deg", nm, re.I):
            nm = colname(cy_) or nm
        nm = nm or f"{name} #{k + 1}"
        x, y = X[m], Y[m]
        o = np.argsort(x); x, y = x[o], y[o]
        keep = np.r_[True, np.diff(x) > 1e-9]
        x, y = x[keep], y[keep]
        is_sim = np.mean(y == 0) > 0.05 or bool(re.search(r"sim|calc|\bcal\b|model", nm, re.I))
        out.append(Series(nm, x, y, is_sim, bool((y < 0).any())))
    if not out:
        raise ValueError("no (2θ, intensity) columns found")
    return out


def load_xrd(path):
    with open(path, encoding="utf-8-sig", errors="replace") as fh:
        return parse_xrd_text(fh.read(), name=re.sub(r"\.[^.]+$", "", str(path).split("/")[-1].split("\\")[-1]))


# ───────────────────────── CIF ─────────────────────────
@dataclass
class Structure:
    cell: tuple                      # a, b, c, alpha, beta, gamma
    elements: np.ndarray             # element symbols
    frac: np.ndarray                 # (n, 3) fractional coordinates
    occ: np.ndarray                  # occupancies
    name: str = "structure"


def _parse_symop(op):
    """'x-y+1/2' style operator → (3×3 matrix, translation) without eval."""
    rows, trans = [], []
    parts = op.lower().replace(" ", "").split(",")
    if len(parts) != 3:
        raise ValueError(f"bad symmetry operator: {op}")
    for p in parts:
        row, t = [0.0, 0.0, 0.0], Fraction(0)
        for term in re.findall(r"[+-]?[^+-]+", p):
            sign = -1 if term.startswith("-") else 1
            term = term.lstrip("+-")
            m = re.fullmatch(r"(\d*\.?\d*(?:/\d+)?)\*?([xyz])", term)
            if m:
                coef = float(Fraction(m.group(1))) if m.group(1) else 1.0
                row["xyz".index(m.group(2))] += sign * coef
            else:
                t += sign * Fraction(term)
        rows.append(row); trans.append(float(t))
    return np.array(rows), np.array(trans)


def read_cif(path_or_text):
    """Minimal CIF reader: cell, symmetry operators, atom sites (type, fract xyz, occupancy)."""
    text = path_or_text
    if "\n" not in str(path_or_text):
        with open(path_or_text, encoding="utf-8", errors="replace") as fh:
            text = fh.read()
    toks, lines, i = [], text.splitlines(), 0
    while i < len(lines):
        ln = lines[i]
        if ln.startswith(";"):
            buf = [ln[1:]]; i += 1
            while i < len(lines) and not lines[i].startswith(";"):
                buf.append(lines[i]); i += 1
            toks.append("\n".join(buf)); i += 1; continue
        for m in re.finditer(r"'([^']*)'(?=\s|$)|\"([^\"]*)\"(?=\s|$)|(#.*$)|(\S+)", ln):
            if m.group(3) is not None:
                break
            toks.append(next(g for g in (m.group(1), m.group(2), m.group(4)) if g is not None))
        i += 1
    items, loops, i = {}, [], 0
    while i < len(toks):
        t = toks[i]
        if t.lower() == "loop_":
            heads = []; i += 1
            while i < len(toks) and toks[i].startswith("_"):
                heads.append(toks[i].lower()); i += 1
            vals = []
            while i < len(toks) and not toks[i].startswith("_") and toks[i].lower() != "loop_" and not toks[i].lower().startswith("data_"):
                vals.append(toks[i]); i += 1
            loops.append((heads, [vals[r:r + len(heads)] for r in range(0, len(vals) - len(heads) + 1, len(heads))]))
            continue
        if t.startswith("_") and i + 1 < len(toks):
            items[t.lower()] = toks[i + 1]; i += 2; continue
        i += 1
    def num(v):   # handles e.s.d.s in parentheses and fractions such as 1/2
        t = re.sub(r"\(.*\)", "", str(v)).strip()
        return float(Fraction(t)) if "/" in t else float(t)
    cell = tuple(num(items[k]) for k in ("_cell_length_a", "_cell_length_b", "_cell_length_c",
                                          "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma"))
    ops = ["x,y,z"]
    for heads, rows in loops:
        for key in ("_symmetry_equiv_pos_as_xyz", "_space_group_symop_operation_xyz"):
            if key in heads:
                ops = [r[heads.index(key)] for r in rows]
    ops = [_parse_symop(o) for o in ops]
    site_loop = next(((h, r) for h, r in loops if "_atom_site_fract_x" in h), None)
    if site_loop is None:
        raise ValueError("no _atom_site_fract_x loop in CIF")
    H, rows = site_loop
    ix, iy, iz = (H.index(f"_atom_site_fract_{k}") for k in "xyz")
    it = H.index("_atom_site_type_symbol") if "_atom_site_type_symbol" in H else H.index("_atom_site_label")
    io_ = H.index("_atom_site_occupancy") if "_atom_site_occupancy" in H else None
    els, frs, occs = [], [], []
    for r in rows:
        el = (re.match(r"[A-Z][a-z]?", r[it]) or re.match(r".", "C")).group(0)
        occ = num(r[io_]) if io_ is not None and r[io_] not in (".", "?") else 1.0
        p = np.array([num(r[ix]), num(r[iy]), num(r[iz])])
        seen = []
        for R, t in ops:
            q = (R @ p + t) % 1.0
            if any(np.all(np.minimum(np.abs(q - s), 1 - np.abs(q - s)) < 1e-3) for s in seen):
                continue
            seen.append(q); els.append(el); frs.append(q); occs.append(occ)
    name = items.get("_chemical_name_common") or items.get("_chemical_formula_sum") or "CIF"
    return Structure(cell, np.array(els), np.array(frs), np.array(occs), name)


def structure_reflections(st: Structure, lam=LAM_CU, tt_max=40.0, B=2.0):
    """All reflections to tt_max: returns hkl (n,3), d, 2θ, I = |F|²·Lp·exp(−2Bs²)."""
    a, b, c = st.cell[:3]
    Gs = inv_metric(*st.cell)
    dmin = lam / (2 * math.sin(tt_max * D2R / 2))
    H, K, L = (int(v / dmin) + 1 for v in (a, b, c))
    h, k, l = np.meshgrid(np.arange(-H, H + 1), np.arange(-K, K + 1), np.arange(-L, L + 1), indexing="ij")
    hkl = np.stack([h.ravel(), k.ravel(), l.ravel()], 1)
    hkl = hkl[np.any(hkl != 0, 1)]
    d = 1 / np.sqrt(np.einsum("ni,ij,nj->n", hkl, Gs, hkl))
    hkl, d = hkl[d >= dmin], d[d >= dmin]
    s = 1 / (2 * d)
    F = np.zeros(len(hkl), complex)
    for el in np.unique(st.elements):
        m = st.elements == el
        F += f0(el, s) * (np.exp(2j * np.pi * hkl @ st.frac[m].T) * st.occ[m]).sum(1)
    tt = two_theta(d, lam)
    I = np.abs(F) ** 2 * np.exp(-2 * B * s**2) * lp_factor(tt)
    return hkl, d, tt, I


def simulate_pattern(st: Structure, x, lam=LAM_CU, fwhm=0.1, eta=0.5, B=2.0, lp_across=True):
    """Powder pattern of a structure on grid x (pseudo-Voigt, constant FWHM).

    lp_across=True lets the Lorentz–polarization factor vary across each profile, as it does in a
    measurement; this matters only for broad peaks below ~5°.
    """
    x = np.asarray(x, float)
    _, _, tt, I = structure_reflections(st, lam, x.max() + 2, B)
    y = np.zeros_like(x)
    lpx = lp_factor(x) if lp_across else None
    for t, i in zip(tt, I):
        z = 4 * (x - t) ** 2 / fwhm**2
        v = eta / (1 + z) + (1 - eta) * np.exp(-math.log(2) * z)
        if lpx is not None:
            v = v * lpx / lp_factor(t)
            v *= _trapz(eta / (1 + z) + (1 - eta) * np.exp(-math.log(2) * z), x) / max(_trapz(v, x), 1e-12)
        y += i * v
    return y


# ───────────────────────── peaks & indexing ─────────────────────────
@dataclass
class Peak:
    tt: float
    d: float
    fwhm: float
    height: float
    area: float
    size_nm: float


def find_peaks(x, y, rel_prominence=0.012, smooth_deg=0.08, lam=LAM_CU, instrumental=0.05):
    x, y = np.asarray(x, float), np.asarray(y, float)
    step = (x[-1] - x[0]) / (len(x) - 1)
    w = max(1, int(round(smooth_deg / step)) | 1)
    ys = np.convolve(y, np.ones(w) / w, mode="same") if w > 1 else y
    prom = rel_prominence * (ys.max() - ys.min())
    idx, props = _find_peaks(ys, prominence=prom, width=1, rel_height=0.5)
    out = []
    for k, i in enumerate(idx):
        if not (0 < i < len(x) - 1) or x[i] < x[0] + 0.3:
            continue
        a0, a1, a2 = ys[i - 1], ys[i], ys[i + 1]
        den = a0 - 2 * a1 + a2
        xc = x[i] + (0.5 * (a0 - a2) / den if den else 0) * step
        fwhm = props["widths"][k] * step
        base = ys[i] - props["prominences"][k]
        lo, hi = max(0, int(i - 1.5 * props["widths"][k])), min(len(x), int(i + 1.5 * props["widths"][k]) + 1)
        area = float(np.clip(ys[lo:hi] - base, 0, None).sum() * step)
        out.append(Peak(float(xc), float(d_spacing(xc, lam)), float(fwhm), float(ys[i]), area,
                        scherrer(xc, fwhm, lam, instrumental)))
    return out


def _inv_d2(lattice, h, k, l, a, c):
    inplane = 4 / 3 * (h*h + h*k + k*k) / a**2 if lattice == "hex" else (h*h + k*k) / a**2
    return inplane + (l*l / c**2 if c else 0.0)


def index_peaks(peaks, lattice="hex", lam=LAM_CU):
    """Index a peak list on a hexagonal/tetragonal lattice. Returns (a, c, [(peak, hkl, rel_err)])."""
    low = sorted((p for p in peaks if p.tt < 18), key=lambda p: p.tt)
    if not low:
        return None
    d1 = low[0].d
    a = 2 * d1 / math.sqrt(3) if lattice == "hex" else d1
    stack = sorted((p for p in peaks if 18 < p.tt < 32), key=lambda p: -p.area)
    c = stack[0].d if stack else 0.0
    cand = [(h, k, l) for h in range(10) for k in range(h + 1) for l in range(3 if c else 1) if h or k or l]
    assign = []
    for _ in range(4):
        assign = []
        for p in peaks:
            errs = [abs(1 / math.sqrt(_inv_d2(lattice, *m, a, c)) - p.d) / p.d for m in cand]
            j = int(np.argmin(errs)); assign.append((p, cand[j], errs[j]))
        f = [(4/3*(h*h+h*k+k*k) if lattice == "hex" else h*h+k*k, 1/p.d**2) for p, (h, k, l), e in assign if e < 0.02 and l == 0]
        if f:
            F, Y = np.array(f).T
            a = 1 / math.sqrt((F @ Y) / (F @ F))
    return a, c, assign


# ───────────────────────── Pawley ─────────────────────────
@dataclass
class Reflection:
    label: str
    tt: float
    d: float
    mult: int
    intensity: float
    esd: float
    fwhm: float
    is00l: bool


@dataclass
class PawleyResult:
    lattice: str
    a: float
    c: float | None
    zero: float
    esd: dict
    Rwp: float
    Rp: float
    DW: float
    esd_inflation: float
    size_inplane_nm: float | None
    size_stack_nm: float | None
    x: np.ndarray
    y: np.ndarray
    calc: np.ndarray
    background: np.ndarray
    reflections: list
    params: dict = field(default_factory=dict)
    unit_weights: bool = True

    def relative_intensities(self):
        ref = next((r for r in self.reflections if not r.is00l), self.reflections[0])
        return [(r.label, 100 * r.intensity / ref.intensity, 100 * r.esd / ref.intensity if np.isfinite(r.esd) else np.nan)
                for r in self.reflections]

    def summary(self):
        e = lambda k, n=3: f" ± {self.esd[k]:.{n}f}" if k in self.esd and np.isfinite(self.esd[k]) else ""
        lines = [f"Pawley ({'hexagonal' if self.lattice == 'hex' else 'tetragonal'}), "
                 f"{self.x[0]:.2f}–{self.x[-1]:.2f}°, {len(self.x)} points, {'unit' if self.unit_weights else '1/I'} weights",
                 f"  a = {self.a:.4f}{e('a', 4)} Å" + (f", c = {self.c:.4f}{e('c', 4)} Å" if self.c else ""),
                 f"  zero = {self.zero:+.4f}{e('zero', 4)}°   Rwp = {100*self.Rwp:.2f} %   Rp = {100*self.Rp:.2f} %"
                 f"   DW = {self.DW:.2f} (e.s.d. × {self.esd_inflation:.1f})",
                 f"  crystallite size: in-plane {self.size_inplane_nm:.1f} nm" +
                 (f", stacking {self.size_stack_nm:.1f} nm" if self.size_stack_nm else ""),
                 "  hkl      2θ(°)     d(Å)    I/Iref(%)   ±σ"]
        for (lab, pct, sg), r in zip(self.relative_intensities(), self.reflections):
            lines.append(f"  {lab:6s} {r.tt:8.3f} {r.d:8.3f} {pct:10.2f} {sg:7.2f}")
        return "\n".join(lines)


def _groups(lattice, a, c, use_c, lam, tmin, tmax):
    dmin = lam / (2 * math.sin(min(tmax, 170) * D2R / 2))
    H, L = int(a / dmin) + 2, (int(c / dmin) + 1 if use_c and c else 0)
    table = {}
    for h in range(-H, H + 1):
        for k in range(-H, H + 1):
            for l in range(0, L + 1):
                if not (h or k or l) or (l == 0 and (h < 0 or (h == 0 and k < 0))):
                    continue
                q = _inv_d2(lattice, h, k, l, a, c if use_c else 0)
                if 1 / math.sqrt(q) < dmin:
                    continue
                table.setdefault(round(q * 1e9), []).append((h, k, l))
    out = []
    for key, mem in table.items():
        rep = mem[0]
        d = 1 / math.sqrt(_inv_d2(lattice, *rep, a, c if use_c else 0))
        tt = float(two_theta(d, lam))
        if tmin - 1 <= tt <= tmax + 1:
            best = next((m for m in mem if m[0] >= m[1] >= 0 and m[2] >= 0), rep)
            label = "".join(str(v) if v >= 0 else f"-{-v}" for v in best)
            out.append(dict(rep=rep, mult=len(mem), label=label, is00l=all(m[0] == 0 and m[1] == 0 for m in mem), tt=tt))
    return sorted(out, key=lambda g: g["tt"])


def _tch(Hg, Hl):
    H = (Hg**5 + 2.69269*Hg**4*Hl + 2.42843*Hg**3*Hl**2 + 4.47163*Hg**2*Hl**3 + 0.07842*Hg*Hl**4 + Hl**5) ** 0.2
    r = Hl / H
    return H, np.clip(1.36603*r - 0.47719*r**2 + 0.11116*r**3, 0, 1)


class _Pawley:
    NAMES = ("a", "c", "zero", "G0", "X", "Xc", "asym", "bgt")

    def __init__(self, x, y, w, lattice, a0, c0, use_c, lam, zero, lp_profile, bg_exp, X0, fixed_I=None):
        self.x, self.y, self.sw = x, y, np.sqrt(w)
        self.lattice, self.lam, self.lp_profile, self.bg_exp = lattice, lam, lp_profile, bg_exp
        self.groups = _groups(lattice, a0, c0, use_c, lam, x[0], x[-1])
        self.has00l = bool(use_c and any(g["rep"][2] for g in self.groups))
        self.fixed_I = fixed_I
        c0 = c0 or 3.5
        # value, lower, upper, fixed
        self.P = {"a": [a0, 0.85*a0, 1.15*a0, False], "c": [c0, 0.85*c0, 1.15*c0, not self.has00l],
                  "zero": [0.0 if zero is None else zero, -0.5, 0.5, zero is not None],
                  "G0": [0.05, 0.003, 2.0, False], "X": [X0, 0.003, 6.0, False], "Xc": [2.0, 0.01, 15.0, not self.has00l],
                  "asym": [0.0, -0.6, 0.6, False], "bgt": [0.5, 0.05, 6.0, not bg_exp]}
        self.free = [k for k in self.NAMES if not self.P[k][3]]
        u = 2 * (x - x[0]) / (x[-1] - x[0]) - 1
        self.cheb = np.stack([np.ones_like(u), u, 2*u*u - 1, 4*u**3 - 3*u], 1)
        self.lpx = lp_factor(x) if lp_profile else None

    def theta(self, v=None):
        th = {k: self.P[k][0] for k in self.NAMES}
        if v is not None:
            th.update(zip(self.free, v))
        return th

    def peak_columns(self, th):
        x, cols, meta = self.x, [], []
        for g in self.groups:
            h, k, l = g["rep"]
            d = 1 / math.sqrt(_inv_d2(self.lattice, h, k, l, th["a"], th["c"] if self.has00l else 0))
            t = float(two_theta(d, self.lam)) + th["zero"]
            cphi = min(1.0, abs(l) * d / th["c"]) if self.has00l else 0.0
            Xe = math.sqrt(th["X"]**2 * (1 - cphi**2) + th["Xc"]**2 * cphi**2)
            H, eta = _tch(th["G0"], Xe / math.cos(t * D2R / 2))
            u = x - t
            Hs = np.where(u < 0, H * (1 - th["asym"]), H * (1 + th["asym"]))
            z = 4 * u * u / Hs**2
            v = eta / (1 + z) + (1 - eta) * np.exp(-math.log(2) * z)
            v[np.abs(u) > 14 * H] = 0.0
            if self.lpx is not None:
                v = v * self.lpx / lp_factor(t)
                area = _trapz(v, x)
            else:
                area = (eta * math.pi / 2 + (1 - eta) * math.sqrt(math.pi / (4 * math.log(2)))) * H
            cols.append(v / (area if area > 0 else 1))
            meta.append((t, d, H))
        return np.array(cols).T, meta

    def design(self, th):
        P, meta = self.peak_columns(th)
        if self.fixed_I is not None:
            P = (P @ self.fixed_I)[:, None]
        bg = [np.exp(-(self.x - self.x[0]) / th["bgt"])[:, None]] if self.bg_exp else []
        A = np.hstack([P] + bg + [self.cheb])
        npk = P.shape[1]
        lb = np.r_[np.zeros(npk + len(bg)), -np.inf * np.ones(4)]
        return A, lb, meta, npk, len(bg)

    def solve_linear(self, th):
        A, lb, meta, npk, nbg = self.design(th)
        sol = lsq_linear(A * self.sw[:, None], self.y * self.sw, bounds=(lb, np.inf), method="bvls")
        return A, sol.x, meta, npk, nbg

    def resid(self, v):
        th = self.theta(v)
        A, beta, *_ = self.solve_linear(th)
        return self.sw * (self.y - A @ beta)


def pawley(x, y, lattice="hex", a0=None, c0=None, use_c="auto", tt_range=None, lam=LAM_CU,
           zero=None, lp_profile=True, bg_exp=True, weights="auto", fixed_intensities=None, max_nfev=200):
    """Pawley refinement on a hexagonal ('hex') or tetragonal ('tet') cell.

    zero=None refines the zero shift; a number fixes it. use_c: 'auto' | True | False (include 00l/hkl
    reflections). fixed_intensities: array aligned with the reflection groups → only a scale is refined
    (structure-model fit). Returns PawleyResult.
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    if tt_range:
        m = (x >= tt_range[0]) & (x <= tt_range[1]); x, y = x[m], y[m]
    if len(x) > 1200:                                   # bin: COF peaks are ≥ 0.2° wide
        f = int(math.ceil(len(x) / 1200)); n = len(x) // f * f
        x, y = x[:n].reshape(-1, f).mean(1), y[:n].reshape(-1, f).mean(1)
    pk = find_peaks(x, y, lam=lam)
    if a0 is None or (c0 is None and use_c != False):
        ix = index_peaks(pk, lattice, lam)
        if ix:
            a0 = a0 or ix[0]; c0 = c0 if c0 is not None else (ix[1] or None)
    if a0 is None:
        raise ValueError("cannot guess a0; pass a0=")
    uc = (c0 is not None and c0 > 0 and float(two_theta(c0, lam)) <= x[-1]) if use_c == "auto" else bool(use_c)
    unit = weights == "unit" or (weights == "auto" and (y < 0).any())
    w = np.ones_like(y) if unit else 1 / np.maximum(np.abs(y), 1)
    X0 = max(0.05, 0.7 * pk[0].fwhm) if pk else 0.4
    pr = _Pawley(x, y, w, lattice, a0, c0 or 0, uc, lam, zero, lp_profile and x[0] < 8, bg_exp and x[0] < 6, X0,
                 None if fixed_intensities is None else np.asarray(fixed_intensities, float))
    v0 = [pr.P[k][0] for k in pr.free]
    lo = [pr.P[k][1] for k in pr.free]; hi = [pr.P[k][2] for k in pr.free]
    sol = least_squares(pr.resid, v0, bounds=(lo, hi), x_scale="jac", max_nfev=max_nfev, diff_step=1e-4)
    th = pr.theta(sol.x)
    A, beta, meta, npk, nbg = pr.solve_linear(th)
    calc = A @ beta
    r = pr.sw * (y - calc)
    rss = float(r @ r)
    Rwp = math.sqrt(rss / float((pr.sw * y) @ (pr.sw * y)))
    Rp = float(np.abs(y - calc).sum() / np.abs(y).sum())
    DW = float(np.sum(np.diff(r) ** 2) / rss)
    infl = math.sqrt(2 / max(DW, 1e-3)) if DW < 2 else 1.0
    # covariance: ∂model/∂θ with β fixed, plus the active linear columns
    J, names = [], []
    for k in pr.free:
        i = pr.free.index(k)
        if sol.x[i] <= pr.P[k][1] + 1e-9 or sol.x[i] >= pr.P[k][2] - 1e-9:
            continue
        h = 1e-4 * max(abs(sol.x[i]), 0.05)
        v2 = sol.x.copy(); v2[i] += h
        A2 = pr.design(pr.theta(v2))[0]
        J.append(pr.sw * (A2 @ beta - calc) / h); names.append(k)
    active = [j for j in range(len(beta)) if beta[j] > 0 or j >= npk + nbg]
    for j in active:
        J.append(pr.sw * A[:, j]); names.append(("beta", j))
    J = np.array(J).T
    s2 = rss / max(1, len(y) - J.shape[1])
    try:
        cov = np.linalg.inv(J.T @ J) * s2
        e = np.sqrt(np.abs(np.diag(cov))) * infl
    except np.linalg.LinAlgError:
        e = np.full(J.shape[1], np.nan)
    esd = {n: e[i] for i, n in enumerate(names) if isinstance(n, str)}
    esd_beta = {n[1]: e[i] for i, n in enumerate(names) if not isinstance(n, str)}
    refl = []
    for j, (g, (t, d, H)) in enumerate(zip(pr.groups, meta)):
        if fixed_intensities is None:
            I, sI = float(beta[j]), float(esd_beta.get(j, np.nan))
        else:
            I, sI = float(pr.fixed_I[j] * beta[0]), np.nan
        refl.append(Reflection(g["label"], t, d, g["mult"], I, sI, H, g["is00l"]))
    bg = A[:, npk:] @ beta[npk:]
    size = lambda rr: scherrer(rr.tt, rr.fwhm, lam) if rr else None
    r_in = next((q for q in refl if not q.is00l and q.intensity > 0), None)
    r_st = next((q for q in refl if q.is00l and q.intensity > 0), None)
    return PawleyResult(lattice, th["a"], th["c"] if pr.has00l else None, th["zero"], esd, Rwp, Rp, DW, infl,
                        size(r_in), size(r_st) if pr.has00l else None, x, y, calc, bg, refl,
                        params={k: th[k] for k in pr.NAMES}, unit_weights=unit)


def zero_scan(x, y, zeros=(-0.15, -0.10, -0.05, 0.0, 0.05, 0.10, 0.15), **kw):
    """Refine with the zero shift fixed at each value; returns [(zero, a, c, Rwp)]."""
    out = []
    for z in zeros:
        f = pawley(x, y, zero=z, **kw)
        out.append((z, f.a, f.c, f.Rwp))
    return out


def model_percentages(fit: PawleyResult, model_d, model_I, tol=0.012):
    """Map model reflections (d, I in the model's own cell) onto the Pawley reflections.

    The model d-spacings are scaled so that its strongest low-angle line matches the first
    in-plane Pawley reflection, which absorbs a cell mismatch. Returns % of that reflection.
    """
    ref = next((r for r in fit.reflections if not r.is00l), fit.reflections[0])
    model_d, model_I = np.asarray(model_d, float), np.asarray(model_I, float)
    low = np.argsort(-model_d)
    low = [i for i in low if model_d[i] > 0.7 * ref.d][:4] or list(np.argsort(-model_d)[:1])
    strong = max(low, key=lambda i: model_I[i])
    sc = ref.d / model_d[strong]
    out = np.array([model_I[np.abs(model_d * sc - r.d) / r.d < tol].sum() for r in fit.reflections])
    base = out[fit.reflections.index(ref)] or 1.0
    return 100 * out / base, sc


def compare(fit: PawleyResult, model_pct, name="model", nsigma=3.0):
    """Reflections whose Pawley intensity differs from the model by > nsigma."""
    ref = next((r for r in fit.reflections if not r.is00l), fit.reflections[0])
    bad = []
    for r, pm in zip(fit.reflections, model_pct):
        if r is ref or not np.isfinite(r.esd):
            continue
        pe, sg = 100 * r.intensity / ref.intensity, 100 * r.esd / ref.intensity
        if pe < 0.5 and pm < 0.5:
            continue
        z = (pe - pm) / max(sg, 0.02 * max(pe, pm))
        if abs(z) > nsigma:
            bad.append((r.label, pe, pm, z))
    return bad


# ───────────────────────── report / CLI ─────────────────────────
def analyse(path, exp=None, sim=None, cif=None, lattice="hex", tt_range=None, zero=None, scan=False,
            lam=LAM_CU, B=2.0, plot=False, out_prefix=None, verbose=True):
    """Full workflow used by the CLI; returns a dict with all results."""
    series = load_xrd(path)
    pick = lambda i, want_sim: series[i] if i is not None else next((s for s in series if s.is_sim == want_sim), None)
    E = pick(exp, False) or series[0]
    S = pick(sim, True) if sim != -1 else None
    if S is E:
        S = None
    say = print if verbose else (lambda *a, **k: None)
    say(f"experiment: {E.name}  ({len(E.x)} pts, {E.x[0]:.2f}–{E.x[-1]:.2f}°)" + ("  [background-subtracted]" if E.has_neg else ""))
    peaks = find_peaks(E.x, E.y, lam=lam)
    say("peaks (2θ, d, FWHM, size nm):  " + "; ".join(f"{p.tt:.2f}° {p.d:.2f} Å {p.fwhm:.2f}° {p.size_nm:.1f}" for p in peaks))
    if tt_range is None:
        first = peaks[0] if peaks else None
        lo = max(E.x[0], first.tt - 2.5 * first.fwhm) if first else E.x[0]
        sharp = [p.tt for p in peaks if p.tt < 18 and p.fwhm < 1.5]
        hi = min(E.x[-1], max(12.0, (sharp[-1] if sharp else 9) + 3))
        tt_range = (lo, hi)
    res = {"experiment": E, "peaks": peaks, "tt_range": tt_range}
    simcell = None
    if S is not None:
        spk = find_peaks(S.x, S.y, rel_prominence=0.001, smooth_deg=0, lam=lam)
        ix = index_peaks(spk, lattice, lam)
        if ix:
            simcell = ix[:2]
            say(f"simulation: {S.name}  cell a = {ix[0]:.3f} Å" + (f", c = {ix[1]:.3f} Å" if ix[1] else ""))
        res.update(sim=S, sim_peaks=spk, sim_cell=simcell)
    fit = pawley(E.x, E.y, lattice=lattice, tt_range=tt_range, zero=zero, lam=lam)
    say("\n" + fit.summary())
    res["pawley"] = fit
    if simcell:
        dA = 100 * (fit.a - simcell[0]) / simcell[0]
        say(f"\ncell vs simulation: a_exp − a_sim = {dA:+.1f} %" + ("  (> 2 %: check the simulated cell / zero shift)" if abs(dA) > 2 else ""))
    models = {}
    if S is not None and res.get("sim_peaks"):
        spk = res["sim_peaks"]
        models["simulation"] = model_percentages(fit, [p.d for p in spk], [p.area for p in spk])[0]
    if cif:
        st = read_cif(cif)
        _, d, tt, I = structure_reflections(st, lam, min(60.0, float(two_theta(min(r.d for r in fit.reflections) * 0.8, lam))), B)
        models["CIF"] = model_percentages(fit, d, I)[0]
        res["structure"] = st
    for name, pct in models.items():
        say(f"\n{name} vs Pawley (I/Iref %):  " + "  ".join(f"{r.label} {100*r.intensity/fit.reflections[0].intensity:.1f}/{p:.1f}" for r, p in zip(fit.reflections, pct)))
        bad = compare(fit, pct, name)
        say("  > 3σ: " + (", ".join(f"{b[0]} ({b[1]:.1f} vs {b[2]:.1f} %)" for b in bad) if bad else "none"))
        mf = pawley(E.x, E.y, lattice=lattice, tt_range=tt_range, zero=zero, lam=lam, fixed_intensities=pct)
        say(f"  fit with {name} intensities fixed: Rwp = {100*mf.Rwp:.2f} % (Pawley {100*fit.Rwp:.2f} %)")
        res.setdefault("models", {})[name] = (pct, bad, mf)
    if scan:
        sc = zero_scan(E.x, E.y, lattice=lattice, tt_range=tt_range, lam=lam)
        span = max(q[1] for q in sc) - min(q[1] for q in sc)
        say("\nzero scan:  " + "  ".join(f"{q[0]:+.2f}° → a {q[1]:.2f} Å (Rwp {100*q[3]:.2f} %)" for q in sc))
        if span > 0.3:
            say(f"  a varies by {span:.2f} Å over ±0.15° zero shift → fix the zero with an internal standard (e.g. Si).")
        res["zero_scan"] = sc
    if fit.DW < 0.5:
        say(f"\nnote: Durbin–Watson {fit.DW:.2f} — residuals strongly correlated (smoothed / background-subtracted data); e.s.d.s ×{fit.esd_inflation:.1f}.")
    if out_prefix:
        _write_outputs(res, out_prefix, plot)
        say(f"\nwritten: {out_prefix}_pawley_fit.csv, {out_prefix}_reflections.csv" + (f", {out_prefix}_pawley.png" if plot else ""))
    return res


def _write_outputs(res, prefix, plot):
    fit = res["pawley"]
    with open(f"{prefix}_pawley_fit.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh); w.writerow(["2theta_deg", "obs", "calc", "background", "diff"])
        for row in zip(fit.x, fit.y, fit.calc, fit.background, fit.y - fit.calc):
            w.writerow([f"{v:.6g}" for v in row])
    with open(f"{prefix}_reflections.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        names = list(res.get("models", {}))
        w.writerow(["hkl", "2theta_deg", "d_A", "multiplicity", "I", "I_esd", "I_rel_pct", "I_rel_esd_pct"] + [f"{n}_pct" for n in names])
        for k, (r, (lab, pct, sg)) in enumerate(zip(fit.reflections, fit.relative_intensities())):
            w.writerow([lab, f"{r.tt:.4f}", f"{r.d:.4f}", r.mult, f"{r.intensity:.6g}", f"{r.esd:.4g}", f"{pct:.3f}", f"{sg:.3f}"] +
                       [f"{res['models'][n][0][k]:.3f}" for n in names])
    if plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(fit.x, fit.y, "k.", ms=2.5, label="observed")
        ax.plot(fit.x, fit.calc, color="#c4561f", lw=1.3, label=f"Pawley (Rwp {100*fit.Rwp:.2f} %)")
        ax.plot(fit.x, fit.background, ":", color="gray", lw=1, label="background")
        off = min(fit.y.min(), 0) - 0.08 * np.ptp(fit.y)
        ax.plot(fit.x, fit.y - fit.calc + off, color="gray", lw=0.8, label="difference")
        for r in fit.reflections:
            ax.plot([r.tt, r.tt], [off * 0.55, off * 0.35], color="#1f6f8b", lw=1)
        ax.set_xlabel("2θ (°)"); ax.set_ylabel("intensity"); ax.legend(frameon=False, fontsize=8)
        ax.set_title(f"a = {fit.a:.3f} Å" + (f", c = {fit.c:.3f} Å" if fit.c else "") + f", zero {fit.zero:+.3f}°", loc="left", fontsize=10)
        fig.tight_layout(); fig.savefig(f"{prefix}_pawley.png", dpi=180); plt.close(fig)


def _synthetic(a=30.0, c=3.45, lam=LAM_CU, fwhm=0.6, step=0.02, noise=40.0, seed=1):
    """Honeycomb (hcb) COF of uniform line density — used by the self-test."""
    rng = np.random.default_rng(seed)
    nodeA = np.array([1/3, 2/3]); nbs = np.array([[2/3, 1/3], [-1/3, 1/3], [2/3, 4/3]])
    pts = np.vstack([nodeA + (nb - nodeA) * t for nb in nbs for t in np.linspace(0, 1, 40, endpoint=False)])
    st = Structure((a, a, c, 90, 90, 120), np.array(["C"] * len(pts)), np.c_[pts % 1, np.zeros(len(pts))], np.ones(len(pts)))
    x = np.arange(2, 35 + 1e-9, step)
    y = simulate_pattern(st, x, lam, fwhm, 0.5, B=3)
    y = 8000 * y / y.max() + 9000 * np.exp(-(x - 2) / 0.5) + noise * rng.standard_normal(len(x))
    return x, y, st


def selftest():
    x, y, st = _synthetic()
    fit = pawley(x, y, lattice="hex", tt_range=(2.3, 14), zero=0.0)
    _, d, tt, I = structure_reflections(st, LAM_CU, 16, B=3)
    pct, _ = model_percentages(fit, d, I)
    rel = fit.relative_intensities()
    ok_a = abs(fit.a - 30.0) < 0.05
    dev = [abs(p - m) for (lab, p, s), m in zip(rel, pct) if m > 3]
    ok_I = max(dev) < 1.5 if dev else False
    print(fit.summary())
    print(f"\nself-test: a = {fit.a:.4f} Å (true 30.0000) {'OK' if ok_a else 'FAIL'}; "
          f"max |ΔI| of strong reflections = {max(dev):.2f} % {'OK' if ok_I else 'FAIL'}")
    p = read_cif("data_t\n_cell_length_a 30\n_cell_length_b 30\n_cell_length_c 3.5\n_cell_angle_alpha 90\n_cell_angle_beta 90\n"
                 "_cell_angle_gamma 120\nloop_\n_space_group_symop_operation_xyz\n'x,y,z'\n'-y,x-y,z'\n'-x+y,-x,z'\n'-x,-y,z'\n'y,-x+y,z'\n'x-y,x,z'\n"
                 "loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_occupancy\n"
                 "C1 C 0.3333 0.6667 0 1\nN1 N 0.1000(2) 0.2 0.5 0.5\nS1 S 1/2 0 0 1\n")
    counts = {e: int((p.elements == e).sum()) for e in "CNS"}
    ok_cif = counts == {"C": 2, "N": 6, "S": 3}
    print(f"CIF symmetry expansion (P6): {counts} {'OK' if ok_cif else 'FAIL'}")
    return ok_a and ok_I and ok_cif


def main(argv=None):
    ap = argparse.ArgumentParser(description="COF powder-XRD analysis: peaks, Pawley, simulation/CIF comparison.")
    ap.add_argument("data", nargs="?", help="CSV/TXT/XY file")
    ap.add_argument("--exp", type=int, help="index of the experimental series (default: first non-simulation)")
    ap.add_argument("--sim", type=int, help="index of the simulated series (-1: none; default: auto)")
    ap.add_argument("--cif", help="structure model (CIF) for model intensities")
    ap.add_argument("--lattice", choices=["hex", "tet"], default="hex")
    ap.add_argument("--range", nargs=2, type=float, metavar=("LO", "HI"), help="2θ range for the refinement")
    ap.add_argument("--zero", type=float, help="fix the zero shift (default: refined)")
    ap.add_argument("--scan", action="store_true", help="zero-shift scan (a–zero correlation)")
    ap.add_argument("--lam", type=float, default=LAM_CU, help="wavelength in Å (default Cu Kα1)")
    ap.add_argument("--B", type=float, default=2.0, help="isotropic B (Å²) for CIF intensities")
    ap.add_argument("--out", help="output prefix for CSV (and PNG with --plot)")
    ap.add_argument("--plot", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args(argv)
    if a.selftest:
        return 0 if selftest() else 1
    if not a.data:
        ap.error("data file required (or --selftest)")
    for i, s in enumerate(load_xrd(a.data)):
        print(f"[{i}] {s.name}: {len(s.x)} pts {s.x[0]:.2f}–{s.x[-1]:.2f}°{'  (simulation)' if s.is_sim else ''}")
    out = a.out or (re.sub(r"\.[^.]+$", "", a.data) if a.plot else None)
    analyse(a.data, a.exp, a.sim, a.cif, a.lattice, tuple(a.range) if a.range else None, a.zero, a.scan, a.lam, a.B, a.plot, out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
