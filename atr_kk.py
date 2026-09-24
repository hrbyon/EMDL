"""ATR-FTIR Kramers-Kronig correction following Kim, Cho, Kwak, Anal. Chem. 2024, 96, 15924.

Pipeline (unpolarized, single reflection, theta = 45 deg):
  R (measured, 0-1) -> Rs from R = (Rs + Rs^2)/2  (Abeles: Rp = Rs^2 at 45 deg)
  -> phase Theta_s via KK on ln sqrt(Rs) (Maclaurin integration, Ohta & Ishida 1988)
  -> complex n_hat from r_s = sqrt(Rs) exp(i Theta_s) (Fresnel, s-pol)
  -> corrected absorbance A = 4*pi*k*nu*l/ln10  (= 1.74*pi*k*nu*l, eq 16)

Convention: n_hat = n + i k (k > 0), xi = n_hat cos(phi) = sqrt(n_hat^2 - n0^2 sin^2 theta), Im xi > 0,
r_s = (n0 cos theta - xi) / (n0 cos theta + xi).
"""
import numpy as np

LN10 = np.log(10.0)


# ---------- forward model (used for validation / simulation) ----------
def lorentz_eps(nu, n_inf, oscillators):
    """Complex dielectric function; oscillators = [(nu0, strength S, gamma), ...]."""
    eps = np.full(nu.shape, n_inf ** 2, dtype=complex)
    for nu0, S, g in oscillators:
        eps += S * nu0 ** 2 / (nu0 ** 2 - nu ** 2 - 1j * g * nu)
    return eps


def atr_forward(n_hat, n0=2.4, theta_deg=45.0):
    th = np.deg2rad(theta_deg)
    a = n0 * np.cos(th)
    xi = np.sqrt(n_hat ** 2 - (n0 * np.sin(th)) ** 2 + 0j)
    xi = np.where(xi.imag < 0, -xi, xi)
    rs = (a - xi) / (a + xi)
    Rs = np.abs(rs) ** 2
    return (Rs + Rs ** 2) / 2, Rs, rs


# ---------- inverse (the correction) ----------
def rs_from_unpolarized(R):
    R = np.clip(R, 1e-6, 1.0)
    return (-1.0 + np.sqrt(1.0 + 8.0 * R)) / 2.0


def kk_maclaurin(nu, f):
    """P-integral  (2 nu_i / pi) * P int f(nu') / (nu'^2 - nu_i^2) dnu'  on an even grid (Maclaurin: alternate points)."""
    N = len(nu)
    h = nu[1] - nu[0]
    out = np.empty(N)
    idx = np.arange(N)
    for i in range(N):
        j = idx[(idx - i) % 2 == 1]
        out[i] = (2.0 / np.pi) * 2.0 * h * nu[i] * np.sum(f[j] / (nu[j] ** 2 - nu[i] ** 2))
    return out


def theta_inf(n_inf, n0=2.4, theta_deg=45.0):
    th = np.deg2rad(theta_deg)
    return -2.0 * np.arctan(np.sqrt((n0 * np.sin(th)) ** 2 - n_inf ** 2) / (n0 * np.cos(th)))


def atr_kk_correct(nu, R, n_inf, n0=2.4, theta_deg=45.0, path_um=1.0):
    """nu ascending, evenly spaced (cm-1); R = measured reflectance (0-1, unpolarized).
    Returns dict with n, k, corrected absorbance A_corr, phase."""
    th = np.deg2rad(theta_deg)
    Rs = rs_from_unpolarized(R)
    lnr = 0.5 * np.log(Rs)
    phase = theta_inf(n_inf, n0, theta_deg) - kk_maclaurin(nu, lnr)
    rs = np.sqrt(Rs) * np.exp(1j * phase)
    a = n0 * np.cos(th)
    xi = a * (1 - rs) / (1 + rs)
    nhat = np.sqrt(xi ** 2 + (n0 * np.sin(th)) ** 2)
    nhat = np.where(nhat.imag < 0, -nhat, nhat)
    n, k = nhat.real, np.clip(nhat.imag, 0, None)
    A = 4 * np.pi * k * nu * (path_um * 1e-4) / LN10
    return dict(n=n, k=k, A_corr=A, phase=phase, Rs=Rs)


def partial_correct(nu, R, n_sample=1.3, n0=2.4, theta_deg=45.0, nu_ref=1000.0):
    """'Partially corrected' ATR: pATR normalized by penetration depth for constant n (dp ~ 1/nu)."""
    pATR = -np.log10(np.clip(R, 1e-6, None))
    return pATR * nu / nu_ref


if __name__ == "__main__":
    # validation: simulate a strongly absorbing liquid, invert, compare to the true k
    nu = np.arange(400.0, 4000.0, 0.964)
    osc = [(1180, 0.020, 12), (1380, 0.010, 14), (1085, 0.012, 15), (740, 0.004, 15)]
    for ninf in (1.40, 1.60):
        nhat = np.sqrt(lorentz_eps(nu, ninf, osc))
        R, _, _ = atr_forward(nhat)
        out = atr_kk_correct(nu, R, ninf)
        A_true = 4 * np.pi * nhat.imag * nu * 1e-4 / LN10
        pATR = -np.log10(R)
        for c in (1180, 740):
            m = (nu > c - 40) & (nu < c + 40)
            print(f"n_inf={ninf} band {c}: true peak {nu[m][A_true[m].argmax()]:.1f}, "
                  f"raw pATR peak {nu[m][pATR[m].argmax()]:.1f}, corrected {nu[m][out['A_corr'][m].argmax()]:.1f}; "
                  f"max|k err|/kmax = {np.abs(out['k'][m]-nhat.imag[m]).max()/nhat.imag[m].max():.3f}")
