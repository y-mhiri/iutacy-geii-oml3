"""TP2 -- Série de Fourier, forme complexe. Fonctions fournies, à ne pas
modifier.

Signal : un signal T-périodique, défini par sa fonction f (à valeurs
scalaires), sa période T, et les points de [0, T) où f n'est pas dérivable
(sauts ou changements de pente). Ces points servent uniquement à guider
l'intégration numérique -- ils n'ont pas besoin d'être exhaustifs.
"""

from dataclasses import dataclass, field
from typing import Callable

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad


@dataclass
class Signal:
    f: Callable[[float], float]
    T: float
    name: str
    breakpoints: list = field(default_factory=list)


# ─── Signaux ────────────────────────────────────────────────────────

def make_square(T=1.0, A=1.0):
    """Signal carré k(t) : +A sur [0,T/2), -A sur [T/2,T). Le même signal
    que dans le TP1 (a_n=0, b_n=4A/(n*pi) pour n impair)."""
    def f(t):
        return A if (t % T) < T / 2 else -A
    return Signal(f=f, T=T, name="carré", breakpoints=[0.0, T / 2])


def make_levels(N=4, T=1.0, A=1.0, seed=0):
    """Signal en escalier à N paliers de largeur égale T/N, de hauteurs
    v_0,...,v_{N-1} tirées indépendamment dans [-A, A] (sans lien de
    symétrie particulier entre les paliers)."""
    rng = np.random.default_rng(seed)
    values = rng.uniform(-A, A, N)

    def f(t):
        k = min(int((t % T) // (T / N)), N - 1)
        return float(values[k])

    breakpoints = [k * T / N for k in range(N)]
    return Signal(f=f, T=T, name=f"{N} paliers", breakpoints=breakpoints)


def shift_signal(signal, phi, name=None):
    """g(t) = signal.f(t + phi), même période. Les points de rupture sont
    translatés en conséquence."""
    def g(t):
        return signal.f(t + phi)
    bp = sorted((p - phi) % signal.T for p in signal.breakpoints)
    return Signal(f=g, T=signal.T, name=name or f"{signal.name} décalé",
                  breakpoints=bp)


def combine_signals(sig_a, sig_b, op, name=None):
    """h(t) = op(sig_a.f(t), sig_b.f(t)) -- par exemple op = lambda a, b: a + b,
    ou op = lambda a, b: a - b. Les deux signaux doivent avoir la même
    période."""
    assert sig_a.T == sig_b.T
    T = sig_a.T

    def h(t):
        return op(sig_a.f(t), sig_b.f(t))
    bp = sorted(set([p % T for p in sig_a.breakpoints] +
                     [p % T for p in sig_b.breakpoints]))
    return Signal(f=h, T=T, name=name or "combinaison", breakpoints=bp)


# ─── Coefficients de Fourier ────────────────────────────────────────

def fourier_coeffs(signal, N):
    """Calcule a0, (a_n)_1..N, (b_n)_1..N par intégration numérique des
    formules du cours (scipy.integrate.quad), sans passer par un
    échantillonnage régulier ni par une FFT. Identique au TP1."""
    T, f = signal.T, signal.f
    w0 = 2 * np.pi / T
    interior = sorted(p for p in signal.breakpoints if 0 < p < T)
    bounds = [0.0] + interior + [T]
    pts = interior or None

    def oscillatory_integral(weight, wvar):
        return sum(
            quad(f, lo, hi, weight=weight, wvar=wvar, limit=200)[0]
            for lo, hi in zip(bounds[:-1], bounds[1:])
        )

    a0 = quad(f, 0, T, points=pts, limit=200)[0] / T
    an = np.zeros(N)
    bn = np.zeros(N)
    for n in range(1, N + 1):
        an[n - 1] = 2 / T * oscillatory_integral("cos", n * w0)
        bn[n - 1] = 2 / T * oscillatory_integral("sin", n * w0)
    return a0, an, bn


def to_complex_coeffs(a0, an, bn):
    """Construit les coefficients complexes (c_n) pour n = -N..N à partir
    de a0, (a_n), (b_n) (coefficients réels, cf. TD1), via
        c_0 = a0,   c_n = (a_n - j b_n)/2  (n>0),   c_{-n} = conj(c_n).
    Retourne un dict {n: c_n}.
    """
    N = len(an)
    c = {0: complex(a0, 0.0)}
    for n in range(1, N + 1):
        c[n] = (an[n - 1] - 1j * bn[n - 1]) / 2
        c[-n] = np.conj(c[n])
    return c


def reconstruct_complex(t, T, c, orders=None):
    """Évalue sum_n c_n exp(j n w0 t) sur le vecteur t, pour n dans
    `orders` (par défaut, toutes les clés de c). Retourne un tableau
    complexe : pour un signal réel, la partie imaginaire doit rester
    numériquement négligeable si `orders` respecte la symétrie
    c_{-n} = conj(c_n).
    """
    w0 = 2 * np.pi / T
    t = np.asarray(t, dtype=float)
    orders = c.keys() if orders is None else orders
    s = np.zeros_like(t, dtype=complex)
    for n in orders:
        s = s + c[n] * np.exp(1j * n * w0 * t)
    return s


# ─── Tracés ──────────────────────────────────────────────────────────

def plot_signal(signal, t=None, tmin=None, tmax=None, npts=2000,
                 title=None, ax=None, xlim=None, ylim=None):
    """Trace signal.f sur une grille temporelle.

    Deux façons de préciser la grille, au choix :
    - passer directement un tableau de temps via `t` ;
    - ou donner `tmin`, `tmax` (et éventuellement `npts`, la densité de
      points) : la grille `np.linspace(tmin, tmax, npts)` est alors
      construite pour vous.

    Pour "zoomer" sur une zone, resserrer (tmin, tmax) -- en augmentant
    `npts` si la zone observée est petite -- plutôt que recadrer après
    coup avec `xlim` sur une grille large et grossière.
    """
    ax = ax or plt.gca()
    if t is None:
        if tmin is None or tmax is None:
            raise ValueError(
                "plot_signal : fournir soit `t`, soit `tmin` et `tmax`."
            )
        t = np.linspace(tmin, tmax, npts)
    y = [signal.f(ti) for ti in t]
    ax.plot(t, y)
    ax.set_xlabel("t")
    ax.set_title(title or signal.name)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    return ax


def plot_spectrum_complex(c, title=None, ax=None, xlim=None, ylim=None,
                           freq=False, T=None):
    """Trace |c_n| en fonction du rang n (spectre bilatéral, n négatifs et
    positifs), ou de la fréquence n/T si freq=True (fournir T dans ce cas).
    """
    ax = ax or plt.gca()
    ns = np.array(sorted(c.keys()))
    amp = np.array([abs(c[n]) for n in ns])
    x = ns / T if freq else ns
    ax.stem(x, amp)
    ax.set_xlabel("fréquence n/T" if freq else "rang n")
    ax.set_ylabel("|c_n|")
    ax.set_title(title or "Spectre bilatéral")
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    return ax


def plot_reconstruction(t, reconstructions, original=None, title=None,
                         ax=None, xlim=None, ylim=None, show_imag=False):
    """reconstructions : dict {étiquette: valeurs} (réelles ou complexes)
    à superposer. Si une reconstruction est complexe, sa partie réelle est
    tracée (et, si show_imag=True, sa partie imaginaire aussi, en
    pointillés) : utile pour observer ce qui se passe quand la symétrie
    c_{-n} = conj(c_n) est brisée.
    """
    ax = ax or plt.gca()
    if original is not None:
        ax.plot(t, original, color="black", lw=2, label="signal original",
                 alpha=0.6)
    for label, y in reconstructions.items():
        y = np.asarray(y)
        if np.iscomplexobj(y):
            ax.plot(t, y.real, label=f"{label} (Re)")
            if show_imag:
                ax.plot(t, y.imag, "--", label=f"{label} (Im)")
        else:
            ax.plot(t, y, label=label)
    ax.set_xlabel("t")
    ax.legend()
    ax.set_title(title or "Reconstruction")
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    return ax
