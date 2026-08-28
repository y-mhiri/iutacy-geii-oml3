"""TP1 -- Série de Fourier. Fonctions fournies, à ne pas modifier.

Signal : un signal T-périodique, défini par sa fonction f (à valeurs
scalaires), sa période T, et les points de [0, T) où f n'est pas dérivable
(sauts ou changements de pente). Ces points servent uniquement à guider
l'intégration numérique -- ils n'ont pas besoin d'être exhaustifs.
"""

import warnings
from dataclasses import dataclass, field
from typing import Callable

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from scipy.integrate import IntegrationWarning


@dataclass
class Signal:
    f: Callable[[float], float]
    T: float
    name: str
    breakpoints: list = field(default_factory=list)


# ─── Signaux ────────────────────────────────────────────────────────

def make_square(T=1.0, A=1.0):
    """Signal carré : +A sur [0,T/2), -A sur [T/2,T)."""
    def f(t):
        return A if (t % T) < T / 2 else -A
    return Signal(f=f, T=T, name="carré", breakpoints=[0.0, T / 2])


def make_sawtooth(T=1.0):
    """Signal dents de scie : f(t) = t sur [0,T). Une seule discontinuité
    par période, exactement à la jonction t=0/T (contrairement au carré,
    qui en a une au milieu et une à la jonction)."""
    def f(t):
        return t % T
    return Signal(f=f, T=T, name="dents de scie", breakpoints=[0.0])


def make_random_walk(T=1.0, n_points=25, amplitude=1.0, seed=42):
    """Signal en dents irrégulières : n_points sommets tirés au hasard,
    reliés par des segments de droite, refermés en fin de période (donc
    continu). Nombre de sommets fixé -> nombre fini d'extrema : le signal
    respecte les conditions de Dirichlet malgré son allure erratique.
    """
    rng = np.random.default_rng(seed)
    xs = np.sort(rng.uniform(0, T, n_points))
    xs[0] = 0.0
    ys = amplitude * rng.uniform(-1, 1, n_points)
    ys[-1] = ys[0]  # boucle la période : pas de saut en t=0/T

    def f(t):
        return float(np.interp(t % T, xs, ys))
    return Signal(f=f, T=T, name="marche aléatoire", breakpoints=list(xs))


def make_exotique(T=1.0):
    """Signal oscillant infiniment vite au voisinage de t=0 (sin(2π/t)).
    Continu partout ailleurs, mais nombre INFINI d'extrema par période :
    ne respecte pas les conditions de Dirichlet, malgré une allure
    globalement discrète sur un tracé classique.
    """
    def f(t):
        tt = t % T
        if tt == 0.0:
            return 0.0
        return np.sin(2 * np.pi / tt)
    return Signal(f=f, T=T, name="exotique", breakpoints=[])


# ─── Coefficients de Fourier ────────────────────────────────────────

def fourier_coeffs(signal, N):
    """Calcule a0, (a_n)ـ1..N, (b_n)_1..N par intégration numérique des
    formules du cours (scipy.integrate.quad), sans passer par un
    échantillonnage régulier ni par une FFT.

    Les intégrales a_n, b_n sont oscillatoires (facteur cos/sin(n w0 t)) :
    on utilise le mode `weight='cos'/'sin'` de quad, prévu pour ce cas et
    fiable même pour n grand (un quad "nu" perd en précision dès que n
    dépasse ~200, faute de subdivisions). Ce mode n'accepte pas de points
    de rupture : quand le signal en a, on découpe l'intégrale à la main
    sur chaque morceau et on additionne.
    """
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

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", IntegrationWarning)
        a0 = quad(f, 0, T, points=pts, limit=200)[0] / T
        an = np.zeros(N)
        bn = np.zeros(N)
        for n in range(1, N + 1):
            an[n - 1] = 2 / T * oscillatory_integral("cos", n * w0)
            bn[n - 1] = 2 / T * oscillatory_integral("sin", n * w0)

    if caught:
        print(f"[{signal.name}] le calcul numérique peine à converger "
              f"({len(caught)} intégrale(s) difficile(s)) -- signe que ce "
              f"signal a un comportement inhabituel quelque part sur sa période.")

    return a0, an, bn


def parseval_energy(a0, an, bn):
    """a0^2 + (1/2) sum_n (a_n^2 + b_n^2) -- calculée à partir des seuls
    coefficients, sans toucher au signal lui-même."""
    return a0**2 + 0.5 * np.sum(an**2 + bn**2)


def signal_energy(signal):
    """(1/T) intégrale sur une période de f(t)^2 -- calculée directement
    à partir du signal, sans passer par ses coefficients de Fourier."""
    T, f = signal.T, signal.f
    pts = [p for p in signal.breakpoints if 0 < p < T] or None
    return quad(lambda t: f(t) ** 2, 0, T, points=pts, limit=200)[0] / T


def reconstruct(t, T, a0, an, bn, orders=None, include_dc=True):
    """Évalue a0 + sum_n [a_n cos(n w0 t) + b_n sin(n w0 t)] sur le
    vecteur t. `orders` : rangs à inclure (par défaut, tous ceux fournis).
    """
    w0 = 2 * np.pi / T
    t = np.asarray(t, dtype=float)
    orders = range(1, len(an) + 1) if orders is None else orders
    s = np.full_like(t, a0 if include_dc else 0.0)
    for n in orders:
        s = s + an[n - 1] * np.cos(n * w0 * t) + bn[n - 1] * np.sin(n * w0 * t)
    return s


# ─── Tracés ──────────────────────────────────────────────────────────

def plot_signal(signal, t, title=None, ax=None, xlim=None, ylim=None):
    """xlim, ylim : ex. (0, 0.1) -- permet de "zoomer" sans tracé interactif :
    changez ces bornes et ré-exécutez la cellule."""
    ax = ax or plt.gca()
    y = [signal.f(ti) for ti in t]
    ax.plot(t, y)
    ax.set_xlabel("t")
    ax.set_title(title or signal.name)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    return ax


def plot_spectrum(an, bn, title=None, ax=None, xlim=None, ylim=None):
    ax = ax or plt.gca()
    n = np.arange(1, len(an) + 1)
    amplitude = np.sqrt(an**2 + bn**2)
    ax.stem(n, amplitude)
    ax.set_xlabel("rang n")
    ax.set_ylabel("|A_n|")
    ax.set_title(title or "Spectre d'amplitude")
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    return ax


def plot_reconstruction(t, reconstructions, original=None, title=None, ax=None,
                         xlim=None, ylim=None):
    """reconstructions : dict {étiquette: valeurs} à superposer.
    xlim, ylim : ex. (0.4, 0.6) -- pour zoomer sur une zone (un saut, un
    dépassement de Gibbs...) : changez ces bornes et ré-exécutez la cellule.
    """
    ax = ax or plt.gca()
    if original is not None:
        ax.plot(t, original, color="black", lw=2, label="signal original", alpha=0.6)
    for label, y in reconstructions.items():
        ax.plot(t, y, label=label)
    ax.set_xlabel("t")
    ax.legend()
    ax.set_title(title or "Reconstruction")
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    return ax
