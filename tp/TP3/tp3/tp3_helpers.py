"""TP3 -- Transformée de Fourier. Fonctions fournies, à ne pas modifier.

Signal : un signal défini par sa fonction f (à valeurs scalaires), son
support (tmin, tmax) en dehors duquel f vaut 0, et éventuellement une
liste de segments (sous-intervalles du support) qui accélèrent et
fiabilisent l'intégration numérique quand f est une somme de morceaux
séparés (ex. plusieurs impulsions). Par défaut, segments = [support].
"""

from dataclasses import dataclass, field
from typing import Callable

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad


@dataclass
class Signal:
    f: Callable[[float], float]
    support: tuple
    name: str
    segments: list = field(default_factory=list)


# ─── Signaux ────────────────────────────────────────────────────────

def make_door(theta=1.0, A=1.0):
    """Porte de largeur theta, hauteur A, centrée en 0."""
    def f(t):
        return A if abs(t) < theta / 2 else 0.0
    return Signal(f=f, support=(-theta / 2, theta / 2), name="porte")


def make_beat(f1=1.0, f2=1.2, T=3.0, A=1.0):
    """Somme de deux cosinus de fréquences proches f1, f2, observée sur
    une fenêtre temporelle de durée T (signal nul en dehors)."""
    def f(t):
        return A * (np.cos(2 * np.pi * f1 * t) + np.cos(2 * np.pi * f2 * t)) \
            if abs(t) < T / 2 else 0.0
    return Signal(f=f, support=(-T / 2, T / 2), name=f"battement (T={T})")


# ─── Transformée de Fourier ─────────────────────────────────────────

def fourier_transform(signal, freq, limit=200):
    """Calcule X(freq) = intégrale de x(t) exp(-j 2 pi freq t) sur le
    support de x, par intégration numérique (scipy.integrate.quad, mode
    oscillatoire weight='cos'/'sin', fiable même pour freq grand). Somme
    sur `signal.segments` si fourni (utile quand x est nul entre
    plusieurs morceaux séparés), sinon sur le support entier.
    """
    segs = signal.segments or [signal.support]
    w = 2 * np.pi * freq
    x = signal.f
    re = sum(quad(x, a, b, weight="cos", wvar=w, limit=limit)[0] for a, b in segs)
    im = -sum(quad(x, a, b, weight="sin", wvar=w, limit=limit)[0] for a, b in segs)
    return complex(re, im)


def spectrum(signal, freqs, limit=200):
    """X(f) pour chaque f de `freqs` -- tableau complexe."""
    return np.array([fourier_transform(signal, f, limit=limit) for f in freqs])


def inverse_fourier_transform(X, t, fmin, fmax, limit=200):
    """Calcule x(t) = intégrale de X(f) exp(+j 2 pi f t) sur [fmin, fmax]
    par intégration numérique. X : fonction f -> complexe.
    """
    Xr = lambda f: X(f).real
    Xi = lambda f: X(f).imag
    w = 2 * np.pi * t
    a = quad(Xr, fmin, fmax, weight="cos", wvar=w, limit=limit)[0]
    b = quad(Xi, fmin, fmax, weight="sin", wvar=w, limit=limit)[0]
    c = quad(Xr, fmin, fmax, weight="sin", wvar=w, limit=limit)[0]
    d = quad(Xi, fmin, fmax, weight="cos", wvar=w, limit=limit)[0]
    return complex(a - b, c + d)


def explorer(f0, theta=1.0, A=1.0, montrer_sin=True, tmin=-1.5, tmax=1.5, npts=2000):
    """Pour la porte de largeur theta et hauteur A : trace le produit
    avec cos(2*pi*f0*t) (et avec sin(2*pi*f0*t) si montrer_sin), l'aire
    signée coloriée, et la moyenne -- partie réelle de X(f0) pour le
    produit en cosinus, opposé de la partie imaginaire pour le produit en
    sinus. Affiche aussi la valeur théorique A*theta*sinc(pi*f0*theta),
    à comparer à la valeur numérique.
    """
    porte = make_door(theta=theta, A=A)
    t = np.linspace(tmin, tmax, npts)
    x = np.array([porte.f(ti) for ti in t])
    Xf0 = fourier_transform(porte, f0)
    produits = {"cos": (x * np.cos(2 * np.pi * f0 * t), Xf0.real)}
    if montrer_sin:
        produits["sin"] = (x * np.sin(2 * np.pi * f0 * t), -Xf0.imag)

    fig, axes = plt.subplots(1, len(produits), figsize=(5.5 * len(produits), 4),
                              squeeze=False)
    for ax, nom in zip(axes[0], produits):
        produit, moyenne = produits[nom]
        ax.plot(t, x, color="gray", alpha=0.3)
        ax.plot(t, produit, color="C0")
        ax.fill_between(t, produit, 0, where=(produit >= 0), color="C0", alpha=0.3)
        ax.fill_between(t, produit, 0, where=(produit < 0), color="C3", alpha=0.3)
        ax.axhline(moyenne, color="k", ls="--", label=f"moyenne = {moyenne:.4f}")
        ax.set_title(f"porte(t) x {nom}(2 pi f0 t)")
        ax.legend()
    plt.suptitle(f"f0 = {f0}")
    plt.show()

    arg = np.pi * f0 * theta
    theorique = A * theta * (1.0 if arg == 0 else np.sin(arg) / arg)
    print(f"valeur théorique (formule de A1) : {theorique:.4f}")
    print(f"valeur numérique (moyenne en cosinus) : {produits['cos'][1]:.4f}")


# ─── Tracés ──────────────────────────────────────────────────────────

def plot_signal(signal, t=None, tmin=None, tmax=None, npts=2000,
                 title=None, ax=None, xlim=None, ylim=None):
    """Trace signal.f sur une grille temporelle (voir TP1/TP2 : passer
    `t`, ou `tmin`/`tmax`/`npts`)."""
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


def plot_spectrum(freqs, values, title=None, ax=None, xlim=None, ylim=None,
                   kind="abs"):
    """Trace le spectre en fonction de la fréquence. `values` : tableau
    complexe (ex. sortie de `spectrum`). `kind` : "abs" (module, par
    défaut), "re" ou "im"."""
    ax = ax or plt.gca()
    values = np.asarray(values)
    y = {"abs": np.abs, "re": np.real, "im": np.imag}[kind](values)
    ax.plot(freqs, y)
    ax.set_xlabel("fréquence f")
    ax.set_ylabel({"abs": "|X(f)|", "re": "Re[X(f)]", "im": "Im[X(f)]"}[kind])
    ax.set_title(title or "Spectre")
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    return ax


def plot_reconstruction(t, reconstructions, original=None, title=None,
                         ax=None, xlim=None, ylim=None):
    """reconstructions : dict {étiquette: valeurs} (réelles ou complexes,
    partie réelle tracée) à superposer à `original`."""
    ax = ax or plt.gca()
    if original is not None:
        ax.plot(t, original, color="black", lw=2, label="signal original",
                 alpha=0.6)
    for label, y in reconstructions.items():
        y = np.asarray(y)
        if np.iscomplexobj(y):
            y = y.real
        ax.plot(t, y, label=label)
    ax.set_xlabel("t")
    ax.legend()
    ax.set_title(title or "Reconstruction")
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    return ax
