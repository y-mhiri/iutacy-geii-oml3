"""TP4 -- Convolution (domaine continu). Fonctions fournies, à ne pas modifier.

Signal : un signal défini par sa fonction f (à valeurs scalaires), son
support (tmin, tmax) en dehors duquel f vaut 0. Même convention que
TP3.
"""

from dataclasses import dataclass
from typing import Callable

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from math import sin, pi


@dataclass
class Signal:
    f: Callable[[float], float]
    support: tuple
    name: str


# ─── Signaux ────────────────────────────────────────────────────────

def make_door(theta=1.0, A=1.0, t0=0.0):
    """Porte de largeur theta, hauteur A, centrée en t0."""
    def f(t):
        return A if abs(t - t0) < theta / 2 else 0.0
    return Signal(f=f, support=(t0 - theta / 2, t0 + theta / 2), name="porte")


def make_triangle(theta=1.0, A=1.0, t0=0.0):
    """Triangle de demi-largeur theta (support de largeur 2*theta),
    hauteur A, centré en t0."""
    def f(t):
        u = t - t0
        return A * (1 - abs(u) / theta) if abs(u) < theta else 0.0
    return Signal(f=f, support=(t0 - theta, t0 + theta), name="triangle")

def make_sinc(theta=1.0, t0=0.0):
    """Sinus cardinal de demi-largeur theta (support de largeur 2*theta),
    hauteur theta, centré en t0."""
    def f(t):
        u = t - t0
        return theta * sin(pi * theta * u ) / (pi * theta * u)  if abs(u) > 1e-16 else theta
    return Signal(f=f, support=(-100,100), name="sinc")

def make_sinc_abs(theta=1.0, t0=0.0):
    """ Carré d'un sinus cardinal de demi-largeur theta (support de largeur 2*theta),
    hauteur theta, centré en t0."""
    def f(t):
        u = t - t0
        return abs(theta * sin(pi * theta * u) / (pi * theta * u) ) if abs(u) > 1e-16 else theta
    return Signal(f=f, support=(-100,100), name="sinc abs")


def reverse(signal):
    """Retourné temporel : renvoie le Signal x(-t)."""
    a, b = signal.support
    return Signal(f=lambda t: signal.f(-t), support=(-b, -a),
                  name=f"{signal.name} retourné")


# ─── Convolution ──────────────────────────────────────────────────

def convolve(signal, t, limit=200):
    """Calcule (x*y)(t) = intégrale de x(tau) y(t-tau) dtau, pour
    `signal` = (x, y) un couple de Signal, sur le support de x."""
    x, y = signal
    a, b = x.support
    return quad(lambda tau: x.f(tau) * y.f(t - tau), a, b, limit=limit)[0]


def shifted_copies(signal, positions, weights=None):
    """Fonction f telle que f(t) = somme des w_k * signal.f(t - t_k),
    une copie décalée de `signal` par position de `positions` -- c'est
    exactement x convolué avec une somme pondérée d'impulsions de
    Dirac aux instants `positions`, poids `weights` (1 par défaut)."""
    if weights is None:
        weights = [1.0] * len(positions)

    def f(t):
        return sum(w * signal.f(t - t0) for t0, w in zip(positions, weights))
    return f


def moving_average(x_func, t, theta, limit=200):
    """Calcule (x*porte_theta)(t)/theta = moyenne de x sur
    [t-theta/2, t+theta/2]. `x_func` : fonction t -> valeur (pas
    besoin d'un Signal borné, x_func peut être défini sur tout R)."""
    val = quad(x_func, t - theta / 2, t + theta / 2, limit=limit)[0]
    return val / theta


def autocorrelation(signal, tau, tmin=None, tmax=None, limit=200):
    """Calcule R_x(tau) = intégrale de x(t) x(t+tau) dt, sur
    [tmin, tmax] (par défaut, un peu plus large que le support de
    `signal`)."""
    a, b = signal.support
    if tmin is None:
        tmin = a - abs(tau)
    if tmax is None:
        tmax = b + abs(tau)
    return quad(lambda t: signal.f(t) * signal.f(t + tau), tmin, tmax, limit=limit)[0]


# ─── Tracés ──────────────────────────────────────────────────────────

def plot_function(f, tmin, tmax, npts=2000, ax=None, title=None, xlabel="t"):
    """Trace une fonction t -> f(t) quelconque sur [tmin, tmax]."""
    ax = ax or plt.gca()
    t = np.linspace(tmin, tmax, npts)
    y = [f(ti) for ti in t]
    ax.plot(t, y)
    ax.set_xlabel(xlabel)
    ax.set_title(title or "")
    return ax


def plot_signal(signal, tmin=None, tmax=None, npts=2000, ax=None, title=None):
    """Trace signal.f sur [tmin, tmax] (par défaut, le support de signal)."""
    if tmin is None or tmax is None:
        a, b = signal.support
        tmin = tmin if tmin is not None else a
        tmax = tmax if tmax is not None else b
    return plot_function(signal.f, tmin, tmax, npts=npts, ax=ax,
                          title=title or signal.name)
