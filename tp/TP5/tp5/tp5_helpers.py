"""TP5 -- Transformée de Fourier Discrète. Fonctions fournies, à ne pas modifier.

Un signal discret est ici simplement un tableau numpy x de longueur N,
x[n] pour n=0,...,N-1 -- pas de classe dédiée, contrairement aux TP
précédents : il n'y a plus d'intégrale à calculer, seulement des sommes
finies.
"""

import numpy as np
import matplotlib.pyplot as plt


def tfd(x):
    """Calcule la TFD de x : X[k] = somme_{n=0}^{N-1} x[n] exp(-2j pi k n / N),
    pour k=0,...,N-1. Implémentation directe (pas de FFT), pour que le
    calcul corresponde exactement à la formule du cours."""
    x = np.asarray(x, dtype=complex)
    N = len(x)
    k = np.arange(N).reshape(-1, 1)
    n = np.arange(N).reshape(1, -1)
    W = np.exp(-2j * np.pi * k * n / N)
    return W @ x


def tfd_k(x, k):
    """Calcule X[k] = somme_{n=0}^{N-1} x[n] exp(-2j pi k n / N) pour un
    tableau k quelconque -- la même formule que tfd, mais évaluée pour
    les rangs k fournis, y compris en dehors de 0,...,N-1."""
    x = np.asarray(x, dtype=complex)
    N = len(x)
    k = np.asarray(k).reshape(-1, 1)
    n = np.arange(N).reshape(1, -1)
    W = np.exp(-2j * np.pi * k * n / N)
    return (W @ x).flatten()


def plot_tfd_k(k, X, kind="abs", ax=None, title=None):
    """Comme plot_tfd, mais pour des rangs k quelconques (fournis
    explicitement, pas forcément 0,...,N-1) -- à utiliser avec tfd_k."""
    ax = ax or plt.gca()
    y = {"abs": np.abs, "re": np.real, "im": np.imag}[kind](X)
    markerline, stemlines, baseline = ax.stem(k, y)
    plt.setp(baseline, visible=False)
    ax.set_xlabel("rang k")
    ax.set_ylabel({"abs": "|X[k]|", "re": "Re[X[k]]", "im": "Im[X[k]]"}[kind])
    ax.set_title(title or "Spectre (TFD)")
    return ax


def plot_signal_discret(x, ax=None, title=None):
    """Trace x[n] en bâtons, n=0,...,N-1."""
    ax = ax or plt.gca()
    N = len(x)
    markerline, stemlines, baseline = ax.stem(np.arange(N), np.real(x))
    plt.setp(baseline, visible=False)
    ax.set_xlabel("n")
    ax.set_ylabel("x[n]")
    ax.set_title(title or "Signal discret")
    return ax


def plot_tfd(X, kind="abs", ax=None, title=None):
    """Trace le spectre en bâtons, k=0,...,N-1. `kind` : "abs" (module,
    par défaut), "re" ou "im"."""
    ax = ax or plt.gca()
    N = len(X)
    y = {"abs": np.abs, "re": np.real, "im": np.imag}[kind](X)
    markerline, stemlines, baseline = ax.stem(np.arange(N), y)
    plt.setp(baseline, visible=False)
    ax.set_xlabel("rang k")
    ax.set_ylabel({"abs": "|X[k]|", "re": "Re[X[k]]", "im": "Im[X[k]]"}[kind])
    ax.set_title(title or "Spectre (TFD)")
    return ax


def make_square_discret(nu0, N):
    """Signal carré discret x[n] = +1 si n*nu0 mod 1 < 1/2, -1 sinon,
    n=0,...,N-1 -- un créneau de rapport cyclique 50%, fréquence
    normalisée nu0 (pas besoin que nu0 soit un multiple de 1/N).
    Construit à partir de la phase n*nu0 mod 1 (et non de sign(sin(...)))
    pour éviter tout problème d'arrondi exactement aux instants de
    transition."""
    n = np.arange(N)
    phase = (nu0 * n) % 1.0
    return np.where(phase < 0.5, 1.0, -1.0)


def zero_pad(x, N):
    """Complète x par des zéros, avant et après, jusqu'à la longueur N
    (le signal d'origine reste centré dans la fenêtre agrandie)."""
    pad_total = N - len(x)
    before = pad_total // 2
    after = pad_total - before
    return np.concatenate([np.zeros(before), x, np.zeros(after)])
