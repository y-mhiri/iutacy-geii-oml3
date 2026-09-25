# %% [markdown]
# # TP3 --- Transformée de Fourier
#
# **Ce notebook est un compagnon d'exécution, pas le sujet.** Le sujet de
# référence -- questions complètes, contexte, rappels -- est **`TP3.pdf`**.
# Gardez-le ouvert à côté : les cellules ci-dessous sont repérées par leur
# numéro de question dans le PDF (A1, A2, B1, ...) et n'en reprennent que
# l'essentiel.
#
# Vos réponses vont sur **feuille**.

# %%
# Sur Google Colab, seul ce fichier .ipynb est importé (pas le reste du
# dépôt) : on récupère tp3_helpers.py depuis GitHub avant de l'importer.
# En local (Jupyter, VS Code...), ce fichier est déjà à côté du notebook :
# rien ne se passe.
try:
    import google.colab
    IN_COLAB = True
except ImportError:
    IN_COLAB = False

if IN_COLAB:
    import urllib.request
    urllib.request.urlretrieve(
        "https://raw.githubusercontent.com/y-mhiri/iutacy-geii-oml3/main/tp/TP3/notebook/tp3_helpers.py",
        "tp3_helpers.py",
    )

# %%
import numpy as np
import matplotlib.pyplot as plt
from tp3_helpers import (
    make_door, make_beat, explorer,
    fourier_transform, spectrum, inverse_fourier_transform,
    plot_signal, plot_spectrum, plot_reconstruction,
)

# %% [markdown]
# ## A --- Valeur moyenne et exponentielle complexe

# %% [markdown]
# ### A1 : sur papier

# %% [markdown]
# ### A2 : porte, valeur moyenne en cosinus et en sinus

# %%
explorer(f0=1)

# %% [markdown]
# ### A3 : tracer le spectre point par point

# %%
for f0 in [0.1, 1.0, 1.43, 2.0]:  # <- les quatre valeurs demandées
    explorer(f0, montrer_sin=False)

# %% [markdown]
# ## B --- La transformée de Fourier

# %% [markdown]
# ### B1 : sur papier

# %% [markdown]
# ### B2 : spectre de la porte, plusieurs largeurs

# %%
freqs = np.linspace(-8, 8, 800)
for theta in [0.5, 1.0, 2.0]:  # <- modifiez ces largeurs et ré-exécutez
    k = make_door(theta=theta, A=1.0)
    X = spectrum(k, freqs)
    plt.plot(freqs, np.abs(X), label=f"theta={theta}")
plt.legend()
plt.xlabel("f")
plt.ylabel("|X(f)|")
plt.show()

# %% [markdown]
# ### B3 : porte rétrécie, aire conservée

# %%
freqs = np.linspace(-8, 8, 800)
for theta in [1.0, 0.3, 0.05]:  # <- modifiez ces largeurs et ré-exécutez
    p = make_door(theta=theta, A=1.0 / theta)  # aire = 1, cf. impulsion de Dirac
    X = spectrum(p, freqs)
    plt.plot(freqs, X.real, label=f"theta={theta}")
plt.legend()
plt.xlabel("f")
plt.ylabel("X(f)")
plt.show()

# %% [markdown]
# ### B4 : sur papier

# %% [markdown]
# ### B5 : sur papier

# %% [markdown]
# ### B6 : troncature fréquentielle, reconstruction

# %%
theta = 1.0
porte = make_door(theta=theta, A=1.0)


def sinc(x):
    return 1.0 if x == 0 else np.sin(x) / x


def X_theorique(f):
    # X(f) de la porte, formule de B1
    return complex(theta * sinc(np.pi * f * theta), 0.0)


def X_tronque(f, W):
    return X_theorique(f) if abs(f) < W / 2 else 0j


t = np.linspace(-1.5, 1.5, 300)
original = [porte.f(ti) for ti in t]

for W in [4, 10, 40]:  # <- modifiez ces largeurs de bande et ré-exécutez
    y = [inverse_fourier_transform(lambda f: X_tronque(f, W), ti, -W / 2, W / 2)
         for ti in t]
    plot_reconstruction(t, {f"W={W}": y}, original=original,
                         title="reconstruction depuis le spectre tronqué")
    plt.show()

# %% [markdown]
# ## C --- Résolution fréquentielle

# %% [markdown]
# ### C1 : le signal, durée courte

# %%
f1, f2, T = 1.0, 1.2, 3.0  # <- changez T (et ré-exécutez C1, C2)
beat = make_beat(f1, f2, T)
plot_signal(beat, tmin=-T / 2, tmax=T / 2, npts=2000)
plt.show()

# %% [markdown]
# ### C2 : spectre, durée courte

# %%
freqs = np.linspace(0.7, 1.5, 400)
X = spectrum(beat, freqs)
plot_spectrum(freqs, X, title=f"Spectre, T={T}")
plt.show()

# %% [markdown]
# ### C3 : durée longue

# %%
T = 8.0  # <- changez T (et ré-exécutez cette cellule)
beat = make_beat(f1, f2, T)
plot_signal(beat, tmin=-T / 2, tmax=T / 2, npts=2000)
plt.show()

freqs = np.linspace(0.7, 1.5, 400)
X = spectrum(beat, freqs)
plot_spectrum(freqs, X, title=f"Spectre, T={T}")
plt.show()

# %% [markdown]
# ### C4 : comparer plusieurs durées

# %%
freqs = np.linspace(0.7, 1.5, 400)
for T in [3, 5, 8, 15]:  # <- durées d'observation à comparer
    beat = make_beat(f1, f2, T)
    X = spectrum(beat, freqs)
    plot_spectrum(freqs, X, title=f"T = {T}")
    plt.show()

# %% [markdown]
# ### C5 : sur papier
