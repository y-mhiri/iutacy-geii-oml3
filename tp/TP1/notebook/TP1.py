# %% [markdown]
# # TP1 --- Série de Fourier
#
# **Ce notebook est un compagnon d'exécution, pas le sujet.** Le sujet de
# référence -- questions complètes, contexte, rappels -- est **`TP1.pdf`**.
# Gardez-le ouvert à côté : les cellules ci-dessous sont repérées par leur
# numéro de question dans le PDF (A1, A2, B1, ...) et n'en reprennent que
# l'essentiel.
#
# Toutes les fonctions dont vous avez besoin sont déjà écrites dans
# `tp1_helpers.py`. Vous n'avez pas à les modifier : votre travail est de
# faire tourner les cellules, changer les valeurs indiquées, lire les
# tracés, et répondre aux questions du PDF.

# %%
import numpy as np
import matplotlib.pyplot as plt
from tp1_helpers import (
    make_square, make_random_walk, make_exotique, make_sawtooth,
    fourier_coeffs, reconstruct, parseval_energy, signal_energy,
    plot_signal, plot_spectrum, plot_reconstruction,
)

# %% [markdown]
# ## A --- Première manipulation de séries trigonométriques
#
# ### A1-A2 : sinus et cosinus

# %%
f = 1.0
t = np.linspace(-10, 10, 2000)
plt.plot(t, np.sin(2 * np.pi * f * t), label="sin")
plt.plot(t, np.cos(2 * np.pi * f * t), label="cos")
plt.xlabel("t")
plt.legend()
plt.show()

# %% [markdown]
# _Réponse A2 :_

# %% [markdown]
# ### A3 : même chose à 440 Hz, puis zoom
#
# On trace sur une fenêtre large, puis on zoome avec `xlim` -- pas besoin de
# refaire le tracé, changez juste les bornes et ré-exécutez.

# %%
f = 440.0
t = np.linspace(-0.02, 0.02, 4000)
plt.plot(t, np.sin(2 * np.pi * f * t), label="sin")
plt.plot(t, np.cos(2 * np.pi * f * t), label="cos")
plt.xlabel("t")
plt.legend()
plt.xlim(-0.002, 0.002)  # <- zoom sur ~une période : modifiez ces bornes et ré-exécutez
plt.show()

# %% [markdown]
# _Réponse A3 :_

# %% [markdown]
# ### A4 : une série qui ne vient d'aucun signal
#
# $T=1$, $a_n=b_n=1$ pour tout $n$. Choisissez une valeur de $N$, regardez,
# recommencez avec une autre valeur.

# %%
T = 1.0
t = np.linspace(0, T, 2000)
N = 5  # <- modifiez cette valeur (essayez 1, 2, 5, 20, 100) et ré-exécutez
an = np.ones(N)
bn = np.ones(N)
y = reconstruct(t, T, a0=0.0, an=an, bn=bn)
plot_reconstruction(t, {f"N={N}": y}, title="a_n = b_n = 1")
plt.show()

# %% [markdown]
# Cas extrême, pour finir :

# %%
N = 100_000
an = np.ones(N)
bn = np.ones(N)
y = reconstruct(t, T, a0=0.0, an=an, bn=bn)  # quelques secondes
plt.plot(t, y)
plt.title(f"N={N}")
plt.show()

# %% [markdown]
# _Réponse A4 :_

# %% [markdown]
# ## B --- Analyse et synthèse d'un signal périodique par série de Fourier
#
# ### B1 : le signal carré

# %%
carre = make_square(T=1.0, A=1.0)
t = np.linspace(-10, 10, 4000)
plot_signal(carre, t, xlim=(-2, 2))  # <- modifiez xlim pour zoomer/dézoomer
plt.show()

# %% [markdown]
# ### B2 : coefficients de Fourier

# %%
a0, an, bn = fourier_coeffs(carre, N=100)
print("a0 =", a0)
print("b1..b5 =", bn[:5])

# %% [markdown]
# _Réponse B2 :_

# %% [markdown]
# ### B3 : spectre

# %%
plot_spectrum(an, bn)
plt.show()

# %% [markdown]
# ### B4 : reconstruction
#
# Choisissez une valeur de $N$, regardez, recommencez avec une autre valeur.

# %%
t = np.linspace(-2, 2, 2000)
original = [carre.f(ti) for ti in t]
N = 10  # <- modifiez cette valeur (essayez 2, 5, 10, 100) et ré-exécutez
y = reconstruct(t, carre.T, a0, an, bn, orders=range(1, N + 1))
plot_reconstruction(t, {f"N={N}": y}, original=original,
                     xlim=(-0.2, 0.2))  # <- zoomez près d'un saut pour observer Gibbs
plt.show()

# %% [markdown]
# _Réponse B4 :_

# %% [markdown]
# ### B5 : reconstruction avec des rangs élevés seulement

# %%
n_min, n_max = 50, 100  # <- modifiez ces bornes et ré-exécutez
y_partiel = reconstruct(t, carre.T, a0, an, bn,
                         orders=range(n_min, n_max + 1), include_dc=False)
plot_reconstruction(t, {f"rangs {n_min}-{n_max}": y_partiel}, original=original)
plt.show()

# %% [markdown]
# _Réponse B5 :_

# %% [markdown]
# ## C --- Deux signaux à l'allure trompeuse
#
# On traite les deux signaux l'un après l'autre, en entier, plutôt qu'en
# parallèle.

# %% [markdown]
# ### C1-C2 : la marche aléatoire

# %%
marche = make_random_walk(T=1.0, n_points=25, seed=42)
t = np.linspace(0, 3, 3000)
plot_signal(marche, t)  # <- ajoutez xlim=(...) pour zoomer sur un segment
plt.show()

# %% [markdown]
# _Réponse C1 :_

# %%
t = np.linspace(0, 1, 2000)
original = [marche.f(ti) for ti in t]
a0, an, bn = fourier_coeffs(marche, N=60)
N = 20  # <- modifiez cette valeur (essayez 5, 20, 60) et ré-exécutez
y = reconstruct(t, marche.T, a0, an, bn, orders=range(1, N + 1))
plot_reconstruction(t, {f"N={N}": y}, original=original, title=marche.name)
plt.show()

# %% [markdown]
# _Réponse C2 :_

# %% [markdown]
# ### C3-C4 : le signal exotique

# %%
exotique = make_exotique(T=1.0)
t = np.linspace(0, 3, 3000)
plot_signal(exotique, t)  # <- essayez xlim=(0, 0.05) pour zoomer près de t=0
plt.show()

# %% [markdown]
# _Réponse C3 :_

# %%
t = np.linspace(0, 1, 2000)
original = [exotique.f(ti) for ti in t]
a0, an, bn = fourier_coeffs(exotique, N=60)
N = 20  # <- modifiez cette valeur (essayez 5, 20, 60) et ré-exécutez
y = reconstruct(t, exotique.T, a0, an, bn, orders=range(1, N + 1))
plot_reconstruction(t, {f"N={N}": y}, original=original, title=exotique.name)
plt.show()

# %% [markdown]
# _Réponse C4 :_

# %% [markdown]
# ## D (optionnel, si le temps le permet) --- Le signal dents de scie
#
# ### D1 : afficher le signal

# %%
scie = make_sawtooth(T=1.0)
t = np.linspace(-3, 3, 4000)
plot_signal(scie, t)
plt.show()

# %% [markdown]
# _Réponse D1 :_

# %% [markdown]
# ### D2 : coefficients, reconstruction, zoom sur la jonction

# %%
a0, an, bn = fourier_coeffs(scie, N=100)
t = np.linspace(-0.5, 1.5, 2000)
original = [scie.f(ti) for ti in t]
N = 10  # <- modifiez cette valeur et ré-exécutez
y = reconstruct(t, scie.T, a0, an, bn, orders=range(1, N + 1))
plot_reconstruction(t, {f"N={N}": y}, original=original,
                     xlim=(-0.2, 0.2))  # <- zoomez sur la jonction t=0
plt.show()

# %% [markdown]
# _Réponse D2 :_

# %% [markdown]
# ### D3 : spectre, comparaison avec le signal carré

# %%
plot_spectrum(an, bn, title="Spectre -- dents de scie")
plt.show()

# %% [markdown]
# _Réponse D3 :_ (comparer à B3)

# %% [markdown]
# ### D4 : somme du spectre vs. énergie du signal

# %%
print("énergie du signal (intégrale) :", signal_energy(scie))
for N in [1, 5, 20, 100]:
    a0, an, bn = fourier_coeffs(scie, N)
    print(f"N={N:4d}  somme du spectre (Parseval) = {parseval_energy(a0, an, bn):.6f}")

# %% [markdown]
# _Réponse D4 :_
