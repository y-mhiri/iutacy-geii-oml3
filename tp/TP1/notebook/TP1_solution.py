# %% [markdown]
# # TP1 --- Série de Fourier (corrigé)
#
# Compagnon de `TP1.pdf` (sujet de référence). Cellules repérées par
# numéro de question (A1, A2, B1, ...).

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
# _Réponse A1 :_ $T=1/f=1\,$s, donc 20 périodes sur $[-10,10]$.
#
# _Réponse A2 :_ $\cos(t)=\sin(t+\pi/2)$ : cos est en avance d'un quart de
# période sur sin. Sur le graphe, cos vaut son maximum en $t=0$ tandis que
# sin y est nul (et croissant) -- c'est ce qui permet de les distinguer.

# %% [markdown]
# ### A3 : même chose à 440 Hz, puis zoom

# %%
f = 440.0
t = np.linspace(-0.02, 0.02, 4000)
plt.plot(t, np.sin(2 * np.pi * f * t), label="sin")
plt.plot(t, np.cos(2 * np.pi * f * t), label="cos")
plt.xlabel("t")
plt.legend()
plt.xlim(-0.002, 0.002)
plt.show()

# %% [markdown]
# _Réponse A3 :_ $T=1/440\approx 2{,}27\,$ms. Le décalage reste un quart de
# période (proportion inchangée), donc environ $0{,}57\,$ms en valeur
# absolue : seule sa durée en secondes diminue avec $f$.

# %% [markdown]
# ### A4 : une série qui ne vient d'aucun signal

# %%
T = 1.0
t = np.linspace(0, T, 2000)
N = 5
an = np.ones(N)
bn = np.ones(N)
y = reconstruct(t, T, a0=0.0, an=an, bn=bn)
plot_reconstruction(t, {f"N={N}": y}, title="a_n = b_n = 1")
plt.show()

# %%
N = 100_000
an = np.ones(N)
bn = np.ones(N)
y = reconstruct(t, T, a0=0.0, an=an, bn=bn)
plt.plot(t, y)
plt.title(f"N={N}")
plt.show()

# %% [markdown]
# _Réponse A4 :_ La période ne change pas (chaque terme est de période
# $\le 1$, donc la somme reste $1$-périodique quel que soit $N$). En
# revanche la série ne converge pas : en $t=0$ chaque terme vaut 1, donc la
# somme partielle vaut exactement $N$ -- elle diverge (visible
# numériquement : la valeur en $t=0$ atteint $10^5$ pour $N=10^5$).
# Ailleurs, la somme oscille sans se stabiliser. C'est attendu : rien ne
# garantit qu'une série trigonométrique quelconque converge -- c'est le
# rôle des conditions de Dirichlet (section suivante) que de le garantir.

# %% [markdown]
# ## B --- Analyse et synthèse d'un signal périodique par série de Fourier
#
# ### B1 : le signal carré

# %%
carre = make_square(T=1.0, A=1.0)
t = np.linspace(-10, 10, 4000)
plot_signal(carre, t, xlim=(-2, 2))
plt.show()

# %% [markdown]
# ### B2 : coefficients de Fourier

# %%
a0, an, bn = fourier_coeffs(carre, N=100)
print("a0 =", a0)
print("b1..b5 =", bn[:5])
print("théorie b_n = 4/(n*pi), n impair :", [4 / (n * np.pi) for n in range(1, 6, 2)])

# %% [markdown]
# _Réponse B2 :_ D'après le cours (TD1), $a_n=0$ pour tout $n$, et
# $b_n=4/(n\pi)$ si $n$ impair, $0$ si $n$ pair. Les valeurs numériques
# concordent (aux erreurs d'arrondi près sur $a_n$, $\sim 10^{-17}$).

# %% [markdown]
# ### B3 : spectre

# %%
plot_spectrum(an, bn)
plt.show()

# %% [markdown]
# ### B4 : reconstruction

# %%
t = np.linspace(-2, 2, 2000)
original = [carre.f(ti) for ti in t]
N = 10
y = reconstruct(t, carre.T, a0, an, bn, orders=range(1, N + 1))
plot_reconstruction(t, {f"N={N}": y}, original=original, xlim=(-0.2, 0.2))
plt.show()

# %% [markdown]
# _Réponse B4 :_ Plus $N$ augmente, plus la reconstruction se rapproche du
# signal carré loin des sauts. Mais au voisinage de chaque discontinuité,
# un dépassement d'environ 9% de l'amplitude du saut persiste quel que soit
# $N$ (bien visible en zoomant avec `xlim` près d'un saut) : il se resserre
# sans jamais disparaître. C'est le phénomène de Gibbs.

# %% [markdown]
# ### B5 : reconstruction avec des rangs élevés seulement

# %%
n_min, n_max = 50, 100
y_partiel = reconstruct(t, carre.T, a0, an, bn,
                         orders=range(n_min, n_max + 1), include_dc=False)
plot_reconstruction(t, {f"rangs {n_min}-{n_max}": y_partiel}, original=original)
plt.show()

# %% [markdown]
# _Réponse B5 :_ Le résultat n'a plus l'allure d'un carré : c'est une
# oscillation rapide de faible amplitude (les $b_n$ de rang 50-100 valent
# $4/(n\pi)\approx 0{,}01$--$0{,}025$, bien plus petits que l'amplitude 1
# d'origine). Les harmoniques basses portent la forme globale du signal ;
# les hautes n'ajoutent que du détail fin, avec une énergie bien plus
# faible.

# %% [markdown]
# ## C --- Deux signaux à l'allure trompeuse

# %% [markdown]
# ### C1-C2 : la marche aléatoire

# %%
marche = make_random_walk(T=1.0, n_points=25, seed=42)
t = np.linspace(0, 3, 3000)
plot_signal(marche, t)
plt.show()

# %% [markdown]
# _Réponse C1 :_ Oui, périodique par construction. Il respecte Dirichlet
# malgré son allure erratique : construit à partir d'un nombre **fixe** de
# sommets (25) reliés par des segments de droite -- donc un nombre fini
# d'extrema, aucune discontinuité (la dernière valeur est recalée sur la
# première), et borné.

# %%
t = np.linspace(0, 1, 2000)
original = [marche.f(ti) for ti in t]
a0, an, bn = fourier_coeffs(marche, N=60)
N = 20
y = reconstruct(t, marche.T, a0, an, bn, orders=range(1, N + 1))
plot_reconstruction(t, {f"N={N}": y}, original=original, title=marche.name)
plt.show()

# %% [markdown]
# _Réponse C2 :_ La reconstruction converge proprement vers le signal
# quand $N$ augmente (cohérent avec Dirichlet respecté, réponse C1).

# %% [markdown]
# ### C3-C4 : le signal exotique

# %%
exotique = make_exotique(T=1.0)
t = np.linspace(0, 3, 3000)
plot_signal(exotique, t, xlim=(0, 0.05))
plt.show()

# %% [markdown]
# _Réponse C3 :_ Périodique, mais il oscille infiniment vite au voisinage
# de $t=0$ dans chaque période ($\sin(2\pi/t)$ a une infinité d'extrema qui
# s'accumulent quand $t\to0^+$, bien visible en zoomant avec `xlim`). Il ne
# respecte donc **pas** les conditions de Dirichlet (nombre fini
# d'extrema), bien que son tracé à l'échelle globale ait l'air raisonnable.

# %%
t = np.linspace(0, 1, 2000)
original = [exotique.f(ti) for ti in t]
a0, an, bn = fourier_coeffs(exotique, N=60)
N = 20
y = reconstruct(t, exotique.T, a0, an, bn, orders=range(1, N + 1))
plot_reconstruction(t, {f"N={N}": y}, original=original, title=exotique.name)
plt.show()

# %% [markdown]
# _Réponse C4 :_ Le calcul des coefficients affiche lui-même un
# avertissement (« le calcul numérique peine à converger »), et la
# reconstruction reste mauvaise près de $t=0$ même à $N=60$ : cohérent avec
# la non-convergence attendue pour un signal qui ne vérifie pas Dirichlet
# (réponse C3).

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
# _Réponse D1 :_ Une seule discontinuité par période, exactement à la
# jonction $t=0/T$ (contrairement au carré, qui en a une au milieu de la
# période *et* une à la jonction). Le signal monte linéairement puis
# retombe brutalement à chaque période.

# %% [markdown]
# ### D2 : coefficients, reconstruction, zoom sur la jonction

# %%
a0, an, bn = fourier_coeffs(scie, N=100)
t = np.linspace(-0.5, 1.5, 2000)
original = [scie.f(ti) for ti in t]
N = 10
y = reconstruct(t, scie.T, a0, an, bn, orders=range(1, N + 1))
plot_reconstruction(t, {f"N={N}": y}, original=original, xlim=(-0.2, 0.2))
plt.show()

# %% [markdown]
# _Réponse D2 :_ Même phénomène que pour le carré : un dépassement
# d'environ 9% de l'amplitude du saut apparaît juste avant/après la
# jonction et ne disparaît pas quand $N$ augmente -- Gibbs, ici localisé au
# seul endroit où le signal est discontinu.

# %% [markdown]
# ### D3 : spectre, comparaison avec le signal carré

# %%
plot_spectrum(an, bn, title="Spectre -- dents de scie")
plt.show()

# %% [markdown]
# _Réponse D3 :_ Ici $b_n=-1/(n\pi)$ pour **tout** $n$ (spectre dense),
# contre $b_n=4/(n\pi)$ sur les rangs **impairs seulement** pour le carré
# (spectre creux, amplitude plus grande par raie). Mais dans les deux cas
# l'enveloppe décroît en $1/n$ : c'est attendu, puisque les deux signaux
# ont une discontinuité de première espèce -- c'est précisément cette
# décroissance lente qui est responsable du phénomène de Gibbs dans les
# deux cas.

# %% [markdown]
# ### D4 : somme du spectre vs. énergie du signal

# %%
print("énergie du signal (intégrale) :", signal_energy(scie))
for N in [1, 5, 20, 100]:
    a0, an, bn = fourier_coeffs(scie, N)
    print(f"N={N:4d}  somme du spectre (Parseval) = {parseval_energy(a0, an, bn):.6f}")

# %% [markdown]
# _Réponse D4 :_ La somme partielle $a_0^2+\frac12\sum(a_n^2+b_n^2)$ se
# rapproche, quand $N$ augmente, de la même valeur que l'intégrale directe
# $\frac1T\int_0^T f(t)^2\,\d t$ (ici $\approx 0{,}333$ dans les deux cas).
# Autrement dit, l'"énergie" du signal peut se lire soit directement sur
# $f$, soit sur son spectre. On retrouvera ce résultat plus tard dans le
# cours sous le nom de théorème de Parseval.
