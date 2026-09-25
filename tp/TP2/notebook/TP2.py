# %% [markdown]
# # TP2 --- Série de Fourier, forme complexe
#
# **Ce notebook est un compagnon d'exécution, pas le sujet.** Le sujet de
# référence -- questions complètes, contexte, rappels -- est **`TP2.pdf`**.
# Gardez-le ouvert à côté : les cellules ci-dessous sont repérées par leur
# numéro de question dans le PDF (A1, A2, B1, ...) et n'en reprennent que
# l'essentiel.
#
# Contrairement au TP1, vos réponses vont sur **feuille** (calculs sur
# papier, réponses aux questions) : ce notebook ne sert qu'à exécuter et
# observer.

# %%
# Sur Google Colab, seul ce fichier .ipynb est importé (pas le reste du
# dépôt) : on récupère tp2_helpers.py depuis GitHub avant de l'importer.
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
        "https://raw.githubusercontent.com/y-mhiri/iutacy-geii-oml3/main/tp/TP2/notebook/tp2_helpers.py",
        "tp2_helpers.py",
    )

# %%
import numpy as np
import matplotlib.pyplot as plt
from tp2_helpers import (
    make_square, make_levels, shift_signal, combine_signals,
    fourier_coeffs, to_complex_coeffs, reconstruct_complex,
    plot_signal, plot_spectrum_complex, plot_reconstruction,
)

# %% [markdown]
# ## A --- Forme complexe
#
# ### A1 : nombres complexes en Python

# %%
z = 2 + 3j
print("z =", z)
print("Re(z) =", z.real, " Im(z) =", z.imag)

x = 2 + 0j  # un réel, écrit comme un complexe de partie imaginaire nulle
print("x =", x, " Im(x) =", x.imag)

# %%
for theta, label in [(0, "0"), (2 * np.pi, "2*pi"), (4 * np.pi, "4*pi"),
                      (np.pi, "pi"), (3 * np.pi, "3*pi"),
                      (np.pi / 2, "pi/2"), (np.pi / 3, "pi/3")]:
    val = np.exp(1j * theta)
    print(f"exp(j*{label:6s}) = {val:.4f}")

# %% [markdown]
# ### A3 : vérification numérique

# %%
carre = make_square(T=1.0, A=1.0)
a0, an, bn = fourier_coeffs(carre, N=10)
c = to_complex_coeffs(a0, an, bn)
for n in range(-5, 6):
    print(f"c_{n:+d} = {c[n]:.4f}")

# %% [markdown]
# ### A4 : spectre bilatéral

# %%
plot_spectrum_complex(c, title="Spectre bilatéral -- carré")
plt.show()

# %%
plot_spectrum_complex(c, freq=True, T=carre.T,
                       title="Spectre bilatéral -- carré (en fréquence)")
plt.show()

# %% [markdown]
# ### A6 : coefficients négatifs retirés

# %%
c_asym = {n: (v if n >= 0 else 0j) for n, v in c.items()}
t = np.linspace(-1, 1, 1000)
y_sym = reconstruct_complex(t, carre.T, c)
y_asym = reconstruct_complex(t, carre.T, c_asym)
original = [carre.f(ti) for ti in t]

plot_reconstruction(t, {"symétrique (tous les n)": y_sym.real,
                         "asymétrique (n >= 0 seulement)": y_asym.real},
                     original=original, title="Reconstruction (partie réelle)")
plt.show()

plot_reconstruction(t, {"symétrique (tous les n)": y_sym.imag,
                         "asymétrique (n >= 0 seulement)": y_asym.imag},
                     title="Reconstruction (partie imaginaire)")
plt.show()

# %% [markdown]
# ## B --- Décalage temporel
#
# ### B1-B2 : observer le décalage

# %%
f = make_square(T=1.0, A=1.0)
phi = 0.15  # <- modifiez cette valeur et ré-exécutez
g = shift_signal(f, phi)
plot_signal(f, tmin=-1, tmax=1, npts=2000, title="f")
plt.show()
plot_signal(g, tmin=-1, tmax=1, npts=2000, title=f"g(t) = f(t + {phi})")
plt.show()

# %% [markdown]
# ### B4 : vérification numérique

# %%
a0f, anf, bnf = fourier_coeffs(f, N=10)
cf = to_complex_coeffs(a0f, anf, bnf)
a0g, ang, bng = fourier_coeffs(g, N=10)
dg = to_complex_coeffs(a0g, ang, bng)

w0 = 2 * np.pi / f.T
dg_predit = {n: cf[n] * np.exp(1j * n * w0 * phi) for n in cf}

for n in range(1, 6):
    print(f"n={n}  d_n (numérique) = {dg[n]:.4f}   "
          f"d_n (prédit depuis c_n) = {dg_predit[n]:.4f}")

# %% [markdown]
# ### B5 : comparaison des spectres

# %%
plot_spectrum_complex(cf, title="Spectre de f")
plt.show()
plot_spectrum_complex(dg, title="Spectre de g = f(t+phi)")
plt.show()

# %% [markdown]
# ## C --- Superposition d'un signal et de son décalage
#
# ### C1 : le signal à N paliers

# %%
N = 4
f = make_levels(N=N, T=1.0, A=1.0, seed=3)
plot_signal(f, tmin=-1, tmax=2, npts=2000)
plt.show()

# %% [markdown]
# ### C2 : somme avec un décalage d'une demi-période

# %%
phi = f.T / 2
g = shift_signal(f, phi)
h = combine_signals(f, g, lambda a, b: a + b, name="f(t)+f(t+T/2)")
plot_signal(h, tmin=-1, tmax=2, npts=2000)
plt.show()

# %% [markdown]
# ### C4 : vérification numérique

# %%
Ncoef = 8
a0f, anf, bnf = fourier_coeffs(f, N=Ncoef)
cf = to_complex_coeffs(a0f, anf, bnf)
a0h, anh, bnh = fourier_coeffs(h, N=Ncoef)
ch = to_complex_coeffs(a0h, anh, bnh)

w0 = 2 * np.pi / f.T
for n in range(1, Ncoef + 1):
    predit = abs(cf[n]) * abs(1 + np.exp(1j * n * w0 * phi))
    print(f"n={n}  |c_n[h]| numérique = {abs(ch[n]):.5f}   "
          f"prédit = {predit:.5f}")

# %% [markdown]
# ### C5 : comparaison des spectres

# %%
plot_spectrum_complex(cf, title="Spectre de f")
plt.show()
plot_spectrum_complex(ch, title="Spectre de h = f(t)+f(t+T/2)")
plt.show()

# %% [markdown]
# ### C6 : recommencer avec un décalage adapté

# %%
phi2 = f.T / 4  # moitié de la période apparente de h (T/2)
h2 = combine_signals(h, shift_signal(h, phi2), lambda a, b: a + b, name="h2")

plot_signal(h2, tmin=-1, tmax=2, npts=2000)
plt.show()

a0h2, anh2, bnh2 = fourier_coeffs(h2, N=Ncoef)
ch2 = to_complex_coeffs(a0h2, anh2, bnh2)
plot_spectrum_complex(ch2, title="Spectre de h2 = h(t)+h(t+T/4)")
plt.show()

print("h2 est-il constant ? valeurs à quelques instants :")
for t in [0.05, 0.2, 0.5, 0.83]:
    print(f"  t={t:.2f}  h2={h2.f(t):.5f}")
