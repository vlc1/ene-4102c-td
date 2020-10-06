### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 0162aa9a-f363-11ea-192a-8350869f18df
md"""
Versions [Pluto](https://github.com/vlc1/ene-4102c-td/blob/master/td1/part2.jl) et [Jupyter](https://vlc1.github.io/ene-4102c/td1/part2.ipynb) de ce notebook, et correction au format [Pluto](https://github.com/vlc1/ene-4102c-td/blob/solution/td1/part2.jl).

"""

# ╔═╡ d3f9af78-f35e-11ea-0812-1b1069364fba
md"""
# Exercice 1

Établir pour chacune des fonctions proposées ci-dessous un développement limité en ``0`` à l'ordre ``n``.

| $f$ | $n$ |
|:--------:|:-----:|
| ``x \mapsto \frac{\ln \left ( 1 + x \right )}{1 + x}`` | ``3`` |
| ``x \mapsto \frac{\ln \left [ \cosh \left ( x \right ) \right ]}{x \ln \left ( 1 + x \right )}`` | ``2`` |
| ``x \mapsto \frac{\cos \left ( x \right ) - 1}{\ln \left ( 1 + x \right ) \sinh \left ( x \right )}`` | ``3`` |

"""

# ╔═╡ 07859f52-f360-11ea-0540-893a89f816ff
md"""
# Exercice 2

Déterminer la limite
```math
\lim_{x \to 0} \frac{\sinh \left ( x \right )}{\sin \left ( x \right )}.
```

"""

# ╔═╡ 2937f110-f366-11ea-10b1-1d324a07f056
md"""
# Exercice 3

1. Adapter l'exemple du développement limité de la fonction exponentielle aux fonctions ``\cosh`` et ``\sinh``.
1. Utiliser la bibliothèque `Plots.jl` pour visualiser ces développements limités.

"""

# ╔═╡ Cell order:
# ╟─0162aa9a-f363-11ea-192a-8350869f18df
# ╟─d3f9af78-f35e-11ea-0812-1b1069364fba
# ╟─07859f52-f360-11ea-0540-893a89f816ff
# ╟─2937f110-f366-11ea-10b1-1d324a07f056
