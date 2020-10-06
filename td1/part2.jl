### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ ecbbc1ca-07a9-11eb-012c-edb5b379dda1
using Plots

# ╔═╡ 0162aa9a-f363-11ea-192a-8350869f18df
md"""
Versions [Pluto](https://github.com/vlc1/ene-4102c-td/blob/master/td1/part2.jl) et [Jupyter](https://vlc1.github.io/ene-4102c/td1/part2.ipynb) de ce notebook.

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

# ╔═╡ 0ecdebd0-07a8-11eb-2504-c9c77f451a58
md"""
## Correction

```math
\frac{\ln \left ( 1 + x \right )}{1 + x} = x - \frac{3 x ^ 2}{2} + \frac{11 x ^ 3}{6} + \mathcal{O} \left ( x ^ 4 \right )
```
```math
\frac{\ln \left [ \cosh \left ( x \right ) \right ]}{x \ln \left ( 1 + x \right )} = \frac{1}{2} + \frac{x}{4} - \frac{x ^ 2}{8} + \mathcal{O} \left ( x ^ 3 \right )
```
```math
\frac{\cos \left ( x \right ) - 1}{\ln \left ( 1 + x \right ) \sinh \left ( x \right )} = -\frac{1}{2} - \frac{x}{4} + \frac{x ^ 2}{6} + \frac{x ^ 3}{24} + \mathcal{O} \left ( x ^ 4 \right )
```

"""

# ╔═╡ 07859f52-f360-11ea-0540-893a89f816ff
md"""
# Exercice 2

Déterminer la limite
```math
\lim_{x \to 0} \frac{\sinh \left ( x \right )}{\sin \left ( x \right )}.
```

"""

# ╔═╡ 7250c606-07a9-11eb-1bcd-a7e371200e58
md"""
## Correction

```math
\lim_{x \to 0} \frac{\sinh \left ( x \right )}{\sin \left ( x \right )} = 1
```

"""

# ╔═╡ 2937f110-f366-11ea-10b1-1d324a07f056
md"""
# Exercice 3

1. Adapter l'exemple du développement limité de la fonction exponentielle aux fonctions ``\cosh`` et ``\sinh``.
1. Utiliser la bibliothèque `Plots.jl` pour visualiser ces développements limités.

"""

# ╔═╡ 9a139d6c-07a9-11eb-2f62-ed69c34cac69
md"""
## Correction

"""

# ╔═╡ a088258c-07a9-11eb-0859-f135b5379d27
function mycos(x, n = 0)
	y = zero(x)
	for i in 0:n
		y += (-1) ^ i * x ^ 2i / factorial(2i)
	end
	y
end

# ╔═╡ f0aaa9e0-07a9-11eb-0b99-f355ceeed754
begin
	local fig
	fig = plot(cos, lw = 2, xlim = (-π, π), ylim = (-1.5, 1.5))
	for n in 0:4
		plot!(fig, x -> mycos(x, n), lw = 2, xlim = (-π, π), label = "DL(0, $(2n))")
	end
	fig
end

# ╔═╡ c79e1fdc-07a9-11eb-34ff-d3cdf9916161
function mysin(x, n = 0)
	y = zero(x)
	for i in 0:n
		y += (-1) ^ i * x ^ (2i + 1) / factorial(2i + 1)
	end
	y
end

# ╔═╡ a74be92a-07aa-11eb-1716-9b44c98ae491
begin
	local fig
	fig = plot(sin, lw = 2, xlim = (-π, π), ylim = (-1.5, 1.5))
	for n in 0:4
		plot!(fig, x -> mysin(x, n), lw = 2, xlim = (-π, π), label = "DL(0, $(2n + 1))")
	end
	fig
end

# ╔═╡ Cell order:
# ╟─0162aa9a-f363-11ea-192a-8350869f18df
# ╟─d3f9af78-f35e-11ea-0812-1b1069364fba
# ╟─0ecdebd0-07a8-11eb-2504-c9c77f451a58
# ╟─07859f52-f360-11ea-0540-893a89f816ff
# ╟─7250c606-07a9-11eb-1bcd-a7e371200e58
# ╟─2937f110-f366-11ea-10b1-1d324a07f056
# ╟─9a139d6c-07a9-11eb-2f62-ed69c34cac69
# ╠═ecbbc1ca-07a9-11eb-012c-edb5b379dda1
# ╠═a088258c-07a9-11eb-0859-f135b5379d27
# ╠═f0aaa9e0-07a9-11eb-0b99-f355ceeed754
# ╠═c79e1fdc-07a9-11eb-34ff-d3cdf9916161
# ╠═a74be92a-07aa-11eb-1716-9b44c98ae491
