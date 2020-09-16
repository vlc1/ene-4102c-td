### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 24926e0a-f7ec-11ea-0be8-c90fc7e13813
using Zygote, Plots

# ╔═╡ 4ae18622-f7ec-11ea-2f71-d5b166ff50fb
md"""
# Recherche de la racine d'une fonction

Nous avons vu lors de la séance 2 que les éléments ``y_1 \simeq y \left ( t_1 \right )``, ``y_2 \simeq y \left ( t_2 \right )``... de la solution numérique du problème de Cauchy
```math
\left \{ \begin{aligned}
\dot{y} \left ( t \right ) & = f \left [ t, y \left ( t \right ) \right ], \\
y \left ( 0 \right ) & = y_0
\end{aligned} \right .
```
par les méthodes à pas unique sont définis implicitement, c'est à dire comme racine de fonctions.

Dans le cas du schéma implicite d'Euler, par exemple, l'élément ``y_{n + 1}`` (``n = 0, \ldots N - 1``) est la racine de l'équation
```math
F_n \left ( x \right ) = 0
```
où (``\tau`` dénote le pas de temps)
```math
F_n \colon x \mapsto x - y_n - \tau f \left ( t_{n + 1}, x \right ).
```

L'objectif de cette première partie est de présenter comment résoudre cette équation numériquement dans le cas scalaire (``\forall n = 1, \ldots N``, ``y_n \in \mathbb{K}`` où ``\mathbb{K} = \mathbb{R}, \mathbb{C}``...).

## Différentiation algorithmique

"""

# ╔═╡ 63796c5e-f7ec-11ea-2ce7-11ff767d897f
md"""
## Méthode de Newton

"""

# ╔═╡ b847075e-f7f2-11ea-2995-090a48903ddc
function newton(f, x, p...)
    r = f(x, p...)
    while abs(r) > √eps(r)
        x -= r / gradient(f, x, p...)[1]
        r = f(x, p...)
    end
    x, r
end

# ╔═╡ 2b44a7a2-f7ee-11ea-3599-c7ac149c74ea
f(x) = x ^ 2 - 3x

# ╔═╡ 37e809a4-f7ee-11ea-1cbe-2b6b97b9fb28
gradient(f, π)

# ╔═╡ c7ae0958-f7ee-11ea-0844-256a0910a9b4
newton(f, π)

# ╔═╡ d0bb311e-f7f4-11ea-23c2-6df6cbd3b862
g(x, y, τ) = x - y + τ * x

# ╔═╡ e6f0ce26-f7f4-11ea-04dc-edb855a30bf0
newton(g, 1., 1., 0.125)

# ╔═╡ 7dc345ee-f7ec-11ea-138f-a1c6b81e0260
md"""
# Modèle et solution exacte

"""

# ╔═╡ 58162516-f7ec-11ea-3095-85bde6d71604
model(t, q, λ = -1) = λ * q

# ╔═╡ 3ab81476-f7f8-11ea-3633-09930c9cdffe
solution(t, λ = -1) = exp(λ * t)

# ╔═╡ 8a4674da-f7ec-11ea-2faf-4b332a41d7fc
md"""
# Schéma numérique

"""

# ╔═╡ a71a5b6c-f7ec-11ea-3881-cd909df3a068
explicit(f, x, y, τ, t) = x - y - τ * f(t, y)

# ╔═╡ aad46496-f7ec-11ea-3d35-5bf3365ec5e3
implicit(f, x, y, τ, t) = x - y - τ * f(t + τ, x)

# ╔═╡ ad6a1be2-f7ec-11ea-376b-3d1b8ad005c0
midpoint(f, x, y, τ, t) = x - y - τ * f(t + τ / 2, (x + y) / 2)

# ╔═╡ b417acac-f7ec-11ea-00cc-a1921aa7b239
trapezoidal(f, x, y, τ, t) = x - y - τ * (f(t, y) + f(t + τ, x)) / 2

# ╔═╡ b5f58c6a-f7ec-11ea-34eb-7beed843a51c
function rk2(f, x, y, τ, t)
    p = y + τ * f(t, y) / 2
    x - y - τ * f(t + τ / 2, p)
end

# ╔═╡ bdcbf0be-f7ec-11ea-3d8e-93cf936c03ff
function rk4(f, x, y, τ, t)
    k₁ = τ * f(t, y)
    k₂ = τ * f(t + τ / 2, y + k₁ / 2)
    k₃ = τ * f(t + τ / 2, y + k₂ / 2)
    k₄ = τ * f(t + τ, y + k₃)
    x - y - (k₁ + 2k₂ + 2k₃ + k₄) / 6
end

# ╔═╡ c87b46c0-f7ec-11ea-2918-d306ffd1c2bd
md"""
# Time stepping

"""

# ╔═╡ 2cd31668-f7ed-11ea-1b9b-d5c6ae89db19
struct Problem{F, G}
    scheme::F
    model::G
end

# ╔═╡ 2f44266c-f7ed-11ea-0add-37905cfaf54c
(this::Problem)(var...) = this.scheme(this.model, var...)

# ╔═╡ c5638ff8-f7ec-11ea-270f-cd22f82a2a23
function integrate(problem, y, t, τ, n)
    T = [t]
    Y = [y]

    for i in 1:n
        y, _ = newton(problem, y, y, τ, t)
        t += τ
        
        push!(Y, y)
        push!(T, t)
    end

    T, Y
end

# ╔═╡ 11baa86e-f7ed-11ea-3439-519c8a8fc61e
md"""
# Premier exemple

"""

# ╔═╡ 19b1b2e2-f7ed-11ea-2e66-c7391d4f7443
problem = Problem(implicit, model)

# ╔═╡ 48205f2a-f7ed-11ea-1574-f1ed0b1d7034
y, t = 1.0, 0.

# ╔═╡ 78044ff6-f7f7-11ea-1501-a354bed24082
τ, n = 0.1, 10

# ╔═╡ 6f6dac7c-f7f5-11ea-3726-45c4888b6d8a
T, Y = integrate(problem, y, t, τ, n)

# ╔═╡ ba191574-f7f5-11ea-29a5-639e8c561661
begin
	a, b = T[begin], T[end]
	c, d = minimum(Y), maximum(Y)
	fig = scatter(T, Y, xlims = (a, b), ylims = (1.1c - 0.1d, -0.1c + 1.1d))
	plot!(fig, a:(b - a) / 100:b, solution)
end

# ╔═╡ Cell order:
# ╠═24926e0a-f7ec-11ea-0be8-c90fc7e13813
# ╟─4ae18622-f7ec-11ea-2f71-d5b166ff50fb
# ╟─63796c5e-f7ec-11ea-2ce7-11ff767d897f
# ╠═b847075e-f7f2-11ea-2995-090a48903ddc
# ╠═2b44a7a2-f7ee-11ea-3599-c7ac149c74ea
# ╠═37e809a4-f7ee-11ea-1cbe-2b6b97b9fb28
# ╠═c7ae0958-f7ee-11ea-0844-256a0910a9b4
# ╠═d0bb311e-f7f4-11ea-23c2-6df6cbd3b862
# ╠═e6f0ce26-f7f4-11ea-04dc-edb855a30bf0
# ╟─7dc345ee-f7ec-11ea-138f-a1c6b81e0260
# ╠═58162516-f7ec-11ea-3095-85bde6d71604
# ╠═3ab81476-f7f8-11ea-3633-09930c9cdffe
# ╟─8a4674da-f7ec-11ea-2faf-4b332a41d7fc
# ╠═a71a5b6c-f7ec-11ea-3881-cd909df3a068
# ╠═aad46496-f7ec-11ea-3d35-5bf3365ec5e3
# ╠═ad6a1be2-f7ec-11ea-376b-3d1b8ad005c0
# ╠═b417acac-f7ec-11ea-00cc-a1921aa7b239
# ╠═b5f58c6a-f7ec-11ea-34eb-7beed843a51c
# ╠═bdcbf0be-f7ec-11ea-3d8e-93cf936c03ff
# ╟─c87b46c0-f7ec-11ea-2918-d306ffd1c2bd
# ╠═2cd31668-f7ed-11ea-1b9b-d5c6ae89db19
# ╠═2f44266c-f7ed-11ea-0add-37905cfaf54c
# ╠═c5638ff8-f7ec-11ea-270f-cd22f82a2a23
# ╟─11baa86e-f7ed-11ea-3439-519c8a8fc61e
# ╠═19b1b2e2-f7ed-11ea-2e66-c7391d4f7443
# ╠═48205f2a-f7ed-11ea-1574-f1ed0b1d7034
# ╠═78044ff6-f7f7-11ea-1501-a354bed24082
# ╠═6f6dac7c-f7f5-11ea-3726-45c4888b6d8a
# ╠═ba191574-f7f5-11ea-29a5-639e8c561661
