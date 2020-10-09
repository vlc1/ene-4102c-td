### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 24926e0a-f7ec-11ea-0be8-c90fc7e13813
using Zygote, Plots

# ╔═╡ d61e1ea6-f924-11ea-00dc-794c93177d22
md"""
Versions [Pluto](https://github.com/vlc1/ene-4102c-td/blob/master/td2.jl) et [Jupyter](https://vlc1.github.io/ene-4102c/td2/sujet.ipynb) de ce notebook.

"""

# ╔═╡ 4ae18622-f7ec-11ea-2f71-d5b166ff50fb
md"""
# Recherche de la racine d'une fonction

Nous avons vu lors de la séance 3 que les éléments
```math
\begin{aligned}
y_1 & \simeq y \left ( t_1 \right ), \\
y_2 & \simeq y \left ( t_2 \right ), \\
& \ldots
\end{aligned}
```
de la solution numérique du problème de Cauchy
```math
\left \{ \begin{aligned}
\dot{y} \left ( t \right ) & = f \left [ t, y \left ( t \right ) \right ], \\
y \left ( 0 \right ) & = y_0
\end{aligned} \right .
```
par les méthodes à pas unique sont définis implicitement, c'est à dire comme racines de fonctions.

Dans le cas du schéma implicite d'Euler, par exemple, l'élément ``y_{n + 1}`` (``n = 0, \ldots N - 1``) est la racine de l'équation
```math
F_n \left ( x \right ) = 0
```
où (``\tau`` dénote le pas de temps)
```math
F_n \colon x \mapsto x - y_n - \tau f \left ( t_{n + 1}, x \right ).
```

L'objectif de cette première partie est de présenter, dans le cas scalaire (``\forall n = 1, \ldots N``, ``y_n \in \mathbb{K}`` où ``\mathbb{K} = \mathbb{R}, \mathbb{C}``...), une méthode de résolution numérique de ces équations : la méthode de Newton.

## Différentiation algorithmique

Avant toute chose, la méthode de Newton requiert l'évaluation du gradient d'une fonction. Considérons la fonction de deux variables
```math
f \colon \left ( x, y \right ) \mapsto x ^ 2 e ^ {-y}.
```

Cette fonction étant relativement simple, nous pouvons bien sûr la dériver par rapport à ``x`` ou ``y``, pour en obtenir le gradient
```math
\nabla f \colon \left ( x, y \right ) \mapsto \left ( 2x e ^ {-y}, -x ^ 2 e ^ {-y} \right ).
```

Bien que tout à fait acceptable, on préfère souvent à cette approche "à la main" une autre approche, dite **différentiation algorithmique**. L'idée est que l'utilisateur se contente de fournir la fonction `f` seule, le code faisant le reste pour en calculer le gradient (ou une valeur approchée).

Cette étape peut être réalisée grâce à plusieurs techniques, notamment par différences finies, *dual numbers*, *source-to-source transformation*... Nous nous contenterons ici d'utiliser le package `Zygote.jl`. Dans le cas de la fonction précédente, le gradient au point ``\left ( 1, 2 \right )`` est tout simplement évalué comme suit.

"""

# ╔═╡ 24bed800-f858-11ea-095b-45829493711d
gradient((x, y) -> x ^ 2 * exp(-y), 1., 2.)

# ╔═╡ 36288188-f858-11ea-0bcb-158f5d177506
md"""
**Question** -- Soit la fonction
```math
\left ( x, y, z \right ) \mapsto \sin \left ( x \right ) \cos \left ( y - z \right ).
```
En dériver le gradient et comparer la valeur au point ``\left ( \pi / 2, \pi / 3, \pi / 4 \right )`` avec celle obtenue par la fonction `gradient` de `Zygote.jl`.

"""

# ╔═╡ 63796c5e-f7ec-11ea-2ce7-11ff767d897f
md"""
## Méthode de Newton

À partir d'une estimation initiale ``x _ 0`` de la solution ``x _ \star`` de l'équation
```math
f \left ( x _ \star \right ) = 0
```
(où ``f`` dénote une fonction ``\mathcal{C} ^ 1``), la méthode de Newton produit une séquence ``\left ( x _ k \right )`` qui, sous certaines conditions, tend vers la solution ``x _ \star``.

La séquence est définie comme suit :
```math
\forall k \in \mathbb{N}, \quad x _ {k + 1} = x _ k - \frac{f \left ( x _ k \right )}{f' \left ( x _ k \right )}.
```

"""

# ╔═╡ b847075e-f7f2-11ea-2995-090a48903ddc
function newton(f, x, p...)
    r = f(x, p...)
    while abs(r) > √eps(r)
        x -= r / first(gradient(f, x, p...))
        r = f(x, p...)
    end
    x, r
end

# ╔═╡ a0fb1520-f8b7-11ea-1cc2-41cd1c29c12f
md"""
Pour une interprétation graphique de la méthode, voir l'animation suivante (source : [wikipedia](https://en.wikipedia.org/wiki/Newton%27s_method)).

![Méthode de Newton](https://upload.wikimedia.org/wikipedia/commons/e/e0/NewtonIteration_Ani.gif)

On remarquera la nécessité d'évaluer la dérivée de ``f`` à chaque itération. La fonction suivante implémente cet algorithme. Elle cherche la racine par rapport au premier argument de la fonction ``f``, dénoté ici ``x``. ``p`` représente ici un ou plusieurs paramètres (implémentée ici par une **fonction variadique**, c'est à dire une fonction qui accepte un nombre arbitraire d'arguments.

**Exemple** -- Soit la famille de fonctions
```math
f_\alpha \colon x \mapsto x ^ \alpha - 2.
```
Cherchons sa racine dans le cas ``\alpha = 3`` à partir de l'estimation initiale ``x = 2``.

"""

# ╔═╡ 7aa3438e-f8b7-11ea-1c1e-b937e725e9fa
newton((x, α) -> x ^ α - 2, 2., 3)

# ╔═╡ b9da1834-f8b7-11ea-016b-6dd917e044a3
md"""
Cette valeur est évidemment à comparer à la racine cubique de ``2``, soit
```julia
julia> cbrt(2.)
1.2599210498948732
```

On note enfin que la fonction `newton` retourne un `Tuple` de deux éléments :

1. Le premier est la dernière estimation de ``x _ \star``,
1. Le second est la valeur de la fonction en ce point.

"""

# ╔═╡ 7dc345ee-f7ec-11ea-138f-a1c6b81e0260
md"""
# Modèle et solution exacte

On se concentre pour l'instant sur le modèle linéaire homogène (décroissance radioactive), pour lequel le second membre de l'EDO s'écrit
```math
f \colon \left ( t, y \right ) \mapsto \lambda y.
```

La solution exacte s'écrit alors sous la forme :
```math
y \colon t \mapsto \exp \left ( \lambda t \right ) y_0.
```

Ces deux fonctions sont implémentées dans les cellules suivantes.

"""

# ╔═╡ 58162516-f7ec-11ea-3095-85bde6d71604
linear(t, y, λ = -1) = λ * y

# ╔═╡ 3ab81476-f7f8-11ea-3633-09930c9cdffe
solution(t, λ = -1, y₀ = 1., t₀ = 0.) = exp(λ * (t - t₀)) * y₀

# ╔═╡ 8a4674da-f7ec-11ea-2faf-4b332a41d7fc
md"""
# Schéma numérique

On rappelle que lors de la séance 3, quatre schémas numériques ont été présentés :
```math
y_{n + 1} - y_n - \tau f \left ( t_n, y_n \right ) = 0 \quad \text{(Euler explicite)},
```
```math
y_{n + 1} - y_n - \tau f \left ( t_{n + 1}, y_{n + 1} \right ) = 0 \quad \text{(Euler implicite)},
```
```math
y_{n + 1} - y_n - \tau \frac{f \left ( t_n, y_n \right ) + f \left ( t_{n + 1}, y_{n + 1} \right )}{2} = 0 \quad \text{(Méthode des trapèzes)},
```
```math
y_{n + 1} - y_n - \tau f \left ( \frac{t_n + t_{n + 1}}{2}, \frac{y_n + y_{n + 1}}{2} \right ) = 0 \quad \text{(Méthode du point milieu)}.
```

**Question** -- En suivant l'exemple suivant (`explicit`), implémenter les fonctions `implicit`, `trapezoidal` et `midpoint` dont la racine est ``y_{n + 1}``. On préservera l'ordre des paramètres, au nombre de 3, à savoir

* `y` -- la solution précédente, ``y _ n`` ;
* `τ` -- le pas de temps, ``\tau`` ;
* `t` -- l'instant précédent, ``t _ n``.

"""

# ╔═╡ a71a5b6c-f7ec-11ea-3881-cd909df3a068
explicit(f, x, y, τ, t) = x - y - τ * f(t, y)

# ╔═╡ 0b0ec174-f9bb-11ea-2165-b156550d1c7e
implicit(f, x, y, τ, t) = x - y - τ * f(t + τ, x)

# ╔═╡ 17a3d276-f9bb-11ea-179a-75e3d70c266a
trapezoidal(f, x, y, τ, t) = x - y - τ * (f(t, y) + f(t + τ, y)) / 2

# ╔═╡ 2f1b501e-f9bb-11ea-1e93-43ad4f994902
midpoint(f, x, y, τ, t) = x - y - τ * f(t + τ / 2, (x + y) / 2)

# ╔═╡ 605783fa-f8bc-11ea-293a-77b05130e32c
md"""
Il existe bien sûr un grand nombre de schémas à pas unique, ceux de Runge-Kutta étant parmi les plus connus, notamment le schéma explicite d'ordre 2
```math
\begin{aligned}
y _ * & = y _ n + \frac{\tau}{2} f \left ( t_n, y_n \right ), \\
y _ {n + 1} & = y _ n + \tau f \left ( t_n + \frac{\tau}{2}, y _ * \right )
\end{aligned}
```
et celui d'ordre 4
```math
\begin{aligned}
k _ 1 & = \tau f \left ( t _ n, y _ n \right ), \\
k _ 2 & = \tau f \left ( t _ n + \frac{\tau}{2}, y_n + \frac{k _ 1}{2} \right ), \\
k _ 3 & = \tau f \left ( t _ n + \frac{\tau}{2}, y_n + \frac{k _ 2}{2} \right ), \\
k _ 4 & = \tau f \left ( t _ n + \tau, y_n + k _ 3 \right )
\end{aligned}
```
et enfin
```math
y _ {n + 1} = y _ n + \frac{k _ 1 + 2 k _ 2 + 2 k _ 3 + k _ 4}{6}.
```

**Question** -- À partir de la fonction `rk2` suivante, implémenter la fonction `rk4`.

"""

# ╔═╡ b5f58c6a-f7ec-11ea-34eb-7beed843a51c
function rk2(f, x, y, τ, t)
    z = y + τ * f(t, y) / 2
    x - y - τ * f(t + τ / 2, z)
end

# ╔═╡ a78344e4-f9ac-11ea-34f0-8fe50a08cb15
function rk4(f, x, y, τ, t)
	k₁ = τ * f(t, y)
	k₂ = τ * f(t + τ / 2, y + k₁ / 2)
	k₃ = τ * f(t + τ / 2, y + k₂ / 2)
	k₄ = τ * f(t + τ, y + k₃)
	x - y - (k₁ + 2k₂ + 2k₃ + k₄) / 6
end

# ╔═╡ c87b46c0-f7ec-11ea-2918-d306ffd1c2bd
md"""
# Intégration temporelle

Il reste à présent à assembler les différentes briques, à savoir les modèles et les schémas implémentés. Pour ce faire, on définit un nouveau **composite type**, que l'on nommera `Problem` (la première lettre du nom des types est, par convention, en majuscule).

"""

# ╔═╡ 2cd31668-f7ed-11ea-1b9b-d5c6ae89db19
struct Problem{F, G}
    scheme::F
    model::G
end

# ╔═╡ aec5fd4a-f8bf-11ea-1cb2-77441c10648f
md"""
On souhaite enfin que les instances de ce nouveau type puissent être appelées comme une fonction de plusieurs arguments. Ces instances seront ensuite passées comme premier argument à la fonction `newton`. On dira que la cellule suivante rend les objets de type `Problem` **callable**.

"""

# ╔═╡ 2f44266c-f7ed-11ea-0add-37905cfaf54c
(this::Problem)(var...) = this.scheme(this.model, var...)

# ╔═╡ d7b451b4-f8bf-11ea-07ac-2df552a0d944
md"""
Il ne reste plus qu'à implémenter la bouche d'intégration temporelle. Étant donnés

* Une condition et un instant initiaux, `y` et `t`,
* Un pas de temps et un nombre d'itérations, `τ` et `n`,

la fonction `integrate` retourne un `Tuple` de deux vecteurs, le premier contenant les instants
```math
t_0 \quad t_1 \quad \cdots \quad t_N
```
et le second la solution numérique, à savoir
```math
y_0 \quad y_1 \quad \cdots \quad y_N.
```

"""

# ╔═╡ c5638ff8-f7ec-11ea-270f-cd22f82a2a23
function integrate(problem, y, t, τ, s)
    T = [t]
    Y = [y]

    while t < (1 - √eps(t)) * s
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

Les cellules suivantes illustrent la résolution du problème linéaire par la méthode explicit d'Euler.

"""

# ╔═╡ 19b1b2e2-f7ed-11ea-2e66-c7391d4f7443
problem = Problem(explicit, linear)

# ╔═╡ 48205f2a-f7ed-11ea-1574-f1ed0b1d7034
y, t = 1.0, 0.

# ╔═╡ 78044ff6-f7f7-11ea-1501-a354bed24082
τ, s = 0.25, 1.

# ╔═╡ 6f6dac7c-f7f5-11ea-3726-45c4888b6d8a
T, Y = integrate(problem, y, t, τ, s)

# ╔═╡ f08317e2-f8c0-11ea-11d8-7d486eb7e057
md"""
La solution numérique peut être visualisée et comparée à la solution exacte grâce aux commandes suivantes.

"""

# ╔═╡ ba191574-f7f5-11ea-29a5-639e8c561661
begin
	a, b = first(T), last(T)
	c, d = minimum(Y), maximum(Y)
	fig = scatter(T, Y, xlims = (a, b), ylims = (1.1c - 0.1d, -0.1c + 1.1d))
	plot!(fig, a:(b - a) / 100:b, solution)
end

# ╔═╡ 107eb128-f8c1-11ea-0ffc-c3d57849b589
md"""
**Question** -- Augmenter progressivement le pas de temps d'intégration de `Problem(explicit, linear)`. Que remarque t'on ?

**Question** -- Tester différents schémas avec différents pas de temps. Qu'observe t'on ?

**Question** -- Calculer l'erreur à un temps donné (par exemple, ``T = 1``) en fonction du pas en temps pour chacun des schémas implémentés. Visualiser ces données sur un graphique en échelle logarithmique.

"""

# ╔═╡ d4637d80-f8c1-11ea-1f7f-df462373ca2d
md"""
# Au delà du cas linéaire

Tout l'intérêt de la méthode de Newton est qu'on peut bien sûr s'affranchir de la linéarité du modèle (même si, pour rester simple, ce notebook se limite au cas scalaire).

**Question** -- Utiliser ou imaginer un modèle non-linéaire et en visualiser la solution.

"""

# ╔═╡ 91a712b2-f8bd-11ea-3b8c-1bfbd521d29a
nonlinear(t, y, λ = -1) = λ * y ^ 2 / (1 + t)

# ╔═╡ 46222752-f91a-11ea-372e-2dde47c81add
md"""
# Au delà du cas scalaire

On se propose à présent d'utiliser une bibliothèque existante, `DifferentialEquations.jl`, pour résoudre l'[équation de prédation de Lotka-Volterra](https://fr.wikipedia.org/wiki/%C3%89quations_de_pr%C3%A9dation_de_Lotka-Volterra) :

> En mathématiques, les équations de prédation de Lotka-Volterra, que l'on désigne aussi sous le terme de "modèle proie-prédateur", sont un couple d'équations différentielles non-linéaires du premier ordre, et sont couramment utilisées pour décrire la dynamique de systèmes biologiques dans lesquels un prédateur et sa proie interagissent. Elles ont été proposées indépendamment par Alfred James Lotka en 1925 et Vito Volterra en 1926.

Le système d'équations s'écrit :
```math
\left \{ \begin{aligned}
\dot{x} \left ( t \right ) & = x \left ( t \right ) \left [ \alpha - \beta y \left ( t \right ) \right ], \\
\dot{y} \left ( t \right ) & = y \left ( t \right ) \left [ \delta x \left ( t \right ) - \gamma \right ]
\end{aligned} \right .
```
où

* ``t`` est le temps ;
* ``x \left ( t \right )`` est l'effectif des proies à l'instant ``t`` ;
* ``y \left ( t \right )`` est l'effectif des prédateurs à l'instant ``t``.

Les paramètres suivants enfin caractérisent les interactions entre les deux espèces :
```math
\alpha = 0.1, \quad \beta = 0.003, \quad \gamma = 0.06, \quad \delta = 0.0012.
```

**Question** -- Résoudre l'équation de Lotka-Volterra avec la bibliothèque `DifferenialEquations.jl` et de sa [documentation](https://diffeq.sciml.ai/stable/tutorials/ode_example/).

"""

# ╔═╡ Cell order:
# ╟─d61e1ea6-f924-11ea-00dc-794c93177d22
# ╠═24926e0a-f7ec-11ea-0be8-c90fc7e13813
# ╟─4ae18622-f7ec-11ea-2f71-d5b166ff50fb
# ╠═24bed800-f858-11ea-095b-45829493711d
# ╟─36288188-f858-11ea-0bcb-158f5d177506
# ╟─63796c5e-f7ec-11ea-2ce7-11ff767d897f
# ╠═b847075e-f7f2-11ea-2995-090a48903ddc
# ╟─a0fb1520-f8b7-11ea-1cc2-41cd1c29c12f
# ╠═7aa3438e-f8b7-11ea-1c1e-b937e725e9fa
# ╟─b9da1834-f8b7-11ea-016b-6dd917e044a3
# ╟─7dc345ee-f7ec-11ea-138f-a1c6b81e0260
# ╠═58162516-f7ec-11ea-3095-85bde6d71604
# ╠═3ab81476-f7f8-11ea-3633-09930c9cdffe
# ╟─8a4674da-f7ec-11ea-2faf-4b332a41d7fc
# ╠═a71a5b6c-f7ec-11ea-3881-cd909df3a068
# ╠═0b0ec174-f9bb-11ea-2165-b156550d1c7e
# ╠═17a3d276-f9bb-11ea-179a-75e3d70c266a
# ╠═2f1b501e-f9bb-11ea-1e93-43ad4f994902
# ╟─605783fa-f8bc-11ea-293a-77b05130e32c
# ╠═b5f58c6a-f7ec-11ea-34eb-7beed843a51c
# ╠═a78344e4-f9ac-11ea-34f0-8fe50a08cb15
# ╟─c87b46c0-f7ec-11ea-2918-d306ffd1c2bd
# ╠═2cd31668-f7ed-11ea-1b9b-d5c6ae89db19
# ╟─aec5fd4a-f8bf-11ea-1cb2-77441c10648f
# ╠═2f44266c-f7ed-11ea-0add-37905cfaf54c
# ╟─d7b451b4-f8bf-11ea-07ac-2df552a0d944
# ╠═c5638ff8-f7ec-11ea-270f-cd22f82a2a23
# ╟─11baa86e-f7ed-11ea-3439-519c8a8fc61e
# ╠═19b1b2e2-f7ed-11ea-2e66-c7391d4f7443
# ╠═48205f2a-f7ed-11ea-1574-f1ed0b1d7034
# ╠═78044ff6-f7f7-11ea-1501-a354bed24082
# ╠═6f6dac7c-f7f5-11ea-3726-45c4888b6d8a
# ╟─f08317e2-f8c0-11ea-11d8-7d486eb7e057
# ╠═ba191574-f7f5-11ea-29a5-639e8c561661
# ╟─107eb128-f8c1-11ea-0ffc-c3d57849b589
# ╟─d4637d80-f8c1-11ea-1f7f-df462373ca2d
# ╠═91a712b2-f8bd-11ea-3b8c-1bfbd521d29a
# ╟─46222752-f91a-11ea-372e-2dde47c81add
