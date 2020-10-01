### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 6f564db6-031e-11eb-10df-dfdb07cacd7d
using DifferentialEquations

# ╔═╡ 2aea5cfc-03bf-11eb-3922-7b05c68e795a
using LinearAlgebra

# ╔═╡ 9da78850-03c4-11eb-011c-975615ab606e
using Plots

# ╔═╡ 389f51d8-0318-11eb-1ada-19cfbeb88947
md"""
# Résolution du problème de Blasius

On se propose de résoudre le problème de Blasius de deux façons :

1. Par la méthode de tir (*shooting*), implémentée dans la bibliothèque `DifferentialEquations.jl`,
1. En exploitant une symétrie du problème.

Le problème s'écrit :
```math
y''' \left ( x \right ) + y \left ( x \right ) y'' \left ( x \right ) = 0
```
où la variable dépendante ``y`` vérifie les conditions aux limites
```math
\left \{ \begin{aligned}
y \left ( 0 \right ) & = 0, \\
y' \left ( 0 \right ) & = 0, \\
y' \left ( \infty \right ) & = 1.
\end{aligned} \right .
```

L'une des difficultés rencontrées lors de la résolution numérique de ce problème est qu'il est défini sur un domaine semi-infini. On se propose ici d'utiliser la méthode dite *domain truncation*, qui consiste à résoudre l'équation sur un domaine fini mais suffisamment grand pour ne pas influencer la solution numérique, ce qui doit bien sûr être vérifié *a posteriori*.

## Problème aux limites et méthode de tir

1. Quel est le degré de cette EDO ? La réécrire sous la forme d'un système d'EDO d'ordre 1.
2. À partir de ce [tutoriel](https://diffeq.sciml.ai/stable/tutorials/bvp_example/), résoudre l'équation de Blasius sur le domaine ``\left [ 0, 10000 \right ]``. On pourra utiliser les fonctions suivantes.

```julia
function blasius!(z, y, p, x)
	z[1] = y[2]
	z[2] = y[3]
	z[3] = -y[1] * y[3]
end
```
```julia
function bc!(r, y, p, x)
	r[1] = y[begin][1]
	r[2] = y[begin][2]
	r[3] = y[end][2] - 1
end
```

## Problème de Cauchy et symétrie

3. Montrer que l'équation de Blasius et les conditions en ``0`` sont inchangées par la transformation
```math
\begin{aligned}
\overline{y} & \leftarrow c y, \\
\overline{x} & \leftarrow x / c.
\end{aligned}
```

On suppose maintenant que ``\overline{y}`` est soumis à la même équation que ``y``, à l'exception de la condition en $\infty$ qui est remplacée par :
```math
\overline{y}'' \left ( 0 \right ) = 1.
```
On définit alors
```math
\alpha = \overline{y}' \left ( \infty \right ).
```

4. Montrer que le choix de la constante ``c = \sqrt{\alpha}`` mène à
```math
y' \left ( \infty \right ) = 1.
```
5. Résoudre ``\overline{y}`` numériquement à partir de ce [tutoriel](https://diffeq.sciml.ai/stable/tutorials/ode_example/) et en déduire la solution de Blasius intégrant le problème numériquement une seule fois.

## Vérification

6. Vérifier qu'augmenter la taille du domaine d'intégration n'affecte pas la solution numérique.

"""

# ╔═╡ 27071d3a-031f-11eb-337e-ef188b1ff882
u₀ = [0., 0., 1.]

# ╔═╡ 5ba72e6c-0322-11eb-0216-05672a01f4c5
blasius!(z, y, p, x) = z .= [y[2], y[3], -y[1] * y[3]]

# ╔═╡ 3426ac9e-0324-11eb-2238-5d7522a221db
function bc!(r, y, p, x)
	r[1] = y[begin][1]
	r[2] = y[begin][2]
	r[3] = y[end][2] - 1
end

# ╔═╡ 59010f78-0324-11eb-0053-01a5b40429df
bvp = BVProblem(blasius!, bc!, [0., 0., 1.], (0., 100.))

# ╔═╡ 8351f7d8-0324-11eb-3015-d1cf06eff4be
sol = solve(bvp, Shooting(Tsit5()))

# ╔═╡ e304b5b2-0324-11eb-0b67-b5b1c5a3ce36
md"""
# Conduction stationnaire dans une barre

On se propose de résoudre l'équation
```math
0 = y'' \left ( x \right ) + f \left ( x \right )
```
avec les conditions
```math
\left \{ \begin{aligned}
y ' \left ( 0 \right ) & = g_0, \\
y \left ( 1 \right ) & = y_1.
\end{aligned} \right .
```

``f`` représente une source (par exemple un dépot de chaleur par laser).

1. Soit ``n`` le nombre de points. Quelle valeur de ``h`` utiliser pour obtenir une discrétisation d'ordre 2 ?
1. Discrétiser l'équation.

"""

# ╔═╡ 70e982b8-03c0-11eb-283d-67e9be295243
f(x) = -one(x)

# ╔═╡ b4d14eac-03c0-11eb-31c4-331908f6e6c7
y₁ = 0.5

# ╔═╡ b8183b66-03c0-11eb-22f2-b1055db79a0f
g₀ = 0.

# ╔═╡ cd38ee84-03be-11eb-13f6-fb65706e7293
n = 8

# ╔═╡ d69f1fac-03be-11eb-18f3-f7a77ce68ca0
h = 1 / (1 / √3 + n)

# ╔═╡ 872d5e9c-03bf-11eb-1ca2-272d8a149d11
H = [i == 1 ? h / √3 : h for i in 1:n + 1]

# ╔═╡ d888a030-03bf-11eb-17b4-2b88fb01c7f5
X = [sum(H[1:i-1]) for i in 2:n + 1]

# ╔═╡ 2603c732-03bf-11eb-3c01-a3b6d7e04ef6
begin
	A = Tridiagonal(zeros.((n - 1, n, n - 1))...)

	A[1, 1], A[1, 2] = 6 / (3 + 2 * √3) / h ^ 2, -6 / (3 + 2 * √3) / h ^ 2

	for i in 2:n - 1
		A[i, i - 1], A[i, i], A[i, i + 1] = -1 / h ^ 2, 2 / h ^ 2, - 1/ h ^ 2
	end

	A[n, n - 1], A[n, n] = -1 / h ^ 2, 2 / h ^ 2
end

# ╔═╡ 5d9b193a-03c0-11eb-2b91-3f838aaf7c4a
begin
	B = f.(X)
	B[1] -= 6g₀ / (3 + 2 * √3) / h
	B[n] += y₁ / h ^ 2
end

# ╔═╡ 6377b8b2-03bf-11eb-138d-0700381fcceb
A \ B

# ╔═╡ bac0d47a-03c2-11eb-3531-79e6d80379be
X .^ 2 ./ 2

# ╔═╡ a1c7017c-03c4-11eb-1425-6b42db5f8878


# ╔═╡ Cell order:
# ╟─389f51d8-0318-11eb-1ada-19cfbeb88947
# ╠═6f564db6-031e-11eb-10df-dfdb07cacd7d
# ╠═27071d3a-031f-11eb-337e-ef188b1ff882
# ╠═5ba72e6c-0322-11eb-0216-05672a01f4c5
# ╠═3426ac9e-0324-11eb-2238-5d7522a221db
# ╠═59010f78-0324-11eb-0053-01a5b40429df
# ╠═8351f7d8-0324-11eb-3015-d1cf06eff4be
# ╟─e304b5b2-0324-11eb-0b67-b5b1c5a3ce36
# ╠═2aea5cfc-03bf-11eb-3922-7b05c68e795a
# ╠═70e982b8-03c0-11eb-283d-67e9be295243
# ╠═b4d14eac-03c0-11eb-31c4-331908f6e6c7
# ╠═b8183b66-03c0-11eb-22f2-b1055db79a0f
# ╠═cd38ee84-03be-11eb-13f6-fb65706e7293
# ╠═d69f1fac-03be-11eb-18f3-f7a77ce68ca0
# ╠═872d5e9c-03bf-11eb-1ca2-272d8a149d11
# ╠═d888a030-03bf-11eb-17b4-2b88fb01c7f5
# ╠═2603c732-03bf-11eb-3c01-a3b6d7e04ef6
# ╠═5d9b193a-03c0-11eb-2b91-3f838aaf7c4a
# ╠═6377b8b2-03bf-11eb-138d-0700381fcceb
# ╠═bac0d47a-03c2-11eb-3531-79e6d80379be
# ╠═9da78850-03c4-11eb-011c-975615ab606e
# ╠═a1c7017c-03c4-11eb-1425-6b42db5f8878
