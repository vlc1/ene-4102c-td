### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 2aea5cfc-03bf-11eb-3922-7b05c68e795a
using DifferentialEquations, LinearAlgebra, Plots, LsqFit

# ╔═╡ 5d7fafec-0422-11eb-1062-710a53eb59f0
md"""
Versions [Pluto](https://github.com/vlc1/ene-4102c-td/blob/master/td4.jl) et [Jupyter](https://vlc1.github.io/ene-4102c/td4.ipynb) de ce notebook.

"""

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

1. Soit ``n`` le nombre de points. Quelle valeur de ``h`` utiliser pour obtenir une discrétisation d'ordre 2 ? Implémenter la fonction `spacing` qui retourne ``h`` en fonction de ``n``, et enfin `mesh` qui retourne un tableau contenant les abscisses de chacune des inconnues.

"""

# ╔═╡ 8bd6afce-04b5-11eb-39a4-817b45ff8848
phi() = 1 / √3

# ╔═╡ 65969a00-0419-11eb-3bd9-41931e3e837a
spacing(n) = 1 / (n + phi())

# ╔═╡ a8cfd4a4-0418-11eb-3d7d-257307b0d2bd
mesh(n) = [spacing(n) * (phi() + (j - 1)) for j in 1:n]

# ╔═╡ 103e3098-041b-11eb-2a89-855f41985558
md"""
2. Discrétiser l'équation, et l'implémenter à l'aide de deux fonctions, `laplacian` et `rhs`, invoquée depuis la fonction `numerical` ci-dessous.

"""

# ╔═╡ 2603c732-03bf-11eb-3c01-a3b6d7e04ef6
function laplacian(n)
	h = spacing(n)

	A = Tridiagonal(zeros.((n - 1, n, n - 1))...)

	# gauche
	A[1, 1] = 1 / (phi() + 1 / 2) / h ^ 2
	A[1, 2] = -1 / (phi() + 1 / 2) / h ^ 2

	# intérieur
	for j in 2:n - 1
		A[j, j - 1] = -1 / h ^ 2
		A[j, j] = 2 / h ^ 2
		A[j, j + 1] = -1 / h ^ 2
	end

	# droite
	A[n, n - 1] = -1 / h ^ 2
	A[n, n] = 2 / h ^ 2

	A
end

# ╔═╡ 5d9b193a-03c0-11eb-2b91-3f838aaf7c4a
# rhs = right hand side = second membre
function rhs(n, f, g₀, y₁)
	h, X = spacing(n), mesh(n)

	# source
	B = f.(X)

	# conditions aux limites
	B[begin] -= g₀ / (phi() + 1 / 2) / h
	B[end] += y₁ / h ^ 2

	B
end

# ╔═╡ a2321e98-041d-11eb-1a10-297553e4aee7
function numerical(n, f, g₀, y₁)
	A, B = laplacian(n), rhs(n, f, g₀, y₁)
	A \ B
end

# ╔═╡ 867ddb82-041b-11eb-3e88-43190945ccab
md"""
3. Définir la solution analytique du problème suivant (fonction `analytical`) :

"""

# ╔═╡ 70e982b8-03c0-11eb-283d-67e9be295243
begin
	f(x) = π ^ 2 * sinpi(x)
	g₀ = π
	y₁ = 0
	f, g₀, y₁
end

# ╔═╡ ea42b9e4-041b-11eb-02cb-c33e4890adb5
analytical(x) = sinpi(x)

# ╔═╡ 444fa398-041c-11eb-3365-8d70eed19687
md"""
4. Résoudre le problème pour ``n = 2``, ``4``... et visualiser les solutions numérique et analytique. Que constatez-vous ?

"""

# ╔═╡ fa72e35a-041d-11eb-0512-5b63e47199b5
N = [2 ^ i for i in 4:6]

# ╔═╡ a1c7017c-03c4-11eb-1425-6b42db5f8878
begin
	fig = plot()
	ϵ = Float64[]
	for n in N
		X, Y = mesh(n), numerical(n, f, g₀, y₁)
		push!(ϵ, norm(numerical(n, f, g₀, y₁) .- analytical.(X), Inf))
		scatter!(fig, X, Y, label = "n = $n")
	end
	plot!(fig, analytical, xlim = (0, 1), label = "analytical")
end

# ╔═╡ f6aa25e8-041e-11eb-307b-7b6c8411573c
md"""
5. Construire un autre problème en modifiant `f`, `g₀`, `y₁` et `analytical` de sorte à que la solution ne soit plus un polynôme d'ordre 2.
6. Calculer les erreurs ``L_1``, ``L_2`` et ``L_\infty`` pour plusieurs maillages et normes. Enfin, estimer l'ordre de convergence pour chacune des normes par la méthode des moindres carrés.

"""

# ╔═╡ Cell order:
# ╠═2aea5cfc-03bf-11eb-3922-7b05c68e795a
# ╟─5d7fafec-0422-11eb-1062-710a53eb59f0
# ╟─389f51d8-0318-11eb-1ada-19cfbeb88947
# ╟─e304b5b2-0324-11eb-0b67-b5b1c5a3ce36
# ╠═8bd6afce-04b5-11eb-39a4-817b45ff8848
# ╠═65969a00-0419-11eb-3bd9-41931e3e837a
# ╠═a8cfd4a4-0418-11eb-3d7d-257307b0d2bd
# ╟─103e3098-041b-11eb-2a89-855f41985558
# ╠═a2321e98-041d-11eb-1a10-297553e4aee7
# ╠═2603c732-03bf-11eb-3c01-a3b6d7e04ef6
# ╠═5d9b193a-03c0-11eb-2b91-3f838aaf7c4a
# ╟─867ddb82-041b-11eb-3e88-43190945ccab
# ╠═70e982b8-03c0-11eb-283d-67e9be295243
# ╠═ea42b9e4-041b-11eb-02cb-c33e4890adb5
# ╟─444fa398-041c-11eb-3365-8d70eed19687
# ╠═fa72e35a-041d-11eb-0512-5b63e47199b5
# ╠═a1c7017c-03c4-11eb-1425-6b42db5f8878
# ╟─f6aa25e8-041e-11eb-307b-7b6c8411573c
