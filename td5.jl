### A Pluto.jl notebook ###
# v0.12.3

using Markdown
using InteractiveUtils

# ╔═╡ 028e9c9a-08ac-11eb-0f5e-01125e8da26f
begin
	using Zygote, StaticArrays, LinearAlgebra, Kronecker, LsqFit
	import Zygote.hessian
end

# ╔═╡ faaba6b8-0a08-11eb-391e-89ae75354ccf
md"""
Versions [Pluto](https://github.com/vlc1/ene-4102c-td/blob/master/td5.jl) et [Jupyter](https://vlc1.github.io/ene-4102c/td5/sujet.ipynb) de ce notebook.

"""

# ╔═╡ 4a689a06-08a3-11eb-3461-fd2e89b7d99f
md"""
# Method of Manufactured Solutions

On s'intéresse ici à la résolution du problème de conduction dans une plaque ``\left ( x, y \right ) \in \left [ 0, 1 \right ] ^ 2`` en régime stationnaire. On se propose donc de résoudre le problème suivant par la méthode des différences finies :
```math
0 = \Delta \theta \left ( x, y \right ) + \omega \left ( x, y \right )
```
où le Laplacien en coordonnées Cartésiennes, c'est à dire
```math
\Delta \theta = \left . \frac{\partial ^ 2 \theta}{\partial x ^ 2} \right \vert _ y +  \left . \frac{\partial ^ 2 \theta}{\partial y ^ 2} \right \vert _ x.
```

Le champs de température ``\theta`` sera soumis aux conditions aux limites suivantes :
```math
\left \{ \begin{aligned}
\partial_x \theta \left ( 0, y \right ) & = g_1 \left ( y \right ), \\
\theta \left ( 1, y \right ) & = \theta_1 \left ( y \right ),
\end{aligned} \right . \quad \mathrm{and} \quad \left \{ \begin{aligned}
\partial_y \theta \left ( x, 0 \right ) & = g_2 \left ( x \right ), \\
\theta \left ( x, 1 \right ) & = \theta_2 \left ( x \right ).
\end{aligned} \right .
```

On se propose d'employer la *Method of Manufactured Solutions* afin de vérifier l'implémentation. Étant donnée une solution ``\theta`` choisie à l'avance, l'idée est de

* Déterminer ``\omega``, ``g_1``, ``g_2``, ``\theta_1`` et ``\theta_2`` à partir de cette solution,
* Résoudre numériquement le problème et comparer solution numérique à la solution choisie afin d'obtenir l'ordre de la méthode.

1. Calculer  ``\omega``, ``g_1``, ``g_2``, ``\theta_1`` et ``\theta_2`` correspondant à la solution
```math
\theta \colon \left ( x, y \right ) \mapsto \sin \left ( \pi x \right ) \cos \left (\pi y \right ).
```
2. Les fonctions `Δ`, `left`, `bottom`, `right` et `top` prennent ``\theta`` comme argument et retournent respectivement ``\omega``, ``g_1``, ``g_2``, ``\theta_1`` et ``\theta_2`` automatiquement grâce à la bibliothèque d'*Algorithmic Differentiation* `Zygote.jl`. Vérifier graphiquement que ces fonctions retournent bien le résultat attendu.

"""

# ╔═╡ fed3b534-08d0-11eb-3442-65c59de227c2
# Question 2
begin
	θ(x, y) = sinpi(x) * cospi(y)
	ω(x, y) = -2π ^ 2 * sinpi(x) * cospi(y)
	g₁(y) = π * cospi(y)
	g₂(x) = zero(x)
	θ₁(y) = zero(y)
	θ₂(x) = -sinpi(x)
end

# ╔═╡ 90130a3a-08a5-11eb-0c26-ff20028c13a0
# Question 2 - NE PAS MODIFIER
begin
	θ(x::AbstractVector) = θ(x[1], x[2])
	Δ(f, x, y) = first(hessian(f, SVector(x, y))) + last(hessian(f, SVector(x, y)))
	left(f, y) = first(gradient(f, zero(y), y))
	bottom(f, x) = last(gradient(f, x, zero(x)))
	right(f, y) = f(one(y), y)
	top(f, x) = f(x, one(x))
end

# ╔═╡ 9b479164-0a03-11eb-3c93-a1e929c5fd2e
md"""
3. Les fonctions `phi`, `spacing` et `mesh` du TD4 sont rappelées ci-dessous. Les utiliser pour implémenter la fonction `exact` qui évalue ``θ`` à chaque point du maillage.

"""

# ╔═╡ ecc8d690-08a2-11eb-04ab-275138e29f23
# Question 3 - NE PAS MODIFIER
phi() = 1 / √3

# ╔═╡ 2a88081a-0a02-11eb-1f89-adc0728ee515
# Question 3 - NE PAS MODIFIER
spacing(n) = 1 / (n + phi())

# ╔═╡ 90791206-0a02-11eb-1d91-19f5454706a7
# Question 5 - NE PAS MODIFIER
mesh(n) = [spacing(n) * (phi() + (j - 1)) for j in 1:n]

# ╔═╡ b669d1c8-0a03-11eb-0d79-d7232b76bc79
# Question 3
function exact(f, n)
	x, y = mesh.(n)

	map(Tuple.(CartesianIndices(n))) do (i, j)
		f(x[i], y[j])
	end
end

# ╔═╡ 62510d90-0a01-11eb-12a7-a73a23bd3917
md"""
4. La fonction `laplacian` du TD4 est rappelée plus bas. S'en servir pour implémenter la discrétisation du laplacien en deux dimensions (voir le [produit de Kronecker](https://fr.wikipedia.org/wiki/Produit_de_Kronecker) implémenté par la fonction `kron` de `Kronecker.jl`).

"""

# ╔═╡ 747cca1e-0a02-11eb-100c-73e14b996048
# Question 4
function laplacian(n::NTuple{2, Int})
	id = Diagonal.(fill.(-1.0, n))
	fd = laplacian.(n)
	kron(id[2], fd[1]) + kron(fd[2], id[1])
end

# ╔═╡ 8e3fb30a-0a00-11eb-3d1a-cf1d430f9524
# Question 4 - NE PAS MODIFIER
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

# ╔═╡ fc6b1422-0a01-11eb-2db4-6f24865d02d4
md"""
5. Se servir de `mesh`, `Δ`, `left`, `bottom`, `right` et `top` pour implémenter la fonction `rhs` qui assemble le second membre.

"""

# ╔═╡ 4ef4ed6a-08d3-11eb-3fc3-4963a265809a
# Question 5
function rhs(f, n)
	h = spacing.(n)
	x, y = mesh.(n)

	# source
	b = map(Tuple.(CartesianIndices(n))) do (i, j)
		Δ(f, x[i], y[j])
	end

	# boundary conditions
	b[1, :] .+= left.(f, y) / (phi() + 1 / 2) / h[1]
	b[end, :] .-= right.(f, y) / h[1] ^ 2
	b[:, 1] .+= bottom.(f, x) / (phi() + 1 / 2) / h[2]
	b[:, end] .-= top.(f, x) / h[2] ^ 2

	b
end

# ╔═╡ 39436fd2-0a04-11eb-21b0-f9be513a8477
md"""
6. La fonction `numerical` définie plus bas résoudre numériquement le problème. Implémenter la fonction `error` ci-dessous qui calcule la norme de la différence entre la solution numérique et la solution analytique. Vérifier que celle-ci est nulle lorsque ``\theta`` est une fonction quadratique.

"""

# ╔═╡ f8e7f1de-08db-11eb-3e1b-4177c637d838
# Question 6 - NE PAS MODIFIER
function numerical(f, n)
	A = laplacian(n)
	b = reshape(rhs(f, n), prod(n))

	x = A \ b

	reshape(x, n...)
end

# ╔═╡ d4660498-0a04-11eb-2f37-7310b0a79832
# Question 6
error(f, n) = norm(numerical(f, n) - exact(f, n))

# ╔═╡ f9b732a4-0a05-11eb-3f07-df7a037b4699
md"""
7. Pour la solution de la question 1, vérifier le schéma est bien d'ordre deux.

"""

# ╔═╡ 52280146-0a68-11eb-0f38-abc19c00370a
begin
	n = [(2 ^ i, 2 ^ i) for i in 4:6]
	ϵ, h = spacing.(first.(n)), error.(θ, n)
end

# ╔═╡ Cell order:
# ╟─faaba6b8-0a08-11eb-391e-89ae75354ccf
# ╠═028e9c9a-08ac-11eb-0f5e-01125e8da26f
# ╟─4a689a06-08a3-11eb-3461-fd2e89b7d99f
# ╠═fed3b534-08d0-11eb-3442-65c59de227c2
# ╠═90130a3a-08a5-11eb-0c26-ff20028c13a0
# ╟─9b479164-0a03-11eb-3c93-a1e929c5fd2e
# ╠═b669d1c8-0a03-11eb-0d79-d7232b76bc79
# ╠═ecc8d690-08a2-11eb-04ab-275138e29f23
# ╠═2a88081a-0a02-11eb-1f89-adc0728ee515
# ╠═90791206-0a02-11eb-1d91-19f5454706a7
# ╟─62510d90-0a01-11eb-12a7-a73a23bd3917
# ╠═747cca1e-0a02-11eb-100c-73e14b996048
# ╠═8e3fb30a-0a00-11eb-3d1a-cf1d430f9524
# ╟─fc6b1422-0a01-11eb-2db4-6f24865d02d4
# ╠═4ef4ed6a-08d3-11eb-3fc3-4963a265809a
# ╟─39436fd2-0a04-11eb-21b0-f9be513a8477
# ╠═d4660498-0a04-11eb-2f37-7310b0a79832
# ╠═f8e7f1de-08db-11eb-3e1b-4177c637d838
# ╟─f9b732a4-0a05-11eb-3f07-df7a037b4699
# ╠═52280146-0a68-11eb-0f38-abc19c00370a
