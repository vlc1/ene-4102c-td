### A Pluto.jl notebook ###
# v0.12.4

using Markdown
using InteractiveUtils

# ╔═╡ 028e9c9a-08ac-11eb-0f5e-01125e8da26f
using Plots, LinearAlgebra, Kronecker, PlutoUI

# ╔═╡ faaba6b8-0a08-11eb-391e-89ae75354ccf
md"""
Versions [Pluto](https://github.com/vlc1/ene-4102c-td/blob/master/td6.jl) et [Jupyter](https://vlc1.github.io/ene-4102c/td6/sujet.ipynb) de ce notebook.

"""

# ╔═╡ 4a689a06-08a3-11eb-3461-fd2e89b7d99f
md"""
# Conduction dans une plaque en régime instationnaire

On reprend ici la résolution du problème de la conduction dans une plaque ``\left ( x, y \right ) \in \left [ 0, 1 \right ] ^ 2``, à ceci près que cette fois-ci le régime est supposé instationnaire. On se propose donc de résoudre numériquement le problème suivant :
```math
\rho c_p \left . \frac{\partial \theta}{\partial t} \right \vert _ {x, y} = \lambda \Delta \theta \left ( t, x, y \right ) + \omega \left ( t, x, y \right )
```
où on écrit encore une fois le Laplacien en coordonnées Cartésiennes, c'est à dire
```math
\Delta \theta = \left . \frac{\partial ^ 2 \theta}{\partial x ^ 2} \right \vert _ y +  \left . \frac{\partial ^ 2 \theta}{\partial y ^ 2} \right \vert _ x.
```

Le champs de température ``\theta`` sera soumis aux conditions aux limites suivantes (noter la dépendance temporelle) :
```math
\left \{ \begin{aligned}
\partial_x \theta \left ( t, 0, y \right ) & = g_1 \left ( t, y \right ), \\
\theta \left ( t, 1, y \right ) & = \theta_1 \left ( t, y \right ),
\end{aligned} \right . \quad \mathrm{and} \quad \left \{ \begin{aligned}
\partial_y \theta \left ( t, x, 0 \right ) & = g_2 \left ( t, x \right ), \\
\theta \left ( t, x, 1 \right ) & = \theta_2 \left ( t, x \right ).
\end{aligned} \right .
```

À ces conditions aux limites, il convient d'ajouter la condition initiale
```math
\theta \left ( 0, x, y \right ) = \theta_0 \left ( x, y \right ).
```

On se propose de résoudre ce problème par

1. La méthode des différences finies pour ce qui est des dérivées spatiales,
1. Les méthodes d'Euler (explicite et implicites) pour ce qui est de la dérivée temporelle.

Pour ce faire, on utilisera les codes ci-dessous, modifiés à partir de ceux des séances précédentes. On commencera par implémenter les fonctions `ω`, `g₁`, `g₂`, `θ₁`, `θ₂` et `θ₀` à la fin de ce notebook.

"""

# ╔═╡ ecc8d690-08a2-11eb-04ab-275138e29f23
# NE PAS MODIFIER
ϕ() = 1 / √3

# ╔═╡ 2a88081a-0a02-11eb-1f89-adc0728ee515
# NE PAS MODIFIER
spacing(n) = 1 / (n + ϕ())

# ╔═╡ 90791206-0a02-11eb-1d91-19f5454706a7
# NE PAS MODIFIER
mesh(n) = [spacing(n) * (ϕ() + (j - 1)) for j in 1:n]

# ╔═╡ 8e3fb30a-0a00-11eb-3d1a-cf1d430f9524
# NE PAS MODIFIER
function laplacian(n)
	h = spacing(n)

	A = Tridiagonal(zeros.((n - 1, n, n - 1))...)

	# gauche
	A[1, 1] = 1 / (ϕ() + 1 / 2) / h ^ 2
	A[1, 2] = -1 / (ϕ() + 1 / 2) / h ^ 2

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

# ╔═╡ 747cca1e-0a02-11eb-100c-73e14b996048
# NE PAS MODIFIER
function laplacian(n::NTuple{2, Int})
	id = Diagonal.(fill.(-1.0, n))
	fd = laplacian.(n)
	kron(id[2], fd[1]) + kron(fd[2], id[1])
end

# ╔═╡ 4ef4ed6a-08d3-11eb-3fc3-4963a265809a
# NE PAS MODIFIER
function rhs(ω, g₁, g₂, θ₁, θ₂, n, t)
	h = spacing.(n)
	x, y = mesh.(n)

	# source
	b = map(Tuple.(CartesianIndices(n))) do (i, j)
		ω(t, x[i], y[j])
	end

	# boundary conditions
	b[1, :] .-= g₁.(t, y) / (ϕ() + 1 / 2) / h[1]
	b[end, :] .+= θ₁.(t, y) / h[1] ^ 2
	b[:, 1] .-= g₂.(t, x) / (ϕ() + 1 / 2) / h[2]
	b[:, end] .+= θ₂.(t, x) / h[2] ^ 2

	reshape(b, prod(n))
end

# ╔═╡ 34048188-0fb2-11eb-3c97-795e8bc6ea22
# NE PAS MODIFIER
function initial(θ₀, n)
	x, y = mesh.(n)

	θ = map(Tuple.(CartesianIndices(n))) do (i, j)
		θ₀(x[i], y[j])
	end

	reshape(θ, prod(n))
end

# ╔═╡ e58f6d58-0fe6-11eb-1221-95bd801b1623
# MODIFIER
function integrate(ω, g₁, g₂, θ₁, θ₂, z, t, τ, n, m)
	T, Z = [t], [z]

	for i in 1:m
		z = z
		t += τ

		push!(Z, z)
		push!(T, t)
	end

	T, Z
end

# ╔═╡ 949fe75a-0f8f-11eb-0371-e3f6cdc9856c
# MODIFIER
begin
	ω(t, x, y) = zero(x)
	g₁(t, y) = zero(y)
	g₂(t, x) = zero(x)
	θ₁(t, y) = zero(y)
	θ₂(t, x) = zero(x)
	θ₀(x, y) = zero(x)
end

# ╔═╡ Cell order:
# ╟─faaba6b8-0a08-11eb-391e-89ae75354ccf
# ╠═028e9c9a-08ac-11eb-0f5e-01125e8da26f
# ╟─4a689a06-08a3-11eb-3461-fd2e89b7d99f
# ╠═ecc8d690-08a2-11eb-04ab-275138e29f23
# ╠═2a88081a-0a02-11eb-1f89-adc0728ee515
# ╠═90791206-0a02-11eb-1d91-19f5454706a7
# ╠═8e3fb30a-0a00-11eb-3d1a-cf1d430f9524
# ╠═747cca1e-0a02-11eb-100c-73e14b996048
# ╠═4ef4ed6a-08d3-11eb-3fc3-4963a265809a
# ╠═34048188-0fb2-11eb-3c97-795e8bc6ea22
# ╠═e58f6d58-0fe6-11eb-1221-95bd801b1623
# ╠═949fe75a-0f8f-11eb-0371-e3f6cdc9856c
