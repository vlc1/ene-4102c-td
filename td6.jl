### A Pluto.jl notebook ###
# v0.12.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

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

# ╔═╡ 264df660-0fe8-11eb-0466-8daf53d0c2e1
# NE PAS MODIFIER
begin
	struct System{I,S,L,R,B,T}
		initial::I
		source::S
		left::L
		right::R
		bottom::B
		top::T
	end

	for field in fieldnames(System)
		eval(:($field(prob::System) = prob.$field))
	end
end

# ╔═╡ 4ef4ed6a-08d3-11eb-3fc3-4963a265809a
# NE PAS MODIFIER
function rhs(sys, n, t)
	h = spacing.(n)
	x, y = mesh.(n)

	# source
	b = map(Tuple.(CartesianIndices(n))) do (i, j)
		source(sys)(t, x[i], y[j])
	end

	# boundary conditions
	b[1, :] .-= left(sys).(t, y) / (ϕ() + 1 / 2) / h[1]
	b[end, :] .+= right(sys).(t, y) / h[1] ^ 2
	b[:, 1] .-= bottom(sys).(t, x) / (ϕ() + 1 / 2) / h[2]
	b[:, end] .+= top(sys).(t, x) / h[2] ^ 2

	reshape(b, prod(n))
end

# ╔═╡ 34048188-0fb2-11eb-3c97-795e8bc6ea22
# NE PAS MODIFIER
function initial(sys, n)
	x, y = mesh.(n)

	θ = map(Tuple.(CartesianIndices(n))) do (i, j)
		initial(sys)(x[i], y[j])
	end

	reshape(θ, prod(n))
end

# ╔═╡ e58f6d58-0fe6-11eb-1221-95bd801b1623
# MODIFIER
function integrate(sys, τ, m, n)
	t, z = zero(τ), initial(sys, n)

	T, Z = [t], [z]

	A = laplacian(n)
	I = UniformScaling(1.0)

	for i in 1:m
		b = rhs(sys, n, t + τ)

		z = (I - τ * A) \ (z + τ * b)
		t += τ

		push!(Z, z)
		push!(T, t)
	end

	T, Z
end

# ╔═╡ 949fe75a-0f8f-11eb-0371-e3f6cdc9856c
# MODIFIER
begin
	θ₀(x, y) = zero(x * y)
	ω(t, x, y) = zero(x)
	g₁(t, y) = y
	g₂(t, x) = x
	θ₁(t, y) = y
	θ₂(t, x) = x
end

# ╔═╡ 49543c2a-0feb-11eb-261f-8baf0b90c4ff
begin
	τ, m, n = 0.001, 10, (16, 16)
	sys = System(θ₀, ω, g₁, θ₁, g₂, θ₂)
	T, Z = integrate(sys, τ, m, n)
end

# ╔═╡ 82796ac2-0ff1-11eb-13ea-af8c0badf657
axes(Z)

# ╔═╡ 62aadebc-0ff1-11eb-3c2a-838b2c3a234b
@bind step Slider(axes(Z, 1))

# ╔═╡ 8ff89a58-0ff1-11eb-1d25-89abb1d7af84
begin
	x, y = mesh.(n)
	heatmap(x, y, reshape(Z[step], n...))
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
# ╠═264df660-0fe8-11eb-0466-8daf53d0c2e1
# ╠═4ef4ed6a-08d3-11eb-3fc3-4963a265809a
# ╠═34048188-0fb2-11eb-3c97-795e8bc6ea22
# ╠═e58f6d58-0fe6-11eb-1221-95bd801b1623
# ╠═949fe75a-0f8f-11eb-0371-e3f6cdc9856c
# ╠═49543c2a-0feb-11eb-261f-8baf0b90c4ff
# ╠═82796ac2-0ff1-11eb-13ea-af8c0badf657
# ╠═62aadebc-0ff1-11eb-3c2a-838b2c3a234b
# ╠═8ff89a58-0ff1-11eb-1d25-89abb1d7af84
