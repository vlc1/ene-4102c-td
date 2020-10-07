### A Pluto.jl notebook ###
# v0.12.0

using Markdown
using InteractiveUtils

# ╔═╡ 028e9c9a-08ac-11eb-0f5e-01125e8da26f
begin
	using Zygote, StaticArrays, LinearAlgebra, Kronecker
	import Zygote.hessian
end

# ╔═╡ 4b159bea-08ab-11eb-1a40-a768dfbc2051
md"""
# Method of manufactured solution

"""

# ╔═╡ 4a689a06-08a3-11eb-3461-fd2e89b7d99f
md"""
# Conduction dans une plaque en régime stationnaire

On se propose de résoudre numériquement le problème suivant en coordonnée cartésienne :
```math
0 = \Delta \theta \left ( x, y \right ) + f \left ( x, y \right ),
```
où le champs de température ``\theta`` est soumis aux conditions aux limites suivantes :
```math
\left \{ \begin{aligned}
\partial_x \theta \left ( 0, y \right ) & = g_1 \left ( y \right ), \\
\theta \left ( 1, y \right ) & = \theta_1 \left ( y \right ),
\end{aligned} \right . \quad \mathrm{and} \quad \left \{ \begin{aligned}
\partial_y \theta \left ( x, 0 \right ) & = g_2 \left ( x \right ), \\
\theta \left ( x, 1 \right ) & = \theta_2 \left ( x \right ).
\end{aligned} \right .
```

"""

# ╔═╡ ecc8d690-08a2-11eb-04ab-275138e29f23
begin
	phi() = 1 / √3
	spacing(n) = 1 / (n + phi())
	mesh(n) = [spacing(n) * (phi() + (j - 1)) for j in 1:n]
	#=
	function numerical(n, f, g₀, y₁)
		A, B = laplacian(n), rhs(n, f, g₀, y₁)
		A \ B
	end
	=#
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
end

# ╔═╡ fed3b534-08d0-11eb-3442-65c59de227c2
begin
	θ(x, y) = x ^ 2 / 2 + y ^ 2 / 2 + sin(y) * cos(x)
	θ(x::AbstractVector) = θ(x[1], x[2])
end

# ╔═╡ 90130a3a-08a5-11eb-0c26-ff20028c13a0
begin
	Δ(f) = (x, y) -> first(hessian(f, SVector(x, y))) + last(hessian(f, SVector(x, y)))
	left(f) = y -> first(first(gradient(f, SVector(zero(y), y))))
	bottom(f) = x -> last(first(gradient(f, SVector(x, zero(x)))))
	right(f) = y -> f(one(y), y)
	top(f) = x -> f(x, one(x))
end

# ╔═╡ cd4c9688-08a4-11eb-1dce-7b6bdf78f402
n = (4, 8)

# ╔═╡ f5570038-08da-11eb-25a4-ab003a6a7e16
h = spacing.(n)

# ╔═╡ 37d5bf6a-08d3-11eb-2676-b946778d1eea
x, y = mesh.(n)

# ╔═╡ 4ef4ed6a-08d3-11eb-3fc3-4963a265809a
begin
	# source
	b = map(Tuple.(CartesianIndices(n))) do (i, j)
		Δ(θ)(x[i], y[j])
	end

	# boundary conditions
	b[1, :] .+= left(θ).(y) / (phi() + 1 / 2) / h[1]
	b[end, :] .-= right(θ).(y) / h[1] ^ 2
	b[:, 1] .+= bottom(θ).(x) / (phi() + 1 / 2) / h[2]
	b[:, end] .-= top(θ).(x) / h[2] ^ 2
end

# ╔═╡ 8c5f6f3e-08da-11eb-3546-f9e75010c2dc
begin
	id = Diagonal.(fill.(-1.0, n))
	fd = laplacian.(n)
	A = kron(id[2], fd[1]) + kron(fd[2], id[1])
end

# ╔═╡ f8e7f1de-08db-11eb-3e1b-4177c637d838
numerical = reshape(A \ reshape(b, prod(n)), n...)

# ╔═╡ f1510d38-08dc-11eb-256b-5709d1bef09e
exact = map(Tuple.(CartesianIndices(n))) do (i, j)
	θ(x[i], y[j])
end

# ╔═╡ 7db5ab58-08dd-11eb-3f77-2f628e94ca40
norm(numerical - exact)

# ╔═╡ Cell order:
# ╠═028e9c9a-08ac-11eb-0f5e-01125e8da26f
# ╟─4b159bea-08ab-11eb-1a40-a768dfbc2051
# ╟─4a689a06-08a3-11eb-3461-fd2e89b7d99f
# ╠═ecc8d690-08a2-11eb-04ab-275138e29f23
# ╠═fed3b534-08d0-11eb-3442-65c59de227c2
# ╠═90130a3a-08a5-11eb-0c26-ff20028c13a0
# ╠═cd4c9688-08a4-11eb-1dce-7b6bdf78f402
# ╠═f5570038-08da-11eb-25a4-ab003a6a7e16
# ╠═37d5bf6a-08d3-11eb-2676-b946778d1eea
# ╠═4ef4ed6a-08d3-11eb-3fc3-4963a265809a
# ╠═8c5f6f3e-08da-11eb-3546-f9e75010c2dc
# ╠═f8e7f1de-08db-11eb-3e1b-4177c637d838
# ╠═f1510d38-08dc-11eb-256b-5709d1bef09e
# ╠═7db5ab58-08dd-11eb-3f77-2f628e94ca40
