### A Pluto.jl notebook ###
# v0.12.0

using Markdown
using InteractiveUtils

# ╔═╡ 028e9c9a-08ac-11eb-0f5e-01125e8da26f
begin
	using Zygote, StaticArrays
	import Zygote.hessian
end

# ╔═╡ 61a36274-028f-11eb-06c9-91ed9a2056f7
using LinearAlgebra, Kronecker

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
	function numerical(n, f, g₀, y₁)
		A, B = laplacian(n), rhs(n, f, g₀, y₁)
		A \ B
	end
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
	θ(x, y) = x ^ 2 / 2 + y ^ 2 / 2
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

# ╔═╡ 2b4c801a-08d5-11eb-1b05-3d4d10ba82b8
function (x, y)
	x + y
end

# ╔═╡ cd4c9688-08a4-11eb-1dce-7b6bdf78f402
n = (4, 8)

# ╔═╡ 37d5bf6a-08d3-11eb-2676-b946778d1eea
x, y = mesh.(n)

# ╔═╡ 4ef4ed6a-08d3-11eb-3fc3-4963a265809a
begin
	b = map(Tuple.(CartesianIndices(n))) do (i, j)
		Δ(θ)(x[i], y[j])
	end
end

# ╔═╡ b9a88d00-08d4-11eb-3a12-4de04966d8e1
CartesianIndices(n)

# ╔═╡ c1c967e8-08d4-11eb-1dcf-df221c74b480
Tuple(index)

# ╔═╡ 801d5214-08d4-11eb-2e5d-dbdfb45d357d
#size(b)

# ╔═╡ f30d9adc-08d3-11eb-297d-77c19a5432b1
begin
	foo = [2, 3, 4, 5]
	bar = [6 7 8 9]
	typeof.((foo, bar))
end

# ╔═╡ d57252f8-08a4-11eb-340f-cda87449ccee
laplacian.(n)

# ╔═╡ 4f6452bc-028f-11eb-3146-47029f8d6209
c, d = -1 .* ones.(n .- 1), 2 .* ones.(n)

# ╔═╡ 26953684-0290-11eb-253f-6dc5b199298a
A = Diagonal.(ones.(n))

# ╔═╡ 5e08958a-028f-11eb-2562-79e5d34c1e5a
B = Tridiagonal.(c, d, c)

# ╔═╡ abc7dc54-028f-11eb-0dc3-87a5587fa835
kron(A[1], B[2]) + kron(B[1], A[2])

# ╔═╡ 79e6efae-028f-11eb-2e7a-7170dadf75c6
kron(A...)

# ╔═╡ e7b4f602-028f-11eb-39b2-bd42f7bf0e93
UniformScaling(4)

# ╔═╡ a4d30098-048a-11eb-38e3-5bfdfa4792ac
id = @. Diagonal(ones(n))

# ╔═╡ b7e78fe6-048a-11eb-20f7-31bf91e086d3
der = @. Tridiagonal(ones(n - 1), -2ones(n), ones(n - 1))

# ╔═╡ ea8ab52c-048a-11eb-106e-bd3852a996d6
#kronecker(id[1], id[2], der[3])

# ╔═╡ 110e3124-048b-11eb-2e0e-77dcd19bf0ca
#kronecker(der[1], id[2], id[3])

# ╔═╡ a75408fc-048b-11eb-2c27-7712e459f9f9
#laplacian = kronecker(der[1], id[2], id[3]) + kronecker(id[1], der[2], id[3]) + kronecker(id[1], id[2], der[3])

# ╔═╡ fbe58ed8-048b-11eb-2c27-2588a9488baf
#laplacian = kronecker(der[1], id[2]) + kronecker(id[1], der[2])

# ╔═╡ 4bdeab8e-048c-11eb-2323-59e5c6f74e73
#b = rand(prod(n))

# ╔═╡ a561f8e6-048c-11eb-0206-61e13bc3e06c
#x = laplacian \ b

# ╔═╡ Cell order:
# ╠═028e9c9a-08ac-11eb-0f5e-01125e8da26f
# ╟─4b159bea-08ab-11eb-1a40-a768dfbc2051
# ╟─4a689a06-08a3-11eb-3461-fd2e89b7d99f
# ╠═ecc8d690-08a2-11eb-04ab-275138e29f23
# ╠═fed3b534-08d0-11eb-3442-65c59de227c2
# ╠═90130a3a-08a5-11eb-0c26-ff20028c13a0
# ╠═2b4c801a-08d5-11eb-1b05-3d4d10ba82b8
# ╠═cd4c9688-08a4-11eb-1dce-7b6bdf78f402
# ╠═37d5bf6a-08d3-11eb-2676-b946778d1eea
# ╠═4ef4ed6a-08d3-11eb-3fc3-4963a265809a
# ╠═b9a88d00-08d4-11eb-3a12-4de04966d8e1
# ╠═c1c967e8-08d4-11eb-1dcf-df221c74b480
# ╠═801d5214-08d4-11eb-2e5d-dbdfb45d357d
# ╠═f30d9adc-08d3-11eb-297d-77c19a5432b1
# ╠═d57252f8-08a4-11eb-340f-cda87449ccee
# ╠═61a36274-028f-11eb-06c9-91ed9a2056f7
# ╠═4f6452bc-028f-11eb-3146-47029f8d6209
# ╠═26953684-0290-11eb-253f-6dc5b199298a
# ╠═5e08958a-028f-11eb-2562-79e5d34c1e5a
# ╠═abc7dc54-028f-11eb-0dc3-87a5587fa835
# ╠═79e6efae-028f-11eb-2e7a-7170dadf75c6
# ╠═e7b4f602-028f-11eb-39b2-bd42f7bf0e93
# ╠═a4d30098-048a-11eb-38e3-5bfdfa4792ac
# ╠═b7e78fe6-048a-11eb-20f7-31bf91e086d3
# ╠═ea8ab52c-048a-11eb-106e-bd3852a996d6
# ╠═110e3124-048b-11eb-2e0e-77dcd19bf0ca
# ╠═a75408fc-048b-11eb-2c27-7712e459f9f9
# ╠═fbe58ed8-048b-11eb-2c27-2588a9488baf
# ╠═4bdeab8e-048c-11eb-2323-59e5c6f74e73
# ╠═a561f8e6-048c-11eb-0206-61e13bc3e06c
