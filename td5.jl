### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 028e9c9a-08ac-11eb-0f5e-01125e8da26f
using Zygote

# ╔═╡ 52def07c-08a8-11eb-05fb-9f9caf480190
using BenchmarkTools

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

# ╔═╡ f6a57860-08a4-11eb-2363-2fc5d35b853b


# ╔═╡ 0d4d8be0-08aa-11eb-09e2-3738916f0581
θ(x) = x[1] ^ 3 * x[2] ^ 4

# ╔═╡ 90130a3a-08a5-11eb-0c26-ff20028c13a0
begin
	Δ(θ) = xy -> Zygote.hessian(θ, xy)
	left(θ) = y -> getindex(getindex(gradient(θ, [zero(y), y]), 1), 1)
	bottom(θ) = x -> getindex(getindex(gradient(θ, [x, zero(x)]), 1), 2)
	right(θ) = y -> θ([one(y), y])
	top(θ) = x -> θ([x, one(x)])
end

# ╔═╡ e5dca68a-08a8-11eb-10df-9f537313df55
left(θ)

# ╔═╡ e25ed0f8-08aa-11eb-110f-d32dfc467d40
gradient(θ, [1., 2.])

# ╔═╡ 024e6440-08a7-11eb-06c0-87709bc08bf6
function h(f, x, y)
	g₁(f, x, y) = getindex(gradient(f, x, y), 1)
	g₂(f, x, y) = getindex(gradient(f, x, y), 2)
	getindex(gradient(g₁, x, y), 1) + getindex(gradient(g₂, x, y), 2)
end

# ╔═╡ d8fa560c-08a7-11eb-1551-5712802dc20e
@btime Δ(f)([1., 2.])

# ╔═╡ 439a7b74-08a6-11eb-1b20-7f3941c2978d
h(f, 1., 2.)

# ╔═╡ ee9c5fba-08a6-11eb-016c-25374c0912d8
Δ(1., 2.)

# ╔═╡ cd4c9688-08a4-11eb-1dce-7b6bdf78f402
n = (4, 8)

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
b = rand(prod(n))

# ╔═╡ a561f8e6-048c-11eb-0206-61e13bc3e06c
#x = laplacian \ b

# ╔═╡ Cell order:
# ╠═028e9c9a-08ac-11eb-0f5e-01125e8da26f
# ╟─4b159bea-08ab-11eb-1a40-a768dfbc2051
# ╠═4a689a06-08a3-11eb-3461-fd2e89b7d99f
# ╠═ecc8d690-08a2-11eb-04ab-275138e29f23
# ╠═f6a57860-08a4-11eb-2363-2fc5d35b853b
# ╠═0d4d8be0-08aa-11eb-09e2-3738916f0581
# ╠═90130a3a-08a5-11eb-0c26-ff20028c13a0
# ╠═e5dca68a-08a8-11eb-10df-9f537313df55
# ╠═e25ed0f8-08aa-11eb-110f-d32dfc467d40
# ╠═024e6440-08a7-11eb-06c0-87709bc08bf6
# ╠═d8fa560c-08a7-11eb-1551-5712802dc20e
# ╠═52def07c-08a8-11eb-05fb-9f9caf480190
# ╠═439a7b74-08a6-11eb-1b20-7f3941c2978d
# ╠═ee9c5fba-08a6-11eb-016c-25374c0912d8
# ╠═cd4c9688-08a4-11eb-1dce-7b6bdf78f402
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
