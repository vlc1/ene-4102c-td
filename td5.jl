### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 61a36274-028f-11eb-06c9-91ed9a2056f7
using LinearAlgebra, Kronecker

# ╔═╡ 4020216e-028f-11eb-147f-a3093365e0a5
n = (6, 4)

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
laplacian = kronecker(der[1], id[2]) + kronecker(id[1], der[2])

# ╔═╡ 4bdeab8e-048c-11eb-2323-59e5c6f74e73
b = rand(prod(n))

# ╔═╡ a561f8e6-048c-11eb-0206-61e13bc3e06c
x = laplacian \ b

# ╔═╡ Cell order:
# ╠═61a36274-028f-11eb-06c9-91ed9a2056f7
# ╠═4020216e-028f-11eb-147f-a3093365e0a5
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
