### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 61a36274-028f-11eb-06c9-91ed9a2056f7
using LinearAlgebra

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

# ╔═╡ Cell order:
# ╠═61a36274-028f-11eb-06c9-91ed9a2056f7
# ╠═4020216e-028f-11eb-147f-a3093365e0a5
# ╠═4f6452bc-028f-11eb-3146-47029f8d6209
# ╠═26953684-0290-11eb-253f-6dc5b199298a
# ╠═5e08958a-028f-11eb-2562-79e5d34c1e5a
# ╠═abc7dc54-028f-11eb-0dc3-87a5587fa835
# ╠═79e6efae-028f-11eb-2e7a-7170dadf75c6
# ╠═e7b4f602-028f-11eb-39b2-bd42f7bf0e93
