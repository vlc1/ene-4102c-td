### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ e8196eec-fe3a-11ea-2e52-6736a89588a2
using Plots, Zygote, LsqFit, Roots

# ╔═╡ e926a51e-ff0b-11ea-03e3-5369116e7cae
md"""
Versions [Pluto](https://github.com/vlc1/ene-4102c-td/blob/master/td3.jl) et [Jupyter](https://vlc1.github.io/ene-4102c/td3.ipynb) de ce notebook.

"""

# ╔═╡ 2740f61e-fe3e-11ea-33a6-251b4bf2e579
md"""
# Domaine de stabilité

$(begin
	local fig
	fig = contour(
		-3:0.01:.5,
		-3:0.01:3,
		(x, y) -> (1 + x) ^ 2 + y ^ 2,
		levels = [1],
		aspect_ratio = :equal,
		color = :blue,
		lw = 2,
		title = "Domaines de stabilité des méthodes EE et RK2"
	)
	contour!(
		fig,
		-3:0.01:.5,
		-3:0.01:3,
		(x, y) -> (1 + x + (x ^ 2 - y ^ 2) / 2) ^ 2 + (y * (1 + x)) ^ 2,
		levels = [1],
		aspect_ratio = :equal,
		color = :red,
		lw = 2
	)
end)

On a vu en cours que les domaines de stabilité des schémas explicit d'Euler et RK2 sont respectivement
```math
\left \{ z \in \mathbb{C} \left  / \left \vert 1 + z \right \vert \le 1 \right . \right \}
```
et
```math
\left \{ z \in \mathbb{C} \left  / \left \vert 1 + z + z ^ 2 / 2 \right \vert \le 1 \right . \right \}.
```

1. Définir la raison ``\overline{\sigma} \left ( z \right )`` (``z = \lambda \tau``) pour la méthode de RK4.
1. Visualiser le domaine de stabilité de RK4 grâce à la fonction `contour` de la bibliothèque `Plots.jl`.
1. Dans le cas ``\lambda \in \mathbb{R} ^ -``, trouver le pas de temps maximal de RK4 (en fonction de ``\lambda``). On utilisera la fonction `find_zero` de la bibliothèque `Roots.jl` (à installer).
1. Vérifier qu'au delà de ce pas de temps, cette méthode devient instable.

"""

# ╔═╡ ca7f0732-fe46-11ea-2f6a-cd77ccce2e3a
md"""
# Ordre de convergence

On rappelle que l'erreur d'un schéma à l'instant $t _ n$, définie par
```math
\epsilon = y_n - y \left ( t_n \right ),
```
peut s'écrire
```math
\epsilon = C \tau ^ p
```
où ``p`` dénote l'**ordre** de la méthode.

1. Montrer que ``\ln \epsilon`` est une fonction affine de ``\ln \tau``.
1. En réutilisant le code du TD2 (voir ci-dessous), pour chaque schéma, calculer l'erreur à un temps donné (par exemple, ``s = 1``) en utilisant plusieurs pas de temps. On pourra s'inspirer du code suivant :
```julia
Δ = [0.125 / 2 ^ i for i in 5:-1:1]
methods = [explicit, implicit, midpoint, rk2, rk4]
errors = Dict(method => Float64[] for method in methods)
for method in methods
	pb = Problem(method, linear)
	for δ in Δ
		T, Y = integrate(pb, y, t, δ, s)
		append!(errors[method], abs(Y[end] - solution(T[end])))
	end
end
```
3. Visualiser l'erreur en fonction du pas de temps grâce aux fonctions `scatter`/`scatter!` en échelle logarithmique (`scale = :log`).
4. Utiliser le logarithme de ces valeurs ainsi que la méthode des moindres carrés, implémentée par la fonction `curve_fit` de `LsqFit.jl` (à installer), pour estimer les paramètres de la fonction suivante : `e(τ, p) = p[1] .+ p[2] .* τ`.

"""

# ╔═╡ ae942a96-feff-11ea-030c-5f9dfb04c404
function newton(f, x, p...)
    r = f(x, p...)
    while abs(r) > √eps(r)
        x -= r / first(gradient(f, x, p...))
        r = f(x, p...)
    end
    x, r
end

# ╔═╡ ba74296c-feff-11ea-2254-c97edc9bdd20
struct Problem{F, G}
    scheme::F
    model::G
end

# ╔═╡ ceb61f98-feff-11ea-3cab-93a03f609676
(this::Problem)(var...) = this.scheme(this.model, var...)

# ╔═╡ efd147b6-feff-11ea-12d0-c3e2be6b5bd5
function integrate(problem, y, t, τ, s)
    T = [t]
    Y = [y]

	while t < (1 - √eps(t)) * s
        y, _ = newton(problem, y, y, τ, t)
        t += τ
        
        push!(Y, y)
        push!(T, t)
    end

    T, Y
end

# ╔═╡ 78bfaf38-ff00-11ea-1293-ed87bb52db18
explicit(f, x, y, τ, t) = x - y - τ * f(t, y)

# ╔═╡ 7a85c9ae-ff00-11ea-00b3-71c88caa094b
implicit(f, x, y, τ, t) = x - y - τ * f(t + τ, x)

# ╔═╡ 7e893fcc-ff00-11ea-05fe-35dd7bd3c73c
trapezoidal(f, x, y, τ, t) = x - y - τ * (f(t, y) + f(t + τ, y)) / 2

# ╔═╡ 8214cf30-ff00-11ea-2143-b10461e0e230
midpoint(f, x, y, τ, t) = x - y - τ * f(t + τ / 2, (x + y) / 2)

# ╔═╡ 85646d94-ff00-11ea-124c-619e03f113d9
function rk2(f, x, y, τ, t)
    z = y + τ * f(t, y) / 2
    x - y - τ * f(t + τ / 2, z)
end

# ╔═╡ 893b7186-ff00-11ea-0689-a725c0e27ec6
function rk4(f, x, y, τ, t)
	k₁ = τ * f(t, y)
	k₂ = τ * f(t + τ / 2, y + k₁ / 2)
	k₃ = τ * f(t + τ / 2, y + k₂ / 2)
	k₄ = τ * f(t + τ, y + k₃)
	x - y - (k₁ + 2k₂ + 2k₃ + k₄) / 6
end

# ╔═╡ 95617548-ff00-11ea-2d7f-7fafc8887d98
linear(t, y, λ = -1) = λ * y

# ╔═╡ 970e31b2-ff00-11ea-2c5d-efc0a71dabfd
solution(t, λ = -1, y₀ = 1., t₀ = 0.) = exp(λ * (t - t₀)) * y₀

# ╔═╡ 6dcdad44-ff00-11ea-1e93-7f4f6d70e55d
problem = Problem(rk2, linear)

# ╔═╡ b5c87ac0-ff00-11ea-04a5-c95672fe2a61
y, t = 1.0, 0.

# ╔═╡ b85216ea-ff00-11ea-266d-ebcccf7327c8
τ, s = 0.1, 1.

# ╔═╡ c1a56c4a-ff00-11ea-0298-717ac761f22d
T, Y = integrate(problem, y, t, τ, s)

# ╔═╡ Cell order:
# ╟─e926a51e-ff0b-11ea-03e3-5369116e7cae
# ╠═e8196eec-fe3a-11ea-2e52-6736a89588a2
# ╟─2740f61e-fe3e-11ea-33a6-251b4bf2e579
# ╟─ca7f0732-fe46-11ea-2f6a-cd77ccce2e3a
# ╠═ae942a96-feff-11ea-030c-5f9dfb04c404
# ╠═ba74296c-feff-11ea-2254-c97edc9bdd20
# ╠═ceb61f98-feff-11ea-3cab-93a03f609676
# ╠═efd147b6-feff-11ea-12d0-c3e2be6b5bd5
# ╠═78bfaf38-ff00-11ea-1293-ed87bb52db18
# ╠═7a85c9ae-ff00-11ea-00b3-71c88caa094b
# ╠═7e893fcc-ff00-11ea-05fe-35dd7bd3c73c
# ╠═8214cf30-ff00-11ea-2143-b10461e0e230
# ╠═85646d94-ff00-11ea-124c-619e03f113d9
# ╠═893b7186-ff00-11ea-0689-a725c0e27ec6
# ╠═95617548-ff00-11ea-2d7f-7fafc8887d98
# ╠═970e31b2-ff00-11ea-2c5d-efc0a71dabfd
# ╠═6dcdad44-ff00-11ea-1e93-7f4f6d70e55d
# ╠═b5c87ac0-ff00-11ea-04a5-c95672fe2a61
# ╠═b85216ea-ff00-11ea-266d-ebcccf7327c8
# ╠═c1a56c4a-ff00-11ea-0298-717ac761f22d
