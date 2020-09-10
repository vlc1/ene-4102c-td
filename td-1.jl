### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 9df1a0de-f33e-11ea-276a-5175cf87a652
md"""
# Configuration

## Étape 1 : installation de Julia 1.5.1

Depuis le site <https://julialang.org/downloads/>, télécharger la dernière version stable, à savoir Julia 1.5.1 en sélectionnant le bon système d'exploitation (Linux x86, Mac, Windows...).

Lancer enfin l'installation.

## Étape 2 : installation et démarrage de Julia

Le programme "Julia 1.5.1" doit maintenant être installé sur votre ordinateur. Pour le démarrer, chercher l'exécutable et lancer l'installation. Vous devriez voir apparaitre la console Julia, aussi appelée *Julia REPL*.
"""

# ╔═╡ 8e0646ee-f347-11ea-22f6-efa2515444b3
md"""
# Première utilisation

Pour une première prise en main de `Pluto`,

1. modifier et exécuter les cellules suivantes,
2. ajouter quelques cellules...

"""

# ╔═╡ dfb0b044-f347-11ea-1dee-2d161f8aa0c8
a = 1 + 3

# ╔═╡ e8c0f304-f347-11ea-20c2-7f3c5b9327f7
b = 3a

# ╔═╡ bfe43eb2-f347-11ea-257a-e5ebe3e91eb9
function myexp(x, n = 1)
	y = zero(x)
	for i in 0:n
		y += x ^ i / factorial(i)
	end
	y
end

# ╔═╡ ddea46f4-f347-11ea-3714-b5212290ca66
myexp(1., 4), exp(1.), Float64(ℯ)

# ╔═╡ 0bb68dfc-f347-11ea-1d1b-37be8fcd1cff
md"""
# Configuration et utilisation de `Pluto`

## Étape 1 : installation de la bibliothèque

L'étape suivante est l'installation de la bibliothèque [`Pluto`](https://github.com/fonsp/Pluto.jl) grâce au *package manager*.

Pour ce faire, ouvrir la console puis entrer `]` (fermeture de crochet) pour accéder au mode `Pkg` :
```julia
julia> ]
(@v1.5) pkg>
```

La line devient bleue et le *prompt* change. Pour installer `Pluto`, exécuter la commande suivante :
```julia
(@v1.5) pkg> add Pluto
```

L'installation devrait se terminer après quelques minutes.

## Étape 2 : installation du navigateur

L'exécution des notebooks `Pluto` nécessite un navigateur. Pour un rendu optimal, utiliser [Mozilla Firefox](https://www.mozilla.org/firefox/) ou [Google Chrome](https://www.google.com/chrome/).

## Étape 3 : démarrage

Depuis la console Julia, exécuter les commandes suivantes :
```julia
julia> using Pluto
julia> Pluto.run()
```

La console vous propose à présent de vous rendre à l'URL <http://localhost:1234/> (ou quelque chose comme ça) qu'il vous suffit d'ouvrir dans votre navigateur.

La page d'accueil de `Pluto` devrait maintenant s'ouvrir dans votre navigateur.

## Étape 4 : ouverture d'un notebook

Il s'agit maintenant d'ouvrir d'un notebook (par exemple, celui dont ce document est un rendu HTML, et qui se trouve à l'adresse <https://github.com/vlc1/ene-4102c-td/blob/master/td-1.jl>). Il vous suffit simplement de copier-coller cette adresse dans le champs *Open from file* et de cliquer sur *Open* !
"""

# ╔═╡ Cell order:
# ╟─9df1a0de-f33e-11ea-276a-5175cf87a652
# ╟─8e0646ee-f347-11ea-22f6-efa2515444b3
# ╠═dfb0b044-f347-11ea-1dee-2d161f8aa0c8
# ╠═e8c0f304-f347-11ea-20c2-7f3c5b9327f7
# ╠═bfe43eb2-f347-11ea-257a-e5ebe3e91eb9
# ╠═ddea46f4-f347-11ea-3714-b5212290ca66
# ╟─0bb68dfc-f347-11ea-1d1b-37be8fcd1cff
