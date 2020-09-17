### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ c0db5a3e-f39c-11ea-12fa-3b2f78b7b2b7
md"""
Ce notebook est à retrouver [ici](https://github.com/vlc1/ene-4102c-td/blob/master/td1-part3.jl).

"""

# ╔═╡ 03cb7526-f392-11ea-123a-797739ba60de
md"""
# Exercice 1 -- Manipulation de matrices

1. Définir le vecteur $u = \left [ \begin{array}{cccccc} 0 & 1 & 2 & 3 & \cdots & 49 & 50 \end{array} \right ]$. Quelle est sa taille ?
1. Définir le vecteur $v$ contenant les cinq premiers éléments de $u$, et le vecteur $w$ contenant les cinq premiers et les cinq derniers éléments de $u$.
1. Définir la matrice
$$M =
\left [ \begin{array}{ccccccc}
1 & 2 & 3 & \cdots & 8 & 9 & 10 \\
11 & 12 & 13 & \cdots & 18 & 19 & 20 \\
21 & 22 & 23 & \cdots & 28 & 29 & 30
\end{array}
\right ].$$
4. Extraire de $M$ les matrices
$$N =
\left [ \begin{array}{cc}
1 & 2 \\
11 & 12 \\
21 & 22
\end{array} \right ], \quad P =
\left [ \begin{array}{ccc}
8 & 9 & 10 \\
18 & 19 & 20 \\
28 & 29 & 30
\end{array} \right ] \quad \mathrm{et} \quad Q =
\left [ \begin{array}{cc}
3 & 7 \\
23 & 27 \end{array} \right ].$$
5. Extraire de la matrice $M$ la matrice $R$ obtenue en prenant une colonne sur deux.
1. Définir les matrices $A = \left [ \begin{array}{ccccc} 2 & 4 & 6 & 8 & \cdots & 100 \end{array} \right ]$ et $B = \left [ \begin{array}{ccccc} -1 & -3 & -5 & \cdots & -99 \end{array} \right ]$ puis le vecteur $x = \left [ \begin{array}{ccccccccc} -1 & 2 & -3 & 4 & -5 & 6 & \cdots & -99 & 100 \end{array} \right ]$.

"""

# ╔═╡ ca8cb49a-f39c-11ea-307e-97825e864b0f
md"""
# Exercice 2 -- Matrices et systèmes linéaires

1. Écrire une fonction, n'utilisant aucune boucle (`for`, `while`...) qui prend comme paramètre un entier $n$ et qui construit la matrice suivante :
$$\left [ \begin{array}{ccccccc}
1 & 1 & 0 & \cdots & 0 & 0 & 0 \\
\frac{1}{n} & 2 & \frac{n-1}{n} & \cdots & 0 & 0 & 0 \\
0 & \frac{2}{n} & 3 & \cdots & 0 & 0 & 0 \\
\vdots & \vdots & \vdots & \ddots & \vdots & \vdots & \vdots \\
0 & 0 & 0 & \cdots & n - 1 & \frac{2}{n} & 0 \\
0 & 0 & 0 & \cdots & \frac{n - 1}{n} & n & \frac{1}{n} \\
0 & 0 & 0 & \cdots & 0 & 1 & n + 1
\end{array} \right ].$$
"""

# ╔═╡ 46497466-f39e-11ea-2634-cd5985460790
Markdown.MD(Markdown.Admonition("hint", "Indice", [md"Renseignez-vous sur le type `Tridiagonal` de la bibliothèque standart `LinearAlgebra` grâce au *Live docs*."]))

# ╔═╡ 7fb06e78-f39f-11ea-38d9-ab51b94af7c5
md"""
2. Résoudre numériquement le système
$$\left \{ \begin{aligned}
x + 2y + 3z + 4t & = 1, \\
2x + 3y + 4z + t & = -2, \\
-2x + 4y -5z + 2t & = 0, \\
8x + y - z + 3t & = 1.
\end{aligned} \right .$$

"""

# ╔═╡ Cell order:
# ╟─c0db5a3e-f39c-11ea-12fa-3b2f78b7b2b7
# ╟─03cb7526-f392-11ea-123a-797739ba60de
# ╟─ca8cb49a-f39c-11ea-307e-97825e864b0f
# ╟─46497466-f39e-11ea-2634-cd5985460790
# ╟─7fb06e78-f39f-11ea-38d9-ab51b94af7c5
