### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ c0db5a3e-f39c-11ea-12fa-3b2f78b7b2b7
md"""
Versions [Pluto](https://github.com/vlc1/ene-4102c-td/blob/master/td1/part3.jl) et [Jupyter](https://vlc1.github.io/ene-4102c/td1/part3.ipynb) de ce notebook.

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

# ╔═╡ f54ffbe4-07ac-11eb-1ac9-7d4528dc1dd8
md"""
## Correction

"""

# ╔═╡ fb4c298e-07ac-11eb-2982-4d3322766983
# Question 1
begin
	U = collect(0:50)
	U, length(U)
end

# ╔═╡ 08cc94cc-07ad-11eb-38d8-3bb9b7bd0c67
# Question 2
begin
	V = U[begin:begin + 4]
	W = [U[begin:begin + 4]; U[end - 4:end]]
	V, W
end

# ╔═╡ 12209582-07ad-11eb-2d27-ad6a96675c39
# Question 3
M = transpose(reshape(1:30, 10, 3))

# ╔═╡ 60bc41aa-07ad-11eb-32fc-e5087f1353ae
# Question 4
begin
	N = M[:, begin:begin + 1]
	P = M[:, end - 2:end]
	Q = M[begin:2:begin + 2, begin + 2:4:begin + 6]
	N, P, Q
end

# ╔═╡ e9822b76-07ad-11eb-1f95-b5e39f3a7615
# Question 5
R = M[:, begin:2:end]

# ╔═╡ 1c1e59ec-07ae-11eb-009f-4d280ccbf6dd
# Question 6
begin
	A = collect(-1:-2:-99)
	B = collect(2:2:100)
	C = collect(Iterators.flatten(zip(A, B)))
	A, B, C
end

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
# ╟─f54ffbe4-07ac-11eb-1ac9-7d4528dc1dd8
# ╠═fb4c298e-07ac-11eb-2982-4d3322766983
# ╠═08cc94cc-07ad-11eb-38d8-3bb9b7bd0c67
# ╠═12209582-07ad-11eb-2d27-ad6a96675c39
# ╠═60bc41aa-07ad-11eb-32fc-e5087f1353ae
# ╠═e9822b76-07ad-11eb-1f95-b5e39f3a7615
# ╠═1c1e59ec-07ae-11eb-009f-4d280ccbf6dd
# ╟─ca8cb49a-f39c-11ea-307e-97825e864b0f
# ╟─46497466-f39e-11ea-2634-cd5985460790
# ╟─7fb06e78-f39f-11ea-38d9-ab51b94af7c5
