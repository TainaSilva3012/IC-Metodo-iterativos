### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 084d9850-281d-46ad-9ec3-5c391dd11a40
begin
using LinearAlgebra, Random, BenchmarkTools, PlutoUI, SparseArrays, DataFrames
    Random.seed!(3)
    n=111
    r=200
    A₀ = randn(n,r)
    b₀ = randn(n)
  #Primeiro precisamos fazer uma função que normalize as linhas da matriz A
end 

# ╔═╡ 32e51ae3-b6e1-4a3b-ad86-83c809005756
md"""
## Ferramentas
"""

# ╔═╡ e80dedf1-b93b-451f-9b28-d19230c25397
md"""
**Lema 1:**  Um espaço afim $v + \mathcal{M}$, em que $dim(\mathcal{M}) = m-1$ é dito um hiperplano. A $i$-ésima equação de um sistema linear $Ax=b$, $A \in \mathbb{R}^{m \times n}$, isto é, $A_{i*}x = b_i$ é um hiperplano em $\mathbb{R}^{m}$. Desta forma, a solução do sistema se encontra na intersecção dos $m$ hiperplanos definidos pelas linhas de $A$. É possível fazer as seguintes afirmações:

i) Dado um escalar $\beta$ e um vetor não-nulo $u$, o conjunto $\mathcal{H}= \{ x \mid u^{T}x= \beta\}$ é um hiperplano em $\mathbb{R}^{m}$.

ii) A projeção ortogonal de $b \in \mathbb{R}^{m}$ sobre $\mathcal{H}$ é $p=b-(\frac{u^{T}b - \beta}{u^{T}u})u$.

**Demonstração:** Os pontos em $\mathcal{H}$ são apenas soluções para um sistema linear $u^{T}x=\beta$. Usando o fato de que a solução geral de um sistema linear é uma solução particular mais a solução geral de qualquer sistema linear é a uma solução particular mais a solução geral da equação homogênea associada, assim, obtemos

$\mathcal{H} = \frac{\beta u}{u^{T} u} + N(u^{T})= \frac{\beta u}{u^{T} u} + [R(u)]^{\perp} = \frac{\beta u}{u^{T} u} + u^{\perp},$

em que $u^{\perp}$ denota o complemeto ortogonal do espaço unidimensional abrangido pelo vetor $u$. Assim $\mathcal{H}= v + \mathcal{M}$, em que $v = \frac{\beta u}{u}$ e $\mathcal{M} = u^{\perp}$, em que pelo teorema do núcleo e da imagem, $dim(\mathcal{M})=dim(u^{\perp})=n-1$, visto que $u$ é unidimensional.
"""

# ╔═╡ 28652520-d768-11eb-1af1-b5b2542039de
md"""
## Método das projeções de Kaczmarz
"""

# ╔═╡ 002d6935-5c8b-4e64-aa73-f9e5c31ef1c2
md"""
 O método foi criado pelo polonês Stefan Kaczmarz no ano de 1937. O algoritmo utiliza, a cada passo, uma equação do sistema $Ax=b$, de forma ordenada.

"O método de Kaczmarz clássico, também conhecido como Técnica de Reconstrução Algébrica (ART - Algebraic Reconstruction Technique), é um algoritmo iterativo para solucionar sistemas do tipo $Ax = b$ que age realizando operações nas linhas da matriz $A$" (ESTACIO, 2014, p.34).

"Utilizando uma equação da matriz $A$ por vez e atualizando a solução aproximada corrente $(x_{k})$ a cada iteração o método converge para a solução deste sistema linear caso ele seja consistente. E, caso o sistema seja inconsistente, o algoritmo converge até uma região próxima a interseção entre os hiperplanos e, não necessariamente, a uma solução de mínimos quadrados" (ESTACIO, 2014, p.34).

A solução do sitema linear não-singular

$\begin{bmatrix}
a_{11} & a_{12}\\
a_{21} & a_{22}  
\end{bmatrix}
\cdot 
\begin{bmatrix}
x_{1}\\
x_{2}
\end{bmatrix}
=
\begin{bmatrix}
b_{1}\\
b_{2}
\end{bmatrix}$

é a intersecção de dois hiperplanos (neste caso duas retas) definidos por

$\mathcal{H}_{1} = \{(x_{1}, x_{2}) \mid a_{11} \cdot x_{1} + a_{12} \cdot x_{2} = b_{1}\}$

$\text{ e }$ 

$\mathcal{H}_{2} = \{(x_{1}, x_{2}) \mid a_{21} \cdot x_{1} + a_{22} \cdot x_{2} = b_{2}\}.$

É visualmente evidente que, começando com um ponto arbitrário $P_{0}$ e realizando uma squência de projeções ortogonais alternadas em $\mathcal{H}_{1}$ e $\mathcal{H}_{2}$ como representado na figura abaixo, a sequência resultante de projeções $\{p_{1}, p_{2}, p_{3}, p_{4}, \ldots \}$ converge para $\mathcal{H}_{1} \cap \mathcal{H}_{2}$, que é a solução de $Ax=b$.

"""

# ╔═╡ 6562f2a7-9e3a-4c89-91aa-3d26e1d7f5d4
md"""

![Imagem](https://i.ibb.co/BBRs4cs/kaczmarzs.jpg)

"""

# ╔═╡ 63954ace-60fd-4938-b4e5-17a80051d427
md"""
Essa idéia pode ser generalizada usando o Lema 1.
Para um sistema consistente $A_{n \times r}x=b$ em que $rank(A)=r$, escale as linhas de modo que $\| A_{i*}\|_{2} = 1$ para todo $i$, em que $\mathcal{H}_{i}=\{ x \mid A_{i*}x = b_{i}\}$ é o hiperplano definido pela i-ésima equação.
Comece com um vetor arbitrário $p_{0} \in \mathbb{R}^{r \times 1}$, e sucessivamente realizar projeções ortogonais em cada hiperplano para gerar a seguinte sequência:

$\begin{align}
p_{1} =& p_{0}-(A_{1*}p_{0}-b_{1})(A_{1*})^{T}\\
p_{2} =& p_{1}-(A_{2*}p_{1}-b_{2})(A_{2*})^{T}\\
&\vdots\\
p_{n} =& p_{n-1}-(A_{n*}p_{n-1}-b_{n})(A_{n*})^{T}
\end{align}$

Quando todos os $n$ hiperplanos tiverem sido usados, continue repetindo o processo. Por exemplo, na segunda passagem projete $p_{n}$ em $\mathcal{H}_{1}$; então projete $p_{n+1}$ em $\mathcal{H}_{2}$, etc. Para um $p_{0}$ arbitrário, toda a sequência de Kaczmarz é gerada executando o seguinte loop duplo:
"""

# ╔═╡ 59729285-eca1-4ffb-95f9-f10f06053e4e
md"""
For $k = 0, 1, 2, 3, \ldots$

For $i = 1, 2, \ldots, n$

$p_{kn+i}=p_{kn+i-1}-(A_{i*}p_{kn+i-1}-b_{i})(A_{i*})^{T}$
"""

# ╔═╡ f80e4298-a076-4139-aca7-c49fc257bf8b
md"""
Provaremos que a sequência de Kaczmarz converge para a solução de $Ax = b$ mostrando

$\|p_{kn+i}-x\|_{2}^{2}=\|p_{kn+i-1}-x\|_{2}^{2}-(A_{i*}p_{kn+i-1}-b_{i})^{2}$ 
"""

# ╔═╡ acbde7f4-6bd9-4aef-b1c5-d8c1269f01f4
md"""
**Demonstação:** Por conveniência defina $\beta=A_{i*}p_{kn+i-1}-b_{i}$ de forma que, de acordo com a sequência de passos do algoritmo temos $p_{kn+i}=p_{kn+i-1}-\beta(A_{i*})^{T}$. 
Também utilizaremos o seguinte fato $A_{i*}(p_{kn+i-1}-x)=A_{i*}p_{kn+i-1}-b_{i}=\beta$, junto com $\| A_{i*} \|^{2}_{2} = 1$ de forma a escrever

$\begin{align}
\|p_{kn+i}-x \|^{2}_{2} &= \| p_{kn+i-1} - \beta(A_{i*})^{T} - x \|^{2}_{2} \\
&=\| (p_{kn+i-1} - x) - \beta(A_{i*})^{T} \|^{2}_{2} \\
&= (p_{kn+i-1} - x)^{T} (p_{kn+i-1} - x) - 2 \beta A_{i*} (p_{kn+i-1} - x) + \beta^{2}A_{i*}(A_{i*})^{T}\\
&= \|p_{kn+i-1}-x \|^{2}_{2} - \beta^{2}.
\end{align}$

Consequentemente, observamos que $\|p_{kn+i}-x \|_{2} \leq \|p_{kn+i-1}-x \|_{2}$. A igualdade vale se, e somente se, $\beta = 0$, ou de maneira equivalente, se, e somente se, $p_{kn+i-1} \in \mathcal{H}_{i-1} \cap \mathcal{H}_{i}$. Conclui-se que a sequência de normas $\|p_{kn+i}-x \|_{2}$ é monótona decrescente, e desta forma possui um valor limitante. Isto implica que sequência de $\beta$'s definida anteriormente deve tender a 0, e desta forma a sequência de $p_{kn+i}$'s converge para $x$.
"""

# ╔═╡ 040be751-a13b-4e08-99d5-50fc8a3bbeea
md"""
**Implementação Computacional**
"""

# ╔═╡ cb71b9b2-6d8c-4dc8-9016-e3be6ebb465a
md"""
Vamos agora programar o método de Kaczmarz.
"""

# ╔═╡ 28d7ddef-0276-48ed-9c88-2de9c3eb5e15
  function normalize_matriz(A₀,b₀)
      n,r = size(A₀)
      for i=1:n
           β = norm(A₀[i,:])
           A₀[i,:] = A₀[i,:]/β
           b₀[i] = b₀[i]/β
       end
       return A₀,b₀
    end

# ╔═╡ 387cb3ba-34d1-4c19-9c13-5843297a0eb9
A,b = normalize_matriz(A₀,b₀) #Teste da função que normaliza as linhas de A e recalcula o vetor b

# ╔═╡ 3b3c923d-5e9f-489a-9ae7-e9b3041370d6
# Então, precisamos fazer uma função de projeção
function proj(p₀, u, β)
	return  p₀ - ((dot(u, p₀)- β)/dot(u, u))*u
end

# ╔═╡ 8a050838-4131-48bd-a143-e9d5977c8c84
#Por fim fazemos o algoritmo de kaczmars utilizando a função de normalização e a de projeção
function Kaczmarz(A::Matrix, b::Vector; itmax::Int=10000, ε::Float64= 1e-3)
    n,r = size(A) 
    pₖ = ones(r)
    k = 1 
    A,b = normalize_matriz(A,b)
    tol = norm(b- A*pₖ)
    while k <= itmax && tol> ε
      for i=1:n
          pₖ= proj(pₖ, A[i,:], b[i])
      end
        tol = norm(b- A*pₖ)
        k += 1
    end
    return pₖ, tol, k
end

# ╔═╡ 476969f9-27c4-49c1-84a2-d38e81c1437f
Kaczmarz(A, b, ε=1e-6) #Teste da função que implementa o método iterativo de Kaczmarz

# ╔═╡ 55aa1871-59c9-40c8-a556-ef9bab302043
md"""
####  Método de reflexão de Cimmino
"""

# ╔═╡ 98aa94e9-7922-4450-84ba-ae644f0780ff
md"""
O método de Cimmino é um  método iterativo para solução de sistemas lineares, foi publicado em 1938 pelo matemático italiano Gianfranco Cimmino (1908-1989), estudante de Mauro Picone. O algoritmo usa simultaneamente todas as equações do sistema linear a cada passo iterativo.

Segundo Benzi (2004, p.6), o método data de pelo menos 1932 e é possível que Cimmino nunca se preocupou em publicá-lo se não pela insistência de Picone. É provável que Cimmino estava ocupado com sua pesquisa [...] e relutante a dedicar seu tempo para escrever sobre um tópico considerado de menor importância.

Os métodos de Cimmino e Kaczmarz são proximamente relacionados. O algoritmo de Cimmino se provou ser melhor formatado para computadores paralelos, enquanto o método de Kaczmarz tende a convergir mais rapidamente" (BENZI, 2004, p.8).

O matemático italiano Gianfranco Cimmino usou a seguinte observação elementar para construir um algoritmo iterativo para resolver sistemas lineares. Para um sistema $2 \times 2$ tem $Ax = b$, sejam $\mathcal{H}_{1}$ e $\mathcal{H}_{2}$ as duas retas (hiperplanos) definidas pelas duas equações. Para um palpite arbitrário $r_{0}$, seja $r_{1}$ o reflexão de $r_{0}$ sobre a linha $\mathcal{H}_{1}$, e seja $r_{2}$ o reflexão de $r_{0}$ sobre a linha $\mathcal{H}_{2}$. Conforme ilustrado na abaixo, os três pontos $r_{0}$, $r_{1}$ e $r_{2}$ encontram-se em um círculo cujo centro é $\mathcal{H}_{1} \cap \mathcal{H}_{2}$ (a solução do sistema).

"""

# ╔═╡ dc992285-74e8-4da1-83d0-a143ee146b90
md"""
![Imagem](https://i.ibb.co/CMrVJfc/cimmino.jpg)
"""

# ╔═╡ d3bb2863-e421-4dae-8819-5a69f6016730
md"""
O valor médio $m = \frac{(r_{1} + r_{2})}{2}$ está estritamente dentro do círculo, então $m$ é uma melhor aproximação da solução do que $r_{0}$. É visualmente evidente que iteração produz uma sequência que converge para a solução de $Ax = b$. Provaremos isso em geral usando o esquema a seguir.

a) Para um escalar $\beta$ e um vetor $u \in \mathbb{R}^{n}$ tal que $\|u\|_{2} = 1$, considere o hiperplano $\mathcal{H}= \{x \mid u^{T}x = \beta\}$. Usaremos a **Definição 2** para mostrar que a reflexão de um vetor $b$ sobre $\mathcal{H}$ é $r=b-2(u^{T}b - \beta)u$.

b) Para um sistema $Ax = b$ em que as linhas de $A \in R^{n \times r}$ foram escalonadas de forma que $\|A_{i∗}\|^{2} = 1$ para todo $i$, seja $\mathcal{H}_{i}=\{x \mid A_{i∗}x=b_{i}\}$ o hiperplano definido pela $i$-ésima equação. Se $r_{0} \in \mathbb{R}^{r \times 1}$ é um vetor arbitrário, e se $r_{i}$ é a reflexão de $r_{0}$ sobre $\mathcal{H}_{i}$, explicaremos por que o valor médio das reflexões $\{r_{1}, r_{2}, \ldots, r_{n}\}$ é $m=r_{0} - (\frac{2}{n})A^{T}\varepsilon$, em que $\varepsilon = Ar_{0}-b$.

c) A iteração da parte b produz $m_{k}=m_{k−1}-\frac{2}{n}A^{T}\varepsilon_{k − 1}$, em que $\varepsilon_{k−1} = Am_{k−1}-b$. Mostraremos que se $A$ é não-singular, e se $x=A^{−1}b$, então $x-m_{k}=(I-\frac{2}{n}A^{T}A)^{k}(x-m_{0})$. 
**Nota:** Pode ser provado que se $(I-\frac{2}{n} A^{T}A)^{k} \rightarrow 0$ quando $k \rightarrow \infty$, então $m_{k} \rightarrow x$ para qualquer $m_{0}$ escolhido. Na verdade, $m_{k}$ converge mesmo se o  posto $A$ for deficiente - se consistente, ele converge para uma solução e, se inconsistente, o limite é uma solução de mínimos quadrados. Método de Cimmino também funciona com ponderado significa. Se $W = diag (w_{1}, w_{2}, \ldots, w_{n})$, em que $w_{i}> 0$ e $w_{i}=1$, então $m_{k}=m_{k−1}-\omega A^{T}W \varepsilon_{k−1}$ é uma sequência convergente em que $0 < \omega < 2$ é um “relaxamento parâmetro” que pode ser ajustado para alterar a taxa de convergência. 


**Demonstração:**

a) Como já definido anteriormente temos $\mathcal{H} = v + u^{\perp}$, em que $v = \beta u$. Então, subtraindo $v$ de todos os termos no hiperplano $\mathcal{H}$ e também da solução $b$ transladamos o hiperplano para a origem. Como mostrado na figura abaixo, este processo move $\mathcal{H}$ para $u^{\perp}$ e translada $b$ para $b-v$ e $r$ para $r-v$.
Utilizando reflexões de Householder, sabemos que a reflexão de $b-v$ sobre $u^{\perp}$ é $r -v = R(b-v)=(I -2uu^{T})(b-v)= b+ (\beta - 2u^{T}b)u$.
Portanto, a reflexão de $b$ sobre $\mathcal{H}$ é $r = R(b-v)+v=b-2(u^{T}b- \beta )u$.

b) Da parte anterior, conclui-se que a reflexão de $r_{0}$ sobre $\mathcal{H}_{i}$ é $r_{i} = r_{0} - 2(A_{i*}r_{0}-b_{i})(A_{i*}^{T})$, e portanto, o valor médio de todas as reflexões $\{r_{1}, r_{2}, \ldots, r_{n}\}$ é

$\begin{align}
m &= \frac{1}{n}\sum_{i=1}^{n}r_{i} = \frac{1}{n}\sum_{i=1}^{n}(r_{0}-2(A_{i*}r_{0}-b_{i})(A_{i*}^{T}))\\
&= r_{0} - \frac{2}{n}\sum_{i=1}^{n}(A_{i*}r_{0}-b_{i})(A_{i*})^{T}\\
&= r_{0} - \frac{2}{n} A^{T}(Ar_{0} - b) = r_{0} - \frac{2}{n}A^{T} \varepsilon.
\end{align}$

c) Primeiro note que

$\begin{align}
x - m_{k} &= x - m_{k-1} + \frac{2}{n}A^{T} \varepsilon_{k-1}\\
&= x - m_{k-1} + \frac{2}{n}A^{T}(Am_{k-1} - b)\\
&= x - m_{k-1} + \frac{2}{n}A^{T}(Am_{k-1} - Ax)\\
&= x - m_{k-1} + \frac{2}{n}A^{T}A(m_{k-1} - x)\\
&= (I - \frac{2}{n}A^{T}A)(x - m_{k-1})
\end{align}$

e então usando substituições sucessivas se conclui que $x - m_{k} = (I - \frac{2}{n}A^{T}A)^{k}(x-m_{0})$.
"""

# ╔═╡ d2a457d1-228f-408d-8c54-fe594610737c
md"""
**Implementação Computacional**
"""

# ╔═╡ ea96100b-14b9-4ebf-aa7d-12d8d4563f38
# Então, precisamos fazer uma função de reflexão
function reflexao(p₀, u, β)
	return  2*proj(p₀, u, β) - p₀
end
   

# ╔═╡ 1d054b6c-20e1-42f9-a174-7c20b95353a3
function Cimmino(A::Matrix, b::Vector; itmax::Int=10000, ε::Float64= 1e-3)
	n,r = size(A)
	xₖ = ones(r)
	k = 1
	A,b = normalize_matriz(A,b)
	tol = norm(b- A*xₖ)
 		while k <= itmax && tol> ε 
			rₖ = zeros(r)
			for i=1:n
				rₖ += reflexao(xₖ, A[i,:], b[i])
			end
			xₖ= rₖ/n
			tol = norm(b- A*xₖ)
			k += 1
		end
	return xₖ, tol, k
end

# ╔═╡ 3e4d9b19-c458-429a-9660-ab9c4275e964
Cimmino(A, b, ε=1e-3)

# ╔═╡ 45ea236e-3371-44f4-ba45-118cd95aea51
Kaczmarz(A, b, ε=1e-3)

# ╔═╡ 3cf130fe-0a03-4c94-b869-665c84ae12c2
md"""
**Referências**

BENZI, Michele. **Gianfranco Cimmino’s contributions to Numerical Mathematics.** Atlanta, Georgia, 2005.

CERQUEIRA, Ana Cecília Sanches. **Um estudo sobre sequências e séries.** 2013. Diss. (Mestrado) – Instituto de Geociências e Ciências Exatas, Universidade Estadual Paulista Júlio de Mesquita Filho, Rio Claro.

ESTÁCIO, Leonardo Bravo. **Sobre a escolha da relaxação e ordenação dasprojeções no método de Kaczmarz com ênfase em implementações altamente paralelas e aplicações em reconstrução tomográfica.** 2014. Diss. (Mestrado) –Instituto de Ciências Matemáticas e de Computação, Universidade de São Paulo, São Carlos.

MEYER, Carl D. **Matrix analysis and Applied linear algebra.** Philadelphia:SIAM, 2000. p. 869.
"""

# ╔═╡ 8cef3582-ba2f-46bb-90f7-dd4be4fa11c9
md"""
## CAV
"""

# ╔═╡ ad6bf2b9-7851-4afb-9478-475ef1cf6e22
md"""
Cav é basicamente Cimmino mas com a uma normalização diferente
"""

# ╔═╡ 55633e53-5e31-43a2-a60a-d297c2334ddf
  function norma_S_matriz(A₀, b₀; S::Vector=[]) #Caso não seja fornecida a E, cacula norma Euclidiana
    n,r = size(A₀)
	if isempty(S) 
		S = ones(r)
	end   
	for i=1:n
		βᵢ = A₀[i,:]'*S[i]*A₀[i,:]
        A₀[i,:] /=βᵢ
        b₀[i] /= βᵢ
       end
       return A₀, b₀
    end

# ╔═╡ c2880646-720f-4d31-9f4a-15b55a897305
let
	L = [ 3. 4;
	  	8  7]
	R = [15.; 10]
	K =[1/2.,  1/3]
	norma_S_matriz(L,R,S=K)
end

# ╔═╡ 8206aa80-e5f7-4361-8e4c-a5be2deb8385
#Ver como só permitir E de size rxr para não comprometer as operações

# ╔═╡ 4bc789a8-f4db-4050-8e76-828163c0a2ed
spa = sprand(100, 100, 0.05)

# ╔═╡ 472b35f1-f277-4d8c-988c-e81aa3726568
vet_spa = rand(100)

# ╔═╡ b4cb7e0e-6b7c-485b-87ad-220da542898d
function Esparsidade(A::SparseMatrixCSC)
	S = Int[]
	for col in eachcol(A)
		push!(S, nnz(col))
	end
	return S
end

# ╔═╡ 750ce1cf-e6b4-45af-868a-71fabf479421
S = Esparsidade(spa)

# ╔═╡ a16eae89-22a8-4728-9b81-ab4265ffbb51
function CAV_sparse(A::SparseMatrixCSC, b::Vector; itmax::Int=10000, ε::Float64= 1e-3)
	n,r = size(A)
	xₖ = ones(r)
	k = 0
	E = Esparsidade(A)
	A,b = norma_S_matriz(A₀, b₀; S=E) 
	tol = norm(b- A*xₖ)
 		while k <= itmax && tol> ε 
			rₖ = zeros(r)
			for i=1:n
				rₖ += reflexao(xₖ, A[i,:], b[i])
			end
			xₖ= rₖ/n
			tol = norm(b- A*xₖ)
			k += 1
		end
	return xₖ, tol, k
end

# ╔═╡ f6b32c1c-17ed-4238-ba46-01d3a21c5319
CAV_sparse(spa, vet_spa)

# ╔═╡ 95e67977-9703-4eba-bede-f71bf7576e10
begin
#df = DataFrame(Matriz=Matrix[], Vetor=Vector[], Solução=Float64[], Tolerância=Int64[], Iterações=Int64[]);
#push!(df, [xₖ, tol, k])
end

# ╔═╡ 2d28fae0-5c11-46d3-bdde-fb43dcd99bf0
function CAV(A::Matrix, b::Vector; itmax::Int=10000, ε::Float64= 1e-3)
	if issparse(A) == true 
		CAV_sparse(A, b; itmax=10000, ε=1e-3)
	else 
		Cimmino(A, b; itmax=10000, ε=1e-3)
	end
	return xₖ, tol, k
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[compat]
BenchmarkTools = "~1.1.0"
DataFrames = "~1.2.0"
PlutoUI = "~0.7.9"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Statistics", "UUIDs"]
git-tree-sha1 = "ffabdf5297c9038973a0a3724132aa269f38c448"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.1.0"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "dc7dedc2c2aa9faf59a55c622760a25cbefbe941"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.31.0"

[[Crayons]]
git-tree-sha1 = "3f71217b538d7aaee0b69ab47d9b7724ca8afa0d"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.0.4"

[[DataAPI]]
git-tree-sha1 = "ee400abb2298bd13bfc3df1c412ed228061a2385"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.7.0"

[[DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "1dadfca11c0e08e03ab15b63aaeda55266754bad"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.2.0"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "4437b64df1e0adccc3e5d1adbc3ac741095e4677"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.9"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[InvertedIndices]]
deps = ["Test"]
git-tree-sha1 = "15732c475062348b0165684ffe28e85ea8396afc"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.0.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "81690084b6198a2e1da36fcfda16eeca9f9f24e4"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.1"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "4ea90bd5d3985ae1f9a908bd4500ae88921c5ce7"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.0"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "c8abc88faa3f7a3950832ac5d6e690881590d6dc"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "1.1.0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlutoUI]]
deps = ["Base64", "Dates", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "Suppressor"]
git-tree-sha1 = "44e225d5837e2a2345e69a1d1e01ac2443ff9fcb"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.9"

[[PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "cde4ce9d6f33219465b55162811d8de8139c0414"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.2.1"

[[PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "0d1245a357cc61c8cd61934c07447aa569ff22e6"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.1.0"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Reexport]]
git-tree-sha1 = "5f6c21241f0f655da3952fd60aa18477cf96c220"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.1.0"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[Suppressor]]
git-tree-sha1 = "a819d77f31f83e5792a76081eee1ea6342ab8787"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.0"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "8ed4a3ea724dac32670b062be3ef1c1de6773ae8"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.4.4"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─32e51ae3-b6e1-4a3b-ad86-83c809005756
# ╟─e80dedf1-b93b-451f-9b28-d19230c25397
# ╟─28652520-d768-11eb-1af1-b5b2542039de
# ╟─002d6935-5c8b-4e64-aa73-f9e5c31ef1c2
# ╟─6562f2a7-9e3a-4c89-91aa-3d26e1d7f5d4
# ╟─63954ace-60fd-4938-b4e5-17a80051d427
# ╟─59729285-eca1-4ffb-95f9-f10f06053e4e
# ╟─f80e4298-a076-4139-aca7-c49fc257bf8b
# ╟─acbde7f4-6bd9-4aef-b1c5-d8c1269f01f4
# ╟─040be751-a13b-4e08-99d5-50fc8a3bbeea
# ╟─cb71b9b2-6d8c-4dc8-9016-e3be6ebb465a
# ╠═084d9850-281d-46ad-9ec3-5c391dd11a40
# ╠═28d7ddef-0276-48ed-9c88-2de9c3eb5e15
# ╠═387cb3ba-34d1-4c19-9c13-5843297a0eb9
# ╠═3b3c923d-5e9f-489a-9ae7-e9b3041370d6
# ╠═8a050838-4131-48bd-a143-e9d5977c8c84
# ╠═476969f9-27c4-49c1-84a2-d38e81c1437f
# ╟─55aa1871-59c9-40c8-a556-ef9bab302043
# ╟─98aa94e9-7922-4450-84ba-ae644f0780ff
# ╟─dc992285-74e8-4da1-83d0-a143ee146b90
# ╟─d3bb2863-e421-4dae-8819-5a69f6016730
# ╟─d2a457d1-228f-408d-8c54-fe594610737c
# ╠═ea96100b-14b9-4ebf-aa7d-12d8d4563f38
# ╠═1d054b6c-20e1-42f9-a174-7c20b95353a3
# ╠═3e4d9b19-c458-429a-9660-ab9c4275e964
# ╠═45ea236e-3371-44f4-ba45-118cd95aea51
# ╟─3cf130fe-0a03-4c94-b869-665c84ae12c2
# ╟─8cef3582-ba2f-46bb-90f7-dd4be4fa11c9
# ╟─ad6bf2b9-7851-4afb-9478-475ef1cf6e22
# ╠═55633e53-5e31-43a2-a60a-d297c2334ddf
# ╠═c2880646-720f-4d31-9f4a-15b55a897305
# ╠═8206aa80-e5f7-4361-8e4c-a5be2deb8385
# ╠═4bc789a8-f4db-4050-8e76-828163c0a2ed
# ╠═472b35f1-f277-4d8c-988c-e81aa3726568
# ╠═b4cb7e0e-6b7c-485b-87ad-220da542898d
# ╠═750ce1cf-e6b4-45af-868a-71fabf479421
# ╠═a16eae89-22a8-4728-9b81-ab4265ffbb51
# ╠═f6b32c1c-17ed-4238-ba46-01d3a21c5319
# ╠═95e67977-9703-4eba-bede-f71bf7576e10
# ╠═2d28fae0-5c11-46d3-bdde-fb43dcd99bf0
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
