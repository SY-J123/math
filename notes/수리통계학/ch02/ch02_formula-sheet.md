# Chapter 2. Random Variables

## 2.1 Definition of a Random Variable

**Definition 2.1.1 (Def 2.1) — random variable**

A random variable is a real-valued function $X:S\to\mathbb{R}$ defined on the sample space, assigning the real number $X(s)$ to each outcome $s\in S$.

**Example 2.1.1 — three coin tosses, number of heads**

$X=$ number of heads: $X(HHH)=3$; $X(HHT)=X(HTH)=X(THH)=2$; $X(HTT)=X(THT)=X(TTH)=1$; $X(TTT)=0$. With $P(\text{head})=p$ and independent tosses,

$$P(X=0)=(1-p)^3,\quad P(X=1)=3p(1-p)^2,\quad P(X=2)=3p^2(1-p),\quad P(X=3)=p^3$$

**Definition 2.1.2 — discrete and continuous random variables**

- **discrete**: the possible values are finite or countably infinite (e.g. $\{0,1,2,\dots\}$)
- **continuous**: the possible values fill an interval (uncountable), e.g. a lifetime or a length

## 2.2 Probability Density Function and Cumulative Distribution Function

**Definition 2.2.1 — probability density function (pmf, pdf)**

$f(x)$ describes the probability that $X$ takes each value, and must satisfy

1. $f(x_i)\ge0$ (discrete), $f(x)\ge0$ (continuous)
2. total mass $1$: $\ \sum_i f(x_i)=1$ (discrete), $\int_{-\infty}^{\infty}f(x)dx=1$ (continuous)

In the continuous case $P(X=a)=0$, and probability is measured as area:

$$P(X\in A)=\begin{cases}\sum_{x_i\in A}f(x_i)&\text{discrete}\\ \int_A f(x)dx&\text{continuous}\end{cases}\qquad P(a\le X\le b)=\int_a^b f(x)dx \tag{2.1}$$

**Example 2.2.1 (Ex 2.1) — checking a density**

$f(x)=\dfrac{x^2}{9}$ on $(0,3)$ is a pdf, since $f\ge0$ and

$$\int_0^3\frac{x^2}{9}dx=\Big[\frac{x^3}{27}\Big]_0^3=1,\qquad P(0<X<2)=\int_0^2\frac{x^2}{9}dx=\frac{8}{27}$$

**Definition 2.2.2 (Def 2.2) — cumulative distribution function**

$$F(x)=P(X\le x),\qquad -\infty<x<\infty$$

Written $X\sim f(x)$ or $X\sim F(x)$.

**Theorem 2.2.1 (Thm 2.1) — characterization of a cdf**

$F$ is the cdf of some random variable if and only if

1. $F$ is non-decreasing: $a<b\Rightarrow F(a)\le F(b)$
2. $\lim_{x\to-\infty}F(x)=0$ and $\lim_{x\to\infty}F(x)=1$
3. $F$ is right-continuous: $\lim_{h\to0^+}F(x+h)=F(x)$

From the density,

$$F(x)=\sum_{x_i\le x}f(x_i)\ \ \text{(discrete)} \tag{2.2}$$

$$F(x)=\int_{-\infty}^{x}f(t)dt\ \ \text{(continuous)} \tag{2.3}$$

**Example 2.2.2 (Ex 2.2, 2.3) — obtaining a cdf**

- three tosses with $p=\tfrac12$, $X=$ number of heads: pmf $\big(\tfrac18,\tfrac38,\tfrac38,\tfrac18\big)$ for $x=0,1,2,3$, and the cdf is the step function jumping to $\tfrac18,\tfrac48,\tfrac78,1$
- $f(x)=c(1+x)^{-2}$, $x>0$: $\int_0^\infty c(1+x)^{-2}dx=c=1$, so $F(x)=\int_0^x(1+t)^{-2}dt=1-\dfrac1{1+x}$ and $P(1<X<2)=F(2)-F(1)=\tfrac23-\tfrac12=\tfrac16$

**Theorem 2.2.2 (Thm 2.2) — pdf from cdf**

For a continuous random variable, wherever $F$ is differentiable,

$$f(x)=\frac{d}{dx}F(x)$$

Proof. Differentiate $F(x)=\int_{-\infty}^x f(t)dt$ (by the fundamental theorem of calculus). Integration takes the density to the cdf, differentiation takes the cdf back to the density. $\blacksquare$

**Example 2.2.3 (Ex 2.4) — differentiating a cdf**

$F(x)=1-e^{-x/\lambda}$, $x>0$, gives $f(x)=F'(x)=\dfrac1\lambda e^{-x/\lambda}$, the exponential density.

## 2.3 Joint and Conditional Distributions

### 2.3.1 Joint distribution

**Definition 2.3.1 (Def 2.3) — joint density**

$f_{X,Y}(x,y)\ge0$ with total mass $1$, and for any region $A$

$$P\big((X,Y)\in A\big)=\begin{cases}\sum\sum_{(x,y)\in A}f_{X,Y}(x,y)&\text{discrete}\\ \iint_{A}f_{X,Y}(x,y)\,dx\,dy&\text{continuous}\end{cases}$$

**Definition 2.3.2 (Def 2.4) — joint cdf**

$$F(x_1,\dots,x_k)=P(X_1\le x_1,\dots,X_k\le x_k)$$

$$F(x_1,\dots,x_k)=\int_{-\infty}^{x_1}\!\cdots\int_{-\infty}^{x_k}f(t_1,\dots,t_k)\,dt_k\cdots dt_1,\qquad f=\frac{\partial^k}{\partial x_1\cdots\partial x_k}F \tag{2.6}$$

**Example 2.3.1 (Ex 2.6) — probability of a region, choosing the limits**

$f_{X,Y}(x,y)=xy\,e^{-(x+y)}$, $x,y>0$; find $P(X\ge2Y)$. Fix the outer variable and solve the inequality $x\ge2y$ for the inner one: with $x$ outside, $0\le y\le x/2$ and $x$ runs from $0$ to $\infty$:

$$P(X\ge2Y)=\int_0^{\infty}\!\int_0^{x/2}xy\,e^{-(x+y)}\,dy\,dx$$

$$\int_0^{x/2}ye^{-y}dy=1-\Big(\frac x2+1\Big)e^{-x/2}$$

$$P(X\ge2Y)=\int_0^\infty xe^{-x}\Big[1-\Big(\tfrac x2+1\Big)e^{-x/2}\Big]dx=1-\frac{20}{27}=\frac{7}{27}$$

Reversing the order gives $\int_0^\infty\!\int_{2y}^\infty xy\,e^{-(x+y)}dx\,dy$, with the same value.

**Theorem 2.3.1 (Thm 2.3) — characterization of a bivariate cdf**

$F(x_1,x_2)$ is a bivariate joint cdf if and only if

1. $F$ is non-decreasing and right-continuous in each variable
2. $F(-\infty,x_2)=F(x_1,-\infty)=0$ and $F(\infty,\infty)=1$
3. for all $a<b$, $c<d$: $\ F(b,d)-F(b,c)-F(a,d)+F(a,c)\ge0$ (rectangle probabilities are non-negative)

**Definition 2.3.3 — support**

$\{x:f(x)>0\}$, or $\{(x,y):f_{X,Y}(x,y)>0\}$ in two variables. Outside the support $f=0$, so only the part of the integration region inside the support contributes; probabilities and cdfs are therefore split into cases according to whether the point lies inside the support.

If the support is not a product (rectangle) set — e.g. $0<x<y$ — the variables are necessarily dependent. If it is a rectangle and the density factors as (function of $x$) $\times$ (function of $y$), they are independent.

**Example 2.3.2 (Ex 2.8) — joint cdf on a triangular support**

$f_{X,Y}(x,y)=8xy$ on the triangle $0<x<y<1$, and $F(x,y)=\iint_{u\le x,\,v\le y}f_{X,Y}(u,v)\,du\,dv$.

Above the diagonal ($x\le y$): $u$ from $0$ to $x$, then $v$ from $u$ to $y$:

$$F(x,y)=\int_0^x\!\int_u^y 8uv\,dv\,du=\int_0^x 4u(y^2-u^2)du=2x^2y^2-x^4$$

Below the diagonal ($x>y$): $x$ lies outside the support, so it cannot serve as an upper limit; take the limit from the support instead ($y$ in place of $x$), and $x$ drops out, leaving the marginal cdf of $Y$:

$$F(x,y)=\int_0^{y}\!\int_u^{y}8uv\,dv\,du=y^4=F_Y(y)$$

Only $y$ past the support ($0<x\le1<y$): replace the upper limit $y$ by the support edge $1$, and $y$ drops out, leaving the marginal cdf of $X$:

$$F(x,y)=\int_0^{x}\!\int_u^{1}8uv\,dv\,du=\int_0^x 4u(1-u^2)du=2x^2-x^4=F_X(x)$$

The cases are the possible combinations of where $x$ and $y$ sit relative to the support (before it, inside it, past it):

$$F(x,y)=\begin{cases}0,& x\le0\ \text{or}\ y\le0\\ 2x^2y^2-x^4,& 0<x\le y\le1\\ y^4,& 0<y<x,\ y\le1\\ 2x^2-x^4,& 0<x\le1<y\\ 1,& x\ge1\ \text{and}\ y\ge1\end{cases}$$

The pieces agree on the boundaries: both middle cases give $x^4$ at $x=y$, and all give $1$ at $x=y=1$.

**Method 2.3.1 — integration limits in higher dimensions**

With four or five variables the support cannot be drawn; the rule "if a limit is outside the support, take it from the support" becomes a $\min$:

- upper limit $=\min$(the wall you imposed, the limit the support allows)
- lower limit $=$ the floor set by the support (the preceding variable, or $0$)

For the support $0<X<Y<Z<1$, stacking the limits from the inside out,

$$F(x,y,z)=\int_0^{\min(x,y,z,1)}\!\int_{X}^{\min(y,z,1)}\!\int_{Y}^{\min(z,1)}f(X,Y,Z)\,dZ\,dY\,dX$$

### 2.3.2 Marginal distribution

**Theorem 2.3.2 (Thm 2.4) — marginal density**

Sum (integrate) out the other variables:

$$f_X(x)=\begin{cases}\sum_{y}f_{X,Y}(x,y)&\text{discrete}\\ \int_{-\infty}^{\infty}f_{X,Y}(x,y)\,dy&\text{continuous}\end{cases}$$

and symmetrically for $f_Y$.

Proof (continuous). The marginal cdf is

$$F_X(x)=P(X\le x)=P(X\le x,\ Y<\infty)=\int_{-\infty}^{x}\!\int_{-\infty}^{\infty}f_{X,Y}(t,y)\,dy\,dt$$

Differentiating with respect to $x$ gives $f_X(x)=\int_{-\infty}^{\infty}f_{X,Y}(x,y)\,dy$ (by Theorem 2.2.2). $\blacksquare$

**Example 2.3.3 (Ex 2.10) — marginal density**

$f_{X,Y}(x,y)=xy\,e^{-(x+y)}$, $x,y>0$:

$$f_X(x)=\int_0^\infty xye^{-(x+y)}dy=xe^{-x}\int_0^\infty ye^{-y}dy=xe^{-x},\qquad f_Y(y)=ye^{-y}$$

### 2.3.3 Conditional distribution

**Definition 2.3.4 (Def 2.5) — conditional density**

For $x$ with $f_X(x)>0$,

$$f_{Y\mid X}(y\mid x)=\frac{f_{X,Y}(x,y)}{f_X(x)}$$

and for $k>2$ variables

$$f(x_{i+1},\dots,x_k\mid x_1,\dots,x_i)=\frac{f(x_1,\dots,x_k)}{f(x_1,\dots,x_i)} \tag{2.7}$$

**Example 2.3.4 (Ex 2.12) — conditional density**

$f_{X,Y}(x,y)=x^2e^{-x(y+1)}$, $x,y>0$. The marginal is $f_X(x)=\int_0^\infty x^2e^{-x(y+1)}dy=xe^{-x}$, so

$$f_{Y\mid X}(y\mid x)=\frac{x^2e^{-x(y+1)}}{xe^{-x}}=xe^{-xy}$$

i.e. $Y\mid X=x\sim\mathrm{EXP}(1/x)$.

### 2.3.4 Independent random variables

**Definition 2.3.5 (Def 2.6) — independence**

$X,Y$ are independent when $P(X\in A,\ Y\in B)=P(X\in A)P(Y\in B)$ for all events $A,B$.

**Theorem 2.3.3 (Thm 2.5) — factorization criterion**

$X,Y$ are independent if and only if

$$f_{X,Y}(x,y)=f_X(x)f_Y(y)\qquad\text{for all }x,y$$

Proof (continuous case). For any intervals $A$ and $B$, independence says

$$\int_B\!\int_A f_{X,Y}(x,y)\,dx\,dy=P(X\in A,\ Y\in B)=P(X\in A)P(Y\in B)=\int_A f_X(x)dx\int_B f_Y(y)dy$$

$$=\int_B\!\int_A f_X(x)f_Y(y)\,dx\,dy$$

Since the two integrals agree for all intervals $A,B$, the integrands agree: $f_{X,Y}(x,y)=f_X(x)f_Y(y)$. The discrete case is the same with $\int$ replaced by $\sum$. $\blacksquare$

If $n$ random variables are independent with common density $f$ (a random sample), then

$$f(x_1,\dots,x_n)=f(x_1)f(x_2)\cdots f(x_n) \tag{2.8}$$

**Example 2.3.5 (Ex 2.13) — testing independence**

Joint pmf $f_{X,Y}(x,y)=\dfrac{2^{x+y}e^{-4}}{x!\,y!}$, $x,y=0,1,2,\dots$. The marginal is

$$f_X(x)=\sum_{y}\frac{2^{x+y}e^{-4}}{x!\,y!}=\frac{2^xe^{-4}}{x!}\sum_{y=0}^{\infty}\frac{2^y}{y!}=\frac{2^xe^{-4}}{x!}e^{2}=\frac{2^xe^{-2}}{x!}\quad\text{(by the exponential series)}$$

and likewise $f_Y(y)=\dfrac{2^ye^{-2}}{y!}$. Since $f_{X,Y}=f_Xf_Y$, $X$ and $Y$ are independent, each $\mathrm{POI}(2)$.

By contrast $f_{X,Y}=\tfrac12$ on $0\le x\le y\le2$ does not factor — the support is a triangle — so those variables are dependent.

## 2.4 Expectation

**Definition 2.4.1 (Def 2.7) — expectation (mean)**

When $E(|X|)<\infty$,

$$E(X)=\mu=\begin{cases}\sum_x x f(x)&\text{discrete}\\ \int_{-\infty}^{\infty}xf(x)dx&\text{continuous}\end{cases}$$

and for a function of $X$,

$$E[g(X)]=\sum_x g(x)f(x)\ \text{(discrete)},\qquad E[g(X)]=\int_{-\infty}^{\infty}g(x)f(x)dx\ \text{(continuous)} \tag{2.9, 2.10}$$

The two ways of computing $E(Y)$ for $Y=g(X)$ always agree, so the distribution of $Y$ need not be found first (law of the unconscious statistician, LOTUS):

$$E(Y)=\sum_y y f_Y(y)=\sum_x g(x)f_X(x)=E[g(X)]$$

**Example 2.4.1 (Ex 2.15, 2.16) — computing an expectation**

- three tosses, $X=$ number of heads: $E(X)=0\cdot\tfrac18+1\cdot\tfrac38+2\cdot\tfrac38+3\cdot\tfrac18=\tfrac32$
- $f(x)=\tfrac{x^2}{9}$ on $(0,3)$: $E(X)=\int_0^3 x\cdot\tfrac{x^2}{9}dx=\big[\tfrac{x^4}{36}\big]_0^3=\tfrac94$

### 2.4.1 Properties of expectation

**Theorem 2.4.1 (Thm 2.6) — linearity**

For constants $a,b,c$ and random variables $X,Y$: $E(c)=c$ and

$$E(aX+b)=aE(X)+b,\qquad E(aX+bY)=aE(X)+bE(Y),\qquad E\Big(\sum_i a_iX_i\Big)=\sum_i a_iE(X_i)$$

Proof (continuous). One variable, using $\int f=1$:

$$E(aX+b)=\int(ax+b)f(x)dx=a\int xf(x)dx+b\int f(x)dx=aE(X)+b$$

Two variables, integrating the joint density (marginalization):

$$E(aX+bY)=\iint(ax+by)f_{X,Y}\,dx\,dy=a\iint xf_{X,Y}\,dx\,dy+b\iint yf_{X,Y}\,dx\,dy=aE(X)+bE(Y)$$

The discrete case is identical with sums. $\blacksquare$

**Theorem 2.4.2 (Thm 2.7) — expectation of a product under independence**

If $X,Y$ are independent then $E(XY)=E(X)E(Y)$, and more generally $E[g(X)h(Y)]=E[g(X)]E[h(Y)]$.

Proof. Independence gives $f_{X,Y}=f_Xf_Y$ (by Theorem 2.3.3), so

$$E(XY)=\iint xy\,f_X(x)f_Y(y)\,dx\,dy=\Big(\int xf_Xdx\Big)\Big(\int yf_Ydy\Big)=E(X)E(Y)$$

The same separation applies to $g(X)h(Y)$. $\blacksquare$

The converse fails: $E(XY)=E(X)E(Y)$ (zero covariance, uncorrelated) does not imply independence. Counterexample: $U\sim U[0,2\pi)$ with $X=\cos U$, $Y=\sin U$ gives $E(X)=E(Y)=0$ and $E(XY)=\tfrac12E(\sin2U)=0$, so $\mathrm{Cov}(X,Y)=0$, yet $X^2+Y^2=1$ ties them together completely. Independence implies uncorrelatedness, not conversely.

### 2.4.2 Variance and covariance

**Definition 2.4.2 (Def 2.8) — variance, standard deviation**

$$\mathrm{Var}(X)=\sigma^2=E\big[(X-\mu)^2\big]=E(X^2)-[E(X)]^2,\qquad \sigma_X=\sqrt{\mathrm{Var}(X)}$$

The standardized variable $Z=\dfrac{X-\mu}{\sigma}$ has mean $0$ and variance $1$. (2.13)

**Theorem 2.4.3 (Thm 2.8) — properties of variance**

1. $\mathrm{Var}(aX+b)=a^2\mathrm{Var}(X)$
2. if $X_1,\dots,X_n$ are independent, $\mathrm{Var}\big(\sum_i X_i\big)=\sum_i\mathrm{Var}(X_i)$
3. $\mathrm{Var}(c)=0$
4. if $X_1,X_2$ are independent, $\mathrm{Var}(X_1-X_2)=\mathrm{Var}(X_1)+\mathrm{Var}(X_2)$

Proof. **(1)** $E(aX+b)=a\mu+b$, so

$$\mathrm{Var}(aX+b)=E\big[(aX+b-a\mu-b)^2\big]=E\big[a^2(X-\mu)^2\big]=a^2\mathrm{Var}(X)$$

**(2)** By Theorem 2.4.5, $\mathrm{Var}(\sum_i X_i)=\sum_i\mathrm{Var}(X_i)+2\sum_{j<k}\mathrm{Cov}(X_j,X_k)$, and independence makes every covariance $0$.

**(3)** Put $a=0$ in (1).

**(4)** Write $X_1-X_2=X_1+(-1)X_2$; by (1) $\mathrm{Var}(-X_2)=(-1)^2\mathrm{Var}(X_2)=\mathrm{Var}(X_2)$ — the coefficient is squared — so by (2)

$$\mathrm{Var}(X_1-X_2)=\mathrm{Var}(X_1)+\mathrm{Var}(X_2)\qquad\blacksquare$$

**Definition 2.4.3 (Def 2.9) — covariance**

$$\mathrm{Cov}(X,Y)=E\big[(X-\mu_X)(Y-\mu_Y)\big]=E(XY)-E(X)E(Y) \tag{2.14}$$

$\mathrm{Cov}(X,X)=\mathrm{Var}(X)$, and independence implies $\mathrm{Cov}(X,Y)=0$ (not conversely).

**Theorem 2.4.4 (Thm 2.9) — covariance under linear transformation**

$$\mathrm{Cov}(aX+b,\ cY+d)=ac\,\mathrm{Cov}(X,Y)$$

Proof. The means are $a\mu_X+b$ and $c\mu_Y+d$, so the shifts cancel inside the deviations:

$$\mathrm{Cov}(aX+b,cY+d)=E\big[(aX-a\mu_X)(cY-c\mu_Y)\big]=E\big[ac(X-\mu_X)(Y-\mu_Y)\big]=ac\,\mathrm{Cov}(X,Y)$$

$\blacksquare$

**Theorem 2.4.5 (Thm 2.10) — variance of a sum**

$$\mathrm{Var}\Big(\sum_{i=1}^{n}X_i\Big)=\sum_{i=1}^{n}\mathrm{Var}(X_i)+2\sum_{j<k}\mathrm{Cov}(X_j,X_k)$$

Proof (deviation form). With $S=\sum_i X_i$ and $E(S)=\sum_i\mu_i$,

$$\mathrm{Var}(S)=E\Big[\Big(\sum_i(X_i-\mu_i)\Big)^2\Big]=E\Big[\sum_i(X_i-\mu_i)^2+\sum_{j\ne k}(X_j-\mu_j)(X_k-\mu_k)\Big]$$

$$=\sum_i\mathrm{Var}(X_i)+2\sum_{j<k}\mathrm{Cov}(X_j,X_k)$$

Proof (expansion form). Expanding $E[(\sum X_i)^2]-[E(\sum X_i)]^2$,

$$E\Big[\sum_i X_i^2+2\sum_{j<k}X_jX_k\Big]-\Big[\sum_i E(X_i)\Big]^2$$

$$=\Big[\sum_i E(X_i^2)+2\sum_{j<k}E(X_jX_k)\Big]-\Big[\sum_i\{E(X_i)\}^2+2\sum_{j<k}E(X_j)E(X_k)\Big]$$

$$=\sum_i\big\{E(X_i^2)-[E(X_i)]^2\big\}+2\sum_{j<k}\big\{E(X_jX_k)-E(X_j)E(X_k)\big\}=\sum_i\mathrm{Var}(X_i)+2\sum_{j<k}\mathrm{Cov}(X_j,X_k)$$

$\blacksquare$

### 2.4.3 Conditional expectation

**Definition 2.4.4 (Def 2.10) — conditional expectation**

$$E(Y\mid X=x)=\sum_y y\,f_{Y\mid X}(y\mid x)\ \text{(discrete)},\qquad E(Y\mid X=x)=\int_{-\infty}^{\infty}y\,f_{Y\mid X}(y\mid x)\,dy\ \text{(continuous)}$$

**Example 2.4.2 (Ex 2.24) — uniform on a triangle**

$f_{X,Y}(x,y)=\tfrac12$ on $0\le x\le y\le2$. The marginals are

$$f_X(x)=\int_x^2\tfrac12dy=\frac{2-x}{2}\ (0\le x\le2),\qquad f_Y(y)=\int_0^y\tfrac12dx=\frac y2\ (0\le y\le2)$$

$$f_{Y\mid X}(y\mid x)=\frac{1/2}{(2-x)/2}=\frac1{2-x}\quad(x\le y\le2)\ \Rightarrow\ Y\mid X=x\sim U(x,2)$$

$$E(Y\mid X=x)=\int_x^2 y\cdot\frac1{2-x}dy=\frac1{2-x}\cdot\frac{4-x^2}{2}=\frac{2+x}{2}$$

the midpoint of $[x,2]$, as expected for a uniform conditional distribution.

**Theorem 2.4.6 (Thm 2.11) — law of total expectation**

$$E\big[E(Y\mid X)\big]=E(Y)$$

Proof (continuous).

$$E\big[E(Y\mid X)\big]=\int E(Y\mid x)f_X(x)dx=\int\Big(\int y\,\frac{f_{X,Y}(x,y)}{f_X(x)}dy\Big)f_X(x)dx$$

$$=\iint y\,f_{X,Y}(x,y)\,dy\,dx=\int y\,f_Y(y)dy=E(Y)\qquad\blacksquare$$

Note that $E(Y\mid x)$ is a function of $x$ alone ($y$ has been integrated out), so the outer integral is just the expectation of a function of $X$.

**Theorem 2.4.7 (Thm 2.12) — conditional expectation under independence**

If $X,Y$ are independent then $E(Y\mid x)=E(Y)$ and $E(X\mid y)=E(X)$.

Proof. Independence gives $f_{Y\mid X}(y\mid x)=\dfrac{f_X(x)f_Y(y)}{f_X(x)}=f_Y(y)$, so $E(Y\mid x)=\int yf_Y(y)dy=E(Y)$. $\blacksquare$

**Definition 2.4.5 (Def 2.11) — conditional variance**

$$\mathrm{Var}(Y\mid x)=E\big[(Y-E(Y\mid x))^2\mid x\big]=E(Y^2\mid x)-[E(Y\mid x)]^2 \tag{2.15}$$

**Example 2.4.3 (Ex 2.26) — conditional variance**

$f_{Y\mid X}(y\mid x)=xe^{-xy}$, $y>0$, i.e. $Y\mid X=x\sim\mathrm{EXP}(1/x)$, so $E(Y\mid x)=\tfrac1x$ and $E(Y^2\mid x)=\tfrac2{x^2}$:

$$\mathrm{Var}(Y\mid x)=\frac2{x^2}-\frac1{x^2}=\frac1{x^2}$$

**Theorem 2.4.8 (Thm 2.13) — law of total variance**

$$\mathrm{Var}(Y)=E\big[\mathrm{Var}(Y\mid X)\big]+\mathrm{Var}\big[E(Y\mid X)\big]$$

Proof. Taking the expectation of $\mathrm{Var}(Y\mid X)=E(Y^2\mid X)-[E(Y\mid X)]^2$,

$$E\big[\mathrm{Var}(Y\mid X)\big]=E(Y^2)-E\big[(E(Y\mid X))^2\big]\quad\text{(by Theorem 2.4.6)}$$

$$\mathrm{Var}\big[E(Y\mid X)\big]=E\big[(E(Y\mid X))^2\big]-[E(Y)]^2$$

Adding gives $E(Y^2)-[E(Y)]^2=\mathrm{Var}(Y)$. $\blacksquare$

**Corollary 2.4.1 — consequences**

1. $E[\mathrm{Var}(Y\mid X)]\le\mathrm{Var}(Y)$, since $\mathrm{Var}[E(Y\mid X)]\ge0$: using the information in $X$ can only reduce the variance of $Y$
2. if $X,Y$ are independent then $E(Y\mid X)=E(Y)$ is constant, so $\mathrm{Var}[E(Y\mid X)]=0$ and $E[\mathrm{Var}(Y\mid X)]=\mathrm{Var}(Y)$ — conditioning reduces nothing

### 2.4.4 Probability inequalities

**Theorem 2.4.9 (Thm 2.14) — Markov's inequality**

For $u(X)\ge0$ and $c>0$,

$$P\big[u(X)\ge c\big]\le\frac{E[u(X)]}{c} \tag{2.16}$$

Proof (continuous). With $A=\{x:u(x)\ge c\}$,

$$E[u(X)]=\int u(x)f(x)dx\ \ge\ \int_A u(x)f(x)dx\ \ge\ \int_A c\,f(x)dx=c\,P[u(X)\ge c]$$

Divide by $c$. $\blacksquare$

**Theorem 2.4.10 (Thm 2.15) — Chebyshev's inequality**

$$P\big(|X-\mu|\ge k\sigma\big)\le\frac1{k^2},\qquad P\big(|X-\mu|<k\sigma\big)\ge1-\frac1{k^2} \tag{2.17}$$

Proof. Put $u(X)=(X-\mu)^2$ and $c=k^2\sigma^2$ in Markov's inequality:

$$P\big[(X-\mu)^2\ge k^2\sigma^2\big]\le\frac{E[(X-\mu)^2]}{k^2\sigma^2}=\frac{\sigma^2}{k^2\sigma^2}=\frac1{k^2}$$

and the left side is $P(|X-\mu|\ge k\sigma)$. $\blacksquare$

**Example 2.4.4 (Ex 2.27) — applying Chebyshev**

$\mu=25$, $\sigma^2=16$ ($\sigma=4$). The event $17\le X\le33$ is $|X-25|<8=2\sigma$, so with $k=2$

$$P(17\le X\le33)\ge1-\frac1{2^2}=0.75$$

whatever the shape of the distribution.

**Theorem 2.4.11 (Thm 2.16) — Cauchy–Schwarz inequality**

$$[E(XY)]^2\le E(X^2)E(Y^2)$$

with equality if and only if $P(Y=cX)=1$ for some constant $c$.

Proof. For any real $t$, $(tX-Y)^2\ge0$, so its expectation is non-negative:

$$0\le E\big[(tX-Y)^2\big]=E(X^2)t^2-2E(XY)t+E(Y^2)=At^2-2Bt+C$$

with $A=E(X^2)$, $B=E(XY)$, $C=E(Y^2)$ constants in $t$ (by linearity). An upward parabola that is non-negative for every $t$ cannot cross the axis twice, so its discriminant is $\le0$:

$$\frac D4=B^2-AC\le0\ \Rightarrow\ B^2\le AC$$

Substituting back gives $[E(XY)]^2\le E(X^2)E(Y^2)$. $\blacksquare$

Equality means $D=0$: the parabola touches the axis, so $E[(t_0X-Y)^2]=0$ for some $t_0$, which forces $Y=t_0X$ with probability $1$ — perfect linear dependence (after centering, $\rho=\pm1$).

## 2.5 Discrete Distributions

**Definition 2.5.1 — Bernoulli distribution**

$X\sim\mathrm{Ber}(p)$ with $f(1)=p$, $f(0)=q=1-p$:

$$E(X)=p,\qquad \mathrm{Var}(X)=p(1-p)=pq \tag{2.18}$$

since $E(X)=1\cdot p+0\cdot q=p$ and $E(X^2)=p$, so $\mathrm{Var}=p-p^2=pq$.

**Definition 2.5.2 — indicator function**

$$I_A(\omega)=\begin{cases}1,&\omega\in A\\ 0,&\omega\notin A\end{cases}$$

Since $P(I_A=1)=P(A)$, $I_A\sim\mathrm{Ber}(P(A))$ — the simplest discrete random variable.

**Theorem 2.5.1 — properties of indicators**

1. idempotent: $I_A^2=I_A$, and $I_A^k=I_A$ for $k\ge1$
2. complement: $I_{A^c}=1-I_A$
3. intersection is a product: $I_{A\cap B}=I_AI_B=\min(I_A,I_B)$
4. union by inclusion–exclusion: $I_{A\cup B}=I_A+I_B-I_AI_B=\max(I_A,I_B)$
5. monotone: $A\subseteq B\iff I_A\le I_B$

**Theorem 2.5.2 — moments of indicators**

Expectation is probability — the bridge that turns probabilities into integrals and sums:

$$E(I_A)=1\cdot P(A)+0\cdot P(A^c)=P(A)$$

By idempotence $E(I_A^2)=E(I_A)=P(A)$, so

$$\mathrm{Var}(I_A)=E(I_A^2)-[E(I_A)]^2=P(A)-P(A)^2=P(A)P(A^c)$$

and since $I_AI_B=I_{A\cap B}$,

$$\mathrm{Cov}(I_A,I_B)=P(A\cap B)-P(A)P(B)$$

so independence of $A,B$ is equivalent to $\mathrm{Cov}(I_A,I_B)=0$.

**Method 2.5.1 — the indicator method**

Split a count into a sum of indicators; linearity alone then gives the mean, with no independence needed:

$$X=\sum_{i=1}^n I_{A_i}\ \Longrightarrow\ E(X)=\sum_{i=1}^n E(I_{A_i})=\sum_{i=1}^n P(A_i)$$

The binomial distribution is of this form with $A_i=\{$success on trial $i\}$ and $P(A_i)=p$, giving $E(X)=np$. The formula holds even when the $A_i$ are dependent, which is the strength of the method.

**Definition 2.5.3 — binomial distribution**

$X\sim B(n,p)$, the number of successes in $n$ independent Bernoulli trials, $X=\sum_{i=1}^n X_i$:

$$f(x)=\binom{n}{x}p^xq^{n-x},\quad x=0,1,\dots,n \tag{2.20}$$

$$E(X)=np,\qquad \mathrm{Var}(X)=npq \tag{2.21, 2.22}$$

Proof. $X=\sum X_i$ with independent $X_i\sim\mathrm{Ber}(p)$, so $E(X)=\sum E(X_i)=np$ (by linearity) and $\mathrm{Var}(X)=\sum\mathrm{Var}(X_i)=npq$ (by independence, Theorem 2.4.3(2)). $\blacksquare$

The binomial counts successes under sampling with replacement (or from an infinite population); without replacement the success probability changes from trial to trial and the count is hypergeometric.

**Definition 2.5.4 — multinomial distribution**

Each of $n$ independent trials falls in one of $k$ categories with probabilities $p_1,\dots,p_k$, $\sum_i p_i=1$; with $X_i$ the count in category $i$ and $\sum_i X_i=n$,

$$f(x_1,\dots,x_k)=\frac{n!}{x_1!\,x_2!\cdots x_k!}p_1^{x_1}p_2^{x_2}\cdots p_k^{x_k},\qquad \sum_{i=1}^k x_i=n$$

For $k=2$ this is the binomial with $(X_1,X_2)=(X,n-X)$; the coefficient is the multinomial coefficient of Theorem 1.4.4.

**Theorem 2.5.3 — moments of the multinomial**

1. the marginals are binomial: $X_i\sim B(n,p_i)$, so $E(X_i)=np_i$ and $\mathrm{Var}(X_i)=np_i(1-p_i)$
2. the covariances are negative: $\mathrm{Cov}(X_i,X_j)=-np_ip_j$ for $i\ne j$

Proof. **(1)** Treating category $i$ as success and all others as failure gives $n$ independent trials with success probability $p_i$, so $X_i\sim B(n,p_i)$.

**(2)** Let $I_{im}=1$ when trial $m$ falls in category $i$, so $X_i=\sum_{m=1}^n I_{im}$. A single trial belongs to one category only, so for $i\ne j$ we have $I_{im}I_{jm}=0$ and hence $E(I_{im}I_{jm})=0$:

$$\mathrm{Cov}(I_{im},I_{jm})=0-p_ip_j=-p_ip_j$$

Different trials are independent, so their covariances vanish, and

$$\mathrm{Cov}(X_i,X_j)=\sum_{m=1}^n\mathrm{Cov}(I_{im},I_{jm})=-np_ip_j\qquad\blacksquare$$

The negative sign is natural: $\sum_i X_i=n$ is fixed, so a large count in one category forces smaller counts elsewhere.

**Example 2.5.1 (Ex 2.31) — multinomial covariance via conditional expectation**

$(X_1,X_2,X_3)\sim\mathrm{MULT}(n,p_1,p_2,p_3)$. The marginals are $X_1\sim B(n,p_1)$, $X_2\sim B(n,p_2)$, and the conditional distribution is

$$X_2\mid X_1=x_1\ \sim\ B\Big(n-x_1,\ \frac{p_2}{1-p_1}\Big)$$

since once $X_1=x_1$ is fixed, the remaining $n-x_1$ trials fall in category $2$ with probability $p_2/(1-p_1)$. Then, $X_1$ being constant given $X_1$,

$$E(X_1X_2)=E\big[X_1E(X_2\mid X_1)\big]=E\Big[X_1(n-X_1)\frac{p_2}{1-p_1}\Big]=\frac{p_2}{1-p_1}\big[nE(X_1)-E(X_1^2)\big]$$

(by Theorem 2.4.6). With $E(X_1)=np_1$ and $E(X_1^2)=np_1(1-p_1)+n^2p_1^2$, the bracket is $n(n-1)p_1(1-p_1)$, so

$$E(X_1X_2)=\frac{p_2}{1-p_1}n(n-1)p_1(1-p_1)=n(n-1)p_1p_2$$

$$\mathrm{Cov}(X_1,X_2)=E(X_1X_2)-E(X_1)E(X_2)=n(n-1)p_1p_2-n^2p_1p_2=-np_1p_2$$

the same result as by the indicator method.

**Definition 2.5.5 — Poisson distribution**

$X\sim\mathrm{POI}(\lambda)$, the number of rare events in a unit interval with mean rate $\lambda$:

$$f(x)=\frac{e^{-\lambda}\lambda^x}{x!},\quad x=0,1,2,\dots \tag{2.24}$$

$$E(X)=\lambda,\qquad \mathrm{Var}(X)=\lambda \tag{2.25}$$

$$E(X)=\sum_{x\ge1}x\frac{e^{-\lambda}\lambda^x}{x!}=\lambda\sum_{x\ge1}\frac{e^{-\lambda}\lambda^{x-1}}{(x-1)!}=\lambda$$

and from $E[X(X-1)]=\lambda^2$ we get $E(X^2)=\lambda^2+\lambda$, hence $\mathrm{Var}(X)=\lambda$.

**Theorem 2.5.4 — derivation from the Poisson process**

Model the count of a rare event $A$ by three assumptions on a short interval of length $h$:

1. one occurrence has probability proportional to length: $\lim_{h\to0}P(1\text{ occurrence})/h=\lambda$
2. two or more occurrences are negligible: $\lim_{h\to0}P(\ge2)/h=0$, written $o(h)$
3. counts in disjoint intervals are independent

Split $[0,t]$ into $n$ subintervals of length $h=t/n$. Each is a Bernoulli trial (one occurrence with probability $\approx\lambda t/n$, or none), and the subintervals are independent by (3), so the count is approximately $B(n,\lambda t/n)$:

$$P(X(t)=x)=\binom nx\Big(\frac{\lambda t}{n}\Big)^x\Big(1-\frac{\lambda t}{n}\Big)^{n-x}=\frac{(\lambda t)^x}{x!}\cdot\frac{n(n-1)\cdots(n-x+1)}{n^x}\Big(1-\frac{\lambda t}{n}\Big)^{n}\Big(1-\frac{\lambda t}{n}\Big)^{-x}$$

The three factors tend to $1$, $e^{-\lambda t}$ and $1$ (by the exponential limit), so

$$P(X(t)=x)\ \longrightarrow\ e^{-\lambda t}\frac{(\lambda t)^x}{x!},\qquad x=0,1,2,\dots$$

that is, the count on $[0,t]$ is $\mathrm{POI}(\lambda t)$, and $t=1$ gives (2.24).

Why the "two or more" term disappears: by assumption (2) a single subinterval has probability $o(t/n)$ of two or more occurrences, and over all $n$ subintervals

$$P\Big(\bigcup_{i=1}^n\{\ge2\text{ in subinterval }i\}\Big)\le\sum_{i=1}^n o\Big(\frac tn\Big)=n\cdot o\Big(\frac tn\Big)=t\cdot\frac{o(t/n)}{t/n}\ \longrightarrow\ 0$$

(by Boole's inequality, Theorem 1.2.4, and $o(t/n)\big/(t/n)\to0$). In the limit each subinterval is a $0$-or-$1$ Bernoulli trial.

**Theorem 2.5.5 (Thm 2.17) — Poisson approximation to the binomial**

Holding $np=\lambda$ fixed as $n\to\infty$ and $p\to0$,

$$\binom nx p^x(1-p)^{n-x}\ \longrightarrow\ \frac{e^{-\lambda}\lambda^x}{x!}$$

Proof. Substituting $p=\lambda/n$,

$$\binom nx p^x(1-p)^{n-x}=\frac{n!}{x!(n-x)!}\frac{\lambda^x}{n^x}\Big(1-\frac\lambda n\Big)^{n-x}=\frac{\lambda^x}{x!}\cdot\frac{n!}{(n-x)!\,n^x}\Big(1-\frac\lambda n\Big)^{n}\Big(1-\frac\lambda n\Big)^{-x}$$

where the three factors tend to $1$, $e^{-\lambda}$ and $1$, giving $e^{-\lambda}\lambda^x/x!$. $\blacksquare$

For large $n$ and small $p$ (roughly $n\ge20$, $p\le0.05$, or $n\ge100$, $np\le10$) the binomial coefficient is awkward and $\mathrm{POI}(np)$ is used instead.

**Definition 2.5.6 — geometric distribution**

$X\sim\mathrm{GEO}(p)$, the number of trials up to and including the first success:

$$f(x)=p(1-p)^{x-1},\quad x=1,2,\dots \tag{2.27}$$

$$E(X)=\frac1p,\qquad \mathrm{Var}(X)=\frac{1-p}{p^2}=\frac{q}{p^2} \tag{2.28, 2.29}$$

**Theorem 2.5.6 (Thm 2.18) — memorylessness**

$$P(X>j+k\mid X>j)=P(X>k)$$

Proof. First the tail probability; factoring out $(1-p)^k$ leaves the total probability of a geometric distribution, which is $1$:

$$P(X>k)=\sum_{j=k+1}^{\infty}p(1-p)^{j-1}=(1-p)^k\sum_{l=1}^{\infty}p(1-p)^{l-1}=(1-p)^k$$

Since $\{X>j+k\}\subset\{X>j\}$, the intersection of the two events is $\{X>j+k\}$, so

$$P(X>j+k\mid X>j)=\frac{P(X>j+k)}{P(X>j)}=\frac{(1-p)^{j+k}}{(1-p)^{j}}=(1-p)^k=P(X>k)\qquad\blacksquare$$

Equivalently $P(X>j+k)=P(X>j)P(X>k)$.

**Definition 2.5.7 — negative binomial distribution**

$X\sim\mathrm{NB}(r,p)$, the number of trials up to and including the $r$-th success:

$$f(x)=\binom{x-1}{r-1}p^r(1-p)^{x-r},\quad x=r,r+1,\dots \tag{2.30}$$

$$E(X)=\frac rp,\qquad \mathrm{Var}(X)=\frac{rq}{p^2} \tag{2.31, 2.32}$$

Since $X=\sum_{i=1}^r X_i$ is a sum of $r$ independent geometric variables, the mean and variance are $r$ times the geometric ones.

**Theorem 2.5.7 — relation to the binomial**

For a constant $n\ge r$, with $Y\sim B(n,p)$ the number of successes in $n$ trials,

$$P(X\le n)=P(Y\ge r)=\sum_{k=r}^{n}\binom nk p^k(1-p)^{n-k}$$

The two events coincide: "the $r$-th success occurs within the first $n$ trials" is exactly "the first $n$ trials already contain at least $r$ successes". The negative binomial stops at the $r$-th success and records the trial number; the binomial runs all $n$ trials and counts successes.

**Example 2.5.2 (Ex 2.35) — negative binomial (best-of-seven series)**

Best of seven (first to four wins), each game won with $p=0.7$; winning in exactly five games means $X=5$ with $r=4$:

$$P(X=5)=\binom43(0.7)^4(0.3)^1=0.3108$$

**Definition 2.5.8 — hypergeometric distribution**

From $r$ red and $w$ white balls ($r+w=N$), draw $n$ at random without replacement; $X$ is the number of red balls drawn:

$$f(x)=\frac{\binom rx\binom{N-r}{n-x}}{\binom Nn},\qquad x=\max(0,n-N+r),\dots,\min(r,n) \tag{2.33}$$

written $X\sim\mathrm{HYP}(n,N,r)$.

$$E(X)=n\frac rN,\qquad \mathrm{Var}(X)=n\frac rN\cdot\frac{N-r}{N}\cdot\frac{N-n}{N-1} \tag{2.34, 2.35}$$

Reading the formula: $\binom Nn$ counts all draws, $\binom rx$ the ways of taking $x$ of the $r$ red balls, $\binom{N-r}{n-x}$ the ways of taking the rest from the white ones. The support is bounded above by $\min(r,n)$ and below by $\max(0,n-N+r)$, since at least $n-(N-r)$ of the drawn balls must be red.

Proof (mean and variance). **(1)** Using the absorption identities $x\binom rx=r\binom{r-1}{x-1}$ and $\binom Nn=\frac Nn\binom{N-1}{n-1}$, the $x=0$ term drops out and

$$E(X)=\sum_{x=0}^{n}x\frac{\binom rx\binom{N-r}{n-x}}{\binom Nn}=n\frac rN\sum_{x=1}^{n}\frac{\binom{r-1}{x-1}\binom{N-r}{n-x}}{\binom{N-1}{n-1}}=n\frac rN\sum_{y=0}^{n-1}\frac{\binom{r-1}{y}\binom{(N-1)-(r-1)}{(n-1)-y}}{\binom{N-1}{n-1}}=n\frac rN$$

since with $y=x-1$ the last sum is the total probability of $\mathrm{HYP}(n-1,N-1,r-1)$, namely $1$.

**(2)** Similarly, using $x(x-1)\binom rx=r(r-1)\binom{r-2}{x-2}$,

$$E[X(X-1)]=n(n-1)\frac{r(r-1)}{N(N-1)}\sum_{y=0}^{n-2}\frac{\binom{r-2}{y}\binom{(N-2)-(r-2)}{(n-2)-y}}{\binom{N-2}{n-2}}=n(n-1)\frac{r(r-1)}{N(N-1)}$$

and from $\mathrm{Var}(X)=E[X(X-1)]+E(X)-[E(X)]^2$,

$$\mathrm{Var}(X)=n(n-1)\frac{r(r-1)}{N(N-1)}+n\frac rN-n^2\frac{r^2}{N^2}=n\frac rN\Big[(n-1)\frac{r-1}{N-1}+1-\frac{nr}{N}\Big]=\frac{nr}{N}\cdot\frac{(N-r)(N-n)}{N(N-1)}$$

$\blacksquare$

**Theorem 2.5.8 — finite population correction and the binomial limit**

With $p=r/N$,

$$\mathrm{Var}(X)=np(1-p)\cdot\frac{N-n}{N-1}$$

the binomial variance $np(1-p)$ times the finite population correction $\frac{N-n}{N-1}$ (equal to $1$ when $n=1$, and shrinking the variance as $n$ grows). Moreover, as $N\to\infty$ and $r\to\infty$ with $r/N\to p$,

$$\frac{\binom rx\binom{N-r}{n-x}}{\binom Nn}\ \longrightarrow\ \binom nx p^x(1-p)^{n-x}$$

so for large $N$ and $r$ the hypergeometric may be approximated by the binomial with $p=r/N$.

Proof. **(1)** Writing all three binomial coefficients as factorials,

$$P(X=x)=\frac{\dfrac{r!}{x!(r-x)!}\cdot\dfrac{(N-r)!}{(n-x)!(N-r-n+x)!}}{\dfrac{N!}{n!(N-n)!}}=\binom nx\cdot\frac{r!}{(r-x)!}\cdot\frac{(N-r)!}{(N-r-n+x)!}\cdot\frac{(N-n)!}{N!}$$

collecting $x!,(n-x)!,n!$ into $\binom nx$. The remaining three ratios are products of consecutive integers, with $x+(n-x)=n$ factors in the numerator and $n$ in the denominator:

$$=\binom nx\cdot\frac{r(r-1)\cdots(r-x+1)\cdot(N-r)(N-r-1)\cdots(N-r-n+x+1)}{N(N-1)\cdots(N-n+1)}$$

**(2)** Divide numerator and denominator by $N^n=N^x\cdot N^{n-x}$, one $N$ per factor.

**(3)** As $N\to\infty$ with $n,x$ fixed and $k/N\to0$ for any constant $k$:

- denominator ($n$ factors): $1,1-\tfrac1N,\dots,1-\tfrac{n-1}N$ all tend to $1$
- first part of the numerator ($x$ factors): $\dfrac{r-k}{N}=\dfrac rN-\dfrac kN\to p$, giving $p^x$
- second part ($n-x$ factors): $\dfrac{N-r-k}{N}=1-\dfrac rN-\dfrac kN\to1-p$, giving $(1-p)^{n-x}$

Hence $P(X=x)\to\binom nx p^x(1-p)^{n-x}$. $\blacksquare$

Intuitively, when $N$ is large removing one ball barely changes the proportion $r/N\approx p$, so drawing without replacement behaves like drawing with replacement.

**Example 2.5.3 (Ex 2.36) — hypergeometric versus binomial approximation**

A box of $100$ contains $30$ defectives; $10$ are drawn without replacement, so $X\sim\mathrm{HYP}(10,100,30)$:

$$P(X\le2)=\sum_{x\le2}\frac{\binom{30}{x}\binom{70}{10-x}}{\binom{100}{10}}=0.372857$$

The binomial approximation with $p=30/100=0.3$ gives

$$P(X\le2)=\sum_{x\le2}\binom{10}{x}(0.3)^x(0.7)^{10-x}=0.382783$$

reasonably close, since $n/N$ is small.

## 2.6 Continuous Distributions

**Definition 2.6.1 — uniform distribution**

$X\sim U(a,b)$:

$$f(x)=\frac1{b-a},\ a\le x\le b;\qquad F(x)=\frac{x-a}{b-a},\ a\le x\le b \tag{2.36}$$

$$E(X)=\frac{a+b}{2},\qquad \mathrm{Var}(X)=\frac{(b-a)^2}{12}$$

**Definition 2.6.2 — exponential distribution**

$X\sim\mathrm{EXP}(\theta)$, the waiting time between events of a Poisson process, with scale parameter $\theta=1/\lambda$:

$$f(x)=\frac1\theta e^{-x/\theta},\ x>0;\qquad F(x)=1-e^{-x/\theta}$$

$$E(X)=\theta,\qquad \mathrm{Var}(X)=\theta^2 \tag{2.41}$$

**Theorem 2.6.1 — derivation as a waiting time**

In a Poisson process of rate $\lambda$ the count over time $t$ is $Y_t\sim\mathrm{POI}(\lambda t)$. For the waiting time $X$ to the first event, "$X>t$" is the same event as "no occurrence in $[0,t]$":

$$P(X>t)=P(Y_t=0)=\frac{e^{-\lambda t}(\lambda t)^0}{0!}=e^{-\lambda t}$$

$$F_X(t)=1-e^{-\lambda t},\qquad f_X(t)=F_X'(t)=\lambda e^{-\lambda t},\ t>0$$

With $\theta=1/\lambda$ this is the exponential density: the gaps between events of a Poisson process are exponential.

**Theorem 2.6.2 (Thm 2.19) — memorylessness**

$$P(X>a+t\mid X>a)=P(X>t)$$

Proof. Using $P(X>x)=e^{-x/\theta}$, and $\{X>a+t\}\subset\{X>a\}$ so that the intersection is $\{X>a+t\}$,

$$P(X>a+t\mid X>a)=\frac{P(X>a+t)}{P(X>a)}=\frac{e^{-(a+t)/\theta}}{e^{-a/\theta}}=e^{-t/\theta}=P(X>t)\qquad\blacksquare$$

If $X$ is the lifetime of a component, the chance of surviving a further $t$ having already lasted $a$ equals that of a new component: the age is not remembered. The exponential is the only continuous distribution with this property.

**Definition 2.6.3 — normal distribution**

$X\sim N(\mu,\sigma^2)$:

$$f(x)=\frac1{\sqrt{2\pi}\,\sigma}\exp\Big[-\frac{(x-\mu)^2}{2\sigma^2}\Big],\quad-\infty<x<\infty \tag{2.42}$$

$$E(X)=\mu,\qquad \mathrm{Var}(X)=\sigma^2 \tag{2.44}$$

bell-shaped and symmetric about $x=\mu$.

**Theorem 2.6.3 — derivation from the target model (Herschel–Maxwell)**

For a shot at the centre of a target with coordinates $(X,Y)$, assume only

1. iid: the horizontal and vertical errors $X,Y$ are independent with the same density $f$
2. rotational symmetry: the density depends only on the distance $r=\sqrt{x^2+y^2}$

Independence gives the joint density $f(x)f(y)$; symmetry makes it a function $h$ of $x^2+y^2$:

$$f(x)f(y)=h(x^2+y^2)$$

Taking logarithms with $L=\ln f$ gives $L(x)+L(y)=\ln h(x^2+y^2)$, and differentiating with respect to $x$ and to $y$,

$$\frac{L'(x)}{2x}=\frac{L'(y)}{2y}\qquad\text{for all }x,y$$

The left side depends only on $x$ and the right only on $y$, so both must equal a constant, say $-k$. Then $L'(x)=-2kx$, so $L(x)=-kx^2+c$ and $f(x)=Ce^{-kx^2}$. The two constants follow from two conditions.

Total integral $1$, by the Gaussian integral $\int_{-\infty}^\infty e^{-kx^2}dx=\sqrt{\pi/k}$:

$$C\sqrt{\frac\pi k}=1\ \Rightarrow\ C=\sqrt{\frac k\pi}$$

Variance $\sigma^2$, using $\int_{-\infty}^\infty x^2e^{-kx^2}dx=\frac1{2k}\sqrt{\pi/k}$:

$$\sigma^2=\int_{-\infty}^\infty x^2Ce^{-kx^2}dx=\frac1{2k}\ \Rightarrow\ k=\frac1{2\sigma^2}$$

so $C=\sqrt{k/\pi}=\dfrac1{\sqrt{2\pi}\,\sigma}$ and

$$f(x)=\frac1{\sqrt{2\pi}\,\sigma}e^{-x^2/2\sigma^2}\ \sim\ N(0,\sigma^2)$$

Shifting the centre to $\mu$ gives (2.42).

**Theorem 2.6.4 — linear transformation of a normal**

If $X\sim N(\mu,\sigma^2)$ and $Y=aX+b$ with $a\ne0$, then

$$Y\sim N\big(a\mu+b,\ a^2\sigma^2\big) \tag{2.43}$$

In particular $a=1/\sigma$, $b=-\mu/\sigma$ gives $Z=\dfrac{X-\mu}{\sigma}\sim N(0,1)$.

Proof. Take $a>0$ and write the cdf of $Y$ as a probability for $X$:

$$F_Y(y)=P(aX+b\le y)=P\Big(X\le\frac{y-b}{a}\Big)=\int_{-\infty}^{(y-b)/a}\frac1{\sqrt{2\pi}\,\sigma}e^{-(x-\mu)^2/2\sigma^2}dx$$

Substituting $t=ax+b$ (so $dx=\tfrac1a dt$ and the upper limit becomes $t=y$),

$$F_Y(y)=\int_{-\infty}^{y}\frac1{\sqrt{2\pi}\,a\sigma}\exp\Big\{-\frac{[t-(a\mu+b)]^2}{2a^2\sigma^2}\Big\}dt$$

$$f_Y(y)=\frac{d}{dy}F_Y(y)=\frac1{\sqrt{2\pi}\,a\sigma}\exp\Big\{-\frac{[y-(a\mu+b)]^2}{2a^2\sigma^2}\Big\}$$

the $N(a\mu+b,a^2\sigma^2)$ density. For $a<0$ replace $a$ by $|a|$. $\blacksquare$

**Theorem 2.6.5 (Thm 2.20) — standardization**

$Z=\dfrac{X-\mu}{\sigma}\sim N(0,1)$, and with $\phi,\Phi$ the standard normal pdf and cdf,

$$F(x)=\Phi\Big(\frac{x-\mu}{\sigma}\Big),\qquad \Phi(-z)=1-\Phi(z)$$

Proof (normalization of the normal density). The density integrates to $1$ provided $\int_{-\infty}^\infty e^{-x^2/2}dx=\sqrt{2\pi}$. This integral $I$ has no elementary antiderivative, so square it and use polar coordinates:

$$I^2=\Big(\int_{-\infty}^\infty e^{-x^2/2}dx\Big)\Big(\int_{-\infty}^\infty e^{-y^2/2}dy\Big)=\iint_{\mathbb R^2}e^{-(x^2+y^2)/2}dx\,dy$$

With $x=r\cos\theta$, $y=r\sin\theta$ and $dx\,dy=r\,dr\,d\theta$, the factor $r$ makes the inner integral elementary:

$$I^2=\int_0^{2\pi}\!\int_0^\infty e^{-r^2/2}r\,dr\,d\theta=\int_0^{2\pi}1\,d\theta=2\pi\ \Rightarrow\ I=\sqrt{2\pi}$$

so $\int\frac1{\sqrt{2\pi}}e^{-x^2/2}dx=1$, and the general case follows by the substitution $z=(x-\mu)/\sigma$. This is the origin of the $\sqrt{2\pi}$ in the normalizing constant. $\blacksquare$

**Definition 2.6.4 — standard normal cdf**

$$\Phi(z)=P(Z\le z)=\int_{-\infty}^{z}\phi(t)dt=\int_{-\infty}^{z}\frac1{\sqrt{2\pi}}e^{-t^2/2}dt$$

There is no closed elementary form, so values come from tables or numerical computation; $\Phi'(z)=\phi(z)$.

**Theorem 2.6.6 — properties of $\Phi$**

1. $\Phi(-\infty)=0$, $\Phi(\infty)=1$, and $\Phi$ is continuous and increasing
2. symmetry: $\Phi(-z)=1-\Phi(z)$, in particular $\Phi(0)=\tfrac12$
3. interval probabilities: $P(a\le Z\le b)=\Phi(b)-\Phi(a)$, $\ P(Z>z)=1-\Phi(z)$
4. for $X\sim N(\mu,\sigma^2)$: $\ P(a\le X\le b)=\Phi\big(\frac{b-\mu}{\sigma}\big)-\Phi\big(\frac{a-\mu}{\sigma}\big)$

**Example 2.6.1 (Ex 2.38) — normal probability by standardization**

$X\sim N(25,16)$:

$$P(20\le X\le35)=\Phi\Big(\frac{35-25}{4}\Big)-\Phi\Big(\frac{20-25}{4}\Big)=\Phi(2.5)-\Phi(-1.25)\approx0.8882$$

**Theorem 2.6.7 — mean and variance of the normal**

Using $\phi(z)=\frac1{\sqrt{2\pi}}e^{-z^2/2}$ and the substitution $x=\mu+\sigma z$ ($dx=\sigma dz$), three integrals recur:

$$\int_{-\infty}^{\infty}\phi(z)dz=1,\qquad \int_{-\infty}^{\infty}z\phi(z)dz=0\ \text{(odd integrand)},\qquad \int_{-\infty}^{\infty}z^2\phi(z)dz=1$$

Proof. The third follows by integration by parts with $u=z$, $dv=ze^{-z^2/2}dz$:

$$\int_{-\infty}^{\infty}z^2\phi(z)dz=\frac1{\sqrt{2\pi}}\Big\{\big[-ze^{-z^2/2}\big]_{-\infty}^{\infty}+\int_{-\infty}^{\infty}e^{-z^2/2}dz\Big\}=\frac1{\sqrt{2\pi}}\big(0+\sqrt{2\pi}\big)=1$$

the boundary term vanishing since $ze^{-z^2/2}\to0$, and the remaining integral being the Gaussian integral. Then

$$E(X)=\int_{-\infty}^{\infty}(\mu+\sigma z)\phi(z)dz=\mu\cdot1+\sigma\cdot0=\mu$$

$$E(X^2)=\int_{-\infty}^{\infty}(\mu^2+2\mu\sigma z+\sigma^2z^2)\phi(z)dz=\mu^2+0+\sigma^2$$

$$\mathrm{Var}(X)=E(X^2)-[E(X)]^2=(\mu^2+\sigma^2)-\mu^2=\sigma^2\qquad\blacksquare$$

**Definition 2.6.5 — bivariate normal distribution**

With parameters $\mu_X,\mu_Y,\sigma_X,\sigma_Y,\rho$ ($\sigma_X,\sigma_Y>0$, $|\rho|<1$),

$$f(x,y)=\frac1{2\pi\sigma_X\sigma_Y\sqrt{1-\rho^2}}\exp\left\{-\frac1{2(1-\rho^2)}\Big[\Big(\frac{x-\mu_X}{\sigma_X}\Big)^2-2\rho\frac{x-\mu_X}{\sigma_X}\frac{y-\mu_Y}{\sigma_Y}+\Big(\frac{y-\mu_Y}{\sigma_Y}\Big)^2\Big]\right\} \tag{2.47}$$

**Theorem 2.6.8 — marginals of the bivariate normal**

$X\sim N(\mu_X,\sigma_X^2)$ and $Y\sim N(\mu_Y,\sigma_Y^2)$: the marginals are univariate normal, with $\rho$ gone.

Proof. Standardize with $u=\dfrac{x-\mu_X}{\sigma_X}$, $v=\dfrac{y-\mu_Y}{\sigma_Y}$ ($dy=\sigma_Ydv$) and complete the square in $v$:

$$u^2-2\rho uv+v^2=(v-\rho u)^2+u^2(1-\rho^2)$$

The part free of $v$ leaves the integral:

$$f_X(x)=\frac1{2\pi\sigma_X\sqrt{1-\rho^2}}e^{-u^2/2}\int_{-\infty}^{\infty}\exp\Big\{-\frac{(v-\rho u)^2}{2(1-\rho^2)}\Big\}dv$$

The remaining integral is normal in form with variance $1-\rho^2$, hence equal to $\sqrt{2\pi}\sqrt{1-\rho^2}$ (by the Gaussian integral), and $\sqrt{1-\rho^2}$ cancels:

$$f_X(x)=\frac1{\sqrt{2\pi}\,\sigma_X}e^{-u^2/2}=\frac1{\sqrt{2\pi}\,\sigma_X}\exp\Big\{-\frac{(x-\mu_X)^2}{2\sigma_X^2}\Big\}\ \sim\ N(\mu_X,\sigma_X^2)$$

and symmetrically for $Y$. $\blacksquare$

The converse fails: normal marginals do not imply a bivariate normal joint distribution.

**Theorem 2.6.9 (Thm 2.21, 2.22) — moments and conditional distribution**

$\mathrm{Cov}(X,Y)=\rho\sigma_X\sigma_Y$, the conditional distributions are normal, and

$$E(Y\mid x)=\mu_Y+\rho\frac{\sigma_Y}{\sigma_X}(x-\mu_X),\qquad \mathrm{Var}(Y\mid x)=\sigma_Y^2(1-\rho^2) \tag{2.49}$$

In the bivariate normal, $\rho=0$ is equivalent to independence.

Proof (conditional distribution). Divide the joint density by the marginal, writing both exponents over the same denominator:

$$f(x,y)=\frac1{2\pi\sigma_X\sigma_Y\sqrt{1-\rho^2}}\exp\Big\{-\frac{(v-\rho u)^2+u^2(1-\rho^2)}{2(1-\rho^2)}\Big\}$$

$$f_X(x)=\frac1{\sqrt{2\pi}\,\sigma_X}\exp\Big\{-\frac{u^2(1-\rho^2)}{2(1-\rho^2)}\Big\}$$

The constants leave $\dfrac1{\sqrt{2\pi}\,\sigma_Y\sqrt{1-\rho^2}}$ (with $\sigma_X$ cancelling), and in the exponent the term $u^2(1-\rho^2)$ cancels exactly:

$$f_{Y\mid X}(y\mid x)=\frac1{\sqrt{2\pi}\big(\sigma_Y\sqrt{1-\rho^2}\big)}\exp\Big\{-\frac{(v-\rho u)^2}{2(1-\rho^2)}\Big\}$$

already a univariate normal. Undoing the standardization,

$$v-\rho u=\frac{y-\mu_Y}{\sigma_Y}-\rho\frac{x-\mu_X}{\sigma_X}=\frac1{\sigma_Y}\Big[y-\Big(\mu_Y+\rho\frac{\sigma_Y}{\sigma_X}(x-\mu_X)\Big)\Big]$$

$$f_{Y\mid X}(y\mid x)=\frac1{\sqrt{2\pi}\,\sigma_Y\sqrt{1-\rho^2}}\exp\left\{-\frac{\big[y-\big(\mu_Y+\rho\frac{\sigma_Y}{\sigma_X}(x-\mu_X)\big)\big]^2}{2\sigma_Y^2(1-\rho^2)}\right\}$$

and reading off the mean and variance gives (2.49). $\blacksquare$

The conditional mean is linear in $x$ (the population regression line) and the conditional variance $\sigma_Y^2(1-\rho^2)$ is smaller than $\sigma_Y^2$ regardless of $x$: knowing $X$ shrinks the uncertainty in $Y$ by the factor $1-\rho^2$.

**Theorem 2.6.10 — joint mgf of the bivariate normal**

$$M(t_1,t_2)=E\big[e^{t_1X+t_2Y}\big]=\exp\Big\{\mu_Xt_1+\mu_Yt_2+\tfrac12\big(\sigma_X^2t_1^2+2\rho\sigma_X\sigma_Yt_1t_2+\sigma_Y^2t_2^2\big)\Big\} \tag{2.48}$$

obtained by writing $X=\mu_X+\sigma_XZ_1$ and $Y=\mu_Y+\sigma_Y(\rho Z_1+\sqrt{1-\rho^2}Z_2)$ with independent standard normals and using $E[e^{aZ}]=e^{a^2/2}$.

Proof (moments from the mgf). Put $g(t_1,t_2)=\mu_Xt_1+\mu_Yt_2+\tfrac12(\sigma_X^2t_1^2+2\rho\sigma_X\sigma_Yt_1t_2+\sigma_Y^2t_2^2)$, so $M=e^g$ and

$$g_{t_1}=\mu_X+\sigma_X^2t_1+\rho\sigma_X\sigma_Yt_2,\quad g_{t_2}=\mu_Y+\sigma_Y^2t_2+\rho\sigma_X\sigma_Yt_1,\quad g_{t_1t_1}=\sigma_X^2,\ g_{t_2t_2}=\sigma_Y^2,\ g_{t_1t_2}=\rho\sigma_X\sigma_Y$$

Differentiate and set $(t_1,t_2)=(0,0)$, noting $M(0,0)=1$ (by Theorem 2.8.1).

$$E(X)=\frac{\partial M}{\partial t_1}\Big|_0=M(0,0)g_{t_1}(0,0)=\mu_X,\qquad E(Y)=\mu_Y$$

$$\frac{\partial^2M}{\partial t_1^2}=Mg_{t_1}^2+Mg_{t_1t_1}\ \Rightarrow\ E(X^2)=\mu_X^2+\sigma_X^2\ \Rightarrow\ \mathrm{Var}(X)=\sigma_X^2$$

$$\frac{\partial^2M}{\partial t_1\partial t_2}=Mg_{t_1}g_{t_2}+Mg_{t_1t_2}\ \Rightarrow\ E(XY)=\mu_X\mu_Y+\rho\sigma_X\sigma_Y\ \Rightarrow\ \mathrm{Cov}(X,Y)=\rho\sigma_X\sigma_Y$$

so $\rho=\dfrac{\mathrm{Cov}(X,Y)}{\sigma_X\sigma_Y}$ is the correlation coefficient. $\blacksquare$

Setting $t_2=0$ gives $M(t_1,0)=\exp\{\mu_Xt_1+\tfrac12\sigma_X^2t_1^2\}$, the mgf of $N(\mu_X,\sigma_X^2)$ (by uniqueness, Theorem 2.8.3); and if $\rho=0$ then (2.48) factors as $M(t_1,0)M(0,t_2)$, so $X,Y$ are independent (by Theorem 2.8.5) — the reason $\rho=0$ means independence here.

**Definition 2.6.6 (Def 2.13) — gamma function**

$$\Gamma(k)=\int_0^\infty t^{k-1}e^{-t}dt\qquad(k>0)$$

**Theorem 2.6.11 — properties of the gamma function**

1. recursion: $\Gamma(k)=(k-1)\Gamma(k-1)$
2. factorials: $\Gamma(n)=(n-1)!$ for natural $n$
3. half-integer: $\Gamma(\tfrac12)=\sqrt\pi$

Proof. **(1)** Integrate by parts with $u=t^{k-1}$, $dv=e^{-t}dt$:

$$\Gamma(k)=\big[-t^{k-1}e^{-t}\big]_0^\infty+(k-1)\int_0^\infty t^{k-2}e^{-t}dt=(k-1)\Gamma(k-1)$$

the boundary term vanishing because $e^{-t}$ decays faster than any polynomial.

**(2)** $\Gamma(1)=\int_0^\infty e^{-t}dt=1$, and iterating (1) gives $\Gamma(n)=(n-1)!\,\Gamma(1)=(n-1)!$.

**(3)** Substitute $t=\tfrac{u^2}{2}$ ($dt=u\,du$, $t^{-1/2}=\sqrt2/u$):

$$\Gamma(\tfrac12)=\int_0^\infty\frac{\sqrt2}{u}e^{-u^2/2}u\,du=\sqrt2\int_0^\infty e^{-u^2/2}du=\sqrt2\cdot\frac{\sqrt{2\pi}}{2}=\sqrt\pi$$

(by the Gaussian integral). $\blacksquare$

The gamma function extends the factorial to real arguments, $\Gamma(n+1)=n!$.

**Definition 2.6.7 — gamma distribution**

$X\sim\mathrm{GAM}(k,\theta)$ with shape $k$ and scale $\theta$:

$$f(x)=\frac1{\theta^k\Gamma(k)}x^{k-1}e^{-x/\theta},\ x>0 \tag{2.50}$$

$$E(X)=k\theta,\qquad \mathrm{Var}(X)=k\theta^2$$

$k=1$ gives the exponential; integer $k$ gives the sum of $k$ independent exponentials (Erlang).

**Theorem 2.6.12 — sum of independent exponentials**

If $X_1,\dots,X_k$ are independent $\mathrm{EXP}(\theta)$ then $S_k=X_1+\cdots+X_k\sim\mathrm{GAM}(k,\theta)$.

Proof (induction on $k$, by convolution, Theorem 2.7.6). For $k=1$, the density $\tfrac1\theta e^{-x/\theta}=\dfrac1{\theta^1\Gamma(1)}x^{1-1}e^{-x/\theta}$ is that of $\mathrm{GAM}(1,\theta)$.

Assume $f_{S_k}(s)=\dfrac1{\theta^k\Gamma(k)}s^{k-1}e^{-s/\theta}$ and add an independent $X_{k+1}\sim\mathrm{EXP}(\theta)$. Both densities vanish off $(0,\infty)$, so the convolution runs over $0<s<z$:

$$f_{S_{k+1}}(z)=\int_0^z f_{S_k}(s)f_{X_{k+1}}(z-s)ds=\int_0^z\frac{s^{k-1}e^{-s/\theta}}{\theta^k\Gamma(k)}\cdot\frac{e^{-(z-s)/\theta}}{\theta}ds=\frac{e^{-z/\theta}}{\theta^{k+1}\Gamma(k)}\int_0^z s^{k-1}ds=\frac{e^{-z/\theta}}{\theta^{k+1}\Gamma(k)}\cdot\frac{z^k}{k}$$

the key step being $e^{-s/\theta}e^{-(z-s)/\theta}=e^{-z/\theta}$, free of $s$. With $k\Gamma(k)=\Gamma(k+1)$,

$$f_{S_{k+1}}(z)=\frac1{\theta^{k+1}\Gamma(k+1)}z^{(k+1)-1}e^{-z/\theta}$$

the $\mathrm{GAM}(k+1,\theta)$ density. $\blacksquare$

The moments follow from linearity: $E(S_k)=k\theta$ and $\mathrm{Var}(S_k)=k\theta^2$.

**Definition 2.6.8 (Def 2.14) — beta function**

$$B(a,b)=\int_0^1 x^{a-1}(1-x)^{b-1}dx=\frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}$$

**Definition 2.6.9 — beta distribution**

$X\sim\mathrm{BETA}(a,b)$:

$$f(x)=\frac1{B(a,b)}x^{a-1}(1-x)^{b-1},\ 0<x<1$$

$$E(X)=\frac{a}{a+b},\qquad \mathrm{Var}(X)=\frac{ab}{(a+b)^2(a+b+1)}$$

Supported on $[0,1]$, it models proportions and probabilities; $a=b=1$ gives $U(0,1)$.

**Definition 2.6.10 — other continuous distributions**

- Weibull: $f(x)=\dfrac{\beta}{\theta^\beta}x^{\beta-1}e^{-(x/\theta)^\beta}$ (lifetimes, reliability; $\beta=2$ is Rayleigh)
- Pareto: $f(x)=\dfrac k\theta(1+x/\theta)^{-(k+1)}$ (heavy tail)
- Cauchy: $f(x)=\dfrac1{\pi\beta[1+((x-\alpha)/\beta)^2]}$ (no mean)

together with the lognormal, double exponential and logistic distributions.

## 2.7 Functions of Random Variables

### 2.7.1 The cdf method

**Theorem 2.7.1 — discrete transformation**

The event $\{Y=y\}$ for $Y=g(X)$ is the event that $X$ lies in the preimage $g^{-1}(y)=\{x:g(x)=y\}$, and the events $\{X=x\}$ for those $x$ are mutually exclusive, so

$$f_Y(y)=P\Big(\bigcup_{x\in g^{-1}(y)}\{X=x\}\Big)=\sum_{\{x:\,g(x)=y\}}f_X(x)\quad\text{(by countable additivity)}$$

For example $Y=X^2$ has $g^{-1}(4)=\{-2,2\}$, so $f_Y(4)=f_X(-2)+f_X(2)$. In the continuous case one instead writes $F_Y(y)=P(g(X)\le y)$ as a probability for $X$ and differentiates.

**Theorem 2.7.2 (Thm 2.23) — monotone transformation**

If $g$ is one-to-one (monotone) and $X\sim f_X$, then $Y=g(X)$ has density

$$f_Y(y)=f_X\big(g^{-1}(y)\big)\Big|\frac{d}{dy}g^{-1}(y)\Big|$$

Proof. For increasing $g$, $F_Y(y)=P(g(X)\le y)=P(X\le g^{-1}(y))=F_X(g^{-1}(y))$, and differentiating gives $f_Y(y)=f_X(g^{-1}(y))\frac{d}{dy}g^{-1}(y)$. For decreasing $g$, $F_Y(y)=1-F_X(g^{-1}(y))$ and the derivative changes sign; the two cases combine into the absolute value. $\blacksquare$

**Theorem 2.7.3 — probability integral transform**

If the cdf $F$ of a continuous random variable is continuous and strictly increasing, then feeding $X$ into its own cdf gives

$$U=F(X)\ \sim\ U(0,1)$$

Proof. For $0\le u\le1$, monotonicity gives $F(X)\le u\iff X\le F^{-1}(u)$, so

$$P(U\le u)=P\big(X\le F^{-1}(u)\big)=F\big(F^{-1}(u)\big)=u$$

i.e. $F_U(u)=u$ on $[0,1]$. $\blacksquare$

Conversely, if $U\sim U(0,1)$ then $X=F^{-1}(U)$ has cdf $F$, since $P(X\le x)=P(U\le F(x))=F(x)$ — the inverse transform method for generating random numbers from any distribution.

**Example 2.7.1 (Ex 2.40) — $Y=X^2$, a many-to-one transformation**

$g(x)=x^2$ is not one-to-one: $g^{-1}(y)=\{-\sqrt y,\sqrt y\}$ for $y>0$, so Theorem 2.7.2 does not apply and we go through the cdf:

$$F_Y(y)=P(X^2\le y)=P(-\sqrt y\le X\le\sqrt y)=F_X(\sqrt y)-F_X(-\sqrt y)$$

Differentiating (chain rule, $\frac{d}{dy}\sqrt y=\frac1{2\sqrt y}$), the two preimages contribute additively:

$$f_Y(y)=f_X(\sqrt y)\frac1{2\sqrt y}-f_X(-\sqrt y)\Big(-\frac1{2\sqrt y}\Big)=\frac1{2\sqrt y}\big[f_X(\sqrt y)+f_X(-\sqrt y)\big],\quad y>0$$

the continuous analogue of summing over the preimage in Theorem 2.7.1. For $X\sim N(0,1)$ the density is even, so

$$f_Y(y)=\frac1{\sqrt y}\cdot\frac1{\sqrt{2\pi}}e^{-y/2}=\frac1{\sqrt{2\pi}}y^{-1/2}e^{-y/2},\quad y>0$$

which is the $\mathrm{GAM}(\tfrac12,2)$ density, i.e. $\chi^2(1)$: the square of a standard normal is chi-squared with one degree of freedom.

**Example 2.7.2 (Ex 2.42) — sum of two independent uniforms**

$X,Y\sim U(0,1)$ independent, $Z=X+Y\in[0,2]$, by convolution (Theorem 2.7.6):

$$f_Z(z)=\int_{-\infty}^{\infty}f_X(x)f_Y(z-x)dx$$

The integrand is $1$ only when $0\le x\le1$ and $0\le z-x\le1$, i.e. $x\in[\max(0,z-1),\min(1,z)]$, so the range of integration depends on $z$:

$$0\le z\le1:\ f_Z(z)=\int_0^z 1\,dx=z,\qquad 1<z\le2:\ f_Z(z)=\int_{z-1}^1 1\,dx=2-z$$

$$f_Z(z)=\begin{cases}z,&0<z\le1\\ 2-z,&1<z\le2\\ 0,&\text{otherwise}\end{cases}$$

a triangular density peaking at $z=1$.

**Method 2.7.1 — integration limits in a convolution**

The range of integration is the overlap of the two supports: in $f_Z(z)=\int f_X(x)f_Y(z-x)dx$ the integrand is non-zero only where $x$ and $z-x$ both lie in their supports, which for $X,Y\in[0,1]$ means

$$x\in\big[\max(0,z-1),\ \min(1,z)\big]$$

The density changes form at the values of $z$ where these expressions switch. To find them mechanically: substitute every combination of the endpoints of the supports (the corners of the support box) into the sum, then deduplicate and sort. For $Z=X+Y$ with $X,Y\in[0,1]$ the corner sums $\{0,1,1,2\}$ give breakpoints $0,1,2$.

**Example 2.7.3 — breakpoints for a trivariate sum**

$X\in[0,1]$, $Y\in[0,2]$, $Z\in[0,3]$ and $W=X+Y+Z$. The $2^3=8$ corners give

$$W\in\{0,1,2,3,3,4,5,6\}$$

so after deduplicating and sorting the breakpoints are $0,1,2,3,4,5,6$ and $f_W$ has a different form on each of the six intervals. The routine is the same for $n$ variables: substitute the $2^n$ corners, deduplicate, sort, then integrate over the overlap on each interval.

### 2.7.2 Joint transformations

**Theorem 2.7.4 (Thm 2.24) — discrete joint transformation**

For a discrete vector $Y=g(X)$,

$$f_{Y}(y)=\sum_{\{x:\,g(x)=y\}}f_{X}(x)$$

and if $g$ is one-to-one the sum has a single term:

$$f_{Y}(y)=f_{X}\big(g^{-1}(y)\big)$$

with no Jacobian factor, unlike the continuous case.

Proof. The event $\{Y=y\}$ is the event that $X$ takes one of the values in the preimage $g^{-1}(y)$, and those events are mutually exclusive, so the probabilities add (by countable additivity):

$$f_{Y}(y)=P\Big(\bigcup_{x\in g^{-1}(y)}\{X=x\}\Big)=\sum_{\{x:\,g(x)=y\}}f_{X}(x)$$

A one-to-one $g$ leaves one term. In the discrete case probability is mass attached to points, and a one-to-one map merely relabels the points; no stretching factor $|J|$ arises. $\blacksquare$

**Theorem 2.7.5 (Thm 2.25) — continuous joint transformation, Jacobian**

If $Y=g(X)$ is one-to-one and the inverse $x_i=h_i(y)$ has continuous partial derivatives,

$$f_{Y}(y)=f_{X}\big(h_1(y),\dots,h_n(y)\big)|J|,\qquad J=\det\Big(\frac{\partial(x_1,\dots,x_n)}{\partial(y_1,\dots,y_n)}\Big)$$

**Example 2.7.4 (Ex 2.45) — sum and difference of independent exponentials**

$X_1,X_2\sim\mathrm{EXP}(1)$ independent, $Y_1=X_1+X_2$, $Y_2=X_1-X_2$. Independence gives

$$f_{X_1,X_2}(x_1,x_2)=e^{-x_1}e^{-x_2}=e^{-(x_1+x_2)},\qquad x_1>0,\ x_2>0$$

(1) Inverse transformation:

$$x_1=\frac{y_1+y_2}{2},\qquad x_2=\frac{y_1-y_2}{2}$$

(2) Support: $x_1>0$ and $x_2>0$ become $y_1+y_2>0$ and $y_1-y_2>0$, i.e. $y_1>|y_2|$. The support must be carried over too, since it fixes the later limits of integration.

(3) Jacobian:

$$J=\det\begin{pmatrix}\tfrac12&\tfrac12\\ \tfrac12&-\tfrac12\end{pmatrix}=-\tfrac12\ \Rightarrow\ |J|=\tfrac12$$

(4) Joint density; the exponent is $x_1+x_2=y_1$, so $y_2$ disappears:

$$f_{Y_1,Y_2}(y_1,y_2)=e^{-y_1}\cdot\tfrac12=\tfrac12 e^{-y_1},\qquad y_1>|y_2|$$

(5) Marginals, minding the support $y_1>|y_2|$:

$$f_{Y_1}(y_1)=\int_{-y_1}^{y_1}\tfrac12 e^{-y_1}dy_2=y_1e^{-y_1},\quad y_1>0$$

which is $\mathrm{GAM}(2,1)$, agreeing with Theorem 2.6.12, and

$$f_{Y_2}(y_2)=\int_{|y_2|}^{\infty}\tfrac12 e^{-y_1}dy_1=\tfrac12 e^{-|y_2|},\quad-\infty<y_2<\infty$$

the double exponential (Laplace) density. The product of the marginals differs from the joint density — the support is a wedge $y_1>|y_2|$, not a rectangle — so $Y_1,Y_2$ are not independent.

**Example 2.7.5 (Ex 2.46) — ratio of two standard normals is Cauchy**

$X_1,X_2\sim N(0,1)$ independent, $Y_2=X_1/X_2$. A single ratio is not a two-dimensional one-to-one map, so add the auxiliary variable $Y_1=X_2$.

(1) Inverse: $x_2=y_1$, $x_1=y_1y_2$.

(2) Jacobian:

$$J=\det\begin{bmatrix}y_2&y_1\\ 1&0\end{bmatrix}=-y_1,\qquad|J|=|y_1|$$

(3) Joint density from $f_{X_1,X_2}(x_1,x_2)=\dfrac1{2\pi}e^{-(x_1^2+x_2^2)/2}$ (by Theorem 2.7.5):

$$f_{Y_1,Y_2}(y_1,y_2)=\frac1{2\pi}e^{-(y_1^2y_2^2+y_1^2)/2}|y_1|=\frac{|y_1|}{2\pi}e^{-y_1^2(1+y_2^2)/2}$$

(4) Integrate out the auxiliary variable; the integrand is even in $y_1$, and with $a=\dfrac{1+y_2^2}{2}$ we use $\int_0^\infty y_1e^{-ay_1^2}dy_1=\dfrac1{2a}$:

$$f_{Y_2}(y_2)=\frac{2}{2\pi}\cdot\frac1{2a}=\frac1{\pi(1+y_2^2)},\qquad-\infty<y_2<\infty$$

the Cauchy density — heavy-tailed, with no mean.

**Theorem 2.7.6 (Thm 2.26) — convolution**

If $X,Y$ are independent, $Z=X+Y$ has

$$f_Z(z)=\begin{cases}\sum_{x}f_X(x)f_Y(z-x)&\text{discrete}\\ \int_{-\infty}^{\infty}f_X(x)f_Y(z-x)dx&\text{continuous}\end{cases}$$

Proof (continuous).

$$F_Z(z)=P(X+Y\le z)=\iint_{x+y\le z}f_X(x)f_Y(y)\,dy\,dx=\int_{-\infty}^{\infty}f_X(x)\Big(\int_{-\infty}^{z-x}f_Y(y)dy\Big)dx$$

Differentiating with respect to $z$ turns the inner integral into $f_Y(z-x)$. $\blacksquare$

Proof (via conditional expectation). Since $E[I_A]=P(A)$ (by Theorem 2.5.2), probabilities become expectations and the law of total expectation applies:

$$F_Z(z)=P(X+Y\le z)=E\big[I(X+Y\le z)\big]=E\Big[E\big[I(X+Y\le z)\mid X=x\big]\Big]=E\big[P(Y\le z-x)\big]$$

$$=E\big[F_Y(z-x)\big]=\int_{-\infty}^{\infty}F_Y(z-x)f_X(x)dx$$

(by Theorem 2.4.6 and independence). Differentiating in $z$ gives $f_Z(z)=\int f_Y(z-x)f_X(x)dx$. $\blacksquare$

**Example 2.7.6 (Ex 2.47) — $Z=X+Y$ by two routes**

$X,Y\sim U(0,1)$ independent.

Route 1, joint transformation. Add the auxiliary $U=X$, so $(x,y)\mapsto(z,u)=(x+y,x)$ with inverse $x=u$, $y=z-u$ and

$$J=\det\begin{bmatrix}0&1\\ 1&-1\end{bmatrix}=-1,\qquad|J|=1$$

Since $f_{X,Y}=1$ on the unit square, $f_{Z,U}(z,u)=1$ for $0\le u\le1$, $0\le z-u\le1$; integrating $u$ over $[\max(0,z-1),\min(1,z)]$,

$$f_Z(z)=\begin{cases}\int_0^z du=z,&0<z\le1\\ \int_{z-1}^1 du=2-z,&1<z\le2\\ 0,&\text{otherwise}\end{cases}$$

Route 2, convolution with indicators. Putting $f_X=f_Y=I_{(0,1)}$ into Theorem 2.7.6, and using $I_{(0,1)}(z-x)=1\iff z-1<x<z$,

$$f_Z(z)=\int_{-\infty}^{\infty}I_{(0,1)}(x)I_{(0,1)}(z-x)dx=z\,I_{(0,1)}(z)+(2-z)I_{[1,2)}(z)$$

the same triangular density. Convolution is the special case of a joint transformation in which the auxiliary variable $U=X$ is integrated out, with $|J|=1$.

## 2.8 Moment Generating Function

**Definition 2.8.1 (Def 2.15) — moments**

- $r$-th moment: $\mu_r'=E(X^r)$
- $r$-th central moment: $\mu_r=E[(X-\mu)^r]$
- $\mu_1'=\mu$, the mean
- $\mu_2=\sigma^2$, the variance
- $\mu_3/\sigma^3$, the skewness
- $\mu_4/\sigma^4$, the kurtosis

**Definition 2.8.2 (Def 2.16) — moment generating function**

If the expectation exists for $-h<t<h$,

$$M_X(t)=E\big(e^{tX}\big)$$

**Theorem 2.8.1 (Thm 2.27) — moments from the mgf**

$$M_X^{(r)}(0)=E(X^r)$$

Proof. By the exponential series,

$$M_X(t)=E(e^{tX})=E\Big(\sum_{r\ge0}\frac{(tX)^r}{r!}\Big)=\sum_{r\ge0}\frac{t^r}{r!}E(X^r)$$

Differentiating $r$ times and setting $t=0$ leaves only the term $\frac{t^r}{r!}E(X^r)$, so $M_X^{(r)}(0)=E(X^r)$; e.g. $M'(0)=E(X)$ and $M''(0)=E(X^2)$. $\blacksquare$

**Theorem 2.8.2 — mgf of the standard distributions**

- $B(n,p)$: $(q+pe^t)^n$
- $\mathrm{POI}(\lambda)$: $e^{\lambda(e^t-1)}$
- $\mathrm{EXP}(\theta)$: $(1-\theta t)^{-1}$, $t<\tfrac1\theta$
- $\mathrm{GAM}(k,\theta)$: $(1-\theta t)^{-k}$, $t<\tfrac1\theta$
- $N(\mu,\sigma^2)$: $\exp\big(\mu t+\tfrac12\sigma^2t^2\big)$

Proof. Substitute each distribution into $M_X(t)=E(e^{tX})$.

$$B(n,p):\ \sum_{x=0}^{n}e^{tx}\binom nx p^xq^{n-x}=\sum_x\binom nx(pe^t)^xq^{n-x}=(q+pe^t)^n\quad\text{(by the binomial theorem)}$$

$$\mathrm{POI}(\lambda):\ \sum_{x\ge0}e^{tx}\frac{e^{-\lambda}\lambda^x}{x!}=e^{-\lambda}\sum_x\frac{(\lambda e^t)^x}{x!}=e^{-\lambda}e^{\lambda e^t}=e^{\lambda(e^t-1)}\quad\text{(by the exponential series)}$$

$$\mathrm{EXP}(\theta):\ \int_0^\infty e^{tx}\frac1\theta e^{-x/\theta}dx=\frac1\theta\cdot\frac1{1/\theta-t}=(1-\theta t)^{-1}$$

For $\mathrm{GAM}(k,\theta)$, collecting the exponent as $-(1/\theta-t)x$ and using the gamma integral gives $(1-\theta t)^{-k}$. For the normal, completing the square in the exponent $tx-\frac{(x-\mu)^2}{2\sigma^2}$ gives $-\frac{(x-\mu-\sigma^2t)^2}{2\sigma^2}+\mu t+\frac12\sigma^2t^2$, and the first part integrates to $1$. $\blacksquare$

Proof (normal, by standardization). For $Z\sim N(0,1)$, completing the square $z^2-2tz=(z-t)^2-t^2$,

$$M_Z(t)=\frac1{\sqrt{2\pi}}\int_{-\infty}^{\infty}e^{tz-z^2/2}dz=e^{t^2/2}\cdot\frac1{\sqrt{2\pi}}\int_{-\infty}^{\infty}e^{-(z-t)^2/2}dz=e^{t^2/2}$$

the remaining integral being the total mass of $N(t,1)$. Then with $X=\mu+\sigma Z$ (by Theorem 2.8.4(1)),

$$M_X(t)=e^{\mu t}M_Z(\sigma t)=e^{\mu t}e^{(\sigma t)^2/2}=\exp\big(\mu t+\tfrac12\sigma^2t^2\big)\qquad\blacksquare$$

**Example 2.8.1 (Ex 2.52) — mean and variance from the mgf**

$B(n,p)$: $M(t)=(q+pe^t)^n$ and $M'(t)=n(q+pe^t)^{n-1}pe^t$, so $M'(0)=np=E(X)$; computing $E(X^2)$ from $M''(0)$ gives $\mathrm{Var}(X)=np(1-p)$.

**Theorem 2.8.3 (Thm 2.28) — uniqueness**

If two random variables have the same mgf on an interval, they have the same distribution: the mgf determines the distribution uniquely.

This gives the standard identification argument: compute an mgf, match it to a known one, and conclude the distribution — the main tool for finding distributions of sums and transformations.

**Theorem 2.8.4 (Thm 2.29) — linear transformation and independent sums**

1. $Y=aX+b\ \Rightarrow\ M_Y(t)=e^{bt}M_X(at)$
2. $X_1,\dots,X_n$ independent, $Y=\sum_i X_i\ \Rightarrow\ M_Y(t)=\prod_{i=1}^{n}M_{X_i}(t)$

Proof. **(1)** $M_Y(t)=E[e^{t(aX+b)}]=e^{bt}E[e^{(at)X}]=e^{bt}M_X(at)$.

**(2)** By independence the expectation of the product separates (by Theorem 2.4.2):

$$M_Y(t)=E\big[e^{t\sum X_i}\big]=E\Big[\prod_i e^{tX_i}\Big]=\prod_i E[e^{tX_i}]=\prod_i M_{X_i}(t)\qquad\blacksquare$$

**Example 2.8.2 (Ex 2.56–2.59) — distributions of independent sums**

The mgf of a sum is a product, so uniqueness identifies the distribution at once:

- $n$ independent $\mathrm{Ber}(p)$: $(q+pe^t)^n\ \Rightarrow\ B(n,p)$
- independent $\mathrm{POI}(\lambda_i)$: $e^{(\sum\lambda_i)(e^t-1)}\ \Rightarrow\ \mathrm{POI}(\sum\lambda_i)$
- $n$ independent $\mathrm{EXP}(\theta)$: $(1-\theta t)^{-n}\ \Rightarrow\ \mathrm{GAM}(n,\theta)$
- independent normals, $\sum a_iX_i$ with $X_i\sim N(\mu_i,\sigma_i^2)$: $\ N\big(\sum a_i\mu_i,\ \sum a_i^2\sigma_i^2\big)$

**Definition 2.8.3 (Def 2.17) — joint mgf**

$$M_{X}(t_1,\dots,t_k)=E\Big[\exp\Big(\sum_{i=1}^{k}t_iX_i\Big)\Big]$$

Setting all other arguments to $0$ gives the marginal mgf:

$$M_{X_i}(t_i)=M_{X}(0,\dots,0,t_i,0,\dots,0) \tag{2.62}$$

**Theorem 2.8.5 (Thm 2.30) — independence and the joint mgf**

$X,Y$ are independent if and only if the joint mgf factors:

$$M_{X,Y}(t_1,t_2)=M_X(t_1)M_Y(t_2)$$

## 2.9 Probability Generating Function

**Definition 2.9.1 (Def 2.18) — probability generating function**

For $X$ taking values $0,1,2,\dots$ with $p_j=P(X=j)$, and $|s|\le1$,

$$G_X(s)=E\big(s^X\big)=\sum_{j=0}^{\infty}p_js^j$$

the mgf with $e^t$ replaced by $s$. Its main properties:

1. mean: $E(X)=G_X'(1)$, and in general the factorial moments $E[X(X-1)\cdots(X-k+1)]=G_X^{(k)}(1)$
2. independent sums: $X,Y$ independent and $Z=X+Y$ give $G_Z(s)=G_X(s)G_Y(s)$
3. recovery of probabilities: $P(X=k)=\dfrac1{k!}G_X^{(k)}(0)$

**Theorem 2.9.1 — recovery of the probabilities**

Proof. Expand the pgf as a series:

$$G_X(s)=p_0+p_1s+p_2s^2+p_3s^3+p_4s^4+\cdots$$

Each differentiation brings the exponent down as a factor:

$$G_X'(s)=p_1+2p_2s+3p_3s^2+\cdots,\quad G_X''(s)=2p_2+6p_3s+\cdots,\quad G_X'''(s)=6p_3+24p_4s+\cdots$$

The constant term of the third derivative is not $p_3$ alone but $3!\,p_3$, the exponents $3,2,1$ having come down in turn. Setting $s=0$ kills every term containing $s$:

$$G_X'''(0)=3!\,p_3\ \Longrightarrow\ p_3=\frac1{3!}G_X'''(0)$$

and after $k$ differentiations the coefficient of $s^k$ has picked up $k(k-1)\cdots1=k!$:

$$G_X^{(k)}(0)=k!\,p_k\ \Longrightarrow\ p_k=\frac1{k!}G_X^{(k)}(0)$$

The factor $1/k!$ removes exactly the $k!$ generated by differentiating. $\blacksquare$

**Example 2.9.1 (Ex 2.61) — pgf of the binomial and the Poisson**

- binomial $B(n,p)$: one Bernoulli trial has pgf $q+ps$, so by property 2 $\ G_X(s)=(q+ps)^n$, and $G_X'(1)=np=E(X)$
- Poisson $\mathrm{POI}(\lambda)$: $\ G_X(s)=\sum_{j\ge0}\dfrac{e^{-\lambda}\lambda^j}{j!}s^j=e^{-\lambda}e^{\lambda s}=e^{\lambda(s-1)}$, and $G_X'(1)=\lambda=E(X)$

For several variables one uses the joint pgf $G_{X_1,\dots,X_r}(s_1,\dots,s_r)=E\big[\prod_i s_i^{X_i}\big]$.
