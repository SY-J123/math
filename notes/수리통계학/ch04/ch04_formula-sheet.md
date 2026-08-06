# Chapter 4. Estimation

## 4.1 Method of Moments

**Definition 4.1.1 (Def 4.1) — statistic, estimator, estimate**

- **statistic**: a function $T(X)=T(X_1,\dots,X_n)$ of the random sample containing no unknown parameter
- **estimator**: a statistic $T(X)$ (a random variable) used to estimate $g(\theta)$
- **estimate**: the value $T(x)=T(x_1,\dots,x_n)$ (a number) obtained by substituting the observations $X_i=x_i$

An estimator of $\theta$ is written $\hat\theta$ or $\hat\theta_n$, the subscript emphasizing that the sample size used is $n$.

**Example 4.1.1 (Ex 4.1) — estimators of the normal mean and variance**

For a random sample from $N(\mu,\sigma^2)$ the sample mean and sample variance may be used:

$$\hat\mu=T_1(X)=\overline{X}_n=\frac1n\sum_{i=1}^n X_i,\qquad
\hat\sigma^2=T_2(X)=S_n^2=\frac{1}{n-1}\sum_{i=1}^n(X_i-\overline{X}_n)^2$$

**Definition 4.1.2 — population moment, sample moment**

- $r$-th **population moment**: $\ \mu_r'=E(X^r)$
- $r$-th **sample moment**: $\ m_r'=\dfrac1n\sum_{i=1}^n X_i^{\,r}$

**Method 4.1.1 — moment estimator**

Each population moment is a function of the parameters, $\mu_j'=\mu_j'(\theta_1,\dots,\theta_k)$. Equating the $k$ lowest-order moments to the corresponding sample moments,

$$m_j'=\mu_j'(\theta_1,\theta_2,\dots,\theta_k),\qquad j=1,2,\dots,k$$

and solving gives the **moment estimator** $(\hat\theta_1,\dots,\hat\theta_k)$.

The justification is the law of large numbers: $m_r'\xrightarrow{\ p\ }\mu_r'$, so replacing population moments by sample moments and solving yields an estimate close to the true value. The method is comparatively simple, but the solution need not be unique.

**Example 4.1.2 (Ex 4.2) — normal $N(\mu,\sigma^2)$**

$\mu_1'=\mu$ and $\mu_2'=E(X^2)=\sigma^2+\mu^2$, so the moment equations are $m_1'=\mu$ and $m_2'=\sigma^2+\mu^2$:

$$\hat\mu=\frac1n\sum_{i=1}^n x_i=\overline{X}_n,\qquad
\hat\sigma^2=m_2'-\hat\mu^2=\frac1n\sum x_i^2-\overline{X}_n^2=\frac1n\sum_{i=1}^n(x_i-\overline{X}_n)^2$$

Note that this $\hat\sigma^2$ divides by $n$, not by $n-1$ as $S_n^2$ does; the moment estimator is therefore **biased**, with $E(\hat\sigma^2)=\frac{n-1}{n}\sigma^2$.

**Example 4.1.3 (Ex 4.3) — gamma $\mathrm{GAM}(k,\theta)$**

$\mu_1'=k\theta$ and $\mu_2'=k\theta^2+k^2\theta^2$ (since $\operatorname{Var}(X)=k\theta^2$). From $m_2'-m_1'^2=k\theta^2$ together with $\theta=\dfrac{k\theta^2}{k\theta}$ and $k=\dfrac{(k\theta)^2}{k\theta^2}$,

$$\hat\theta=\frac{\frac1n\sum x_i^2-\overline{X}_n^2}{\overline{X}_n}=\frac{\sum_{i=1}^n(x_i-\overline{X}_n)^2}{n\,\overline{X}_n},\qquad
\hat k=\frac{\overline{X}_n}{\hat\theta}=\frac{n\,\overline{X}_n^2}{\sum_{i=1}^n(x_i-\overline{X}_n)^2}$$

## 4.2 Maximum Likelihood Estimation

Maximum likelihood chooses the parameter that makes the observed data most plausible. If a Bernoulli success probability is known to be either $0.1$ or $0.9$ and $10$ observations sum to $8$, everyone picks $0.9$ — because a sum of $8$ is more probable when $0.9$ is true.

**Definition 4.2.1 — likelihood function**

The joint density $f(x_1,\dots,x_n;\theta)$ is normally read as a function of $(x_1,\dots,x_n)$ for fixed $\theta$. Read instead as a function of $\theta\in\Omega$ for the observed $X_1=x_1,\dots,X_n=x_n$, it is the **likelihood function**

$$L(\theta)=L(\theta\,;x_1,\dots,x_n)=f(x_1,\dots,x_n;\theta)$$

For one fixed $\theta$ and one fixed $(x_1,\dots,x_n)$ the density and the likelihood take the same value, but they are different functions — just as $y=x^2+z$ may be read as $f(x)=x^2+z$ or as $g(z)=x^2+z$. In particular the likelihood is a function of $\theta$ and is **not** a density.

**Theorem 4.2.1 — likelihood of an independent sample**

If $X_1,\dots,X_n$ are independent with densities $f_i(x_i;\theta)$, the joint density is a product, so

$$L(\theta\,;x_1,\dots,x_n)=\prod_{i=1}^n f_i(x_i;\theta)$$

and for a random sample from $f(x;\theta)$, where $f_i=f$ for all $i$,

$$L(\theta\,;x_1,\dots,x_n)=\prod_{i=1}^n f(x_i;\theta) \tag{4.3}$$

**Example 4.2.1 — how a likelihood is formed (three coin tosses)**

Let $\theta$ be the probability of heads, and let the outcome be head, tail, head.

1. Joint density, both data and parameter unknown:
   $$f(x_1,x_2,x_3;\theta)=\theta^{x_1}(1-\theta)^{1-x_1}\cdot\theta^{x_2}(1-\theta)^{1-x_2}\cdot\theta^{x_3}(1-\theta)^{1-x_3}$$
2. Substituting the observations $x_1=1,\ x_2=0,\ x_3=1$ gives $\theta$, $1-\theta$, $\theta$ respectively.
3. The $x$'s are gone and only $\theta$ remains: $\ L(\theta)=\theta\cdot(1-\theta)\cdot\theta=\theta^2(1-\theta)$.
4. $\log L(\theta)=2\log\theta+\log(1-\theta)$, and
   $$\frac{d}{d\theta}\log L(\theta)=\frac2\theta-\frac1{1-\theta}=0\ \Longrightarrow\ 2(1-\theta)=\theta\ \Longrightarrow\ \hat\theta=\frac23$$
   with $\frac{d^2}{d\theta^2}\log L=-\frac2{\theta^2}-\frac1{(1-\theta)^2}<0$, so this is the maximum.

The answer $\hat\theta=2/3$ is exactly the sample proportion of heads.

**Definition 4.2.2 (Def 4.2) — maximum likelihood estimator**

If $\hat\theta=\hat\theta(x_1,\dots,x_n)\in\Omega$ **maximizes** $L(\theta\,;x_1,\dots,x_n)$, then $\hat\theta=\hat\theta(X_1,\dots,X_n)$ is the **maximum likelihood estimator** (MLE) of $\theta$ — the parameter value that makes the data actually observed most probable.

**Example 4.2.2 (Ex 4.4) — reading a likelihood table, $B(5,p)$**

For $Y=\sum_{i=1}^5 X_i\sim B(5,p)$, tabulate $f_Y(y;p)$ over a grid of $p$. The likelihood is a function of $p$ for fixed $Y=y$, so the table is read **along a row**: for $Y=3$ the row is largest at $p=0.6$ ($0.34560$), so $\hat p=0.6$; for $Y=1$ it is largest at $p=0.2$. Larger $Y$ giving larger $\hat p$ is what one would expect of a success probability.

**Method 4.2.1 — maximizing the log-likelihood**

Maximizing $L$ is equivalent to maximizing the **log-likelihood**

$$\log L(\theta\,;x_1,\dots,x_n)=\log\prod_{i=1}^n f_i(x_i;\theta)=\sum_{i=1}^n\log f_i(x_i;\theta) \tag{4.4}$$

which in many cases reduces to solving $\dfrac{d}{d\theta}\log L(\theta\,;x_1,\dots,x_n)=0$.

Two reasons the logarithm helps:

- $\log$ is strictly increasing, so it does not move the maximizing $\theta$
- the likelihood of a random sample is a **product**, and taking logs turns it into a **sum**, which is far easier to differentiate

A root of $\frac{d}{d\theta}\log L=0$ need not be a maximum; this must be checked, e.g. by the second derivative. Here $\log$ always means the natural logarithm.

**Example 4.2.3 (Ex 4.5) — exponential $\mathrm{EXP}(\theta)$**

For $x_i>0$,

$$L(\theta)=\Big(\frac1\theta\Big)^{\!n}\exp\Big(-\sum_{i=1}^n x_i/\theta\Big),\qquad
\log L(\theta)=-n\log\theta-\sum_{i=1}^n x_i/\theta$$

$$\frac{d}{d\theta}\log L(\theta)=-\frac n\theta+\frac{\sum_{i=1}^n x_i}{\theta^2}=0\ \Longrightarrow\ \hat\theta=\overline{X}_n$$

and the second derivative at $\bar x_n$ is negative, so $\overline{X}_n$ is the MLE.

**Example 4.2.4 (Ex 4.6) — Poisson $\mathrm{POI}(\lambda)$**

$$L(\lambda)=\frac{\exp(-n\lambda)\lambda^{\sum x_i}}{\prod_{i=1}^n x_i!},\qquad
\log L(\lambda)=-n\lambda+\sum_{i=1}^n x_i\log\lambda-\log\Big(\prod_{i=1}^n x_i!\Big)$$

$$\frac{d}{d\lambda}\log L(\lambda)=-n+\sum_{i=1}^n\frac{x_i}{\lambda}=0\ \Longrightarrow\ \hat\lambda=\bar x_n$$

$$\frac{d^2}{d\lambda^2}\log L(\lambda)\Big|_{\lambda=\bar x_n}=-\frac{\sum x_i}{\lambda^2}\Big|_{\lambda=\bar x_n}=-\frac{n}{\bar x_n}<0$$

so the MLE is $\overline{X}_n$, assuming $\bar x_n>0$ (for $\bar x_n=0$ see Example 4.2.8).

Sometimes the second derivative need not be checked: here $\log L\to-\infty$ as $\lambda\to0$ and as $\lambda\to\infty$, and there is exactly one critical point, so that point must be the maximum.

**Method 4.2.2 — vector parameter**

For $\theta=(\theta_1,\dots,\theta_k)$, solve the $k$ simultaneous equations

$$\frac{\partial}{\partial\theta_i}\log L(\theta_1,\dots,\theta_k)=0,\qquad i=1,\dots,k$$

and check that the solution maximizes the log-likelihood.

**Example 4.2.5 (Ex 4.7) — normal $N(\mu,\sigma^2)$**

$$L(\mu,\sigma^2)=\frac{1}{(2\pi\sigma^2)^{n/2}}\exp\Big[-\sum_{i=1}^n(x_i-\mu)^2/2\sigma^2\Big],\qquad
\log L=-\frac n2\log(2\pi\sigma^2)-\frac{\sum_{i=1}^n(x_i-\mu)^2}{2\sigma^2}$$

$$\frac{\partial}{\partial\mu}\log L=\sum_{i=1}^n\frac{x_i-\mu}{\sigma^2}=0,\qquad
\frac{\partial}{\partial\sigma^2}\log L=-\frac{n}{2\sigma^2}+\frac{\sum_{i=1}^n(x_i-\mu)^2}{2\sigma^4}=0$$

$$(\hat\mu,\ \hat\sigma^2)=\Big(\overline{X}_n,\ \sum_{i=1}^n(X_i-\overline{X}_n)^2/n\Big)$$

Verification of the maximum without the $2\times2$ second-derivative matrix: from the identity $\sum(x_i-\mu)^2=\sum(x_i-\bar x_n)^2+n(\bar x_n-\mu)^2$, any $\mu\ne\bar x_n$ makes $\sum(x_i-\mu)^2$ strictly larger, so $\log L(\mu,\sigma^2)$ is maximized at $\mu=\bar x_n$ **for every** $\sigma^2$. Substituting $\mu=\bar x_n$ leaves a function of $\sigma^2$ alone, maximized at $n^{-1}\sum(x_i-\bar x_n)^2$.

**Example 4.2.6 (Ex 4.8) — logistic distribution, no closed form**

$$f(x;\theta)=\frac{\exp\{-(x-\theta)\}}{(1+\exp\{-(x-\theta)\})^2},\qquad
\log L(\theta)=n\theta-n\bar x_n-2\sum_{i=1}^n\log\big(1+\exp(-(x_i-\theta))\big)$$

$$\frac{d}{d\theta}\log L(\theta)=n-2\sum_{i=1}^n\frac{\exp(-(x_i-\theta))}{1+\exp(-(x_i-\theta))}=0$$

has no simple closed-form root, but the derivative is strictly decreasing in $\theta$ (the second derivative is negative) and tends to $n$ as $\theta\to-\infty$ and to $-n$ as $\theta\to+\infty$, so the root is **unique** and can be found numerically.

**Example 4.2.7 (Ex 4.9) — uniform $U(0,\theta)$: maximizing by inspection**

$$L(\theta\,;x_1,\dots,x_n)=\prod_{i=1}^n f(x_i;\theta)=
\begin{cases}\dfrac1{\theta^n},&0\le x_i\le\theta\ \text{for all }i\\[4pt]0,&\text{otherwise}\end{cases}$$

This is not continuous in $\theta$, hence not differentiable, so inspect the shape instead. The condition "$0\le x_i\le\theta$ for all $i$" is the condition $\theta\ge x_{(n)}$, and $\theta^{-n}$ is decreasing, so $L$ jumps up at $\theta=x_{(n)}$ and decreases thereafter: the MLE is the sample maximum $X_{(n)}$.

**Example 4.2.8 (Ex 4.10) — the MLE need not exist (Poisson with $\bar x_n=0$)**

In Example 4.2.4 the MLE is $\overline{X}_n$ provided $\bar x_n>0$. If the observed $\bar x_n$ is $0$, then $\frac{d}{d\lambda}\log L=-n<0$ for every $\lambda$, and since the parameter space $(0,\infty)$ is **open**, no maximizer exists.

**Example 4.2.9 (Ex 4.11) — the MLE need not be unique ($U[\theta-1,\theta+1]$)**

$$L(\theta)=\prod_{i=1}^n\frac12 I_{[\theta-1,\theta+1]}(x_i)
=2^{-n}I_{[\theta-1\le x_{(1)},\ x_{(n)}\le\theta+1]}
=2^{-n}I_{[x_{(n)}-1,\ x_{(1)}+1]}(\theta)$$

attains its maximum $2^{-n}$ at **every** point of $[x_{(n)}-1,\ x_{(1)}+1]$, so every point of $[X_{(n)}-1,\ X_{(1)}+1]$ is an MLE.

**Example 4.2.10 (Ex 4.12) — Hardy–Weinberg model**

$$(N_1,N_2,N_3)\sim\mathrm{MULT}(n,p_1,p_2,p_3),\qquad
p_1=\theta^2,\ \ p_2=2\theta(1-\theta),\ \ p_3=(1-\theta)^2,\ \ 0<\theta<1$$

$$L(\theta)=\frac{n!}{n_1!n_2!n_3!}\theta^{2n_1}\big(2\theta(1-\theta)\big)^{n_2}(1-\theta)^{2n_3}$$

$$\log L(\theta)=\log\frac{n!}{n_1!n_2!n_3!}+2n_1\log\theta+n_2\log2\theta+n_2\log(1-\theta)+2n_3\log(1-\theta)$$

$$\frac{d}{d\theta}\log L(\theta)=\frac{2n_1}{\theta}+\frac{n_2}{\theta}-\frac{n_2}{1-\theta}-\frac{2n_3}{1-\theta}=0
\ \Longrightarrow\ \theta=\frac{2n_1+n_2}{2n}$$

$$\frac{d^2}{d\theta^2}\log L(\theta)=-\frac{2n_1+n_2}{\theta^2}-\frac{n_2+2n_3}{(1-\theta)^2}<0$$

so for $n_1\ne n$ and $n_3\ne n$ the MLE is $\hat\theta=\dfrac{2N_1+N_2}{2n}$, where $n=N_1+N_2+N_3$.

**Theorem 4.2.2 (Thm 4.1) — invariance of the MLE**

If $\hat\theta_n$ is the MLE of $\theta$, then for any function $g$,

$$g(\hat\theta_n)\ \text{is the MLE of}\ g(\theta)$$

Proof (for one-to-one $g$; the property holds generally). Put $\eta=g(\theta)$, let $L$ be the likelihood in the parametrization $\theta$ and $L^*$ the likelihood in the parametrization $\eta$. Then

$$L^*(g(\theta))=L^*(\eta)=\prod_{i=1}^n f(x_i;g^{-1}(\eta))=L(g^{-1}(\eta))=L(\theta)$$

so the two parametrizations have corresponding maximizers:

$$\max_{\eta\in g(\Omega)}L^*(\eta)=\max_{\eta\in g(\Omega)}L(g^{-1}(\eta))=\max_{\theta\in\Omega}L(\theta)$$

$$L^*(\hat\eta)=\max_{\eta\in g(\Omega)}L^*(\eta)=\max_{\theta\in\Omega}L(\theta)=L(\hat\theta_n)=L^*(g(\hat\theta_n))$$

hence $\hat\eta=g(\hat\theta_n)$. $\blacksquare$

Practical use: obtaining an MLE generally requires differentiation, so parametrize in whatever way makes the differentiation easy, then transform by invariance. The property extends to vector parameters: if $\hat\theta=(\hat\theta_1,\dots,\hat\theta_k)$ is the MLE of $\theta$, then $g(\hat\theta_1,\dots,\hat\theta_k)$ is the MLE of $g(\theta_1,\dots,\theta_k)$.

**Example 4.2.11 (Ex 4.13) — the standard deviation $\sigma$**

From $\hat\sigma^2=\sum_{i=1}^n(X_i-\overline{X}_n)^2/n$ (Example 4.2.5), invariance gives $\hat\sigma=\sqrt{\hat\sigma^2}$.

**Example 4.2.12 (Ex 4.14) — $P(X\ge1)$ for the exponential**

$\log L(\lambda)=-n\log\lambda-\sum x_i/\lambda$ gives $\hat\lambda=\overline{X}_n$, and for $g(\lambda)=P(X\ge1)=\exp(-1/\lambda)$,

$$g(\hat\lambda)=\exp\big(-1/\overline{X}_n\big)$$

## 4.3 Criteria for Estimation

An estimator is a function of the sample, hence a random variable with a distribution of its own; inspecting that distribution reveals the character of the estimator.

**Example 4.3.1 (Ex 4.15) — distribution of the MLE $X_{(n)}$ for $U(0,\theta)$**

With $\theta=1$ and $n=5$,

$$F_{X_{(5)}}(x)=P(X_{(5)}\le x)=P(X_1\le x,\dots,X_5\le x)=\{P(X\le x)\}^5=\{F_X(x)\}^5=x^5,\qquad 0<x<1$$

$$f(x)=5x^4,\qquad 0<x<1$$

The second equality holds because the sample maximum is below $x$ exactly when all observations are, the third because the $X_i$ are iid. Probability concentrates near $\theta=1$, but the distribution is **asymmetric** about $\theta$.

**Example 4.3.2 (Ex 4.16) — distribution of $\overline{X}_n$**

For a random sample from $N(\mu,\sigma^2)$, the estimator $\overline{X}_n$ of $\mu$ has distribution $N(\mu,\sigma^2/n)$ — centered at the parameter, **symmetric** about it, and concentrated near it.

**Definition 4.3.1 — loss function, mean square error**

We want $T(X)$ to be "close" to $g(\theta)$; closeness depends on the criterion, e.g.

- **absolute error**: $\ |T(X)-g(\theta)|$
- **squared error**: $\ (T(X)-g(\theta))^2$
- **weighted squared error**: $\ (T(X)-g(\theta))^2w(\theta)$

These **loss functions** are random variables, so no estimator can minimize them pointwise; their **expectations** are used instead. The usual criterion is the **mean square error**

$$\mathrm{MSE}=E\big(T(X)-g(\theta)\big)^2 \tag{4.5}$$

The MSE is a **function of $\theta$**, so an estimator minimizing it need not exist. Comparing two estimators, one may have smaller MSE for every $\theta$ (then it is preferred), or the ordering may reverse with $\theta$ (then the comparison is inconclusive).

No estimator minimizes the MSE **uniformly** in $\theta\in\Omega$. For instance $T(X)\equiv\theta_0$, ignoring the sample entirely, has $\mathrm{MSE}=0$ when the true value happens to be $\theta_0$, and nothing can beat that; yet its MSE grows without bound as $\theta$ moves away from $\theta_0$. Restricting attention to **unbiased** estimators excludes such degenerate rules.

**Definition 4.3.2 (Def 4.3) — bias, unbiased estimator**

For an estimator $T(X)$ of $g(\theta)$:

- **bias**: $\ E[T(X)]-g(\theta)$
- $T(X)$ is **unbiased** if $E[T(X)]=g(\theta)$, i.e. the bias is $0$
- the variance of the estimator is $E\big[T(X)-E[T(X)]\big]^2$

An unbiased estimator has expectation equal to the target, so it neither systematically over- nor under-estimates.

**Example 4.3.3 (Ex 4.17) — two unbiased estimators of $\mu$**

For a random sample of size $10$ from $N(\mu,\sigma_0^2)$ with $\sigma_0^2$ known, both $T_1(X)=\overline{X}_{10}$ and $T_2(X)=(X_1+X_2)/2$ satisfy $E[T_i(X)]=\mu$.

**Example 4.3.4 (Ex 4.18) — $X_{(n)}$ for $U(0,\theta)$: bias and its correction**

Since every observation lies in $[0,\theta]$, the sample maximum can never exceed $\theta$ and its expectation must fall short of it:

$$E(X_{(n)})=\int_0^\theta t\cdot\frac{nt^{n-1}}{\theta^n}dt=\frac{n\,t^{n+1}}{(n+1)\theta^n}\Big|_0^\theta=\frac{n}{n+1}\theta$$

so $X_{(n)}$ is **not** unbiased. Rescaling,

$$T_1(X)=\Big(\frac{n+1}{n}\Big)X_{(n)}\qquad\Longrightarrow\qquad
E(T_1(X))=\frac{n+1}{n}\cdot\frac{n}{n+1}\theta=\theta$$

**Theorem 4.3.1 (Thm 4.2) — MSE = variance + bias$^2$**

$$\mathrm{MSE}=\operatorname{Var}(T(X))+(\text{bias})^2$$

Proof. Add and subtract $E(T(X))$ and expand:

$$\begin{aligned}
\mathrm{MSE}&=E\big[T(X)-g(\theta)\big]^2\\
&=E\big[\{T(X)-E(T(X))\}+\{E(T(X))-g(\theta)\}\big]^2\\
&=E\{T(X)-E(T(X))\}^2+\{E(T(X))-g(\theta)\}^2+2E\{T(X)-E(T(X))\}\{E(T(X))-g(\theta)\}\\
&=\operatorname{Var}(T(X))+(\text{bias})^2
\end{aligned}$$

The cross term vanishes: $\{E(T(X))-g(\theta)\}$ is a constant and comes out of the expectation, leaving $E\{T(X)-E(T(X))\}=0$. $\blacksquare$

Estimating under the MSE criterion therefore trades off **variance** against **bias**, and a good estimator needs both small. Within the class of **unbiased** estimators the bias term is $0$, so minimizing the MSE means minimizing the **variance** — uniformly in $\theta$ if possible.

**Definition 4.3.3 (Def 4.4) — relative efficiency**

For unbiased estimators $T_1(X),T_2(X)$ of $g(\theta)$, the **relative efficiency** of $T_1$ with respect to $T_2$ is the variance ratio

$$\frac{\operatorname{Var}[T_2(X)]}{\operatorname{Var}[T_1(X)]}$$

which exceeds $1$ exactly when $T_1$ has the smaller variance and is therefore the better estimator under the MSE criterion.

**Example 4.3.5 (Ex 4.19) — continuing Example 4.3.3**

$$\operatorname{Var}(T_1(X))=\frac{\sigma_0^2}{10},\qquad \operatorname{Var}(T_2(X))=\frac{\sigma_0^2}{2}$$

The efficiency of $T_2$ relative to $T_1$ is $2/10=0.2$, so $T_2$ is the less efficient estimator.

**Example 4.3.6 (Ex 4.20) — $U(0,\theta)$: maximum-based vs minimum-based**

The unbiased estimator $T_1(X)=\big(\frac{n+1}{n}\big)X_{(n)}$ of Example 4.3.4 has density

$$f_{T_1}(t)=\Big(\frac{n}{n+1}\Big)^{\!n}\times\frac{n\,t^{n-1}}{\theta^n},\qquad 0<t<\Big(\frac{n+1}{n}\Big)\theta$$

For the sample minimum, the density (3.13) is $f_{X_{(1)}}(t;\theta)=\frac n\theta\big(1-\frac t\theta\big)^{n-1}$ on $0<t<\theta$, so

$$E(X_{(1)})=\int_0^\theta\frac{nt}{\theta}\Big(1-\frac t\theta\Big)^{\!n-1}dt=\frac{\theta}{n+1}
\qquad\Longrightarrow\qquad T_2(X)=(n+1)X_{(1)}\ \text{is unbiased}$$

$$E(T_1^2)=\int_0^{(n+1)\theta/n}\Big(\frac{n}{n+1}\Big)^{\!n}\frac{n\,t^{n+1}}{\theta^n}dt=\frac{(n+1)^2}{n(n+2)}\theta^2$$

$$\operatorname{Var}(T_1(X))=\frac{(n+1)^2}{n(n+2)}\theta^2-\theta^2=\frac{\theta^2}{n(n+2)},
\qquad \operatorname{Var}(T_2(X))=\frac{n\theta^2}{n+2}$$

$$\frac{\operatorname{Var}(T_2(X))}{\operatorname{Var}(T_1(X))}=\frac{n\theta^2/(n+2)}{\theta^2/n(n+2)}=n^2$$

so $T_1$ is better for $n>1$, overwhelmingly so for large $n$.

**Remark 4.3.1 — relative efficiency as a ratio of sample sizes**

If $\operatorname{Var}(T_1(X))=\nu_1/n$ and $\operatorname{Var}(T_2(X))=\nu_2/n$, the efficiency of $T_1$ relative to $T_2$ is $\nu_2/\nu_1$. Letting $n_1,n_2$ be the sample sizes at which the two methods attain the **same** variance, $\nu_1/n_1=\nu_2/n_2$ gives $\nu_2/\nu_1=n_2/n_1$: relative efficiency is the ratio of sample sizes needed for equal precision.

**Definition 4.3.4 — $100\%$ efficiency**

$T^*(X)$ is said to have $100\%$ efficiency if, for **every** unbiased estimator $T(X)$ of $g(\theta)$,

$$\frac{\operatorname{Var}(T^*(X))}{\operatorname{Var}(T(X))}\le1 \tag{4.6}$$

Finding such an estimator — a minimum variance unbiased estimator — is the subject of Section 4.4.

## 4.4 Minimum Variance Unbiased Estimator

**Definition 4.4.1 (Def 4.5) — minimum variance unbiased estimator (MVUE)**

$T^*(X)$ is the **MVUE** of $g(\theta)$ if

1. $T^*(X)$ is unbiased for $g(\theta)$, i.e. $E[T^*(X)]=g(\theta)$
2. $\operatorname{Var}(T^*(X))\le\operatorname{Var}(T(X))$ for **every** unbiased estimator $T(X)$ of $g(\theta)$

Two routes to an MVUE: the **Cramér–Rao lower bound** (a lower bound on the variance of any unbiased estimator — check whether some estimator attains it), and **complete sufficient statistics** via the **Rao–Blackwell** and **Lehmann–Scheffé** theorems.

### 4.4.1 Cramér–Rao lower bound

**Definition 4.4.2 — Fisher's information**

For a density $f(x;\theta)$,

$$I(\theta)=\operatorname{Var}\Big(\frac{\partial}{\partial\theta}\log f(X;\theta)\Big)
=E\left[\Big(\frac{\partial}{\partial\theta}\log f(X;\theta)\Big)^{\!2}\right] \tag{4.7}$$

the expected squared rate of change of the log-density in $\theta$. The larger $I(\theta)$, the more one observation says about the parameter.

**Theorem 4.4.1 — Fisher information = variance of the score**

Call $\frac{\partial}{\partial\theta}\log f(X;\theta)$ the **score**. Then the two expressions in (4.7) agree.

Proof. By definition of variance,

$$\operatorname{Var}\Big(\frac{\partial}{\partial\theta}\log f(X;\theta)\Big)
=E\left[\Big(\frac{\partial}{\partial\theta}\log f\Big)^{\!2}\right]-\left[E\Big(\frac{\partial}{\partial\theta}\log f\Big)\right]^{2}$$

The second term vanishes: since $\frac{\partial}{\partial\theta}\log f=\dfrac{\partial f/\partial\theta}{f}$, the denominator cancels against the density, and interchanging differentiation and integration,

$$E\Big(\frac{\partial}{\partial\theta}\log f(X;\theta)\Big)
=\int\frac{\frac{\partial}{\partial\theta}f(x;\theta)}{f(x;\theta)}f(x;\theta)dx
=\frac{\partial}{\partial\theta}\int f(x;\theta)dx=\frac{\partial}{\partial\theta}(1)=0 \qquad\blacksquare$$

**Theorem 4.4.2 — second-derivative form**

$$I(\theta)=E\left[\Big(\frac{\partial}{\partial\theta}\log f\Big)^{\!2}\right]
=-E\left[\frac{\partial^2}{\partial\theta^2}\log f(X;\theta)\right]$$

Proof. By the quotient rule,

$$\frac{\partial^2}{\partial\theta^2}\log f(X;\theta)
=\frac{\partial}{\partial\theta}\left(\frac{\frac{\partial}{\partial\theta}f}{f}\right)
=\frac{\big(\frac{\partial^2}{\partial\theta^2}f\big)f-\big(\frac{\partial}{\partial\theta}f\big)^2}{f^2}
=\frac{\frac{\partial^2}{\partial\theta^2}f}{f}-\Big(\frac{\partial}{\partial\theta}\log f\Big)^{\!2}$$

Taking expectations, the first term vanishes:

$$E\left[\frac{\frac{\partial^2}{\partial\theta^2}f}{f}\right]
=\int\frac{\partial^2}{\partial\theta^2}f(x;\theta)dx
=\frac{\partial^2}{\partial\theta^2}\int f(x;\theta)dx=\frac{\partial^2}{\partial\theta^2}(1)=0 \qquad\blacksquare$$

**Remark 4.4.1 — what Fisher information measures**

It is the variance of the **slope** of the log-likelihood, not of the likelihood itself. A large score variance means the likelihood reacts sharply to small changes in $\theta$, i.e. it is **sharply peaked** near the true value — and the more peaked, the more precisely the data pin the parameter down.

**Example 4.4.1 (Ex 4.21) — normal**

For $X\sim N(\mu,\sigma_0^2)$ with $\sigma_0^2$ known, $f(x;\mu)=(\sqrt{2\pi}\sigma_0)^{-1}\exp[-(x-\mu)^2/2\sigma_0^2]$, so

$$I(\mu)=E\big[(\partial/\partial\mu)\log f(X;\mu)\big]^2=E\left[\frac{X-\mu}{\sigma_0^2}\right]^{\!2}=\frac1{\sigma_0^2}$$

The smaller the variance, the larger the information: observations lie closer to $\mu$, so each says more about it. If $\sigma_1^2<\sigma_2^2$, use $X\sim N(\mu,\sigma_1^2)$ rather than $Y\sim N(\mu,\sigma_2^2)$ to estimate $\mu$.

**Example 4.4.2 (Ex 4.22) — Bernoulli**

$\log f_X(X;p)=X\log p+(1-X)\log(1-p)$, so

$$\frac{\partial\log f_X(X;p)}{\partial p}=\frac Xp-\frac{1-X}{1-p}=\frac{X-p}{p(1-p)},
\qquad I(p)=\frac{E(X-p)^2}{p^2(1-p)^2}=\frac1{p(1-p)}$$

the reciprocal of the Bernoulli variance; the information is smallest at $p=1/2$.

**Theorem 4.4.3 (Thm 4.3) — information inequality (Cramér–Rao lower bound)**

Assume

1. distinct parameters give distinct densities: $\theta\ne\theta'\Rightarrow f(x;\theta)\ne f(x;\theta')$
2. $A=\{x:f(x;\theta)>0\}$ does not depend on $\theta$, and $\log f(x;\theta)$ is twice differentiable in $\theta$ with continuous derivatives for all $x\in A$, $\theta\in\Omega$
3. differentiation and integration may be interchanged whenever $E(T(X))<\infty$

If $X_1,\dots,X_n$ is a random sample from $f(x;\theta)$ with $\operatorname{Var}(T(X))<\infty$, $E(T(X))=g(\theta)$ and $0<I(\theta)<\infty$ for all $\theta\in\Omega$, then $g$ is differentiable and

$$\operatorname{Var}(T(X))\ \ge\ \frac{[g'(\theta)]^2}{n\,I(\theta)} \tag{4.8}$$

In particular $\operatorname{Var}(T(X))\ge\dfrac1{n\,I(\theta)}$ when $g(\theta)=\theta$.

Proof. All integrals are $n$-fold over $(-\infty,\infty)$, and $T=T(x_1,\dots,x_n)$.

**(1) Start from unbiasedness.** $g(\theta)=E(T(X))$, and by assumption (3),

$$g'(\theta)=\frac{\partial}{\partial\theta}\int\!\cdots\!\int T\prod_{i=1}^n f(x_i;\theta)dx_1\cdots dx_n
=\int\!\cdots\!\int T\,\frac{\partial}{\partial\theta}\left[\prod_{i=1}^n f(x_i;\theta)\right]dx_1\cdots dx_n$$

**(2) The log-derivative trick.** From $\frac{\partial}{\partial\theta}\log\prod f=\dfrac{\partial\prod f/\partial\theta}{\prod f}$ we get $\frac{\partial}{\partial\theta}\prod f=\big[\frac{\partial}{\partial\theta}\log\prod f\big]\prod f$, which restores the integrand to the form $(\cdots)\times\prod f$:

$$g'(\theta)=\int\!\cdots\!\int T\left[\frac{\partial}{\partial\theta}\log\prod_{i=1}^n f(x_i;\theta)\right]\prod_{i=1}^n f\,dx_1\cdots dx_n$$

**(3) Subtract $g(\theta)$ times zero.** Differentiating $\int\!\cdots\!\int\prod f=1$ gives

$$\int\!\cdots\!\int\left[\frac{\partial}{\partial\theta}\log\prod_{i=1}^n f\right]\prod_{i=1}^n f\,dx_1\cdots dx_n=0$$

so subtracting $g(\theta)$ times this expression changes nothing:

$$g'(\theta)=E\left[\big(T(X)-g(\theta)\big)\Big(\frac{\partial}{\partial\theta}\log\prod_{i=1}^n f(X_i;\theta)\Big)\right]$$

**(4) Cauchy–Schwarz.** With $U=T(X)-g(\theta)$ and $V=\frac{\partial}{\partial\theta}\log\prod f$ in $[E(UV)]^2\le E(U^2)E(V^2)$,

$$[g'(\theta)]^2\le E\big[\{T(X)-g(\theta)\}^2\big]\cdot
E\left[\Big(\frac{\partial}{\partial\theta}\log\prod_{i=1}^n f(X_i;\theta)\Big)^{\!2}\right]
=\operatorname{Var}(T(X))\cdot E\left[\Big(\frac{\partial}{\partial\theta}\log\prod_{i=1}^n f\Big)^{\!2}\right]$$

using that $T$ is unbiased, so $E[\{T-g(\theta)\}^2]$ is its variance.

**(5) Each score has mean zero.** With $S_i=\frac{\partial}{\partial\theta}\log f(X_i;\theta)$, the same interchange as in (1) and (3) gives

$$E(S_i)=\int\Big[\frac{\partial}{\partial\theta}\log f(x;\theta)\Big]f(x;\theta)dx
=\frac{\partial}{\partial\theta}\int f(x;\theta)dx=\frac{\partial}{\partial\theta}(1)=0$$

**(6) The right-hand expectation is $n\,I(\theta)$.** The logarithm turns the product into a sum, $\frac{\partial}{\partial\theta}\log\prod_{i=1}^n f(X_i;\theta)=\sum_{i=1}^n S_i$, and expanding the square,

$$E\left[\Big(\sum_{i=1}^n S_i\Big)^{\!2}\right]=\sum_{i=1}^n\sum_{j=1}^n E(S_iS_j)$$

The cross terms vanish, since $X_i\perp X_j$ for $i\ne j$ gives $E(S_iS_j)=E(S_i)E(S_j)=0$ by (5); only the $n$ diagonal terms survive:

$$E\left[\Big(\sum_{i=1}^n S_i\Big)^{\!2}\right]=n\,E\left[\Big(\frac{\partial}{\partial\theta}\log f(X;\theta)\Big)^{\!2}\right]=n\,I(\theta)$$

**(7) Conclusion.** Substituting (6) into (4) gives $[g'(\theta)]^2\le\operatorname{Var}(T(X))\cdot n\,I(\theta)$, and dividing by $n\,I(\theta)>0$ gives (4.8). $\blacksquare$

An unbiased estimator with $\operatorname{Var}(T(X))=1/nI(\theta)$ therefore has the smallest variance attainable and is the MVUE of $\theta$.

**Example 4.4.3 (Ex 4.23) — the Poisson sample mean is the MVUE of $\lambda$**

$$I(\lambda)=E\left[\frac{\partial}{\partial\lambda}\log f(X;\lambda)\right]^2
=E\left[-1+\frac X\lambda\right]^2=E\left[\frac{X-\lambda}{\lambda}\right]^2
=\frac{\operatorname{Var}(X)}{\lambda^2}=\frac1\lambda$$

and $\operatorname{Var}(\overline{X}_n)=\lambda/n$ attains the bound in (4.8).

**Example 4.4.4 (Ex 4.24) — the normal sample mean is the MVUE of $\mu$**

$I(\mu)=1/\sigma^2$ (Example 4.4.1) and $\operatorname{Var}(\overline{X}_n)=\sigma^2/n$, which attains the bound.

**Remark 4.4.2 — limitation of the Cramér–Rao method**

No estimator need attain the bound. When the bound is too small to be attained, the method cannot identify the MVUE — it only certifies that no smaller variance is possible. The second route, complete sufficiency, handles such cases.

### 4.4.2 Sufficient statistics

**Definition 4.4.3 (Def 4.6) — (jointly) sufficient statistic**

Let $X=(X_1,\dots,X_n)$ have joint density $f(x_1,\dots,x_n;\theta_1,\dots,\theta_k)$ and let $S(X)=(S_1(X),\dots,S_l(X))$ be a vector of $l$ statistics. $S$ is a **jointly sufficient statistic** for $\theta=(\theta_1,\dots,\theta_k)$ if the conditional distribution of $(X_1,\dots,X_n)\mid S(X)$ **does not depend on $\theta$**. For $l=1$, $S(X)$ is a **sufficient statistic** for $\theta$.

**Remark 4.4.3 — data reduction**

Once the value of $S$ is given, the conditional distribution of the sample carries no information about the parameter: **all information about $\theta$ is contained in $S(X)$**. The sample splits into "information about $\theta$ ($S$)" plus "a part unrelated to $\theta$", so a large data set may be replaced by the sufficient statistic with no loss of estimation efficiency — **data reduction**.

Concretely: $1000$ customers visit a shop and the mean spend per customer is wanted.

- the raw data are $X_1=x_1,\dots,X_{1000}=x_{1000}$, a list of $1000$ numbers recording who spent what
- the **total** $S=\sum_{i=1}^{1000}X_i$, a single number, is all that is needed to estimate the mean
- "what the third customer spent", "in what order they paid" contributes nothing to that estimate
- a $1000$-dimensional data set is thus reduced to a $1$-dimensional one with no loss of information about the parameter

**Example 4.4.5 (Ex 4.25) — Bernoulli, by direct computation of the conditional distribution**

Let $S=\sum_{i=1}^n X_i$ for a random sample from Bernoulli($p$). Since $S\sim B(n,p)$ we have $f_S(s;p)=\binom ns p^s(1-p)^{n-s}$, so for $\sum x_i=s$,

$$f(x_1,\dots,x_n\mid s)=\frac{f(x_1,\dots,x_n;p)}{f_S(s;p)}
=\frac{p^s(1-p)^{n-s}}{\binom ns p^s(1-p)^{n-s}}=\frac{1}{\binom ns}$$

which is free of $p$: $S=\sum X_i$ is sufficient for $p$.

**Remark 4.4.4 — what is discarded**

Toss a coin $5$ times and observe $(1,0,1,1,0)$.

- **original information** $X_1=x_1,\dots,X_n=x_n$: a record of $n$ outcomes *including their order* ("first head, second tail, third head, …"); the joint density is the product of these ordered outcome probabilities
- **reduced information** $S=\sum X_i$: the single number "how many heads in total", here $(1,0,1,1,0)\mapsto S=3$
- **discarded information**: *which* tosses were heads — of no use in learning $p$

The result $f(x_1,\dots,x_n\mid S=s)=1/\binom ns$ says: given that there were $s$ heads in total, the probability that the data were specifically $(1,0,1,1,0)$ is just the uniform probability of one of $\binom ns$ arrangements, with no information about $p$ left in it. The order is redundant for estimation, so reducing the list of $n$ values to the single total loses nothing.

**Remark 4.4.5 — why the parameter cancels**

By definition $f(x_1,\dots,x_n\mid s)=\dfrac{f(x_1,\dots,x_n;\theta)}{f_S(s;\theta)}$, so read numerator and denominator separately:

- **numerator** $f(x_1,\dots,x_n;\theta)$ — the whole sample: information about $\theta$ together with the rest of the sample information
- **denominator** $f_S(s;\theta)$ — the summary $S$: if $S$ is sufficient it has absorbed **all** of the part connected to $\theta$
- **the division** — the $\theta$-factor in the denominator matches the one in the numerator exactly, so it cancels wholesale and no $\theta$ survives

In Example 4.4.5, $p^s(1-p)^{n-s}$ cancels and only $1/\binom ns$ remains. Conversely, if $\theta$ fails to cancel, then $S$ has not captured all the information about $\theta$ and is not sufficient. The factorization theorem checks the structure "numerator $=$ ($\theta$-part) $\times$ ($\theta$-free part)" without computing any conditional probability.

**Theorem 4.4.4 (Thm 4.4) — factorization theorem (Neyman–Fisher)**

$S(X)=(S_1(X),\dots,S_k(X))$ is jointly sufficient **if and only if** the joint density factors into a function $g$ of $s$ and $\theta$ only and a function $h$ of $(x_1,\dots,x_n)$ only:

$$f(x_1,x_2,\dots,x_n;\theta)=g(s(x);\theta)\times h(x_1,x_2,\dots,x_n)$$

In one line: if the **only channel** through which $\theta$ touches the data is $S$, then $S$ is sufficient — and the theorem tests that condition in multiplicative form.

**Notation.** The rule $s(\cdot)$ is fixed once and for all, and $\theta$ never enters its definition.

- $S=s(X_1,\dots,X_n)$ — a random variable, with its own distribution $f_S(s;\theta)$
- $s=s(x_1,\dots,x_n)$ — the value attached to the observed data $x$; that is, $s$ is a **function of $x$**

Below, $x=(x_1,\dots,x_n)$ and $A_s:=s^{-1}(\{s\})=\{x:s(x)=s\}$.

Proof (discrete case). Since $s(x)$ is determined by $x$, the events $\{X=x\}$ and $\{X=x,\ S=s(x)\}$ coincide, so $P(A\cap B)=P(B)P(A\mid B)$ gives

$$f(x;\theta)=f_S\big(s(x);\theta\big)\cdot f_{X\mid S}\big(x\mid s(x)\big) \tag{$\ast$}$$

($\Rightarrow$) If $S$ is sufficient, then by definition $f_{X\mid S}$ contains no $\theta$. Putting

$$g(s;\theta):=f_S(s;\theta),\qquad h(x):=f_{X\mid S}\big(x\mid s(x)\big)$$

in $(\ast)$ gives $f(x;\theta)=g(s(x);\theta)h(x)$. Here $h$ is the composite $x\mapsto s(x)\mapsto f_{X\mid S}(x\mid s(x))$, and since $s$ is a function of $x$, $h$ is a function of $x$ alone; by sufficiency it contains no $\theta$.

($\Leftarrow$) Conversely let $f(x;\theta)=g(s(x);\theta)h(x)$. Since $\{S=s\}=\{X\in A_s\}$ and $g(s(x);\theta)=g(s;\theta)$ is constant on $A_s$,

$$f_S(s;\theta)=\sum_{x\in A_s}f(x;\theta)=g(s;\theta)\sum_{x\in A_s}h(x)=g(s;\theta)c(s),
\qquad c(s):=\sum_{x\in A_s}h(x)$$

Hence for $x\in A_s$, where $\{X=x,S=s\}=\{X=x\}$,

$$f_{X\mid S}(x\mid s)=\frac{f(x;\theta)}{f_S(s;\theta)}=\frac{g(s;\theta)h(x)}{g(s;\theta)c(s)}=\frac{h(x)}{c(s)}$$

and it is $0$ for $x\notin A_s$. Neither $h$ nor $c$ contains $\theta$, so the conditional distribution does not depend on $\theta$. $\blacksquare$

**Remark 4.4.6 — reading $s$ as $s(x)$**

$f_{X\mid S}(x\mid s)$ looks like a function of the two variables $x$ and $s$, but the $s$ in the condition is not supplied from outside — it is $s(x)$. The only free variable is $x$, and

$$h(x)=f_{X\mid S}\big(x\mid s(x)\big),\qquad f_{X\mid S}\big(x\mid s(x)\big)=\frac{h(x)}{c(s(x))}$$

are both functions of $x$ alone. What disappears is not $s$ but $\theta$; $s$ survives inside the formula through $x$.

**Remark 4.4.7 — $g$ and $h$ are not unique**

The theorem asserts that such a factorization **exists**, and $g=f_S$, $h=f_{X\mid S}$ is merely the candidate exhibited in ($\Rightarrow$). In applications any factorization meeting the conditions ($g$ in $s,\theta$ only; $h$ in $x$ only) will do — moving a constant from $g$ into $h$ gives another one.

**Remark 4.4.8 — the level set $A_s$ concretely**

Bernoulli with $n=3$ and $s=2$: $A_2=\{(1,1,0),(1,0,1),(0,1,1)\}$ with $|A_2|=\binom32$. All three have the same probability $p^2(1-p)$, so $f_S(2;p)=\binom32p^2(1-p)$, and with $h\equiv1$, $c(2)=\binom32$ and $f_{X\mid S}(x\mid2)=1/\binom32$ — matching Example 4.4.5. So $A_s$ is the level set of data sharing one summary value, and $c(s)$ is the total of $h$ over it.

**Remark 4.4.9 — sufficiency seen through level sets**

$S$ partitions the sample space into level sets $A_s$.

- sufficiency means: **within one bundle, whichever data arose, there is nothing further to learn about $\theta$**
- indeed the within-bundle distribution $h(x)/c(s)$ contains no $\theta$ — the arrangement inside a bundle is pure randomness unrelated to $\theta$
- so one may forget which data arose and remember only the bundle label $s$, with no loss of information about the parameter

**Example 4.4.6 (Ex 4.26) — Bernoulli by factorization**

From $f(x_1,\dots,x_n;p)=p^{\sum x_i}(1-p)^{n-\sum x_i}\prod_{i=1}^n I_{\{0,1\}}(x_i)$, take

$$g(s;p)=p^s(1-p)^{n-s}\ (\text{the parameter part}),\qquad
h(x_1,\dots,x_n)=\prod_{i=1}^n I_{\{0,1\}}(x_i)\ (\text{the data part})$$

so $f=g(s;p)h(x)$ and $\sum_{i=1}^n X_i$ is sufficient.

**Example 4.4.7 (Ex 4.27, 4.28) — uniform $U(\theta_1,\theta_2)$: minimum and maximum**

$$f(x_1,\dots,x_n;\theta_1,\theta_2)=\prod_{i=1}^n\frac{1}{\theta_2-\theta_1}I_{(\theta_1,\theta_2)}(x_i)
=\frac{1}{(\theta_2-\theta_1)^n}I_{(\theta_1,\,x_{(n)}]}(x_{(1)})I_{[x_{(1)},\,\theta_2)}(x_{(n)})$$

since all $x_i$ lying in the interval is equivalent to the sample minimum exceeding $\theta_1$ and the sample maximum falling below $\theta_2$. With $h\equiv1$, the factorization theorem makes $(X_{(1)},X_{(n)})$ jointly sufficient for $(\theta_1,\theta_2)$.

**Example 4.4.8 (Ex 4.29) — normal: $(\sum X_i,\sum X_i^2)$**

$$f(x_1,\dots,x_n;\mu,\sigma^2)=(\sqrt{2\pi}\sigma)^{-n}
\exp\left[-\frac{\sum_{i=1}^n x_i^2-2\mu\sum_{i=1}^n x_i+n\mu^2}{2\sigma^2}\right]$$

so $S_1=\sum_{i=1}^n X_i$ and $S_2=\sum_{i=1}^n X_i^2$ are jointly sufficient, and so are $\overline{X}_n$ and $\sum_{i=1}^n(X_i-\overline{X}_n)^2$, being in one-to-one correspondence with $(S_1,S_2)$.

**Remark 4.4.10 — sufficient statistics are not unique**

The order statistics $(X_{(1)},\dots,X_{(n)})$ are always jointly sufficient: $n!$ samples map to one ordered vector, each with conditional probability $1/n!$, free of $\theta$. Any **one-to-one function** of a sufficient statistic is again sufficient. Among all of them, one that reduces the data as far as possible while retaining all information about the parameter is called a **minimal sufficient statistic**.

### 4.4.3 Rao–Blackwell theorem

**Theorem 4.4.5 (Thm 4.5) — Rao–Blackwell**

Let $S$ be sufficient and let $T(X)$ be an unbiased estimator of $g(\theta)$. Then $\delta(S)=E(T(X)\mid S)$ is also unbiased for $g(\theta)$ and, for every $\theta$,

$$\operatorname{Var}(\delta(S))\ \le\ \operatorname{Var}(T(X))$$

Conditioning an unbiased estimator on a sufficient statistic therefore **improves** it.

Proof.

**(1) $\delta(S)$ is an estimator.** By **sufficiency** the conditional distribution does not depend on $\theta$, so $\delta(S)=E(T(X)\mid S)$ contains no unknown parameter and is a statistic.

**(2) Unbiasedness.** By the double expectation theorem,

$$E(\delta(S))=E\big[E(T(X)\mid S)\big]=E(T(X))=g(\theta)$$

**(3) Variance.** By the total variance formula (Theorem 2.13) with $Y=T(X)$ and $X=S$,

$$\operatorname{Var}(T(X))=E\big(\operatorname{Var}(T(X)\mid S)\big)+\operatorname{Var}\big(E(T(X)\mid S)\big)\ \ge\ \operatorname{Var}(\delta(S))$$

because the first term is nonnegative. Equality holds when $\operatorname{Var}(T(X)\mid S)=0$, i.e. $T(X)=\delta(S)$. $\blacksquare$

**Example 4.4.9 (Ex 4.30) — Rao–Blackwellizing $X_1$ in the Bernoulli case**

$E(X_1)=p$, so $X_1$ is unbiased for $p$. Conditioning on the sufficient statistic $S=\sum_{i=1}^n X_i$, and using

$$P(X_1=1\mid S=s)=\frac{p\binom{n-1}{s-1}p^{s-1}(1-p)^{n-s}}{\binom ns p^s(1-p)^{n-s}}=\frac sn$$

$$E\Big(X_1\Big|\sum_{i=1}^n X_i=s\Big)=0\cdot\Big(1-\frac sn\Big)+1\cdot\frac sn=\frac sn$$

so $E(T(X)\mid S)=S/n=\overline{X}_n$, and the variance indeed drops from $\operatorname{Var}(X_1)=p(1-p)$ to $\operatorname{Var}(\overline{X}_n)=p(1-p)/n$.

Rao–Blackwell improves any unbiased estimator into a function of a sufficient statistic. But is the result the MVUE? Different starting estimators may give different $\delta$, so a tool guaranteeing **uniqueness** is needed: completeness.

### 4.4.4 Completeness

**Definition 4.4.4 (Def 4.7) — complete statistic**

A statistic $S=S(X_1,\dots,X_n)$ is **complete** if the only function $g$ with

$$E\big(g(S)\big)=0\quad\text{for all }\theta\in\Omega$$

is $g(\cdot)\equiv0$. A statistic that is both sufficient and complete is a **complete sufficient statistic**.

**Example 4.4.10 (Ex 4.31) — a statistic that is not complete: $X_1-X_2$**

In a Bernoulli($p$) sample, $E(X_1-X_2)=0$ for every $p$, yet $X_1-X_2$ is not identically $0$: the nonzero function $g(s)=s$ satisfies $E(g(X_1-X_2))=0$, so $X_1-X_2$ is **not** complete.

**Example 4.4.11 (Ex 4.32) — checking completeness: binomial and uniform**

**(1) Binomial.** For $T=\sum_{i=1}^n X_i\sim B(n,p)$,

$$E[g(T)]=\sum_{t=0}^n g(t)\binom nt p^t(1-p)^{n-t}
=(1-p)^n\sum_{t=0}^n g(t)\binom nt\Big(\frac{p}{1-p}\Big)^{\!t}=0,\qquad 0<p<1$$

Holding for all $p$ means a polynomial in $\frac{p}{1-p}>0$ vanishes identically, so all coefficients vanish and $g(t)\equiv0$: $T$ is complete.

**(2) Uniform.** For $U(0,\theta)$ the sample maximum has density $f_{X_{(n)}}(y)=ny^{n-1}/\theta^n$, so

$$E[g(X_{(n)})]=\frac{n}{\theta^n}\int_0^\theta g(y)y^{n-1}dy=0\ \ (\forall\theta>0)
\iff \int_0^\theta g(y)y^{n-1}dy=0\ \ (\forall\theta>0)$$

Differentiating in $\theta$ gives $g(\theta)\theta^{n-1}=0$, hence $g(\cdot)\equiv0$. Since $X_{(n)}$ is also sufficient (Example 4.4.7), it is a **complete sufficient statistic** for $\theta$.

### 4.4.5 Lehmann–Scheffé theorem

**Theorem 4.4.6 (Thm 4.6) — Lehmann–Scheffé**

Let $S$ be a **complete sufficient** statistic for $\theta$ and let $T(X)$ be unbiased for $g(\theta)$. Then

$$\delta(S)=E\big(T(X)\mid S\big)$$

is the **unique MVUE** of $g(\theta)$.

Proof. By the double expectation theorem $\delta(S)$ is unbiased for $g(\theta)$. Let $\delta^*(S)$ be any unbiased estimator that is a **function of $S$**. Then for all $\theta$,

$$E\big[\delta(S)-\delta^*(S)\big]=g(\theta)-g(\theta)=0$$

and $\delta(S)-\delta^*(S)$ is a function of $S$, so **completeness** gives $P[\delta(S)=\delta^*(S)]=1$: there is essentially only one unbiased estimator among the functions of $S$.

Now let $\eta(X_1,\dots,X_n)$ be an arbitrary unbiased estimator, not necessarily a function of $S$. By Rao–Blackwell, $E(\eta\mid S)$ is unbiased with variance no larger than that of $\eta$, and being a function of $S$ it coincides with $\delta(S)$ by the uniqueness just shown. Hence $\delta(S)$ has the smallest variance among all unbiased estimators. $\blacksquare$

**Remark 4.4.11 — strategy for finding an MVUE**

Theorem 4.4.6 also says: if $S$ is complete sufficient and $E(S(X))=g(\theta)$, then $S(X)$ is the unique MVUE of $g(\theta)$. So find a complete sufficient statistic for $\theta$, then either

1. find an **unbiased** function of it, or
2. take a convenient unbiased estimator and **condition** it on the complete sufficient statistic.

**Example 4.4.12 (Ex 4.33) — MVUE for $U(0,\theta)$**

$X_{(n)}$ is complete sufficient (Examples 4.4.7 and 4.4.11) and $E[(n+1)X_{(n)}/n]=\theta$ (Example 4.3.4), so

$$\frac{(n+1)X_{(n)}}{n}$$

is the unique MVUE of $\theta$.

### 4.4.6 Exponential family

Complete sufficiency is very useful for point estimation, but proving completeness directly can require heavy computation. For the **exponential family** it is automatic, which covers many of the commonly used distributions.

**Definition 4.4.5 (Def 4.8) — exponential family**

If for suitable functions $a,b,c_i,t_i\ (i=1,\dots,k)$ the density can be written

$$f(x;\theta)=a(\theta)\,b(x)\,\exp\left[\sum_{i=1}^k c_i(\theta)t_i(x)\right],\qquad -\infty<x<\infty$$

with $\theta=(\theta_1,\dots,\theta_k)$, it belongs to the **exponential family** with $k$ parameters. Here $a$ depends on the parameter only, $b$ on the data only, the exponent is a sum of products of a parameter function $c_i$ and a data function $t_i$, and the set $\{x:f(x;\theta)>0\}$ must not depend on $\theta$.

**Example 4.4.13 (Ex 4.34) — Bernoulli, Poisson and normal are exponential families**

**(1) Bernoulli.** $f(x;p)=p^x(1-p)^{1-x}=(1-p)\exp\big[x\log\frac{p}{1-p}\big]$ for $x=0,1$, so $a(p)=1-p$, $b(x)=1$, $c(p)=\log\frac{p}{1-p}$, $t(x)=x$.

**(2) Poisson.** $f(x;\lambda)=\frac{e^{-\lambda}\lambda^x}{x!}=e^{-\lambda}\frac{1}{x!}\exp(x\log\lambda)$, so $a(\lambda)=e^{-\lambda}$, $b(x)=\frac1{x!}$, $c(\lambda)=\log\lambda$, $t(x)=x$.

**(3) Normal.**

$$f(x;\mu,\sigma^2)=\frac{1}{\sqrt{2\pi\sigma^2}}\exp\Big(-\frac{\mu^2}{2\sigma^2}\Big)
\exp\Big(-\frac{1}{2\sigma^2}x^2+\frac{\mu}{\sigma^2}x\Big)$$

so $a(\mu,\sigma^2)=\frac{1}{\sqrt{2\pi\sigma^2}}e^{-\mu^2/2\sigma^2}$, $b(x)=1$, $c_1=-\frac1{2\sigma^2},\ t_1(x)=x^2$, $c_2=\frac{\mu}{\sigma^2},\ t_2(x)=x$ — a two-parameter exponential family.

The exponential, gamma, beta and negative binomial distributions also belong. $U(0,\theta)$ and the shifted exponential $f(x;\theta)=\exp(-(x-\theta)),\ x>\theta$ do **not**, because the region where $f(x;\theta)>0$ depends on $\theta$.

**Theorem 4.4.7 (Thm 4.7) — complete sufficient statistic of an exponential family**

For a random sample $X_1,\dots,X_n$ from an exponential family $f(x;\theta)=a(\theta)b(x)\exp\big[\sum_{i=1}^k c_i(\theta)t_i(x)\big]$, the statistics

$$S_1=\sum_{i=1}^n t_1(X_i),\quad\dots,\quad S_k=\sum_{i=1}^n t_k(X_i)$$

are **jointly complete sufficient** for $\theta_1,\dots,\theta_k$. (The proof is omitted.)

Consequently any unbiased estimator of a function of $(\theta_1,\dots,\theta_k)$ that is built as a function of $(S_1,\dots,S_k)$ is the MVUE. A one-to-one function of a complete sufficient statistic is again complete sufficient.

**Example 4.4.14 (Ex 4.35) — Bernoulli: MVUE of $p(1-p)$**

$S=\sum_{i=1}^n X_i$ is complete sufficient (Example 4.4.13, Theorem 4.4.7). Trying the function $\overline{X}_n(1-\overline{X}_n)$ of it,

$$E\big[\overline{X}_n(1-\overline{X}_n)\big]=E(\overline{X}_n)-E(\overline{X}_n^2)
=p-\big[p^2+\operatorname{Var}(\overline{X}_n)\big]=p-p^2-\frac{p(1-p)}{n}=p(1-p)\Big[1-\frac1n\Big]$$

so after rescaling,

$$\frac{\overline{X}_n(1-\overline{X}_n)}{1-\frac1n}=\frac{n\,\overline{X}_n(1-\overline{X}_n)}{n-1}$$

is an unbiased function of a complete sufficient statistic — the MVUE of $p(1-p)$.

**Example 4.4.15 (Ex 4.36) — Poisson: MVUE of $\lambda$**

$$\prod_{i=1}^n f(x_i;\lambda)=\frac{\exp(-n\lambda)\lambda^{\sum_{i=1}^n x_i}}{\prod_{i=1}^n x_i!}$$

so $\sum_{i=1}^n X_i$ is complete sufficient, as is its one-to-one function $\overline{X}_n$. Since $E(\overline{X}_n)=\lambda$, the sample mean is the MVUE — an unbiased function $g(T)$ of a complete sufficient statistic $T$ is the UMVUE.

**Example 4.4.16 (Ex 4.37) — normal: MVUE of $(\mu,\sigma^2)$**

$$f(x_1,\dots,x_n;\mu,\sigma^2)=(2\pi\sigma^2)^{-n/2}
\exp\left[-\sum_{i=1}^n x_i^2/2\sigma^2+(\mu/\sigma^2)\sum_{i=1}^n x_i-n\mu^2/2\sigma^2\right]$$

so $S_1=\sum_{i=1}^n X_i$ and $S_2=\sum_{i=1}^n X_i^2$ are jointly complete sufficient, and so is their one-to-one function $(\overline{X}_n,S_n^2)$. Since $E(\overline{X}_n)=\mu$ and $E(S_n^2)=\sigma^2$, the sample mean and sample variance are the MVUEs of $\mu$ and $\sigma^2$.

### 4.4.7 Ancillary statistics and Basu's theorem

**Definition 4.4.6 — ancillary statistic**

A statistic $T$ whose **distribution does not depend on $\theta$**, even though the density does, is called **ancillary** for $\theta$.

**Theorem 4.4.8 (Thm 4.8) — Basu's theorem**

If $S=(S_1,\dots,S_k)$ is jointly **complete sufficient** for $\theta=(\theta_1,\dots,\theta_k)$ and $T$ is ancillary for $\theta$, then $S$ and $T$ are **independent**.

Proof (discrete case). Let $f_T(t)$, $f_S(s;\theta)$ and $f(t\mid s)$ be the densities of $T$, of $S$, and of $T$ given $S=s$. By the law of total probability $f_T(t)=\sum_s f(t\mid s)f_S(s;\theta)$, so

$$E_\theta\big[f_T(t)-f(t\mid S)\big]=f_T(t)-\sum_s f(t\mid s)f_S(s;\theta)=f_T(t)-f_T(t)=0$$

for all $\theta$. Now $f_T(t)-f(t\mid S)$ is a function of $S$: by **ancillarity** $f_T(t)$ is free of $\theta$, and by **sufficiency** $f(t\mid S)$ is free of $\theta$ as well. Hence **completeness** of $S$ forces $f(t\mid s)=f_T(t)$ for all $s$ — the conditional distribution equals the marginal, so $S$ and $T$ are independent. $\blacksquare$

**Example 4.4.17 (Ex 4.38) — $\overline{X}_n\perp\hat\sigma^2$: a second proof of Theorem 3.3.2(1)**

For a random sample from $N(\mu,\sigma^2)$ put $\hat\mu=\overline{X}_n$ and $\hat\sigma^2=\sum_{i=1}^n(X_i-\overline{X}_n)^2/n$. **For fixed $\sigma^2$**, the sample mean is complete sufficient for $\mu$ (Example 4.4.16), while $n\hat\sigma^2/\sigma^2\sim\chi^2(n-1)$ does not depend on $\mu$, so $\hat\sigma^2$ is **ancillary** for $\mu$. Basu's theorem then gives independence of $\overline{X}_n$ and $\hat\sigma^2$ — a different proof from the orthogonal-transformation argument of Theorem 3.3.2(1).
