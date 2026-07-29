# Chapter 1. Probability Theory

## 1.1 Sample Space and Events

**Definition 1.1.1 — experiment, sample space, event**

- **experiment**: a process carried out to obtain an outcome
- **sample space** $S$: the set of all possible outcomes of the experiment
- **event**: a subset of $S$

**Example 1.1.1 (Ex 1.1–1.4) — sample spaces and events**

- three tosses of a coin: $S=\{HHH,HHT,HTH,HTT,THH,THT,TTH,TTT\}$ ($8$ outcomes); "exactly two heads" is $A=\{HHT,HTH,THH\}$
- same experiment, counting heads only: $S=\{0,1,2,3\}$; "at least one head" is $A=\{1,2,3\}$
- toss until the first head: $S=\{H,TH,TTH,TTTH,\dots\}$, an infinite set; "first head on the second toss" is $A=\{TH\}$
- lifetime of a machine in hours: $S=\{x\mid 0\le x<\infty\}$, continuous; "runs at least $20$ hours" is $A=\{x\mid x\ge20\}$

**Definition 1.1.2 (Def 1.1) — intersection, union**

- **intersection** $A\cap B$: the event of the outcomes belonging to both $A$ and $B$
- **union** $A\cup B$: the event of the outcomes belonging to $A$ or $B$

**Definition 1.1.3 (Def 1.2) — mutually exclusive**

$A$ and $B$ are mutually exclusive when they share no outcome:

$$A\cap B=\varnothing$$

**Definition 1.1.4 (Def 1.3) — complement**

$A^c$ is the event of all outcomes of $S$ not in $A$.

**Example 1.1.2 (Ex 1.5–1.7) — operations on events**

- one die, $S=\{1,2,3,4,5,6\}$, $A=\{2,4,6\}$ (even), $B=\{3,6\}$ (multiples of $3$): $\ A\cap B=\{6\}$, $A\cup B=\{2,3,4,6\}$
- one coin toss: head and tail cannot occur together, so they are mutually exclusive
- one die: the complement of "odd" is "even"

**Theorem 1.1.1 — De Morgan's laws**

$$(A\cup B)^c=A^c\cap B^c,\qquad (A\cap B)^c=A^c\cup B^c$$

**Theorem 1.1.2 — distributive laws**

$$A\cap(B\cup C)=(A\cap B)\cup(A\cap C),\qquad A\cup(B\cap C)=(A\cup B)\cap(A\cup C)$$

Proof. For both laws, show that $x$ belongs to the left side if and only if it belongs to the right side (equivalently, check that the two sides describe the same region of a Venn diagram). $\blacksquare$

## 1.2 Definition of Probability

**Definition 1.2.1 — relative frequency definition (Von Mises)**

If $A$ occurs $m$ times in $n$ repetitions of the experiment, its relative frequency is $m/n$, and

$$P(A)=\lim_{n\to\infty}\frac{m}{n}$$

**Definition 1.2.2 — axioms of probability (Kolmogorov)**

A function $P$ on the events of $S$ is a probability if

1. $P(A)\ge0$ for every event $A$
2. $P(S)=1$
3. ($\sigma$-additivity) for mutually exclusive $A_1,A_2,\dots$ with $A_i\cap A_j=\varnothing$, $i\ne j$: $\ P\big(\bigcup_{i=1}^{\infty}A_i\big)=\sum_{i=1}^{\infty}P(A_i)$

**Definition 1.2.3 — probability space**

The triple $(S,\mathcal F,P)$, where $\mathcal F$ is the family of events to which probabilities are assigned. $\mathcal F$ is closed under complement, union and intersection, so $P(A^c)$, $P(A\cup B)$ and the like are always defined and the axioms apply to them.

**Theorem 1.2.1 (Thm 1.1) — basic properties of probability**

1. $P(A^c)=1-P(A)$
2. $P(\varnothing)=0$
3. $A\subset B\ \Rightarrow\ P(A)\le P(B)$
4. $P(A\cup B)=P(A)+P(B)-P(A\cap B)$

Proof. **(1)** $S=A\cup A^c$ with $A\cap A^c=\varnothing$, so

$$1=P(S)=P(A)+P(A^c)\quad\text{(by axioms 2 and 3)}$$

**(2)** $\varnothing=S^c$, so

$$P(\varnothing)=P(S^c)=1-P(S)=0\quad\text{(by (1))}$$

**(3)** $A\subset B$ gives the disjoint decomposition $B=A\cup(B\cap A^c)$, so

$$P(B)=P(A)+P(B\cap A^c)\ \ge\ P(A)\quad\text{(by axiom 3 and }P(B\cap A^c)\ge0\text{, axiom 1)}$$

**(4)** Both $A\cup B=A\cup(B\cap A^c)$ and $B=(A\cap B)\cup(B\cap A^c)$ are disjoint decompositions, so

$$P(A\cup B)=P(A)+P(B\cap A^c),\qquad P(B\cap A^c)=P(B)-P(A\cap B)$$

$$\therefore\ P(A\cup B)=P(A)+P(B)-P(A\cap B)\qquad\blacksquare$$

**Corollary 1.2.1 — consequences of Theorem 1.2.1**

- $0\le P(A)\le1$, since $\varnothing\subset A\subset S$ gives $0=P(\varnothing)\le P(A)\le P(S)=1$ (by monotonicity (3))
- probability of a difference: $\ P(A\cap B^c)=P(A)-P(A\cap B)$, from the disjoint decomposition $A=(A\cap B)\cup(A\cap B^c)$
- Boole's inequality for two events: $\ P(A\cup B)\le P(A)+P(B)$ (by (4) and $P(A\cap B)\ge0$)

**Theorem 1.2.2 — inclusion–exclusion for three events**

$$P(A\cup B\cup C)=P(A)+P(B)+P(C)-P(A\cap B)-P(A\cap C)-P(B\cap C)+P(A\cap B\cap C)$$

Proof (Ex 1.11). Treat $A\cup B$ as a single event and apply Theorem 1.2.1(4) twice:

$$P((A\cup B)\cup C)=P(A\cup B)+P(C)-P((A\cup B)\cap C)$$

$$(A\cup B)\cap C=(A\cap C)\cup(B\cap C)\quad\text{(by the distributive law, Theorem 1.1.2)}$$

$$P((A\cup B)\cap C)=P(A\cap C)+P(B\cap C)-P(A\cap B\cap C)\quad\text{(by Theorem 1.2.1(4))}$$

Substituting $P(A\cup B)=P(A)+P(B)-P(A\cap B)$ and the last line gives the claim. $\blacksquare$

**Theorem 1.2.3 — inclusion–exclusion for $n$ events**

$$P\Big(\bigcup_{i=1}^{n}A_i\Big)=\sum_{i}P(A_i)-\sum_{i<j}P(A_i\cap A_j)+\sum_{i<j<k}P(A_i\cap A_j\cap A_k)-\cdots+(-1)^{n-1}P\Big(\bigcap_{i=1}^{n}A_i\Big)$$

Intersections of an odd number of events enter with $+$ and of an even number with $-$; the term collecting all $m$-fold intersections carries the sign $(-1)^{m-1}$.

Proof (by induction on $n$). The case $n=2$ is Theorem 1.2.1(4) and $n=3$ is Theorem 1.2.2. Assume the formula for $n-1$ events and group the first $n-1$ into $B=A_1\cup\cdots\cup A_{n-1}$, leaving two events:

$$P(A_1\cup\cdots\cup A_n)=P(B\cup A_n)=P(B)+P(A_n)-P(B\cap A_n)\quad\text{(by Theorem 1.2.1(4))}$$

Both remaining pieces are unions of $n-1$ events, so the induction hypothesis expands them:

- $P(B)=P(A_1\cup\cdots\cup A_{n-1})$
- $B\cap A_n=(A_1\cap A_n)\cup\cdots\cup(A_{n-1}\cap A_n)$ (by the distributive law)

Checking the signs: an index set $S$ not containing $n$ comes from $P(B)$ with sign $(-1)^{|S|-1}$; $S=\{n\}$ comes from $+P(A_n)$; and $S=T\cup\{n\}$ with $|S|\ge2$ comes from $-P(B\cap A_n)$ with sign $-(-1)^{|T|-1}=(-1)^{|S|-1}$, since $|S|=|T|+1$. Every $S$ therefore appears exactly once with sign $(-1)^{|S|-1}$. $\blacksquare$

**Theorem 1.2.4 — Boole's inequality**

$$P\Big(\bigcup_{i=1}^{\infty}A_i\Big)\le\sum_{i=1}^{\infty}P(A_i)$$

Proof. Strip off the overlaps to make the events disjoint: $B_1=A_1$ and $B_k=A_k\cap(A_1\cup\cdots\cup A_{k-1})^c$. Then the $B_k$ are mutually exclusive, $B_k\subset A_k$, and $\bigcup_i B_i=\bigcup_i A_i$, so

$$P\Big(\bigcup_i A_i\Big)=P\Big(\bigcup_i B_i\Big)=\sum_i P(B_i)\le\sum_i P(A_i)\quad\text{(by axiom 3, then monotonicity of }B_i\subset A_i)$$

$\blacksquare$

**Theorem 1.2.5 — Bonferroni's inequality**

$$P\Big(\bigcap_{i=1}^{k}A_i\Big)\ge 1-\sum_{i=1}^{k}P(A_i^c)$$

Proof.

$$\Big(\bigcap_{i}A_i\Big)^c=\bigcup_{i}A_i^c\quad\text{(by De Morgan, Theorem 1.1.1)}$$

$$P\Big(\bigcap_{i}A_i\Big)=1-P\Big(\bigcup_{i}A_i^c\Big)\ \ge\ 1-\sum_{i}P(A_i^c)\quad\text{(by Theorem 1.2.4)}$$

$\blacksquare$

When each $A_i$ is nearly certain (each $P(A_i^c)$ small), Bonferroni's inequality bounds from below the probability that they hold simultaneously; the Bonferroni correction for simultaneous tests rests on it.

## 1.3 Conditional Probability

**Definition 1.3.1 (Def 1.4) — conditional probability**

For $P(B)>0$, the probability of $A$ given that $B$ has occurred is

$$P(A\mid B)=\frac{P(A\cap B)}{P(B)}$$

**Example 1.3.1 (Ex 1.12) — computing a conditional probability**

Two tosses, $S=\{HH,HT,TH,TT\}$, $H=\{HH,HT\}$ (first toss a head), $A=\{HH\}$:

$$P(A\mid H)=\frac{P(A\cap H)}{P(H)}=\frac{1/4}{1/2}=\frac12$$

**Theorem 1.3.1 — $P(\,\cdot\mid B)$ is itself a probability**

With $B$ fixed, $P(\,\cdot\mid B)$ satisfies the three axioms, so all properties of Theorem 1.2.1 carry over:

1. $P(A\mid B)\ge0$, $\ P(S\mid B)=P(B\mid B)=1$
2. $P(A^c\mid B)=1-P(A\mid B)$, $\ P(\varnothing\mid B)=0$
3. $A_1\subset A_2\ \Rightarrow\ P(A_1\mid B)\le P(A_2\mid B)$
4. $P(A_1\cup A_2\mid B)=P(A_1\mid B)+P(A_2\mid B)-P(A_1\cap A_2\mid B)$

Proof. For additivity, $A_1\cap A_2=\varnothing$ implies $(A_1\cap B)\cap(A_2\cap B)=\varnothing$, so

$$P(A_1\cup A_2\mid B)=\frac{P[(A_1\cup A_2)\cap B]}{P(B)}=\frac{P[(A_1\cap B)\cup(A_2\cap B)]}{P(B)}=\frac{P(A_1\cap B)+P(A_2\cap B)}{P(B)}=P(A_1\mid B)+P(A_2\mid B)$$

(by axiom 3). Also $P(A\mid B)\ge0$ and $P(S\mid B)=P(S\cap B)/P(B)=P(B)/P(B)=1$ follow from the definition. All three axioms hold, so the remaining properties follow by the proof of Theorem 1.2.1. $\blacksquare$

**Theorem 1.3.2 — multiplication rule**

For $P(B)>0$, and in general for $k$ events,

$$P(A\cap B)=P(B)\,P(A\mid B)$$

$$P(A_1\cap A_2\cap\cdots\cap A_k)=P(A_1)\,P(A_2\mid A_1)\,P(A_3\mid A_1\cap A_2)\cdots P(A_k\mid A_1\cap\cdots\cap A_{k-1})$$

**Example 1.3.2 (Ex 1.15) — urn, sampling without replacement**

$1$ red and $9$ white balls; white on each of the first three draws and red on the fourth:

$$P(W_1\cap W_2\cap W_3\cap R_4)=\frac{9}{10}\cdot\frac{8}{9}\cdot\frac{7}{8}\cdot\frac{1}{7}=\frac{1}{10}\quad\text{(by Theorem 1.3.2)}$$

Each draw's outcome becomes the condition for the next, so numerator and denominator shrink in turn.

**Definition 1.3.2 — partition**

$B_1,\dots,B_k$ partition $S$ when they are mutually exclusive ($B_i\cap B_j=\varnothing$, $i\ne j$) and $\bigcup_{i=1}^k B_i=S$.

**Theorem 1.3.3 (Thm 1.2) — law of total probability**

If $B_1,\dots,B_k$ partition $S$, then for any event $A$

$$P(A)=\sum_{i=1}^{k}P(A\cap B_i)=\sum_{i=1}^{k}P(B_i)\,P(A\mid B_i)$$

Proof. Since $\bigcup B_i=S$,

$$A=A\cap S=A\cap\Big(\bigcup_{i=1}^k B_i\Big)=\bigcup_{i=1}^k(A\cap B_i)$$

and the $A\cap B_i$ are mutually exclusive because the $B_i$ are, so

$$P(A)=\sum_{i=1}^{k}P(A\cap B_i)=\sum_{i=1}^{k}P(B_i)\,P(A\mid B_i)\quad\text{(by axiom 3, then Theorem 1.3.2)}$$

$\blacksquare$

**Theorem 1.3.4 (Thm 1.3) — Bayes' theorem**

If $B_1,\dots,B_k$ partition $S$ and $A$ has occurred, the probability that it came from $B_j$ is

$$P(B_j\mid A)=\frac{P(A\cap B_j)}{P(A)}=\frac{P(B_j)\,P(A\mid B_j)}{\sum_{i=1}^{k}P(B_i)\,P(A\mid B_i)}$$

Proof. Start from the definition of conditional probability, then replace the numerator by the multiplication rule (Theorem 1.3.2) and the denominator by the law of total probability (Theorem 1.3.3). $\blacksquare$

**Example 1.3.3 (Ex 1.16) — total probability**

Three production lines make $50\%,30\%,20\%$ of the output with defect rates $2\%,5\%,10\%$:

$$P(D)=0.5(0.02)+0.3(0.05)+0.2(0.10)=0.045\quad\text{(by Theorem 1.3.3)}$$

**Example 1.3.4 (Ex 1.17) — Bayes' theorem**

Prevalence $P(A)=0.01$, sensitivity $P(+\mid A)=0.95$, false positive rate $P(+\mid A^c)=0.10$:

$$P(A\mid+)=\frac{0.95(0.01)}{0.95(0.01)+0.10(0.99)}\approx0.088\quad\text{(by Theorem 1.3.4)}$$

**Definition 1.3.3 (Def 1.5) — independence of two events**

$A$ and $B$ are independent when

$$P(A\cap B)=P(A)\,P(B)$$

For $P(B)>0$ this is equivalent to $P(A\mid B)=P(A)$: knowing $B$ does not change the probability of $A$.

**Definition 1.3.4 (Def 1.6) — mutual independence**

$A_1,\dots,A_n$ are mutually independent when, for every choice of $k$ indices $i_1,\dots,i_k$ out of $\{1,\dots,n\}$ and every $k=2,\dots,n$,

$$P(A_{i_1}\cap A_{i_2}\cap\cdots\cap A_{i_k})=P(A_{i_1})P(A_{i_2})\cdots P(A_{i_k})$$

This is stronger than pairwise independence, which checks only $k=2$.

**Example 1.3.5 (Ex 1.19, 1.22) — independence, pairwise independence**

- one card from $52$, $K$ = king, $D$ = diamond: $\ P(K\cap D)=\frac1{52}=\frac4{52}\cdot\frac{13}{52}=P(K)P(D)$, so $K$ and $D$ are independent
- pairwise but not mutually independent: $\ P(A\cap B\cap C)=\frac28=\frac14\ne P(A)P(B)P(C)=\frac18$, although all pairs factor

## 1.4 Counting

**Theorem 1.4.1 — multiplication rule for counting**

If $A$ can occur in $m$ ways and $B$ in $n$ ways, then $A$ followed by $B$ can occur in $m\times n$ ways.

**Definition 1.4.1 — permutation**

The number of ways of drawing $r$ items from $n$ distinct items with order:

$${}_nP_r=n(n-1)(n-2)\cdots(n-r+1)=\frac{n!}{(n-r)!},\qquad 0!=1$$

**Definition 1.4.2 — combination**

The number of ways of drawing $r$ from $n$ distinct items without order; each combination corresponds to $r!$ orderings, so

$$\binom{n}{r}={}_nC_r=\frac{{}_nP_r}{r!}=\frac{n!}{r!\,(n-r)!},\qquad \binom{n}{r}=\binom{n}{n-r}$$

**Theorem 1.4.2 — Pascal's identity**

$$\binom{n}{r}=\binom{n-1}{r-1}+\binom{n-1}{r},\qquad 1\le r\le n-1$$

Proof (combinatorial). Fix one element $a$ and split the selections according to whether $a$ is drawn.

- $a$ included: choose the remaining $r-1$ from the other $n-1$ items, $\binom{n-1}{r-1}$ ways
- $a$ excluded: choose all $r$ from the other $n-1$ items, $\binom{n-1}{r}$ ways

The two cases are mutually exclusive and together exhaust all selections, so the counts add. $\blacksquare$

This identity is the rule behind Pascal's triangle: each entry is the sum of the two above it.

**Theorem 1.4.3 — binomial theorem**

$$(a+b)^n=\sum_{k=0}^{n}\binom{n}{k}a^{k}b^{n-k}$$

Proof. Expanding the $n$ factors of $(a+b)(a+b)\cdots(a+b)$ picks $a$ or $b$ from each factor. A term $a^kb^{n-k}$ arises by choosing which $k$ of the $n$ factors contribute $a$, and there are $\binom nk$ such choices. $\blacksquare$

Putting $a=b=1$ gives $\sum_{k=0}^n\binom nk=2^n$, the number of subsets of an $n$-element set.

**Theorem 1.4.4 — multinomial theorem**

With the multinomial coefficient, for $r_1+\cdots+r_k=n$,

$$\binom{n}{r_1\,r_2\,\cdots\,r_k}=\frac{n!}{r_1!\,r_2!\cdots r_k!}$$

$$(a_1+a_2+\cdots+a_k)^n=\sum_{r_1+\cdots+r_k=n}\binom{n}{r_1\,r_2\,\cdots\,r_k}a_1^{r_1}a_2^{r_2}\cdots a_k^{r_k}$$

Proof. Each of the $n$ factors contributes one of $a_1,\dots,a_k$. The term $a_1^{r_1}\cdots a_k^{r_k}$ arises as often as the $n$ factors can be split into groups of sizes $r_1,\dots,r_k$, and that count is the multinomial coefficient. The binomial theorem is the case $k=2$. $\blacksquare$

**Definition 1.4.3 — sampling with and without replacement**

- **with replacement**: a drawn item may be drawn again
- **without replacement**: a drawn item is removed from the pool

The count also depends on whether order is taken into account, so both choices must be fixed before counting.

**Example 1.4.1 (Ex 1.25–1.28) — counting**

- $n$ balls into $n$ baskets at random: $n^n$ arrangements in all, $n!$ with one ball per basket, so the probability is $n!/n^n$
- $k$ balls ($n>k$) with at most one per basket: $\ {}_nP_k/n^k$
- $5$ white and $10$ black balls, $3$ drawn without replacement: counted by combinations, e.g. all black with probability $\binom{10}{3}\big/\binom{15}{3}$
- $5$ cards from $52$, a straight: $\ \dfrac{10\,(4^5-4)}{\binom{52}{5}}\approx0.0039$
