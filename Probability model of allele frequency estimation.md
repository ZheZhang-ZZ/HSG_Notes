### Probability model of allele frequency estimation

#### Conditional independence assumptions

Example 1.1 is a good start point to understand the idea of Bayesian statistics and basic assumptions of conditional independence. Here I will try to derive formula (1.1), which is very important for the following EM algorithm section.

Consider that we have sequenced $n$ individuals at a specific locus, and the sequece data can be represented as ${\boldsymbol \rm d}=\{d_1,\dots,d_n\}$. Here we have a unknown parameter vector ($\boldsymbol \theta$) comprised of allele frequency $f$, sequencing error $\epsilon$, and we can construct a probability model based on joint probability of $\{{\boldsymbol \rm d}, \boldsymbol \theta\}$. Based on Bayes formula, we can get the following formula

$$
\mathbb{P}({\boldsymbol \rm d}, \boldsymbol \theta)=\mathbb{P}({\boldsymbol \rm d}|\boldsymbol \theta)\mathbb{P}(\boldsymbol \theta) \tag{1}
$$

where $\boldsymbol \theta = \{f, \epsilon\}$, substitute this into formula (1), we get

$$
\mathbb{P}({\boldsymbol \rm d}, f, \epsilon)=\mathbb{P}({\boldsymbol \rm d}|f, \epsilon)\mathbb{P}(f, \epsilon) \tag{2}
$$

Formula (2) is just the first part of formula (1.1) in the book that as was shown below:

$$
\mathbb{P}({\boldsymbol \rm d}, f, \epsilon)=\mathbb{P}({\boldsymbol \rm d}|f, \epsilon)\mathbb{P}(f, \epsilon)=\left[\prod_{i=1}^{n}\sum_{g_{i}}\mathbb{P}(d_i|g_i,\epsilon)\mathbb{P}(g_i|f)\right]\mathbb{P}(f, \epsilon)
$$

It is a little bit confusing that how to convert $\mathbb{P}({\boldsymbol \rm d}|f, \epsilon)$ to $\left[\prod_{i=1}^{n}\sum_{g_{i}}\mathbb{P}(d_i|g_i,\epsilon)\mathbb{P}(g_i|f)\right]$, and we will derive it as follows:

Firstly, because the sequecing for different individuals are independent, so we can get the following formula based on the product role of probabilities:

$$
\mathbb{P}({\boldsymbol \rm d}|f, \epsilon)=\prod_{i=1}^{n}\mathbb{P}(d_{i}|f, \epsilon) \tag{3}
$$

For each individual $i$, we can now put a new variable $g_i$, which is the unobserved genotype for this individual, in $\mathbb{P}(d_{i}|f, \epsilon)$, and then we can just get rid of this unobservable variable by integral, and this trick is just how we get the marginal distribution of a variable (eg. consider two variables $\{{\boldsymbol \rm x},{\boldsymbol \rm y}\}$ and we can get the marginal probability of ${\boldsymbol \rm x}$ by the integral computation: $\mathbb{P}({\boldsymbol \rm x})=\int \mathbb{P}({\boldsymbol \rm x}, {\boldsymbol \rm y}){\mathrm{d}{\boldsymbol \rm x}}$). As integral is for continuous variable and $g$ is actually a discrete varaible, so "integral" for $g$ is just summation. For one individual, the formula which can be presented below:

$$
\begin{aligned}
\mathbb{P}(d_{i}|f, \epsilon)&=\sum_{g_i}\mathbb{P}(d_{i},g_i|f, \epsilon) \\
&=\sum_{g_i}\mathbb{P}(d_{i}|g_i,f, \epsilon)\mathbb{P}(g_i|f,\epsilon)
\end{aligned} \tag{4}
$$

Formula (4) is still not the same as the part after summation symbol in the square brackets in formula (1.1). In fact, there is a trick here called **conditonal independence** assumed here. When $g_i$ has been known $d_i$ and $f$ is independent, which means that $\mathbb{P}(d_{i}|g_i,f, \epsilon)=\mathbb{P}(d_{i}|g_i,\epsilon)$, and $\mathbb{P}(g_i|f,\epsilon)=\mathbb{P}(g_i|f)$, which is due to that $g_i$ and $\epsilon$ is independent based on the assumption of the scenario. Therefore, formula can be converted to the following formula:

$$
\sum_{g_i}\mathbb{P}(d_{i}|g_i,f, \epsilon)\mathbb{P}(g_i|f,\epsilon) = \sum_{g_i}\mathbb{P}(d_{i}|g_i,\epsilon)\mathbb{P}(g_i|f) \tag{5}
$$

Now, substitute (5) into (3), we get:

$$
\mathbb{P}({\boldsymbol \rm d}|f, \epsilon)=\prod_{i=1}^{n}\sum_{g_i}\mathbb{P}(d_{i}|g_i,\epsilon)\mathbb{P}(g_i|f) \tag{6}
$$

Again, substitute (6) into (2), we get the same formula (1.1) in the book as follows:

$$
\mathbb{P}({\boldsymbol \rm d}, f, \epsilon)=\mathbb{P}({\boldsymbol \rm d}|f, \epsilon)\mathbb{P}(f, \epsilon)=\left[\prod_{i=1}^{n}\sum_{g_{i}}\mathbb{P}(d_i|g_i,\epsilon)\mathbb{P}(g_i|f)\right]\mathbb{P}(f, \epsilon) \tag{7}
$$

If we look into the equation (7) from Bayesian view, we can get that $\mathbb{P}(f, \epsilon)$ is the **prior probability** of $\boldsymbol \theta$, $\left[\prod_{i=1}^{n}\sum_{g_{i}}\mathbb{P}(d_i|g_i,\epsilon)\mathbb{P}(g_i|f)\right]$ is called **likelihood** and $\mathbb{P}({\boldsymbol \rm d}, f, \epsilon)$ is the **posterior probability**. 

#### Probability Distributions

For the likelihood in equation (7), we need to specified two probability funcitons, $\mathbb{P}(d_i|g_i,\epsilon)$ and $\mathbb{P}(g_i|f)$. Firstly, if we known $f$，then we can estimate the probability of  genotype using H-W equilibrium:

$$
\mathbb{P}(g_i|f)=\tbinom{2}{g_i}f^{g_i}(1-f)^{g_i}=\left\{
\begin{aligned}
&(1-f)^2  &\text{ if } g_i=0 \text{,}\\
&2f(1-f)  &\text{ if } g_i=1 \text{,}\\
&f^2      &\text{ if } g_i=2 \text{,}
\end{aligned}
\right. \tag{8}
$$
where you can see that the piecewise function was converted to a single-line equation nicely. The function of $\mathbb{P}(d_i|g_i,\epsilon)$ can be represented below:
$$
\mathbb{P}(d_i|g_i,\epsilon)=\prod_{j=1}^{n_i}\mathbb{P}(b_{ij}|g_i,\epsilon) \tag{9}
$$
where  $b_{ij}$ is the *j*th read for *i*th individual. We can obtain this equation because we assume that reads for each individual were independent from each other, and $\mathbb{P}(b_{ij}|g_i,\epsilon)$ can be further converted to equation (10):
$$
\mathbb{P}(b_{ij}|g_i,\epsilon)=\left\{
\begin{aligned}
&(1-\epsilon)  &\text{ if } g_i=0, b_{ij}=\text{A} \text{ or } g_i=2, b_{ij}=\text{D} \text{,}\\
&\frac{1}{2}  &\text{ if } g_i=1 \text{,}\\
&\epsilon      &\text{ if } g_i=0, b_{ij}=\text{D} \text{ or } g_i=2, b_{ij}=\text{A} \text{,}
\end{aligned}
\right. \tag{10}
$$
wehre **A** is the ancestral allele (which can just be denoted as "REF" allele in VCF file), and **D** is the derived allele (which can just be denoted as "ALT" allele in VCF file). For example, if we get a sequence result of first indiviual as $d_1=\{\text{A, D, A}\}$, and if the acual genotype of individual is "0", this means the sequence result of **D** is an error, because we cannot get a derive allele from two ancestral alleles. Therefore, we can calculate the genotype likelihood for this individual as $\mathbb{P}(d_1|g_1=0,\epsilon)=(1-\epsilon)*\epsilon*(1-\epsilon)=(1-\epsilon)^2\epsilon$

There is actually a prior prbability $\mathbb{P}(f, \epsilon)$ in equation 7. However, if we don't want to estimate parameters using Bayesian inference, we can ignore it at present.

### An introduction to EM algorithm

In this section, we will firstly make a simple introduction to EM algorithm. Please note that we will not dive into the whole mathematical theory, but highlight some key concepts to help you understand this classic algorithm and to make a foundation for its application in estimation of allele frequencies.

#### Definition

EM (expectation-maximization) algorithm was proposed to handle the maximum likelihood estimation (MLE) of parametres in likelihood function with latent variables (which is genotype in our example). Below is the formula of EM:
$$
\boldsymbol \theta^{(t+1)}=\mathop{\arg\max}_{\boldsymbol \theta}\int_Z\left[\log(\mathbb{P}(\boldsymbol {X,Z|\theta}))\right]\mathbb{P}(\boldsymbol {Z|X,\theta^{(t)}}) \tag{11}
$$
Let's explain a little bit of equation 11. Firstly, EM algorithm is a iteration procedure, so $t$ means current iteration and $t+1$ means the next iteration. Our aim is to find MLE of parameters for the next round of iteration until the procedure converges. Secondly,  $\mathbb{P}(\boldsymbol {Z|X,\theta^{(t)}})$ is called posterior probability of $\boldsymbol {\rm Z}$, where $\boldsymbol {\rm Z}$ is the latent variable, $\boldsymbol {\rm X}$ is the observaed data, and $\boldsymbol {\rm \theta}$ is the unknown parameters. $\log(\mathbb{P}(\boldsymbol {X,Z|\theta}))$ is called log likelihood of complete data, and $\int_Z\left[\log(\mathbb{P}(\boldsymbol {X,Z|\theta}))\right]\mathbb{P}(\boldsymbol {Z|X,\theta^{(t)}})$ is the expectation of $\log(\mathbb{P}(\boldsymbol {X,Z|\theta}))$ in the space of $\boldsymbol {\rm Z}$. Why is it a expectation?

We all know that for a random variable of $\boldsymbol {\rm Z}$, the expectation of $\boldsymbol {\rm Z}$ is $\mathbb{E}\left[\boldsymbol {\rm Z}\right]=\sum{z\mathbb{P}(z)}$, when $\boldsymbol {\rm Z}$ is a discrete variable and $\mathbb{E}\left[\boldsymbol {\rm Z}\right]=\int{z\mathbb{P}(z)}\mathrm{d}z$, when $\boldsymbol {\rm Z}$ is a continuous variable. The regulation here is that the integral of a random variable multiplies the pdf (pmf) is actually the expectation. Therefore, $\int_Z\left[\log(\mathbb{P}(\boldsymbol {X,Z|\theta}))\right]\mathbb{P}(\boldsymbol {Z|X,\theta^{(t)}})$ can be regarded as $\mathbb{E}_{\boldsymbol Z}\left[\log(\mathbb{P}(\boldsymbol {X,Z|\theta}))\right]$, which is referred as E-step in EM algorithm. Meanwhile, $\mathop{\arg\max}_{\boldsymbol \theta}\mathbb{E}_{\boldsymbol Z}\left[\log(\mathbb{P}(\boldsymbol {X,Z|\theta}))\right]$ is called M-step because it is used to maximize the expectation.

#### Proof of convergence

In the book, the authors listed the simple proof of the convergence of EM algorithm, which is not clear. Here we will give the details of the proof:

The convergence of EM means that the next round of log likelihood of the data given parameters are always larger than the current those of current round, i.e., $\log{\mathbb{P}(\boldsymbol X|\boldsymbol \theta^{(t+1)})}>\log{\mathbb{P}(\boldsymbol X|\boldsymbol \theta^{(t)})}$.

Based on Bayes rule, we get:
$$
\log{\mathbb{P}(\boldsymbol x|\boldsymbol\theta)}=\log{\mathbb{P}(\boldsymbol {x,z}|\boldsymbol\theta)} - \log{\mathbb{P}(\boldsymbol z|\boldsymbol{x,\theta)}} \tag{12}
$$
Now let us calculate the Expectation for the left and righ hand side based on $\boldsymbol {Z|X,\theta^{(t)}}$, for the left hand side:
$$
\begin{aligned}
\mathbb{E}_{z|x,\theta^{(t)}}\left[\log\mathbb{P}(\boldsymbol {x|\theta})\right]&=\int_{\boldsymbol z}\left[\log(\mathbb{P}(\boldsymbol {x|\theta}))\right]\mathbb{P}(\boldsymbol {z|x,\theta^{(t)}})\mathrm{d}\boldsymbol z \\
&=\log(\mathbb{P}(\boldsymbol {x|\theta}))\int_{\boldsymbol z}\mathbb{P}(\boldsymbol {z|x,\theta^{(t)}})\mathrm{d}\boldsymbol z \\
&=\log(\mathbb{P}(\boldsymbol {x|\theta})) \times 1 \\
&= \log(\mathbb{P}(\boldsymbol {x|\theta}))
\end{aligned} \tag{13}
$$
The first line of eq.13 is the definition of Expectation. The second line holds because $\log(\mathbb{P}(\boldsymbol {x|\theta}))$ is not related to $\boldsymbol Z$ and can be taken outside the $\int$ directly. The third line holds because $\int_{\boldsymbol z}\mathbb{P}(\boldsymbol {z|x,\theta^{(t)}})\mathrm{d}\boldsymbol z$ is just the whole area of a probability density function (pdf), which is 1 from the definition. 

For the right hand side, we get:
$$
\begin{aligned}
\mathbb{E}_{z|x,\theta^{(t)}}\left[\log\mathbb{P}(\boldsymbol {x,z|\theta})\right]-\mathbb{E}_{z|x,\theta^{(t)}}\left[\log\mathbb{P}(\boldsymbol {z|x,\theta})\right]&=\int_{\boldsymbol z}\left[\log(\mathbb{P}(\boldsymbol {x,z|\theta}))\right]\mathbb{P}(\boldsymbol {z|x,\theta^{(t)}})\mathrm{d}\boldsymbol z - \\&\int_{\boldsymbol z}\left[\log(\mathbb{P}(\boldsymbol {z|x,\theta}))\right]\mathbb{P}(\boldsymbol {z|x,\theta^{(t)}})\mathrm{d}\boldsymbol z\\
\end{aligned} \tag{14}
$$
Let's define $Q(\boldsymbol {\theta,\theta^{(t)}})=\int_{\boldsymbol z}\left[\log(\mathbb{P}(\boldsymbol {x,z|\theta}))\right]\mathbb{P}(\boldsymbol {z|x,\theta^{(t)}})\mathrm{d}\boldsymbol z$, and $H(\boldsymbol {\theta,\theta^{(t)}})=\int_{\boldsymbol z}\left[\log(\mathbb{P}(\boldsymbol {z|x,\theta}))\right]\mathbb{P}(\boldsymbol {z|x,\theta^{(t)}})\mathrm{d}\boldsymbol z$

Now, let us recall the aim of this proof, i.e., $\log{\mathbb{P}(\boldsymbol X|\boldsymbol \theta^{(t+1)})}>\log{\mathbb{P}(\boldsymbol X|\boldsymbol \theta^{(t)})}$. If we can prove that $Q(\boldsymbol {\theta^{(t+1)},\theta^{(t)}})\geq Q(\boldsymbol {\theta,\theta^{(t)}}) \text{ && } H(\boldsymbol {\theta^{(t+1)},\theta^{(t)}})\leq H(\boldsymbol {\theta,\theta^{(t)}})$, then it means $\log{\mathbb{P}(\boldsymbol X|\boldsymbol \theta^{(t+1)})}>\log{\mathbb{P}(\boldsymbol X|\boldsymbol \theta^{(t)})}$ holds.

Firstly, we can find that $Q$ function is acutually the expectation in eq. 11 (EM algorithm) that we want to maximize. Now that we are maximizing $Q$ function during iteration, so the next round of $Q$ is must be larger or equal than the current round, which means that $Q(\boldsymbol {\theta^{(t+1)},\theta^{(t)}})\geq Q(\boldsymbol {\theta,\theta^{(t)}})$ holds. 

Secondly, we need to prove $H(\boldsymbol {\theta^{(t+1)},\theta^{(t)}})\leq H(\boldsymbol {\theta,\theta^{(t)}})$, i.e., $H(\boldsymbol {\theta^{(t+1)},\theta^{(t)}}) - H(\boldsymbol {\theta^{(t)},\theta^{(t)}}) \leq0$:
$$
\begin{aligned}
&H(\boldsymbol {\theta^{(t+1)},\theta^{(t)}}) - H(\boldsymbol {\theta^{(t)},\theta^{(t)}}) \\
&=\int_{\boldsymbol z}\left[\log(\mathbb{P}(\boldsymbol {z|x,\theta^{(t+1)}}))\right]\mathbb{P}(\boldsymbol {z|x,\theta^{(t)}})\mathrm{d}\boldsymbol z-\int_{\boldsymbol z}\left[\log(\mathbb{P}(\boldsymbol {z|x,\theta^{(t)}}))\right]\mathbb{P}(\boldsymbol {z|x,\theta^{(t)}})\mathrm{d}\boldsymbol z \\
&=\int_{\boldsymbol z}\mathbb{P}(\boldsymbol {z|x,\theta^{(t)}})\left[\frac{\log(\mathbb{P}(\boldsymbol {z|x,\theta^{(t+1)}}))}{\log(\mathbb{P}(\boldsymbol {z|x,\theta^{(t)}}))}\right]\mathrm{d}\boldsymbol z = \mathbb{E}_{z|x,\theta^{(t)}}\left[\frac{\log(\mathbb{P}(\boldsymbol {z|x,\theta^{(t+1)}}))}{\log(\mathbb{P}(\boldsymbol {z|x,\theta^{(t)}}))}\right]\\
&=\mathbb{E}_{z|x,\theta^{(t)}}\log\left[\frac{(\mathbb{P}(\boldsymbol {z|x,\theta^{(t+1)}}))}{(\mathbb{P}(\boldsymbol {z|x,\theta^{(t)}}))}\right]
\end{aligned} \tag{15}
$$
Here we need a famous inequality, called *Jensen's inequality* to help complete the proof. Jensen inequality can be easily illustrated in a figure:

![Jensen equlity](/Users/zhezhang/Desktop/icloud/0-tongbu/2-learning/statistical_genomics_and_bioinformatics/hsg/Jensen_equality.png)

This figure means that for a concave function, e.g., log function, $\mathbb{E}[f(x)]\le f(\mathbb{E}[x])$, which can be observed in the figure above that the line is always below the curve. With Jensen's inequality, we can obtain:
$$
\begin{aligned}
\mathbb{E}[f(x)]=\mathbb{E}_{z|x,\theta^{(t)}}\log\left[\frac{(\mathbb{P}(\boldsymbol {z|x,\theta^{(t+1)}}))}{(\mathbb{P}(\boldsymbol {z|x,\theta^{(t)}}))}\right] \le f(\mathbb{E}[x])&=\log{\mathbb{E}_{z|x,\theta^{(t)}}}\left[\frac{(\mathbb{P}(\boldsymbol {z|x,\theta^{(t+1)}}))}{(\mathbb{P}(\boldsymbol {z|x,\theta^{(t)}}))}\right]\\
&=\log\int_{z}\mathbb{P}(\boldsymbol {z|x,\theta^{(t)}}))\frac{(\mathbb{P}(\boldsymbol {z|x,\theta^{(t+1)}}))}{(\mathbb{P}(\boldsymbol {z|x,\theta^{(t)}}))}\mathrm{d}z\\
&=\log\int_{z}\mathbb{P}(\boldsymbol {z|x,\theta^{(t)}}))\mathrm{d}z\\
&=\log1=0
\end{aligned}
$$
So, $H(\boldsymbol {\theta^{(t+1)},\theta^{(t)}}) - H(\boldsymbol {\theta^{(t)},\theta^{(t)}}) \leq0$ holds, and hence $\log{\mathbb{P}(\boldsymbol X|\boldsymbol \theta^{(t+1)})}>\log{\mathbb{P}(\boldsymbol X|\boldsymbol \theta^{(t)})}$, i.e., EM algorithm can always converge.

#### Steps of EM

Eq.11 has actually illustrated the two steps included in the implementation of EM:

+ E-step: Given the current $\boldsymbol \theta^{(t)}$, calculate the expecatation, i.e., $Q(\boldsymbol {\theta,\theta^{(t)}})=\int_{\boldsymbol z}\left[\log(\mathbb{P}(\boldsymbol {x,z|\theta}))\right]\mathbb{P}(\boldsymbol {z|x,\theta^{(t)}})\mathrm{d}\boldsymbol z$

+ M-step: Optimize $\boldsymbol \theta$ by maximizing the $Q$ function, and then take $\boldsymbol \theta$ of this round to the next round.

### The example of allele frequency estimation using EM

Based on all the knowledge above, we can finally turn to the allele frequency estimation by EM algorithm. Firstly, from eq. 7, we have the likelihood:
$$
\mathcal{L(f,\epsilon)}=\mathbb{P}(\boldsymbol d|f,\epsilon)=\prod_{i=1}^{n}\sum_{g_{i}=0}^2\mathbb{P}(d_i|g_i,\epsilon)\mathbb{P}(g_i|f) \tag{16}
$$
In EM we need the complete data likelihood with latent varible, which is genotype ($\boldsymbol g$) here:
$$
\mathcal{L}_c(f,\epsilon)=\mathbb{P}(\boldsymbol d, \boldsymbol g|f,\epsilon)=\prod_{i=1}^n\mathbb{P}(d_{i}|g_i,\epsilon)\mathbb{P}(g_i|f) \tag{17}
$$
which holds as we have explained in conditional independence section. Convert it to log likelihood:
$$
\mathscr{L}_c(f,\epsilon)=\sum_{i=1}^n(\log\mathbb{P}(d_{i}|g_i,\epsilon)+\log\mathbb{P}(g_i|f)) \tag{18}
$$

+ E-step:

$$
\begin{aligned}
Q(\boldsymbol {\theta,\theta^{(t)}})&=\int_{\boldsymbol g}\left[\log(\mathbb{P}(\boldsymbol {d,g|\theta}))\right]\mathbb{P}(\boldsymbol {g|d,\theta^{(t)}})\mathrm{d}\boldsymbol g\\
&=\sum_{g=1}^2\sum_{i=1}^n(\log\mathbb{P}(d_{i}|g_i,\epsilon)+\log\mathbb{P}(g_i|f))\mathbb{P}(g_i|d_i,f^{(t)},\epsilon^{(t)}) 

\end{aligned} \tag{19}
$$

where
$$
\mathbb{P}(g_i|d_i,f^{(t)},\epsilon^{(t)})=p_{g_i}=\frac{\mathbb{P}(d_{i}|g_i,\epsilon^{(t)})\mathbb{P}(g_i|f^{(t)})}{\sum_{h=0}^2\mathbb{P}(d_{i}|h,\epsilon^{(t)})\mathbb{P}(h|f^{(t)})} \tag{20}
$$
which can also derived from the *Bayes formula* and *law of total probability*. Each part of eq. 20 has already been illustrated in eq. 8~10. 

+ M-step:

The real charm of EM shows in this step, becasuse in eq.19 we can find $f$ and $\epsilon$ are separated in different terms, which means we can calculate the derivations of them separately:

