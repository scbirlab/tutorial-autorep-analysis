## Transcript dynamics in growing cells

Let $M_i$ be the total number of reporter transcripts from strain $i$ in the pool, and $n_i$ the 
number of cells. Transcripts are made at rate $E_i$ per cell and decay at rate $d$:

$$
\frac{dM_i}{dt} = E_i n_i - d M_i
$$

This is awkward to work with because $M_i$ also grows whenever the population grows. We want the 
transcript level **per cell**, $m_i = M_i / n_i$. Differentiating $M_i / n_i$ with the quotient rule, and substituting $\frac{dM_i}{dt}$ from above, 
$$
\frac{d}{dt} \frac{M_i}{n_i} = \frac{n_i \left( E_i n_i - d M_i \right) - M_i \frac{dn_i}{dt}}{n_i^2} = E_i - d \, m_i - \left( \frac{1}{n_i} \frac{dn_i}{dt} \right) m_i
$$

The bracketed term is the instantaneous per-capita growth rate of the strain. Call it $r_i$:

$$
r_i = \frac{1}{n_i} \frac{dn_i}{dt} = \frac{d \log n_i}{dt}
$$

Then

$$
\frac{dm_i}{dt} = E_i - (d + r_i) m_i
$$

The extra term $r_i m_i$ is dilution. Every division splits a cell's transcripts between the two 
daughters, so growth removes transcripts per cell in the same way decay does.

- $m_i(t)$: reporter transcripts per cell in strain $i$.
- $d$: mRNA decay rate. For bacterial mRNA, half-lives of 2–5 minutes correspond to $d \approx$ 0.15–0.35 per minute.
- $r_i(t)$: instantaneous per-capita growth rate of strain $i$. In exponential phase it equals the intrinsic growth rate $w_i$ of the previous tutorial. In the Lotka–Volterra model it is $w_i$ multiplied by 
$\left(1 - \frac{\Sigma_j n_j}{K}\right)$, a factor shared by all strains. It can be measured from the DNA counts as $d \log n_i / dt$, with the spike-in putting $n_i(t)$ on an absolute scale (see the previous tutorial).

This equation has a closed-form solution if $r_i$ is constant over the relaxation. If $m_i$ starts at $m_i(0)$,

$$
m_i(t) = m_i^* + \left( m_i(0) - m_i^* \right) e^{-(d + r_i) t}, \qquad m_i^* = \frac{E_i}{d + r_i}
$$

The transcript level relaxes to its steady state $m_i^*$ at rate $d + r_i$. With a decay rate 
of 0.15–0.35 per minute, this takes a few minutes.

### Steady state

Samples taken tens of minutes or hours after induction are on the plateau, so we can use 
$m_i = m_i^*$ (the quasi-steady state, QSS). This holds when $r_i(t)$ changes slowly compared with the relaxation rate $d + r_i$, so that $m_i$ tracks its instantaneous steady state $E_i / (d + r_i(t))$. Below, $r_i$ means its value at the sampling time. Compare strain $i$ to a reference strain 0 with no repression 
($K_0 = 0$, e.g. a non-targeting sgRNA):

$$
\frac{m_i}{m_0} = \frac{d + r_0}{d + r_i} \cdot \frac{1}{1 + K_i R}
$$

The first factor is a growth term, evaluated with the instantaneous growth rates at the sampling time. It is 1 when $r_i = r_0$, and then

$$
\frac{m_i}{m_0} = \frac{1}{1 + K_i R}
$$

So transcript level per cell falls with sgRNA strength, with a known shape, and no time 
dependence under steady state (late timepoint) consitions. 

We assume the growth term is 1 for the analysis below. This holds near carrying capacity, where every strain 
grows slowly and similarly. It fails for guides that have a strong effect on an essential chromosomal gene. Those strains 
grow slower, dilute their transcripts less, and look weaker than they are. 

It might seem better to fit the time course after induction. But dCas9 
does not switch on instantly: $R(t)$ rises over a timescale comparable to or slower than mRNA relaxation. 
The early transient therefore mostly reports how fast dCas9 accumulates, not $d$ or $K_i$. 
Fitting it with a model that assumes a step change in $R$ misattributes the rates. 
For ranking, use late timepoints only.
