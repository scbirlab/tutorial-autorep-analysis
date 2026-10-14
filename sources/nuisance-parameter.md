### Removing the sample-level nuisance parameter

The factor $\frac{k_s}{\phi_s}$ is the RNA-to-DNA depth ratio of sample $s$ relative to $t = 0$. It is 
one unknown number per sample, and is shared by every guide. We get rid of it the same way as 
in the previous tutorial: take a reference.

For the non-targeting reference guides ($K_0 R = 0$, and $m_0$ constant),

$$
A_0 = \frac{k_s}{\phi_s}
$$

Before induction, dCas9 is off, so $R = 0$ and every guide is at steady state with 
$m_i(0) = \frac{E_{max}}{d + r_i}$. If the guides do not change growth before induction 
($r_i = r_0$ at $t = 0$), then $m_i(0) = m_0(0)$ and the ratio of the two normalised ratios 
is the ratio of late steady states:

$$
\frac{A_i}{A_0} = \frac{m_i}{m_0} = \frac{d + r_0}{d + r_i} \cdot \frac{1}{1 + K_i R}
$$

Plotted against the reference term across late samples, each guide gives a straight line with gradient 1 and intercept $-\log(1 + K_i R)$. The ingredients are:

- Dividing RNA by DNA removes cell number and population growth.
- Dividing by each guide's own $t = 0$ removes per-guide capture efficiency and mRNA stability.
- Dividing by the pooled non-targeting guides removes the per-sample depth factor $\frac{k_s}{\phi_s}$.
- Late timepoints are on the transcript plateau, and the growth term is 1.

What remains is $-\log(1 + K_i R)$. Averaging the left side minus the reference term over the late samples gives this value per guide, $\hat{b}_i$. The most negative $\hat{b}_i$ is the strongest guide, so $-\hat{b}_i = \log(1 + K_i R)$ ranks the guides. The transient 
after induction, the mRNA decay rate $d$, and $E_{max}$ are not needed.
