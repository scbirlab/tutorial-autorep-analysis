## Read counts from next-generation sequencing

But we don't measure $m_i$ directly. Each sample is a random sampling of the pool, followed 
by molecular biology handling and uneven sequencing depth. Each guide also has its own capture and 
amplification efficiency, which differs between strains because the 5'UTR sequence differs.

For one guide in sample $s$, relative to its own $t = 0$ sample,

$$
\frac{c_i(t)}{c_i(0)} = \phi_s \frac{n_i(t)}{n_i(0)}
\qquad
\frac{x_i(t)}{x_i(0)} = k_s \frac{M_i(t)}{M_i(0)} = k_s \frac{m_i(t)}{m_i(0)} \frac{n_i(t)}{n_i(0)}
$$

Variables:
- $c_i(t)$: gDNA read count of strain $i$ at time $t$.
- $x_i(t)$: RNA read count of strain $i$ at time $t$.
- $\phi_s$: ratio of DNA sampling depth at time $t$ to that at time 0 for sample $s$.
- $k_s$: the same for RNA-seq.

Dividing by $t = 0$ cancels each guide's constant capture efficiency and any per-guide mRNA stability 
effect that is the same before and after induction. Dividing the RNA ratio by the DNA ratio 
cancels the number of cells, $\frac{n_i(t)}{n_i(0)}$, which contains all the population dynamics. 
Define

$$
A_i = \frac{x_i(t) / c_i(t)}{x_i(0) / c_i(0)} = \frac{k_s}{\phi_s} \frac{m_i(t)}{m_i(0)}
$$

$A_i$ is the normalised RNA/DNA ratio. It is proportional to transcripts per cell, 
relative to before induction.
