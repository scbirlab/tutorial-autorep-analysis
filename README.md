---
title: Tutorial - Autorep sgRNA strength from pooled CRISPRi-seq and RNA-seq
emoji: 🧮
colorFrom: green
colorTo: blue
sdk: gradio
sdk_version: 6.30.0
app_file: app.py
pinned: true
license: mit
short_description: Inferring sgRNA strength from NGS autorep experiments
tags:
    - biology
    - sequencing
    - pooled-screen
    - fitness
    - gradio
---

# Tutorial – Autorep sgRNA strength from pooled CRISPRi-seq and RNA-seq

[![Open in Spaces](https://huggingface.co/datasets/huggingface/badges/resolve/main/open-in-hf-spaces-md-dark.svg)](https://huggingface.co/spaces/scbirlab/tutorial-autorep-analysis)

This is the repository for the interactive tutorial, acccessed [here](https://huggingface.co/spaces/scbirlab/tutorial-autorep-analysis). 
Non-interactive text from the tutorial is below.

## The autorep experiment

CRISPRi sgRNAs differ in how strongly they repress their target, and the difference depends 
on the spacer sequence. How do we rank sgRNAs by strength, in a pool?

In the autorep design, every strain carries a fluorescent reporter (mWasabi) with the sgRNA's 
own target site in its 5'UTR. Each strain has one sgRNA, so the spacer is also the strain's barcode. 
The strains are grown together in a pool. When dCas9 is induced, each sgRNA represses its own 
reporter, and we measure two things by sequencing:

- $c_i(t)$: gDNA spacer-amplicon counts for strain $i$ at time $t$. These count **cells**.
- $x_i(t)$: RNA-seq counts of the barcoded reporter transcript for strain $i$ at time $t$. 
These count **transcripts**.

A dCas9 bound at the target site blocks transcription. If a single site is at equilibrium 
occupancy, the fraction of time the promoter is free (unrepressed) is

$$
\frac{1}{1 + K_i R}
$$

- $K_i$: affinity of the dCas9–sgRNA complex for its site. This is the sgRNA strength we want to rank.
- $R$: concentration of free dCas9. It is the same for every strain.

So production of reporter transcript from strain $i$ is

$$
E_i = \frac{E_{max}}{1 + K_i R}
$$

where $E_{max}$ is production with no dCas9 bound. Because $R$ is shared, we can only identify 
the dimensionless product $K_i R$. That is enough for ranking, since multiplying every $K_i$ 
by the same constant does not change the order.

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
