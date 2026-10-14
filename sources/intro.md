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