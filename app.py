from functools import partial
import os

from carabiner.mpl import add_legend, grid, colorblind_palette
import gradio as gr
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.stats import multivariate_hypergeom, nbinom

# Set the default color cycle
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(
    color=("lightgrey", "dimgrey") + colorblind_palette(name="vibrant_0")[1:],
)

SEED: int = 42
MAX_TIME: float = 5.
SOURCES_DIR: str = "sources"

E_MAX: float = 2.  # transcript rate per cell, with no dCas9 bound
GROWTH_RATE: float = 1.  # intrinsic growth of non-targeting strains
PSEUDOCOUNT: float = 1.  # only to avoid log(0)

N_NT: int = 2  # non-targeting guides, pooled as the reference
STRENGTH_DEFAULTS = [1., 5., 20.]
FITNESS_DEFAULTS = [.5, 2., .2]
N_GUIDES: int = len(STRENGTH_DEFAULTS)
STRAIN_LABELS = (
    ["Non-targeting"]
    + ["_none"] * (N_NT - 1)
    + [f"Guide {i + 1}" for i in range(N_GUIDES)]
)


def inject_markdown(filename):
    with open(os.path.join(SOURCES_DIR, filename), 'r') as f:
        md = f.read()
    return gr.Markdown(
        md,
        latex_delimiters=[
            {"left": "$$", "right": "$$", "display": True},
            {"left": "$", "right": "$", "display": False},
        ],
    )


def colors():
    return mpl.rcParams['axes.prop_cycle'].by_key()['color']


def dynamics(t, y, w, K, KR, d, k_R):
    """Cells and per-cell transcripts in a pool.

    y   : cells per strain, then transcripts per cell per strain
    w   : intrinsic growth rate per strain
    K   : carrying capacity of the pool
    KR  : K_i * R at full dCas9 induction, per strain
    d   : mRNA decay rate
    k_R : rate at which dCas9 accumulates after induction

    """
    n_strains = len(w)
    n, m = y[:n_strains], y[n_strains:]
    r = w * (1. - np.sum(n) / K)  # instantaneous per-capita growth
    occupancy = KR * (1. - np.exp(-k_R * t))  # K_i R(t)
    dn = r * n
    dm = E_MAX / (1. + occupancy) - (d + r) * m
    return np.concatenate([dn, dm])


def calculate_pool(
    inoculum, inoculum_var, carrying_capacity, mrna_decay, dcas9_rate,
    *guide_params,
    n_timepoints=100,
):
    strength, fitness = guide_params[:N_GUIDES], guide_params[N_GUIDES:]
    w = GROWTH_RATE * np.asarray([1.] * N_NT + list(fitness))
    KR = np.asarray([0.] * N_NT + list(strength))
    K = inoculum * carrying_capacity

    inoculum_var = inoculum + inoculum_var * np.square(inoculum)
    p = inoculum / inoculum_var
    n = (inoculum ** 2.) / (inoculum_var - inoculum)
    n0 = nbinom.rvs(n, p, size=len(w), random_state=SEED).astype(float)
    # before induction, dCas9 is off and growth is unaffected by the guides
    m0 = E_MAX / (mrna_decay + GROWTH_RATE * (1. - np.sum(n0) / K))

    t = np.linspace(0., MAX_TIME, num=int(n_timepoints))
    ode_solution = solve_ivp(
        dynamics,
        t_span=(0., MAX_TIME),
        t_eval=t,
        y0=np.concatenate([n0, np.full(len(w), m0)]),
        args=(w, K, KR, mrna_decay, dcas9_rate),
        rtol=1e-8,
        atol=1e-8,
    )
    cells = ode_solution.y[:len(w)]
    transcripts_per_cell = ode_solution.y[len(w):]
    r = w[:,None] * (1. - np.sum(cells, axis=0, keepdims=True) / K)
    occupancy = KR[:,None] * (1. - np.exp(-dcas9_rate * t[None]))
    qss = E_MAX / (1. + occupancy) / (mrna_decay + r)
    return t, cells, transcripts_per_cell, qss


def plotter_t(x, growths, scatter=False, dashed=None, vline=None, **kwargs):
    fig, axes = grid(aspect_ratio=1.5)
    plotter_f = partial(axes.scatter, s=5.) if scatter else axes.plot
    for i, y in enumerate(growths):
        plotter_f(
            x.flatten(), y.flatten(),
            label=STRAIN_LABELS[i],
        )
    if dashed is not None:
        for i, y in enumerate(dashed):
            axes.plot(
                x.flatten(), y.flatten(),
                color=colors()[i],
                linestyle="--",
                label="_none",
            )
        axes.plot([], [], color="black", linestyle="--", label="Steady-state prediction")
    if vline is not None:
        axes.axvline(vline, color="grey", linestyle=":", label="_none")
    axes.set(
        xlabel="Time",
        yscale="log",
        **kwargs,
    )
    add_legend(axes)
    return fig


def plotter_ref(x, growths, scatter=False, fitlines=None, text=None, **kwargs):
    fig, axes = grid(aspect_ratio=1.5 if text is None else 1.7)
    plotter_f = partial(axes.scatter, s=5.) if scatter else axes.plot
    for i, y in enumerate(growths):
        plotter_f(
            x.flatten(), y.flatten(),
            label=STRAIN_LABELS[i],
        )
    if fitlines is not None:
        fit_x, intercepts = fitlines
        for i, b in enumerate(intercepts.flatten()):
            # gradient fixed at 1 on log-log axes
            axes.plot(
                fit_x.flatten(), fit_x.flatten() * np.exp(b),
                color=colors()[i],
                label="_none",
            )
    if text is not None:
        axes.text(
            1.05, 0.,
            text,
            fontsize=10,
            va="top",
            transform=axes.transAxes,
        )
    axes.set(
        xscale="log",
        yscale="log",
        **kwargs,
    )
    axes.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    axes.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    add_legend(axes)
    return fig


def plotter_estimates(true_strength, estimated_strength, text=None, **kwargs):
    fig, axes = grid(aspect_ratio=1.5 if text is None else 1.7)
    for i, (true, estimated) in enumerate(zip(true_strength, estimated_strength)):
        axes.scatter(
            [true], [estimated],
            s=40.,
            color=colors()[i + N_NT],
            label=STRAIN_LABELS[i + N_NT],
        )
    lims = [min(true_strength) / 2., max(true_strength) * 2.]
    axes.plot(lims, lims, color="lightgrey", linestyle="--", zorder=0, label="_none")
    if text is not None:
        axes.text(
            1.05, 0.,
            text,
            fontsize=10,
            va="top",
            transform=axes.transAxes,
        )
    axes.set(
        xscale="log",
        yscale="log",
        **kwargs,
    )
    add_legend(axes)
    return fig


def pool_plotter(
    inoculum, inoculum_var, carrying_capacity, mrna_decay, dcas9_rate, *guide_params,
):
    t, cells, transcripts_per_cell, qss = calculate_pool(
        inoculum, inoculum_var, carrying_capacity, mrna_decay, dcas9_rate, *guide_params,
    )
    initial = transcripts_per_cell[:,:1]
    return [
        plotter_t(
            t,
            cells,
            ylabel="Number of cells per strain",
        ),
        plotter_t(
            t,
            transcripts_per_cell / initial,
            dashed=qss / initial,
            ylabel="Transcripts per cell, relative to $t=0$",
        ),
    ]


def reads_sampler(population, sample_frac, seq_depth, reps, variance, capture=None, seed=SEED):
    samples = []
    for i, timepoint_pop in enumerate(np.split(population.astype(int), population.shape[-1], axis=-1)):
        sample_size = np.floor(timepoint_pop.sum() * sample_frac).astype(int)
        samples.append(
            multivariate_hypergeom.rvs(
                m=timepoint_pop.flatten(),
                n=sample_size,
                size=reps,
                random_state=seed + i,
            ).T
        )
    samples = np.stack(samples, axis=-2)
    if capture is not None:
        samples = samples * capture[:,None,None]
    read_means = np.floor(seq_depth * samples.shape[0] * samples / samples.sum(axis=0, keepdims=True))
    # negative binomial with variance = mean + variance * mean ** 2
    n = 1. / variance
    p = 1. / (1. + variance * read_means)
    return nbinom.rvs(n, p, random_state=seed)


def capture_efficiency(n_strains, spread, seed):
    """Per-guide capture and amplification efficiency, constant across samples."""
    return np.exp(np.random.default_rng(seed).normal(scale=spread, size=n_strains))


def calculate_reads(
    sample_frac, seq_reps, seq_depth, read_var, capture_spread,
    inoculum, inoculum_var, carrying_capacity, mrna_decay, dcas9_rate, *guide_params,
):
    t, cells, transcripts_per_cell, qss = calculate_pool(
        inoculum, inoculum_var, carrying_capacity, mrna_decay, dcas9_rate, *guide_params,
        n_timepoints=10,
    )
    n_strains = cells.shape[0]
    dna_counts = reads_sampler(
        cells,
        sample_frac, seq_depth, seq_reps, read_var,
        capture=capture_efficiency(n_strains, capture_spread, SEED + 1),
        seed=SEED + 1000,
    )
    rna_counts = reads_sampler(
        cells * transcripts_per_cell,
        sample_frac, seq_depth, seq_reps, read_var,
        capture=capture_efficiency(n_strains, capture_spread, SEED + 2),
        seed=SEED + 2000,
    )
    return t, dna_counts, rna_counts


def normalised_ratio(dna_counts, rna_counts):
    """A_i = (x_i / c_i) / (x_i(0) / c_i(0)), with t=0 counts summed over replicates."""
    x0 = rna_counts[:,0].sum(axis=-1)
    c0 = dna_counts[:,0].sum(axis=-1)
    log_r0 = np.log((x0 + PSEUDOCOUNT) / (c0 + PSEUDOCOUNT))
    log_A = (
        np.log((rna_counts + PSEUDOCOUNT) / (dna_counts + PSEUDOCOUNT))
        - log_r0[:,None,None]
    )
    return log_A, x0, c0


def autorep_fitter(t, dna_counts, rna_counts, late_time):
    log_A, x0, c0 = normalised_ratio(dna_counts, rna_counts)
    late_start = max(1, int(np.flatnonzero(t >= late_time)[0]) if np.any(t >= late_time) else len(t) - 1)
    late = slice(late_start, None)
    # pooled non-targeting reference for each late sample
    nt_ratio = rna_counts[:N_NT,late].sum(axis=0) / dna_counts[:N_NT,late].sum(axis=0)
    log_A0 = np.log(nt_ratio / (x0[:N_NT].sum() / c0[:N_NT].sum()))
    y = log_A[:,late] - log_A0[None]
    b = y.mean(axis=(1, 2))
    return log_A[:,late], log_A0, b, t[late_start]


def reads_plotter(
    sample_frac, seq_reps, seq_depth, read_var, capture_spread,
    inoculum, inoculum_var, carrying_capacity, mrna_decay, dcas9_rate, *guide_params,
):
    t, dna_counts, rna_counts = calculate_reads(
        sample_frac, seq_reps, seq_depth, read_var, capture_spread,
        inoculum, inoculum_var, carrying_capacity, mrna_decay, dcas9_rate, *guide_params,
    )
    log_A, x0, c0 = normalised_ratio(dna_counts, rna_counts)
    t_tiled = np.tile(t[:,None], (1, seq_reps))
    return pool_plotter(
        inoculum, inoculum_var, carrying_capacity, mrna_decay, dcas9_rate, *guide_params,
    ) + [
        plotter_t(
            t_tiled,
            dna_counts,
            scatter=True,
            ylabel="gDNA read counts per strain",
        ),
        plotter_t(
            t_tiled,
            rna_counts,
            scatter=True,
            ylabel="RNA read counts per strain",
        ),
        plotter_t(
            t_tiled,
            np.exp(log_A),
            scatter=True,
            ylabel="$A_i = \\frac{x_i(t)}{c_i(t)} / \\frac{x_i(0)}{c_i(0)}$",
        ),
    ]


def fit_plotter(
    late_time,
    sample_frac, seq_reps, seq_depth, read_var, capture_spread,
    inoculum, inoculum_var, carrying_capacity, mrna_decay, dcas9_rate, *guide_params,
):
    t, dna_counts, rna_counts = calculate_reads(
        sample_frac, seq_reps, seq_depth, read_var, capture_spread,
        inoculum, inoculum_var, carrying_capacity, mrna_decay, dcas9_rate, *guide_params,
    )
    log_A_late, log_A0, b, late_start_time = autorep_fitter(t, dna_counts, rna_counts, late_time)
    true_strength = list(guide_params[:N_GUIDES])
    estimated_strength = np.expm1(-b[N_NT:])
    plot_text = "\n".join(
        f"Guide {i + 1}: $\\hat{{b}}={bi:.2f}$"
        for i, bi in enumerate(b[N_NT:])
    ) + f"\nNon-targeting: $\\hat{{b}}={b[:N_NT].mean():.2f}$"
    rank_recovered = np.array_equal(
        np.argsort(estimated_strength),
        np.argsort(true_strength),
    )
    estimates_text = "\n".join(
        f"Guide {i + 1}: true {true:.1f}, estimated {estimated:.1f}"
        for i, (true, estimated) in enumerate(zip(true_strength, estimated_strength))
    ) + f"\nRank order recovered: {'yes' if rank_recovered else 'no'}"

    reference = np.exp(log_A0)
    fit_x = np.geomspace(reference.min(), reference.max(), num=20)
    return reads_plotter(
        sample_frac, seq_reps, seq_depth, read_var, capture_spread,
        inoculum, inoculum_var, carrying_capacity, mrna_decay, dcas9_rate, *guide_params,
    ) + [
        plotter_ref(
            reference,
            np.exp(log_A_late),
            scatter=True,
            fitlines=(fit_x, b),
            text=plot_text,
            xticks=np.round(np.geomspace(reference.min(), reference.max(), num=3), 2),
            xlabel="$A_0$, pooled non-targeting",
            ylabel="$A_i$",
        ),
        plotter_estimates(
            true_strength,
            estimated_strength,
            text=estimates_text,
            xlabel="True $K_iR$",
            ylabel="Estimated $K_iR$",
        ),
    ]


with gr.Blocks() as demo:
    inject_markdown("header.md")
    # Autorep and transcript dynamics
    inject_markdown("intro.md")
    inject_markdown("transcript-dynamics.md")
    with gr.Row():
        guide_strength = [
            gr.Slider(
                .1, 100.,
                step=.1,
                value=k,
                label=f"Strength $K R$, guide {i + 1}",
            )
            for i, k in enumerate(STRENGTH_DEFAULTS)
        ]
    with gr.Row():
        guide_fitness = [
            gr.Slider(
                0., 2.,
                step=.05,
                value=w,
                label=f"Relative growth rate, guide {i + 1}",
            )
            for i, w in enumerate(FITNESS_DEFAULTS)
        ]
    with gr.Row():
        mrna_decay = gr.Slider(
            .05, 100.,
            step=.01,
            value=2.,
            label="mRNA decay rate",
        )
        dcas9_rate = gr.Slider(
            .005, 1.,
            step=.005,
            value=.8,
            label="dCas9 accumulation rate",
        )
    with gr.Row():
        inoculum = gr.Slider(
            10, 1_000_000,
            step=10,
            value=1000,
            label="Average inoculum per strain",
        )
        inoculum_var = gr.Slider(
            .001, 1.,
            step=.001,
            value=.001,
            label="Inoculum variance between strains",
        )
        carrying_capacity = gr.Slider(
            N_NT + N_GUIDES + 1, 10_000,
            step=1,
            value=10,
            label="Total carrying capacity (x inoculum)",
        )
    plot_pool = gr.Button("Plot growth and transcript levels")
    pool_curves_n = gr.Plot(label="Growth vs time", format="png")
    pool_curves_m = gr.Plot(label="Transcripts per cell vs time", format="png")
    pool_curves = [pool_curves_n, pool_curves_m]
    pool_inputs = [
        inoculum, inoculum_var, carrying_capacity,
        mrna_decay, dcas9_rate,
        *guide_strength, *guide_fitness,
    ]

    # Read counts
    inject_markdown("read-counts.md")
    with gr.Row():
        sample_frac = gr.Slider(
            .001, 1., step=.001,
            value=.1,
            label="Fraction of population per sample",
        )
        seq_reps = gr.Slider(
            1, 10,
            step=1,
            value=3,
            label="Technical replicates",
        )
        seq_depth = gr.Slider(
            10, 10_000,
            step=10,
            value=10_000,
            label="Average reads per strain per sample",
        )
    with gr.Row():
        read_var = gr.Slider(
            .001, 1.,
            step=.001,
            value=.001,
            label="Sequencing variance",
        )
        capture_spread = gr.Slider(
            0., 1.,
            step=.05,
            value=.5,
            label="Per-guide capture variability (log-scale SD)",
        )
    plot_reads = gr.Button("Plot read counts")
    read_curves_dna = gr.Plot(label="gDNA read counts vs time", format="png")
    read_curves_rna = gr.Plot(label="RNA read counts vs time", format="png")
    read_curves_ratio = gr.Plot(label="Normalised RNA/DNA ratio vs time", format="png")
    read_curves = [read_curves_dna, read_curves_rna, read_curves_ratio]
    read_inputs = [sample_frac, seq_reps, seq_depth, read_var, capture_spread]

    # Fitting
    inject_markdown("nuisance-parameter.md")
    late_time = gr.Slider(
        0., MAX_TIME,
        step=1.,
        value=0.,
        label="Use samples after (minutes)",
    )
    plot_fit = gr.Button("Estimate sgRNA strength")
    fit_curves_ref = gr.Plot(label="Normalised ratio vs non-targeting reference", format="png")
    fit_curves_estimates = gr.Plot(label="Estimated vs true strength", format="png")
    fit_curves = [fit_curves_ref, fit_curves_estimates]

    # Events
    plot_pool.click(
        fn=pool_plotter,
        inputs=pool_inputs,
        outputs=pool_curves,
    )
    plot_reads.click(
        fn=reads_plotter,
        inputs=read_inputs + pool_inputs,
        outputs=pool_curves + read_curves,
    )
    plot_fit.click(
        fn=fit_plotter,
        inputs=[late_time] + read_inputs + pool_inputs,
        outputs=pool_curves + read_curves + fit_curves,
    )
if __name__ == "__main__":
    demo.queue() 
    demo.launch(share=True)