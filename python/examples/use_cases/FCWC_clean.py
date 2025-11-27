import numpy as np
import pandas as pd
from scipy.stats import beta
from plotnine import (
    ggplot, aes, geom_line, geom_point, facet_wrap,
    scale_color_manual, labs, theme_minimal, theme, element_text
)
import matplotlib
matplotlib.use("TkAgg")
# matplotlib.use("Qt5Agg")
from matplotlib import pyplot as plt


# Reproducibility
np.random.seed(2025)

n_sims = 100000

# Beta helper: from mean m and pseudo-count k
def beta_from_mean_k(m, k):
    shape1 = m * k
    shape2 = (1 - m) * k
    return shape1, shape2


### Fixed components from authors ####

# Confession prevalence P(C | I) ~ 0.44
pC_par = beta_from_mean_k(0.44, 100)

# Overall wrongful conviction rate P(W) ~ 0.03
pW_overall_par_base = beta_from_mean_k(0.03, 100)

# Share of wrongful convictions involving false confession P(F | W) ~ 0.15
pF_givenW_par = beta_from_mean_k(0.15, 100)

# Confession-specific wrongful conviction P(W | C,I) ~ 0.06 (higher than overall rate)
pW_givenC_fixed = 0.06
pW_givenC_par_fixed = beta_from_mean_k(pW_givenC_fixed, 100)


###############################################
# PANEL A — vary P(W | C,I) from 0.04 to 0.10 #
###############################################

pW_givenC_means = np.arange(0.04, 0.101, 0.01)

def run_panelA(mean_pW_givenC):

    # Shared components
    pC = beta.rvs(*pC_par, size=n_sims)
    pW_overall = beta.rvs(*pW_overall_par_base, size=n_sims)
    pF_givenW = beta.rvs(*pF_givenW_par, size=n_sims)

    # Scenario-specific P(W | C, I)
    pW_givenC_par = beta_from_mean_k(mean_pW_givenC, 100)
    pW_givenC = beta.rvs(*pW_givenC_par, size=n_sims)

    # Authors' factorisation
    fcwc_auth = pC * pW_overall * pF_givenW

    # Corrected factorisation: depends on P(W | C,I)
    fcwc_corr = pC * pW_givenC * pF_givenW

    return pd.DataFrame({
        "panel": ["Panel A: Vary P(W | C,I)"],
        "x_param": [mean_pW_givenC],
        "authors_mean": [fcwc_auth.mean()],
        "corrected_mean": [fcwc_corr.mean()]
    })


panelA = pd.concat([run_panelA(m) for m in pW_givenC_means],
                   ignore_index=True)


#########################################
# PANEL B — vary P(W) from 0.01 to 0.05 #
#########################################

pW_overall_means = np.arange(0.01, 0.051, 0.01)

def run_panelB(mean_pW_overall):

    # Shared components
    pC = beta.rvs(*pC_par, size=n_sims)
    pF_givenW = beta.rvs(*pF_givenW_par, size=n_sims)

    # Scenario-specific P(W)
    pW_overall_par = beta_from_mean_k(mean_pW_overall, 100)
    pW_overall = beta.rvs(*pW_overall_par, size=n_sims)

    # Fixed P(W | C,I)
    pW_givenC = beta.rvs(*pW_givenC_par_fixed, size=n_sims)

    # Authors' factorisation (depends on P(W))
    fcwc_auth = pC * pW_overall * pF_givenW

    # Corrected factorisation (depends on P(W | C,I))
    fcwc_corr = pC * pW_givenC * pF_givenW

    return pd.DataFrame({
        "panel": ["Panel B: Vary P(W)"],
        "x_param": [mean_pW_overall],
        "authors_mean": [fcwc_auth.mean()],
        "corrected_mean": [fcwc_corr.mean()]
    })


panelB = pd.concat([run_panelB(m) for m in pW_overall_means],
                   ignore_index=True)


#####################
# Combine & prepare #
#####################

plot_data = pd.concat([panelA, panelB], ignore_index=True)

plot_data["x_pct"] = 100 * plot_data["x_param"]
plot_data["authors_mean_pct"] = 100 * plot_data["authors_mean"]
plot_data["corrected_mean_pct"] = 100 * plot_data["corrected_mean"]

plot_data["panel"] = pd.Categorical(
    plot_data["panel"],
    categories=[
        "Panel A: Vary P(W | C,I)",
        "Panel B: Vary P(W)"
    ],
    ordered=True
)


###########
# PLOT    #
###########

plot = (
    ggplot(plot_data, aes("x_pct"))
    + geom_line(aes(y="authors_mean_pct",
                    color='"Authors’ factorisation"'),
                size=1.2)
    + geom_point(aes(y="authors_mean_pct",
                     color='"Authors’ factorisation"'),
                 size=2)
    + geom_line(aes(y="corrected_mean_pct",
                    color='"Corrected factorisation"'),
                size=1.2)
    + geom_point(aes(y="corrected_mean_pct",
                     color='"Corrected factorisation"'),
                 size=2)
    + facet_wrap("~panel", scales="free_x")
    + scale_color_manual(
        values={
            "Authors’ factorisation": "#D55E00",
            "Corrected factorisation": "#0072B2"
        },
        name="Model"
    )
    + labs(
        title="FCWC Base-Rate Estimates Under Two Factorisations",
        subtitle="Using authors’ approximate mean inputs: "
                 "P(C) ≈ 0.44, P(W) ≈ 0.03, P(F | W) ≈ 0.15",
        x="Assumed mean (%)",
        y="Mean FCWC base rate (%)",
        caption=(
            "Simulated using 100,000 Monte Carlo draws per scenario.\n"
            "In Panel A, only the corrected factorisation depends on P(W | C,I).\n"
            "In Panel B, only the authors’ factorisation depends on P(W)."
        )
    )
    + theme_minimal(base_size=14)
    + theme(
        strip_text=element_text(weight="bold", size=12),
        legend_position="bottom",
        plot_title=element_text(weight="bold", size=18),
        plot_subtitle=element_text(size=11)
    )
)

print(plot)
plot.show()
