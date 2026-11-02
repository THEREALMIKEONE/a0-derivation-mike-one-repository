# HMG as one slice of the scale axis: scripts

Mike One, 29 October 2026. These scripts produce every new number in the paper *Hyperconical gravity applied to the renormalization scale*. Each was run on public data with the pass rule written in the script before running. The outputs of those runs are in `outputs/`.

S is counted as in *A scale-dependent vacuum*: S = ln(R/R_floor), with R a radius and R_floor = sqrt(2) l_P.

| Number in the note | Script | Output |
| --- | --- | --- |
| HMG's slice formulas on every layer; vacuum at our layer 0.694 rho_crit; the vacuum curve equals gravity's coupling (l_P/R)^2; galactic scale 1.06e-10 m/s^2, gamma0 = 12.3, angle 0.4624 pi | `hmg_layers.py` | `outputs/hmg_layers.txt` |
| Late-type environment parameter s = 0.295 from the vacuum clock (fitted 0.329 +- 0.023) | `late_type_s.py` | `outputs/late_type_s.txt` |
| Baryonic Tully-Fisher on SPARC, observed - predicted = -0.02 +- 0.02 dex | `btfr_sparc.py` | `outputs/btfr_sparc.txt` |
| Does the vacuum follow the slice: p = 0.011 +- 0.005 (Pantheon+, DESI DR2 BAO, Planck distance priors); 145 high-z discs with family distances, p = 0.33 +- 0.20; the crossing (omega_J t0 = 0.986); g at z = 1, 2 | `family_cosmology.py` | `outputs/family_cosmology.txt` |
| p = 0.011 +- 0.005 from a fine profile in p | `p_profile.py` | `outputs/p_profile.txt` |
| Vacuum fraction in the same fit, 0.697 +- 0.004; HMG's 0.6937 at delta chi^2 = 0.94 | `vacuum_fraction_profile.py` | `outputs/vacuum_fraction_profile.txt` |
| The slice and the stack (Section 2, Table 2): gravity (l_P/R)^2 and the slice vacuum on every layer; alpha_s = 0.224 at the bottom-quark scale and 0.313 at the tau scale; 1/alpha_em = 136.7 at the electron | `stack_model.py` | `outputs/stack_model.txt`, `outputs/stack_model_rows.json` |
| alpha_s at one and two loops, independent check of the strong medium | `stacked_running.py` | `outputs/stacked_running.txt` |
| The (S, t) sheet, Figure 1: the edge S = ln(ct/R_floor) from -0.35 (Planck time) to 139.90 (today); the dark-energy layer from p = 0.011 | `sheet_figure.py` | `outputs/fig_sheet.pdf`, `outputs/sheet_figure.txt` |
| Section 6, KiDS-1000 isolated lenses: late-type 1.01 g_dag and blue 0.85 g_dag (ln = +0.01 +- 0.11, -0.16 +- 0.12); early 2.42 g_dag and red 2.71 g_dag; inner vs outer halves 2.8 and 5.0 sigma; mass corrections remove 26% and 9%; added neighbour mass beats a changed clock by delta chi^2 = 130 and 198 (amplitudes 0.47 and 0.48); MICE control | `kids_lens_tests.py` | `outputs/kids_lens_tests.txt` |

## The three maps

| Number in the paper | Script | Output |
| --- | --- | --- |
| Map 2 (Fig. 2): all four couplings from the floor to our layer; floor values 0.50, 0.019, 0.020, 0.030; gravity = EM at S = 2.0 | `map2_forces.py` | `outputs/map2_forces.txt`, `outputs/fig_map2_forces.pdf` |
| Map 3 (Table 5, Fig. 4): layer stability checks; riding the edge excluded at ~1900 sigma; strong vs electron |Delta S| < 2.6e-7 | `map3_layers_in_time.py`, `map3_film.py` | `outputs/map3_layers_in_time.txt`, `outputs/fig_map3_film.pdf` |
| Map 1 (Fig. 5): edge and active layer, standard vs HMG; 67 e-fold gap; nucleosynthesis at 34 yr vs ~2 min | `map1_edge_and_heat.py` | `outputs/map1_edge_and_heat.txt`, `outputs/fig_map1.pdf` |
| No phantom crossing: w_eff(z) of the p-family (w0 = -0.997) against the DESI DR2 CPL fit | `top_layer_w.py` | `outputs/top_layer_w.txt` |
| Fig. 1 with the standard-history edge c/H and HMG's ct | `sheet_figure.py` | `outputs/fig_sheet.pdf` |

The measurements in map 3 are quoted from their papers (sources in the script header); no data download is needed for these scripts.

## Running

    pip install numpy scipy
    python hmg_layers.py
    python stack_model.py                                     # no data needed
    python stacked_running.py                                 # no data needed
    python sheet_figure.py                                    # no data needed
    python kids_lens_tests.py <lensing_RAR_KiDS checkout>     # github.com/robertmonjo/lensing_RAR_KiDS
    python late_type_s.py <lensing_RAR_KiDS checkout>          # github.com/robertmonjo/lensing_RAR_KiDS
    python btfr_sparc.py <data dir with SPARC_Lelli2016c.mrt and sparc/>
    python family_cosmology.py <GM_physical_plane_paper> <GM_physical_plane_paper>
    python p_profile.py <GM_physical_plane_paper> <GM_physical_plane_paper>
    python vacuum_fraction_profile.py <GM_physical_plane_paper> <GM_physical_plane_paper>

`GM_physical_plane_paper` is the test suite of *A scale-dependent vacuum* (github.com/THEREALMIKEONE/A-scale-dependent-vacuum-Repository). Run its `run_all.py` once first so that the data are downloaded into its `data/` folder; the last two scripts reuse its likelihood code (Pantheon+, DESI DR2, Planck distance priors) and its high-redshift samples. The cosmology fits take 10 to 30 minutes each.

## Notes

- In `family_cosmology.py` the vacuum scales as rho_vac ~ H^(2p) with energy exchange between vacuum and matter; p = 0 is a constant vacuum, p = 1 keeps Omega_L constant at every age. Its rule W1 tests the earlier postulate Omega_vac = 2/3 and fails it (delta chi^2 = 59); the proposal uses the fitted value instead.
- The profile in p inside `family_cosmology.py` is coarse (steps of 0.025) and reports +0.004 +- 0.004; the free fit gives +0.011, and the fine profile in `p_profile.py` (steps of 0.005) gives the quoted p = 0.011 +- 0.005.
- AI-assisted tools were used to write the code. The author checked it and takes responsibility for it.
