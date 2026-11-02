# The galactic acceleration scale as the Jeans rate of the vacuum

Mike One · Independent researcher, Italy

> g† = c·√(4πGρ_Λ)/2π = c²·√(Λ/2)/2π = 1.056 × 10⁻¹⁰ m s⁻²

A short letter proposing that the acceleration scale of the galactic radial acceleration relation is the Jeans rate of the vacuum energy density. A ground-state postulate plus scale invariance gives the observed form of the relation with no free parameter. Every test below was run with its pass rule fixed in advance, and no parameter was fitted to the data tested.

## Contents

- `The_galactic_acceleration_scale_as_the_Jeans_rate_of_the_vacuum_Mike_One.pdf`: the letter
- `letter_source.zip`: its LaTeX source (A&A Letters format, with `aa.cls`)
- `research_note.pdf`: a longer note with the full derivation, every test and the search history
- Scripts:
  - `jeans_test.py`: galaxy–galaxy weak-lensing circular velocities (Mistele et al. 2024)
  - `jeans_test_late.py`: KiDS-1000 lensing RAR, late-type and blue lenses (Brouwer et al. 2021)
  - `sparc_vacuum_gdag.py`: SPARC test of g† = c²√(Λ/3)/2π, the commonly quoted form, against cH₀/2π
  - `sparc_likelihood_scan.py`: SPARC ΔBIC as a function of g†; output in `dbic_curve.json`

- Paper 5, *Hyperconical gravity applied to the renormalization scale*, and its scripts now have their own repository: https://github.com/THEREALMIKEONE/Hyperconical-gravity-renormalization-scale-Repository

## Results

| Dataset | Jeans rate (this work) | cH(z)/2π | c²√(Λ/3)/2π |
|---|---|---|---|
| SPARC, 153 galaxies (ΔBIC) | −7.8 | −7.8 | +53 |
| Galaxy–galaxy lensing (Δχ²) | +3.2 | +0.9 | +19.0 |
| KiDS late/blue lensing (Δχ²) | +1.2 | +0.2 | +13.8 |
| High-z kinematics, R = +0.33 ± 0.20 | 1.7σ | excluded, 3.4σ | 1.7σ |

Pass rules: ΔBIC < 10, Δχ² < 4, within 2σ.

## Rerunning

```
pip install numpy scipy sympy

# lensing data (public compilation by R. Monjo)
git clone https://github.com/robertmonjo/lensing_RAR_KiDS
python jeans_test.py      lensing_RAR_KiDS/data/external_reference/mistele_gg_wl.txt
python jeans_test_late.py lensing_RAR_KiDS/data
```

The SPARC scripts reuse the unchanged likelihood of the companion paper's test suite (`GM_physical_plane_paper.zip`, at https://github.com/THEREALMIKEONE/A-scale-dependent-vacuum-Repository). Put the SPARC rotation curves in `<run dir>/data/sparc/` and `SPARC_Lelli2016c.mrt` in `<run dir>/data/`, then:

```
python sparc_vacuum_gdag.py   <GM_physical_plane_paper dir> <run dir>
python sparc_likelihood_scan.py <GM_physical_plane_paper dir> <run dir> dbic_curve.json
```

The high-redshift result is the redshift test (cell p08) of the same suite.

## Predictions

1. The scale is constant for late-type lenses at z ≳ 1. The expansion form predicts 1.8 times larger at z = 1.
2. It follows the vacuum density: with Ω_m h² fixed by the CMB, H₀ = 73 instead of 67.4 km/s/Mpc would raise it by about 12%, against 8% for cH₀/2π.
3. If dark energy evolves, g† ∝ √ρ_DE(z).

## Use of AI tools

AI-assisted tools were used to write the analysis code and to draft and edit the text. The author directed the work, checked all code, results and text, and takes full responsibility for the content.

## Acknowledgments

The lensing data are from the public compilation of R. Monjo (lensing_RAR_KiDS).
