# Map 3: the film. Does each layer of the sheet hold still while the edge moves?
# Mike One, 30 October 2026. Rules written before running.
#
# Two readings of a layer through time:
#   STILL: each layer keeps its place on S; its couplings do not change with time. Prediction: every drift = 0.
#   EDGE:  the layers ride with the universe's edge, S_layer(t) = S_layer,0 + ln(t/t0), like p = 1 for the vacuum.
# What can be seen: only dimensionless numbers. A shift of ALL layers together, floor included, is a change of
# units and invisible. What is visible is a layer moving relative to another layer:
#   (a) particle layers relative to the floor (Planck scale): G m_p^2 / (hbar c) = (m_p / M_P)^2 changes as
#       exp(-2 dS) -> EDGE predicts (dG/dt)/G = -2/t0 today (sign: G m_p^2 falls as the layers climb).
#   (b) the strong layer relative to the electron layer: mu = m_p / m_e.
#   (c) the charged medium between the Z and the electron: alpha at low energy.
#   (d) our vacuum layer: p (Section 7 of paper 5).
# PASS for STILL: each measured drift within 2 sigma of 0. EDGE is excluded where its prediction is more than
# 5 sigma from the measurement.
import math

GYR = 3.156e16; t0 = 13.8   # Gyr
rows = []
def add(layer, what, meas, err, where, source, edge_pred=None):
    z = abs(meas) / err
    rows.append((layer, what, meas, err, where, source, z, edge_pred))

# (a) gravity vs particle layers: lunar laser ranging, today
add("particles vs floor (gravity)", "dG/dt / G  [1/yr]", 7.1e-14, 7.6e-14, "today (LLR, 48 yr)",
    "Hofmann & Mueller 2018, CQG 35, 035015", edge_pred=-2 / (t0 * 1e9))
# (b) strong layer vs electron layer
add("strong vs electron (mu = m_p/m_e)", "Delta mu / mu", -1.0e-7, math.hypot(0.8e-7, 1.0e-7), "z = 0.89, 7.5 Gyr ago",
    "Bagdonaite et al. 2013, PRL 111, 231101 (methanol, PKS 1830-211)")
# (c) charged medium: alpha
add("charged medium (alpha)", "d alpha/dt / alpha  [1/yr]", 1.8e-19, 2.5e-19, "today (Yb+ optical clocks)",
    "Filzinger et al. 2023, arXiv:2301.03433")
add("charged medium (alpha)", "Delta alpha / alpha", 1.3e-6, math.hypot(1.3e-6, 0.4e-6), "z = 1.15, ~8 Gyr ago",
    "Murphy et al. 2022, A&A 658, A123 (ESPRESSO, HE 0515-4414)")
add("charged medium (alpha)", "Delta alpha / alpha", 3.6e-3, 3.7e-3, "z ~ 1100, 13.8 Gyr ago (CMB)",
    "Planck Collaboration 2015, A&A 580, A22")
add("electron layer (m_e)", "Delta m_e / m_e", 4e-3, 11e-3, "z ~ 1100 (CMB)", "Planck Collaboration 2015, A&A 580, A22")
# (d) our vacuum layer
add("our layer (vacuum)", "p  (0 = still, 1 = rides the edge)", 0.011, 0.005, "z < ~2 (SN + BAO + CMB)",
    "paper 5, Section 7")

print(f"{'layer':36s} {'quantity':36s} {'measured':>22s} {'from 0':>7s}  STILL")
allpass = True
for layer, what, m, e, where, src, z, edge in rows:
    ok = z < 2
    allpass &= ok
    print(f"{layer:36s} {what:36s} {m:+10.2e} +- {e:8.2e} {z:6.1f}s  {'pass' if ok else 'not zero'}   [{where}; {src}]")
print()
llr = rows[0]
zedge = abs(llr[7] - llr[2]) / llr[3]
print(f"EDGE reading, particle layers riding the edge: predicts dG/dt / G = {llr[7]:+.2e} /yr; "
      f"LLR measures {llr[2]:+.1e} +- {llr[3]:.1e} -> {zedge:.0f} sigma away: {'EXCLUDED' if zedge > 5 else 'not excluded'}")
# relative motion of the strong and electron layers allowed by mu: m_p ~ Lambda_QCD ~ exp(-S_conf); m_e ~ exp(-S_e)
dS_lim = 2 * rows[1][3]
print(f"strong layer vs electron layer: |Delta S| < {dS_lim:.1e} over 7.5 Gyr (2 sigma, from Delta mu / mu); "
      f"riding the edge would move them by ln(13.8/6.3) = {math.log(13.8/6.3):.2f} if they moved differently")
print(f"our vacuum layer: p = 0.011 +- 0.005 is 2.2 sigma from 0 and 198 sigma from 1 (a slow drift, not a ride)")
