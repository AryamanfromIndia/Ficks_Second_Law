# Fick's Second Law: Non-Steady-State Diffusion in Solids

An interactive Streamlit web app that solves Fick's second law for a semi-infinite solid whose surface is held at a constant concentration. Enter a temperature, the concentrations and a point (x, Cx) inside the material. The app finds the time at which that concentration is reached and plots the full concentration profile at that time.

All concentrations are in **weight percent (wt %)**.

## Theory

For a semi-infinite solid with constant surface concentration Cs and uniform initial concentration C₀, the solution of Fick's second law is

```
(Cx − C₀) / (Cs − C₀) = 1 − erf( x / (2√(Dt)) )
```

where the Gaussian error function is

```
erf(z) = (2/√π) ∫₀ᶻ e^(−y²) dy
```

and the diffusion coefficient follows an Arrhenius relation:

```
D = D₀ · exp( −Q_d / RT )
```

| Symbol | Meaning | Unit |
|---|---|---|
| D₀ | Pre-exponential factor | m²/s |
| Q_d | Activation energy for diffusion | J/mol |
| R | Gas constant (8.314) | J/(mol·K) |
| T | Absolute temperature | K |
| C₀ | Initial concentration in the solid | wt % |
| Cs | Constant surface concentration | wt % |
| Cx | Concentration at depth x | wt % |

**Solving for time.** Rearranging gives

```
erf(z) = (Cs − Cx) / (Cs − C₀),   z = x / (2√(Dt))
z = erf⁻¹[ (Cs − Cx) / (Cs − C₀) ]
t = x² / (4 D z²)
```

The app computes z with the inverse error function (`scipy.special.erfinv`), then t.

## Features

- Built-in data table of D₀ and Q_d in three categories:
  - Interstitial diffusion (C and N in α-Fe and γ-Fe)
  - Self-diffusion (Fe, Cu, Al, Mg, Zn, Mo, Ni)
  - Interdiffusion, vacancy mechanism (Zn/Cu, Cu/Zn, Cu/Al, Mg/Al, Cu/Ni, Ni/Cu)
- Option to enter your own D₀ and Q_d.
- Temperature in °C or K.
- Depth in nm, µm, mm, cm or m.
- Time given in minutes; very large or small values use powers of 10.
- Step-by-step working shown in the app.
- Plot of Cx vs x with Cs, C₀ and your (x, Cx) point marked.
- Download of the concentration profile as CSV.
- Input checks: Cx must lie between C₀ and Cs, x must be positive, and Cs must differ from C₀.
- Warning if α-Fe is used above 912 °C or γ-Fe below 912 °C.

## Installation

Requires Python 3.9 or newer.

```bash
pip install streamlit numpy scipy pandas matplotlib
```

## Running the app

```bash
streamlit run fick_app.py
```

If `streamlit` is not recognised, use:

```bash
python3 -m streamlit run fick_app.py
```

The app opens in your browser at http://localhost:8501. Stop it with `Ctrl + C` in the terminal.

> **Tip:** avoid apostrophes or special characters in the project folder name (for example `Fick's_LawProject`), since they can break the terminal command.

## How to use

1. In the sidebar, choose the diffusion type and the material pair, or tick the box to enter D₀ and Q_d manually.
2. Enter the temperature.
3. Enter the initial concentration C₀ and the surface concentration Cs (wt %).
4. Enter the depth x and the concentration Cx you want at that depth.
5. Read the time, D and z from the result cards, and see the profile plot below.

### Worked example

Carburizing of steel: carbon in γ-Fe (FCC) at 927 °C, with C₀ = 0.25 wt %, Cs = 1.20 wt %, and Cx = 0.80 wt % at x = 0.50 mm.

- D ≈ 8.3 × 10⁻¹² m²/s
- t ≈ 815 min (about 13.6 h)

## Project structure

```
fick_app.py   # the Streamlit app (data, calculation, plot)
README.md     # this file
```

## Assumptions and limitations

- Semi-infinite solid, with diffusion in one direction only.
- D is constant, independent of concentration and position.
- Surface concentration Cs is constant from t = 0.
- Initial concentration C₀ is uniform.
- Diffusion data are taken from a textbook table, and real values vary with composition and purity.
- Iron has two sets of coefficients because of its phase change at 912 °C.

## Dependencies

| Package | Use |
|---|---|
| streamlit | Web interface |
| numpy | Numerical calculations |
| scipy | `erf` and `erfinv` |
| pandas | Data tables and CSV export |
| matplotlib | Concentration profile plot |
