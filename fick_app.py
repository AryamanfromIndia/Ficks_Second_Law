"""
Fick's Second Law - Non-Steady-State Diffusion in Solids
Run with:  streamlit run fick_app.py
Requires:  pip install streamlit numpy scipy pandas matplotlib
"""
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
from scipy.special import erf, erfinv

R = 8.314  # J/(mol K)

# ---------------------------------------------------------------- DATA
# (Diffusing species, Host metal, D0 [m^2/s], Qd [J/mol])
DATA = {
    "Interstitial Diffusion": [
        ("C", "Fe (α or BCC)", 1.1e-6, 87_400),
        ("C", "Fe (γ or FCC)", 2.3e-5, 148_000),
        ("N", "Fe (α or BCC)", 5.0e-7, 77_000),
        ("N", "Fe (γ or FCC)", 9.1e-5, 168_000),
    ],
    "Self-Diffusion": [
        ("Fe", "Fe (α or BCC)", 2.8e-4, 251_000),
        ("Fe", "Fe (γ or FCC)", 5.0e-5, 284_000),
        ("Cu", "Cu (FCC)", 2.5e-5, 200_000),
        ("Al", "Al (FCC)", 2.3e-4, 144_000),
        ("Mg", "Mg (HCP)", 1.5e-4, 136_000),
        ("Zn", "Zn (HCP)", 1.5e-5, 94_000),
        ("Mo", "Mo (BCC)", 1.8e-4, 461_000),
        ("Ni", "Ni (FCC)", 1.9e-4, 285_000),
    ],
    "Interdiffusion (Vacancy)": [
        ("Zn", "Cu (FCC)", 2.4e-5, 189_000),
        ("Cu", "Zn (HCP)", 2.1e-4, 124_000),
        ("Cu", "Al (FCC)", 6.5e-5, 136_000),
        ("Mg", "Al (FCC)", 1.2e-4, 130_000),
        ("Cu", "Ni (FCC)", 2.7e-5, 256_000),
        ("Ni", "Cu (FCC)", 1.9e-4, 230_000),
    ],
}


def make_df(rows):
    return pd.DataFrame(
        {
            "Diffusing Species": [r[0] for r in rows],
            "Host Metal": [r[1] for r in rows],
            "D₀ (m²/s)": [sci(r[2], 1) for r in rows],
            "Q_d (J/mol)": [f"{r[3]:,}" for r in rows],
        }
    )


_SUP = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")


def _split(x, d=2):
    """x -> (mantissa, exponent) with d decimals."""
    if x == 0:
        return 0.0, 0
    e = int(np.floor(np.log10(abs(x))))
    m = round(x / 10**e, d)
    if abs(m) >= 10:
        m /= 10
        e += 1
    return m, e


def sci(x, d=2):
    """Plain text: 2.30 × 10⁻⁵"""
    m, e = _split(x, d)
    return f"{m:.{d}f} × 10{str(e).translate(_SUP)}"


def sci_tex(x, d=2):
    """LaTeX: 2.30\\times10^{-5}"""
    m, e = _split(x, d)
    return f"{m:.{d}f}\\times10^{{{e}}}"


def fmt_time(t):
    """Seconds -> minutes string."""
    m = t / 60
    if m >= 1e6 or m < 1e-2:
        return f"{sci(m, 3)} min"
    return f"{m:,.2f}".rstrip("0").rstrip(".") + " min"


# ---------------------------------------------------------------- UI
st.set_page_config(page_title="Fick's 2nd Law", page_icon="🧪", layout="wide")
st.title("🧪 Fick's Second Law — Non-Steady-State Diffusion in Solids")
st.caption(
    "Semi-infinite solid, constant surface concentration:  "
    "(Cx − C₀)/(Cs − C₀) = 1 − erf( x / 2√(Dt) ),   D = D₀ exp(−Q_d / RT)"
)

with st.expander("📋 Table of diffusion data", expanded=False):
    for cat, rows in DATA.items():
        st.markdown(f"**{cat}**")
        st.dataframe(make_df(rows), hide_index=True, use_container_width=True)
    st.caption(
        "Iron has two sets of coefficients: BCC α-Fe is stable below 912 °C, "
        "FCC γ-Fe above 912 °C."
    )

# ---- sidebar inputs
st.sidebar.header("1 · Material system")
category = st.sidebar.selectbox("Diffusion type", list(DATA.keys()))
rows = DATA[category]
labels = [f"{r[0]} in {r[1]}" for r in rows]
use_custom = st.sidebar.checkbox("Enter D₀ and Q_d manually")
if use_custom:
    st.sidebar.write("D₀ (m²/s) = a × 10ᵇ")
    _a = st.sidebar.number_input("a (mantissa)", value=1.0, min_value=0.0001, format="%.3f")
    _b = st.sidebar.number_input("b (power of 10)", value=-5, step=1)
    D0 = _a * 10.0 ** _b
    Qd = st.sidebar.number_input("Q_d (J/mol)", value=150_000.0, step=1000.0)
    host = "custom"
else:
    choice = st.sidebar.selectbox("Diffusing species in host", labels)
    sp, host, D0, Qd = rows[labels.index(choice)]
    st.sidebar.write(f"**D₀** = {sci(D0)} m²/s")
    st.sidebar.write(f"**Q_d** = {Qd:,} J/mol")

st.sidebar.header("2 · Temperature")
unit = st.sidebar.radio("Unit", ["°C", "K"], horizontal=True)
T_in = st.sidebar.number_input(
    f"Temperature ({unit})", value=927.0 if unit == "°C" else 1200.0, step=10.0
)
T = T_in + 273.15 if unit == "°C" else T_in

st.sidebar.header("3 · Concentrations (wt %)")
C0 = st.sidebar.number_input("Initial concentration C₀", value=0.25, min_value=0.0, format="%.4f")
Cs = st.sidebar.number_input("Surface concentration Cs (constant)", value=1.20, min_value=0.0, format="%.4f")

st.sidebar.header("4 · Target point")
x_unit = st.sidebar.selectbox("Depth unit", ["mm", "µm", "m", "cm", "nm"], index=0)
x_factor = {"m": 1.0, "cm": 1e-2, "mm": 1e-3, "µm": 1e-6, "nm": 1e-9}[x_unit]
x_val = st.sidebar.number_input(f"Depth x ({x_unit})", value=0.50, min_value=0.0, format="%.4f")
Cx = st.sidebar.number_input("Concentration at x, Cx", value=0.80, min_value=0.0, format="%.4f")
x = x_val * x_factor

# ---------------------------------------------------------------- CALC
if T <= 0:
    st.error("Temperature must be above 0 K.")
    st.stop()

if "Fe" in host and not use_custom:
    if "α" in host and T - 273.15 > 912:
        st.warning("α-Fe (BCC) is not stable above 912 °C — consider the γ-Fe (FCC) data.")
    if "γ" in host and T - 273.15 < 912:
        st.warning("γ-Fe (FCC) is not stable below 912 °C — consider the α-Fe (BCC) data.")

D = D0 * np.exp(-Qd / (R * T))

errors = []
if Cs == C0:
    errors.append("Cs and C₀ must be different.")
else:
    frac = (Cs - Cx) / (Cs - C0)  # = erf(z)
    if not (0 < frac < 1):
        lo, hi = sorted([C0, Cs])
        errors.append(
            f"Cx must lie strictly between C₀ and Cs (between {lo:g} and {hi:g} wt %) "
            "for a finite, positive time."
        )
if x <= 0:
    errors.append("Depth x must be greater than 0.")
if errors:
    for e in errors:
        st.error(e)
    st.stop()

z = erfinv(frac)                 # x / (2 sqrt(Dt))
t = x**2 / (4 * D * z**2)        # seconds

# ---------------------------------------------------------------- RESULTS
c1, c2, c3, c4 = st.columns(4)
c1.metric("Diffusion coefficient D", f"{sci(D, 3)} m²/s")
c2.metric("z = x / 2√(Dt)", f"{z:.4f}")
c3.metric("Time t", fmt_time(t))
c4.metric("Diffusion length √(Dt)", f"{np.sqrt(D*t)*1e3:.4g} mm")
st.success(
    f"Cx = {Cx:g} wt % is reached at x = {x_val:g} {x_unit} after "
    f"**t = {fmt_time(t)}**  ({sci(t, 3)} s) at {T-273.15:.1f} °C ({T:.1f} K)."
)

with st.expander("Show working"):
    st.latex(r"D = D_0 \exp\!\left(-\frac{Q_d}{RT}\right)"
             rf" = {sci_tex(D0)}\exp\!\left(-\frac{{{Qd:,.0f}}}{{8.314\times{T:.1f}}}\right) = {sci_tex(D,3)}\ \mathrm{{m^2/s}}")
    st.latex(r"\mathrm{erf}(z)=\frac{2}{\sqrt{\pi}}\int_{0}^{z}e^{-y^{2}}\,dy")
    st.latex(r"\frac{C_x-C_0}{C_s-C_0}=1-\mathrm{erf}\!\left(\frac{x}{2\sqrt{Dt}}\right)"
             r"\;\Rightarrow\;\mathrm{erf}(z)=\frac{C_s-C_x}{C_s-C_0}"
             rf"={frac:.5f}\;\Rightarrow\; z={z:.5f}")
    st.latex(r"t=\frac{x^2}{4Dz^2}"
             rf"=\frac{{({sci_tex(x,3)})^2}}{{4({sci_tex(D,3)})({z:.5f})^2}}={sci_tex(t,3)}\ \mathrm{{s}}")

# ---------------------------------------------------------------- PLOT
x_max = max(4 * x, 3.5 * np.sqrt(D * t) * 2)  # show whole profile incl. target
xs = np.linspace(0, x_max, 600)
Cprof = Cs - (Cs - C0) * erf(xs / (2 * np.sqrt(D * t)))

fig, ax = plt.subplots(figsize=(9, 5))
xs_plot = xs / x_factor
ax.plot(xs_plot, Cprof, lw=2.5, color="#1f4e9c", label=f"Profile at t = {fmt_time(t)}")
ax.axhline(Cs, ls="--", color="tab:red", lw=1, label=f"Cs = {Cs:g}")
ax.axhline(C0, ls="--", color="tab:green", lw=1, label=f"C₀ = {C0:g}")
ax.plot([x_val], [Cx], "o", color="black", ms=9, zorder=5, label=f"(x, Cx) = ({x_val:g} {x_unit}, {Cx:g})")
ax.vlines(x_val, min(C0, Cs), Cx, colors="gray", ls=":")
ax.hlines(Cx, 0, x_val, colors="gray", ls=":")
ax.set_xlabel(f"Depth x ({x_unit})")
ax.set_ylabel("Concentration Cx (wt %)")
ax.set_title(f"Concentration profile — {'custom' if use_custom else choice}, {T-273.15:.0f} °C")
ax.set_xlim(0, x_max / x_factor)
ax.grid(alpha=0.3)
ax.legend()
st.pyplot(fig)

# download profile
prof_df = pd.DataFrame({f"x ({x_unit})": xs_plot, "Cx (wt %)": Cprof})
st.download_button("⬇ Download profile as CSV", prof_df.to_csv(index=False), "profile.csv", "text/csv")