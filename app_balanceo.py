"""
Aplicación Streamlit: modelo de la bicicleta con balanceo, transferencia de carga
y límites de vuelco.

Ejecutar con:
    streamlit run app_balanceo.py

Necesita modelo_balanceo.py en la misma carpeta.
"""
import math

import numpy as np
import plotly.graph_objects as go
import streamlit as st
from plotly.subplots import make_subplots

from modelo_balanceo import (G, PRESETS, RELACION_DIR, Vehiculo, carriles_dlc,
                             comprueba_dlc, delta_03g, esquinas, simula, y_ref_dlc)

st.set_page_config(page_title="Balanceo y vuelco", layout="wide")

COL_A = "#2A2275"     # con balanceo (azul)
COL_B = "#7A1712"     # sin balanceo / plano (granate)
COL_C = "#D4910F"     # oro
GRIS = "#8A8A8A"
gr = np.degrees

# ======================================================================
#  Barra lateral
# ======================================================================
CAMPOS = ["m", "Iz", "Ix", "a", "b", "h", "hr", "tf", "tr", "Cf", "Cr", "muf", "mur", "C",
          "Kphi", "repf", "zeta"]


def aplica_preset():
    p = PRESETS[st.session_state["preset"]]
    for k in CAMPOS:
        val = p[k]
        if k in ("Cf", "Cr", "Kphi"):
            val = val / 1e3
        if k == "repf":
            val = int(round(val * 100))
        st.session_state[k] = val


if "m" not in st.session_state:
    st.session_state["preset"] = list(PRESETS)[0]
    aplica_preset()

with st.sidebar:
    st.header("Vehículo")
    st.selectbox("Vehículo predefinido", list(PRESETS), key="preset", on_change=aplica_preset)

    with st.expander("Masa y geometría", expanded=True):
        st.number_input("Masa, m (kg)", 300.0, 40000.0, step=50.0, key="m")
        st.number_input("Inercia de guiñada, I_z (kg·m²)", 200.0, 200000.0, step=100.0, key="Iz")
        st.number_input("Inercia de balanceo, I_x (kg·m²)", 50.0, 50000.0, step=50.0, key="Ix")
        st.number_input("c.d.g. – eje delantero, a (m)", 0.3, 6.0, step=0.05, key="a")
        st.number_input("c.d.g. – eje trasero, b (m)", 0.3, 6.0, step=0.05, key="b")
        st.number_input("Altura del c.d.g., h (m)", 0.2, 3.0, step=0.05, key="h")
        st.number_input("Altura del eje de balanceo, h_r (m)", 0.0, 2.0, step=0.05, key="hr")
        st.number_input("Vía delantera, t_f (m)", 0.8, 3.0, step=0.05, key="tf")
        st.number_input("Vía trasera, t_r (m)", 0.8, 3.0, step=0.05, key="tr")

    with st.expander("Neumáticos (por eje)", expanded=False):
        st.number_input("Rigidez delantera, C_αf (kN/rad)", 10.0, 2000.0, step=5.0, key="Cf")
        st.number_input("Rigidez trasera, C_αr (kN/rad)", 10.0, 2000.0, step=5.0, key="Cr")
        st.slider("Adherencia delantera, μ_f", 0.2, 1.5, step=0.05, key="muf")
        st.slider("Adherencia trasera, μ_r", 0.2, 1.5, step=0.05, key="mur")
        st.slider("Factor de forma de la fórmula mágica, C", 1.1, 1.9, step=0.05, key="C")
        st.checkbox("Sensibilidad a la carga del neumático", value=True, key="sens",
                    help="Si está activada, la transferencia de carga reduce la fuerza lateral "
                         "máxima del eje. Si no, el eje no la nota.")

    with st.expander("Suspensión", expanded=True):
        st.number_input("Rigidez total al balanceo, K_φ (kN·m/rad)", 10.0, 5000.0, step=5.0,
                        key="Kphi")
        st.slider("Reparto delantero de K_φ (%)", 20, 80, step=5, key="repf")
        st.slider("Amortiguamiento relativo del balanceo, ζ_φ", 0.05, 1.0, step=0.05, key="zeta")

ss = st.session_state
PAR = dict(m=ss.m, Iz=ss.Iz, Ix=ss.Ix, a=ss.a, b=ss.b, h=ss.h, hr=min(ss.hr, ss.h - 0.05),
           tf=ss.tf, tr=ss.tr, Cf=ss.Cf * 1e3, Cr=ss.Cr * 1e3, muf=ss.muf, mur=ss.mur, C=ss.C,
           Kphi=ss.Kphi * 1e3, repf=ss.repf / 100, zeta=ss.zeta, sens=ss.sens)
v = Vehiculo(**PAR)
CLAVE = tuple(sorted(PAR.items()))

if v.Kneta <= 0:
    st.error("La rigidez al balanceo es menor que el momento desestabilizador de la gravedad "
             "(K_φ ≤ m g h'): la carrocería vuelca por sí sola. Aumenta K_φ.")
    st.stop()


@st.cache_data(show_spinner=False)
def barrido(clave, balanceo):
    return Vehiculo(**dict(clave)).barrido_estacionario(balanceo)


@st.cache_data(show_spinner="Simulando el ensayo…")
def ensayo(clave, tipo, opciones, balanceo):
    return simula(Vehiculo(**dict(clave)), tipo, dict(opciones), balanceo)


sb = barrido(CLAVE, True)
y_lim = float(sb["y"][-1]) if len(sb["y"]) else 0.0
lf, lr = v.levantamiento_eje()

# ======================================================================
#  Cabecera
# ======================================================================
st.title("Dinámica lateral con balanceo: límites de vuelco")
st.caption("Modelo de la bicicleta con tres grados de libertad (β, ω, φ), transferencia "
           "lateral de carga en cada eje y neumáticos no lineales rueda a rueda.")

c = st.columns(6)
c[0].metric("SSF = t/(2h)", f"{v.ssf:.2f} g", help="Factor de estabilidad estática: vuelco de "
            "un sólido rígido, sin suspensión.")
c[1].metric("Umbral de vuelco con balanceo", f"{v.umbral_vuelco:.2f} g")
c[2].metric("Gradiente de balanceo", f"{gr(v.grad_balanceo):.1f} °/g")
c[3].metric("Primera rueda que se levanta",
            f"{min(lf, lr):.2f} g", "eje delantero" if lf < lr else "eje trasero",
            delta_color="off")
c[4].metric("Límite en giro estacionario", f"{y_lim:.2f} g")
if sb["causa"] == "vuelco":
    c[5].metric("¿Qué ocurre primero?", "Vuelca", "antes de deslizar", delta_color="inverse")
elif sb["causa"] == "limite_tras":
    c[5].metric("¿Qué ocurre primero?", "Derrapa", "sobreviraje límite", delta_color="off")
else:
    c[5].metric("¿Qué ocurre primero?", "Desliza", "subviraje límite", delta_color="off")


def figura(alto=420):
    fig = go.Figure()
    fig.update_layout(height=alto, margin=dict(l=10, r=10, t=40, b=10), template="plotly_white",
                      legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0))
    return fig


def vista_trasera(phi, Fz, ay, titulo="", alto=360):
    """Vista desde detrás: carrocería girada phi alrededor del centro de balanceo."""
    t = (v.tf + v.tr) / 2
    rr = 0.32
    fig = figura(alto)
    fig.add_shape(type="line", x0=-t, x1=t, y0=0, y1=0, line=dict(color="black", width=2))

    def rota(y, z):
        # punto (y, z) de la carrocería respecto al centro de balanceo; x_img = -y
        y2 = y * math.cos(phi) - (z - v.hr) * math.sin(phi)
        z2 = y * math.sin(phi) + (z - v.hr) * math.cos(phi) + v.hr
        return -y2, z2

    Fl, Fr = Fz[0] + Fz[2], Fz[1] + Fz[3]       # izquierda, derecha (suma de los dos ejes)
    for lado, F in ((1, Fl), (-1, Fr)):          # lado +1 = izquierda del vehículo
        xc = -lado * t / 2
        lev = 0.06 if F <= 1.0 else 0.0
        fig.add_shape(type="rect", x0=xc - 0.11, x1=xc + 0.11, y0=lev, y1=lev + 2 * rr,
                      fillcolor="#3A3A3A", line=dict(width=0))
        L_ = 0.45 * F / (v.m * G / 2)
        if F > 1.0:
            fig.add_annotation(x=xc, y=0, ax=xc, ay=-L_, xref="x", yref="y", axref="x", ayref="y",
                               showarrow=True, arrowhead=2, arrowwidth=2.5, arrowcolor=COL_B)
        fig.add_annotation(x=xc, y=-L_ - 0.12 if F > 1 else -0.12,
                           text=f"{F / 1e3:.1f} kN", showarrow=False, font=dict(color=COL_B))
    # carrocería
    ancho = v.ancho / 2
    z0, z1 = max(0.18, v.h - 0.40), v.h + 0.80
    pts = [rota(-ancho, z0), rota(ancho, z0), rota(ancho, z1), rota(-ancho, z1), rota(-ancho, z0)]
    fig.add_scatter(x=[p[0] for p in pts], y=[p[1] for p in pts], mode="lines", fill="toself",
                    fillcolor="rgba(42,34,117,0.15)", line=dict(color=COL_A, width=2),
                    hoverinfo="skip", showlegend=False)
    cx, cz = rota(0.0, v.h)
    fig.add_scatter(x=[cx], y=[cz], mode="markers", marker=dict(size=12, color="black"),
                    name="c.d.g.", showlegend=False)
    fig.add_scatter(x=[0], y=[v.hr], mode="markers", marker=dict(size=8, symbol="x", color=COL_C),
                    name="centro de balanceo", showlegend=False)
    if abs(ay) > 0.01:
        sx = np.sign(ay)   # a_y > 0 (giro a izquierdas): la inercia apunta a la derecha de la imagen
        fig.add_annotation(x=cx + 0.6 * sx, y=cz, ax=cx, ay=cz, xref="x", yref="y",
                           axref="x", ayref="y", showarrow=True, arrowhead=2, arrowwidth=2.5,
                           arrowcolor=COL_C)
        fig.add_annotation(x=cx + 0.6 * sx, y=cz + 0.16, text="<b>m·a<sub>y</sub></b>",
                           showarrow=False, font=dict(color="#8A5A00", size=13))
    fig.update_xaxes(range=[-t - 0.6, t + 0.6], visible=False)
    fig.update_yaxes(range=[-0.75, v.h + 1.1], visible=False, scaleanchor="x",
                     scaleratio=1)
    fig.update_layout(title=dict(text=titulo, font=dict(size=13)))
    return fig


tabs = st.tabs(["Vuelco cuasi-estático", "Giro estacionario", "Ensayos", "Ecuaciones"])

# ======================================================================
#  1. Vuelco cuasi-estático
# ======================================================================
with tabs[0]:
    col1, col2 = st.columns([3, 2])
    with col1:
        st.subheader("Cargas en las ruedas frente a la aceleración lateral")
        ys, Fzs = [], []
        for y in np.arange(0, 1.6, 0.0025):
            ay = y * G
            phi = v.m * v.hp * ay / v.Kneta
            Fz, mg = v.cargas(phi, 0.0, v.m * v.b / v.L * ay, v.m * v.a / v.L * ay)
            ys.append(y)
            Fzs.append(Fz)
            if mg <= 0:
                break
        ys, Fzs = np.array(ys), np.array(Fzs) / 1e3
        fig = figura(430)
        for k, (nm, col, dash) in enumerate((("Delantera izq. (interior)", COL_A, "solid"),
                                              ("Delantera dcha. (exterior)", COL_A, "dash"),
                                              ("Trasera izq. (interior)", COL_C, "solid"),
                                              ("Trasera dcha. (exterior)", COL_C, "dash"))):
            fig.add_scatter(x=ys, y=Fzs[:, k], mode="lines", name=nm,
                            line=dict(color=col, dash=dash, width=2.5))
        fig.add_vline(x=y_lim, line_dash="dot", line_color=COL_B,
                      annotation_text=f"límite de adherencia {y_lim:.2f} g",
                      annotation_position="top left")
        fig.add_vline(x=v.umbral_vuelco, line_dash="dot", line_color="black",
                      annotation_text=f"vuelco {v.umbral_vuelco:.2f} g",
                      annotation_position="top right")
        fig.add_vline(x=v.ssf, line_dash="dot", line_color=GRIS,
                      annotation_text=f"SSF {v.ssf:.2f}", annotation_position="bottom right")
        fig.update_xaxes(title="a_y / g", range=[0, max(v.ssf, y_lim) * 1.08])
        fig.update_yaxes(title="F_z (kN)", rangemode="tozero")
        st.plotly_chart(fig, width="stretch")
        st.caption("Giro a izquierdas. Se ha prolongado la curva más allá del límite de "
                   "adherencia (línea granate) para ver dónde estaría el vuelco si los "
                   "neumáticos aguantaran: en la realidad, el vehículo solo vuelca por esta vía "
                   "si el vuelco llega antes que el deslizamiento.")
    with col2:
        y_v = min(y_lim, v.umbral_vuelco)
        phi_v = v.m * v.hp * y_v * G / v.Kneta
        Fz_v, _ = v.cargas(phi_v, 0.0, v.m * v.b / v.L * y_v * G, v.m * v.a / v.L * y_v * G)
        st.plotly_chart(vista_trasera(phi_v, Fz_v, y_v * G,
                                      f"Vista trasera a {y_v:.2f} g (φ = {gr(phi_v):.1f}°)"),
                        width="stretch")
        st.markdown("**Umbral de vuelco**")
        st.latex(r"\frac{a_{y,v}}{g} = \frac{t}{2h}\,\frac{1}{1 + R_\varphi\,\dfrac{h - h_r}{h}}")
        st.latex(r"R_\varphi = \frac{\varphi}{a_y/g} = \frac{m g (h-h_r)}{K_\varphi - m g (h-h_r)}")
        st.markdown(f"Con estos datos: $t/(2h)$ = {v.ssf:.3f}, $R_\\varphi$ = "
                    f"{v.grad_balanceo:.4f} rad ({gr(v.grad_balanceo):.2f}°) por $g$ → "
                    f"**{v.umbral_vuelco:.3f} g**. El balanceo desplaza el c.d.g. hacia fuera "
                    f"y reduce el umbral un {(1 - v.umbral_vuelco / v.ssf) * 100:.0f} %.")

    st.subheader("¿Cuánto puede subir el centro de gravedad?")
    hs = np.linspace(0.3, 2.0, 200)
    umb = []
    for hh in hs:
        vv = Vehiculo(**{**PAR, "h": hh, "hr": min(PAR["hr"], hh - 0.05)})
        umb.append(vv.umbral_vuelco if vv.Kneta > 0 else np.nan)
    umb = np.array(umb)
    mu_min = min(v.muf, v.mur)
    fig = figura(360)
    fig.add_scatter(x=hs, y=umb, mode="lines", name="Umbral de vuelco con balanceo",
                    line=dict(color=COL_A, width=3))
    fig.add_scatter(x=hs, y=v.t_eq / (2 * hs), mode="lines", name="SSF (sólido rígido)",
                    line=dict(color=GRIS, dash="dash"))
    fig.add_hline(y=mu_min, line_dash="dot", line_color=COL_B,
                  annotation_text=f"adherencia μ = {mu_min:.2f}")
    fig.add_vline(x=v.h, line_dash="dot", line_color="black", annotation_text="h actual")
    ok = np.where(umb < mu_min)[0]
    if len(ok):
        h_crit = hs[ok[0]]
        fig.add_vrect(x0=h_crit, x1=hs[-1], fillcolor="rgba(122,23,18,0.07)", line_width=0)
        st.caption(f"Por encima de h ≈ {h_crit:.2f} m (zona sombreada) el umbral de vuelco "
                   f"queda por debajo de la adherencia: el vehículo puede volcar en llano sin "
                   f"tropezar con nada.")
    fig.update_xaxes(title="Altura del c.d.g., h (m)")
    fig.update_yaxes(title="a_y / g", range=[0, 2])
    st.plotly_chart(fig, width="stretch")

# ======================================================================
#  2. Giro estacionario
# ======================================================================
with tabs[1]:
    sp = barrido(CLAVE, False)
    causas = {"limite_del": "se satura el eje delantero (subviraje límite)",
              "limite_tras": "se satura el eje trasero (sobreviraje límite)",
              "vuelco": "se levantan las dos ruedas interiores (vuelco)"}
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("Diagrama de maniobrabilidad")
        fig = figura(430)
        fig.add_scatter(x=gr(sb["af"] - sb["ar"]), y=sb["y"], mode="lines",
                        name="Con balanceo y transferencia", line=dict(color=COL_A, width=3))
        fig.add_scatter(x=gr(sp["af"] - sp["ar"]), y=sp["y"], mode="lines",
                        name="Plano (sin transferencia)", line=dict(color=COL_B, width=2, dash="dash"))
        yl = np.linspace(0, 1.0, 2)
        fig.add_scatter(x=gr(v.Kus) * yl, y=yl, mode="lines", name="Lineal",
                        line=dict(color=GRIS, width=1.5, dash="dot"))
        fig.update_xaxes(title="α_f − α_r = δ − L/R (grados)",
                         range=[min(-2, gr(sb["af"] - sb["ar"]).min() - 0.5), 10])
        fig.update_yaxes(title="a_y / g", range=[0, max(sb["y"].max(), sp["y"].max()) * 1.1])
        st.plotly_chart(fig, width="stretch")
        st.markdown(f"Con balanceo el límite es **{y_lim:.2f} g**: {causas.get(sb['causa'], '—')}. "
                    f"Sin transferencia: {sp['y'][-1]:.2f} g.")
    with col2:
        st.subheader("Ángulo de balanceo")
        fig = figura(430)
        fig.add_scatter(x=sb["y"], y=gr(sb["phi"]), mode="lines", name="φ",
                        line=dict(color=COL_A, width=3))
        Fz = np.array(sb["Fz"])
        fig.update_xaxes(title="a_y / g")
        fig.update_yaxes(title="φ (grados)")
        st.plotly_chart(fig, width="stretch")
        if len(Fz):
            lev = [(nm, sb["y"][np.argmax(Fz[:, k] <= 1.0)]) for k, nm in
                   ((0, "delantera interior"), (2, "trasera interior")) if (Fz[:, k] <= 1.0).any()]
            if lev:
                st.markdown("Ruedas que se levantan antes del límite: " +
                            ", ".join(f"**{n}** a {y:.2f} g" for n, y in lev) + ".")
    st.info("Prueba a cambiar el **reparto delantero de K_φ** en la barra lateral: con más "
            "rigidez delante, el eje delantero soporta más transferencia, pierde adherencia y el "
            "vehículo se vuelve más subvirador. Es el efecto de las barras estabilizadoras. "
            "Requiere la sensibilidad a la carga activada.")

# ======================================================================
#  3. Ensayos
# ======================================================================
def series(res, plano=None):
    fig = make_subplots(rows=3, cols=2, shared_xaxes=True, vertical_spacing=0.07,
                        subplot_titles=("Velocidad de guiñada ω (°/s)", "Deriva del vehículo β (°)",
                                        "Aceleración lateral a_y / g", "Balanceo φ (°)",
                                        "Transferencia LTR = (F_z,dcha − F_z,izq)/F_z,eje",
                                        "Ángulo en la rueda δ (°)"))
    for r, nm, col, dash in ((res, "Con balanceo", COL_A, "solid"),
                             (plano, "Plano", COL_B, "dash")):
        if r is None:
            continue
        kw = dict(mode="lines", line=dict(color=col, dash=dash, width=2.2), legendgroup=nm)
        fig.add_scatter(x=r["t"], y=gr(r["omega"]), row=1, col=1, name=nm, **kw)
        fig.add_scatter(x=r["t"], y=gr(r["beta"]), row=1, col=2, name=nm, showlegend=False, **kw)
        fig.add_scatter(x=r["t"], y=r["ay"] / G, row=2, col=1, name=nm, showlegend=False, **kw)
        fig.add_scatter(x=r["t"], y=gr(r["delta"]), row=3, col=2, name=nm, showlegend=False, **kw)
    fig.add_scatter(x=res["t"], y=gr(res["phi"]), row=2, col=2, mode="lines", showlegend=False,
                    line=dict(color=COL_A, width=2.2))
    fig.add_scatter(x=res["t"], y=res["ltr_f"], row=3, col=1, mode="lines", name="LTR delantero",
                    line=dict(color=COL_A, width=2.2))
    fig.add_scatter(x=res["t"], y=res["ltr_r"], row=3, col=1, mode="lines", name="LTR trasero",
                    line=dict(color=COL_C, width=2.2))
    for yv in (1, -1):
        fig.add_hline(y=yv, line_dash="dot", line_color=COL_B, row=3, col=1)
    fig.update_xaxes(title_text="t (s)", row=3, col=1)
    fig.update_xaxes(title_text="t (s)", row=3, col=2)
    fig.update_layout(height=820, margin=dict(l=10, r=10, t=60, b=10), template="plotly_white",
                      legend=dict(orientation="h", yanchor="bottom", y=1.04, x=0))
    return fig


def planta(res, plano=None, n=12, dlc=False, R=None, alto=380):
    fig = figura(alto)
    if dlc:
        for x0, x1, yc, anc in carriles_dlc(v.ancho):
            xs = np.arange(x0, x1 + 0.01, 2.5)
            for lado in (-1, 1):
                fig.add_scatter(x=xs, y=np.full_like(xs, yc + lado * anc / 2), mode="markers",
                                showlegend=False, hoverinfo="skip",
                                marker=dict(symbol="triangle-up", size=8, color="#E8740C"))
        xr = np.linspace(-20, 140, 300)
        fig.add_scatter(x=xr, y=y_ref_dlc(xr), mode="lines", name="Referencia",
                        line=dict(color=GRIS, width=1, dash="dot"))
    if R is not None:
        th = np.linspace(0, 2 * np.pi, 300)
        fig.add_scatter(x=R * np.sin(th), y=R - R * np.cos(th), mode="lines", name="Círculo",
                        line=dict(color=COL_C, width=1.5, dash="dot"))
    idx = np.linspace(0, len(res["t"]) - 1, n).astype(int)
    ex, ey = esquinas(res["X"][idx], res["Y"][idx], res["psi"][idx], v.a, v.b, v.ancho)
    for k in range(len(idx)):
        o = [0, 1, 3, 2, 0]
        fig.add_scatter(x=ex[o, k], y=ey[o, k], mode="lines", fill="toself", showlegend=False,
                        fillcolor="rgba(42,34,117,0.10)", line=dict(color=COL_A, width=1),
                        hoverinfo="skip")
    fig.add_scatter(x=res["X"], y=res["Y"], mode="lines", name="Con balanceo",
                    line=dict(color=COL_A, width=2.5))
    if plano is not None:
        fig.add_scatter(x=plano["X"], y=plano["Y"], mode="lines", name="Plano",
                        line=dict(color=COL_B, width=1.5, dash="dash"))
    fig.update_xaxes(title="X (m)")
    if dlc:
        fig.update_xaxes(range=[-20, 140])
        fig.update_yaxes(title="Y (m)", range=[-4, 8])
    else:
        fig.update_yaxes(title="Y (m)", scaleanchor="x", scaleratio=1)
    return fig


def resumen(res):
    m = st.columns(5)
    fin = {"ok": "Sin vuelco", "vuelco": "Vuelco", "trompo": "Trompo"}[res["fin"]]
    m[0].metric("Resultado", fin, f"en t = {res['t'][-1]:.2f} s" if res["fin"] != "ok" else None,
                delta_color="inverse")
    m[1].metric("|LTR| máx. delantero", f"{np.abs(res['ltr_f']).max():.2f}")
    m[2].metric("|LTR| máx. trasero", f"{np.abs(res['ltr_r']).max():.2f}")
    m[3].metric("|φ| máx.", f"{gr(np.abs(res['phi']).max()):.1f}°")
    m[4].metric("|a_y| máx.", f"{np.abs(res['ay']).max() / G:.2f} g")
    st.caption("LTR (load transfer ratio): 0 sin transferencia, ±1 rueda interior levantada. "
               "El ensayo se detiene si se levantan las dos ruedas interiores (inicio del vuelco) "
               "o si |β| > 60° (trompo).")


def instante(res):
    t = res["t"]
    k0 = int(np.argmax(np.abs(res["ltr_f"]) + np.abs(res["ltr_r"])))
    tt = st.slider("Instante para la vista trasera (s)", float(t[0]), float(t[-1]), float(t[k0]),
                   step=0.01, key=f"inst_{len(t)}_{t[-1]:.3f}")
    k = int(np.argmin(np.abs(t - tt)))
    return st.plotly_chart(vista_trasera(res["phi"][k], res["Fz"][k], res["ay"][k],
                                         f"t = {t[k]:.2f} s · φ = {gr(res['phi'][k]):.1f}° · "
                                         f"a_y = {res['ay'][k] / G:.2f} g"), width="stretch")


with tabs[2]:
    tipo = st.radio("Ensayo", ["Radio constante", "Doble cambio de carril",
                               "Escalón de volante (J-turn)", "Fishhook (NHTSA)"], horizontal=True)
    comparar = st.checkbox("Comparar con el modelo plano (sin balanceo ni transferencia)",
                           value=True)

    if tipo == "Radio constante":
        st.caption("Tipo ISO 4138: un conductor sigue un círculo mientras la velocidad sube "
                   "lentamente; se ve hasta dónde llega la aceleración lateral y si el vehículo "
                   "desliza o vuelca.")
        cc = st.columns(4)
        R = cc[0].slider("Radio, R (m)", 20, 300, 60, step=10)
        V0 = cc[1].slider("Velocidad inicial (km/h)", 10, 80, 30, step=5)
        V1 = cc[2].slider("Velocidad final (km/h)", 40, 200, 110, step=5)
        Tr = cc[3].slider("Duración de la rampa (s)", 20, 120, 60, step=10)
        op = (("R", float(R)), ("V0", float(V0)), ("V1", float(max(V1, V0 + 5))), ("Tr", float(Tr)))
        res = ensayo(CLAVE, "radio", op, True)
        pl = ensayo(CLAVE, "radio", op, False) if comparar else None
        resumen(res)
        col1, col2 = st.columns([3, 2])
        with col1:
            st.plotly_chart(planta(res, pl, R=float(R), alto=430), width="stretch")
        with col2:
            instante(res)
        with st.expander("Series temporales"):
            st.plotly_chart(series(res, pl), width="stretch")

    elif tipo == "Doble cambio de carril":
        st.caption("Geometría inspirada en ISO 3888-1, con la anchura de carril calculada a "
                   "partir de la del vehículo. Conductor con anticipación de la curvatura.")
        Vd = st.slider("Velocidad (km/h)", 30, 120, 70, step=5)
        op = (("V", float(Vd)), ("desp", 3.5))
        res = ensayo(CLAVE, "dlc", op, True)
        pl = ensayo(CLAVE, "dlc", op, False) if comparar else None
        resumen(res)
        marg = comprueba_dlc(res, v)
        st.markdown("Margen mínimo a los conos: " + " · ".join(
            f"{n}: **{'—' if np.isnan(m_) else f'{m_ * 100:+.0f} cm'}**"
            for n, m_ in zip(("entrada", "carril desplazado", "salida"), marg)))
        st.plotly_chart(planta(res, pl, n=14, dlc=True), width="stretch")
        col1, col2 = st.columns([3, 2])
        with col1:
            st.plotly_chart(series(res, pl), width="stretch")
        with col2:
            instante(res)

    elif tipo == "Escalón de volante (J-turn)":
        st.caption("Giro brusco del volante a velocidad constante y mantenido (J-turn de la "
                   "NHTSA: 1000 °/s en el volante).")
        cc = st.columns(3)
        Ve = cc[0].slider("Velocidad (km/h)", 30, 150, 80, step=5)
        amp = cc[1].slider("δ en la rueda (grados)", 0.5, 20.0, 4.0, step=0.5)
        T = cc[2].slider("Tiempo simulado (s)", 2.0, 10.0, 5.0, step=0.5)
        tr = math.radians(amp) / (math.radians(1000.0) / RELACION_DIR)
        st.caption(f"Volante: {amp * RELACION_DIR:.0f}° (relación {RELACION_DIR:.0f}:1). "
                   f"Tiempo de subida: {tr:.2f} s.")
        op = (("V", float(Ve)), ("delta", float(amp)), ("tr", tr), ("T", float(T)))
        res = ensayo(CLAVE, "escalon", op, True)
        pl = ensayo(CLAVE, "escalon", op, False) if comparar else None
        resumen(res)
        col1, col2 = st.columns([3, 2])
        with col1:
            st.plotly_chart(series(res, pl), width="stretch")
        with col2:
            instante(res)
            st.plotly_chart(planta(res, pl, alto=360), width="stretch")

    else:
        st.caption("Fishhook de la NHTSA: giro a 720 °/s en el volante hasta A, espera de "
                   "0,25 s, contravolante hasta −A, se mantiene 3 s y se vuelve a cero. "
                   "A = 6,5 veces el ángulo que da 0,3 g en giro estacionario. Es la maniobra "
                   "con la que se ensaya el vuelco dinámico.")
        cc = st.columns(3)
        Vf = cc[0].slider("Velocidad de entrada (km/h)", 40, 120, 65, step=5)
        fac = cc[1].slider("Factor sobre δ(0,3 g)", 2.0, 10.0, 6.5, step=0.5)
        A = fac * delta_03g(v, Vf / 3.6)
        cc[2].metric("A en la rueda / en el volante",
                     f"{gr(A):.1f}° / {gr(A) * RELACION_DIR:.0f}°")
        op = (("V", float(Vf)), ("delta", float(gr(A))))
        res = ensayo(CLAVE, "fishhook", op, True)
        pl = ensayo(CLAVE, "fishhook", op, False) if comparar else None
        resumen(res)
        if st.checkbox("Buscar la velocidad mínima de vuelco (de 5 en 5 km/h)"):
            vmin = None
            with st.spinner("Barriendo velocidades…"):
                for Vk in range(40, 125, 5):
                    Ak = fac * delta_03g(v, Vk / 3.6)
                    r = ensayo(CLAVE, "fishhook", (("V", float(Vk)), ("delta", float(gr(Ak)))), True)
                    if r["fin"] == "vuelco":
                        vmin = Vk
                        break
            if vmin is None:
                st.success("No vuelca en el fishhook entre 40 y 120 km/h.")
            else:
                st.error(f"Vuelca a partir de **{vmin} km/h**.")
        col1, col2 = st.columns([3, 2])
        with col1:
            st.plotly_chart(series(res, pl), width="stretch")
        with col2:
            instante(res)
            st.plotly_chart(planta(res, pl, alto=360), width="stretch")

# ======================================================================
#  4. Ecuaciones
# ======================================================================
with tabs[3]:
    st.markdown("**Grados de libertad:** deriva β = V_y/V_x, velocidad de guiñada ω y balanceo "
                "φ de la masa suspendida alrededor del eje de balanceo (altura $h_r$). "
                "$h' = h - h_r$, $I_{xr} = I_x + m h'^2$.")
    st.latex(r"m\,(a_y - h'\ddot\varphi) = F_{yf} + F_{yr} \qquad a_y = V_x(\dot\beta + \omega)")
    st.latex(r"I_z\,\dot\omega = a\,F_{yf} - b\,F_{yr}")
    st.latex(r"I_{xr}\,\ddot\varphi = m h' a_y + m g h'\varphi - K_\varphi\,\varphi - C_\varphi\,\dot\varphi")
    st.markdown("**Transferencia de carga** en cada eje $i$ (parte elástica + parte del centro "
                "de balanceo):")
    st.latex(r"\Delta F_{z,i} = \frac{K_{\varphi i}\,\varphi + C_{\varphi i}\,\dot\varphi + F_{yi}\,h_r}{t_i}"
             r"\qquad F_{z,i}^{\,ext/int} = \frac{F_{z,i}}{2} \pm \Delta F_{z,i}")
    st.markdown("**Neumáticos:** fórmula mágica en cada rueda, con la misma deriva en las dos "
                "ruedas del eje. Con sensibilidad a la carga ($r = F_z/F_{z0}$):")
    st.latex(r"C_\alpha = C_{\alpha 0}\,r\,(1{,}5 - 0{,}5\,r) \qquad \mu = \mu_0\,[1 - 0{,}1\,(r-1)]")
    st.markdown("**Vuelco cuasi-estático**")
    st.latex(r"\frac{a_{y,v}}{g} = \frac{t}{2h}\,\frac{1}{1 + R_\varphi\,h'/h}\qquad "
             r"R_\varphi = \frac{m g h'}{K_\varphi - m g h'}")
    st.markdown("**Hipótesis y límites:** toda la masa es suspendida; eje de balanceo "
                "horizontal; carrocería rígida a torsión (cuando se levanta una rueda, el otro "
                "eje absorbe el momento que falta); sin vuelco por tropiezo (bordillo, terreno "
                "blando). El modelo deja de ser válido cuando se levantan las dos ruedas "
                "interiores: ahí empieza el vuelco.")
