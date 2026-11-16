"""
Modelo de la bicicleta con balanceo y transferencia lateral de carga.

Grados de libertad: deriva del vehículo beta, velocidad de guiñada omega y
balanceo phi de la masa suspendida (con su velocidad p = dphi/dt).

Hipótesis
  - Velocidad de avance V_x impuesta (constante o en rampa lenta).
  - Toda la masa se considera suspendida (sin masas no suspendidas).
  - La masa suspendida gira alrededor de un eje de balanceo horizontal, a una
    altura h_r sobre el suelo; su c.d.g. está a h' = h - h_r por encima de él.
  - Las fuerzas laterales de los neumáticos actúan en los centros de balanceo.
  - Cinemática de la bicicleta: las dos ruedas de un eje tienen la misma deriva,
    pero cada una tiene su propia carga vertical y su propia fuerza lateral.
  - Neumático: fórmula mágica por rueda (E = 0) con sensibilidad a la carga
    opcional; con carga nula la rueda no da fuerza.
  - Cuando una rueda interior se descarga del todo, el momento que ya no puede
    absorber su eje pasa al otro eje (carrocería rígida a torsión). Cuando las
    dos ruedas interiores se levantan, el modelo deja de ser válido: inicio del
    vuelco.

Convenios (los mismos que en las clases)
  - Ejes del vehículo: x hacia delante, y hacia la izquierda, z hacia arriba.
  - beta = V_y / V_x;  omega: velocidad de guiñada;  derivas positivas en un giro
    a izquierdas, F_y = C_alpha * alpha.
  - phi > 0: la carrocería se inclina hacia la derecha (hacia fuera en un giro a
    izquierdas).
  - Transferencia de carga dF > 0: la rueda derecha (exterior en un giro a
    izquierdas) gana carga.

Ecuaciones (masa suspendida m, I_xr = I_x + m h'^2 respecto al eje de balanceo)
    m (a_y - h' phi'') = F_yf + F_yr,              a_y = V_x (beta' + omega)
    I_z omega'      = a F_yf - b F_yr
    I_xr phi''      = m h' a_y + m g h' phi - K_phi phi - C_phi phi'
Transferencia en cada eje i (f, r):
    dF_i = (K_phi,i phi + C_phi,i phi' + F_yi h_r) / t_i
"""
from dataclasses import dataclass
import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

G = 9.81


@dataclass(frozen=True)
class Vehiculo:
    m: float = 1500.0        # masa (kg)
    Iz: float = 2500.0       # inercia de guiñada (kg m^2)
    Ix: float = 500.0        # inercia de balanceo respecto a su c.d.g. (kg m^2)
    a: float = 1.2           # c.d.g. - eje delantero (m)
    b: float = 1.5           # c.d.g. - eje trasero (m)
    h: float = 0.55          # altura del c.d.g. (m)
    hr: float = 0.10         # altura del eje de balanceo (m)
    tf: float = 1.55         # vía delantera (m)
    tr: float = 1.55         # vía trasera (m)
    Cf: float = 80e3         # rigidez a deriva del eje delantero (N/rad)
    Cr: float = 90e3         # rigidez a deriva del eje trasero (N/rad)
    muf: float = 0.90        # adherencia del eje delantero
    mur: float = 1.00        # adherencia del eje trasero
    C: float = 1.5           # factor de forma de la fórmula mágica
    Kphi: float = 115e3      # rigidez total al balanceo (N m/rad)
    repf: float = 0.60       # fracción delantera de la rigidez al balanceo
    zeta: float = 0.30       # amortiguamiento relativo del balanceo
    sens: bool = True        # sensibilidad a la carga del neumático

    # ------------------------------------------------------------------
    #  Magnitudes derivadas
    # ------------------------------------------------------------------
    @property
    def L(self):
        return self.a + self.b

    @property
    def Fzf(self):
        return self.m * G * self.b / self.L

    @property
    def Fzr(self):
        return self.m * G * self.a / self.L

    @property
    def hp(self):
        return self.h - self.hr

    @property
    def Ixr(self):
        return self.Ix + self.m * self.hp ** 2

    @property
    def Kf(self):
        return self.repf * self.Kphi

    @property
    def Kr(self):
        return (1 - self.repf) * self.Kphi

    @property
    def Kneta(self):
        """Rigidez al balanceo neta de la gravedad (N m/rad)."""
        return self.Kphi - self.m * G * self.hp

    @property
    def Cphi(self):
        return 2 * self.zeta * math.sqrt(max(self.Kneta, 1.0) * self.Ixr)

    @property
    def Kus(self):
        """Gradiente de subviraje lineal (rad por g)."""
        return self.Fzf / self.Cf - self.Fzr / self.Cr

    @property
    def ancho(self):
        return max(self.tf, self.tr) + 0.25

    # ------------------------------------------------------------------
    #  Vuelco cuasi-estático
    # ------------------------------------------------------------------
    @property
    def grad_balanceo(self):
        """Gradiente de balanceo R_phi = phi / (a_y/g) en rad por g."""
        return self.m * G * self.hp / self.Kneta if self.Kneta > 0 else math.inf

    @property
    def t_eq(self):
        """Vía equivalente del vehículo completo."""
        return (self.b * self.tf + self.a * self.tr) / self.L

    @property
    def ssf(self):
        """Factor de estabilidad estática t/(2h): vuelco de un sólido rígido."""
        return self.t_eq / (2 * self.h)

    @property
    def umbral_vuelco(self):
        """a_y/g de vuelco con suspensión: SSF / (1 + R_phi h'/h)."""
        return self.ssf / (1 + self.grad_balanceo * self.hp / self.h)

    def levantamiento_eje(self):
        """a_y/g a la que se descarga la rueda interior de cada eje (cuasi-estático)."""
        R = self.grad_balanceo / G                     # rad por m/s^2
        ayf = (self.Fzf / 2) * self.tf / (self.Kf * R + self.m * self.b / self.L * self.hr)
        ayr = (self.Fzr / 2) * self.tr / (self.Kr * R + self.m * self.a / self.L * self.hr)
        return ayf / G, ayr / G

    # ------------------------------------------------------------------
    #  Neumáticos
    # ------------------------------------------------------------------
    def fy_rueda(self, alpha, Fz, Fz0, Ca0, mu0):
        """Fuerza lateral de una rueda (fórmula mágica, E = 0)."""
        if Fz <= 0.0:
            return 0.0
        r = Fz / Fz0
        if self.sens:
            Ca = Ca0 * r * (1.5 - 0.5 * r)               # rigidez a deriva cóncava con la carga
            mu = mu0 * (1.0 - 0.10 * (r - 1.0))          # la adherencia baja con la carga
        else:
            Ca, mu = Ca0 * r, mu0                        # proporcional: el eje no nota la transferencia
        if Ca <= 0.0:
            return 0.0
        B = Ca / (self.C * mu * Fz)
        return mu * Fz * math.sin(self.C * math.atan(B * alpha))

    def fy_eje(self, eje, alpha, Fz_izq, Fz_dch):
        if eje == "f":
            Fz0, Ca0, mu0 = self.Fzf / 2, self.Cf / 2, self.muf
        else:
            Fz0, Ca0, mu0 = self.Fzr / 2, self.Cr / 2, self.mur
        return (self.fy_rueda(alpha, Fz_izq, Fz0, Ca0, mu0)
                + self.fy_rueda(alpha, Fz_dch, Fz0, Ca0, mu0))

    # ------------------------------------------------------------------
    #  Transferencia de carga
    # ------------------------------------------------------------------
    def cargas(self, phi, p, Fyf, Fyr):
        """Cargas de las cuatro ruedas (izq. y dcha. de cada eje) y margen de vuelco.

        margen = 1 - |momento de vuelco| / momento máximo que pueden absorber las
        ruedas; margen <= 0: se levantan las dos ruedas interiores."""
        cf = self.Cphi * self.repf
        cr = self.Cphi * (1 - self.repf)
        dFf = (self.Kf * phi + cf * p + Fyf * self.hr) / self.tf
        dFr = (self.Kr * phi + cr * p + Fyr * self.hr) / self.tr
        limf, limr = self.Fzf / 2, self.Fzr / 2
        M = self.tf * dFf + self.tr * dFr
        Mmax = self.tf * limf + self.tr * limr
        margen = 1.0 - abs(M) / Mmax
        if abs(dFf) > limf:
            dFf = math.copysign(limf, dFf)
            dFr = (M - self.tf * dFf) / self.tr
        if abs(dFr) > limr:
            dFr = math.copysign(limr, dFr)
            dFf = (M - self.tr * dFr) / self.tf
            dFf = max(-limf, min(limf, dFf))
        return (limf - dFf, limf + dFf, limr - dFr, limr + dFr), margen

    # ------------------------------------------------------------------
    #  Dinámica
    # ------------------------------------------------------------------
    def evalua(self, x, delta, V, balanceo=True):
        """Derivadas y magnitudes auxiliares para el estado x = (beta, omega, phi, p)."""
        beta, omega, phi, p = x
        af = delta - beta - self.a * omega / V
        ar = self.b * omega / V - beta
        if balanceo:
            # punto fijo: las cargas dependen de las fuerzas en los centros de balanceo
            Fz, margen = self.cargas(phi, p, 0.0, 0.0)
            for _ in range(3):
                Fyf = self.fy_eje("f", af, Fz[0], Fz[1])
                Fyr = self.fy_eje("r", ar, Fz[2], Fz[3])
                Fz, margen = self.cargas(phi, p, Fyf, Fyr)
            Fyf = self.fy_eje("f", af, Fz[0], Fz[1])
            Fyr = self.fy_eje("r", ar, Fz[2], Fz[3])
            m, hp, Ixr = self.m, self.hp, self.Ixr
            # [ m V      -m h' ] [beta' ]   [ Fyf + Fyr - m V omega                       ]
            # [ -m h' V   I_xr ] [phi'' ] = [ m h' V omega + m g h' phi - K phi - C p     ]
            r1 = Fyf + Fyr - m * V * omega
            r2 = m * hp * V * omega + (m * G * hp - self.Kphi) * phi - self.Cphi * p
            det = m * V * Ixr - (m * hp) ** 2 * V
            db = (r1 * Ixr + m * hp * r2) / det
            dp = (m * V * r2 + m * hp * V * r1) / det
        else:
            Fz = (self.Fzf / 2, self.Fzf / 2, self.Fzr / 2, self.Fzr / 2)
            margen = 1.0
            Fyf = self.fy_eje("f", af, Fz[0], Fz[1])
            Fyr = self.fy_eje("r", ar, Fz[2], Fz[3])
            db = (Fyf + Fyr) / (self.m * V) - omega
            dp = 0.0
        dw = (self.a * Fyf - self.b * Fyr) / self.Iz
        ay = V * (db + omega)
        return (db, dw, p if balanceo else 0.0, dp), dict(af=af, ar=ar, Fyf=Fyf, Fyr=Fyr,
                                                         Fz=Fz, margen=margen, ay=ay)

    # ------------------------------------------------------------------
    #  Giro estacionario (cuasi-estático)
    # ------------------------------------------------------------------
    def _deriva_para(self, eje, F, Fz_izq, Fz_dch):
        """Deriva necesaria para que el eje dé la fuerza F (o None si no puede)."""
        al = np.radians(np.linspace(0, 30, 301))
        Fy = np.array([self.fy_eje(eje, x, Fz_izq, Fz_dch) for x in al])
        i = int(np.argmax(Fy))
        if F > Fy[i]:
            return None
        if F <= 0:
            return 0.0
        return brentq(lambda x: self.fy_eje(eje, x, Fz_izq, Fz_dch) - F, 0.0, al[i])

    def estacionario(self, ay, balanceo=True):
        """Estado de giro estacionario para una aceleración lateral ay (m/s^2).

        Devuelve dict con alpha_f, alpha_r, phi, Fz y 'estado':
        'ok', 'limite_del', 'limite_tras' o 'vuelco'."""
        Fyf = self.m * self.b / self.L * ay
        Fyr = self.m * self.a / self.L * ay
        if balanceo:
            phi = self.m * self.hp * ay / self.Kneta if self.Kneta > 0 else math.inf
            Fz, margen = self.cargas(phi, 0.0, Fyf, Fyr)
            if margen <= 0:
                return dict(estado="vuelco", phi=phi, Fz=Fz)
        else:
            phi = 0.0
            Fz = (self.Fzf / 2, self.Fzf / 2, self.Fzr / 2, self.Fzr / 2)
        af = self._deriva_para("f", Fyf, Fz[0], Fz[1])
        ar = self._deriva_para("r", Fyr, Fz[2], Fz[3])
        if af is None:
            return dict(estado="limite_del", phi=phi, Fz=Fz)
        if ar is None:
            return dict(estado="limite_tras", phi=phi, Fz=Fz)
        return dict(estado="ok", af=af, ar=ar, phi=phi, Fz=Fz)

    def barrido_estacionario(self, balanceo=True, paso=0.005, ymax=1.6):
        """Recorre a_y/g hasta el límite. Devuelve arrays y la causa del límite."""
        ys, afs, ars, phis, Fzs = [], [], [], [], []
        causa = None
        for y in np.arange(0.0, ymax, paso):
            e = self.estacionario(y * G, balanceo)
            if e["estado"] != "ok":
                causa = e["estado"]
                break
            ys.append(y); afs.append(e["af"]); ars.append(e["ar"])
            phis.append(e["phi"]); Fzs.append(e["Fz"])
        return dict(y=np.array(ys), af=np.array(afs), ar=np.array(ars), phi=np.array(phis),
                    Fz=np.array(Fzs), causa=causa)


# ======================================================================
#  Geometría del doble cambio de carril (inspirada en ISO 3888-1)
# ======================================================================
VOLADIZO = 0.9
DELTA_MAX = math.radians(25.0)
TAU_DIR = 0.1
RELACION_DIR = 16.0          # relación de desmultiplicación de la dirección


def carriles_dlc(w, desp=3.5):
    return [(0.0, 15.0, 0.0, 1.1 * w + 0.25),
            (45.0, 70.0, desp, 1.2 * w + 0.25),
            (95.0, 125.0, 0.0, 1.3 * w + 0.25)]


def y_ref_dlc(x, desp=3.5):
    x = np.asarray(x, dtype=float)
    s1 = np.clip((x - 15.0) / 30.0, 0, 1)
    s2 = np.clip((x - 70.0) / 25.0, 0, 1)
    return desp * ((1 - np.cos(np.pi * s1)) - (1 - np.cos(np.pi * s2))) / 2


def curv_ref_dlc(x, desp=3.5):
    x = np.asarray(x, dtype=float)
    d1 = np.zeros_like(x)
    d2 = np.zeros_like(x)
    for x0, l, sg in ((15.0, 30.0, 1), (70.0, 25.0, -1)):
        m = (x > x0) & (x < x0 + l)
        u = np.pi * (x - x0) / l
        d1 += np.where(m, sg * desp / 2 * np.sin(u) * np.pi / l, 0.0)
        d2 += np.where(m, sg * desp / 2 * np.cos(u) * (np.pi / l) ** 2, 0.0)
    return d2 / (1 + d1 ** 2) ** 1.5


def esquinas(X, Y, psi, a, b, w, vol=VOLADIZO):
    lx = np.array([a + vol, a + vol, -(b + vol), -(b + vol)])[:, None]
    ly = np.array([w / 2, -w / 2, w / 2, -w / 2])[:, None]
    c, s = np.cos(psi)[None, :], np.sin(psi)[None, :]
    return X[None, :] + lx * c - ly * s, Y[None, :] + lx * s + ly * c


def comprueba_dlc(res, v, desp=3.5):
    ex, ey = esquinas(res["X"], res["Y"], res["psi"], v.a, v.b, v.ancho)
    margenes = []
    for x0, x1, yc, anc in carriles_dlc(v.ancho, desp):
        dentro = (ex >= x0) & (ex <= x1)
        if not dentro.any():
            margenes.append(np.nan)
            continue
        margenes.append(float((anc / 2 - np.abs(ey[dentro] - yc)).min()))
    return margenes


# ======================================================================
#  Señales de volante en lazo abierto
# ======================================================================
def delta_fishhook(t, A, t0=0.5, vel=math.radians(720.0) / RELACION_DIR, espera=0.25):
    """Fishhook de la NHTSA: giro hasta A, espera, contravolante hasta -A durante 3 s
    y vuelta a cero. Velocidad de giro 720 grados/s en el volante."""
    t1 = t0 + A / vel
    t2 = t1 + espera
    t3 = t2 + 2 * A / vel
    t4 = t3 + 3.0
    t5 = t4 + A / vel
    if t < t0:
        return 0.0
    if t < t1:
        return vel * (t - t0)
    if t < t2:
        return A
    if t < t3:
        return A - vel * (t - t2)
    if t < t4:
        return -A
    if t < t5:
        return -A + vel * (t - t4)
    return 0.0


def duracion_fishhook(A, t0=0.5, vel=math.radians(720.0) / RELACION_DIR, espera=0.25):
    return t0 + 4 * A / vel + espera + 3.0 + 1.0


def delta_03g(v, V):
    """Ángulo en la rueda para 0,3 g en giro estacionario lineal a velocidad V."""
    ay = 0.3 * G
    return v.L * ay / V ** 2 + v.Kus * 0.3


# ======================================================================
#  Simulación de ensayos
# ======================================================================
def simula(v, tipo, opciones, balanceo=True):
    """
    tipo: 'radio', 'dlc', 'escalon' o 'fishhook'.
    opciones: dict (o tupla de pares) con
      radio:    R (m), V0, V1 (km/h), Tr (s)
      dlc:      V (km/h), desp (m)
      escalon:  V (km/h), delta (grados en la rueda), tr (s, tiempo de subida), T (s)
      fishhook: V (km/h), delta (grados en la rueda)
    Se detiene si |beta| > 60 grados (trompo) o si se levantan las dos ruedas
    interiores (vuelco).
    """
    o = dict(opciones)
    L, Kus = v.L, v.Kus
    t0 = 0.5
    if tipo == "radio":
        R, V0, V1, Tr = o["R"], o["V0"] / 3.6, o["V1"] / 3.6, o["Tr"]
        T = Tr
        Vfun = lambda t: V0 + (V1 - V0) * min(t / Tr, 1.0)
        ay0 = V0 ** 2 / R
        phi0 = (v.m * v.hp * ay0 / v.Kneta) if balanceo else 0.0
        x0 = [0.0, V0 / R, phi0, 0.0, 0.0, 0.0, 0.0, 0.0, L / R + Kus * ay0 / G]
    else:
        Vc = o["V"] / 3.6
        Vfun = lambda t: Vc
        x0 = [0.0] * 9
        if tipo == "dlc":
            x0[5] = -20.0
            T = (140.0 + 20.0) / Vc
        elif tipo == "escalon":
            T = o["T"]
        else:
            T = duracion_fishhook(math.radians(o["delta"]), t0)

    def delta_mando(t, s):
        beta, omega, phi, p, psi, X, Y, I, _ = s
        V = Vfun(t)
        if tipo == "escalon":
            return math.radians(o["delta"]) * min(max((t - t0) / max(o["tr"], 1e-3), 0.0), 1.0)
        if tipo == "fishhook":
            return delta_fishhook(t, math.radians(o["delta"]), t0)
        if tipo == "radio":
            d = max(V * 1.0, 4.0)
            Px, Py = X + d * math.cos(psi), Y + d * math.sin(psi)
            e = math.hypot(Px, Py - o["R"]) - o["R"]
            dl = 2 * L * e / d ** 2 + Kus * V ** 2 / (G * o["R"]) + 0.01 * I
        else:
            des = o.get("desp", 3.5)
            k = float(curv_ref_dlc(X + V * 0.25, des))
            e_y = Y - float(y_ref_dlc(X, des))
            pend = float(y_ref_dlc(X + 0.05, des) - y_ref_dlc(X - 0.05, des)) / 0.1
            e_c = psi + beta - math.atan(pend)
            dl = (L + Kus * V ** 2 / G) * k - 0.01 * e_y - 0.5 * e_c
        return max(-DELTA_MAX, min(DELTA_MAX, dl))

    conductor = tipo in ("radio", "dlc")

    def delta_real(t, s):
        return s[8] if conductor else delta_mando(t, s)

    def rhs(t, s):
        beta, omega, phi, p, psi, X, Y, I, dlt = s
        V = Vfun(t)
        (db, dw, dphi, dp), _ = v.evalua((beta, omega, phi, p), delta_real(t, s), V, balanceo)
        Vy = V * beta
        e_rad = math.hypot(X, Y - o["R"]) - o["R"] if tipo == "radio" else 0.0
        ddl = (delta_mando(t, s) - dlt) / TAU_DIR if conductor else 0.0
        return [db, dw, dphi, dp, omega, V * math.cos(psi) - Vy * math.sin(psi),
                V * math.sin(psi) + Vy * math.cos(psi), e_rad, ddl]

    def trompo(t, s):
        return math.radians(60) - abs(s[0])
    trompo.terminal = True

    def vuelco(t, s):
        _, aux = v.evalua(s[:4], delta_real(t, s), Vfun(t), balanceo)
        return aux["margen"] + 1e-9 if balanceo else 1.0
    vuelco.terminal = True

    sol = solve_ivp(rhs, (0, T), x0, max_step=0.01, rtol=1e-6, atol=1e-8,
                    events=[trompo, vuelco])
    t = sol.t
    beta, omega, phi, p, psi, X, Y, I, _ = sol.y
    Vt = np.array([Vfun(tt) for tt in t])
    dl = np.array([delta_real(tt, s) for tt, s in zip(t, sol.y.T)])
    aux = [v.evalua(s[:4], d, V, balanceo)[1] for s, d, V in zip(sol.y.T, dl, Vt)]
    Fz = np.array([a["Fz"] for a in aux])
    ltr_f = (Fz[:, 1] - Fz[:, 0]) / (Fz[:, 1] + Fz[:, 0])
    ltr_r = (Fz[:, 3] - Fz[:, 2]) / (Fz[:, 3] + Fz[:, 2])
    fin = "ok"
    if sol.status == 1:
        fin = "trompo" if len(sol.t_events[0]) else "vuelco"
    return dict(t=t, beta=beta, omega=omega, phi=phi, p=p, psi=psi, X=X, Y=Y, delta=dl,
                V=Vt, ay=np.array([a["ay"] for a in aux]), af=np.array([a["af"] for a in aux]),
                ar=np.array([a["ar"] for a in aux]), Fz=Fz, ltr_f=ltr_f, ltr_r=ltr_r,
                margen=np.array([a["margen"] for a in aux]), fin=fin)


# ======================================================================
#  Vehículos predefinidos
# ======================================================================
PRESETS = {
    "Turismo": dict(m=1500.0, Iz=2500.0, Ix=500.0, a=1.20, b=1.50, h=0.55, hr=0.10,
                    tf=1.55, tr=1.55, Cf=80e3, Cr=90e3, muf=0.90, mur=1.00, C=1.5,
                    Kphi=115e3, repf=0.60, zeta=0.30),
    "Todoterreno (SUV)": dict(m=2100.0, Iz=4000.0, Ix=900.0, a=1.35, b=1.50, h=0.75, hr=0.25,
                              tf=1.60, tr=1.60, Cf=110e3, Cr=120e3, muf=0.90, mur=1.00, C=1.5,
                              Kphi=128e3, repf=0.60, zeta=0.30),
    "Furgoneta cargada": dict(m=3500.0, Iz=7000.0, Ix=2000.0, a=1.60, b=1.90, h=1.05, hr=0.40,
                              tf=1.70, tr=1.70, Cf=150e3, Cr=170e3, muf=0.90, mur=1.00, C=1.5,
                              Kphi=235e3, repf=0.55, zeta=0.30),
}
