# Dinámica lateral con balanceo: límites de vuelco

Este documento explica el modelo de `modelo_balanceo.py` y su aplicación interactiva `app_balanceo.py` (Streamlit). Amplía el modelo de la bicicleta no lineal visto en clase con dos efectos:

- el **balanceo** de la carrocería;
- la **transferencia lateral de carga** en cada eje.

Con ellos se puede estudiar qué le ocurre antes a un vehículo en una curva: que deslice o que vuelque.

---

## 1. Convenios y notación

| Símbolo | Significado |
|---|---|
| $V_x$, $V_y$ | Velocidades longitudinal y lateral del c.d.g. (ejes del vehículo) |
| $\beta = V_y/V_x$ | Ángulo de deriva del vehículo |
| $\omega$ | Velocidad de guiñada |
| $\varphi$ | Ángulo de balanceo de la masa suspendida (positivo hacia el exterior en un giro a izquierdas) |
| $\delta$ | Ángulo de giro de la rueda delantera |
| $a$, $b$, $L = a+b$ | Distancias del c.d.g. a los ejes delantero y trasero, y batalla |
| $h$ | Altura del c.d.g. |
| $h_r$ | Altura del eje de balanceo |
| $h' = h - h_r$ | Brazo de balanceo |
| $t_f$, $t_r$ | Vías delantera y trasera |
| $C_{\alpha f}$, $C_{\alpha r}$ | Rigidez de deriva de cada eje (las dos ruedas juntas) |
| $\mu_f$, $\mu_r$ | Coeficientes de adherencia |
| $K_\varphi = K_{\varphi f} + K_{\varphi r}$ | Rigidez total al balanceo (muelles y barras estabilizadoras) |
| $C_\varphi$ | Amortiguamiento del balanceo |
| $W_f = m g\,b/L$, $W_r = m g\,a/L$ | Cargas estáticas en los ejes |

El convenio de signos es el de las clases:

- la fuerza lateral se escribe $F_y = C_\alpha\,\alpha$;
- en un giro a izquierdas las derivas son positivas.

Las derivas de los ejes son:

$$\alpha_f = \delta - \beta - \frac{a\,\omega}{V_x}, \qquad \alpha_r = \frac{b\,\omega}{V_x} - \beta$$

En un giro a izquierdas, la rueda **izquierda** es la **interior** (se descarga) y la **derecha** es la **exterior** (se carga).

---

## 2. Hipótesis del modelo

1. El vehículo es una bicicleta con tres grados de libertad: $\beta$, $\omega$ y $\varphi$. La velocidad $V_x$ es constante o una función impuesta del tiempo.
2. Toda la masa se considera **suspendida**: no se separan masas no suspendidas.
3. La carrocería gira alrededor de un **eje de balanceo horizontal** a la altura $h_r$. Esa altura se obtiene de los centros de balanceo de las dos suspensiones.
4. Las dos ruedas de un eje tienen la misma deriva, pero cada una tiene su propia carga vertical. Por eso la fuerza del eje depende de la transferencia de carga.
5. La carrocería es rígida a torsión. Si una rueda interior se levanta, ese eje ya no puede transmitir más momento, y el que le falta lo absorbe el otro eje.
6. No se modela el **vuelco por tropiezo** (bordillo, terreno blando, llanta que se clava). Solo se considera el vuelco en llano, provocado por la aceleración lateral.

---

## 3. Ecuaciones del movimiento

**Fuerza lateral.** La masa suspendida, al balancear, desplaza lateralmente su c.d.g. Por eso aparece el término $-h'\ddot\varphi$:

$$m\,(a_y - h'\ddot\varphi) = F_{yf} + F_{yr}, \qquad a_y = V_x\,(\dot\beta + \omega)$$

**Guiñada:**

$$I_z\,\dot\omega = a\,F_{yf} - b\,F_{yr}$$

**Balanceo**, con los momentos tomados respecto al eje de balanceo:

$$I_{xr}\,\ddot\varphi = m\,h'\,a_y + m\,g\,h'\,\varphi - K_\varphi\,\varphi - C_\varphi\,\dot\varphi, \qquad I_{xr} = I_x + m\,h'^2$$

Los términos del lado derecho son:

- $m h' a_y$: momento de la fuerza de inercia lateral;
- $m g h' \varphi$: momento del peso desplazado, que es **desestabilizador**;
- $K_\varphi\varphi$ y $C_\varphi\dot\varphi$: reacción de la suspensión.

La rigidez efectiva al balanceo es, por tanto, $K_\varphi - m g h'$. Si esta cantidad es negativa, la carrocería es inestable por sí sola.

**Amortiguamiento.** Se define mediante el amortiguamiento relativo $\zeta_\varphi$:

$$C_\varphi = 2\,\zeta_\varphi\sqrt{(K_\varphi - m g h')\,I_{xr}}$$

**Resolución.** Las ecuaciones de fuerza lateral y de balanceo están acopladas a través de $\dot\beta$ y $\ddot\varphi$. En cada paso de integración se resuelven como un sistema lineal de $2\times2$.

---

## 4. Transferencia lateral de carga

En cada eje $i$, la carga se reparte entre la rueda exterior y la interior así:

$$F_{z,i}^{\,ext} = \frac{F_{z,i}}{2} + \Delta F_{z,i}, \qquad F_{z,i}^{\,int} = \frac{F_{z,i}}{2} - \Delta F_{z,i}$$

La transferencia tiene dos partes:

$$\Delta F_{z,i} = \underbrace{\frac{K_{\varphi i}\,\varphi + C_{\varphi i}\,\dot\varphi}{t_i}}_{\text{elástica (a través de la suspensión)}} + \underbrace{\frac{F_{yi}\,h_r}{t_i}}_{\text{geométrica (centro de balanceo)}}$$

- **Parte elástica.** Pasa por muelles, barras y amortiguadores. Su reparto entre ejes lo fija el **reparto de $K_\varphi$**: $K_{\varphi f} = \text{rep}_f\,K_\varphi$.
- **Parte geométrica.** Pasa directamente por los brazos de la suspensión. Depende de la altura de los centros de balanceo.

**Índice de transferencia (LTR, *load transfer ratio*)** de cada eje:

$$\text{LTR}_i = \frac{F_{z,i}^{\,dcha} - F_{z,i}^{\,izq}}{F_{z,i}^{\,dcha} + F_{z,i}^{\,izq}}$$

Vale 0 sin transferencia y $\pm1$ cuando la rueda interior de ese eje se levanta.

**Saturación.** Cuando $\Delta F_{z,i}$ alcanza $F_{z,i}/2$, la rueda interior de ese eje se ha levantado. Ese eje ya no puede aportar más momento de vuelco, y el exceso pasa al otro eje.

**Vuelco.** El momento total que pueden resistir los dos ejes es:

$$M_{max} = t_f\,\frac{F_{zf}}{2} + t_r\,\frac{F_{zr}}{2}$$

Cuando el momento de vuelco lo alcanza, las **dos ruedas interiores** están en el aire y empieza el vuelco. Las simulaciones se detienen en ese instante, porque a partir de ahí el modelo deja de ser válido.

---

## 5. Neumáticos no lineales con sensibilidad a la carga

En cada rueda se usa la fórmula mágica simplificada:

$$F_y = \mu\,F_z\,\sin\!\big(C\,\arctan(B\,\alpha)\big), \qquad B = \frac{C_\alpha}{C\,\mu\,F_z}$$

Con la **sensibilidad a la carga** activada, y con $r = F_z/F_{z0}$ (la carga dividida por la carga estática de esa rueda):

$$C_\alpha = C_{\alpha 0}\,r\,(1{,}5 - 0{,}5\,r), \qquad \mu = \mu_0\,[1 - 0{,}1\,(r - 1)]$$

**Por qué importa:**

- La fuerza que gana la rueda exterior es **menor** que la que pierde la interior. Por eso, cuanta más transferencia tiene un eje, menos fuerza lateral total puede dar.
- Sin sensibilidad a la carga, $C_\alpha \propto F_z$ y $\mu$ es constante. En ese caso la fuerza del eje no depende de cómo se reparta la carga, y la transferencia no influye en la adherencia. Solo influye en el vuelco.

**Consecuencia práctica (barras estabilizadoras).** Si se aumenta la rigidez al balanceo delantera (más reparto de $K_\varphi$ hacia delante), el eje delantero soporta más transferencia y pierde antes adherencia. El vehículo se vuelve **más subvirador**. Si se aumenta la trasera, se vuelve **más sobrevirador**.

---

## 6. Límites de vuelco cuasi-estáticos

### 6.1 Factor de estabilidad estática (SSF)

Para un sólido rígido, sin suspensión, se toman momentos respecto a la rueda exterior. El vuelco empieza cuando $m\,a_y\,h = m\,g\,t/2$:

$$\text{SSF} = \frac{a_{y,v}}{g} = \frac{t}{2h}$$

Con vías distintas se usa la vía equivalente $t = (b\,t_f + a\,t_r)/L$.

### 6.2 Gradiente de balanceo

En régimen estacionario ($\ddot\varphi = \dot\varphi = 0$):

$$K_\varphi\,\varphi = m\,h'\,(a_y + g\,\varphi) \;\;\Rightarrow\;\; R_\varphi = \frac{\varphi}{a_y/g} = \frac{m\,g\,h'}{K_\varphi - m\,g\,h'} \quad \text{(rad/}g\text{)}$$

### 6.3 Umbral de vuelco con suspensión

Al balancear, el c.d.g. se desplaza hacia el exterior una distancia $h'\varphi$, lo que reduce el brazo estabilizador. El equilibrio de momentos respecto a las ruedas exteriores da:

$$m\,a_y\,h = m\,g\,\Big(\frac{t}{2} - h'\varphi\Big)$$

Sustituyendo $\varphi = R_\varphi\, a_y/g$:

$$\boxed{\;\frac{a_{y,v}}{g} = \frac{t}{2h}\;\frac{1}{1 + R_\varphi\,\dfrac{h'}{h}}\;}$$

El balanceo siempre **reduce** el umbral de vuelco respecto al SSF. La reducción es mayor cuanto más blanda es la suspensión al balanceo y más bajo está el eje de balanceo.

### 6.4 Levantamiento de la primera rueda

La rueda interior del eje $i$ se levanta cuando $\Delta F_{z,i} = F_{z,i}/2$. En estacionario, con $F_{yf} = m\,a_y\,b/L$, $F_{yr} = m\,a_y\,a/L$ y $\varphi = R\,a_y$ ($R$ en rad por m/s²):

$$a_{y,f} = \frac{(F_{zf}/2)\,t_f}{K_{\varphi f}\,R + m\,\dfrac{b}{L}\,h_r}, \qquad a_{y,r} = \frac{(F_{zr}/2)\,t_r}{K_{\varphi r}\,R + m\,\dfrac{a}{L}\,h_r}$$

Normalmente se levanta primero la rueda del eje con mayor rigidez al balanceo. En los turismos suele ser el delantero.

### 6.5 ¿Desliza o vuelca?

Hay que comparar el umbral de vuelco con el límite de adherencia, que es aproximadamente $\mu$ y algo menos si hay sensibilidad a la carga:

- **umbral > adherencia**: el vehículo desliza antes de volcar. Es el caso de los turismos;
- **umbral < adherencia**: el vehículo puede volcar en llano sin tropezar con nada. Es el caso de los vehículos altos y cargados.

---

## 7. Ensayos simulados

| Ensayo | Referencia | Qué se ve |
|---|---|---|
| **Radio constante** | ISO 4138 | Un conductor sigue un círculo mientras la velocidad sube lentamente. Se obtiene la aceleración lateral máxima y si el vehículo desliza o vuelca. |
| **Doble cambio de carril** | Inspirado en ISO 3888-1 | La anchura de los carriles se calcula a partir de la anchura del vehículo. Se comprueba si toca conos y cuánta transferencia aparece en la maniobra. |
| **Escalón de volante (J-turn)** | NHTSA / ISO 7401 | Giro brusco del volante (1000 °/s) a velocidad constante. Muestra la respuesta transitoria y la sobreoscilación del balanceo. |
| **Fishhook** | NHTSA | Volante a 720 °/s hasta una amplitud A, espera de 0,25 s, contravolante hasta −A durante 3 s y vuelta a cero. A = 6,5 veces el ángulo que da 0,3 g en estacionario. Es el ensayo de referencia para el **vuelco dinámico**. |

**Conductores:**

- **Radio constante:** seguimiento de trayectoria (*pure pursuit*) más el término de subviraje $L/R + K_{us}\,a_y/g$, con una pequeña acción integral sobre el error de radio.
- **Doble cambio de carril:** anticipa la curvatura de la trayectoria de referencia 0,25 s y corrige el error lateral y el de orientación. Además, la dirección tiene un retraso de primer orden de 0,1 s.
- **J-turn y fishhook:** ángulo de volante impuesto, en lazo abierto. Se usa una relación de dirección de 16:1.

**Final de un ensayo.** Cada ensayo termina por una de estas causas:

- **ok:** termina el tiempo previsto;
- **vuelco:** se levantan las dos ruedas interiores;
- **trompo:** $|\beta| > 60°$.

### Vuelco dinámico por debajo del umbral estático

En el fishhook, el vehículo puede volcar con una aceleración lateral **inferior** al umbral cuasi-estático, por dos motivos:

1. El balanceo es un sistema de segundo orden poco amortiguado. Ante un cambio brusco, $\varphi$ sobreoscila por encima de su valor estacionario.
2. El contravolante hace que la carrocería pase de un lado a otro. Mientras lo hace, el término $C_\varphi\dot\varphi$ añade transferencia.

---

## 8. Vehículos predefinidos

| Vehículo | SSF | Umbral de vuelco | $R_\varphi$ | Límite en giro estacionario | Fishhook |
|---|---|---|---|---|---|
| Turismo | 1,41 | 1,34 g | 3,5 °/g | 0,85 g, desliza (subviraje) | No vuelca; a 65 km/h, $\lvert\text{LTR}\rvert_{max} \approx 0{,}71$ |
| Todoterreno (SUV) | 1,07 | 1,01 g | 5,0 °/g | 0,82 g, desliza (subviraje) | No vuelca; a 65 km/h, $\lvert\text{LTR}\rvert_{max} \approx 0{,}93$ |
| Furgoneta cargada | 0,81 | 0,76 g | 6,0 °/g | 0,76 g, **vuelca** | Vuelca desde 40 km/h |

Qué muestra cada vehículo:

- **Turismo:** el margen hasta el vuelco es muy amplio. Desliza mucho antes de levantar una rueda.
- **SUV:** el umbral de vuelco está cerca de la adherencia. En maniobras bruscas la transferencia se acerca a 1.
- **Furgoneta cargada:** el umbral de vuelco está por debajo de la adherencia, así que vuelca antes de deslizar. En maniobras dinámicas vuelca incluso con aceleraciones laterales menores que su umbral.

**Sin sensibilidad a la carga:** el límite de adherencia del turismo y del SUV sube a unos 0,90 g, porque la transferencia ya no reduce la fuerza de los ejes. El vuelco de la furgoneta no cambia, ya que depende solo de la geometría y de la suspensión.

---

## 9. Uso de la aplicación

### Barra lateral

- **Vehículo predefinido:** carga los parámetros de uno de los tres vehículos, que después se pueden modificar.
- **Masa y geometría:** $m$, $I_z$, $I_x$, $a$, $b$, $h$, $h_r$, $t_f$, $t_r$.
- **Neumáticos:** $C_{\alpha f}$, $C_{\alpha r}$, $\mu_f$, $\mu_r$, factor de forma $C$ y la casilla de sensibilidad a la carga.
- **Suspensión:** $K_\varphi$, reparto delantero de $K_\varphi$ y amortiguamiento relativo $\zeta_\varphi$.

### Indicadores superiores

- SSF;
- umbral de vuelco con balanceo;
- gradiente de balanceo;
- primera rueda que se levanta;
- límite en giro estacionario;
- qué ocurre primero: desliza, derrapa o vuelca.

### Pestañas

1. **Vuelco cuasi-estático**
   - Cargas en las cuatro ruedas frente a $a_y/g$, con el límite de adherencia, el umbral de vuelco y el SSF marcados.
   - Vista trasera del vehículo en el límite.
   - Curva del umbral de vuelco frente a la altura del c.d.g., con la zona en la que el vehículo vuelca antes de deslizar.
2. **Giro estacionario**
   - Diagrama de maniobrabilidad ($\alpha_f - \alpha_r$ frente a $a_y/g$) del modelo con balanceo y del modelo plano.
   - Ángulo de balanceo frente a $a_y/g$.
3. **Ensayos**
   - Los cuatro ensayos, con comparación opcional con el modelo plano.
   - Series temporales de $\omega$, $\beta$, $a_y$, $\varphi$, LTR y $\delta$.
   - Trayectoria en planta.
   - Vista trasera en el instante que se elija.
   - En el fishhook, búsqueda de la velocidad mínima de vuelco.
4. **Ecuaciones**
   - Resumen del modelo.

### Actividades propuestas

1. Con el **turismo**, aumentar $h$ hasta que la app indique "Vuelca". Comparar el valor obtenido con el que da la fórmula del umbral de vuelco.
2. Cambiar el **reparto delantero de $K_\varphi$** del 40 % al 70 % y observar cómo cambian el diagrama de maniobrabilidad y la rueda que se levanta primero.
3. Repetir la actividad 2 con la sensibilidad a la carga **desactivada** y explicar por qué el efecto desaparece.
4. Con el **SUV**, buscar la velocidad o el factor de amplitud del fishhook con el que empieza a volcar. Probar después a subir $\zeta_\varphi$.
5. Con la **furgoneta**, comprobar en el ensayo de radio constante que vuelca a la aceleración lateral que predice el umbral cuasi-estático (0,76 g).

---

## 10. Limitaciones

- No hay masas no suspendidas, ni cabeceo, ni transferencia longitudinal de carga. La velocidad es impuesta, así que no hay frenada ni tracción.
- El eje de balanceo es horizontal y fijo. La geometría de la suspensión no cambia con el recorrido.
- Las dos ruedas de un eje tienen la misma deriva (no hay convergencia ni caída) y no hay retraso de relajación del neumático.
- La validez termina cuando se levantan las dos ruedas interiores. El vuelco completo, una vez iniciado, no se simula.
- No se incluye el vuelco por tropiezo, que es la causa de la mayoría de los vuelcos reales.

---

## 11. Archivos e instalación

| Archivo | Contenido |
|---|---|
| `modelo_balanceo.py` | Clase `Vehiculo`, cálculo de cargas, régimen estacionario, umbrales de vuelco y simulación de los ensayos |
| `app_balanceo.py` | Interfaz Streamlit |
| `requirements.txt` | `streamlit`, `numpy`, `scipy`, `plotly` |

**Ejecución en local** (en Windows usar `py -m` si `pip` no se reconoce):

```bash
python -m pip install -r requirements.txt
python -m streamlit run app_balanceo.py
```

**Streamlit Community Cloud:**

1. Subir los tres archivos al repositorio de GitHub.
2. Crear la app indicando `app_balanceo.py` como archivo principal.
