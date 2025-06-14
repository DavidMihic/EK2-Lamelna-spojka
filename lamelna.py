from math import sqrt, log, ceil, tan, pi

delta = 0.45
mi = 0.35
mi_a = 0.12
tau_tDI = 230  # N/mm^2
rps = 32
t3 = 0.75  # s
z_k = 40

T_uk = 164.72  # Nm

d_izr = round((5 * T_uk * 1000 / (tau_tDI / 10)) ** (1 / 3), 2)
print(f"d_izr: {d_izr} mm")

d_vr = 40
t = 4.9

# Provjera faktora sigurnosti
b1 = 0.90
b2 = 0.94
K1 = 1.5
beta_kt = 1.9

tau_t = round(T_uk * 1000 / (0.2 * (d_vr - t) ** 3), 2)  # N/mm^2
print(f"tau_t: {tau_t} N/mm^2")
S_post = round(b1 * b2 * tau_tDI / (K1 * beta_kt * tau_t), 2)
S_potr = 1.8
print(f"S_post: {S_post} >= S_potr = {S_potr}", S_post > S_potr)

print()
print("----------------")
print()

# Lamele

R_u_poc = 0.9 * d_vr
print(f"R_u_poc: {R_u_poc} mm")

# Unutarnja 1-100-27
D_eu = 108
D_iu = 75.2
H_u = 32.6

# Vanjska 2-300-31
D_ev = 121
D_iv = 84
H_v = 68.5
H_1v = 61

Ru = (D_iu + 2 * H_u) / 4
Rv = (H_v + H_1v) / 2
R1 = D_eu / 2
R2 = D_iv / 2

print(f"Ru: {Ru} mm")
print(f"Rv: {Rv} mm")
print(f"R1: {R1} mm")
print(f"R2: {R2} mm")

print("R2 > Ru", R2 > Ru)
print("R1 < Rv", R1 < Rv)

Rm = (R1 + R2) / 2
print(f"Rm: {Rm} mm")

print()

v_rel = round(2 * pi * rps * Rm / 1000, 2)  # m/s
v_rel_dop = 20
print(f"v_rel: {v_rel} m/s < v_rel_dop = {v_rel_dop} m/s", v_rel < v_rel_dop)

print()

f3 = round((1 + mi * mi_a * Rm / Rv) * (1 + mi * mi_a * Rm / Ru), 4)
fa = round((1 - mi * mi_a * Rm / Rv) / (1 + mi * mi_a * Rm / Rv), 4)
fi = round((1 - mi * mi_a * Rm / Ru) / (1 + mi * mi_a * Rm / Ru), 4)
print(f"f3: {f3}")
print(f"fa: {fa}")
print(f"fi: {fi}")

print()

n_max1 = round(1 + 2 * log(delta) / log(fa * fi), 2)

c = ceil(n_max1)
n_max = c - (2 if c % 2 else 1)
print(f"n_max = {n_max}")
n = 5  # odabrani broj lamela
z = max(n - 1, 1)
print(f"n:{n}")
print(f"z:{z}")

print()

print(f"suma: {sum((fa * fi) ** q for q in range(int(z // 2)))}")

F_un = round(
    T_uk * 1000 * f3 / (2 * mi * Rm * sum((fa * fi) ** q for q in range(int(z // 2)))),
    2,
)  # N

print(f"Fun: {F_un} N/m^2")

p = round(F_un / (pi * (R1**2 - R2**2)), 3)  # N/mm^2
p_dop = 0.85 * 1  # N/mm^2

print(f"p: {p} N/mm^2 < p_dop = {p_dop} N/mm^2", p <= p_dop)

print()

k = 1.4

T_P = round(mi * k * z * Rm * F_un * 1e-3, 2)
T_NS = round(
    2
    * mi
    * f3**-1
    * Rm
    * F_un
    * sum((fa * fi) ** q for q in range(int(z // 2)))
    * 1e-3,
    2,
)

T_R = round(
    T_uk * k * f3 * z / (2 * sum((fa * fi) ** q for q in range(int(z // 2)))), 2
)  # Nm
print(f"T_R: {T_R} Nm")

Q = round(T_R * pi * rps * t3, 2)  # J/ukljucivanju
print(f"Q: {Q} J/ukljucivanju")
Qz = round(Q * z_k, 2)  # J/h
print(f"Qz: {Qz} J/h")

T_ok = 293.15  # K
T_dop = 723  # K

# Jednostruko ukljucivanje
# m(unutarnja) = 55,4 g
# m(vanjska) = 127 g
m = 0.42  # kg
c = 461  # J/kgK
deltaT = round(Q / (m * c), 2)
T1 = round(deltaT + T_ok, 2)
print(f"T_1: {T1} K <= T_dop = {T_dop} K", T1 <= T_dop)

# Trajno ukljucivanje u radu
R_sr = 64.28  # mm
v_m = round(2 * pi * rps * R_sr / 1000, 2)
alfa_k = round(18800 + 25100 * v_m**1.5, 2)
print(f"alfa_k: {alfa_k}")
As = 0.05675  # m^2
T = round(Qz / (alfa_k * As) + T_ok, 2)
print(f"T: {T} K <= T_dop = {T_dop} K", T <= T_dop)

# Proklizavanje pri preopterecenju
t_pr_dop = round(m * c * (T_dop - T_ok) / (T_uk * 2 * pi * rps), 2)  # s
print(f"t_pr_dop: {t_pr_dop} s")

print()

s_0 = 0.95 * 0.5  # mm
y = 0.98
q_v = 0.086  # mm^3 / Wh
V_0 = round(s_0 * y * z * pi * (R1**2 - R2**2), 2)  # mm^3
print(f"V_0: {V_0} mm^3")
t_h = round(3600 * V_0 / (Qz * q_v), 2)  # h
print(f"t_h: {t_h} h")

print()

# Poluga (E335)
l1 = 50  # mm
l2 = 17.5  # mm
l = 40  # mm
a = 8  # mm
a0 = 4.50  # mm
b = 15  # mm
alfa = 5 * pi / 180  # rad
Re = 320  # MPa
E = 210000  # MPa
P = 2  # mm (hod matice)
fi = 15  # deg

F_uklj = round(F_un / 3, 2)  # N
print(f"F_uklj: {F_uklj} N")

F = round(F_uklj * l2 / l1, 2)  # N
print(f"F: {F} N")

sigma_s = round(6 * F * l / (b * a**2), 2)  # N/mm^2
sigma_dop = round(Re / 3, 1)  # N/mm^2
print(
    f"sigma_s: {sigma_s} N/mm^2 <= sigma_dop = {sigma_dop} N/mm^2", sigma_s <= sigma_dop
)

w = round(
    12
    * F
    / (E * b * tan(alfa) ** 2)
    * (
        -1 / tan(alfa) * log(a0)
        + a0 / (6 * tan(alfa))
        - l * (1 / (l * tan(alfa) + a0) - a0 / (2 * (l * tan(alfa) + a0) ** 2))
        + 1 / tan(alfa) * log(l * tan(alfa) + a0)
        - a0 / (6 * tan(alfa) * (l * tan(alfa) + a0))
    ),
    4,
)  # mm
print(f"w: {w} mm")

w1 = round(w * l2 / l1, 4)  # mm
print(f"w': {w1} mm")

deltaF_uklj = round(F_uklj * (P / w1) * (fi / 360), 2)  # N
print(f"deltaF_uklj: {deltaF_uklj} N")

e = 0.3  # mm (maksimalna zracnost izmedu dvije lamele)
h = e * z  # mm (ukupna maksimalna zracnost izmedu lamela)
h0 = round(h * l1 / l2 + w * (1 + deltaF_uklj / F_uklj), 4)  # mm (ukupni hod poluge)
print(f"h0: {h0} mm")

sigma_s1 = round(sigma_s * (1 + deltaF_uklj / F_uklj), 2)  # N/mm^2
print(
    f"sigma_s1: {sigma_s1} N/mm^2 <= sigma_dop = {sigma_dop} N/mm^2",
    sigma_s1 <= sigma_dop,
)

print()

# Svornjak (S355J0 -> Rm=500 MPa)
l_sv = 26
d_sv = 6

F_sv = round(sqrt(F_uklj**2 + F**2), 2)  # N
print(f"F_sv: {F_sv} N")

sigma_sv = round(0.125 * (F_sv * (l_sv - b)) / (0.1 * d_sv**3), 2)
sigma_sv_dop = 125  # N/mm^2
print(
    f"sigma_sv: {sigma_sv} N/mm^2 <= sigma_sv_dop = {sigma_sv_dop} N/mm^2",
    sigma_sv <= sigma_sv_dop,
)

p_sv_u = round(F_sv / (d_sv * b), 2)  # N/mm^2
p_sv_v = round(F_sv / (d_sv * (l_sv - b)), 2)  # N/mm^2
p_dop = 24  # N/mm^2
print(f"p_sv_u: {p_sv_u} N/mm^2 <= p_dop = {p_dop} N/mm^2", p_sv_u <= p_dop)
print(f"p_sv_v: {p_sv_v} N/mm^2 <= p_dop = {p_dop} N/mm^2", p_sv_v <= p_dop)

tau_a = round(2 * F_sv / (d_sv**2 * pi), 2)  # N/mm^2
tau_a_dop = 72  # N/mm^2
print(f"tau_a: {tau_a} N/mm^2 <= tau_a_dop = {tau_a_dop} N/mm^2", tau_a <= tau_a_dop)

# Ručica
alfa_nos = 40 * pi / 180  # rad
mi_ruc = 0.1
n_pol = 3

F_A = round(F * (tan(alfa_nos) + mi_ruc), 2)
F_A1 = round(n_pol * F_A * (tan(alfa_nos) + mi_ruc), 2)
print(f"F_A: {F_A} N")
print(f"F_A:' {F_A1} N")

R_ruc = 111  # mm (Ortlinghaus)
F_ruc = 150  # N
R = round(F_A1 * R_ruc / F_ruc, 2)
R_dop = 750  # mm
print(f"R: {R} mm <= R_dop = {R_dop} mm", R <= R_dop)
