import numpy as np
import matplotlib.pyplot as plt
import re
import math

# ========== 1. METNİ DOSYADAN OKU ==========
with open("input.txt", "r", encoding="utf8") as f:
    text = f.read()

# ========== 2. SESLİ / SESSİZ ==========
vowels = "aeıioöuüAEIİOÖUÜ"

def is_vowel(c):
    return c in vowels

clean = re.sub(r"[^a-zA-ZçğıöşüÇĞİÖŞÜ]", "", text)

num_vowels = sum(1 for ch in clean if is_vowel(ch))
num_cons   = sum(1 for ch in clean if not is_vowel(ch))
total_chars = len(clean)

p_V = num_vowels / total_chars
p_C = num_cons / total_chars

pi = np.array([p_V, p_C])

print("=== METİNDEN ELDE EDİLEN ORANLAR ===")
print(f"P(V)= {p_V:.4f}, P(C)= {p_C:.4f}\n")

# =================================================
# 3. ASİMETRİK PROPOSAL MATRİSİ ve MH ALGORİTMASI
# =================================================
# q[i,j] = i -> j öneri olasılığı
# i=0: V, i=1: C
q = np.array([
    [0.1, 0.9],   # V -> V %10, V -> C %90
    [0.9, 0.1]    # C -> V %90, C -> C %10
])

def metropolis_hastings(num_steps, pi, q):
    """
    2 durumlu (0=V, 1=C) Metropolis-Hastings.
    pi: hedef dağılım [pi_V, pi_C]
    q : proposal matrisi (asimetrik olabilir)
    """
    chain = np.zeros(num_steps, dtype=int)
    x = 0                       # başlangıç: V
    rng = np.random.default_rng()

    for i in range(num_steps):
        # 1) Proposal: q[x] satırına göre y seç
        y = rng.choice([0, 1], p=q[x])

        # 2) Kabul olasılığı:
        #    alpha = min(1, (pi[y]*q[y,x]) / (pi[x]*q[x,y]))
        if x == y:
            alpha = 1.0  # aynı state'e öneri → kabul et
        else:
            alpha = min(1.0, (pi[y] * q[y, x]) / (pi[x] * q[x, y]))

        u = rng.random()
        if u < alpha:
            x = y  # kabul

        chain[i] = x

    return chain

# ========== 4. MH ÇALIŞTIR ==========
num_steps = 1000
chain = metropolis_hastings(num_steps, pi, q)

p_sim_V = np.mean(chain == 0)
p_sim_C = np.mean(chain == 1)

print("=== METROPOLIS-HASTINGS SONUÇLARI ===")
print(f"P_sim(V) = {p_sim_V:.4f}")
print(f"P_sim(C) = {p_sim_C:.4f}")

print("\nKarşılaştırma:")
print(f"Hedef   P(V)={p_V:.4f}, P(C)={p_C:.4f}")
print(f"MH Sim  P(V)={p_sim_V:.4f}, P(C)={p_sim_C:.4f}\n")

# =======================================================
# 5. NORMAL DAĞILIM (BOX–MULLER) + HISTOGRAM + CDF
# =======================================================
def box_muller(n):
    z = []
    while len(z) < n:
        u1, u2 = np.random.rand(), np.random.rand()
        if u1 == 0:
            continue
        r = math.sqrt(-2 * math.log(u1))
        theta = 2 * math.pi * u2
        z.append(r * math.cos(theta))
        if len(z) < n:
            z.append(r * math.sin(theta))
    return np.array(z)


# =======================================================
# 6. CONVERGENCE (YAKINSAMA) GRAFİĞİ
# =======================================================
running_mean = np.cumsum(chain == 0) / np.arange(1, num_steps + 1)

plt.figure(figsize=(8,5))
plt.plot(running_mean, label="P_sim(V) - Running Mean")
plt.axhline(p_V, color="red", linestyle="--", label="Hedef P(V)")
plt.title("Metropolis-Hastings Yakınsama Grafiği")
plt.xlabel("İterasyon")
plt.ylabel("P(V)")
plt.grid(alpha=0.1)
plt.legend()
plt.show()
