import numpy as np
import re

# ========== 1. METNİ DOSYADAN OKU ==========
with open("input.txt", "r", encoding="utf8") as f:
    text = f.read()

# ========== 2. SESLİ / SESSİZ TESPİTİ ==========
vowels = "aeıioöuüAEIİOÖUÜ"

def is_vowel(c):
    return c in vowels

# Sadece harfleri tut (Türkçe için)
clean = re.sub(r"[^a-zA-ZçğıöşüÇĞİÖŞÜ]", "", text)

# Sesli / sessiz sayıları
num_vowels = sum(1 for ch in clean if is_vowel(ch))
num_cons   = sum(1 for ch in clean if not is_vowel(ch))
total_chars = len(clean)

p_V = num_vowels / total_chars
p_C = num_cons / total_chars

print("=== METİNDEN ELDE EDİLEN HARC ORANLARI ===")
print(f"Sesli (Vowel)   : {num_vowels} (%{p_V*100:.2f})")
print(f"Sessiz (Cons.)  : {num_cons} (%{p_C*100:.2f})")
print(f"Toplam Harf     : {total_chars}")
print(f"Hedef dağılım π = [P(V)={p_V:.4f}, P(C)={p_C:.4f}]\n")

# ========== 3. METROPOLIS HEDEF DAĞILIMI ==========
# State 0: V  (sesli)
# State 1: C  (sessiz)
pi = np.array([p_V, p_C])  # hedef dağılım

# ========== 4. METROPOLIS ALGORİTMASI ==========
def metropolis(num_steps, pi):
    """
    2 durumlu (0=V, 1=C) Metropolis örnekleyici.
    pi: hedef dağılım [pi_V, pi_C]
    """
    states = []
    # Başlangıç durumu: 0 (V) veya 1 (C) - fark etmez
    x = 0
    rng = np.random.default_rng()

    for _ in range(num_steps):
        # Simetrik proposal: diğer durumu öner
        y = 1 - x  # 0 ise 1, 1 ise 0

        # Kabul olasılığı alpha = min(1, pi[y]/pi[x])
        alpha = min(1.0, pi[y] / pi[x])

        u = rng.random()
        if u < alpha:
            x = y  # kabul

        states.append(x)

    return np.array(states)

num_steps = 20
chain = metropolis(num_steps, pi)

# ========== 5. SİMÜLASYON ORANLARI ==========
count_V = np.sum(chain == 0)
count_C = np.sum(chain == 1)

p_sim_V = count_V / num_steps
p_sim_C = count_C / num_steps

print("=== METROPOLIS SİMÜLASYON SONUÇLARI ===")
print(f"Toplam adım            : {num_steps}")
print(f"Simülasyonda V (state 0) sayısı : {count_V} (%{p_sim_V*100:.2f})")
print(f"Simülasyonda C (state 1) sayısı : {count_C} (%{p_sim_C*100:.2f})")

print("\nKarşılaştırma:")
print(f"Hedef   P(V) = {p_V:.4f}, P(C) = {p_C:.4f}")
print(f"Metropolis P_sim(V) = {p_sim_V:.4f}, P_sim(C) = {p_sim_C:.4f}")
