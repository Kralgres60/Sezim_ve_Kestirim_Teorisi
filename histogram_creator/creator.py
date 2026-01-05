import os
import re
import matplotlib.pyplot as plt

# =========================
# AYARLAR
# =========================
INPUT_FILE = "gene_histg.txt"
OUTPUT_DIR = "histograms"

# Klasör yoksa oluştur
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Satır formatı için regex
pattern = re.compile(
    r"X index=(\d+)\s+Silence=(\d+)\s+"
    r"Ade=([\d.]+)\s+Gua=([\d.]+)\s+Thy=([\d.]+)\s+Cyt=([\d.]+)"
)

# =========================
# DOSYAYI OKU
# =========================
with open(INPUT_FILE, "r", encoding="utf-8") as f:
    lines = f.readlines()

# =========================
# HER SATIR İÇİN GRAFİK
# =========================
for line in lines:
    match = pattern.search(line)
    if not match:
        continue

    index, silence, ade, gua, thy, cyt = match.groups()
    values = [float(ade), float(gua), float(thy), float(cyt)]
    labels = ["Ade", "Gua", "Thy", "Cyt"]

    # Grafik
    plt.figure(figsize=(6, 4))
    plt.bar(labels, values)
    plt.title(f"Index {index}")
    plt.xlabel(f"Silence Number : {silence}")
    plt.ylabel("Values")
    plt.grid(axis="y", linestyle="--", alpha=0.6)

    # Kaydet
    filename = f"index_{index}_silence_{silence}.png"
    filepath = os.path.join(OUTPUT_DIR, filename)
    plt.tight_layout()
    plt.savefig(filepath)
    plt.close()

print("Tüm histogramlar başarıyla oluşturuldu.")
