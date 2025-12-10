import pandas as pd
import matplotlib.pyplot as plt

# Excel dosyasını yükle
df = pd.read_excel("histogram.xlsx")

# Sütun isimlerini düzenle
df.columns = ["Silece", "SliceNum", "Ade", "Gua", "Thy", "Cyt"]

# Her Silece (Index) için
for index in df["Silece"].unique():

    # Bu index'e ait tüm satırlar (25, 50, 75, 100)
    subset = df[df["Silece"] == index]

    for _, row in subset.iterrows():

        labels = ["Ade", "Gua", "Thy", "Cyt"]
        values = [row["Ade"], row["Gua"], row["Thy"], row["Cyt"]]

        plt.figure(figsize=(8, 5))
        plt.bar(labels, values, color="steelblue")

        # Başlık
        plt.title(f"Index {int(index)}", fontsize=16)

        # Alt yazı
        plt.xlabel(f"Silece Number : {int(row['SliceNum'])}", fontsize=12)

        plt.ylabel("Values")
        plt.grid(axis='y', linestyle='--', alpha=0.5)

        # Dosya adı
        filename = f"hist_{int(index)}_{int(row['SliceNum'])}.png"

        plt.savefig(filename, dpi=200, bbox_inches="tight")
        plt.close()

        print(f"{filename} kaydedildi.")
