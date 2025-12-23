import subprocess
from pathlib import Path
from collections import Counter
import matplotlib.pyplot as plt

DNA = ["A", "C", "G", "T"]

# ----------------------------
# DOSYA YOLLARI
# ----------------------------
ORIGINAL_TXT = "tsg101_formatted.txt"
MASKED_TXT   = "tsg101_modified.txt"
HOMOLOGS_FASTA = "homologs.fasta"

OUTDIR = Path("out_profile_hmm")
OUTDIR.mkdir(exist_ok=True)

# ----------------------------
# YARDIMCI FONKSİYONLAR
# ----------------------------
def run(cmd):
    print("[CMD]", " ".join(cmd))
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        print(r.stderr)
        raise RuntimeError("Komut başarısız")

def read_txt(path):
    return Path(path).read_text().replace("\n", "").replace("\r", "").upper()

def write_fasta(seq, path, name="query"):
    with open(path, "w") as f:
        f.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")

# ----------------------------
# 1️⃣ MSA
# ----------------------------
msa_fasta = OUTDIR / "msa.fasta"
run(["mafft", "--auto", HOMOLOGS_FASTA])
msa_fasta.write_text(
    subprocess.check_output(["mafft", "--auto", HOMOLOGS_FASTA], text=True)
)

# ----------------------------
# 2️⃣ PROFILE-HMM
# ----------------------------
hmm_path = OUTDIR / "profile.hmm"
run(["hmmbuild", "--dna", str(hmm_path), str(msa_fasta)])

# ----------------------------
# 3️⃣ HMMALIGN
# ----------------------------
masked = read_txt(MASKED_TXT).replace("X", "N")
query_fa = OUTDIR / "query.fasta"
write_fasta(masked, query_fa)

sto_path = OUTDIR / "aligned.sto"
run(["hmmalign", "-o", str(sto_path), str(hmm_path), str(query_fa)])

# ----------------------------
# 4️⃣ HMM DOSYASINDAN MATCH EMİSYONLARI
# ----------------------------
def parse_hmm_emissions(hmm_file):
    emissions = {}
    with open(hmm_file) as f:
        lines = f.readlines()

    i = 0
    k = 0
    while i < len(lines):
        if lines[i].startswith("HMM"):
            i += 2
            continue

        if lines[i].strip() and lines[i].split()[0].isdigit():
            parts = lines[i].split()
            k += 1
            probs = list(map(float, parts[1:5]))
            total = sum(probs)
            emissions[k] = {DNA[j]: probs[j]/total for j in range(4)}
            i += 3
        else:
            i += 1

    return emissions

match_emissions = parse_hmm_emissions(hmm_path)

# ----------------------------
# 5️⃣ STOCKHOLM → MATCH STATE MAP
# ----------------------------
def parse_stockholm(sto):
    seq = ""
    rf  = ""
    with open(sto) as f:
        for line in f:
            if line.startswith("#=GC RF"):
                rf += line.split()[-1]
            elif not line.startswith("#") and not line.startswith("//"):
                parts = line.split()
                if len(parts) == 2:
                    seq += parts[1]
    return seq, rf

aligned_seq, rf = parse_stockholm(sto_path)

pos_to_M = {}
qpos = 0
mk = 0

for i in range(len(aligned_seq)):
    if rf[i] != ".":
        mk += 1
    if aligned_seq[i] not in "-.":
        if rf[i] != ".":
            pos_to_M[qpos] = mk
        qpos += 1

# ----------------------------
# 6️⃣ X DOLDUR
# ----------------------------
original = read_txt(ORIGINAL_TXT)
masked = read_txt(MASKED_TXT)

pred = list(masked)
for i, b in enumerate(pred):
    if b == "X":
        mk = pos_to_M.get(i)
        if mk in match_emissions:
            probs = match_emissions[mk]
            pred[i] = max(probs, key=probs.get)
        else:
            pred[i] = "A"  # fallback

predicted = "".join(pred)
(Path(OUTDIR) / "predicted.txt").write_text(predicted)

# ----------------------------
# 7️⃣ DEĞERLENDİRME
# ----------------------------
correct = 0
total = 0
hist = Counter()

for i, b in enumerate(masked):
    if b == "X":
        total += 1
        hist[(original[i], predicted[i])] += 1
        if original[i] == predicted[i]:
            correct += 1

acc = correct / total

print("Toplam X:", total)
print("Doğru:", correct)
print("Accuracy:", round(acc, 4))

# ----------------------------
# 8️⃣ HİSTOGRAM
# ----------------------------
plt.figure(figsize=(12,5))
plt.bar(
    [f"{a}->{b}" for (a,b) in hist],
    hist.values()
)
plt.xticks(rotation=45)
plt.title("Profile-HMM Posterior (MAP) Tahmin Histogramı")
plt.tight_layout()
plt.savefig(OUTDIR / "histogram.png", dpi=200)
plt.show()
