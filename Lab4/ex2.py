import os
from collections import Counter
from math import log, sqrt
import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter.scrolledtext import ScrolledText
import matplotlib.pyplot as plt

CODON_TO_AA = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}
AA3 = {
    "A":"Ala","R":"Arg","N":"Asn","D":"Asp","C":"Cys","E":"Glu","Q":"Gln","G":"Gly","H":"His",
    "I":"Ile","L":"Leu","K":"Lys","M":"Met","F":"Phe","P":"Pro","S":"Ser","T":"Thr","W":"Trp",
    "Y":"Tyr","V":"Val","*":"Stop"
}
STOPS = {"TAA","TAG","TGA"}

def read_fasta_concat(path: str) -> str:
    if not os.path.exists(path):
        raise FileNotFoundError(f"Not found: {path}")
    seq_lines = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line or line.startswith(">"):
                continue
            seq_lines.append(line.strip())
    seq = "".join(seq_lines).upper()
    seq = "".join(ch for ch in seq if ch in "ACGTU").replace("U","T")
    if not seq:
        raise ValueError(f"No sequence content parsed from: {path}")
    return seq

def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGT","TGCA")
    return seq.translate(comp)[::-1]

def find_orfs_on_strand(seq: str, min_codons: int = 60):
    orfs = []
    n = len(seq)
    for frame in (0,1,2):
        i = frame
        while i <= n - 3:
            if seq[i:i+3] == "ATG":
                j = i
                while j <= n - 3:
                    c = seq[j:j+3]
                    if c in STOPS:
                        cds = seq[i:j]
                        if len(cds) >= 3*min_codons:
                            orfs.append(cds)
                        i = j + 3
                        break
                    j += 3
                else:
                    i = n
            else:
                i += 3
    return orfs

def find_orfs_both_strands(seq: str, min_codons_plus=60, min_codons_minus=45):
    orfs_plus = find_orfs_on_strand(seq, min_codons=min_codons_plus)
    orfs_minus = find_orfs_on_strand(revcomp(seq), min_codons=min_codons_minus)
    return orfs_plus + orfs_minus

def count_codons_in_orfs(orfs):
    cnt = Counter()
    for cds in orfs:
        for i in range(0, len(cds)-2, 3):
            codon = cds[i:i+3]
            if len(codon) == 3 and set(codon) <= set("ACGT"):
                cnt[codon] += 1
    return cnt

def fallback_count_frame0(seq: str):
    cnt = Counter()
    for i in range(0, len(seq)-2, 3):
        codon = seq[i:i+3]
        if len(codon) == 3 and set(codon) <= set("ACGT"):
            cnt[codon] += 1
    return cnt

def top_n(counter: Counter, n=10):
    return sorted(counter.items(), key=lambda kv: (-kv[1], kv[0]))[:n]

def codon_freqs(counter: Counter):
    total = sum(counter.values())
    if total == 0:
        return {k: 0.0 for k in CODON_TO_AA}
    return {k: counter.get(k, 0)/total for k in CODON_TO_AA}

def amino_acid_top3_from_codon_counts(cnt: Counter):
    aa_counts = Counter()
    for codon, c in cnt.items():
        aa = CODON_TO_AA.get(codon)
        if aa and aa != "*":
            aa_counts[aa] += c
    return top_n(aa_counts, 3)

def cosine_similarity(freq_a: dict, freq_b: dict):
    keys = list(CODON_TO_AA.keys())
    a = [freq_a.get(k, 0.0) for k in keys]
    b = [freq_b.get(k, 0.0) for k in keys]
    dot = sum(x*y for x,y in zip(a,b))
    na = (sum(x*x for x in a)) ** 0.5
    nb = (sum(y*y for y in b)) ** 0.5
    return 0.0 if na == 0 or nb == 0 else dot/(na*nb)

def kl_divergence(p: dict, q: dict, eps=1e-12):
    s = 0.0
    for k in CODON_TO_AA:
        pk = max(p.get(k, 0.0), eps)
        qk = max(q.get(k, 0.0), eps)
        s += pk * ( (pk/qk) and ( (pk/qk) != 0) and ( (pk/qk) ).bit_length() )
    return float('nan')

def plot_top10(counter: Counter, title: str, outfile: str):
    top10 = top_n(counter, 10)
    labels = [k for k,_ in top10]
    values = [v for _,v in top10]
    plt.figure()
    plt.bar(labels, values)
    plt.title(title)
    plt.xlabel("Codon")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(outfile, dpi=150)
    plt.show()

def compare_common_top10(cnt_a: Counter, cnt_b: Counter, label_a: str, label_b: str, outfile: str):
    a_top = {k for k,_ in top_n(cnt_a, 10)}
    b_top = {k for k,_ in top_n(cnt_b, 10)}
    common = sorted(a_top & b_top)
    if not common:
        common = sorted(list(a_top))[:5] + sorted(list(b_top))[:5]
    vals_a = [cnt_a.get(c, 0) for c in common]
    vals_b = [cnt_b.get(c, 0) for c in common]
    x = list(range(len(common)))
    width = 0.4
    xa = [i - width/2 for i in x]
    xb = [i + width/2 for i in x]
    plt.figure()
    plt.bar(xa, vals_a, width=width, label=label_a)
    plt.bar(xb, vals_b, width=width, label=label_b)
    plt.xticks(x, common)
    plt.xlabel("Common top codons")
    plt.ylabel("Count")
    plt.title(f"Common top codons: {label_a} vs {label_b}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outfile, dpi=150)
    plt.show()
    return common

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Codon Frequency Comparator — COVID-19 vs Influenza")
        self.geometry("760x520")

        self.path_covid = tk.StringVar()
        self.path_flu = tk.StringVar()

        # Row 0: file pickers
        tk.Label(self, text="COVID-19 FASTA:").grid(row=0, column=0, sticky="w", padx=8, pady=8)
        tk.Entry(self, textvariable=self.path_covid, width=70).grid(row=0, column=1, padx=4, pady=8)
        tk.Button(self, text="Choose...", command=self.choose_covid).grid(row=0, column=2, padx=4, pady=8)

        tk.Label(self, text="Influenza FASTA:").grid(row=1, column=0, sticky="w", padx=8, pady=4)
        tk.Entry(self, textvariable=self.path_flu, width=70).grid(row=1, column=1, padx=4, pady=4)
        tk.Button(self, text="Choose...", command=self.choose_flu).grid(row=1, column=2, padx=4, pady=4)

        tk.Button(self, text="Analyze", command=self.run_analysis).grid(row=2, column=2, padx=8, pady=8, sticky="e")

        tk.Label(self, text="Results:").grid(row=3, column=0, sticky="w", padx=8)
        self.out = ScrolledText(self, width=90, height=20)
        self.out.grid(row=4, column=0, columnspan=3, padx=8, pady=6, sticky="nsew")

        self.grid_rowconfigure(4, weight=1)
        self.grid_columnconfigure(1, weight=1)

    def choose_covid(self):
        path = filedialog.askopenfilename(title="Select COVID-19 FASTA", filetypes=[("FASTA", "*.fa *.fasta *.fna *.fa.gz *.fasta.gz"), ("All files", "*.*")])
        if path:
            self.path_covid.set(path)

    def choose_flu(self):
        path = filedialog.askopenfilename(title="Select Influenza FASTA", filetypes=[("FASTA", "*.fa *.fasta *.fna *.fa.gz *.fasta.gz"), ("All files", "*.*")])
        if path:
            self.path_flu.set(path)

    def log(self, text):
        self.out.insert(tk.END, text + "\n")
        self.out.see(tk.END)
        self.update_idletasks()

    def run_analysis(self):
        self.out.delete("1.0", tk.END)
        covid = self.path_covid.get().strip()
        flu = self.path_flu.get().strip()
        if not covid or not flu:
            messagebox.showerror("Missing files", "Please choose both FASTA files first.")
            return
        try:
            seq_covid = read_fasta_concat(covid)
            seq_flu   = read_fasta_concat(flu)
        except Exception as e:
            messagebox.showerror("Read error", str(e))
            return

        self.log("[1/4] Finding ORFs (both strands)...")
        orfs_covid = find_orfs_both_strands(seq_covid, min_codons_plus=60, min_codons_minus=45)
        orfs_flu   = find_orfs_both_strands(seq_flu,   min_codons_plus=60, min_codons_minus=45)

        if not orfs_covid:
            self.log("  No ORFs found for COVID; using frame-0 fallback.")
            cnt_covid = fallback_count_frame0(seq_covid)
        else:
            cnt_covid = count_codons_in_orfs(orfs_covid)

        if not orfs_flu:
            self.log("  No ORFs found for Influenza; using frame-0 fallback.")
            cnt_flu = fallback_count_frame0(seq_flu)
        else:
            cnt_flu = count_codons_in_orfs(orfs_flu)

        # a) Top 10 COVID
        self.log("[2/4] Plotting Top-10 codons for COVID-19...")
        plot_top10(cnt_covid, "Top 10 codons – SARS-CoV-2", "top10_covid.png")
        self.log("  Saved: top10_covid.png")

        # b) Top 10 Influenza
        self.log("[3/4] Plotting Top-10 codons for Influenza...")
        plot_top10(cnt_flu, "Top 10 codons – Influenza", "top10_influenza.png")
        self.log("  Saved: top10_influenza.png")

        # c) overlap chart
        self.log("[4/4] Comparing common top codons...")
        common = compare_common_top10(cnt_covid, cnt_flu, "SARS-CoV-2", "Influenza", "top_common.png")
        self.log("  Saved: top_common.png")
        if common:
            self.log("Common codons in both Top-10 lists: " + ", ".join(common))

        # d) top 3 AAs per genome
        covid_top3_aa = amino_acid_top3_from_codon_counts(cnt_covid)
        flu_top3_aa   = amino_acid_top3_from_codon_counts(cnt_flu)

        self.log("\nTop 3 amino acids – SARS-CoV-2:")
        for aa, c in covid_top3_aa:
            self.log(f"  {aa} ({AA3.get(aa, aa)}): {c}")

        self.log("\nTop 3 amino acids – Influenza:")
        for aa, c in flu_top3_aa:
            self.log(f"  {aa} ({AA3.get(aa, aa)}): {c}")

        f_covid = codon_freqs(cnt_covid)
        f_flu   = codon_freqs(cnt_flu)
        cos = cosine_similarity(f_covid, f_flu)
        self.log(f"\nCosine similarity (0..1, higher = more similar): {cos:.4f}")

        messagebox.showinfo("Done", "Analysis complete.\nCharts saved: top10_covid.png, top10_influenza.png, top_common.png")

if __name__ == "__main__":
    App().mainloop()
