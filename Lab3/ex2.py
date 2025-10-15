# sliding_tm_gui.py
import math
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import csv
from pathlib import Path

VALID = set("ATGC")

def read_fasta(path: str) -> str:
    seq = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq.append(line.strip())
    s = "".join(seq).upper()
    if not s or any(ch not in VALID for ch in s):
        raise ValueError("FASTA must contain only A,T,G,C (headers starting with '>' are ignored).")
    return s

def tm_wallace(w: str) -> float:
    a, t, g, c = w.count("A"), w.count("T"), w.count("G"), w.count("C")
    return 4*(g+c) + 2*(a+t)

def tm_salt(w: str, na_molar: float = 0.0001) -> float:
    n = len(w)
    gc = (w.count("G") + w.count("C")) / n
    return 81.5 + 16.6*math.log10(na_molar) + 41*gc - 600/n

def sliding_tm(seq: str, k: int, formula: str, na: float):
    if k < 1 or k > len(seq):
        raise ValueError("Window size must be between 1 and the sequence length.")
    f = tm_wallace if formula == "Wallace" else (lambda w: tm_salt(w, na))
    return [(i, seq[i:i+k], round(f(seq[i:i+k]), 2)) for i in range(len(seq) - k + 1)]

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Sliding-Window Tm")
        self.geometry("760x520")
        self.seq_path = None
        self.results = []

        controls = ttk.Frame(self, padding=10)
        controls.pack(fill="x")

        self.file_lbl = ttk.Label(controls, text="No FASTA selected")
        ttk.Button(controls, text="Choose FASTA...", command=self.choose_file).pack(side="right")
        self.file_lbl.pack(side="left", padx=(0,10))

        row2 = ttk.Frame(self, padding=10)
        row2.pack(fill="x")

        ttk.Label(row2, text="Window size:").grid(row=0, column=0, sticky="w")
        self.k_var = tk.StringVar(value="8")
        ttk.Entry(row2, textvariable=self.k_var, width=6).grid(row=0, column=1, padx=(4,15), sticky="w")

        ttk.Label(row2, text="Formula:").grid(row=0, column=2, sticky="w")
        self.formula_var = tk.StringVar(value="Wallace")
        formula_cb = ttk.Combobox(row2, textvariable=self.formula_var, values=["Wallace","Salt-adjusted"], width=14, state="readonly")
        formula_cb.grid(row=0, column=3, padx=(4,15), sticky="w")
        formula_cb.bind("<<ComboboxSelected>>", self.toggle_na)

        ttk.Label(row2, text="[Na⁺] (M):").grid(row=0, column=4, sticky="w")
        self.na_var = tk.StringVar(value="0.0001")
        self.na_entry = ttk.Entry(row2, textvariable=self.na_var, width=8)
        self.na_entry.grid(row=0, column=5, padx=(4,15), sticky="w")

        ttk.Button(row2, text="Compute", command=self.compute).grid(row=0, column=6, padx=(10,0))
        ttk.Button(row2, text="Export CSV", command=self.export_csv).grid(row=0, column=7, padx=(10,0))

        self.info = ttk.Label(self, text="", padding=(10,0))
        self.info.pack(fill="x")

        cols = ("start", "window", "tm")
        self.table = ttk.Treeview(self, columns=cols, show="headings", height=18)
        self.table.heading("start", text="Start index")
        self.table.heading("window", text="Window (k-mer)")
        self.table.heading("tm", text="Tm (°C)")
        self.table.column("start", width=100, anchor="e")
        self.table.column("window", width=420, anchor="w")
        self.table.column("tm", width=100, anchor="e")
        self.table.pack(fill="both", expand=True, padx=10, pady=10)

        self.toggle_na()

    def choose_file(self):
        path = filedialog.askopenfilename(
            title="Choose FASTA file",
            filetypes=[("FASTA files","*.fa *.fasta *.fna *.txt"), ("All files","*.*")]
        )
        if path:
            self.seq_path = path
            self.file_lbl.config(text=Path(path).name)

    def toggle_na(self, *_):
        is_salt = self.formula_var.get() == "Salt-adjusted"
        state = "normal" if is_salt else "disabled"
        self.na_entry.config(state=state)

    def compute(self):
        if not self.seq_path:
            messagebox.showwarning("Missing file", "Please choose a FASTA file.")
            return

        try:
            k = int(self.k_var.get())
        except ValueError:
            messagebox.showerror("Invalid window", "Window size must be an integer.")
            return

        try:
            na = float(self.na_var.get()) if self.formula_var.get()=="Salt-adjusted" else 0.0001
        except ValueError:
            messagebox.showerror("Invalid [Na+]", "Na+ must be a number in mol/L, e.g., 0.0001.")
            return

        try:
            seq = read_fasta(self.seq_path)
            self.results = sliding_tm(seq, k, self.formula_var.get(), na)
        except Exception as e:
            messagebox.showerror("Error", str(e))
            return

        for row in self.table.get_children():
            self.table.delete(row)
        for start, window, tm in self.results[:1000]:
            self.table.insert("", "end", values=(start, window, tm))

        self.info.config(text=f"Sequence length: {len(seq)} | Windows (k={k}): {len(self.results)} | "
                              f"Formula: {self.formula_var.get()}" + (f" | [Na+]={na} M" if self.formula_var.get()=="Salt-adjusted" else ""))

    def export_csv(self):
        if not self.results:
            messagebox.showinfo("No data", "Compute results first.")
            return
        out = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV","*.csv")], title="Save results as")
        if not out:
            return
        with open(out, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["start_index","window","Tm_C"])
            w.writerows(self.results)
        messagebox.showinfo("Saved", f"Exported {len(self.results)} rows to:\n{out}")

if __name__ == "__main__":
    App().mainloop()
