#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import threading
from collections import deque, Counter
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

# Matplotlib (no seaborn; default styles/colors)
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

WINDOW_SIZE = 30
VALID_BASES = {"A", "C", "G", "T"}

class FastaSlidingApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("FASTA Sliding Window Frequencies (k=30)")
        self.geometry("900x650")

        # Top controls
        top = tk.Frame(self)
        top.pack(fill="x", padx=10, pady=8)

        self.path_var = tk.StringVar()
        tk.Label(top, text="FASTA file:").pack(side="left")
        tk.Entry(top, textvariable=self.path_var, width=60).pack(side="left", padx=6)
        tk.Button(top, text="Browse…", command=self.browse_file).pack(side="left", padx=6)
        tk.Button(top, text="Analyze", command=self.start_analysis).pack(side="left")

        # Progress + status
        prog_frame = tk.Frame(self)
        prog_frame.pack(fill="x", padx=10, pady=6)
        self.progress = ttk.Progressbar(prog_frame, mode="determinate")
        self.progress.pack(fill="x")
        self.status_var = tk.StringVar(value="Select a FASTA file and click Analyze.")
        tk.Label(self, textvariable=self.status_var, anchor="w").pack(fill="x", padx=10)

        # Matplotlib figure
        fig = Figure(figsize=(8.8, 4.8))
        self.ax = fig.add_subplot(111)
        self.ax.set_title("Sliding Window Nucleotide Frequencies (window = 30)")
        self.ax.set_xlabel("Window index (end at position i)")
        self.ax.set_ylabel("Relative frequency")

        self.canvas = FigureCanvasTkAgg(fig, master=self)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill="both", expand=True, padx=10, pady=10)

        # Data holders
        self.freq_A = []
        self.freq_C = []
        self.freq_G = []
        self.freq_T = []

    def browse_file(self):
        path = filedialog.askopenfilename(
            title="Select FASTA file",
            filetypes=[("FASTA files", "*.fasta *.fa *.fna *.ffn *.faa *.frn *.fst *.txt"), ("All files", "*.*")]
        )
        if path:
            self.path_var.set(path)

    def start_analysis(self):
        path = self.path_var.get().strip()
        if not path or not os.path.isfile(path):
            messagebox.showerror("Error", "Please select a valid FASTA file.")
            return

        # Reset UI & data
        self.progress["value"] = 0
        self.status_var.set("Analyzing…")
        self.freq_A.clear(); self.freq_C.clear(); self.freq_G.clear(); self.freq_T.clear()
        self.ax.clear()
        self.ax.set_title("Sliding Window Nucleotide Frequencies (window = 30)")
        self.ax.set_xlabel("Window index (end at position i)")
        self.ax.set_ylabel("Relative frequency")
        self.canvas.draw()

        # Run in a background thread to keep UI responsive
        t = threading.Thread(target=self.analyze_file, args=(path,), daemon=True)
        t.start()

    def analyze_file(self, path):
        try:
            file_size = os.path.getsize(path)
        except Exception:
            file_size = 0

        # Sliding window state
        window = deque(maxlen=WINDOW_SIZE)
        counts = Counter()   # counts for A,C,G,T within the current window
        valid_total = 0      # number of valid (A/C/G/T) chars in current window

        bytes_read = 0
        produced_points = 0

        # Stream through file efficiently
        try:
            with open(path, "r", encoding="utf-8", errors="ignore") as f:
                for raw_line in f:
                    bytes_read += len(raw_line.encode("utf-8", errors="ignore"))
                    line = raw_line.strip()
                    if not line or line.startswith(">"):
                        # Skip FASTA headers and blank lines
                        self._update_progress(bytes_read, file_size)
                        continue

                    for ch in line.upper():
                        # Push new char into window
                        old = None
                        if len(window) == WINDOW_SIZE:
                            # About to drop oldest element
                            old = window[0]
                        window.append(ch)

                        # Update counts for leaving char (if any)
                        if old is not None and old in VALID_BASES:
                            counts[old] -= 1
                            if counts[old] == 0:
                                del counts[old]
                            valid_total -= 1

                        # Update counts for incoming char
                        if ch in VALID_BASES:
                            counts[ch] += 1
                            valid_total += 1

                        # Once window is full, record frequencies
                        if len(window) == WINDOW_SIZE:
                            # denominator = count of valid bases in this window
                            denom = max(1, valid_total)
                            self.freq_A.append(counts.get("A", 0) / denom)
                            self.freq_C.append(counts.get("C", 0) / denom)
                            self.freq_G.append(counts.get("G", 0) / denom)
                            self.freq_T.append(counts.get("T", 0) / denom)
                            produced_points += 1

                    self._update_progress(bytes_read, file_size)

        except Exception as e:
            self._set_status(f"Error: {e}")
            self._update_progress(file_size, file_size)
            return

        # Plot only after the vectors are fully computed
        self._plot_results()
        if produced_points == 0:
            self._set_status("Done. (No complete 30-base windows found.)")
        else:
            self._set_status(f"Done. Generated {produced_points} window points.")

    def _plot_results(self):
        # Update the plot on the Tk main thread
        def _update_plot():
            self.ax.clear()
            self.ax.set_title("Sliding Window Nucleotide Frequencies (window = 30)")
            self.ax.set_xlabel("Window index (end at position i)")
            self.ax.set_ylabel("Relative frequency (A,C,G,T)")

            x = range(1, len(self.freq_A) + 1)
            # 4 signals (let matplotlib choose default colors)
            self.ax.plot(x, self.freq_A, label="A")
            self.ax.plot(x, self.freq_C, label="C")
            self.ax.plot(x, self.freq_G, label="G")
            self.ax.plot(x, self.freq_T, label="T")
            self.ax.set_ylim(0.0, 1.0)
            self.ax.legend(loc="upper right")
            self.ax.grid(True, alpha=0.3)

            self.canvas.draw()

        self.after(0, _update_plot)

    def _update_progress(self, bytes_read, file_size):
        # Update progress bar safely from worker thread
        def _upd():
            if file_size > 0:
                val = max(0, min(100, 100 * bytes_read / file_size))
                self.progress["value"] = val
            else:
                # Indeterminate if size is unknown
                self.progress["value"] = 0
            self.update_idletasks()
        self.after(0, _upd)

    def _set_status(self, text):
        def _upd():
            self.status_var.set(text)
        self.after(0, _upd)


if __name__ == "__main__":
    try:
        import matplotlib  # ensure installed
    except ImportError:
        raise SystemExit("Please install matplotlib: pip install matplotlib")
    app = FastaSlidingApp()
    app.mainloop()
