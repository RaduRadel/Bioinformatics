import os
import io
import gzip
import csv
import threading
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from collections import Counter

APP_TITLE = "FASTA DNA Base Percentages (Large-file friendly)"
DNA_ALPHABET = ("A", "C", "G", "T")
DNA_SET = set(DNA_ALPHABET)

class FastaCounterApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title(APP_TITLE)
        self.geometry("750x540")
        self.minsize(650, 500)

        # State
        self.file_path = None
        self.counter = Counter({b: 0 for b in DNA_ALPHABET})
        self.total_symbols = 0
        self._worker = None
        self._cancel = threading.Event()
        self._is_gz = False  # for progress style

        self._build_ui()

    # ---------- UI ----------
    def _build_ui(self):
        top = ttk.Frame(self, padding=12)
        top.pack(fill="x")

        self.pick_btn = ttk.Button(top, text="Choose FASTA file…", command=self.choose_file)
        self.pick_btn.pack(side="left")

        self.file_lbl = ttk.Label(top, text="No file selected")
        self.file_lbl.pack(side="left", padx=10)

        act = ttk.Frame(self, padding=(12, 0, 12, 12))
        act.pack(fill="x")

        self.start_btn = ttk.Button(act, text="Start", command=self.start_count, state="disabled")
        self.start_btn.pack(side="left")

        self.cancel_btn = ttk.Button(act, text="Cancel", command=self.cancel, state="disabled")
        self.cancel_btn.pack(side="left", padx=(8, 0))

        self.export_btn = ttk.Button(act, text="Export CSV…", command=self.export_csv, state="disabled")
        self.export_btn.pack(side="right")

        prog_box = ttk.Frame(self, padding=(12, 0, 12, 12))
        prog_box.pack(fill="x")
        ttk.Label(prog_box, text="Progress").pack(anchor="w")
        self.progress = ttk.Progressbar(prog_box, mode="determinate")
        self.progress.pack(fill="x")
        self.status_lbl = ttk.Label(prog_box, text="")
        self.status_lbl.pack(anchor="w", pady=(4, 0))

        table_box = ttk.Frame(self, padding=(12, 0, 12, 12))
        table_box.pack(fill="both", expand=True)

        cols = ("symbol", "count", "percentage")
        self.tree = ttk.Treeview(table_box, columns=cols, show="headings", height=14)
        for c, w, anchor in (("symbol", 120, "center"), ("count", 180, "e"), ("percentage", 180, "e")):
            self.tree.heading(c, text=c.capitalize())
            self.tree.column(c, width=w, anchor=anchor)
        self.tree.pack(side="left", fill="both", expand=True)

        vsb = ttk.Scrollbar(table_box, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscroll=vsb.set)
        vsb.pack(side="right", fill="y")

        style = ttk.Style(self)
        if "clam" in style.theme_names():
            style.theme_use("clam")

    # ---------- Actions ----------
    def choose_file(self):
        path = filedialog.askopenfilename(
            title="Select FASTA file",
            filetypes=[("FASTA files", "*.fasta *.fa *.fna *.fas *.txt *.gz"),
                       ("All files", "*.*")]
        )
        if not path:
            return
        try:
            size = os.path.getsize(path)
            self.file_lbl.config(text=f"{os.path.basename(path)}  —  {size/(1024*1024):.1f} MB")
        except Exception:
            self.file_lbl.config(text=os.path.basename(path))
        self.file_path = path
        self._is_gz = path.lower().endswith(".gz")
        self.start_btn.config(state="normal")
        self._clear_results()

    def start_count(self):
        if not self.file_path:
            messagebox.showwarning("No file", "Please choose a FASTA file first.")
            return
        self._clear_results()
        self._cancel.clear()
        self.pick_btn.config(state="disabled")
        self.start_btn.config(state="disabled")
        self.cancel_btn.config(state="normal")
        self.export_btn.config(state="disabled")

        # Progress setup
        if self._is_gz:
            # .gz: we cannot know uncompressed size beforehand; show indeterminate progress
            self.progress.config(mode="indeterminate")
            self.progress.start(60)  # ms per move
            self.status_lbl.config(text="Counting A, C, G, T… (decompressing)")
        else:
            self.progress.config(mode="determinate", value=0)
            self.status_lbl.config(text="Counting A, C, G, T…")

        self._worker = threading.Thread(target=self._count_worker, daemon=True)
        self._worker.start()
        self.after(200, self._poll_worker)

    def cancel(self):
        if self._worker and self._worker.is_alive():
            self._cancel.set()
            self.status_lbl.config(text="Cancelling…")
            self.cancel_btn.config(state="disabled")

    def _poll_worker(self):
        if self._worker and self._worker.is_alive():
            self.after(200, self._poll_worker)
        else:
            # Stop indeterminate bar if running
            if str(self.progress.cget("mode")) == "indeterminate":
                self.progress.stop()
                self.progress.config(mode="determinate", value=100)

            self.pick_btn.config(state="normal")
            self.cancel_btn.config(state="disabled")
            self.start_btn.config(state="normal" if self.file_path else "disabled")

            # Always show the four bases, even if zero
            self._populate_table()

            if not self._cancel.is_set():
                if self.total_symbols > 0:
                    self.export_btn.config(state="normal")
                    self.status_lbl.config(text=f"Done. Counted {self.total_symbols:,} A/C/G/T symbols.")
                else:
                    self.status_lbl.config(text="Done, but no A/C/G/T symbols were found in this file.")
            else:
                self.status_lbl.config(text="Cancelled.")

    def _clear_results(self):
        self.counter = Counter({b: 0 for b in DNA_ALPHABET})
        self.total_symbols = 0
        for row in self.tree.get_children():
            self.tree.delete(row)
        self.status_lbl.config(text="")
        self.progress.config(value=0)

    # ---------- File I/O helpers ----------
    def _open_maybe_gz(self, path):
        """
        Return (text_stream, progress_handle_or_none)
        - For .gz: TextIOWrapper over gzip.GzipFile; progress handle is gzip file (compressed bytes).
        - For plain files: text file; progress handle is the same file object (bytes).
        """
        if path.lower().endswith(".gz"):
            gf = gzip.open(path, "rb")
            txt = io.TextIOWrapper(gf, encoding="ascii", errors="ignore", newline="")
            return txt, gf
        else:
            f = open(path, "rt", encoding="ascii", errors="ignore", newline="")
            return f, f

    # ---------- Worker ----------
    def _count_worker(self):
        try:
            file_size = os.path.getsize(self.file_path) if not self._is_gz else None
            processed_bytes = 0

            txt, ph = self._open_maybe_gz(self.file_path)
            try:
                for line in txt:
                    if self._cancel.is_set():
                        return

                    # progress update (only for plain files where size is known)
                    if file_size:
                        try:
                            processed_bytes = ph.tell()
                            pct = min(100.0, (processed_bytes / file_size) * 100.0)
                            self._set_progress(pct, f"Counting… {pct:.1f}%")
                        except Exception:
                            pass

                    # FASTA header?
                    # Some files may have leading spaces; be tolerant.
                    if line.lstrip().startswith(">"):
                        continue

                    s = line.strip().upper()
                    if not s:
                        continue

                    # Map U->T (RNA → DNA) and count only A/C/G/T
                    s = s.replace("U", "T")
                    # Faster than per-char "in set" in pure Python loops for big lines:
                    for ch in s:
                        if ch in DNA_SET:
                            self.counter[ch] += 1
                            self.total_symbols += 1

            finally:
                try:
                    txt.close()
                except Exception:
                    pass

            # Final progress message for .gz
            if self._is_gz:
                self._set_progress(100.0, "Counting… 100.0%")

        except Exception as e:
            self._set_progress(0, f"Error: {e}")

    def _set_progress(self, pct, txt):
        # Safe UI update from worker
        self.after(0, lambda: (self.progress.config(value=pct),
                               self.status_lbl.config(text=txt)))

    # ---------- Results ----------
    def _populate_table(self):
        # Clear first (so repeated runs don't duplicate rows)
        for iid in self.tree.get_children():
            self.tree.delete(iid)

        for base in DNA_ALPHABET:
            cnt = self.counter.get(base, 0)
            pct = (cnt / self.total_symbols) * 100 if self.total_symbols else 0.0
            self.tree.insert("", "end", values=(base, f"{cnt:,}", f"{pct:.4f}%"))

    def export_csv(self):
        if self.total_symbols == 0:
            messagebox.showinfo("Nothing to export", "No A/C/G/T symbols were counted.")
            return
        path = filedialog.asksaveasfilename(
            title="Export CSV",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if not path:
            return
        try:
            with open(path, "w", newline="", encoding="utf-8") as fp:
                w = csv.writer(fp)
                w.writerow(["Symbol", "Count", "Percentage"])
                for iid in self.tree.get_children():
                    row = self.tree.item(iid)["values"]
                    w.writerow(row)
            messagebox.showinfo("Exported", f"Results exported to:\n{path}")
        except Exception as e:
            messagebox.showerror("Export error", str(e))


if __name__ == "__main__":
    app = FastaCounterApp()
    app.mainloop()
