from collections import Counter
from itertools import product

S = "ATTGTCCCAATCTGTTG"
alphabet = set("ATCG")

def kmer_stats(seq, k):
    seq = seq.upper()
    counts = Counter(seq[i:i+k] for i in range(len(seq)-k+1)
                     if set(seq[i:i+k]) <= alphabet)
    total = sum(counts.values())

    data = []
    for km in map(''.join, product("ATCG", repeat=k)):
        c = counts.get(km, 0)
        rel = (c / total) if total else 0.0
        data.append((km, c, rel, rel*100))

    data.sort(key=lambda x: (-x[1], x[0]))
    return data

for k in (2, 3):
    print(f"\n{k}-mers (k={k}):")
    print(f"{'k-mer':>6}  {'count':>5}  {'rel_freq':>10}  {'percent':>9}")
    print("-" * 36)
    for km, c, rel, pct in kmer_stats(S, k):
        print(f"{km:>6}  {c:5d}  {rel:10.4f}  {pct:8.2f}%")
