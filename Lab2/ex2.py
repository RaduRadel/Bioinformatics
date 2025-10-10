S = "ATTGTCCCAATCTGTTG"

def existing_kmers_in_order(seq, k):
    seq = seq.upper()
    n = len(seq)
    total_positions = max(0, n - k + 1)

    seen = set()
    result = []

    for i in range(total_positions):
        km = seq[i:i+k]
        if km not in seen:
            count = 0
            for j in range(total_positions):
                if seq[j:j+k] == km:
                    count += 1
            rel = (count / total_positions) if total_positions > 0 else 0.0
            result.append((km, count, rel, rel * 100))
            seen.add(km)

    return result, total_positions

for k in (2, 3):
    rows, total = existing_kmers_in_order(S, k)
    print(f"\n{k}-nucleotides:")
    print(f"Total scanning windows = {total}")
    print(f"{'k-nucleotide':>6}  {'count':>5}  {'rel_freq':>10}  {'percent':>9}")
    print("-" * 36)
    for km, c, rel, pct in rows:
        print(f"{km:>6}  {c:5d}  {rel:10.4f}  {pct:8.2f}%")
