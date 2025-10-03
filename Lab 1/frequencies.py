from collections import Counter

def alphabet(sequence: str):
    n = len(sequence)
    freq = Counter(sequence)
    result = {symbol: (count / n) * 100 for symbol, count in freq.items()}
    return dict(sorted(result.items()))

if __name__ == "__main__":
    S = "ATTTCGCCGATA"
    alphabet_freqs = alphabet(S)
    for symbol, pct in alphabet_freqs.items():
        print(f"{symbol} = {pct:.2f}%")