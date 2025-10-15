import math

s = input("DNA: ").strip().upper()
valid = set("ATGC")

if not s or any(ch not in valid for ch in s):
    print("Error: DNA sequence must contain only A, T, G, or C.")
    exit()

A, T, G, C = s.count('A'), s.count('T'), s.count('G'), s.count('C')
n = len(s)
gc = (G + C) / n
na = 0.0001

tm_wallace = 4*(G + C) + 2*(A + T)
tm_salt = 81.5 + 16.6*math.log10(na) + 41*gc - 600/n

print(f"Tm (Wallace): {tm_wallace:.1f} °C")
print(f"Tm (Salt-adjusted): {tm_salt:.1f} °C")
