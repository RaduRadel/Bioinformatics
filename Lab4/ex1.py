import sys

GENETIC_CODE = {
    "UUU":"F","UUC":"F","UUA":"L","UUG":"L",
    "UCU":"S","UCC":"S","UCA":"S","UCG":"S",
    "UAU":"Y","UAC":"Y","UAA":"*","UAG":"*",
    "UGU":"C","UGC":"C","UGA":"*","UGG":"W",

    "CUU":"L","CUC":"L","CUA":"L","CUG":"L",
    "CCU":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAU":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "CGU":"R","CGC":"R","CGA":"R","CGG":"R",

    "AUU":"I","AUC":"I","AUA":"I","AUG":"M",
    "ACU":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAU":"N","AAC":"N","AAA":"K","AAG":"K",
    "AGU":"S","AGC":"S","AGA":"R","AGG":"R",

    "GUU":"V","GUC":"V","GUA":"V","GUG":"V",
    "GCU":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAU":"D","GAC":"D","GAA":"E","GAG":"E",
    "GGU":"G","GGC":"G","GGA":"G","GGG":"G",
}

AA3 = {
    "A":"Ala","R":"Arg","N":"Asn","D":"Asp","C":"Cys","E":"Glu","Q":"Gln","G":"Gly","H":"His",
    "I":"Ile","L":"Leu","K":"Lys","M":"Met","F":"Phe","P":"Pro","S":"Ser","T":"Thr","W":"Trp",
    "Y":"Tyr","V":"Val","*":"Stop"
}

def clean_and_validate(seq: str) -> str:
    s = "".join(ch for ch in seq.upper() if not ch.isspace())
    allowed = set("ACGTU")
    if not s or any(ch not in allowed for ch in s):
        bad = sorted({ch for ch in s if ch not in allowed})
        raise ValueError(f"Invalid characters: {''.join(bad) or 'none'}. Allowed: A C G T U.")
    return s.replace("T","U")

def translate_coding_region(rna: str) -> str:
    start = rna.find("AUG")
    if start == -1:
        raise ValueError("No start codon (AUG/ATG) found.")
    aa = []
    for i in range(start, len(rna)-2, 3):
        codon = rna[i:i+3]
        a = GENETIC_CODE.get(codon, "?")
        if a == "?" or len(codon) < 3:
            break
        if a == "*":
            aa.append("*")
            return "".join(aa)
        aa.append(a)
    raise ValueError("No in-frame stop codon (UAA/UAG/UGA) found before sequence end.")

def to_three_letter(p1: str) -> str:
    return "-".join(AA3.get(a, "Xaa") for a in p1)

if __name__ == "__main__":
    seq = input("DNA/RNA sequence:\n> ")
    try:
        rna = clean_and_validate(seq)
        protein = translate_coding_region(rna)
        print("\nProtein (1-letter):", protein)
        print("Protein (3-letter):", to_three_letter(protein))
    except ValueError as e:
        print("Error:", e)