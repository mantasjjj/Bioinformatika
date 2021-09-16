import numpy as np
from Bio import SeqIO


# 1. Pateiktoje sekoje fasta formatu surastu visas start ir stop kodonų poras, tarp kurių nebutu stop kodono (ir tiesioginei sekai ir jos reverse komplementui).
# Susirandam visas imanomas kombinacijas is sekos
def get_triplets(seq):
    frames = []
    frames.append([seq[i:i + 3] for i in range(0, len(seq), 3)])
    frames.append([seq[i:i + 3] for i in range(1, len(seq), 3)])
    frames.append([seq[i:i + 3] for i in range(2, len(seq), 3)])
    frames.append([seq.reverse_complement()[i:i + 3] for i in range(0, len(seq), 3)])
    frames.append([seq.reverse_complement()[i:i + 3] for i in range(1, len(seq), 3)])
    frames.append([seq.reverse_complement()[i:i + 3] for i in range(2, len(seq), 3)])
    return frames


# Is sekos kombinacijos pasiemam kodonus, kurie prasideda ATG ir baigiasi TAA || TAG || TGA
def find_codon(seq_record):
    i = 0
    codon_list = []
    while i < len(seq_record):
        if seq_record[i] == 'ATG':
            start_pos = i
            j = i
            while j < len(seq_record):
                if seq_record[j] == 'TAA' or seq_record[j] == 'TAG' or seq_record[j] == 'TGA':
                    end_pos = j
                    codon_list.append(''.join(str(e) for e in seq_record[start_pos:end_pos + 1]))
                    i = j
                    break
                j += 1
        i += 1
    return codon_list


# Susirandam visus kodonus pagal seku kombinacijas
def find_codons(triplets):
    codons = []
    for triplet in triplets:
        codons.append(find_codon(triplet))
    return codons


# 2. Kiekvienam stop kodonui parinkti toliausiai nuo jo esanti start kodoną (su salyga, kad tarp ju nera kito stop kodono)


# 3. Atfiltruokite visus fragmentus ("tai butu baltymų koduojancios sekos"), kurie trumpesni nei 100 simbolių.
def find_shortest_codons(codons):
    shortest_codons = []
    for codon in codons:
        if len(codon) < 100:
            shortest_codons.append(codon)
    return shortest_codons


if __name__ == '__main__':
    record = SeqIO.read("data/bacterial3.fasta", "fasta")

    triplets = get_triplets(record.seq)

    # 1.
    codons = find_codons(triplets)
    print("All codons in a sequence:")
    codons = np.concatenate(codons)
    print(codons)

    # 3.
    shortest_codons = find_shortest_codons(codons)
    print("Codons that have a length less than 100:")
    # print(shortest_codons)
