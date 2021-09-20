import numpy as np
from Bio import SeqIO
from Bio.Data.CodonTable import unambiguous_dna_by_name


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
def find_farthest_codon_for_each_stop_codon(seq_record):
    i = 0
    codon_list = []
    while i < len(seq_record):
        if seq_record[i] == 'TAA' or seq_record[i] == 'TAG' or seq_record[i] == 'TGA':
            start_pos = i
            j = i + 1
            end_pos = -1
            while j < len(seq_record):
                if seq_record[j] == 'ATG':
                    end_pos = j
                elif seq_record[j] == 'TAA' or seq_record[j] == 'TAG' or seq_record[j] == 'TGA':
                    if end_pos != -1:
                        codon_list.append(''.join(str(e) for e in seq_record[start_pos:end_pos + 1]))
                    i = j
                    break
                j += 1
        i += 1
    return codon_list


def find_farthest_starts(triplets):
    starts = []
    for triplet in triplets:
        starts.append(find_farthest_codon_for_each_stop_codon(triplet))
    return starts


# 3. Atfiltruokite visus fragmentus ("tai butu baltymų koduojancios sekos"), kurie trumpesni nei 100 simbolių.
def find_shortest_codons(codons):
    shortest_codons = []
    for codon in codons:
        if len(codon) < 100:
            shortest_codons.append(codon)
    return shortest_codons


# 4. Parasykite funkcijas, kurios ivertintu kodonu ir dikodonu daznius(visi imanomi kodonai/dikodonai ir ju
# atitinkamas daznis - gali buti nemazai nuliu, jei ju sekoje nerasite).
def find_codon_frequency(codon_list):
    standard_table = unambiguous_dna_by_name["Bacterial"]
    letters = standard_table.nucleotide_alphabet
    alphabet_list = []
    for c1 in letters:
        for c2 in letters:
            for c3 in letters:
                codon = c1 + c2 + c3
                alphabet_list.append(codon)
    codon_alphabet = dict.fromkeys(alphabet_list, 0.0)

    total = 0
    for orf in codon_list:
        codons_in_orf = [orf[i:i + 3] for i in range(0, len(orf), 3)]
        total += len(codons_in_orf)
        for codon in codons_in_orf:
            codon_alphabet[codon] += 1
    for key in codon_alphabet:
        codon_alphabet[key] = (codon_alphabet[key] / total) * 100
    return codon_alphabet


def find_dicodon_frequency(codon_list):
    standard_table = unambiguous_dna_by_name["Bacterial"]
    letters = standard_table.nucleotide_alphabet
    alphabet_list = []
    for c1 in letters:
        for c2 in letters:
            for c3 in letters:
                for c4 in letters:
                    for c5 in letters:
                        for c6 in letters:
                            dicodon = c1 + c2 + c3 + c4 + c5 + c6
                            alphabet_list.append(dicodon)
    dicodon_alphabet = dict.fromkeys(alphabet_list, 0.0)

    total = 0
    for orf in codon_list:
        dicodons_in_orf = [orf[i:i + 6] for i in range(0, len(orf), 6)]
        total += len(dicodons_in_orf)
        for dicodon in dicodons_in_orf:
            if len(dicodon) == 6:
                dicodon_alphabet[dicodon] += 1
    for key in dicodon_alphabet:
        dicodon_alphabet[key] = (dicodon_alphabet[key] / total) * 100
    return dicodon_alphabet


if __name__ == '__main__':
    record = SeqIO.read("data/bacterial4.fasta", "fasta")

    triplets = get_triplets(record.seq)

    # 1.
    codons = find_codons(triplets)
    print("1. All codons in a sequence:")
    codons = np.concatenate(codons)
    # print(codons)

    # 2.
    print("2. Each stop codons farthest start codon")
    farthest_start_codons = find_farthest_starts(triplets)
    # print(farthest_start_codons)

    # 3.
    shortest_codons = find_shortest_codons(codons)
    print("3. Codons that have a length less than 100 symbols:")
    # print(shortest_codons)

    print("4. Codon frequency (%):")
    print(find_codon_frequency(codons))

    print("4. Dicodon frequency (%):")
    print(find_dicodon_frequency(codons))
