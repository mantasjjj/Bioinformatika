import numpy as np
from Bio import SeqIO
from Bio.Data.CodonTable import unambiguous_dna_by_name


alphabet_list = []

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


# Is sekos kombinacijos pasiemam orfus, kurie prasideda ATG ir baigiasi TAA || TAG || TGA
def find_orf(seq_record):
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


def find_orfs(triplets):
    orfs = []
    for triplet in triplets:
        orfs.append(find_orf(triplet))
    return orfs


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
def find_longest_fragments(orfs):
    longest_orfs = []
    for orf in orfs:
        if len(orf) > 100:
            longest_orfs.append(orf)
    return longest_orfs


# 4. Parasykite funkcijas, kurios ivertintu kodonu ir dikodonu daznius(visi imanomi kodonai/dikodonai ir ju
# atitinkamas daznis - gali buti nemazai nuliu, jei ju sekoje nerasite).
def find_codon_frequency(list):
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
    for orf in list:
        codons_in_orf = [orf[i:i + 3] for i in range(0, len(orf), 3)]
        total += len(codons_in_orf)
        for codon in codons_in_orf:
            codon_alphabet[codon] += 1
    for key in codon_alphabet:
        codon_alphabet[key] = (codon_alphabet[key] / total) * 100
    return codon_alphabet


def find_dicodon_frequency(list):
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
    for orf in list:
        dicodons_in_orf = [orf[i:i + 6] for i in range(0, len(orf), 6)]
        total += len(dicodons_in_orf)
        for dicodon in dicodons_in_orf:
            if len(dicodon) == 6:
                dicodon_alphabet[dicodon] += 1
    for key in dicodon_alphabet:
        dicodon_alphabet[key] = (dicodon_alphabet[key] / total) * 100
    return dicodon_alphabet

# 5. Palyginkite kodonu bei dikodonu daznius tarp visu seku (atstumu matrica - kokia formule naudosite/kaip apskaiciuosite - parasykite ataskaitoje).
def print_distance_matrix(frequencies):
    ids = []

    for id in frequencies:
        ids.append(id)

    data = frequencies[ids[0]]

    print("*******", end = ' ')
    for i in ids:
        print(i, end = ' ')

    for i in range(0,8):
        print("")
        print(ids[i], end = ' ')
        for j in range(0, 8):
            total_points = 0
            for k in data:
                points = frequencies[ids[i]][k] - frequencies[ids[j]][k]
                if points < 0:
                    points = points * (-1)
                total_points += points
            print("%.2f" % total_points, end = ' ')


if __name__ == '__main__':
    files = ["bacterial1.fasta", "bacterial2.fasta", "bacterial3.fasta", "bacterial4.fasta", "mamalian1.fasta",
             "mamalian2.fasta", "mamalian3.fasta", "mamalian4.fasta"]

    codon_frequencies = {}
    dicodon_frequencies = {}

    for i in files:
        record = SeqIO.read("data/" + i, "fasta")
        triplets = get_triplets(record.seq)

        print("****************************" + record.id + "***********************************")
        # 1.
        codons = find_orfs(triplets)
        # print("1. All codons in a sequence:")
        codons = np.concatenate(codons)
        # print(codons)

        # 2.
        # print("2. Each stop codons farthest start codon")
        farthest_start_codons = find_farthest_starts(triplets)
        # print(farthest_start_codons)

        # 3.
        longest_fragments = find_longest_fragments(codons)
        # print("3. Fragments that have a length more than 100 symbols:")
        # print(longest_fragments)

        # 4.
        # print("4. Codon frequency (%):")
        codon_frequencies[record.id] = find_codon_frequency(longest_fragments)
        # print(find_codon_frequency(codons))

        # print("4. Dicodon frequency (%):")
        dicodon_frequencies[record.id] = find_dicodon_frequency(longest_fragments)
        # print(find_dicodon_frequency(codons))

    print("Codon frequency matrix:")
    print_distance_matrix(codon_frequencies)

    print("")

    print("Dicodon frequency matrix:")
    print_distance_matrix(dicodon_frequencies)