# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 15:26:03 2024

@author: Pauline
"""

import random
import itertools
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
import subprocess
import tempfile

# Demande les entrées de l'utilisateur
num_sequences = int(input("Nombre de séquences désiré: "))
seq_length = int(input("Longueur des séquences: "))
num_A = int(input("Nombre de A: "))
num_T = int(input("Nombre de T: "))
num_C = int(input("Nombre de C: "))
num_G = int(input("Nombre de G: "))
num_seq = int(input("Nombre de séquences aléatoires à générer (30 000 conseillées) : "))
site_5_prime = input("site de restriction en 5’: ")
site_3_prime = input("site de restriction en 3’: ")
min_Tm = float(input("Température de fusion minimale (Tm): "))
max_Tm = float(input("Température de fusion maximale (Tm): "))
mfe_threshold_RNAfold = float(input("Seuil de MFE pour RNAfold: "))
mfe_threshold_Zipfold = float(input("Seuil de MFE pour Zipfold: "))
mfe_threshold_RNAcofold = float(input("Seuil de MFE pour RNAcofold: "))


# Vérifier que la composition est correcte
assert num_A + num_T + num_C + num_G == seq_length, "La somme des nucléotides doit être égale à la longueur des séquences"

# Générer x séquences aléatoires
sequences = []
for k in range(num_seq):
    seq_list = ['A'] * num_A + ['T'] * num_T + ['C'] * num_C + ['G'] * num_G
    random.shuffle(seq_list)
    sequences.append(''.join(seq_list))

# Ajouter les sites de restriction définis par l’utilisateur
sequences = [site_5_prime + seq + site_3_prime for seq in sequences]

# 1ère sélection : Tm / NE MARCHE PAS / OBSCUR
def calculate_tm(sequence):
    return mt.Tm_NN(Seq(sequence))

sequences = [seq for seq in sequences if min_Tm <= calculate_tm(seq) <= max_Tm]
print(sequences)

# 2ème sélection : RNAfold
def calculate_mfe_rnafold(sequence):
    with tempfile.NamedTemporaryFile(delete=False) as temp_input:
        temp_input.write(sequence.encode())
        temp_input.flush()
        result = subprocess.run(['RNAfold', temp_input.name], capture_output=True, text=True)
    mfe = float(result.stdout.strip().split('\n')[1].split(' ')[-1].strip('()'))
    return mfe

sequences = [seq for seq in sequences if calculate_mfe_rnafold(seq) > mfe_threshold_RNAfold]

# 3ème sélection : Zipfold
def calculate_mfe_zipfold(sequence):
    with tempfile.NamedTemporaryFile(delete=False) as temp_input:
        temp_input.write(sequence.encode())
        temp_input.flush()
        result = subprocess.run(['Zipfold', temp_input.name], capture_output=True, text=True)
    mfe = float(result.stdout.strip().split('\n')[1].split(' ')[-1].strip('()'))
    return mfe

sequences = [seq for seq in sequences if calculate_mfe_zipfold(seq) > mfe_threshold_Zipfold]

# Calculer la MFE avec RNAcofold
def calculate_mfe_rnacofold(sequence1, sequence2):
    with tempfile.NamedTemporaryFile(delete=False) as temp_input:
        temp_input.write((sequence1 + '&' + sequence2).encode())
        temp_input.flush()
        result = subprocess.run(['RNAcofold', temp_input.name], capture_output=True, text=True)
    mfe = float(result.stdout.strip().split('\n')[1].split(' ')[-1].strip('()'))
    return mfe

# Filtrer les séquences par MFE
final_sequences = []
for i, seq1 in enumerate(sequences):
    for seq2 in sequences[i+1:]:
        mfe = calculate_mfe_rnacofold(seq1, seq2)
        if mfe > mfe_threshold_RNAcofold:
            final_sequences.append((seq1, seq2))
            if len(final_sequences) >= num_sequences:
                break
    if len(final_sequences) >= num_sequences:
        break

# Afficher les séquences finales
for seq_pair in final_sequences:
    print(seq_pair)




NIMPORTE QUOI
