from collections import defaultdict
from typing import Dict, List, Tuple

def read_fasta(path: str) -> List[Tuple[str, str]]:
    """
    Lit un fichier FASTA et retourne une liste de (header, sequence).
    - header : ligne sans le '>'
    - sequence : chaîne de nucléotides sans espaces ni retours de ligne
    """
    sequences = []
    header = None
    seq_chunks = []

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue  # ignore les lignes vides

            if line.startswith(">"):
                # Si on avait déjà une séquence en cours, on la sauvegarde
                if header is not None:
                    seq = "".join(seq_chunks).upper()
                    sequences.append((header, seq))
                header = line[1:]  # on enlève le '>'
                seq_chunks = []
            else:
                # ligne de séquence
                seq_chunks.append(line)

    # Ne pas oublier la dernière séquence
    if header is not None:
        seq = "".join(seq_chunks).upper()
        sequences.append((header, seq))

    return sequences


def index_kmers(seq: str, k: int) -> Dict[str, List[int]]:
    """
    Construit un index des k-mers d'une séquence.
    Retourne un dict : kmer -> liste des positions de début (0-based).
    """
    index: Dict[str, List[int]] = defaultdict(list)
    n = len(seq)

    if k > n:
        return {} # pas de k-mers possibles plus long que la sequence

    for i in range(n - k + 1):
        kmer = seq[i:i + k]
        index[kmer].append(i)

    # On peut garder defaultdict ou le convertir en dict normal :
    return dict(index)

if __name__ == "__main__":
    # Lire la séquence inconnue (on suppose 1 seule séquence)
    unknown_list = read_fasta("tRNAs.fasta")
    header_query, seq_query = unknown_list[0]

    k = 11  # valeur pour  PLAST/BLAST
    kmer_index = index_kmers(seq_query, k)

    print(f"Nombre de k-mers distincts (k={k}) :", len(kmer_index))

    # Afficher quelques entrées
    count = 0
    for kmer, positions in kmer_index.items():
        print(f"k-mer {kmer} trouvé aux positions {positions}")
        count += 1
        if count >= 14:
            break