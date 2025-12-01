from collections import defaultdict
from typing import Dict, List, Tuple
from dataclasses import dataclass

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


@dataclass
class HSP:
    db_header: str
    q_start: int
    db_start: int
    length: int
    score: int
    q_end: int
    db_end: int
    aligned_query: str
    aligned_db: str


def find_initial_hsps(query_seq: str, db_sequences, k: int):
    """
    1.3 : Trouve tous les k-mers identiques entre la requête et la base.
    Retourne une liste de HSP initiaux (non étendus).
    """
    hsps = []
    index_q = index_kmers(query_seq, k)

    for db_header, db_seq in db_sequences:
        n = len(db_seq)
        if n < k:
            continue

        for db_pos in range(n - k + 1):
            kmer = db_seq[db_pos:db_pos + k]

            if kmer in index_q:
                for q_pos in index_q[kmer]:
                    # Score initial = k * 5 (matches seulement)
                    score = k * 5

                    hsps.append(
                        HSP(
                            db_header=db_header,
                            q_start=q_pos,
                            db_start=db_pos,
                            length=k,
                            score=score,
                            q_end=q_pos + k - 1,
                            db_end=db_pos + k - 1,
                            aligned_query=query_seq[q_pos:q_pos + k],
                            aligned_db=db_seq[db_pos:db_pos + k],
                        )
                    )

    return hsps

MATCH = 5
MISMATCH = -4

def extend_hsp_greedy(hsp: HSP, query_seq: str, db_seq: str, E: int):
    """
    1.4 : Extension gloutonne dans les deux directions.
    On choisit l'extension qui augmente le plus le score.
    On arrête si la perte dépasse le seuil -E.
    """
    qL = hsp.q_start
    qR = hsp.q_end
    dL = hsp.db_start
    dR = hsp.db_end

    best_score = hsp.score
    current_score = hsp.score
    max_seen = current_score

    # On boucle tant qu'on peut étendre
    while True:
        can_left = qL > 0 and dL > 0
        can_right = qR < len(query_seq) - 1 and dR < len(db_seq) - 1

        if not can_left and not can_right:
            break

        # Score extension gauche
        score_left = None
        if can_left:
            score_left = MATCH if query_seq[qL-1] == db_seq[dL-1] else MISMATCH

        # Score extension droite
        score_right = None
        if can_right:
            score_right = MATCH if query_seq[qR+1] == db_seq[dR+1] else MISMATCH

        # Choisir la meilleure extension
        direction = None
        if score_left is not None and (score_right is None or score_left >= score_right):
            direction = "L"
            delta = score_left
        else:
            direction = "R"
            delta = score_right

        # Tentative d'extension
        new_score = current_score + delta
        
        if max_seen - new_score > E:
            break

        # Appliquer l'extension
        current_score = new_score
        if direction == "L":
            qL -= 1
            dL -= 1
        else:
            qR += 1
            dR += 1

        # Mise à jour du meilleur
        if current_score > max_seen:
            max_seen = current_score

    # Mise à jour du HSP
    hsp.q_start = qL
    hsp.db_start = dL
    hsp.q_end = qR
    hsp.db_end = dR
    hsp.length = (qR - qL + 1)
    hsp.score = max_seen
    hsp.aligned_query = query_seq[qL:qR+1]
    hsp.aligned_db = db_seq[dL:dR+1]

    return hsp





if __name__ == "__main__":

    unknown_list = read_fasta("unknown.fasta")
    db_sequences = read_fasta("tRNAs.fasta")

    k = 11 
    E = 4   

    db_dict = dict(db_sequences)

    for header_query, seq_query in unknown_list:
        print("=" * 60)
        print(f"Requête : {header_query}")
        print(f"Longueur de la séquence requête: {len(seq_query)}")

        kmer_index = index_kmers(seq_query, k)
        print(f"Nombre de k-mers distincts (k={k}):", len(kmer_index))

        # (1.3)
        hsps = find_initial_hsps(seq_query, db_sequences, k=k)
        print(f"Nombre d'HSP initiaux trouvés: {len(hsps)}")

        # (1.4)
        extended = []
        for h in hsps:
            db_seq = db_dict[h.db_header]
            extended_h = extend_hsp_greedy(h, seq_query, db_seq, E)
            extended.append(extended_h)

        extended.sort(key=lambda h: h.score, reverse=True)

        print("\nMeilleurs HSP étendus:")
        for h in extended[:5]:
            print("-" * 40)
            print(">", h.db_header)
            print(f"score = {h.score}")
            print(f"query: [{h.q_start}:{h.q_end}]  db: [{h.db_start}:{h.db_end}]")
            print("Q:", h.aligned_query)
            print("D:", h.aligned_db)
        print()

