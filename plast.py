from collections import defaultdict
from typing import Dict, List, Tuple
from dataclasses import dataclass
import math


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
        return {}  # pas de k-mers possibles plus long que la sequence

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


MATCH = 5
MISMATCH = -4

# constantes pour 1.6
LAMBDA = 0.192
K_CONST = 0.176


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
                    score = k * MATCH

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
            score_left = MATCH if query_seq[qL - 1] == db_seq[dL - 1] else MISMATCH

        # Score extension droite
        score_right = None
        if can_right:
            score_right = MATCH if query_seq[qR + 1] == db_seq[dR + 1] else MISMATCH

        # Choisir la meilleure extension
        if score_left is not None and (score_right is None or score_left >= score_right):
            direction = "L"
            delta = score_left
        else:
            direction = "R"
            delta = score_right

        # Tentative d'extension
        new_score = current_score + delta

        # Vérification du seuil -E
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
    hsp.aligned_query = query_seq[qL:qR + 1]
    hsp.aligned_db = db_seq[dL:dR + 1]

    return hsp


# 1.5

def compute_score(q: str, d: str) -> int:
    """
    Recalcule le score brut d'un alignement sans gap.
    +5 pour un match, -4 pour un mismatch.
    """
    score = 0
    for a, b in zip(q, d):
        score += MATCH if a == b else MISMATCH
    return score


def merge_hsp(
    hsps: List[HSP],
    query_seq: str,
    db_sequences: List[Tuple[str, str]]
) -> List[HSP]:
    """
    1.5 : Fusionne les HSPs qui se chevauchent sur la même séquence de la base,
    et sur la même diagonale (q_start - db_start constant).
    On suppose qu'il n'y a pas de gaps.
    """
    if not hsps:
        return []

    db_dict = dict(db_sequences)

    # grouper par (séquence de la base, diagonale)
    groups: Dict[Tuple[str, int], List[HSP]] = defaultdict(list)
    for h in hsps:
        diag = h.q_start - h.db_start
        groups[(h.db_header, diag)].append(h)

    merged: List[HSP] = []

    for (db_header, diag), group in groups.items():
        group.sort(key=lambda h: h.q_start)
        db_seq = db_dict[db_header]

        cur_q_start = group[0].q_start
        cur_q_end = group[0].q_end

        for h in group[1:]:
            if h.q_start <= cur_q_end:
                # chevauchement ou contiguïté
                cur_q_end = max(cur_q_end, h.q_end)
            else:
                # ferme le bloc courant
                q_start = cur_q_start
                q_end = cur_q_end
                db_start = q_start - diag
                db_end = db_start + (q_end - q_start)

                aligned_q = query_seq[q_start:q_end + 1]
                aligned_d = db_seq[db_start:db_end + 1]
                score = compute_score(aligned_q, aligned_d)

                merged.append(HSP(
                    db_header = db_header,
                    q_start = q_start,
                    db_start = db_start,
                    length = q_end - q_start + 1,
                    score = score,
                    q_end = q_end,
                    db_end = db_end,
                    aligned_query = aligned_q,
                    aligned_db = aligned_d,
                ))

                # nouveau bloc
                cur_q_start = h.q_start
                cur_q_end = h.q_end

        # dernier bloc du groupe
        q_start = cur_q_start
        q_end = cur_q_end
        db_start = q_start - diag
        db_end = db_start + (q_end - q_start)

        aligned_q = query_seq[q_start:q_end + 1]
        aligned_d = db_seq[db_start:db_end + 1]
        score = compute_score(aligned_q, aligned_d)

        merged.append(HSP(
            db_header = db_header,
            q_start = q_start,
            db_start = db_start,
            length = q_end - q_start + 1,
            score = score,
            q_end = q_end,
            db_end = db_end,
            aligned_query = aligned_q,
            aligned_db = aligned_d,
        ))

    return merged


# 1.6

@dataclass
class HSPWithStats:
    hsp: HSP
    bitscore: int
    evalue: float


def bitscore_from_raw(S: int) -> int:

    return int(round((LAMBDA * S - math.log(K_CONST)) / math.log(2)))


def evalue_from_bitscore(B: int, m: int, n: int) -> float:

    return m * n * (2.0 ** (-B))


def total_db_length(db_sequences: List[Tuple[str, str]]) -> int:
    """
    Taille totale m des séquences de la banque.
    """
    return sum(len(seq) for _, seq in db_sequences)


def filter_significant_hsps(
    merged_hsps: List[HSP],
    db_sequences: List[Tuple[str, str]],
    query_length: int,
    ss_threshold: float
) -> List[HSPWithStats]:
    """
    Ajoute bitscore et e-value à chaque HSP, filtre selon le seuil -ss
    et garde un seul HSP (le meilleur) par séquence de la banque.
    """
    if not merged_hsps:
        return []

    m = total_db_length(db_sequences)
    results: List[HSPWithStats] = []

    for h in merged_hsps:
        S = h.score
        B = bitscore_from_raw(S)
        e = evalue_from_bitscore(B, m, query_length)
        results.append(HSPWithStats(hsp=h, bitscore=B, evalue=e))

    # filtrage selon le seuil de significativité
    filtered = [x for x in results if x.evalue <= ss_threshold]

    # on garde le meilleur par séquence de la banque (db_header)
    best_per_db: Dict[str, HSPWithStats] = {}
    for hs in filtered:
        key = hs.hsp.db_header
        if key not in best_per_db:
            best_per_db[key] = hs
        else:
            cur = best_per_db[key]
            if hs.bitscore > cur.bitscore or (
                hs.bitscore == cur.bitscore and hs.evalue < cur.evalue
            ):
                best_per_db[key] = hs

    return list(best_per_db.values())


# 1.7 : affichage des résultats

def print_plast(
    header_query: str,
    seq_query: str,
    db_sequences: List[Tuple[str, str]],
    k: int,
    E: int,
    ss_threshold: float,
):
    """
    Exécute les étapes 1.3–1.7 pour une requête donnée,
    et affiche l'output dans un format proche de l'exemple du TP.
    """
    db_dict = dict(db_sequences)

    print("=" * 60)
    print(f"Requête : {header_query}")
    print(f"Longueur de la séquence requête : {len(seq_query)}")

    # 1.3
    hsps = find_initial_hsps(seq_query, db_sequences, k=k)
    print(f"Nombre d'HSP initiaux trouvés : {len(hsps)}")

    # 1.4
    extended: List[HSP] = []
    for h in hsps:
        db_seq = db_dict[h.db_header]
        extended_h = extend_hsp_greedy(h, seq_query, db_seq, E)
        extended.append(extended_h)

    # 1.5
    merged = merge_hsp(extended, seq_query, db_sequences)

    # 1.6
    significant = filter_significant_hsps(
        merged,
        db_sequences=db_sequences,
        query_length=len(seq_query),
        ss_threshold=ss_threshold,
    )

    # tri par pertinence
    significant.sort(key=lambda x: (-x.bitscore, x.evalue))

    print("\nRésultats significatifs :")

    for hs in significant:
        h = hs.hsp
        print(f">{h.db_header}")
        # Best HSP
        print(
            f"Best HSP score:{h.score:.2f}, "
            f"bitscore:{hs.bitscore:.2f}, "
            f"evalue: {hs.evalue:.2e}"
        )
        # alignement
        print(f"{h.q_start + 1} {h.aligned_query} {h.q_end + 1}")
        print(f"{h.db_start + 1} {h.aligned_db} {h.db_end + 1}")
        print("-" * 40)

    print(f"Total : {len(significant)}\n")


if __name__ == "__main__":
    unknown_list = read_fasta("unknown.fasta")
    db_sequences = read_fasta("tRNAs.fasta")

    k = 11          # taille de la graine
    E = 4           # seuil d'extension
    ss_threshold = 1e-3  # seuil de significativité

    for header_query, seq_query in unknown_list:
        print_plast(
            header_query,
            seq_query,
            db_sequences = db_sequences,
            k = k,
            E = E,
            ss_threshold = ss_threshold,
        )
