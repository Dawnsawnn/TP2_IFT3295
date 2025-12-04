from collections import defaultdict
from typing import Dict, List, Tuple
from dataclasses import dataclass
import math

# lit un fichier FASTA
def read_fasta(path: str) -> List[Tuple[str, str]]:

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

# construit l'index des k-mers
def index_kmers(seq: str, k: int) -> Dict[str, List[int]]:

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

LAMBDA = 0.192
K_CONST = 0.176


# 1.3 trouve les HSP initiaux
def find_initial_hsps(query, db, k):

    hsps = []
    index_q = index_kmers(query, k)

    for db_header, db_seq in db:
        n = len(db_seq)
        if n < k:
            continue

        for db_pos in range(n - k + 1):
            kmer = db_seq[db_pos:db_pos + k]

            if kmer in index_q:
                for q_pos in index_q[kmer]:

                    # calcule le score initial
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
                            aligned_query=query[q_pos:q_pos + k],
                            aligned_db=db_seq[db_pos:db_pos + k],
                        )
                    )

    return hsps

# 1.4 Extension greedy d'un HSP
def extend_hsp_greedy(hsp, query, db, E):
    qL = hsp.q_start
    qR = hsp.q_end
    dL = hsp.db_start
    dR = hsp.db_end

    score = hsp.score
    best = score

    while True:
        left_ok = (qL > 0 and dL > 0)
        right_ok = (qR < len(query)-1 and dR < len(db)-1)

        if not left_ok and not right_ok:
            break

        sL = MATCH if left_ok and query[qL-1] == db[dL-1] else (None if not left_ok else MISMATCH)
        sR = MATCH if right_ok and query[qR+1] == db[dR+1] else (None if not right_ok else MISMATCH)

        if sL is not None and (sR is None or sL >= sR):
            delta = sL
            new_qL, new_qR = qL - 1, qR
            new_dL, new_dR = dL - 1, dR
        else:
            delta = sR
            new_qL, new_qR = qL, qR + 1
            new_dL, new_dR = dL, dR + 1

        new_score = score + delta
        if best - new_score > E:
            break

        qL, qR, dL, dR = new_qL, new_qR, new_dL, new_dR
        score = new_score
        if score > best:
            best = score

    hsp.q_start = qL
    hsp.q_end = qR
    hsp.db_start = dL
    hsp.db_end = dR
    hsp.length = (qR - qL + 1)
    hsp.score = best
    hsp.aligned_query = query[qL:qR+1]
    hsp.aligned_db = db[dL:dR+1]

    return hsp



# 1.5 calcule le score d'un alignement
def compute_score(q: str, d: str) -> int:

    score = 0
    for a, b in zip(q, d):
        score += MATCH if a == b else MISMATCH
    return score

# fusionne les HSPs qui se chevauchent
def merge_hsps(hsps, query, db_all):
    par_seq = defaultdict(list)
    for h in hsps:
        par_seq[h.db_header].append(h)

    fusionnes = []
    db_dict = dict(db_all)

    for nom_seq, groupe in par_seq.items():
        groupe.sort(key=lambda x: x.q_start)
        en_cours = groupe[0]

        for h in groupe[1:]:
            if h.q_start <= en_cours.q_end and h.db_start <= en_cours.db_end:
                en_cours.q_end = max(en_cours.q_end, h.q_end)
                en_cours.db_end = max(en_cours.db_end, h.db_end)

                # recalcule l'alignement
                q = query[en_cours.q_start : en_cours.q_end+1]
                d = db_dict[nom_seq][en_cours.db_start : en_cours.db_end+1]
                en_cours.aligned_query = q
                en_cours.aligned_db = d
                en_cours.score = compute_score(q, d)
            else:
                fusionnes.append(en_cours)
                en_cours = h

        fusionnes.append(en_cours)

    return fusionnes

# 1.6

@dataclass
class HSPResult:
    hsp: HSP
    bitscore: int
    evalue: float


def bitscore(S: int) -> int:

    return int(round((LAMBDA * S - math.log(K_CONST)) / math.log(2)))


def evalue(B: int, m: int, n: int) -> float:

    return m * n * (2.0 ** (-B))


def taille_db(db: List[Tuple[str, str]]) -> int:

    return sum(len(seq) for _, seq in db)

# Filtre les HSP significatifs
def filtrer_hsps(hsps, db_sequences, taille_query, seuil):
    m = sum(len(seq) for _, seq in db_sequences)
    tmp = []

    for h in hsps:
        S = h.score
        B = int(round((0.192 * S - math.log(0.176)) / math.log(2)))
        e = m * taille_query * (2 ** (-B))
        if e <= seuil:
            tmp.append((h, B, e))

    best = {}
    for h, B, e in tmp:
        cle = h.db_header
        if cle not in best:
            best[cle] = (h, B, e)
        else:
            h2, B2, e2 = best[cle]
            if B > B2 or (B == B2 and e < e2):
                best[cle] = (h, B, e)

    return list(best.values())


# 1.7 : affichage des résultats

def print_plast(header_query, seq_query, db, merged, ss_threshold):
    print(40 * "-")
    print(f"Requête : {header_query}")
    print("Longueur de la séquence requête :", len(seq_query))
    print("Nombre d'HSP initiaux trouvés :", len(merged))

    good_hsps = filtrer_hsps(merged, db, len(seq_query), ss_threshold)

    print("\nRésultats significatifs :")
    if not good_hsps:
        print("Aucun HSP significatif.")
        return

    good_hsps.sort(key=lambda x: (-x[1], x[2]))

    count = 0
    for h, B, e in good_hsps:
        print(40 * "-")
        print(">" + h.db_header)
        print(f"# Best HSP score:{h.score:.2f}, bitscore:{B:.2f}, evalue:{e:.2e}")

        print(f"{h.q_start} {h.aligned_query} {h.q_end}")
        print(f"{h.db_start} {h.aligned_db} {h.db_end}")

        count += 1

    print("----------------------------------------")
    print("Total :", count)



if __name__ == "__main__":
    unknown_list = read_fasta("unknown.fasta")
    db = read_fasta("tRNAs.fasta")

    k = 11
    E = 4
    ss_threshold = 1e-3

    for header_query, seq_query in unknown_list:

        initial = find_initial_hsps(seq_query, db, k)

        extended = []
        db_dict = dict(db)
        for h in initial:
            extended.append(extend_hsp_greedy(h, seq_query, db_dict[h.db_header], E))

        merged = merge_hsps(extended, seq_query, db)

        print_plast(
            header_query,
            seq_query,
            db,
            merged,
            ss_threshold
        )


