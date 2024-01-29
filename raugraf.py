from Bio import SeqIO
from tqdm import tqdm

from sys import argv, stdin, stdout, stderr
from collections import Counter, defaultdict
from concurrent.futures import as_completed, ProcessPoolExecutor
from itertools import islice
import argparse
import os

SEB={"A": 0,
     "S": 1, "T": 1,
     "R": 2, "K": 2,
     "H": 3,
     "N": 4, "D": 4,
     "E": 5, "Q": 5,
     "C": 6,
     "P": 7,
     "G": 8,
     "I": 10, "V": 10,
     "L": 11, "M": 11,
     "F": 12, "Y": 12,
     "W": 13,
     "*": 14,
     "X": 15}
NT = {"A": 0, "C": 1, "G": 2, "T":3}



def B(num, bits=10*3):
    return f"{bin(num)[2:]:>0{bits}}"

def mers(seq, k=25, truncate=False, hashed=False, aa=False):
    """Return kmers of SEB-encoded AA sequences (see miniprot paper).

    k = kmer length
    truncate = mask LSD? the coding of AA -> number is done such that similar
               AA's have similar MSB, so masking the LSB gives a more sensitive
               hash.
    """
    bitper = (3 if truncate else 4) if aa else 2
    mask = (1<<(bitper*k)) -1
    h = 0
    l = 0
    i = 0
    while i < (len(seq)-k+1):
        a=seq[i]
        if aa:
            if truncate:
                h = mask & ( (h << bitper) | ((SEB.get(a, 15) & 0xe) >> 1))
            else:
                h = mask & ( (h << bitper) | (SEB.get(a, 15)) ) 
        else:
            try:
                h = mask & ( (h << bitper) | (NT[a]) ) 
            except KeyError:
                l = 0
        i += 1
        if l < k:
            l += 1
        if l == k:
            yield i, hash64(h, mask) if hashed else h

def is_syncmer(n, k=10, s=2, truncate=False, aa=False):
    bitper = (3 if truncate else 4) if aa else 2
    kbit = k*bitper
    sbit = s*bitper
    mask = (1 << (sbit))-1
    mins = mask
    mini = 0
    for i in range(0, kbit-sbit + 1, bitper):
        # scan from right to left for the lowest window of sbit bits
        smer = ((n >> i) & mask)
        if smer < mins:
            mins = smer
            mini = i
    return mini == 0 or mini == (kbit - sbit)

def hash64(key, mask):
    key = (~key + (key << 21)) & mask
    key = key ^ key >> 24
    key = ((key + (key << 3)) + (key << 8)) & mask
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask
    return key


def syncmers(seq, k=9, s=2, truncate=False, hashed=False, aa=False):
    """Return k, s syncmers of SEB-encoded AA sequences (see miniprot paper).

    k = kmer length, aa
    s = submer length, aa
    truncate = mask LSD? the coding of AA -> number is done such that similar
               AA's have similar MSB, so masking the LSB gives a more sensitive
               hash.
    hashed = hash64() kmers to randomise order, or keep lexographic?
    """
    for i, mer in mers(seq, k=k, truncate=truncate, hashed=False, aa=aa):
        if is_syncmer(mer, k=k, s=s, truncate=truncate, aa=aa):
                yield (i, mer)



def noderadius(G, source_node, num_steps, nodes_seen=None):
    if nodes_seen is None:
        nodes_seen = set()
    if num_steps>1:
        neighbors = G[source_node]
        for nbr in neighbors:
            if nbr not in nodes_seen:
                nodes_seen.add(nbr)
                noderadius(G, nbr, num_steps-1, nodes_seen)
    return nodes_seen


def load_seq(seq, k=24, s=4, hashed=False):
    db = defaultdict(set)
    last = None
    for i, mer in syncmers(seq.seq, truncate=False, k=k, s=s, hashed=hashed):
        if last:
            db[last].add(mer)
            db[mer].add(last)
        last = mer
    return db


def batched(iterable, n=1):
    if n < 1:
        raise ValueError('n must be at least one')
    it = iter(iterable)
    while (batch := tuple(islice(it, n))):
        yield batch

def batch_node_radii(G, nodes, num_steps):
    return {n: len(noderadius(G, n, num_steps)) for n in nodes}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", "-t", default=os.cpu_count(), type=int,
            help="Number of compute threads")
    ap.add_argument("--kmer", "-k", default=25, type=int,
            help="Kmer size")
    ap.add_argument("--syncmer", "-s", default=4, type=int,
            help="Syncmers: take only kmers with the N minimal bases at one end")
    ap.add_argument("--steps", "-r", default=25, type=int,
            help="Node radius to traverse when calculating local node complexity")
    ap.add_argument("--output", "-o", default=stdout, type=argparse.FileType("wt"),
            help="Output bedfile")
    ap.add_argument("fasta",
            help="Input fasta file")
    args = ap.parse_args(argv)

    db = defaultdict(set)
    n = 0
    with ProcessPoolExecutor(args.threads) as exc:
        jobs = []
        for seq in SeqIO.parse(args.fasta, "fasta"):
            n+=1
            jobs.append(exc.submit(load_seq, seq, k=args.kmer, s=args.syncmer, hashed=False))
        for job in tqdm(as_completed(jobs), desc="Count chroms"):
            for mer, nbrs in job.result().items():
                db[mer].update(nbrs)
    print(f"DB is {len(db)} syncmers from {n} sequences", file=stderr)

    noderad = {}
    with ProcessPoolExecutor(args.threads) as exc:
        jobs = []
        for nodebatch in batched(db, int(len(db)/args.threads)+1):
            jobs.append(exc.submit(batch_node_radii, db, nodebatch, args.steps))
        for job in tqdm(as_completed(jobs), desc="Compute node radii", total=len(jobs)):
            noderad.update(job.result())

    for seq in tqdm(SeqIO.parse(args.fasta, "fasta"), desc="Iterate graph"):
        left = 0
        for i, mer in syncmers(seq.seq, truncate=False, k=args.kmer, s=args.syncmer, hashed=False):
            right = i+args.kmer
            print(seq.id, left, right, noderad[mer], sep="\t", file=args.output)
            left=right

if __name__ == "__main__":
    main()
