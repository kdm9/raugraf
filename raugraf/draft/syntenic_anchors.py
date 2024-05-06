#!/usr/bin/env python3
from tqdm import tqdm
import natsort

import multiprocessing
from functools import partial
from collections import defaultdict
from itertools import islice
import os
from sys import stderr
import argparse
from pathlib import Path


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('-g', '--graph', required=True, type=Path,
            help = "Input pangenome graph as GFA file")
    ap.add_argument('-d', '--min-depth', type=int,
            help="Only output nodes with >= this depth")
    ap.add_argument('-D', '--max-depth', type=int,
            help="Only output nodes with <= this depth")
    ap.add_argument('-o','--out-tsv', type=Path,  default="/dev/stdout",
            help = "Output local complexity along each path as bedGraph file")

    args = ap.parse_args(argv)

    #G = defaultdict(set)
    nodelength = dict()
    paths = dict()
    node_paths = defaultdict(list)
    print("Load Graph from", args.graph, file=stderr)
    with open(args.graph) as fh:
        for line in tqdm(fh, desc="Parse GFA".ljust(20), unit="lines"):
            L = line.rstrip().split("\t")
            if L[0] == "S":
                nodelength[int(L[1])] = len(L[2])
            if L[0] == "L":
                #left = int(L[1])
                #right = int(L[3])
                #G[left].add(right)
                #G[right].add(left)
                pass
            if L[0] == "P":
                n = L[1]
                P = list(map(lambda x: int(x.rstrip("+-")), L[2].split(',')))
                for node in P:
                    node_paths[node].append(n)
                paths[n] = P
    n_nodes = len(nodelength)
    mean_node_len = sum(nodelength.values())/n_nodes
    print(f"Graph done, {n_nodes} nodes, {len(paths)} paths, mean node length {mean_node_len:0.1f} bp.", file=stderr)

    total_depth = {}
    unique_depth = {}
    for node in tqdm(node_paths, desc="Compute node depths".ljust(20), unit="nodes"):
        np = node_paths[node]
        td = len(np)
        if args.min_depth and td < args.min_depth:
            continue
        if args.max_depth and td > args.max_depth:
            continue
        ud = len(set(np))
        if args.min_depth and ud < args.min_depth:
            continue
        if args.max_depth and ud > args.max_depth:
            continue
        total_depth[node] = td
        unique_depth[node] = ud
    

    with open(args.out_tsv, "w") as ofh:
        print("#path", "left", "right", "node", "total_depth", "unique_depth", sep="\t", file=ofh)
        for path in natsort.natsorted(paths):
            left = 0
            for node in tqdm(paths[path], desc=f"Traverse {path}".ljust(20), unit="nodes"):
                nodel = nodelength[node]
                right = left + nodel
                if node in total_depth:
                    print(path, left, right, node, total_depth[node], unique_depth[node], sep="\t", file=ofh)
                left = right


if __name__ == '__main__':
    main()
