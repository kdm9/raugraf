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


def noderadius(G, source_node, num_steps, nodes_seen=None):
    if nodes_seen is None:
        nodes_seen = set()
    if num_steps>1:
        try:
            neighbors = G[source_node]
        except KeyError:
            return nodes_seen
        for nbr in neighbors:
            if nbr not in nodes_seen:
                nodes_seen.add(nbr)
                noderadius(G, nbr, num_steps-1, nodes_seen)
    return nodes_seen


def batched(iterable, n=1):
    if n < 1:
        raise ValueError('n must be at least one')
    it = iter(iterable)
    while (batch := tuple(islice(it, n))):
        yield batch

def one_batch(nodes, G, jumps):
    return [(n, len(noderadius(G, n, jumps))) for n in nodes]

def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('-g', '--graph', required=True, type=Path,
            help = "Input pangenome graph as GFA file")
    ap.add_argument('-j', '--jumps', type=int, default=25,
            help= "How many nodes to traverse away from each seed node.")
    #ap.add_argument('-t', '--threads', type=int, default=1,
    #        help= "Num Parallel threads")
    ap.add_argument('-n','--node-table', type=Path,
            help = "Output local complexity per node (TSV)")
    ap.add_argument('-b','--node-bed', type=Path,
            help = "Output local complexity along each path as bedGraph file")

    args = ap.parse_args(argv)

    if args.node_table is None and args.node_bed is None:
        ap.error("Must give at least one of --node-table and/or --node-bed")

    G = defaultdict(set)
    nodelength = dict()
    paths = dict()
    print("Load Graph from", args.graph, file=stderr)
    with open(args.graph) as fh:
        for line in tqdm(fh, desc="Parse GFA".ljust(20), unit="lines"):
            L = line.rstrip().split("\t")
            if L[0] == "S":
                nodelength[int(L[1])] = len(L[2])
            if L[0] == "L":
                left = int(L[1])
                right = int(L[3])
                G[left].add(right)
                G[right].add(left)
            if L[0] == "P" and args.node_bed is not None:
                paths[L[1]] = list(map(lambda x: int(x.rstrip("+-")), L[2].split(',')))
    n_nodes = len(nodelength)
    mean_node_len = sum(nodelength.values())/n_nodes
    G = dict(G.items())
    print(f"Graph done, {n_nodes} nodes, {len(paths)} paths, mean node length {mean_node_len:0.1f} bp.", file=stderr)


    print("\nCalculating local complexity over", n_nodes, "nodes", file=stderr)
    results = {}
    for node in tqdm(G, desc="Compute complexity".ljust(20), unit="nodes"):
        results[node]=len(noderadius(G, node, args.jumps))
    #if args.threads == 1:
    #else:
    #    batches = batched(G, int(n_nodes/args.threads) + 1)
    #    with multiprocessing.Pool(processes=args.threads) as pool:
    #        with tqdm(total=n_nodes) as pbar:
    #            batch = partial(one_batch, G=G, jumps=args.jumps)
    #            for res in pool.imap(batch, batches):
    #                pbar.update(len(res))
    #                results.extend(res)
    print("\nComputed local complexity, outputing", file=stderr)

    if args.node_table:
        with open(args.node_table, "w") as ofh:
            print("node_id","node_radius", sep="\t", file=ofh)
            for node, ent in tqdm(sorted(results.items()), desc="To Node Table".ljust(20), unit="nodes"):
                print(node, ent, sep="\t", file=ofh)
    if args.node_bed:
        with open(args.node_bed, "w") as ofh:
            for path in tqdm(natsort.natsorted(paths), desc="Traverse Paths".ljust(20), total=len(paths), unit="paths"):
                left = 0
                for node in paths[path]:
                    nodel = nodelength[node]
                    right = left + nodel
                    print(path, left, right, results.get(node, 0), sep="\t", file=ofh)
                    left = right


if __name__ == '__main__':
    main()
