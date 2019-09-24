#! /usr/bin/env python3

import numpy as np
import argparse
from ALE_parsing import *
import sys
import glob

parser = argparse.ArgumentParser()

parser.add_argument("-f", "--faas", required=True,
                    help="cluster .faas that are to be annotated")
parser.add_argument("-a", "--eggnog_annotation", required=True, type=str,
                    help="eggnog annotations files similar to cognames")
parser.add_argument("-c", "--cog_annotation", required=True, type=str,
                    help="cognames, the official cog annotations")
parser.add_argument("-co", "--cluster_out", required=False, default='cluster_OG.tab',
                    help="output table containing all annotations of all clusters, default 'cluster_OG.tab'")
parser.add_argument("-n", "--ancestor_nodes", required=False, nargs='+', default=[],
                    help="list of ancestor_nodes, ex '-n 83 82 81'")
parser.add_argument("-s", "--species_tree", required=True,
                    help="name of the species tree in 'events.txt'")
parser.add_argument("-e", "--events", required=False, default="events.txt",
                    help="name of the events file, default 'events.txt'")
parser.add_argument("-th", "--threshold", required=False, default=0.5, type=float,
                    help="threshold for considering a gene cluster to be present, default 0.5")
parser.add_argument("-np", "--node_pairs", required=False, nargs='+', default=[],
                    help="list of ancestor node pairs, ex '-np 83 82 82 81'  to get 82to81 and 83to82")
args = parser.parse_args()


def main(args):
    if not args.ancestor_nodes and not args.node_pairs:
        print('either -n or -np has to be given')
        sys.exit()

    annotation = Annotation(args.cog_annotation, args.eggnog_annotation)

    name2cluster = {}
    for f in glob.glob(args.faas + "/*"):
        try:
            c = Cluster(f, annotation)
            name2cluster[c.name] = c
        except KeyError:
            print("could not find ALE results for cluster {}".format(f), file=sys.stderr)

    with open(args.cluster_out, 'w') as outhandle:
        for c in name2cluster.values():
            print(c, file=outhandle)

    events = parse_events(args.events, args.species_tree)

    name2nodes = {}
    for node in set(args.ancestor_nodes + args.node_pairs):
        n = Node(node, annotation, events, name2cluster, args.threshold)
        name2nodes[n.name] = n
        with open("{}_cluster_OG.tab".format(n.name), 'w') as outhandle:
            print(n, file=outhandle)

    if args.node_pairs:
        for i in range(0, len(args.node_pairs), 2):
            gains, losses, transfers = name2nodes[args.node_pairs[i + 1]].gains_losses(name2nodes[args.node_pairs[i + 1]])
            with open("gains_{}to{}.tab".format(args.node_pairs[i], args.node_pairs[i + 1]), 'w') as gainsfile:
                print(gains, file=gainsfile)
            with open("losses_{}to{}.tab".format(args.node_pairs[i], args.node_pairs[i + 1]), 'w') as lossesfile:
                print(losses, file=lossesfile)
            with open("transfers_{}to{}.tab".format(args.node_pairs[i], args.node_pairs[i + 1]), 'w') as transfersfile:
                print(transfers, file=transfersfile)


if __name__ == '__main__':
    main(args)
