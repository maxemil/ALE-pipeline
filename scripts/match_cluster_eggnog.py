#! /usr/bin/env python3

import pandas as pd
import numpy as np
from Bio import SeqIO
import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument("-f", "--faas", required=True, nargs='+',
                    help="cluster .faas that are to be annotated")
parser.add_argument("-a", "--annotations", required=True, nargs='+',
                    help="annotated .faas that are used as a source for annotation")
parser.add_argument("-co", "--cluster_out", required=False, default='cluster_OG.tab',
                    help="output table containing all annotations of all clusters, default cluster_OG.tab")
parser.add_argument("-n", "--ancestor_nodes", required=True, nargs='+',
                    help="list of ancestor_nodes, ex '-n 83 82 81'")
parser.add_argument("-s", "--species_tree", required=True,
                    help="name of the species tree in 'events.txt'")
parser.add_argument("-th", "--threshold", required=False, default=0.5, type=float,
                    help="threshold for considering a gene cluster to be present, default 0.5")
parser.add_argument("-np", "--node_pairs", required=False, nargs='+',
                    help="threshold for considering a gene cluster to be present,ex '-np 83 82 82 81'  to get 82to81 and 83to82")
parser.add_argument("-ae", "--annotation_extension", required=False, default='.faa.emapper.annotations',
                    help="extension fo annotation files, default='.faa.emapper.annotations'")
parser.add_argument("-ce", "--cluster_extension", required=False, default='.bmge.aln.ufboot.clean.ale',
                    help="extension for cluster in events, default='.bmge.aln.ufboot.clean.ale'")
parser.add_argument("--oneline", required=False, default=False, action='store_true',
                    help="flag for formatting, one line per cluster and summary of present COGs")
args = parser.parse_args()


def parse_eggnog_annotations(annotations, annotation_extension):
    ann_dict_cat = {}
    ann_dict_cog = {}
    cog_description = {}
    cog_category = {}
    for ann in annotations:
        species = os.path.basename(ann).replace(annotation_extension, '')
        ann_dict_cat[species] = {}
        ann_dict_cog[species] = {}
        with open(ann) as handle:
            for line in handle:
                if not line.startswith('#'):
                    line = line.split('\t')
                    cog_description[line[9].split('|')[0]] = line[11].strip()
                    cog_category[line[9].split('|')[0]] = ";".join(line[10].split(', '))
                    ann_dict_cog[species][line[0]] = line[9].split('|')[0]
                    ann_dict_cat[species][line[0]] = line[10].split(', ')
    return ann_dict_cat, ann_dict_cog, cog_description, cog_category


def get_OGs_per_cluster(faas, ann_dict_cog):
    cluster_cogs = {}
    for f in faas:
        clst = os.path.basename(f).replace('.faa', '')
        cluster_cogs[clst] = {}
        for rec in SeqIO.parse(f,'fasta'):
            species = rec.id.split('_i_')[0]
            gene = rec.id.split('_i_')[1]
            try:
                cluster_cogs[clst][ann_dict_cog[species][gene]] += 1
            except:
                try:
                    cluster_cogs[clst][ann_dict_cog[species][gene]] = 1
                except:
                    pass
    return cluster_cogs


def get_OG_annotation_for_cluster(cluster_cogs, cog_description, cog_category, outfile):
    with open(outfile, 'w') as out:
        print("\t".join(["cluster", "NOG", '#NOG_annotations',
                         "NOG_category", "NOG_description"]), file=out)
        majority = 0
        maximum = 0
        single = 0
        none = 0
        for key, val in cluster_cogs.items():
            if not val:
                none += 1
                # print("no annotation available for %s" % key)
                print("\t".join([key, "-", "0", "-", "-"]), file=out)
            elif len(val) == 1:
                single += 1
                print("\t".join([key,
                            str(list(val.keys())[0]),
                            str(list(val.values())[0]),
                            cog_category[str(list(val.keys())[0])],
                            cog_category[str(list(val.keys())[0])]]), file=out)
                # print("undoubted COG %s for %s" % (val, key))
            else:
                for cog, count in val.items():
                    print("\t".join([key, cog, str(count), cog_category[cog], cog_description[cog]]), file=out)
                counts = 0
                majority_found = False
                for cog, count in val.items():
                    counts += count
                for cog, count in val.items():
                    if count/counts > .5:
                        # print("majority COG %s in %s" % (cog, key))
                        majority += 1
                        majority_found = True
                if not majority_found:
                    # print("no majority, but this is the most common one %s" %max(val.items(), key=operator.itemgetter(1))[0])
                    maximum += 1
        # print(np.sum([majority, single, none, maximum]))


def add_NOG_drop_cols(df, cluster_cogs, cluster_extension):
    df['clst'] = df.apply(lambda row: row['cluster'].replace(cluster_extension, ""), axis=1)
    del df['cluster']
    del df['species_tree']
    del df['node']

    df['NOG'] = df.apply(lambda row: cluster_cogs[row['clst']] ,axis=1)
    return df

def describe_ancestral_node(cluster_cogs, cog_description, cog_category, ancestor_node, species_tree, events_file, threshold, cluster_extension, oneline):
    def print_internal(nog, line, count, out, nog_cat, nog_desc):
        print("\t".join([line[1]['clst'],
                        nog,
                        str(count),
                        str(line[1]['duplications']),
                        str(line[1]['transfers']),
                        str(line[1]['losses']),
                        str(line[1]['originations']),
                        str(line[1]['copies']),
                        nog_cat,
                        nog_desc]), file=out)

    pd.options.mode.chained_assignment = None
    header = ['species_tree', 'cluster', 'node', 'duplications',
              'transfers', 'losses', 'originations', 'copies']

    df = pd.read_csv(events_file, sep='\t', names=header)
    df_subset = df[df['species_tree'] == species_tree]

    ancestor = df_subset[df_subset['node'] == ancestor_node]

    ancestor = add_NOG_drop_cols(ancestor, cluster_cogs, cluster_extension)

    # drop cluters that are not present in node
    ancestor = ancestor[ancestor.ix[:,4] > threshold]

    with open('%s_cluster_OG.tab' % ancestor_node, 'w') as out:
        print("\t".join(["cluster", "NOG", '#NOG_annotations', "duplications",
                         "transfers", "losses", "originations", "copies",
                         "NOG_category", "NOG_description"]), file=out)
        for line in ancestor.iterrows():
            if oneline:
                nogs = ";".join([nog for nog, count in line[1]['NOG'].items()])
                counts = ";".join([str(count) for nog, count in line[1]['NOG'].items()])
                cog_cats = ";".join(set([cog_category[nog] for nog, count in line[1]['NOG'].items()]))
                cog_descs = ";".join([cog_description[nog] for nog, count in line[1]['NOG'].items()])
                print_internal(nogs,
                               line,
                               counts,
                               out,
                               cog_cats,
                               cog_descs)
            else:
                for nog, count in line[1]['NOG'].items():
                    print_internal(nog,
                                   line,
                                   count,
                                   out,
                                   cog_category[nog],
                                   cog_description[nog])


def differences_nodes(cluster_cogs, cog_description, cog_category, ancestor_node1, ancestor_node2, species_tree, events_file, threshold, cluster_extension):
    pd.options.mode.chained_assignment = None
    header = ['species_tree', 'cluster', 'node', 'duplications',
              'transfers', 'losses', 'originations', 'copies']

    df = pd.read_csv(events_file, sep='\t', names=header)
    df_subset = df[df['species_tree'] == species_tree]

    ancestor1 = df_subset[df_subset['node'] == ancestor_node1]
    ancestor2 = df_subset[df_subset['node'] == ancestor_node2]

    ancestor1 = add_NOG_drop_cols(ancestor1, cluster_cogs, cluster_extension)
    ancestor2 = add_NOG_drop_cols(ancestor2, cluster_cogs, cluster_extension)

    gains_df = pd.DataFrame(columns=ancestor1.columns)
    losses_df = pd.DataFrame(columns=ancestor2.columns)

    for clst in ancestor1['clst'].tolist():
        if (ancestor1[ancestor1['clst'] == clst]['copies'] < 0.5).tolist()[0] and \
        (ancestor2[ancestor2['clst'] == clst]['copies'] >= 0.5).tolist()[0]:
            gains_df = gains_df.append(ancestor2[ancestor2['clst'] == clst])
        elif (ancestor1[ancestor1['clst'] == clst]['copies'] >= 0.5).tolist()[0] and \
        (ancestor2[ancestor2['clst'] == clst]['copies'] < 0.5).tolist()[0]:
            losses_df = losses_df.append(ancestor2[ancestor2['clst'] == clst])

    with open("gains_%sto%s.tab" % (ancestor_node1, ancestor_node2), 'w') as gains,  \
        open("losses_%sto%s.tab" % (ancestor_node1, ancestor_node2), 'w') as losses:
        print("\t".join(["cluster", "type of gain", "NOG cat", "NOGs"]), file=gains)
        print("\t".join(["cluster", "NOG cat", "NOGs"]), file=losses)
        for line in gains_df.iterrows():
            gain_type = "-"
            if line[1]['transfers'] >= 0.5:
                gain_type = "transfer"
            elif line[1]['originations'] >= 0.5:
                gain_type = "origination"
            elif line[1]['transfers'] >= 0.5 and line[1]['originations'] >= 0.5:
                gain_type = "transfer/origination"
            print("\t".join([line[1]['clst'],
                                gain_type,
                                ";".join(set([cog_category[nog] for nog in line[1]['NOG'].keys()])),
                                ";".join(line[1]['NOG'].keys())
                                ]), file=gains)
        for line in losses_df.iterrows():
            print("\t".join([line[1]['clst'],
                                ";".join(set([cog_category[nog] for nog in line[1]['NOG'].keys()])),
                                ";".join(line[1]['NOG'].keys())
                                ]), file=losses)



def main(args):
    ann_dict_cat, ann_dict_cog, cog_description, cog_category = parse_eggnog_annotations(args.annotations, args.annotation_extension)
    cluster_cogs = get_OGs_per_cluster(args.faas, ann_dict_cog)
    get_OG_annotation_for_cluster(cluster_cogs, cog_description, cog_category, args.cluster_out)

    # species_tree = "Alphaproteobacteria_species_recoded_clean"
    events_file = "events.txt"
    for ancestor_node in args.ancestor_nodes:
        describe_ancestral_node(cluster_cogs, cog_description, cog_category, ancestor_node, args.species_tree, events_file, args.threshold, args.cluster_extension, args.oneline)
    if args.node_pairs:
        for i in range(0,len(args.node_pairs),2):
            differences_nodes(cluster_cogs, cog_description, cog_category, args.node_pairs[i], args.node_pairs[i + 1], args.species_tree, events_file, args.threshold, args.cluster_extension)


if __name__ == '__main__':
    main(args)


# for key1, set1 in cluster_cogs.items():
#     if set1:
#         for key2, set2 in cluster_cogs.items():
#             if not key1 == key2:
#                 if set1 == set2:
#                     print(set1)


# faas = !ls cluster_faas/
# cluster_cogs = {}
# cluster_cats = {}
# for f in faas:
#     cluster_cogs[f.replace('.faa', '')] = set()
#     cluster_cats[f.replace('.faa', '')] = set()
#     for rec in SeqIO.parse("cluster_faas/%s" % f,'fasta'):
#         try:
#             cluster_cogs[f.replace('.faa', '')].add(ann_dict_cog[rec.id.split('_i_')[0]][rec.id.split('_i_')[1]])
#         except:
#             print("did not find a cog annotation for %s" % rec.id)
#         try:
#             cluster_cats[f.replace('.faa', '')].add(ann_dict_cat[rec.id.split('_i_')[0]][rec.id.split('_i_')[1]])
#         except:
#             print("did not find a cat annotation for %s" % rec.id)
