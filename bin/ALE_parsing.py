import ete3
from Bio import SeqIO
import pandas as pd
import os
from collections import defaultdict, OrderedDict
import sys


def to_tsv(*args):
    return("\t".join(args))


def to_ssv(*args):
    return("; ".join(args))


class Annotation:
    def __init__(self, cog_annotation, eggnog_annotation):
        self.cog2cat = defaultdict(lambda: "-")
        self.cog2desc = defaultdict(lambda: "-")
        self.parse_cog_annotation(cog_annotation)
        self.parse_eggnog_annotations(eggnog_annotation)

    def parse_cog_annotation(self, cog_file):
        with open(cog_file, 'rU', encoding='windows-1252') as cogfile:
            for line in cogfile:
                if not line.startswith('#'):
                    line = line.strip().split('\t')
                    self.cog2cat[line[0]] = line[1]
                    self.cog2desc[line[0]] = line[2]


    def parse_eggnog_annotations(self, eggnog_file):
        with open(eggnog_file) as nogfile:
            for line in nogfile:
                if not line.startswith('#'):
                    line = line.strip().split('\t')
                    if not line[0] in self.cog2cat.keys():
                        self.cog2desc[line[0]] = line[2]
                        self.cog2cat[line[0]] = line[1]

class Cluster:
    def __init__(self, faa, annotation):
        self.annotation = annotation
        self.faa = faa
        self.name = os.path.basename(faa).split('.')[0]
        self.cog_members = defaultdict(int)
        self.get_cog_members()
        self.cog_members = OrderedDict(sorted(self.cog_members.items(), key=lambda t: t[1], reverse=True))

    def get_cog_members(self):
        for rec in SeqIO.parse(self.faa, 'fasta'):
            gene = rec.id#.split('_i_')[1]
            try:
                self.cog_members[self.annotation.id2cog[gene]] += 1
            except:
                # print("{} not annotated".format(gene), file=sys.stderr)
                pass

    def __str__(self):
        cog_lines = []
        for cog, count in self.cog_members.items():
            cog_lines.append(to_tsv(self.name, cog, str(count), self.annotation.cog2cat[cog], self.annotation.cog2desc[cog]))
        return("\n".join(cog_lines))


def parse_events(events_file, species_tree):
    pd.options.mode.chained_assignment = None
    header = ['species_tree', 'cluster', 'node', 'duplications',
              'transfers', 'losses', 'originations', 'copies']
    all_events = pd.read_csv(events_file, sep='\t', names=header)
    events = all_events[all_events['species_tree'] == species_tree]
    del events['species_tree']
    events['cluster'] = events['cluster'].str.replace('.clean.ufboot.ale', '')
    events.index = pd.MultiIndex.from_arrays([list(events['cluster']), list(events['node'])], names=["cluster", "node"])
    del events['cluster']
    del events['node']
    return events


class Node:
    def __init__(self, name, annotation, events, name2cluster, threshold):
        self.name = name
        self.annotation = annotation
        self.threshold = threshold
        self.events = events
        self.name2cluster = name2cluster
        self.raw_members = events.xs(self.name, level='node')
        self.members = self.raw_members#[self.raw_members["copies"] > self.threshold]
        self.header = ["cluster", "NOG_category", "NOG_description", "duplications",
                                 "transfers", "losses", "originations", "copies"]

    def get_cog_info(self, clst):
        cog_cat = self.annotation.cog2cat[clst]
        cog_desc = self.annotation.cog2desc[clst]
        return (cog_cat, cog_desc)

    def __str__(self):
        lines = []
        lines.append(to_tsv(*self.header))
        for line in self.members.iterrows():
            clst = line[0]
            cog_cat, cog_desc = self.get_cog_info(clst)
            events = [str(i) for i in line[1]]
            clst_line = to_tsv(
                        clst,
                        cog_cat,
                        cog_desc,
                        *events)
            lines.append(clst_line)
        return "\n".join(lines)

    def gains_losses(self, descendant):
        gains_df = pd.DataFrame(columns=self.raw_members.columns)
        losses_df = pd.DataFrame(columns=self.raw_members.columns)
        transfers_df = pd.DataFrame(columns=self.raw_members.columns)
        for line in self.raw_members.iterrows():
            clst = line[0]
            if line[1].loc[['originations']].sum() >= self.threshold:
                gains_df = gains_df.append(self.raw_members.loc[clst])
            elif line[1].loc['losses'] > self.threshold:
                losses_df = losses_df.append(descendant.raw_members.loc[clst])
            if line[1].loc[['transfers']].sum() >= self.threshold:
                transfers_df = transfers_df.append(self.raw_members.loc[clst])
        return self.print_gains(gains_df), self.print_losses(losses_df), self.print_transfers(transfers_df)

    def print_gains(self, gains_df):
        lines = []
        lines.append(to_tsv("cluster", "acquisition_freq", "NOG_category", "NOG_description"))
        for line in gains_df.iterrows():
            clst = line[0]
            cog_cat, cog_desc = self.get_cog_info(clst)
            lines.append(to_tsv(clst,
                                str(line[1].loc['originations']),
                                cog_cat,
                                cog_desc))
        return "\n".join(lines)

    def print_losses(self, losses_df):
        lines = []
        lines.append(to_tsv("cluster", "loss_freq", "NOG_category", "NOG_description"))
        for line in losses_df.iterrows():
            clst = line[0]
            cog_cat, cog_desc = self.get_cog_info(clst)
            lines.append(to_tsv(clst,
                                str(line[1].loc['losses']),
                                cog_cat,
                                cog_desc))
        return "\n".join(lines)


    def print_transfers(self, transfers_df):
        lines = []
        lines.append(to_tsv("cluster", "transfer_freq", "NOG_category", "NOG_description"))
        for line in transfers_df.iterrows():
            clst = line[0]
            cog_cat, cog_desc = self.get_cog_info(clst)
            lines.append(to_tsv(clst,
                                str(line[1].loc['transfers']),
                                cog_cat,
                                cog_desc))
        return "\n".join(lines)

def check_origination(n, clst, events):
    if not (n.is_leaf() or n.is_root()):
    #print(n.name, events.loc['0ZTPB'].loc[[n.up.name, n.name, n.children[0].name, n.children[1].name]].sum())
        try:
            if events.loc[clst].loc[[n.up.name, n.name, n.children[0].name, n.children[1].name], 'originations'].sum() > 0.7:
                return(clst)
        except:
            pass



# pd.options.mode.chained_assignment = None
# header = ['species_tree', 'cluster', 'node', 'duplications','transfers', 'losses', 'originations', 'copies']
# all_events = pd.read_csv("MGIV-eurySel-corrected/events.txt", sep='\t', names=header)
# events = all_events[all_events['species_tree'] == "eurySel.noSA1.chi2prune_clean"]
# del events['species_tree']
# # events['cluster'] = events.apply(lambda row: row['cluster'].split('.')[0], axis=1)
# events['cluster'] = events['cluster'].str.replace('.clean.ufboot.ale', '')
#
# events.index = pd.MultiIndex.from_arrays([list(events['cluster']), list(events['node'])], names=["cluster", "node"])
# del events['cluster']
# del events['node']
#
#
# tree = ete3.PhyloTree('species_trees/eurySel.noSA1.chi2prune_clean.tree', format=1)
#
# #for n in tree.traverse():
# n = tree & "81"
# for clst in set(events.index.get_level_values(0)):
#      res = check_origination(n, clst, events)
#      if res:
#          print(res)
