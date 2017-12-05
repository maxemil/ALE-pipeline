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
    def __init__(self, annotation_files, cog_annotation):
        self.cog_annotation = {}
        self.annotation_files = annotation_files
        self.id2cat = defaultdict(str)
        self.id2cog = defaultdict(str)
        self.cog2cat = defaultdict(str)
        self.cog2desc = defaultdict(str)
        self.parse_cog_annotation(cog_annotation)
        self.parse_eggnog_annotations()

    def parse_cog_annotation(self, cog_file):
        with open(cog_file, 'rU', encoding='windows-1252') as cogfile:
            for line in cogfile:
                if not line.startswith('#'):
                    line = line.strip().split('\t')
                    self.cog_annotation[line[0]] = line[2]

    def parse_eggnog_annotations(self):
        for ann in self.annotation_files:
            with open(ann) as handle:
                for line in handle:
                    if not line.startswith('#'):
                        line = line.split('\t')
                        cog = line[9].split('|')[0]
                        try:
                            desc = self.cog_annotation[cog]
                        except:
                            desc = line[11].strip()
                        seq_id = line[0]
                        cat = line[10].split(', ')
                        self.cog2desc[cog] = desc
                        self.cog2cat[cog] = ";".join(cat)
                        self.id2cog[seq_id] = cog
                        self.id2cat[seq_id] = cat

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


class Events:
    def __init__(self, events_file, species_tree):
        pd.options.mode.chained_assignment = None
        self.header = ['species_tree', 'cluster', 'node', 'duplications',
                  'transfers', 'losses', 'originations', 'copies']
        all_events = pd.read_csv(events_file, sep='\t', names=self.header)
        self.events = all_events[all_events['species_tree'] == species_tree]


class Node:
    def __init__(self, name, annotation, events, name2cluster, threshold):
        self.name = name
        self.annotation = annotation
        self.events = events
        self.name2cluster = name2cluster
        self.raw_members = events.events[events.events['node'] == self.name]
        self.members = self.raw_members[self.raw_members["copies"] > threshold]
        self.members['clst'] = self.members.apply(lambda row: row['cluster'].split('.')[0], axis=1)
        self.raw_members['clst'] = self.raw_members.apply(lambda row: row['cluster'].split('.')[0], axis=1)
        self.header = ["cluster", "Best_NOG", "Best_NOG_category", "Best_NOG_description", "duplications",
                                 "transfers", "losses", "originations", "copies", "NOG", '#NOG_annotations',
                                 "NOG_category", "NOG_description"]

    def get_best_cog_info(self, clst):
        best_cog = list(self.name2cluster[clst].cog_members.keys())[0]
        best_cog_cat = self.annotation.cog2cat[list(self.name2cluster[clst].cog_members.keys())[0]]
        best_cog_desc = self.annotation.cog2desc[list(self.name2cluster[clst].cog_members.keys())[0]]
        return (best_cog, best_cog_cat, best_cog_desc)

    def get_all_cog_info(self, clst):
        cogs = to_ssv(*self.name2cluster[clst].cog_members.keys())
        counts = to_ssv(*[str(i) for i in self.name2cluster[clst].cog_members.values()])
        cats = to_ssv(*set(self.annotation.cog2cat[cog] for cog in self.name2cluster[clst].cog_members.keys()))
        descs = to_ssv(*set(self.annotation.cog2desc[cog] for cog in self.name2cluster[clst].cog_members.keys()))
        return (cogs, counts, cats, descs)

    def __str__(self):
        lines = []
        lines.append(to_tsv(*self.header))
        for line in self.members.iterrows():
            clst = line[1]['clst']
            best_cog, best_cog_cat, best_cog_desc = self.get_best_cog_info(clst)
            events = [str(i) for i in line[1][['duplications','transfers', 'losses', 'originations', 'copies']]]
            cogs, counts, cats, descs = self.get_all_cog_info(clst)
            clst_line = to_tsv(
                        clst,
                        best_cog,
                        best_cog_cat,
                        best_cog_desc,
                        *events,
                        cogs,
                        counts,
                        cats,
                        descs)
            lines.append(clst_line)
        return "\n".join(lines)

    def gains_losses(self, descendant):
        gains_df = pd.DataFrame(columns=self.raw_members.columns)
        losses_df = pd.DataFrame(columns=self.raw_members.columns)
        for clst in self.raw_members['clst']:
            if (self.raw_members[self.raw_members['clst'] == clst]['copies'] < 0.5).tolist()[0] and \
            (descendant.raw_members[descendant.raw_members['clst'] == clst]['copies'] >= 0.5).tolist()[0]:
                gains_df = gains_df.append(descendant.raw_members[descendant.raw_members['clst'] == clst])
            elif (self.raw_members[self.raw_members['clst'] == clst]['copies'] >= 0.5).tolist()[0] and \
            (descendant.raw_members[descendant.raw_members['clst'] == clst]['copies'] < 0.5).tolist()[0]:
                losses_df = losses_df.append(descendant.raw_members[descendant.raw_members['clst'] == clst])
        return self.print_gains(gains_df), self.print_losses(losses_df)

    def print_gains(self, gains_df):
        lines = []
        lines.append(to_tsv("cluster", "type of gain", "Best_NOG", "Best_NOG_category", "Best_NOG_description", "NOG cat", "NOGs", "NOG_description"))
        for line in gains_df.iterrows():
            if line[1]['transfers'] >= 0.5 and line[1]['originations'] >= 0.5:
                gain_type = "transfer/origination"
            elif line[1]['transfers'] >= 0.5:
                gain_type = "transfer"
            elif line[1]['originations'] >= 0.5:
                gain_type = "origination"
            else:
                gain_type = "-"
            clst = line[1]['clst']
            best_cog, best_cog_cat, best_cog_desc = self.get_best_cog_info(clst)
            cogs, counts, cats, descs = self.get_all_cog_info(clst)
            lines.append(to_tsv(clst,
                                gain_type,
                                best_cog,
                                best_cog_cat,
                                best_cog_desc,
                                cogs,
                                cats,
                                descs))
        return "\n".join(lines)

    def print_losses(self, losses_df):
        lines = []
        lines.append(to_tsv("cluster", "Best_NOG", "Best_NOG_category", "Best_NOG_description", "NOG cat", "NOGs", "NOG_description"))
        for line in losses_df.iterrows():
            clst = line[1]['clst']
            best_cog, best_cog_cat, best_cog_desc = self.get_best_cog_info(clst)
            cogs, counts, cats, descs = self.get_all_cog_info(clst)
            lines.append(to_tsv(clst,
                                best_cog,
                                best_cog_cat,
                                best_cog_desc,
                                cogs,
                                cats,
                                descs))
        return "\n".join(lines)
