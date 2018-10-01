import pandas as pd


def to_tsv(*args):
    return("\t".join(args))


def to_ssv(*args):
    return("; ".join(args))

class Counts:
    def __init__(self, count_file):
        pd.options.mode.chained_assignment = None
        self.header = ['cluster', 'node', 'single',
                  'multi', 'losses', 'gain', 'expansion', 'reduction']
        df = pd.read_csv(count_file, sep='\t', comment='#')
        df = df.drop(['Pattern', 'C0/e0,d0,l0,t0'], axis='columns')
        df = df.drop(0, axis='rows')
        nodes = list(set([c.split(':')[0] for c in  df.columns.tolist()]))
        nodes.remove('Family')
        node_dfs = []
        for n in nodes:
            node_df = df.iloc[:,[c.startswith('{}:'.format(n)) for c in df.columns.tolist()]]
            if node_df.shape[1] == 6:
                node_df.columns = ['single', 'multi', 'losses', 'gain', 'expansion', 'reduction']
                node_df['cluster'] = df['Family']
                node_df['node'] = n
                node_dfs.append(node_df)
        self.events = pd.concat(node_dfs)


class Node:
    def __init__(self, name, annotation, events, name2cluster, threshold):
        self.name = name
        self.annotation = annotation
        self.events = events
        self.name2cluster = name2cluster
        self.raw_members = events.events[events.events['node'] == self.name]
        self.members = self.raw_members[((self.raw_members["single"] > threshold) | (self.raw_members["multi"] > threshold)) | ((self.raw_members["gain"] > threshold) | (self.raw_members["losses"] > threshold))]
        self.members['clst'] = self.members.apply(lambda row: row['cluster'].split('.')[0], axis=1)
        self.raw_members['clst'] = self.raw_members.apply(lambda row: row['cluster'].split('.')[0], axis=1)
        self.header = ["cluster", "Best_NOG", "Best_NOG_category", "Best_NOG_description", 'single', 'multi', 'losses', 'gain', 'expansion', 'reduction',
                        "NOG", '#NOG_annotations', "NOG_category", "NOG_description"]

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
            clst = line[1]['cluster']
            best_cog, best_cog_cat, best_cog_desc = self.get_best_cog_info(clst)
            events = [str(i) for i in line[1][['single', 'multi', 'losses', 'gain', 'expansion', 'reduction']]]
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
        for clst in self.raw_members['cluster']:
            if (self.raw_members[self.raw_members['cluster'] == clst]['copies'] < 0.5).tolist()[0] and \
            (descendant.raw_members[descendant.raw_members['cluster'] == clst]['copies'] >= 0.5).tolist()[0]:
                gains_df = gains_df.append(descendant.raw_members[descendant.raw_members['cluster'] == clst])
            elif (self.raw_members[self.raw_members['cluster'] == clst]['copies'] >= 0.5).tolist()[0] and \
            (descendant.raw_members[descendant.raw_members['cluster'] == clst]['copies'] < 0.5).tolist()[0]:
                losses_df = losses_df.append(descendant.raw_members[descendant.raw_members['cluster'] == clst])
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
            clst = line[1]['cluster']
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
            clst = line[1]['cluster']
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
