#! ${params.python3}
import sys
import pandas as pd
import ete3
from ETE3_styles import node_style_basic
from ETE3_Utils import get_ancestor, set_node_style
from numpy import log


def layout(node):
    if node.is_leaf():
        N = ete3.AttrFace("name", fsize=14, fgcolor="black")
        ete3.faces.add_face_to_node(N, node, 0)
    if all([feat in node.features for feat in ['copies', 'transfers', 'duplications', 'losses', 'originations']]):
    # if 'copies' in node.features:
        #C = CircleFace(radius=log(node.weight), color="Red", style="circle")
        copies = ete3.RectFace(width=2, height=node.copies, fgcolor="Red", bgcolor="Red")
        duplications = ete3.RectFace(width=2, height=node.duplications, fgcolor="Green", bgcolor="Green")
        transfers = ete3.RectFace(width=2, height=node.transfers, fgcolor="Blue", bgcolor="Blue")
        losses = ete3.RectFace(width=2, height=node.losses, fgcolor="Orange", bgcolor="Orange")
        originations = ete3.RectFace(width=2, height=node.originations, fgcolor="Black", bgcolor="Black")
        ete3.faces.add_face_to_node(copies, node, 0, position="float")
        ete3.faces.add_face_to_node(duplications, node, 1, position="float")
        ete3.faces.add_face_to_node(losses, node, 2, position="float")
        ete3.faces.add_face_to_node(originations, node, 3, position="float")
        ete3.faces.add_face_to_node(transfers, node, 4, position="float")


def simple_tree_style():
    ts = ete3.TreeStyle()
    ts.layout_fn = layout
    ts.show_leaf_name = False
    ts.scale = 250
    return ts


def add_genome_size(treeO, treeC, species_map, events, filename):
    def scale(x):
        return int((x / 500)) * 2
    events.to_csv(filename.replace('pdf', 'csv'), sep='\t')
    for n in treeO.traverse():
        name = ''
        if n.is_leaf():
            name = species_map[n.name]
        else:
            name = get_ancestor(
                [treeC.get_leaves_by_name(species_map[l.name])[0] for l in n.get_leaves()]).name
        n.add_features(copies=scale(events.loc[name]['copies']))
        n.add_features(transfers=scale(events.loc[name]['transfers']))
        n.add_features(duplications=scale(events.loc[name]['duplications']))
        n.add_features(losses=scale(events.loc[name]['losses']))
        n.add_features(originations=scale(events.loc[name]['originations']))

    ts = simple_tree_style()
    set_node_style(treeO, node_style_basic)
    treeO.ladderize(direction=1)
    treeO.render(filename, tree_style=ts)


def main():
    header = ['species_tree', 'cluster_ID', 'node', 'duplications',
              'transfers', 'losses', 'originations', 'copies']
    df = pd.read_csv("events.txt", sep='\\t', names=header)

    species_map = {line.split()[0]: line.split()[1] for line in open("$species_map")}

    split_df = [df[df['species_tree'] == t] for t in set(df.species_tree)]
    for f in split_df:
        add_genome_size(
            ete3.Tree("$workflow.launchDir/$params.output_trees/%s.tree" %
                      set(f.species_tree).pop().replace('_clean', '_root')),
            ete3.Tree("$workflow.launchDir/$params.output_trees/%s.tree" %
                      set(f.species_tree).pop(), format=1),
            species_map,
            f.groupby('node').sum(),
            set(f.species_tree).pop() + ".pdf")


if __name__ == '__main__':
    main()
