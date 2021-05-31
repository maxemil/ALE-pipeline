#! ${params.python3}
print("${species_tree.baseName}")

import ete3

def get_ancestor(nodes):
    ancestor = nodes[0]
    for n in nodes:
        ancestor = ancestor.get_common_ancestor(n)
    return ancestor

tree = ete3.Tree('$species_tree')

new_root = get_ancestor([tree.get_leaves_by_name(s)[0] for s in $params.outgroup_taxa])
tree.set_outgroup(new_root)
print(tree.write(), file=open("${species_tree.baseName}_root.tree", 'w'))

tree.convert_to_ultrametric()
species_map = {line.split()[0]: line.split()[1] for line in open('map_species.txt')}
for leaf in tree.get_leaves():
    leaf.name = species_map[leaf.name]

with open("${species_tree.baseName}_clean.tree", 'w') as outhandle:
    print(tree.write(), file=outhandle)
