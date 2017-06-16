#! ${params.python3}

with open("$gene") as f, open('events.txt', 'w') as outhandle:
    for line in f:
        if any([line.startswith(x) for x in ['S_terminal_branch', 'S_internal_branch']]):
            line = line.split()
            print("\\t".join(['$species_tree', "${gene.baseName}"] + line[1:]), file=outhandle)
