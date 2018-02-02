#! /usr/bin/env python3

import argparse
import re
from collections import defaultdict


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", required=True,
                    help="list of names to be abbreviated")
parser.add_argument("-o", "--output", default="species_map.txt",
                    help="long species names and their abbreviations")
args = parser.parse_args()


def get_code(name, verbosity=0):
    name_parts = re.split('\.|_|:', name)
    code_base = [s[0:2].upper() for s in name_parts[0:2]]
    code_add = []
    if verbosity > 0 and len(name_parts) > 2:
        code_add = [s[0:verbosity].upper() for s in name_parts[2:]]
    elif verbosity > 0:
        code_base = [s[0:2 + verbosity].upper() for s in name_parts[0:2]]
    return "".join(code_base + code_add)


def main(args):
    verbosity = 0
    names = [line.strip() for line in open(args.input)]
    assert len(names) == len(set(names))
    name2code = {}
    while names:
        code2name = defaultdict(list)
        for n in names:
            code2name[get_code(n, verbosity)].append(n)
        for k,v in code2name.items():
            if len(v) == 1:
                name2code[v[0]] = k
                names.remove(v[0])
        verbosity += 1
        if verbosity > 10:
            raise Exception("some names are too similar, try differentiating better...")
    with open(args.output, 'w') as out:
        for k,v in name2code.items():
            print("{}\t{}".format(k,v), file=out)


if __name__ == '__main__':
    main(args)
