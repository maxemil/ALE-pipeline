#! ${params.python3}

import shutil

def replace(string, patterns):
    for p in patterns:
        string = string.replace(p[0], p[1])
    return string

def main(bootstrap, prefix, separator):
    inputfile = "map_genes.txt"
    replacementfile = "{}.clean.ufboot".format(prefix)

    shutil.copyfile(bootstrap, replacementfile)

    patterns = set()
    for line in open(inputfile):
        line = line.strip().split('\t')
        patterns.add((line[0], line[1]))

    content = open(replacementfile).read()
    content = replace(content, patterns)

    content = content.replace("_","")
    content = content.replace("-","")
    content = content.replace(separator,"_")
    content = content.replace(".","")

    with open(replacementfile, 'w') as outhandle:
        outhandle.write(content)


if __name__ == '__main__':
    main("$bootstrap", "${bootstrap.simpleName}", "$params.separator")
