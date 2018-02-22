#! ${params.python3}

import shutil

def replace(string, patterns):
    for p in patterns:
        string = string.replace(p[0], p[1])
    return string

shutil.copyfile("$bootstrap", "${bootstrap.simpleName}.clean.ufboot")

inputfile = "map_genes.txt"
replacementfile = "${bootstrap.simpleName}.clean.ufboot"

patterns = set()
for line in open(inputfile):
    line = line.strip().split('\t')
    patterns.add((line[0], line[1]))

content = open(replacementfile).read()
content = replace(content, patterns)

content = content.replace("_","")
content = content.replace("-","")
content = content.replace("$params.separator","_")

with open(replacementfile, 'w') as outhandle:
    outhandle.write(content)
