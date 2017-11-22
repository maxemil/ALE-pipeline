#! /usr/bin/env python3
"""replace content in multiple files by patterns in another"""
import argparse
import subprocess


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", required=True,
                    help="pattern, TAB separated")
parser.add_argument("-f", "--files", required=True, nargs='+',
                    help="files that are going to be searched and replaced")
parser.add_argument("-b", "--backup", required=False, default=False,
                    help="switch for creating backup file with extension .bak")

args = parser.parse_args()


def main(inputfile, replacementfiles, backup):
    with open(inputfile) as f:
        patterns = set()
        for line in f:
            line = line.strip().split('\t')
            patterns.add('s/%s/%s/g' % (line[0], line[1]))
        if backup:
            subprocess.call(
                ['sed', '-r', '-i.bak', ";".join(patterns)] + replacementfiles)
        else:
            subprocess.call(
                ['sed', '-r', '-i', ";".join(patterns)] + replacementfiles)

if __name__ == "__main__":
    main(args.input, args.files, args.backup)
