#!/usr/bin/python
# -*- coding: utf-8 -*-
# conda activate python3

import sys
import pandas as pd
def usage():
    print('Usage: python blast_best.py [input_file] [outfile]')


def main():

    df = pd.read_csv(sys.argv[1],sep='\t')
    df = df.sort_values("bit score", ascending=False).drop_duplicates("query", keep='first').reset_index(drop=True)
    df.to_csv(sys.argv[2],sep='\t',index=False)


try:
    main()
except IndexError:
    usage()