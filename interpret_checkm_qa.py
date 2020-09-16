#!/usr/bin/env python

import click
import json
import pandas as pd

@click.command()
@click.option('--input', '-i', type=click.Path(exists=True), help='Absolute path of bin_stats.analyse.tsv')
@click.option('--output', '-o', help='Directory to output csv file')
@click.option('--n', '-n', default='checkm_qa_out', help='Name of output file without extension (optional)')


#'/mnt/c/Users/Eilish McMaster/OneDrive/Documents/2020/Honours/analyses/kea_validation/cp_bin_stats.analyze.tsv'
def workflow(i,o,n):
    checkm_qa = {}
    with open(i) as qa:
        for line in qa:
            mag, items = line.strip().split("\t")
            items = items.replace("'", '"')
            items = json.loads(items)
            checkm_qa[mag] = items

    df = pd.DataFrame(checkm_qa)
    df = pd.DataFrame.transpose(df)
    df.to_csv('%s/%s.csv' %(o,n))

if __name__ == "__main__":
    workflow()
