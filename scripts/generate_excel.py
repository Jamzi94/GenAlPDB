#!/usr/bin/env python3
import pandas as pd, pathlib
root=pathlib.Path(__file__).resolve().parent.parent
summary=root/'logs'/'runtime_summary.csv'
if not summary.exists():
    print('run parse first'); exit()
df=pd.read_csv(summary)
def split(j):
    p=j.split('_')
    return pd.Series(dict(uniprot=p[0],gpu=p[1],cpus=p[2],mem=p[3],partition=p[4]))
df=df.join(df['jobname'].apply(split))
df.to_excel(root/'benchmark_results.xlsx',index=False)
print('Excel ready')
