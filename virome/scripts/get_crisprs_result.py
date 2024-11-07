#-*-coding:utf-8-*-
import pandas as pd
import os
import argparse
parser = argparse.ArgumentParser(description="file merge")
parser.add_argument('--ha', type=str, help='host anno')
parser.add_argument('--sa', type=str, help='scrisprs anno')
parser.add_argument('--i', type=str, help='host scrisprs id')
parser.add_argument('--o', type=str, help='out file')
args = parser.parse_args()
print(args.ha,args.sa,args.i,args.o)
df_id=pd.read_csv(args.i,sep="\t",low_memory=False)
df_id.columns = ['query_ha','query_sa']
print(df_id)

df_ha=pd.read_csv(args.ha,sep="\t",quotechar='"',low_memory=False)
df_ha = df_ha.iloc[:, [0] + list(range(14, df_ha.shape[1]))]
df_ha = df_ha.rename(columns={'query': 'query_ha'})
print(df_ha)

df_sa=pd.read_csv(args.sa,sep="\t",quotechar='"',low_memory=False)
df_sa = df_sa.iloc[:, [0] + list(range(14, df_sa.shape[1]))]
df_sa = df_sa.rename(columns={'query': 'query_sa'})
print(df_sa)

df_merge1 = pd.merge(df_id,df_ha,left_on = 'query_ha',right_on = 'query_ha',how = 'left')
df_merge1=df_merge1.dropna()
df_merge2=pd.merge(df_merge1,df_sa,left_on = 'query_sa',right_on = 'query_sa',how = 'left')
df_merge2=df_merge2.dropna()
print(df_merge2)

df_merge2.to_csv(args.o,sep="\t",index=False)


