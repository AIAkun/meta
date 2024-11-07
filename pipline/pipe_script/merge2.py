#-*-coding:utf-8-*-
import pandas as pd
import os
import argparse
parser = argparse.ArgumentParser(description="file merge")
parser.add_argument('--nt', type=str, help='nt database')
parser.add_argument('--i', type=str, help='input file')
parser.add_argument('--o', type=str, help='out file')
args = parser.parse_args()
print(args.nt,args.i,args.o)
df=pd.read_csv(args.nt,sep="\t",quotechar='"',low_memory=False)
df2=pd.read_csv(args.i,sep="\t",low_memory=False)
df3 = pd.merge(df2,df,left_on = 'subject',right_on = 'subject',how = 'left')
df3.to_csv(args.o,sep="\t",index=False)


