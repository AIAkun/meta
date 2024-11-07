import pandas as pd 
import os
import sys

def metaphlan_taxonomy_stat(metaphlan_path,outdir_path):
    levels=['kingdom','phylum','class','order','family','genus','species','SGB']
    filenames=os.listdir(metaphlan_path)
    data=pd.DataFrame()
    for i in range(len(levels)):
        print(levels[i])
        lens=[]
        samplename=[]
        taxon_data=pd.DataFrame()
        for filename in filenames:
            if filename.endswith(levels[i]+'.txt'):
                file_path=os.path.join(metaphlan_path,filename)
                file=pd.read_table(file_path)
                lens.append(len(file))
                samplename.append(filename.split('_')[0])
                file['source']=filename.split('_')[0]
                taxon_data=pd.concat([taxon_data,file])
        try:
            taxon_data=taxon_data[['taxonomy','taxid','source']].reset_index(drop=True)
            taxonomy=[]
            num=[]
            taxon=[]
            source=[]
            taxid=[]
            for a,b in taxon_data.groupby('taxonomy'):
                taxonomy.append(a)
                taxid.append(list(b['taxid'])[0])
                num.append(len(b))
                source.append(','.join(list(b['source'])))
            taxon_data2=pd.DataFrame({'taxonomy':taxonomy,'taxid':taxid,'num':num,'source':source})
            taxon_data2_path=os.path.join(outdir_path,levels[i]+'_unique.txt')
            taxon_data2.to_csv(taxon_data2_path,sep='\t',index=False)
        except:
            print(taxon_data)
        tmp_data=pd.DataFrame({'sample':samplename,levels[i]:lens})
        print(tmp_data)
        if data.empty is True:
            data=tmp_data
        else:
            data=pd.merge(data,tmp_data,on=['sample'],how='left')
    sample_info_path=os.path.join(outdir_path,'sample_taxonomy.txt')
    data.to_csv(sample_info_path,sep='\t',index=False)
        
def main():
    if len(sys.argv) < 2:
        print("usage: python metaphlan_taxonomy_stat <input_dir> <outdir>")
        sys.exit(1)
    metaphlan_taxonomy_stat(sys.argv[1],sys.argv[2])
 
if __name__ == "__main__":
    main()