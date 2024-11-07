#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import argparse
import os
import pandas as pd
import time

class PipeMetagenome(object):
    """
    @fun: pipeline of metagenome
    """

    def __init__(self):
        # 初始化变量
        super(PipeMetagenome).__init__()
        self.bash_env = "bash"
        self.py_env = "python"
        self.pipe_script_path = '/home/meta/metagenomic/pipeline/pipe_script/'
        self.metaphlan_scripts_path ='/home/meta/anaconda3/envs/metagenomic/lib/python3.7/site-packages/metaphlan/utils/'

        self.fastp_cmd = self.pipe_script_path+"fastp.sh"
        self.ref_remove_cmd = self.pipe_script_path+"ref_remove.sh"
        self.metaphlan_cmd = self.pipe_script_path+"metaphlan.sh"
        self.metaphlan_merge_cmd = self.pipe_script_path+"metaphlan_merge.sh"
        self.metaphlan_merge = self.metaphlan_scripts_path+"merge_metaphlan_tables.py"
        self.hclust2_cmd = self.pipe_script_path+"hclust2.sh"
        self.graphlan_cmd = self.pipe_script_path+"graphlan.sh"
        self.megahit_cmd = self.pipe_script_path+"megahit.sh"
        self.prodigal_cmd = self.pipe_script_path+"prodigal.sh"
        self.cdhit_cmd = self.pipe_script_path+"cd-hit.sh"
        self.emapper_cmd = self.pipe_script_path+"emapper.sh"
        self.count_cmd = self.pipe_script_path+"gene_count.sh"
        self.kegg_cmd = self.pipe_script_path+"kegg/kegg.py"

        self.alpha_cmd = self.pipe_script_path+"alpha.sh"
        self.beta_cmd = self.pipe_script_path+"beta.sh"

        # 新建结果文件夹
        if os.path.exists(args.outDir):
            print("{} outDir already exist, will overwrite!".format(args.outDir))
            self.log = open((args.outDir + "/run.log"),
                            mode="a", encoding="utf-8")
            # quit()
        if not os.path.exists(args.outDir):
            os.mkdir(args.outDir)
            self.log = open((args.outDir + "/run.log"),
                            mode="a", encoding="utf-8")
        chmod_log = 'sudo chmod 777 {0}/run.log'.format(args.outDir)
        os.system(chmod_log)
        print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),'开始本次宏基因组数据分析：',args.outDir,'\n', file=self.log)
        
        self.result_path = os.path.join(args.outDir, "00-result")
        if os.path.exists(self.result_path):
            pass
        else:
            os.makedirs(self.result_path, 0o755)
        self.fastp_path = os.path.join(args.outDir, "01-fastp_trim")
        if os.path.exists(self.fastp_path):
            pass
        else:
            os.makedirs(self.fastp_path, 0o755)
        self.ref_remove_path = os.path.join(args.outDir, "02-ref_remove")
        if os.path.exists(self.ref_remove_path):
            pass
        else:
            os.makedirs(self.ref_remove_path, 0o755)
        self.metaphlan_path = os.path.join(args.outDir, "03-metaphlan")
        if os.path.exists(self.metaphlan_path):
            pass
        else:
            os.makedirs(self.metaphlan_path, 0o755)
        self.megahit_path = os.path.join(args.outDir, "04-megahit")
        if os.path.exists(self.megahit_path):
            pass
        else:
            os.makedirs(self.megahit_path, 0o755)
        self.prodigal_path = os.path.join(args.outDir, "05-prodigal")
        if os.path.exists(self.prodigal_path):
            pass
        else:
            os.makedirs(self.prodigal_path, 0o755)
        self.cdhit_path = os.path.join(args.outDir, "06-cdhit")
        if os.path.exists(self.cdhit_path):
            pass
        else:
            os.makedirs(self.cdhit_path, 0o755)
        self.emapper_path = os.path.join(args.outDir, "07-emapper")
        if os.path.exists(self.emapper_path):
            pass
        else:
            os.makedirs(self.emapper_path, 0o755)
        self.samcount_path = os.path.join(args.outDir, "08-sam_count")
        if os.path.exists(self.samcount_path):
            pass
        else:
            os.makedirs(self.samcount_path, 0o755)
        self.kegg_path = os.path.join(args.outDir, "09-emapper_kegg")
        if os.path.exists(self.kegg_path):
            pass
        else:
            os.makedirs(self.kegg_path, 0o755)

        print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),'分析结果文件夹已新建完成。','\n', file=self.log)
        # 读取fqlist文件
        self.raw_dict = {}
        with open(args.fqList) as f:
            for line in f:
                if line.startswith("#"):
                    pass
                else:
                    line = line.strip().split("\t")
                    self.raw_dict[line[0]] = line[1:]
        # print(self.raw_dict)
        print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),'fqList文件读取完成。','\n', file=self.log)

    def fastp_trim(self):
        """
        @fun: fastp trim
        """
        print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),'fastp trim','\n', file=self.log)
        for k, v in self.raw_dict.items():
            sam_name = k
            sam_r1 = v[0]
            sam_r2 = v[1]
            sam_fastp_path = os.path.join(self.fastp_path, sam_name)
            if os.path.exists(sam_fastp_path):
                pass
            else:
                os.makedirs(sam_fastp_path, 0o755)
            sam_r1_clean = sam_fastp_path + \
                "/{}_clean.1.fastq.gz".format(sam_name)
            sam_r2_clean = sam_fastp_path + \
                "/{}_clean.2.fastq.gz".format(sam_name)
            sam_html = sam_fastp_path + "/{}.html".format(sam_name)
            sam_json = sam_fastp_path + "/{}.json".format(sam_name)
            cmd = "{} {} {} {} {} {} {} {}".format(self.bash_env,
                                                    self.fastp_cmd,
                                                    sam_r1,
                                                    sam_r1_clean,
                                                    sam_r2,
                                                    sam_r2_clean,
                                                    sam_html,
                                                    sam_json)
            os.system(command=cmd)
            print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),cmd,'\n', file=self.log)
            print("Fastp finished: {}, check result dir: {}".format(
                sam_name, sam_fastp_path))
        print("All fastp finished, check dir: {}".format(self.fastp_path))

    def ref_remove(self):
        """
        @run: reference sequence remove
        """
        print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),'reference sequence remove','\n', file=self.log)
        for k, v in self.raw_dict.items():
            sam_name = k
            sam_ref_remove_path = os.path.join(self.ref_remove_path, sam_name)
            if os.path.exists(sam_ref_remove_path):
                os.system("rm -rf {}/*.sam".format(sam_ref_remove_path))
                pass
            else:
                os.makedirs(sam_ref_remove_path, 0o755)
            sam_r1_clean = self.fastp_path + \
                "/{}/{}_clean.1.fastq.gz".format(sam_name, sam_name)
            sam_r2_clean = self.fastp_path + \
                "/{}/{}_clean.2.fastq.gz".format(sam_name, sam_name)
            sam_ref_sam = sam_ref_remove_path + "/{}.sam".format(sam_name)
            sam_ref_sam_log = sam_ref_remove_path + \
                "/{}.mapping.log".format(sam_name)
            sam_r1_unmap = sam_ref_remove_path + \
                "/{}.unmap.1.fastq.gz".format(sam_name)
            sam_r2_unmap = sam_ref_remove_path + \
                "/{}.unmap.2.fastq.gz".format(sam_name)
            sam_unmap_single = sam_ref_remove_path + \
                "/{}.unmap.single.fastq.gz".format(sam_name)
            cmd = "{} {} {} {} {} {} {} {} {} {}".format(self.bash_env,
                                                            self.ref_remove_cmd,
                                                            args.ref,
                                                            sam_r1_clean,
                                                            sam_r2_clean,
                                                            sam_ref_sam,
                                                            sam_ref_sam_log,
                                                            sam_r1_unmap,
                                                            sam_r2_unmap,
                                                            sam_unmap_single)
            os.system(command=cmd)
            os.system("rm -rf {}".format(sam_ref_sam))
            print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),cmd,'\n', file=self.log)
            print("Ref sequence remove finished: {}, check result dir: {}".format(
                sam_name, sam_ref_remove_path))
        print("All ref sequence remove finished, check dir: {}".format(
            self.ref_remove_path))

    def metaphlan(self):
        """
        @fun: metaphlan analysis pipeline
        """
        print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),'metaphlan analysis pipeline','\n', file=self.log)
        for k, v in self.raw_dict.items():
            sam_name = k
            sam_metaphlan_path = os.path.join(self.metaphlan_path, sam_name)
            if os.path.exists(sam_metaphlan_path):
                pass
            else:
                os.makedirs(sam_metaphlan_path, 0o755)
            sam_r1_unmap = self.ref_remove_path + \
                "/{}/{}.unmap.1.fastq.gz".format(sam_name, sam_name)
            sam_r2_unmap = self.ref_remove_path + \
                "/{}/{}.unmap.2.fastq.gz".format(sam_name, sam_name)
            sam_metaphlan_bz2 = sam_metaphlan_path + \
                "/{}_bowtie2.bz2".format(sam_name)
            sam_metaphlan_out = sam_metaphlan_path + \
                "/{}.tsv".format(sam_name)
            cmd = "{} {} {} {} {} {}".format(self.bash_env,
                                                self.metaphlan_cmd,
                                                sam_r1_unmap,
                                                sam_r2_unmap,
                                                sam_metaphlan_bz2,
                                                sam_metaphlan_out)
            os.system(command=cmd)
            print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),cmd,'\n', file=self.log)
            print("Metaphlan finished: {}, check result dir: {}".format(
                sam_name, sam_metaphlan_path))
        metaphlan_merge = self.metaphlan_path + "/00_merged_abundance_table.txt"
        metaphlan_phylum = self.metaphlan_path + "/01_metaphlan_phylum.txt"
        metaphlan_class = self.metaphlan_path + "/02_metaphlan_class.txt"
        metaphlan_order = self.metaphlan_path + "/03_metaphlan_order.txt"
        metaphlan_family = self.metaphlan_path + "/04_metaphlan_family.txt"
        metaphlan_genus = self.metaphlan_path + "/05_metaphlan_genus.txt"
        metaphlan_species = self.metaphlan_path + "/06_metaphlan_species.txt"

        metaphlan_speciesCH = self.metaphlan_path + "/06_metaphlan_speciesCH.txt"
        metaphlan_format = self.metaphlan_path + "/00_metaphlan_format.txt"

        hclust2_png = self.metaphlan_path + '/hclust2.png'
        hclust2_legendpng =self.metaphlan_path + '/hclust2.legend.png'

        graphlantree = self.metaphlan_path + '/graphlantree.txt'
        graphlananno = self.metaphlan_path + '/graphlananno.txt'
        graphlanxml = self.metaphlan_path + '/graphlan.xml'
        graphlanpng = self.metaphlan_path + '/graphlan.png'

        if os.path.exists(metaphlan_merge):
            os.system(
                "cp {}/*.txt {}".format(self.metaphlan_path, self.result_path))
            pass
        else:
            cmd1 = "{} {}/*/*.tsv > {}".format(
                self.metaphlan_merge, self.metaphlan_path, metaphlan_merge)
            os.system(command=cmd1)
            print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),cmd1,'\n', file=self.log)
            cmd2 = "{} {} {} {} {} {} {} {} {}".format(self.bash_env,
                                                       self.metaphlan_merge_cmd,
                                                       metaphlan_merge,
                                                       metaphlan_phylum,
                                                       metaphlan_class,
                                                       metaphlan_order,
                                                       metaphlan_family,
                                                       metaphlan_genus,
                                                       metaphlan_species)
            os.system(command=cmd2)
            print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),cmd2,'\n', file=self.log)
            title_cmd = 'head -n 1 {0}'.format(metaphlan_species)
            title = (os.popen(title_cmd).read().strip())
            metaphlan_speciesCH_cmd = "tail +2 "+metaphlan_species+" |sed 's/s__/ /' |cut -d ' ' -f 2| sed '1i{0}' > {1}".format(title,metaphlan_speciesCH)
            os.system(metaphlan_speciesCH_cmd)
            # species热图
            cmd3 = "{} {} {} {}".format(self.bash_env,
                                            self.hclust2_cmd,
                                            metaphlan_speciesCH,
                                            hclust2_png
                                            )
            os.system(command=cmd3)
            print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),cmd3, '\n',file=self.log)
            metaphlan_format_cmd = "tail -n +2 {0} | cut -f1,2- > {1}".format(metaphlan_merge,metaphlan_format)
            os.system(metaphlan_format_cmd)
            # 树状图
            cmd4 = "{} {} {} {} {} {} {}".format(self.bash_env,
                                self.graphlan_cmd,
                                metaphlan_format,
                                graphlantree,
                                graphlananno,
                                graphlanxml,
                                graphlanpng
                                )
            os.system(command=cmd4)
            print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),cmd4,'\n', file=self.log)
            # alpha多样性
            cmd5 = "{} {} {} {}".format(self.bash_env,
                                self.alpha_cmd,
                                self.result_path,
                                metaphlan_merge                                
                                )
            os.system(command=cmd5)
            print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),cmd5,'\n', file=self.log)
            # beta多样性
            cmd6 = "{} {} {} {}".format(self.bash_env,
                                self.beta_cmd,
                                self.result_path,
                                metaphlan_merge
                                )
            os.system(command=cmd6)
            print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),cmd6,'\n', file=self.log)
            # os.system(
            #     "cp {}/*.* {}".format(self.metaphlan_path, self.result_path))
            
        print("Metaphlan merge finished, check dir: {}".format(self.metaphlan_path))

    def megahit(self):
        """
        @fun: sequence assemble
        """
        print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),'sequence assemble', '\n',file=self.log)
        for k, v in self.raw_dict.items():
            sam_name = k
            sam_megahit_path = os.path.join(self.megahit_path, sam_name)
            if os.path.exists(sam_megahit_path):
                os.system(
                    "rm -rf {}/intermediate_contigs".format(sam_megahit_path))
                continue
            else:
                # os.makedirs(sam_megahit_path,0o755) # The megahit output folder cannot be an existing one!
                sam_r1_unmap = self.ref_remove_path + \
                    "/{}/{}.unmap.1.fastq.gz".format(sam_name, sam_name)
                sam_r2_unmap = self.ref_remove_path + \
                    "/{}/{}.unmap.2.fastq.gz".format(sam_name, sam_name)
                sam_contig = sam_megahit_path + \
                    "/{}.contigs.fa".format(sam_name)
                # remove assembly contigs which < 500bp
                sam_contig_500 = sam_megahit_path + \
                    "/{}.contigs_500.fa".format(sam_name)
                cmd = "{} {} {} {} {} {} {} {}".format(self.bash_env,
                                                        self.megahit_cmd,
                                                        sam_r1_unmap,
                                                        sam_r2_unmap,
                                                        sam_megahit_path,
                                                        sam_name,
                                                        sam_contig,
                                                        sam_contig_500)
                os.system(command=cmd)
                # os.system(
                #     "rm -rf {}/intermediate_contigs".format(sam_megahit_path))
                print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),cmd,'\n', file=self.log)
                print("Megahit finished: {}, check result dir: {}".format(
                    sam_name, sam_megahit_path))
        print("All megahit finished, check dir: {}".format(self.megahit_path))

    def prodigal(self):
        """
        @fun: gene prediction
        """
        print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),'gene prediction', '\n',file=self.log)
        for k, v in self.raw_dict.items():
            sam_name = k
            sam_prodigal_path = os.path.join(self.prodigal_path, sam_name)
            if os.path.exists(sam_prodigal_path):
                pass
            else:
                os.makedirs(sam_prodigal_path, 0o755)
            sam_prot = sam_prodigal_path + "/{}_prot.faa".format(sam_name)
            sam_nucl = sam_prodigal_path + "/{}_nucl.fna".format(sam_name)
            sam_gff = sam_prodigal_path + "/{}_genes.gff".format(sam_name)
            sam_stat = sam_prodigal_path + "/{}.stat".format(sam_name)
            sam_contig_500 = self.megahit_path + \
                "/{}/{}.contigs_500.fa".format(sam_name, sam_name)
            cmd = "{} {} {} {} {} {} {}".format(self.bash_env,
                                                self.prodigal_cmd,
                                                sam_prot,
                                                sam_nucl,
                                                sam_gff,
                                                sam_stat,
                                                sam_contig_500)
            os.system(command=cmd)
            print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),cmd,'\n', file=self.log)
            print("Prodigal finished: {}, check result dir: {}".format(
                sam_name, sam_prodigal_path))
        print("All prodigal finished, check dir: {}".format(self.prodigal_path))

    def cdhit(self):
        """
        @fun: cdhit removes redundancy and non-redundant gene set builds bwa index
        """
        print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),'cdhit removes redundancy and non-redundant gene set builds bwa index', '\n',file=self.log)
        # combine gene sets from individual samples
        prot_cat = self.cdhit_path + "/prot.faa"
        nucl_cat = self.cdhit_path + "/nucl.fna"
        prot_nonerude = self.cdhit_path + "/prot_nonerude.faa"
        nucl_nonerude = self.cdhit_path + "/nucl_nonerude.fna"
        nonerude_list = self.cdhit_path + "/prot_nonerude.list"
        geneset_bwa = self.cdhit_path + "/geneset_bwa"
        geneset_length = self.cdhit_path + "/geneset_length.txt"
        if os.path.exists(geneset_length):
            print("CD-HIT result already exist,regenerate it, check dir: {}".format(self.cdhit_path))
            pass

        cmd1 = "cat {}/*/*_prot.faa > {}".format(
            self.prodigal_path, prot_cat)
        os.system(command=cmd1)
        print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),cmd1,'\n', file=self.log)
        cmd2 = "cat {}/*/*_nucl.fna > {}".format(
            self.prodigal_path, nucl_cat)
        os.system(command=cmd2)
        print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),cmd2,'\n', file=self.log)
        cmd = "{} {} {} {} {} {} {} {} {}".format(self.bash_env,
                                                    self.cdhit_cmd,
                                                    prot_cat,
                                                    prot_nonerude,
                                                    nonerude_list,
                                                    nucl_cat,
                                                    nucl_nonerude,
                                                    geneset_bwa,
                                                    geneset_length)
        os.system(command=cmd)
        print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),cmd,'\n', file=self.log)
        print("CD-HIT finished, check dir: {}".format(self.cdhit_path))

    def emapper(self):
        """
        @fun: function annotation using emapper and generate KO/KEGG information table
        """
        print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),'function annotation using emapper and generate KO/KEGG information table', '\n',file=self.log)
        prot_nonerude = self.cdhit_path + "/prot_nonerude.faa"
        eggnog_profile = self.emapper_path + "/eggnog"
        eggnog_ann = self.emapper_path + "/eggnog.emapper.annotations"
        eggnog_KO = self.emapper_path + "/KEGG_KO.txt"
        eggnog_path = self.emapper_path + "/KEGG_PATHWAY.txt"
        if os.path.exists(eggnog_path):
            print(
                "eggnog-mapper result already exist, regenerate it,check dir: {}".format(self.emapper_path))
            pass
        cmd = "{} {} {} {} {} {} {}".format(self.bash_env,
                                            self.emapper_cmd,
                                            prot_nonerude,
                                            eggnog_profile,
                                            eggnog_ann,
                                            eggnog_KO,
                                            eggnog_path)
        os.system(command=cmd)
        print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),cmd,'\n', file=self.log)
        print("eggnog-mapper finished, check dir: {}".format(self.emapper_path))

    def sam_gene_count(self):
        """
        @fun: gene abundance calculation for each sample and combining all samples
        """
        print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),'gene abundance calculation for each sample and combining all samples', '\n',file=self.log)
        count_list = []
        for k, v in self.raw_dict.items():
            sam_name = k
            sam_count_path = os.path.join(self.samcount_path, sam_name)
            if os.path.exists(sam_count_path):
                pass
            else:
                os.makedirs(sam_count_path, 0o755)
            geneset_bwa = self.cdhit_path + "/geneset_bwa"
            sam_r1_unmap = self.ref_remove_path + \
                "/{}/{}.unmap.1.fastq.gz".format(sam_name, sam_name)
            sam_r2_unmap = self.ref_remove_path + \
                "/{}/{}.unmap.2.fastq.gz".format(sam_name, sam_name)
            sam_count_bam = sam_count_path + \
                "/{}_mapping_geneset.bam".format(sam_name)
            sam_count_txt = sam_count_path + "/{}.count".format(sam_name)
            cmd = "{} {} {} {} {} {} {} {}".format(self.bash_env,
                                                    self.count_cmd,
                                                    geneset_bwa,
                                                    sam_r1_unmap,
                                                    sam_r2_unmap,
                                                    sam_count_bam,
                                                    sam_name,
                                                    sam_count_txt)
            os.system(command=cmd)
            print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),cmd,'\n', file=self.log)
            count_list.append(sam_count_txt)
            print("sample gene count finished, check dir: {}".format(
                self.samcount_path))
        # files = ['DBYX1FB.count', 'DBYX1FB.count.test1', 'DBYX1FB.count.test2']
        if os.path.exists(self.samcount_path + "/merged_file.txt"):
            print(
                "gene abundance table already exist, check file: {}/merged_file.txt".format(self.samcount_path))
            pass
        else:
            df_list = [pd.read_csv(f, sep='\t') for f in count_list]
            df_merged = pd.concat(df_list).groupby(
                'gene', as_index=False).sum()
            df_merged.fillna(0, inplace=True)
            merge_file = self.samcount_path + "/merged_file.txt"
            df_merged.to_csv(merge_file, sep='\t', index=False)
            print(
                "gene abundance calculation finished, check file: {}".format(merge_file))

    def eggnog_kegg(self):
        """
        @fun: obtain the KO/pathway result table
        """
        print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),'obtain the KO/pathway result table', '\n',file=self.log)
        eggnog_KO = self.emapper_path + "/KEGG_KO.txt"
        eggnog_path = self.emapper_path + "/KEGG_PATHWAY.txt"
        merge_file = self.samcount_path + "/merged_file.txt"
        out_ko = self.kegg_path + "/KO_samples.xls"
        out_pathway = self.kegg_path + "/pathway_samples.xls"
        tmp_dir = self.kegg_path + "/tmp"
        if os.path.exists(out_pathway):
            # os.system("cp {}/*.xls {}".format(self.kegg_path, self.result_path))
            print("eggnog_kegg result already exist,regenerate it, check dir: {}".format(
                self.kegg_path))
            pass
        
        cmd = "{} {} -kk {} -kp {} -mt {} -ok {} -op {} -t {}".format(self.py_env,
                                                                        self.kegg_cmd,
                                                                        eggnog_KO,
                                                                        eggnog_path,
                                                                        merge_file,
                                                                        out_ko,
                                                                        out_pathway,
                                                                        tmp_dir)
        os.system(command=cmd)
        os.system("cp {}/*.xls {}".format(self.kegg_path, self.result_path))
        print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),cmd,'\n', file=self.log)
        print("eggnog-kegg finished, check dir: {}".format(self.kegg_path))

    def main(self):
        self.fastp_trim()
        self.ref_remove()
        self.metaphlan()
        self.megahit()
        self.prodigal()
        self.cdhit()
        self.emapper()
        self.sam_gene_count()
        self.eggnog_kegg()
        self.log.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="PipeMetagenome",
                                     usage="\n=================================================================\n"
                                     "python pipe_metagenome.py\n"
                                     "\t--fastq_list fq.list\n"
                                     "\t--output_dir result\n"
                                     "\t--ref ref_bowtie2_index\n"
                                     "ref_bowtie2_index:\n"
                                     "canis: /root/database/Canis_GCF_000002285.5/Canis_GCF_000002285_5\n"
                                     "human: /root/database/hg38_GCF_000001405.40/GCF_000001405.40/hg38\n"
                                     "=================================================================",
                                     description="Pipeline of metagenome")
    # fastq_list must start with "#Sample\tR1\tR2"
    parser.add_argument("-l", "--fastq_list", dest="fqList",
                        required=True, type=str, help="raw fq list")
    parser.add_argument("-o", "--output_dir", dest="outDir",
                        required=True, type=str, help="result output")
    parser.add_argument("-r", "--ref", dest="ref", required=True,
                        type=str, help="ref genome bowtie2 index")
    args = parser.parse_args()
    run = PipeMetagenome()
    run.main()
