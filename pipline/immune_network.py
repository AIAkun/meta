#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import argparse
import os
import pandas as pd
import time

class Pipeline(object):
    """
    @fun: pipeline
    """

    def __init__(self):
        # 初始化变量
        super(Pipeline).__init__()
        self.bash_env = "bash"

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
    
        
        self.result_path = os.path.join(args.outDir, "00-CRISPRCasFinder")
        if os.path.exists(self.result_path):
            pass
        else:
            os.makedirs(self.result_path, 0o755)
        self.fastp_path = os.path.join(args.outDir, "01-MetaCRAST")
        if os.path.exists(self.fastp_path):
            pass
        else:
            os.makedirs(self.fastp_path, 0o755)
        self.ref_remove_path = os.path.join(args.outDir, "02-CRISPR_MapHostDB")
        if os.path.exists(self.ref_remove_path):
            pass
        else:
            os.makedirs(self.ref_remove_path, 0o755)
        self.metaphlan_path = os.path.join(args.outDir, "03-Spacers_MapvMAGs")
        if os.path.exists(self.metaphlan_path):
            pass
        else:
            os.makedirs(self.metaphlan_path, 0o755)
        self.kraken2_path = os.path.join(args.outDir, "03-kraken2")
        if os.path.exists(self.kraken2_path):
            pass
        else:
            os.makedirs(self.kraken2_path, 0o755)
        self.megahit_path = os.path.join(args.outDir, "04-Spacers_MapCRISPR")
        if os.path.exists(self.megahit_path):
            pass
        else:
            os.makedirs(self.megahit_path, 0o755)
        self.prodigal_path = os.path.join(args.outDir, "05-result")
        if os.path.exists(self.prodigal_path):
            pass
        else:
            os.makedirs(self.prodigal_path, 0o755)
        # self.cdhit_path = os.path.join(args.outDir, "06-cdhit")
        # if os.path.exists(self.cdhit_path):
        #     pass
        # else:
        #     os.makedirs(self.cdhit_path, 0o755)
        # self.emapper_path = os.path.join(args.outDir, "07-emapper")
        # if os.path.exists(self.emapper_path):
        #     pass
        # else:
        #     os.makedirs(self.emapper_path, 0o755)
        # self.samcount_path = os.path.join(args.outDir, "08-sam_count")
        # if os.path.exists(self.samcount_path):
        #     pass
        # else:
        #     os.makedirs(self.samcount_path, 0o755)
        # self.kegg_path = os.path.join(args.outDir, "09-emapper_kegg")
        # if os.path.exists(self.kegg_path):
        #     pass
        # else:
        #     os.makedirs(self.kegg_path, 0o755)

        # self.virus_path=os.path.join(args.outDir, "10-virus")
        # if os.path.exists(self.virus_path):
        #     pass
        # else:
        #     os.makedirs(self.virus_path, 0o755)

        print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),'分析结果文件夹已新建完成。','\n', file=self.log)

    def CRISPRCasFinder(self):
        """
        @fun: CRISPRCasFinder
        """
        print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),'CRISPRCasFinder','\n', file=self.log)
        perl /home/meta/metagenomic/softs/CRISPRCasFinder/CRISPRCasFinder.pl -in install_test/sequence.fasta -cas -keep -out CRISPRCasFinder_result 
        # cmd = "{} {} {} {} {} {} {} {}".format(self.bash_env,
        #                                         self.fastp_cmd,
        #                                         sam_r1,
        #                                         sam_r1_clean,
        #                                         sam_r2,
        #                                         sam_r2_clean,
        #                                         sam_html,
        #                                         sam_json)
        # os.system(command=cmd)

    def MetaCRAST(self):
        """
        @fun: MetaCRAST
        """
        print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),'MetaCRAST','\n', file=self.log)

        # cmd = "{} {} {} {} {} {} {} {}".format(self.bash_env,
        #                                         self.fastp_cmd,
        #                                         sam_r1,
        #                                         sam_r1_clean,
        #                                         sam_r2,
        #                                         sam_r2_clean,
        #                                         sam_html,
        #                                         sam_json)
        # os.system(command=cmd)


    def CRISPR_MapHostDB(self):
        """
        @fun: CRISPR_MapHostDB
        """
        print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),'CRISPR_MapHostDB','\n', file=self.log)

        # cmd = "{} {} {} {} {} {} {} {}".format(self.bash_env,
        #                                         self.fastp_cmd,
        #                                         sam_r1,
        #                                         sam_r1_clean,
        #                                         sam_r2,
        #                                         sam_r2_clean,
        #                                         sam_html,
        #                                         sam_json)
        # os.system(command=cmd)
    def Spacers_MapvMAGs(self):
        """
        @fun: Spacers_MapvMAGs
        """
        print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),'Spacers_MapvMAGs','\n', file=self.log)

        # cmd = "{} {} {} {} {} {} {} {}".format(self.bash_env,
        #                                         self.fastp_cmd,
        #                                         sam_r1,
        #                                         sam_r1_clean,
        #                                         sam_r2,
        #                                         sam_r2_clean,
        #                                         sam_html,
        #                                         sam_json)
        # os.system(command=cmd)

    def Spacers_MapCRISPR(self):
        """
        @fun: Spacers_MapCRISPR
        """
        print(time.strftime("%Y-%m-%d, %H:%M:%S\n"),'Spacers_MapCRISPR','\n', file=self.log)

        # cmd = "{} {} {} {} {} {} {} {}".format(self.bash_env,
        #                                         self.fastp_cmd,
        #                                         sam_r1,
        #                                         sam_r1_clean,
        #                                         sam_r2,
        #                                         sam_r2_clean,
        #                                         sam_html,
        #                                         sam_json)
        # os.system(command=cmd)


    def result():
        pass

    def main(self):
        pass
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="",
                                     usage="",
                                     description="")

    parser.add_argument("-m", "--mMAGs", dest="",
                        required=True, type=str, help="")
    parser.add_argument("-v", "--vMAGs", dest="",
                        required=True, type=str, help="")
    # parser.add_argument("-r", "--ref", dest="ref", 
    #                     type=str, help="ref genome bowtie2 index")
    args = parser.parse_args()
    run = Pipeline()
    run.main()
