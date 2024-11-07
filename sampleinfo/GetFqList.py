#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import time
import os

file_path = 'sampleinfo.txt'
data_path = '/media/star/EXTERNAL_USB/吴军华/24-3-7-开始上传/病毒原始下机数据2.17/'

fq_list = 'fq.list'
fq_list_f = open(fq_list,'a')
fq_list_f.write("#Sample\tR1\tR2\n")

fq_gz='fq_gz.txt'
fq_gz_f =open(fq_gz,'a')



with open(file_path, 'r') as f:
	next(f)		# skip first line
	lines = f.readlines()
	for line in lines:
		line=line.strip()
		line_list = line.split('\t')
		run = 'L03'
		number = line_list[1].split('_')[1]
		file_pre = line_list[0]+'_'+run+'_'+number
		
		samplename = line_list[0]+'_'+run+'_'+number+'_'+line_list[3]
		batch = line_list[0]

		R1 = data_path+batch+'/'+run+'/'+file_pre+'_1.fq.gz'
		R2 = data_path+batch+'/'+run+'/'+file_pre+'_2.fq.gz'

		if os.path.exists(R1) and os.path.exists(R2):
			# print(samplename,'\t')
			# print(R1,'\t')
			# print(R2,'\n')
			content = samplename+'\t'+R1+'\t'+R2+'\n'
			fq_list_f.write(content)
			content2 = R1+'\n'+R2+'\n'
			fq_gz_f.write(content2)
		else:
			print(file_pre,'不存在！')

f.close()

