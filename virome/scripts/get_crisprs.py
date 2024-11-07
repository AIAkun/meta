import time


scrisprs_filename='virus_host.crisprs'
out=open('host_scrisprs.id', 'a')
out.write('host_id\tscrisprs_id\n')
out2=open('scrisprs.fasta', 'a')

with open(scrisprs_filename, 'r') as f:
    f=f.readlines()
    for line in f:
        if line.startswith('Sequence'):
            line_list = line.split("'")
            host_name = line_list[1]
            # print(host_name)
            n=1    
        elif line[0].isdigit():
            # print(line)
            line_list = line.split('\t')
            # print(line_list)
            if len(line_list)==5:
                scrisprs_sequence=line_list[3]
                # print(scrisprs_sequence)
                scrisprs_name=host_name+"_"+'S'+str(n)
                # print(scrisprs_name)
                out.write(host_name+"\t"+scrisprs_name+"\n")
                out2.write('>'+scrisprs_name+'\n'+scrisprs_sequence+'\n')
                n=n+1
                # time.sleep(1)