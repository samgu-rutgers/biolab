# Gu lab bioinformatic functions 
# python 3.7
'''
g_pos_to_seq_ce10 (g_chr,g_start,g_end)
revcom(seq)
fa_align(fa_file,ref_seq)
fa_align_bowtie(fa_file,bowtie_ref,barcode_size,end3_trim_size)

dir_path(f,d): return the folder path that contains file f, the path ends with "/"
admera_transfer1(bf,tf)
folder_path_slash(f): add slash to the end of the folder path if it doesn't already have it
admera_sRNAseq_rawFq2fa_IDTlinker1(lib_info_file,report_prefix): parsing Admera sRNA-seq fastq files
cp_rename_chipLib(f): copy admera chipseq fastq.gz files or Smarter RNA-seq libs to demultiplex folder, change name to SG number f: ChIP_lib_info_admera_SG###.csv

admera_RNAseq_rawFq2fa_IDTlinker1(lib_info_file,report_prefix): parsing Admera RNA-seq fastq files, this is specifically for the libraries made with IDT-linker1 
RNAseq_fasta2rawCounts_pdDataframe_geneWise(L): generate gene-wise raw counts in pandas dataframe for RNAseq data
RNAseq_fastq2rawCounts_pdDataframe_geneWise(L): generate gene-wise raw counts for RNAseq data, use fastq files
RNAseq_fastq2rpkm_pdDataframe_geneWise(L): #L: a list of RNA-seq fa file name prefix, ['SG0420_lib11", etc]
get_ce6_cDNA_size(): # return a dict. Key: gene name, value: cDNA size
get_ce6_cDNA_size_v2(ref): # return a dict. Key: gene name, value: cDNA size
sRNAseq_fasta2rawCounts_pdDataframe_geneWise(L): generate gene-wise raw counts in pandas dataframe for sRNAseq data
sRNAseq_fasta2rawCounts_pdDataframe_geneWise_v2(L,ref): generate gene-wise raw counts in pandas dataframe for sRNAseq data

sRNAseq_fasta2rpkm_pdDataframe_geneWise(L)
sRNAseq_fasta2rpkm_pdDataframe_geneWise_v2(L,ref)
fq_align_bowtie(fq_file,bowtie_ref,remove_5_size,remove_3_size): align fastq file with bowtie, return map file path
fq_align_bowtie2_SAM(fq_file_1,fq_file_2,lib_name,bowtie_ref,remove_5_size,remove_3_size): align paired-end fastq file using bowtie2, output is SAM file
##stop using this one: ChIPseq_fastq2rpm_geneWise(L): generate gene-wise rpm values for ChIPseq data. L: a list of ChIPseq fa file names, ['SG0420_lib11", etc] # use ChIPseq_fastq2rpkm_geneWise instead
ChIPseq_fastq2rpkm_geneWise(L): generate gene-wise rpkm values for ChIPseq data. L: a list of ChIPseq fa file names,
ChIPseq_fastq2rawCounts_geneWise(L):generate gene-wise raw counts values for ChIPseq data. L: a list of ChIPseq fa file names, not normalized to gene size and sequencing depth

fq_align_bowtie_SAM

whole_chrom_coverage_ChIP_plot1(name1,name2):show whole chr. coverages for two chip-seq libs; eg (name1,name2)=('SG0820_lib1','SG0820_lib5')
get_total_aligned_number(fq_file,bowtie_ref,remove_5_size,remove_3_size,strand): #10/6/2020 return the total number of bowtie aligned reads. can handle both fa and fastq files, zip or unziped
total_piRNA_count(file_name,barcode_size): count total piRNA reads in a fasta or fastq file, return {'total_reads':c2/f_type,'piRNA_reads':tpc}
sRNA_diff_genomic_feature_counts(l, barcode_size)
sRNAseq2rawCounts_pdDataframe_geneWise(D)
fa_fq_align_bowtie(fq_file,bowtie_ref,remove_5_size,remove_3_size)
pval_beyes(T1,T2,C1,C2): calculate the pval using the Bayesian method
deseq2_sRNA_simple1(lib_list,exp_design,exp_name,working_folder): # generate a R scatter plot and the deseq2 csv file with fold changes and p values
deseq2_RNA_simple1(lib_list,exp_design,exp_name,working_folder):

siRNA_track_plot_cer3(a,l)
siRNA_track_plot_cer3_gag(a,l,m)
siRNA_track_plot_gene_102920(a,l,m):

siRNA_track_plot_oma1_smg1_fusion(l)
RNA_track_plot_oma1_smg1_fusion(l)

MA_plot_sRNAseq(x,y,HL_g)
MA_plot_RNAseq(x,y,HL_g)

bigWig_chip_pairend_ce10(lib_name,end_3_remove) # generate bigwig file using paired end chip-seq sequencing
bigWig_chip_pairend_hg19(lib,end_3_remove)

fastq_to_ucsc_bed_ce6(lib) # generate bed file for RNA-seq or sRNA-seq data
primer3_v1(seq,size_min,size_max)
get_cDNA_fasta(g)

def bigWig_RNAseq_singleEnd_ce10_v1(lib,end_3_remove): # generate bigwig file using single-ended RNA sequencing SG##_lib##_3linker_removed.fastq.gz
'''

import os
import sys
import re
import subprocess
import pandas as pd
import gzip
import numpy as np
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
import mplcursors

# change the following parameters based on your computer setup
HOME=os.getenv("HOME")
#fastq_folder='/home/sam/fastq_demultiplexed/'
fastq_folder='/mnt/da1/fastq_demultiplexed/' # this is where you store your fastq or fasta files. Most functions can search subfolders of fastq_folder for the specified file
#fastq_folder='/mnt/md0/GEO_submissions/'
#fastq_folder='/mnt/md0/raw_fastq_GuLab/'
#fastq_folder='/media/sam/Illumina1/fastq_demultiplexed/'
#fastq_folder='/Users/SG_MBP2012/Box/fastq_demultiplexed/'
#bowtie_folder=HOME+'/Dropbox/bioinformatics/bowtie_mac/'
bowtie_folder=HOME+'/Dropbox/bioinformatics/bowtie_linux/'
bowtie_ref_folder=HOME+'/Dropbox/bioinformatics/bowtie_linux/indexes/'
bowtie2_folder=HOME+'/Dropbox/bioinformatics/bowtie2_linux/'
bowtie2_ref_folder=HOME+'/Dropbox/bioinformatics/bowtie2_linux/indexes/'
map_folder='/mnt/da1/bowtie_map/'

samtoolsfolder='/home/sam/Dropbox/bioinformatics/samtools/bin/'
deeptoolfolder='/home/sam/Dropbox/bioinformatics/deeptool/bin/'
samfolder='/mnt/da1/bowtie_map/'

core_num=12 # use half of the core

def g_pos_to_seq_ce10 (g_chr,g_start,g_end): # get C. elegans genomic DNA sequence, all upper case, based on position (0-based, semi-inclusive), 5/4/2020
    if g_chr not in ('chrI','chrII','chrIII','chrIV','chrV'): sys.exit("check chromosome name.")
    if g_end <= g_start: sys.exit("check start and end positions.")
        
    g_seq_file=HOME+'/Dropbox/bioinformatics/Genome_ce10/ce10.fa'
    FH=open(g_seq_file,'r')
    seq=''
    find_chr=False
    for L in FH:
        if find_chr:
            #check if the positions are longer than the chromosome size:
            if L[0]=='>': sys.exit("target DNA out of chromosome range.")
            c+=1
            if c==0: 
                length_per_line=(len(L.strip()))
                (start_line,start_column) = (int(g_start/length_per_line),g_start%length_per_line)
                (end_line,end_column) = (int(g_end/length_per_line),g_end%length_per_line)
            if c==start_line:
                if start_line<end_line:
                    seq=L.strip().upper()[start_column:]
                if start_line == end_line:
                    seq=L.strip().upper()[start_column:end_column]
                    break
            if c>start_line and c<end_line:
                seq+=L.strip().upper()
            if c==end_line:
                seq+=L.strip()[0:end_column].upper()
                break

        if L[0]=='>' and L.strip()=='>'+g_chr:
            c=-1
            find_chr=True
            
    FH.close()
    return seq
    
def revcom(seq): # return reverse complimentary sequence
    re_seq=seq[::-1]
    revcom_seq=''
    for N in re_seq:
        if N not in ['a','g','c','t','A','G','C','T','n','N']: sys.exit("sequence has letters other than a,t,g,c,A,T,G,C,n,N.")
        if N=='a': revcom_seq+='t'
        if N=='t': revcom_seq+='a'
        if N=='g': revcom_seq+='c'
        if N=='c': revcom_seq+='g'
        if N=='A': revcom_seq+='T'
        if N=='T': revcom_seq+='A'
        if N=='G': revcom_seq+='C'
        if N=='C': revcom_seq+='G'
        if N=='N' or N=='n': revcom_seq+=N
    return revcom_seq

# align fasta file to a reference using bowtie
# return the map file name
def fa_align_bowtie(fa_file,bowtie_ref,barcode_size,end3_trim_size): 
#    print('fa_align_bowtie %s %s' %(fa_file,bowtie_ref))
    bowtie_map_file=map_folder+fa_file.split('/')[-1]+'_'+bowtie_ref+'.map'
    bowtie_prog=bowtie_folder+'bowtie '+bowtie_ref_folder+bowtie_ref+' '+fa_file+' '+bowtie_map_file+' -5 '+str(barcode_size)+' -f -t -v 0 -a -3 '+ str(end3_trim_size) +' -p '+str(core_num)
    if os.path.isfile(bowtie_map_file):
        print ('%s exist.'% bowtie_map_file)
        return bowtie_map_file
    else:
        p=subprocess.run(bowtie_prog, shell=True, check=True,capture_output=True)
        print(p.stderr)
        return bowtie_map_file

def fa_align(fa_file,ref_seq): #align fasta reads to a ref seq, return alignments, using pythons re.finditer, ok speed for light jobs, but slow for millions of reads against large reference
    hits_dict={'sense_hits':[], 'antisense_hits':[]} #keys: "sense_hits" and "antisense_hits", values: lists of [read_name,read_seq,[start positions]] 
    FH=open(fa_file,'r')
    for L in FH:
        if L[0]=='>': seqName=L.strip()[1:]
        else:
            print(L)
            querySeq=L.strip().upper()
            querySeq_revcom=revcom(querySeq)
            sense_hits=[m.start() for m in re.finditer(querySeq,ref_seq)]
            antisense_hits=[m.start() for m in re.finditer(querySeq_revcom,ref_seq)]
            if len(sense_hits)>0:
                hits_dict['sense_hits'].append([seqName,querySeq,sense_hits])
            if len(antisense_hits)>0:
                hits_dict['antisense_hits'].append([seqName,querySeq,antisense_hits])
    FH.close()
    return hits_dict




    
def dir_path(f,d):
# 8/27/2020 return the folder path that contains file f
# if directory d does not have f, exit the program
    for dirpath,dirnames,filenames in os.walk(d):
        if f in filenames:
            return folder_path_slash(dirpath)  
    sys.exit('%s does not have file %s' % (d,f))
    


def admera_transfer1(bf,tf):
# 9/1/2020
# cp raw fastq files from intial downn load folder (bf) to target folder(tf)  
# also compare checksum to make sure the files are the same
    bf=folder_path_slash(bf)
    tf=folder_path_slash(tf)
    fl1=[bf+f for f in os.listdir(bf)]
    for f1 in fl1: 
        #print (f1)
        if os.path.isdir(f1):
            fl2=[f1+'/'+f for f in os.listdir(f1)]
            for f2 in fl2:                
                cmd = 'cp %s %s' %(f2,tf)
                print(cmd)
                os.system(cmd)
                # check checksum values are the same
                #cmd1='md5sum %s' % f2
                #cmd2='md5sum %s' % tf+f2.split('/')[-1]
                #p1=subprocess.run(cmd1, shell=True, check=True,capture_output=True)
                #p2=subprocess.run(cmd2, shell=True, check=True,capture_output=True)
                #print (f2.split('/')[-1])
                #warning='good' if (p1.stdout.split(b' ')[0])==(p1.stdout.split(b' ')[0]) else 'warning: checksum not the same'   
                #print(warning)

def folder_path_slash(f):
# 9/1/2020 add slash to the end of the folder path if it doesn't already have it
    a='' if f[-1]=='/' else '/'
    return f+a

def admera_sRNAseq_rawFq2fa_IDTlinker1(lib_info_file,report_prefix):
# 9/1/2020
# this script is specific to Admera raw fastq data + IDT linker1 ligation based sRNA-seq lib + multiple 4nt inline barcodes
# 150 nt paired end run, but we'll only use R1 reads
# input needs to be  ##.fastq.gz
# For collapsing, only collapse identical reads with indentical barcodes
    # format for the lib_info_file: 
    # line 1: header
    # line 2 and later:
    #         Sample ID (SG number), lib_name (SG number), index, inline barcodes, flowcell ID (the 9 characters before the first '_' of the fastq.gz file), flowcell lane number (the number after the first '_'),
    #         folder of the original fastq, folder of the output1 (raw fastq), folder of output2 (3'linker trmmed fastq), folder of output 3 (collapsed fasta)
    #         ori_file_extension: '_1.fastq', no need to change
    # If multiple barcodes are used for a library, barcodes are seperated by ' ' (a space)
# output
# 1: raw fastq reads (only change the file name for GEO submission)
# 2: fastq reads with 3'linker trimmed but keep the inline barcodes, no collapsing
# 3: fasta reads with 3'linker trimmed but keep the inline barcodes, identicalreads (including the barcode) collapsed, sequence name will tell the number of reads
    final_report='%s_RNAseq_parsing_output.txt' % report_prefix # need to be changed 
    infile_lib_info=open(lib_info_file,'r')
    outfile_report=open(final_report,'w')
    run_length=150 
    min_size=18 #minimal size for RNA seq
    #max_size=47
    linker3='CTGTAGGCACCATCAATC'
    e=0
    for L1 in infile_lib_info:
        e+=1
        if e==1:
            p=re.compile(',')
            m=p.search(L1)
            sep='\t' if m == None else ','
            outfile_report.write(L1.strip()+'\ttotal_fastq_reads_number(millions)\ttotal_fastq_reads_with_barcodes(millions)\tcollapsed_fastq_reads_number(millions)\n')
        else:
            (num_fastq_raw,num_fastq_raw_with_barcodes, num_fasta) = (0.0,0.0, 0.0) # number of fastq reads with barcodes; number of fasta reads, collapsed.
            A1=L1.strip().split(sep)
#            print (A1)
            (sampleID,lib_name,index,barcode,fcID,S_number,ori_folder)= (A1[0],A1[1],A1[2],A1[3].split(' '),A1[4],A1[5],A1[6]) #[1:-1]: remove the quatation marks
            (output1_folder,output2_folder,output3_folder,ori_fastq_extension)=(A1[7],A1[8],A1[9],A1[10])
            file_name=ori_folder+sampleID+'_'+S_number+ori_fastq_extension
            for B in barcode:
                if len(B)!= 4: sys.exit("check barcode info in the input file. Each barcode needs to be followed by a comma and a space, except the last one")
            print (sampleID,barcode,lib_name)
            print (file_name)
            outfile1=open(output1_folder+lib_name+'_raw.fastq','w')
            outfile2=open(output2_folder+lib_name+'_3linker_removed.fastq','w')
            cmd_unzip="gunzip "+file_name+'.gz'
            os.system(cmd_unzip)
            infile=open(file_name,'r')
            c=-1
            for L in infile:
                c+=1
                if c%4==0:
                    lines=['','','','']
                    num_fastq_raw+=1.0
                lines[c%4]=L.strip()
                if c%4==3 and (lines[1][0:4] in barcode):
                    insert=''
                    p1=re.compile(linker3[0:12]) # use the frist 12 nt in the IDT linke 1 to search for the 3' linker
                    m=p1.search(lines[1])
                    if m!=None: # if found perfect match of the 3' linker
                        insert=lines[1][0:m.start()]                   
                    if len(insert)>= min_size:
                        num_fastq_raw_with_barcodes+=1
                        outfile1.write(lines[0]+'\n')
                        outfile2.write(lines[0]+'\n')
                        outfile1.write(lines[1]+'\n')
                        outfile2.write(insert+'\n')
                        outfile1.write(lines[2]+'\n')
                        outfile2.write(lines[2]+'\n')
                        outfile1.write(lines[3]+'\n')
                        outfile2.write(lines[3][0:len(insert)]+'\n') 

            infile.close()
            cmd_zip="gzip "+file_name
            os.system(cmd_zip)
            outfile1.close()
            outfile2.close()
            # output 3
            fastq_count=0
            seq_dict={}
            fastq_file=open(output2_folder+lib_name+'_3linker_removed.fastq','r')
            c=0
            for L2 in fastq_file:
                c+=1
                if c%4==2: 
                    fastq_count+=1
                    seq=L2.strip()
                    if seq in seq_dict: seq_dict[seq]+=1
                    else:seq_dict[seq]=1
            fastq_file.close()
            hist_count={}
            max_count=0
            outfile_fa=open(output3_folder+lib_name+'_3linker_removed_collapsed.fa','w')

            for s in seq_dict:
                outfile_fa.write('>'+lib_name+'_'+str(num_fasta)+'_count_'+str(seq_dict[s])+'\n')
                num_fasta+=1
                outfile_fa.write(s+'\n')
                if seq_dict[s]>max_count: max_count=seq_dict[s]
                if seq_dict[s] in hist_count: hist_count[seq_dict[s]]+=1
                else: hist_count[seq_dict[s]]=1
            outfile_fa.close()
            report_file=open(output3_folder+lib_name+'_collapsed_fasta_report.txt','w')

            report_file.write(lib_name+'\n')
            report_file.write('fastq reads: '+str(fastq_count)+'\n')
            print ('fastq reads: '+str(fastq_count))
            report_file.write('collapsed fasta reads: '+str(num_fasta)+'\n\n')
            print ('collapsed fasta reads: '+str(num_fasta))
            report_file.write('number of identical seq\tnumber of seq\n')
            for k in range(1,max_count+1):
                if k in hist_count:
                    report_file.write(str(k)+'\t'+str(hist_count[k])+'\n')
            report_file.close()
            outfile_report.write('%s\t%.2f\t%.2f\t%.2f\n' % (L1.strip(), num_fastq_raw/1000000, num_fastq_raw_with_barcodes/1000000, num_fasta/1000000))
            cmd_zip="gzip "+output1_folder+lib_name+'_raw.fastq'
            os.system(cmd_zip)
            cmd_zip="gzip "+output2_folder+lib_name+'_3linker_removed.fastq'
            os.system(cmd_zip)
    infile_lib_info.close()
    outfile_report.close()
    
def admera_RNAseq_rawFq2fa_IDTlinker1(lib_info_file,report_prefix):
# 9/15/2020
# this script is specific to Admera raw fastq data + IDT linker1 ligation based RNA-seq lib + multiple 4nt inline barcodes
# 151 nt paired end run, but we'll only use R1 reads
# input needs to be  ##.fastq.gz
# For collapsing, only collapse identical reads with indentical barcodes
    # format for the lib_info_file: 
    # line 1: header
    # line 2 and later:
    #         Sample ID (SG number), lib_name (SG number), index, inline barcodes, flowcell ID (the 9 characters before the first '_' of the fastq.gz file), flowcell lane number (the number after the first '_'),
    #         folder of the original fastq, folder of the output1 (raw fastq), folder of output2 (3'linker trmmed fastq), folder of output 3 (collapsed fasta)
    #         ori_file_extension: '_1.fastq', no need to change
    # If multiple barcodes are used for a library, barcodes are seperated by ' ' (a space)
# output
# 1: raw fastq reads (only change the file name for GEO submission)
# 2: fastq reads with 3'linker trimmed but keep the inline barcodes, no collapsing
# 3: fasta reads with 3'linker trimmed but keep the inline barcodes, identicalreads (including the barcode) collapsed, sequence name will tell the number of reads
    final_report='%s_RNAseq_parsing_output.txt' % report_prefix # need to be changed 
    infile_lib_info=open(lib_info_file,'r')
    outfile_report=open(final_report,'w')
    run_length=151 
    min_size=28 #minimal size for RNA seq
    #max_size=47
    linker3='CTGTAGGCACCATCAATC'
    e=0
    for L1 in infile_lib_info:
        e+=1
        if e==1:
            p=re.compile(',')
            m=p.search(L1)
            sep='\t' if m == None else ','
            outfile_report.write(L1.strip()+'\ttotal_fastq_reads_number(millions)\ttotal_fastq_reads_with_barcodes(millions)\tcollapsed_fastq_reads_number(millions)\n')
        else:
            (num_fastq_raw,num_fastq_raw_with_barcodes, num_fasta) = (0.0,0.0, 0.0) # number of fastq reads with barcodes; number of fasta reads, collapsed.
            A1=L1.strip().split(sep)
#            print (A1)
            (sampleID,lib_name,index,barcode,fcID,S_number,ori_folder)= (A1[0],A1[1],A1[2],A1[3].split(' '),A1[4],A1[5],A1[6]) #[1:-1]: remove the quatation marks
            (output1_folder,output2_folder,output3_folder,ori_fastq_extension)=(A1[7],A1[8],A1[9],A1[10])
            file_name=ori_folder+sampleID+'_'+S_number+ori_fastq_extension
            for B in barcode:
                if len(B)!= 4: sys.exit("check barcode info in the input file. Each barcode needs to be followed by a comma and a space, except the last one")
            print (sampleID,barcode,lib_name)
            print (file_name)
            outfile1=open(output1_folder+lib_name+'_raw.fastq','w')
            outfile2=open(output2_folder+lib_name+'_3linker_removed.fastq','w')
            cmd_unzip="gunzip "+file_name+'.gz'
            os.system(cmd_unzip)
            infile=open(file_name,'r')
            c=-1
            for L in infile:
                c+=1
                if c%4==0:
                    lines=['','','','']
                    num_fastq_raw+=1.0
                lines[c%4]=L.strip()
                if c%4==3 and (lines[1][0:4] in barcode):
                    insert=''
                    p1=re.compile(linker3[0:12]) # use the frist 12 nt in the IDT linke 1 to search for the 3' linker
                    m=p1.search(lines[1])
                    if m!=None: # if found perfect match of the 3' linker
                        insert=lines[1][0:m.start()]                   
                        if len(insert)>= min_size:
                            num_fastq_raw_with_barcodes+=1
                            outfile1.write(lines[0]+'\n')
                            outfile2.write(lines[0]+'\n')
                            outfile1.write(lines[1]+'\n')
                            outfile2.write(insert+'\n')
                            outfile1.write(lines[2]+'\n')
                            outfile2.write(lines[2]+'\n')
                            outfile1.write(lines[3]+'\n')
                            outfile2.write(lines[3][0:len(insert)]+'\n') 
                    else:
                        insert=line[1]
                        num_fastq_raw_with_barcodes+=1
                        outfile1.write(lines[0]+'\n')
                        outfile2.write(lines[0]+'\n')
                        outfile1.write(lines[1]+'\n')
                        outfile2.write(insert+'\n')
                        outfile1.write(lines[2]+'\n')
                        outfile2.write(lines[2]+'\n')
                        outfile1.write(lines[3]+'\n')
                        outfile2.write(lines[3][0:len(insert)]+'\n')                         
                    

            infile.close()
            cmd_zip="gzip "+file_name
            os.system(cmd_zip)
            outfile1.close()
            outfile2.close()
            # output 3
            fastq_count=0
            seq_dict={}
            fastq_file=open(output2_folder+lib_name+'_3linker_removed.fastq','r')
            c=0
            for L2 in fastq_file:
                c+=1
                if c%4==2: 
                    fastq_count+=1
                    seq=L2.strip()
                    if seq in seq_dict: seq_dict[seq]+=1
                    else:seq_dict[seq]=1
            fastq_file.close()
            hist_count={}
            max_count=0
            outfile_fa=open(output3_folder+lib_name+'_3linker_removed_collapsed.fa','w')

            for s in seq_dict:
                outfile_fa.write('>'+lib_name+'_'+str(num_fasta)+'_count_'+str(seq_dict[s])+'\n')
                num_fasta+=1
                outfile_fa.write(s+'\n')
                if seq_dict[s]>max_count: max_count=seq_dict[s]
                if seq_dict[s] in hist_count: hist_count[seq_dict[s]]+=1
                else: hist_count[seq_dict[s]]=1
            outfile_fa.close()
            report_file=open(output3_folder+lib_name+'_collapsed_fasta_report.txt','w')

            report_file.write(lib_name+'\n')
            report_file.write('fastq reads: '+str(fastq_count)+'\n')
            print ('fastq reads: '+str(fastq_count))
            report_file.write('collapsed fasta reads: '+str(num_fasta)+'\n\n')
            print ('collapsed fasta reads: '+str(num_fasta))
            report_file.write('number of identical seq\tnumber of seq\n')
            for k in range(1,max_count+1):
                if k in hist_count:
                    report_file.write(str(k)+'\t'+str(hist_count[k])+'\n')
            report_file.close()
            outfile_report.write('%s\t%.2f\t%.2f\t%.2f\n' % (L1.strip(), num_fastq_raw/1000000, num_fastq_raw_with_barcodes/1000000, num_fasta/1000000))
            cmd_zip="gzip "+output1_folder+lib_name+'_raw.fastq'
            os.system(cmd_zip)
            cmd_zip="gzip "+output2_folder+lib_name+'_3linker_removed.fastq'
            os.system(cmd_zip)
    infile_lib_info.close()
    outfile_report.close()
    
def RNAseq_fasta2rawCounts_pdDataframe_geneWise(L): #L: a list of RNA-seq fa file names, ['SG0420_lib11", etc]
#9/15/2020, generate gene-wise raw counts for RNAseq data
# # input: a list fasta files
# reference: ce6_cDNA_v1
# return: raw counts in pandas dataframe
    gene_dict={'gene_name1':[],
               'gene_name2':[]} 
    genefile=open(HOME+'/Dropbox/bioinformatics/Genome_ce6/ce6_cDNA_v1.txt','r')
    for L1 in genefile:
        A1=(L1.strip()).split()
        (gene_name1,chrom,strand)=(A1[0],A1[1],A1[2])
        gene_name2=A1[10] if len(A1)==11 else gene_name1
        gene_dict['gene_name1'].append(gene_name1)
        gene_dict['gene_name2'].append(gene_name2)
    genefile.close()
    d=pd.DataFrame(gene_dict)
    d.set_index('gene_name1', inplace=True)
    for lib in L:
        print(lib)
        fa_name=lib+'_3linker_removed_collapsed.fa'
        fa_file=(dir_path(fa_name,'/mnt/da1/fastq_demultiplexed/'))+fa_name
        bowtie_map=fa_align_bowtie(fa_file,'ce6_cDNA_v1',4)
        count={}
        for g in gene_dict['gene_name1']:
            count[g]=0
        infile=open(bowtie_map,'r')
        for L in infile:
            A=L.strip().split('\t')
            (strand,g,e)=(A[1],A[2],int(A[-1]))
            if strand == '+':
                count[g]+=1/(1+e)
        infile.close()
        c_list=[]
        for g in gene_dict['gene_name1']:
            c_list.append(int(count[g]))
        d[lib]=c_list
    return d

def RNAseq_fastq2rawCounts_pdDataframe_geneWise(L): #L: a list of RNA-seq fa file name prefixes, ['SG0420_lib11", etc]
#11/10/2020, generate gene-wise raw counts for RNAseq data, 
# 
# # input: a list fastq files, f
# reference: ce6_cDNA_v1
# return: raw counts in pandas dataframe
# 
    #fq_name_appendex='_R1.fastq.gz' # this can be changed depending on the actual name of the fastq files
    #fq_name_appendex='_3link_rm.fastq'
    fq_name_appendex='_3linker_removed.fastq'
    #fq_name_appendex='.fastq.gz'
    nt_rm_5=4
    #nt_rm_3=97 #  use 4 and 97 here for Admera 151 nt sequencing run. This will remove 4 nt from 5'end and align the next 50 nt read 
    nt_rm_3=0 # for 3'liner linker removed RNA-seq
    gene_dict={'gene_name1':[],
               'gene_name2':[]} 
    genefile=open(HOME+'/Dropbox/bioinformatics/Genome_ce6/ce6_cDNA_v1.txt','r')
    for L1 in genefile:
        A1=(L1.strip()).split()
        (gene_name1,chrom,strand)=(A1[0],A1[1],A1[2])
        gene_name2=A1[10] if len(A1)==11 else gene_name1
        gene_dict['gene_name1'].append(gene_name1)
        gene_dict['gene_name2'].append(gene_name2)
    genefile.close()
    d=pd.DataFrame(gene_dict)
    d.set_index('gene_name1', inplace=True)
    for lib in L:
        print(lib)
        fq_name=lib+fq_name_appendex
        fq_file=(dir_path(fq_name,fastq_folder))+fq_name
        bowtie_map=fq_align_bowtie(fq_file,'ce6_cDNA_v1',nt_rm_5,nt_rm_3) 
        count={}
        for g in gene_dict['gene_name1']:
            count[g]=0
        infile=open(bowtie_map,'r')
        for L in infile:
            A=L.strip().split('\t')
            (strand,g,e)=(A[1],A[2],int(A[-1]))
            if strand == '+':
                count[g]+=1/(1+e)
        infile.close()
        c_list=[]
        for g in gene_dict['gene_name1']:
            c_list.append(int(count[g]))
        d[lib]=c_list
    return d

def RNAseq_fastq2rpkm_pdDataframe_geneWise(L): #L: a list of RNA-seq fa file name prefix, ['SG0420_lib11", etc]
#2/15/2021, generate gene-wise rpkm values for RNAseq data, 
# peudo counts for each gene: total ce6 align reads/20000*0.005
# # input: a list fastq files, f
# reference: ce6_cDNA_v1
# return: rpkm counts in pandas dataframe
# 
    #fq_name_appendex='.fastq.gz' # this can be changed depending on the actual name of the fastq files
    fq_name_appendex='_3linker_removed.fastq.gz'
    #fq_name_appendex='.fastq'
    nt_rm_5=4
    #nt_rm_3=97 #  use 4 and 97 here for Admera 151 nt sequencing run. This will remove 4 nt from 5'end and align the next 50 nt read 
    nt_rm_3=0 #  use 0 here for JH DNA core 50nt sequencing run.  
    gene_dict={'gene_name1':[],
               'gene_name2':[]} 
    genefile=open(HOME+'/Dropbox/bioinformatics/Genome_ce6/ce6_cDNA_v1.txt','r')
    cDNA_size=get_ce6_cDNA_size()
    for L1 in genefile:
        A1=(L1.strip()).split()
        (gene_name1,chrom,strand)=(A1[0],A1[1],A1[2])
        gene_name2=A1[10] if len(A1)==11 else gene_name1
        gene_dict['gene_name1'].append(gene_name1)
        gene_dict['gene_name2'].append(gene_name2)
    genefile.close()
    d=pd.DataFrame(gene_dict)
    d.set_index('gene_name1', inplace=True)
    for lib in L:
        print(lib)
        fq_name=lib+fq_name_appendex
        fq_file=(dir_path(fq_name,fastq_folder))+fq_name
        bowtie_map=fq_align_bowtie(fq_file,'ce6_cDNA_v1',nt_rm_5,nt_rm_3) 
        ce6_reads_N=get_total_aligned_number(fq_file,'ce6',nt_rm_5,nt_rm_3,'both')
        count={}
        for g in gene_dict['gene_name1']:
            count[g]=ce6_reads_N/20000*0.005
        infile=open(bowtie_map,'r')
        for L in infile:
            A=L.strip().split('\t')
            (strand,g,e)=(A[1],A[2],int(A[-1]))
            if strand == '+':
                count[g]+=1/(1+e)
        infile.close()
        c_list=[]
        for g in gene_dict['gene_name1']:
            c_list.append(int(count[g])/ce6_reads_N*1000000/cDNA_size[g]*1000)
        d[lib]=c_list
    return d

def sRNAseq_fasta2rawCounts_pdDataframe_geneWise(L): #L: a list of sNA-seq fa file names, ['SG0420_lib11", etc]
#9/15/2020, generate gene-wise raw counts for sRNAseq data
# # input: a list fasta files
# reference: ce6_cDNA_v1
# return: raw counts in pandas dataframe
    (sRNA_size_min,sRNA_size_max)=(20,24)
    gene_dict={'gene_name1':[],
               'gene_name2':[]} 
    genefile=open(HOME+'/Dropbox/bioinformatics/Genome_ce6/ce6_cDNA_v1.txt','r')
    for L1 in genefile:
        A1=(L1.strip()).split()
        (gene_name1,chrom,strand)=(A1[0],A1[1],A1[2])
        gene_name2=A1[10] if len(A1)==11 else gene_name1
        gene_dict['gene_name1'].append(gene_name1)
        gene_dict['gene_name2'].append(gene_name2)
    genefile.close()
    d=pd.DataFrame(gene_dict)
    d.set_index('gene_name1', inplace=True)
    for lib in L:
        print(lib)
        fa_name=lib+'_3linker_removed_collapsed.fa'
        fa_file=(dir_path(fa_name,fastq_folder))+fa_name
        bowtie_map=fa_align_bowtie(fa_file,'ce6_cDNA_v1',4,0)
        count={}
        for g in gene_dict['gene_name1']:
            count[g]=0
        infile=open(bowtie_map,'r')
        for L in infile:
            A=L.strip().split('\t')
            (strand,g,e,size)=(A[1],A[2],int(A[-1]),len(A[4]))
            if strand == '-' and size >= sRNA_size_min and size <= sRNA_size_max:
                count[g]+=1/(1+e)
        infile.close()
        c_list=[]
        for g in gene_dict['gene_name1']:
            c_list.append(int(count[g]))
        d[lib]=c_list
    return d

def sRNAseq_fasta2rawCounts_pdDataframe_geneWise_v2(L,ref): #L: a list of sNA-seq fa file names, ['SG0420_lib11", etc]
# ref: bowtie index
#9/15/2020, generate gene-wise raw counts for sRNAseq data
# # input: a list fasta files
# reference: ce6_cDNA_v1
# return: raw counts in pandas dataframe
    (sRNA_size_min,sRNA_size_max)=(20,24)
    gene_dict={'gene_name1':[],
               'gene_name2':[]} 
    genefile=open(HOME+'/Dropbox/bioinformatics/Genome_ce6/'+ref+'.txt','r')
    for L1 in genefile:
        A1=(L1.strip()).split()
        (gene_name1,chrom,strand)=(A1[0],A1[1],A1[2])
        gene_name2=A1[10] if len(A1)==11 else gene_name1
        gene_dict['gene_name1'].append(gene_name1)
        gene_dict['gene_name2'].append(gene_name2)
    genefile.close()
    d=pd.DataFrame(gene_dict)
    d.set_index('gene_name1', inplace=True)
    for lib in L:
        print(lib)
        fa_name=lib+'_3linker_removed_collapsed.fa'
        fa_file=(dir_path(fa_name,fastq_folder))+fa_name
        bowtie_map=fa_align_bowtie(fa_file,ref,4,0)
        count={}
        for g in gene_dict['gene_name1']:
            count[g]=0
        infile=open(bowtie_map,'r')
        for L in infile:
            A=L.strip().split('\t')
            (strand,g,e,size)=(A[1],A[2],int(A[-1]),len(A[4]))
            if strand == '-' and size >= sRNA_size_min and size <= sRNA_size_max:
                count[g]+=1/(1+e)
        infile.close()
        c_list=[]
        for g in gene_dict['gene_name1']:
            c_list.append(int(count[g]))
        d[lib]=c_list
    return d

def get_ce6_cDNA_size():
    # 2/15/2021, return a dict. Key: gene name, value: cDNA size
    genefile=open(HOME+'/Dropbox/bioinformatics/Genome_ce6/ce6_cDNA_v1.txt','r')
    cDNA_size={}
    for L1 in genefile:
        A1=(L1.strip()).split()
        (gene_name1,chrom,strand)=(A1[0],A1[1],A1[2])
        (exon_start,exon_end)= (A1[8][1:-1].split(','),A1[9][1:-1].split(','))  
        cDNA_size[gene_name1]=0
        for k in range(len(exon_start)-1):
            cDNA_size[gene_name1]+=int(exon_end[k])-int(exon_start[k])
    return cDNA_size

def get_ce6_cDNA_size_v2(ref): #ref: name of the bowtie index
    # 2/15/2021, return a dict. Key: gene name, value: cDNA size
    genefile=open(HOME+'/Dropbox/bioinformatics/Genome_ce6/'+ref+'.txt','r')
    cDNA_size={}
    for L1 in genefile:
        A1=(L1.strip()).split()
        (gene_name1,chrom,strand)=(A1[0],A1[1],A1[2])
        (exon_start,exon_end)= (A1[8][1:-1].split(','),A1[9][1:-1].split(','))  
        cDNA_size[gene_name1]=0
        for k in range(len(exon_start)-1):
            cDNA_size[gene_name1]+=int(exon_end[k])-int(exon_start[k])
    return cDNA_size

def sRNAseq_fasta2rpkm_pdDataframe_geneWise(L): 
#L: a list of sNA-seq fa file prefix, ['SG0420_lib11", etc]
# 02/15/21, generate gene-wise rpkm values for sRNA-seq data
# # input: a list fasta file name prefix: ['SG0420_lib11", etc]
# reference: ce6_cDNA_v1
# return: rpkm in pandas dataframe
    (sRNA_size_min,sRNA_size_max)=(20,24)
    gene_dict={'gene_name1':[],
               'gene_name2':[]} 
    genefile=open(HOME+'/Dropbox/bioinformatics/Genome_ce6/ce6_cDNA_v1.txt','r')
    cDNA_size=get_ce6_cDNA_size()
    for L1 in genefile:
        A1=(L1.strip()).split()
        (gene_name1,chrom,strand)=(A1[0],A1[1],A1[2])
        gene_name2=A1[10] if len(A1)==11 else gene_name1
        gene_dict['gene_name1'].append(gene_name1)
        gene_dict['gene_name2'].append(gene_name2)
    genefile.close()
    d=pd.DataFrame(gene_dict)
    d.set_index('gene_name1', inplace=True)
    for lib in L:
        print(lib)
        fa_name=lib+'_3linker_removed_collapsed.fa'
        fa_file=(dir_path(fa_name,fastq_folder))+fa_name
        ce6_reads_N=get_total_aligned_number(fa_file,'ce6',4,0,'both')
        bowtie_map=fa_align_bowtie(fa_file,'ce6_cDNA_v1',4,0)
        count={}
        for g in gene_dict['gene_name1']:
            count[g]=0
        infile=open(bowtie_map,'r')
        for L in infile:
            A=L.strip().split('\t')
            (strand,g,e,size)=(A[1],A[2],int(A[-1]),len(A[4]))
            if strand == '-' and size >= sRNA_size_min and size <= sRNA_size_max:
                count[g]+=1/(1+e)
        infile.close()
        c_list=[]
        for g in gene_dict['gene_name1']:
            c_list.append(int(count[g])/ce6_reads_N*1000000/cDNA_size[g]*1000)
        d[lib]=c_list
    return d

def sRNAseq_fasta2rpkm_pdDataframe_geneWise_v2(L, ref): 
#L: a list of sNA-seq fa file prefix, ['SG0420_lib11", etc]
#ref: bowtie reference (4/28/20 modification)
# 02/15/21, generate gene-wise rpkm values for sRNA-seq data
# # input: a list fasta file name prefix: ['SG0420_lib11", etc]
# reference: ce6_cDNA_v1
# return: rpkm in pandas dataframe
    (sRNA_size_min,sRNA_size_max)=(20,24)
    gene_dict={'gene_name1':[],
               'gene_name2':[]} 
    genefile=open(HOME+'/Dropbox/bioinformatics/Genome_ce6/'+ref+'.txt','r')
    cDNA_size=get_ce6_cDNA_size_v2(ref)
    for L1 in genefile:
        A1=(L1.strip()).split()
        (gene_name1,chrom,strand)=(A1[0],A1[1],A1[2])
        gene_name2=A1[10] if len(A1)==11 else gene_name1
        gene_dict['gene_name1'].append(gene_name1)
        gene_dict['gene_name2'].append(gene_name2)
    genefile.close()
    d=pd.DataFrame(gene_dict)
    d.set_index('gene_name1', inplace=True)
    for lib in L:
        print(lib)
        fa_name=lib+'_3linker_removed_collapsed.fa'
        fa_file=(dir_path(fa_name,fastq_folder))+fa_name
        ce6_reads_N=get_total_aligned_number(fa_file,'ce6',4,0,'both')
        bowtie_map=fa_align_bowtie(fa_file,ref,4,0)
        count={}
        for g in gene_dict['gene_name1']:
            count[g]=0
        infile=open(bowtie_map,'r')
        for L in infile:
            A=L.strip().split('\t')
            (strand,g,e,size)=(A[1],A[2],int(A[-1]),len(A[4]))
            if strand == '-' and size >= sRNA_size_min and size <= sRNA_size_max:
                count[g]+=1/(1+e)
        infile.close()
        c_list=[]
        for g in gene_dict['gene_name1']:
            c_list.append(int(count[g])/ce6_reads_N*1000000/cDNA_size[g]*1000)
        d[lib]=c_list
    return d

def fq_align_bowtie(fq_file,bowtie_ref,remove_5_size,remove_3_size):
    bowtie_map_file=map_folder+fq_file.split('/')[-1]+'_'+bowtie_ref+'.map'
    bowtie_prog=bowtie_folder+'bowtie '+bowtie_ref_folder+bowtie_ref+' '+fq_file+' '+bowtie_map_file+' -5 '+str(remove_5_size)+' -3 '+str(remove_3_size)+' -q -t -v 0 -a -p '+str(core_num)
#    print(bowtie_prog)
    if os.path.isfile(bowtie_map_file):
        print ('%s exist.'% bowtie_map_file)
        return bowtie_map_file
    else:
        p=subprocess.run(bowtie_prog, shell=True, check=True,capture_output=True)
        print(p.stderr)
        return bowtie_map_file

def fq_align_bowtie_SAM(fq_file,bowtie_ref,remove_5_size,remove_3_size): # return SAM file
    f_type=' -q ' if fq_file.find('.fastq')!=-1 else ' -f '
    bowtie_map_file=map_folder+fq_file.split('/')[-1]+'_'+bowtie_ref+'.sam'
    bowtie_prog=bowtie_folder+'bowtie '+bowtie_ref_folder+bowtie_ref+' '+fq_file+' -S '+bowtie_map_file+' -5 '+str(remove_5_size)+' -3 '+str(remove_3_size)+f_type+' -t -v 0 -a -p '+str(core_num)
    print(bowtie_prog)
    if os.path.isfile(bowtie_map_file):
        print ('%s exist.'% bowtie_map_file)
        return bowtie_map_file
    else:
        p=subprocess.run(bowtie_prog, shell=True, check=True,capture_output=True)
        print(p.stderr)
        return bowtie_map_file

def fq_align_bowtie2_SAM(fq_file_R1, fq_file_R2,lib_name,bowtie_ref,remove_5_size,remove_3_size):
    bowtie_map_file=map_folder+lib_name+'_'+bowtie_ref+'_paired_end.sam'
#    bowtie_prog=bowtie2_folder+'bowtie2 -x '+bowtie2_ref_folder+bowtie_ref+' -U '+fq_file+' -S '+bowtie_map_file+' -5 '+str(remove_5_size)+' -3 '+str(remove_3_size)+' -t -p '+str(core_num)
    bowtie_prog=bowtie2_folder+'bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -x '+bowtie2_ref_folder+bowtie_ref+' -1 '+fq_file_R1+' -2 '+fq_file_R2+' -S '+bowtie_map_file+' -5 '+str(remove_5_size)+' -3 '+str(remove_3_size)+' -t -p '+str(core_num)
#    print(bowtie_prog)
    if os.path.isfile(bowtie_map_file):
        print ('%s exist.'% bowtie_map_file)
        return bowtie_map_file
    else:
        p=subprocess.run(bowtie_prog, shell=True, check=True,capture_output=True)
        print(p.stderr)
        return bowtie_map_file
    
def ChIPseq_fastq2rpm_geneWise(L): #L: a list of ChIPseq fa file names, ['SG0420_lib11", etc]
# Stop using this one because it uses cDNA as bowtie reference
# 9/15/2020, generate gene-wise rpm values for ChIPseq data
# # input: a list fastq R1 files
# reference: ce6_cDNA_v1
# return: rpm in pandas dataframe
    #fq_name_extension='.fastq.gz'
    #remove_3end_size=0 
    fq_name_extension='_R1.fastq.gz'
    remove_3end_size=101
    gene_dict={'gene_name1':[],
               'gene_name2':[]} 
    genefile=open(HOME+'/Dropbox/bioinformatics/Genome_ce6/ce6_cDNA_v1.txt','r')
    for L1 in genefile:
        A1=(L1.strip()).split()
        (gene_name1,chrom,strand)=(A1[0],A1[1],A1[2])
        gene_name2=A1[10] if len(A1)==11 else gene_name1
        gene_dict['gene_name1'].append(gene_name1)
        gene_dict['gene_name2'].append(gene_name2)
    genefile.close()
    d=pd.DataFrame(gene_dict)
    d.set_index('gene_name1', inplace=True)
    for lib in L:
        print(lib)
        fq_name=lib+fq_name_extension
        fq_file=(dir_path(fq_name,fastq_folder))+fq_name
        bowtie_map=fq_align_bowtie(fq_file,'ce6_cDNA_v1',0,0)
        count={}
        for g in gene_dict['gene_name1']:
            count[g]=0
        infile=open(bowtie_map,'r')
        tc=0 # total reads aligned to ce6_cDNA
        for L in infile:
            A=L.strip().split('\t')
            (strand,g,e,size)=(A[1],A[2],int(A[-1]),len(A[4]))
            count[g]+=1/(1+e)
            tc+=1/(1+e)
        infile.close()
        c_list=[]
        for g in gene_dict['gene_name1']:
            c_list.append(round(count[g]/tc*1000000,3))
        d[lib]=c_list
    return d

def ChIPseq_fastq2rpkm_geneWise(L): #L: a list of ChIPseq fa file names, ['SG0420_lib11", etc]
# 6/29/2021 calculate gene by gene rpkm for ChIP-seq libs
# input: a list fastq R1 files
# reference: ce6
# return: rpkm in pandas dataframe
    gene_extension=200 # flanking sequence size to be considered for each each, 
    #fq_name_extension='.fastq.gz'
    remove_3end_size=0 
    fq_name_extension='_R1.fastq.gz'
    #remove_3end_size=101
    gene_dict={'gene_name1':[],
               'gene_name2':[]} 
    gene_info={}
    genefile=open(HOME+'/Dropbox/bioinformatics/Genome_ce6/ce6_cDNA_v1.txt','r')
    for L1 in genefile:
        A1=(L1.strip()).split()
        (gene_name1,chrom,strand,txStart,txEnd)=(A1[0],A1[1],A1[2],int(A1[3]),int(A1[4]) )
        gene_name2=A1[10] if len(A1)==11 else gene_name1
        gene_dict['gene_name1'].append(gene_name1)
        gene_dict['gene_name2'].append(gene_name2)
        gene_info[gene_name1]={'chrom':chrom,
                               'txStart':txStart,
                               'txEnd':txEnd}
    genefile.close()
    d=pd.DataFrame(gene_dict)
    d.set_index('gene_name1', inplace=True)
    for lib in L:
        print(lib)
        fq_name=lib+fq_name_extension
        fq_file=(dir_path(fq_name,fastq_folder))+fq_name
        bowtie_map=fq_align_bowtie(fq_file,'ce6',0,remove_3end_size)
        count={}
        for g in gene_dict['gene_name1']:
            count[g]=0
        infile=open(bowtie_map,'r')
        tc=0 # total reads aligned to ce6
        
        cov={}
        chr_size = {'chrI':15072421,'chrII':15279323,'chrIII':13783681,'chrIV':17493785,'chrV':20919568,'chrX':17718854}
        for chrom in chr_size:
            cov[chrom]={}
            for k in range(-5000,chr_size[chrom]+5000):
                cov[chrom][k]=0

        for L in infile:
            A1=L.strip().split('\t')
            (chromosome,position,nonUnique)=(A1[2],int(A1[3]),int(A1[6]))
            cov[chromosome][position]+=1.0/(1+nonUnique)
            tc+=1/(1+nonUnique)
        infile.close()
        
        c_list=[]
        for g in gene_dict['gene_name1']:
            cov_g=0
            for p in range(gene_info[g]['txStart']-gene_extension,gene_info[g]['txEnd']+gene_extension):
                cov_g+=cov[gene_info[g]['chrom']][p]
            gene_size=gene_info[g]['txEnd']-gene_info[g]['txStart']+2*gene_extension
            c_list.append(round(cov_g/gene_size/tc*1000000*1000,3))
        d[lib]=c_list
    return d

def ChIPseq_fastq2rawCounts_geneWise(L): #L: a list of ChIPseq fa file names, ['SG0420_lib11", etc]
# 6/29/2021 calculate gene by gene raw counts for ChIP-seq libs
# input: a list fastq R1 files
# reference: ce6
# return: raw counts per gene in pandas dataframe
    gene_extension=200 # flanking sequence size to be considered for each each, 
    #fq_name_extension='.fastq.gz'
    #remove_3end_size=0 
    fq_name_extension='_R1.fastq.gz'
    remove_3end_size=101
    gene_dict={'gene_name1':[],
               'gene_name2':[]} 
    gene_info={}
    genefile=open(HOME+'/Dropbox/bioinformatics/Genome_ce6/ce6_cDNA_v1.txt','r')
    for L1 in genefile:
        A1=(L1.strip()).split()
        (gene_name1,chrom,strand,txStart,txEnd)=(A1[0],A1[1],A1[2],int(A1[3]),int(A1[4]) )
        gene_name2=A1[10] if len(A1)==11 else gene_name1
        gene_dict['gene_name1'].append(gene_name1)
        gene_dict['gene_name2'].append(gene_name2)
        gene_info[gene_name1]={'chrom':chrom,
                               'txStart':txStart,
                               'txEnd':txEnd}
    genefile.close()
    d=pd.DataFrame(gene_dict)
    d.set_index('gene_name1', inplace=True)
    for lib in L:
        print(lib)
        fq_name=lib+fq_name_extension
        fq_file=(dir_path(fq_name,fastq_folder))+fq_name
        bowtie_map=fq_align_bowtie(fq_file,'ce6',0,remove_3end_size)
        count={}
        for g in gene_dict['gene_name1']:
            count[g]=0
        infile=open(bowtie_map,'r')

        
        cov={}
        chr_size = {'chrI':15072421,'chrII':15279323,'chrIII':13783681,'chrIV':17493785,'chrV':20919568,'chrX':17718854}
        for chrom in chr_size:
            cov[chrom]={}
            for k in range(-5000,chr_size[chrom]+5000):
                cov[chrom][k]=0

        for L in infile:
            A1=L.strip().split('\t')
            (chromosome,position,nonUnique)=(A1[2],int(A1[3]),int(A1[6]))
            cov[chromosome][position]+=1.0/(1+nonUnique)

        infile.close()
        
        c_list=[]
        for g in gene_dict['gene_name1']:
            cov_g=0
            for p in range(gene_info[g]['txStart']-gene_extension,gene_info[g]['txEnd']+gene_extension):
                cov_g+=cov[gene_info[g]['chrom']][p]
            gene_size=gene_info[g]['txEnd']-gene_info[g]['txStart']+2*gene_extension
            c_list.append(int(cov_g))
        d[lib]=c_list
    return d

def cp_rename_chipLib(f): 
# 9/24/2020 copy admera chipseq fastq.gz files to demultiplex folder, change name to SG number
# also can be used for RNAseq libs that are made by Smarter kit and sequenced by admera 
# needs the ChIP_lib_info_admera_SG###.csv file (f)
                     
    infile_lib_info=open(f,'r')
    e=0
    for L1 in infile_lib_info:
        e+=1
        if e==1:
            p=re.compile(',')
            m=p.search(L1)
            sep='\t' if m == None else ','
        else:
            A=L1.strip().split(sep)
            (sampleID,libName,s_num,rawFolder,demulFolder,L_num)=(A[0],A[1],A[2], A[3],A[4],A[5])
            #fn1="%s%sV1_%s_%s_R1_001.fastq.gz" % (rawFolder,sampleID,s_num,L_num)
            fn1="%s%s_%s_%s_R1_001.fastq.gz" % (rawFolder,sampleID,s_num,L_num)
            fn2="%s%s_R1.fastq.gz" % (demulFolder,libName)
            cmd='cp %s %s' % (fn1,fn2)
            subprocess.run(cmd, shell=True)
            #fn1="%s%sV1_%s_%s_R2_001.fastq.gz" % (rawFolder,sampleID,s_num,L_num)
            fn1="%s%s_%s_%s_R2_001.fastq.gz" % (rawFolder,sampleID,s_num,L_num)
            fn2="%s%s_R2.fastq.gz" % (demulFolder,libName)
            cmd='cp %s %s' % (fn1,fn2)
            print (cmd)
            subprocess.run(cmd, shell=True)
    infile_lib_info.close()

def whole_chrom_coverage_ChIP_plot1(name1,name2):
# 9/25/2020 show whole chr. coverages for two chip-seq libs
# (name1,name2)=('SG0820_lib1','SG0820_lib5')
    (f_n1,f_n2)=(name1+'_R1.fastq.gz',name2+'_R1.fastq.gz')
    (lib1,lib2)=(dir_path(f_n1,fastq_folder)+f_n1,dir_path(f_n2,fastq_folder)+f_n2 )

    (map_file1,map_file2)=(fq_align_bowtie(lib1,'ce6',0,101),fq_align_bowtie(lib2,'ce6',0,101))

    chr_size_chrI_rRNAremoved = {'chrI':15059999,'chrII':15279323,'chrIII':13783681,'chrIV':17493785,'chrV':20919568,'chrX':17718854} 
    chr_size = {'chrI':15072421,'chrII':15279323,'chrIII':13783681,'chrIV':17493785,'chrV':20919568,'chrX':17718854}
    chr_list=['chrI','chrII','chrIII','chrIV','chrV','chrX']
    genome_size=101169 # kb

    (cov1,cov2)=({},{}) #raw counts of perfect alignments of unique and nonunique alignments
    (genome_cov1,genome_cov2)=(0,0)
    for chrom in chr_size:
        (cov1[chrom],cov2[chrom])=({},{})
        for k in range(int(chr_size[chrom]/10000)+1):
            (cov1[chrom][k],cov2[chrom][k])=(0,0)
    infile1=open(map_file1,'r')
    for L1 in infile1:
        A1=L1.split('\t')
        (chromosome,position,strand,nonUnique)=(A1[2],int(A1[3]),A1[1],int(A1[6]))
        genome_cov1+=1.0/(1+nonUnique)
        cov1[chromosome][int(position/10000)]+=1/(1.0+nonUnique)
    infile1.close()

    infile1=open(map_file2,'r')
    for L1 in infile1:
        A1=L1.split('\t')
        (chromosome,position,strand,nonUnique)=(A1[2],int(A1[3]),A1[1],int(A1[6]))
        genome_cov2+=1.0/(1+nonUnique)
        cov2[chromosome][int(position/10000)]+=1/(1.0+nonUnique)
    infile1.close()

    fig, axs = plt.subplots(nrows=6, ncols=1, squeeze=False,figsize=(15,15))
    for r in range(6):
        chrom=chr_list[r]
        cov_d1=pd.DataFrame(cov1[chrom].items(), columns=['pos','cov']).sort_values('pos')
        cov_d2=pd.DataFrame(cov2[chrom].items(), columns=['pos','cov']).sort_values('pos')
        axs[r,0].plot(cov_d1['pos']*10, cov_d1['cov']/genome_cov1*(genome_size/10), label=name1,alpha=0.5)
        axs[r,0].plot(cov_d2['pos']*10, cov_d2['cov']/genome_cov2*(genome_size/10), label=name2,alpha=0.5)
        axs[r,0].set(xlabel='%s (kb)'%chrom, ylabel='coverage', ylim=[-1,10],title=' ')
        axs[r,0].legend(loc="upper left")
        fig.suptitle('Whole chromosome coverage',y=0.9)
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.1, hspace=0.4)
    plt.savefig("whole_chr_coverage_%s_%s.png"% (name1,name2),dpi=72)
    
def get_total_aligned_number(fq_file,bowtie_ref,remove_5_size,remove_3_size,strand):
#10/6/2020 return the total number of bowtie aligned reads
# can handle both fa and fastq files, zip or unziped
# strand: 
# '+': report only + strand alignment number
# '-': report only - strand alignment number
# 'both': report a alignment number regardless strand
# 'sep': report a dictionary of {'+':xxx, '-':xxx}
    fq=' -q ' if fq_file.find('.fastq')!=-1 else ' -f '
    bowtie_map_file=map_folder+'temp.map'
    bowtie_prog=bowtie_folder+'bowtie '+bowtie_ref_folder+bowtie_ref+' '+fq_file+' '+bowtie_map_file+' -5 '+str(remove_5_size)+' -3 '+str(remove_3_size)+ fq +' -t -v 0 -a -p '+str(core_num)
    p=subprocess.run(bowtie_prog, shell=True, check=True,capture_output=True)
    output=p.stderr
    print(output)
    if strand=='both':
        rm_comd='rm '+bowtie_map_file
        p=subprocess.run(rm_comd, shell=True, check=True,capture_output=True)
        return int((output.split(b'\n')[3]).split(b' ')[8])
    else: 
        (c_plus,c_minus)=(0,0)
        FH=open(bowtie_map_file,'r')
        for L1 in FH:
            if L1.split('\t')[1]=='+':
                c_plus+=1/(1+int(L1.strip().split('\t')[-1]))
            else:
                c_minus+=1/(1+int(L1.strip().split('\t')[-1]))
        FH.close()
        rm_comd='rm '+bowtie_map_file
        p=subprocess.run(rm_comd, shell=True, check=True,capture_output=True)
        if strand=='+':
            return int(c_plus)
        if strand=='-':
            return int(c_minus)
        if strand=='sep':
            return {'+':int(c_plus), '-':int(c_minus)}
    
def total_piRNA_count(file_name,barcode_size):
    # count total piRNA reads in a fasta or fastq file
    #  return {'total_reads':c2/f_type,'piRNA_reads':tpc}
    # a bit slow, but more precise than bowtie align, bowtie will only give ~ 50% of what's given from this function
    tpc=0 # total piRNA count
    piRNA_list=[] # store all annotated piRNAs in ce10, grouped in 200-piRNA sub list 
    FH_piRNA=open('/home/sam/Dropbox/bioinformatics/Genome_ce10/piRNA_data_c_elegans.csv','r')
    c=0
    for L1 in FH_piRNA:
        c+=1
        if c>1:
            piRNA_list.append(L1.strip().split(',')[1].replace('U','T'))
    FH_piRNA.close()
    piRNA_set=set(piRNA_list)
    # determine if the file is .gz or not
    gz=True if file_name.find('.gz')!=-1 else False
    if gz:
        FH_reads=gzip.open(file_name,'r')
    else:
        FH_reads=open(file_name,'r')
    f_type=4 if file_name.find('.fastq')!=-1 else 2
    c2=0
    for L1 in FH_reads:
        if gz:
            L1=L1.decode()
        c2+=1
        if c2%f_type==(f_type-2):
            read=L1.strip()[barcode_size:len(L1.strip())]
            if read in piRNA_set:
                tpc+=1
    FH_reads.close()
    return {'total_reads':c2/f_type,'piRNA_reads':tpc}

def sRNA_diff_genomic_feature_counts(l, barcode_size):

# l: lib name eg "SG0321_lib7"
#return a dictionary such as:
#{'ce6_count': 19753094,
#  'miRNA_count': 98724, # miRNA counts / ce6 aligned x 100
#  'cDNA_sense_count': 515551,
#  'cDNA_antisense_count': 9247910,
#  'rRNA_sense_count': 2461346,
#  'rRNA_antisense_count': 7662,
#  'total_reads_count': 25052500.0,
#  'piRNA_count': 26058,
#  'miRNA_percentage': 0.49979005820556516,
#  'cDNA_sense_percentage': 2.6099759359217347,
#  'cDNA_antisense_percentage': 46.81752640877424,
#  'rRNA_sense_percentage': 12.460559343260352,
#  'rRNA_antisense_percentage': 0.038788860114774934,
#  'piRNA_pecentage': 0.13191857437624707,
#  'ce6_percentage': 78.84679772477796}  # ce6 aligned / total reads x 100
#    '''
#    fa_extension='_3linker_removed_collapsed.fa'
    fa_extension='_3linker_removed.fastq.gz'
    fa_folder=dir_path(l+fa_extension,fastq_folder)
    fp=fa_folder+l+fa_extension
    
    c_dict={}
    ce6_count=get_total_aligned_number(fp,'ce6',barcode_size,0,'both')
    c_dict['ce6_count']=ce6_count
    miRNA_count=get_total_aligned_number(fp,'ce_miR_hairpin_mirBase21',barcode_size,0,'+')
    c_dict['miRNA_count']=miRNA_count
    cDNA=get_total_aligned_number(fp,'ce6_cDNA_v1',barcode_size,0,'sep')
    c_dict['cDNA_sense_count']=cDNA['+']
    c_dict['cDNA_antisense_count']=cDNA['-']
    rRNA=get_total_aligned_number(fp,'ce_rRNA',barcode_size,0,'sep')
    c_dict['rRNA_sense_count']=rRNA['+']
    c_dict['rRNA_antisense_count']=rRNA['-']
    piRNA=total_piRNA_count(fp,barcode_size)
    c_dict['total_reads_count']=piRNA['total_reads']
    c_dict['piRNA_count']=piRNA['piRNA_reads']
    (c_dict['miRNA_percentage'],c_dict['cDNA_sense_percentage'],c_dict['cDNA_antisense_percentage'])=(c_dict['miRNA_count']/c_dict['ce6_count']*100,c_dict['cDNA_sense_count']/c_dict['ce6_count']*100,c_dict['cDNA_antisense_count']/c_dict['ce6_count']*100)
    (c_dict['rRNA_sense_percentage'],c_dict['rRNA_antisense_percentage'])=(c_dict['rRNA_sense_count']/c_dict['ce6_count']*100,c_dict['rRNA_antisense_count']/c_dict['ce6_count']*100)
    (c_dict['piRNA_pecentage'],c_dict['ce6_percentage'])=(c_dict['piRNA_count']/c_dict['ce6_count']*100, c_dict['ce6_count']/piRNA['total_reads']*100)
    return c_dict

def sRNAseq2rawCounts_pdDataframe_geneWise(D):
#11/19/2020, generate gene-wise raw counts for sRNAseq data
#input D: a dictionary of sRNA-seq file info, eg {'SG0420_lib11":{'file_appendix':'_3linker_removed_collapsed.fa','barcode_size':4}, fastq, fastq.gz file appendix 
# reference: ce6_cDNA_v1
# return: raw counts in pandas dataframe
    (sRNA_size_min,sRNA_size_max)=(20,24)
    gene_dict={'gene_name1':[],
               'gene_name2':[]} 
    genefile=open(HOME+'/Dropbox/bioinformatics/Genome_ce6/ce6_cDNA_v1.txt','r')
    for L1 in genefile:
        A1=(L1.strip()).split()
        (gene_name1,chrom,strand)=(A1[0],A1[1],A1[2])
        gene_name2=A1[10] if len(A1)==11 else gene_name1
        gene_dict['gene_name1'].append(gene_name1)
        gene_dict['gene_name2'].append(gene_name2)
    genefile.close()
    d=pd.DataFrame(gene_dict)
    d.set_index('gene_name1', inplace=True)
    for lib in D:
        print(lib)
        file_name=dir_path(lib+D[lib]['file_appendix'],'/mnt/da1/fastq_demultiplexed/')+lib+D[lib]['file_appendix']
        bowtie_map=fa_fq_align_bowtie(file_name,'ce6_cDNA_v1',D[lib]['barcode_size'],0)
        count={}
        for g in gene_dict['gene_name1']:
            count[g]=0
        infile=open(bowtie_map,'r')
        for L in infile:
            A=L.strip().split('\t')
            (strand,g,e,size)=(A[1],A[2],int(A[-1]),len(A[4]))
            if strand == '-' and size >= sRNA_size_min and size <= sRNA_size_max:
                count[g]+=1/(1+e)
        infile.close()
        c_list=[]
        for g in gene_dict['gene_name1']:
            c_list.append(int(count[g]))
        d[lib]=c_list
    return d

def fa_fq_align_bowtie(fq_file,bowtie_ref,remove_5_size,remove_3_size):
    # 11/19/2020 This function can align both fastq or fasta/fa file
    f_type=' -q ' if fq_file.find('.fastq')!=-1 else ' -f '
    bowtie_map_file=map_folder+fq_file.split('/')[-1]+'_'+bowtie_ref+'.map'
    bowtie_prog=bowtie_folder+'bowtie '+bowtie_ref_folder+bowtie_ref+' '+fq_file+' '+bowtie_map_file+' -5 '+str(remove_5_size)+' -3 '+str(remove_3_size)+f_type+'-t -v 0 -a -p '+str(core_num)
#    print(bowtie_prog)
    if os.path.isfile(bowtie_map_file):
        print ('%s exist.'% bowtie_map_file)
        return bowtie_map_file
    else:
        p=subprocess.run(bowtie_prog, shell=True, check=True,capture_output=True)
        print(p.stderr)
        return bowtie_map_file

def get_gene_info():
    # 11/20/2020
    # generate a pd dataframe from ce6_cDNA_v1 with the same order of the gene list
    # annotate germline_gene_expresssion, mRNA_size, GRTS, and GRH
    gene_dict={'gene_name1':[],
               'gene_name2':[],
               'mRNA_size':[],
               'germline_expression':[],
               'GRTS':[],
               'GRH':[]} 
    genefile=open(HOME+'/Dropbox/bioinformatics/Genome_ce6/ce6_cDNA_v1.txt','r')
    for L1 in genefile:
        A1=(L1.strip()).split()
        (gene_name1,chrom,strand,tx_start,tx_end)=(A1[0],A1[1],A1[2],int(A1[3]),int(A1[4]))
        gene_name2=A1[10] if len(A1)==11 else gene_name1
        gene_dict['gene_name1'].append(gene_name1)
        gene_dict['gene_name2'].append(gene_name2)
        gene_dict['mRNA_size'].append(tx_end-tx_start)
        gene_dict['germline_expression'].append('no')
        gene_dict['GRTS'].append('no')
        gene_dict['GRH'].append('no')
    genefile.close()
    # annotate germline expressed genes
    # using Kimble 2014 paper
    germline_gene_list_FN='%s/Dropbox/bioinformatics/Genome_ce6/germline_expressed_gene_kimble_2014a.txt'%HOME
    germ_FH=open(germline_gene_list_FN,'r')
    c=0
    germ_gene_dict={'Gender Neutral':[],
                    'Spermatogenic':[],
                    'Oogenic':[]}
    for L1 in germ_FH:
        c+=1
        if c>1:
            A1=L1.strip().split('\t')
            (gn1, gn2, g_type)=(A1[0],A1[1],A1[7])
            germ_gene_dict[g_type].append(gn1)
            germ_gene_dict[g_type].append(gn2)
    germ_FH.close()
    # annotate GRTS or GRH genes        
    # use Julie Ni's 2014 paper
    (GRTS_gene,GRH_gene)=([],[])
    GRTS_file=open('%s/Dropbox/bioinformatics/Genome_ce6/GRTS_GRH/GRTS_gene_ZN2014.txt'%HOME,'r') 
    GRH_file=open('%s/Dropbox/bioinformatics/Genome_ce6/GRTS_GRH/GRH_gene_ZN2014.txt'%HOME,'r')
    for g in GRTS_file:
        GRTS_gene.append(g.strip())
    for g in GRH_file:
        GRH_gene.append(g.strip())    
    GRTS_file.close()
    GRH_file.close()
    for k in range(len(gene_dict['gene_name2'])):
        g=gene_dict['gene_name2'][k]
        for g_type in germ_gene_dict:
            if g in germ_gene_dict[g_type] :
                gene_dict['germline_expression'][k] = g_type           
                break
        if g in  GRTS_gene:
                gene_dict['GRTS'][k] = 'yes'           
        if g in  GRH_gene:
                gene_dict['GRH'][k] = 'yes'          
    d=pd.DataFrame(gene_dict)
    d.set_index('gene_name1',inplace=True)
    return d

def pval_beyes(T1,T2,C1,C2):
    # calculate the pval using the Bayesian method
    # T1, T2: total C. elegans genome aligned reads
    # C1, C2, unnormalized counts per gene in the form of pandas dataframe column
    PC=(C1+C2)/(T1+T2)
    # hypothesis 1:  two samples have the same probability for yielding a read for a given gene 
    # PROB1st = [PC^P1] * [PC^P2] * [(1-PC)^(T1-P1)] * [(1-PC)^(T2-P2)].  
    PROB1st_log=(C1+C2)*np.log(PC)+(T1+T2-(C1+C2))*np.log(1-PC)
    # hypothesis 2:  each sample has a unique probability for yielding a read for a give gene.
    #probablity for individual samples:
    (PC1,PC2)=(C1/T1,C2/T2)
    # PROB2nd = [PC1^P1] * [PC2^P2] * [(1-PC1)^(T1-P1)] *[(1-PC2)^(T2-P2)]
    PROB2nd_log=C1*np.log(PC1)+C2*np.log(PC2)+(T1-C1)*np.log(1-PC1)+(T2-C2)*np.log(1-PC2)
    return(np.exp(PROB1st_log-PROB2nd_log))

def deseq2_sRNA_simple1(lib_list,exp_design,exp_name,working_folder):
# updated on 2/5/2021
# generate a R scatter plot and the deseq2 csv file with fold changes and p values
# user input example
# lib_list=['SG0121_lib%s' % i for i in range(15,19)]
# exp_design={'fem1':['SG0121_lib15','SG0121_lib16'],
#             'glp1':['SG0121_lib17','SG0121_lib18'],}
# exp_name='fem1_glp1_sRNA_013021'
# working_folder='/home/sam/Dropbox/sandbox/2021/202101/013021_fem1vsglp1_sRNA'

    r_fn='deseq2_script_%s.R' % exp_name
    raw_count_csv='rawcounts_%s.csv' % exp_name
    deseq2_result_csv='deseq2_result_%s.csv'  % exp_name
    exp_csv='exp_%s.csv'  % exp_name

    d=sRNAseq_fasta2rawCounts_pdDataframe_geneWise(lib_list)
    d[lib_list].to_csv(raw_count_csv, sep=',')

    # write exp design file
    exp_file=open(exp_csv,'w')
    exp_file.write('lib_name,condition\n')
    for c in exp_design:
        for lib in exp_design[c]:
            exp_file.write('%s,%s\n' %(lib,c))
    exp_file.close()

    # write R script file
    r_file=open(r_fn,'w')
    r_file.write("library('DESeq2')\n")
    r_file.write("library(ggplot2)\n")
    r_file.write("setwd('%s')\n" % working_folder)
    r_file.write("countData <- read.csv('%s',header=TRUE,sep=',')\n" % raw_count_csv)
    r_file.write("metaData <- read.csv('%s',sep=',')\n" % exp_csv)
    r_file.write("dds<-DESeqDataSetFromMatrix(countData=countData, colData=metaData, design=~condition, tidy=TRUE)\n")
    r_file.write("dds <- DESeq(dds)\n")
    r_file.write("res <- results(dds)\n")
    r_file.write("summary(res)\n")
    r_file.write("plotMA(res, ylim=c(-2,2),cex=0.8)\n")
    r_file.write("write.csv(as.data.frame(res),file='%s')\n " % deseq2_result_csv)
    r_file.close()
    r_scriptrun='Rscript %s' % r_fn
    p=subprocess.run(r_scriptrun, shell=True, check=True,capture_output=True)
  
def deseq2_RNA_simple1(lib_list,exp_design,exp_name,working_folder):
# updated on 2/5/2021
# generate a R scatter plot and the deseq2 csv file with fold changes and p values
# user input example
# lib_list=['SG0121_lib%s' % i for i in range(15,19)]
# exp_design={'fem1':['SG0121_lib15','SG0121_lib16'],
#             'glp1':['SG0121_lib17','SG0121_lib18'],}
# exp_name='fem1_glp1_sRNA_013021'
# working_folder='/home/sam/Dropbox/sandbox/2021/202101/013021_fem1vsglp1_sRNA'

    r_fn='deseq2_script_%s.R' % exp_name
    raw_count_csv='rawcounts_%s.csv' % exp_name
    deseq2_result_csv='deseq2_result_%s.csv'  % exp_name
    exp_csv='exp_%s.csv'  % exp_name

    d=RNAseq_fastq2rawCounts_pdDataframe_geneWise(lib_list)
    d[lib_list].to_csv(raw_count_csv, sep=',')

    # write exp design file
    exp_file=open(exp_csv,'w')
    exp_file.write('lib_name,condition\n')
    for c in exp_design:
        for lib in exp_design[c]:
            exp_file.write('%s,%s\n' %(lib,c))
    exp_file.close()

    # write R script file
    r_file=open(r_fn,'w')
    r_file.write("library('DESeq2')\n")
    r_file.write("library(ggplot2)\n")
    r_file.write("setwd('%s')\n" % working_folder)
    r_file.write("countData <- read.csv('%s',header=TRUE,sep=',')\n" % raw_count_csv)
    r_file.write("metaData <- read.csv('%s',sep=',')\n" % exp_csv)
    r_file.write("dds<-DESeqDataSetFromMatrix(countData=countData, colData=metaData, design=~condition, tidy=TRUE)\n")
    r_file.write("dds <- DESeq(dds)\n")
    r_file.write("res <- results(dds)\n")
    r_file.write("summary(res)\n")
    r_file.write("plotMA(res, ylim=c(-2,2),cex=0.8)\n")
    r_file.write("write.csv(as.data.frame(res),file='%s')\n " % deseq2_result_csv)
    r_file.close()
    r_scriptrun='Rscript %s' % r_fn
    p=subprocess.run(r_scriptrun, shell=True, check=True,capture_output=True)
    
def MA_plot_RNAseq(x,y,HL_g):
# 2/18/2021 A simple MA plot using rpkm values for libs x and y
# it can also highligh a list genes with different colors
    d=RNAseq_fastq2rpkm_pdDataframe_geneWise([x,y])
    mean_exp=(d[x]+d[y])/2
    log2_ratio=np.log2(d[y]/d[x])
    plot_title="%s vs %s" % (x,y)
    fig, axs = plt.subplots()
    axs.scatter(mean_exp,log2_ratio,color='lightgray',marker='.')
    for g in HL_g:
        g_s=d['gene_name2']==g
        axs.scatter(mean_exp[g_s],log2_ratio[g_s],marker='o')
    axs.set(xlabel="mean_expression", ylabel=('log2ratio(%s/%s)' % (y,x)), xscale='log',title=plot_title)
    axs.legend(['']+HL_g,loc='upper left')
    plt.savefig(("MAplot_%s_vs_%s.png" % (y,x)),dpi=300)

def MA_plot_sRNAseq(x,y,HL_g):
# 3/5/2021 A simple MA plot using rpkm values for libs x and y
# it can also highligh a list genes with different colors
    d=sRNAseq_fasta2rpkm_pdDataframe_geneWise([x,y])
    mean_exp=(d[x]+d[y])/2
    log2_ratio=np.log2(d[y]/d[x])
    plot_title="%s vs %s" % (x,y)
    fig, axs = plt.subplots()
    axs.scatter(mean_exp,log2_ratio,color='lightgray',marker='.')
    for g in HL_g:
        g_s=d['gene_name2']==g
        axs.scatter(mean_exp[g_s],log2_ratio[g_s],marker='o')
    axs.set(xlabel="mean_expression", ylabel=('log2ratio(%s/%s)' % (y,x)), xscale='log',title=plot_title)
    axs.set_ylim([-6, 6])
    axs.legend(['']+HL_g,loc='upper right')
    plt.savefig(("MAplot_%s_vs_%s.png" % (y,x)),dpi=300)

def siRNA_track_plot_cer3(a,l): 
# 8/27/202 make a VSG plot for various Cer3 alleles, developed for siRNA suppression project
# a: Cer3 allele name, l: library name
# assumes 4nt barcodes
    fa_folder=dir_path(l+'_3linker_removed_collapsed.fa',fastq_folder)
    fa_file_name=fa_folder+l+'_3linker_removed_collapsed.fa'
    barcode_size=4
    
    # VSG options:
    ext_nt = 4 # number (nt) of extension when considering if the two tracks overlap, aka, minimal horizontal gap
    track_width = 2.0 # siRNA track width
    L_per_nt=1 # number of pixel per nt
    line_gap=4.0 #
    gDNA_track_width = 10.0
    track_width=10.0 # track width for sRNA

    bowtie_ref=bowtie_folder+'indexes/'+a
    bowtie_output=l+'_'+a+'_perfect_align.map'

    (sRNA_size_min, sRNA_size_max)=(20,24) # inclusive for both
    # cer3 alleles information, numbers are 1-based offset
    cer3_dict={'cer3_red20_25kb':{'gDNA_size': 25427, 'LTR5_start': 5864, 'LTR5_end': 6287, 'LTR3_start':14579, 'LTR3_end': 15002,
                                  'F58H7_5_start':13325, 'F58H7_5_end':14143, 'oma1_start': 7700, 'oma1_end': 8118,
                                  'MM': [{'pos':7709, 'WT':'T', 'mut':'C'}, {'pos':7739, 'WT':'C', 'mut':'G'},
                                         {'pos':7769, 'WT':'G', 'mut':'A'}, {'pos':7799, 'WT':'A', 'mut':'G'},
                                         {'pos':7829, 'WT':'T', 'mut':'C'}, {'pos':7859, 'WT':'G', 'mut':'A'},
                                         {'pos':7889, 'WT':'C', 'mut':'T'}, {'pos':7919, 'WT':'A', 'mut':'G'},
                                         {'pos':7952, 'WT':'A', 'mut':'C'}, {'pos':7979, 'WT':'G', 'mut':'A'},
                                         {'pos':8009, 'WT':'C', 'mut':'T'}, {'pos':8039, 'WT':'T', 'mut':'C'},
                                         {'pos':8069, 'WT':'T', 'mut':'C'}, {'pos':8099, 'WT':'T', 'mut':'C'}],
                                  'plot_flank': 800},
               'cer3_red40_25kb':{'gDNA_size': 25427, 'LTR5_start': 5864, 'LTR5_end': 6287, 'LTR3_start':14579, 'LTR3_end': 15002,
                                  'F58H7_5_start':13325, 'F58H7_5_end':14143, 'oma1_start': 7700, 'oma1_end': 8118,
                                  'MM': [{'pos':7719, 'WT':'A', 'mut':'G'}, {'pos':7749, 'WT':'A', 'mut':'G'},
                                     {'pos':7779, 'WT':'A', 'mut':'G'}, {'pos':7809, 'WT':'G', 'mut':'A'},
                                     {'pos':7839, 'WT':'C', 'mut':'T'}, {'pos':7869, 'WT':'T', 'mut':'G'},
                                     {'pos':7899, 'WT':'T', 'mut':'C'}, {'pos':7929, 'WT':'G', 'mut':'A'},
                                     {'pos':7959, 'WT':'C', 'mut':'T'}, {'pos':7989, 'WT':'A', 'mut':'G'},
                                     {'pos':8019, 'WT':'T', 'mut':'C'}, {'pos':8049, 'WT':'C', 'mut':'T'},
                                     {'pos':8079, 'WT':'G', 'mut':'C'}, {'pos':8109, 'WT':'A', 'mut':'G'}],
                                  'plot_flank': 800},
              'cer3_red46_25kb':{'gDNA_size': 25507, 'LTR5_start': 5864, 'LTR5_end': 6287, 'LTR3_start':14659, 'LTR3_end': 15082, # Cer3::Cer8 ASG
                                  'F58H7_5_start':13405, 'F58H7_5_end':14223, 'oma1_start': 7700, 'oma1_end': 8198,
                                  'MM': [{'pos':7720, 'WT':'C', 'mut':'G'}, {'pos':7740, 'WT':'A', 'mut':'T'},
                                         {'pos':7760, 'WT':'C', 'mut':'G'}, {'pos':7780, 'WT':'T', 'mut':'A'},
                                         {'pos':7800, 'WT':'T', 'mut':'A'}, {'pos':7820, 'WT':'A', 'mut':'T'},
                                         {'pos':7840, 'WT':'G', 'mut':'C'}, {'pos':7860, 'WT':'G', 'mut':'C'},
                                         {'pos':7880, 'WT':'T', 'mut':'A'}, {'pos':7900, 'WT':'C', 'mut':'G'},
                                         {'pos':7920, 'WT':'C', 'mut':'G'}, {'pos':7940, 'WT':'T', 'mut':'A'},
                                         {'pos':7960, 'WT':'C', 'mut':'G'}, {'pos':7980, 'WT':'A', 'mut':'T'},
                                         {'pos':8000, 'WT':'A', 'mut':'T'}, {'pos':8020, 'WT':'C', 'mut':'G'},
                                         {'pos':8040, 'WT':'A', 'mut':'T'}, {'pos':8060, 'WT':'C', 'mut':'G'},
                                         {'pos':8080, 'WT':'C', 'mut':'G'}, {'pos':8100, 'WT':'A', 'mut':'T'},
                                         {'pos':8120, 'WT':'G', 'mut':'C'}, {'pos':8140, 'WT':'A', 'mut':'T'},
                                         {'pos':8160, 'WT':'C', 'mut':'G'}],
                                  'plot_flank': 800},
               'cer3_red47_25kb':{'gDNA_size': 25473, 'LTR5_start': 5864, 'LTR5_end': 6287, 'LTR3_start':14625, 'LTR3_end': 15048,
                              'F58H7_5_start':13371, 'F58H7_5_end':14189, 'oma1_start': 7700, 'oma1_end': 8164,
                                  'MM': [{'pos':7729, 'WT':'A', 'mut':'G'}, {'pos':7759, 'WT':'T', 'mut':'C'},
                                         {'pos':7789, 'WT':'A', 'mut':'C'}, {'pos':7819, 'WT':'A', 'mut':'G'},
                                         {'pos':7849, 'WT':'T', 'mut':'C'}, {'pos':7879, 'WT':'A', 'mut':'G'},
                                         {'pos':7909, 'WT':'A', 'mut':'G'}, {'pos':7939, 'WT':'T', 'mut':'C'},
                                         {'pos':7969, 'WT':'T', 'mut':'C'}, {'pos':7999, 'WT':'A', 'mut':'G'},
                                         {'pos':8029, 'WT':'C', 'mut':'A'}, {'pos':8059, 'WT':'A', 'mut':'G'},
                                         {'pos':8089, 'WT':'A', 'mut':'G'}, {'pos':8119, 'WT':'C', 'mut':'T'},
                                         {'pos':8149, 'WT':'T', 'mut':'C'}],
                              'plot_flank': 800},
               'cer3_red49_25kb':{'gDNA_size': 25646, 'LTR5_start': 5864, 'LTR5_end': 6287, 'LTR3_start':14798, 'LTR3_end': 15221,
                              'F58H7_5_start':13544, 'F58H7_5_end':14362, 'oma1_start': 7700, 'oma1_end': 8337,
                                  'MM':[],
                              'plot_flank': 800},
    
               'cer3_red52_25kb':{ 'gDNA_size': 25477, 'LTR5_start': 5864, 'LTR5_end': 6287, 'LTR3_start':14629, 'LTR3_end': 15052, 
                                  'F58H7_5_start':13375, 'F58H7_5_end':14193,'oma1_start': 7700, 'oma1_end': 8168,
                                  'MM': [{'pos':7717, 'WT':'T', 'mut':'A'}, {'pos':7747, 'WT':'T', 'mut':'A'}, 
                                         {'pos':7777, 'WT':'C', 'mut':'G'}, {'pos':7807, 'WT':'G', 'mut':'C'},
                                         {'pos':7837, 'WT':'T', 'mut':'A'}, {'pos':7867, 'WT':'G', 'mut':'C'},
                                         {'pos':7897, 'WT':'G', 'mut':'C'}, {'pos':7927, 'WT':'A', 'mut':'T'},
                                         {'pos':7957, 'WT':'G', 'mut':'C'}, {'pos':7987, 'WT':'A', 'mut':'T'},
                                         {'pos':8017, 'WT':'G', 'mut':'C'}, {'pos':8047, 'WT':'C', 'mut':'G'},
                                         {'pos':8077, 'WT':'C', 'mut':'G'}, {'pos':8107, 'WT':'C', 'mut':'G'},
                                         {'pos':8137, 'WT':'A', 'mut':'T'}, {'pos':8167, 'WT':'G', 'mut':'C'}],
                                  'plot_flank': 800},
                'cer3_wt_25kb':{'gDNA_size': 25008, 'LTR5_start': 5864, 'LTR5_end': 6287, 'LTR3_start':14160, 'LTR3_end': 14583, 
                                'F58H7_5_start':12906, 'F58H7_5_end':13724, 'oma1_start': 7463, 'oma1_end': 7929, # oma1_start and oma1_end mark the positions of cer3 sequence inserted in oma-1::cer3(red57)
                                  'MM': [{'pos':7472, 'WT':'A', 'mut':'G'}, {'pos':7502, 'WT':'G', 'mut':'A'}, # WT: red57, mut: wt cer3
                                         {'pos':7532, 'WT':'T', 'mut':'C'}, {'pos':7562, 'WT':'A', 'mut':'G'},
                                         {'pos':7592, 'WT':'T', 'mut':'C'}, {'pos':7622, 'WT':'A', 'mut':'G'},
                                         {'pos':7652, 'WT':'T', 'mut':'C'}, {'pos':7682, 'WT':'A', 'mut':'G'},
                                         {'pos':7712, 'WT':'C', 'mut':'T'}, {'pos':7742, 'WT':'G', 'mut':'A'},
                                         {'pos':7772, 'WT':'A', 'mut':'G'}, {'pos':7802, 'WT':'G', 'mut':'A'},
                                         {'pos':7832, 'WT':'G', 'mut':'A'}, {'pos':7862, 'WT':'G', 'mut':'A'},
                                         {'pos':7892, 'WT':'A', 'mut':'G'}, {'pos':7922, 'WT':'A', 'mut':'G'}],
                                'plot_flank': 800},
               
                'cer8_red35_15kb':{'gDNA_size': 15000, 'LTR5_start': 1501, 'LTR5_end': 2074, 'LTR3_start':12886, 'LTR3_end': 13459, 
                                   'F58H7_5_start':13459, 'F58H7_5_end':13459,'oma1_start': 5507, 'oma1_end': 5928, # F58H7_5 numbers are place holders to avoid error messages
                                  'MM': [{'pos':5519, 'WT':'T', 'mut':'C'}, {'pos':5549, 'WT':'C', 'mut':'G'},
                                         {'pos':5579, 'WT':'G', 'mut':'A'}, {'pos':5609, 'WT':'A', 'mut':'G'},
                                         {'pos':5639, 'WT':'T', 'mut':'C'}, {'pos':5669, 'WT':'G', 'mut':'A'},
                                         {'pos':5699, 'WT':'C', 'mut':'T'}, {'pos':5729, 'WT':'A', 'mut':'G'},
                                         {'pos':5762, 'WT':'A', 'mut':'C'}, {'pos':5789, 'WT':'G', 'mut':'A'},
                                         {'pos':5819, 'WT':'C', 'mut':'T'}, {'pos':5849, 'WT':'T', 'mut':'C'},
                                         {'pos':5879, 'WT':'T', 'mut':'C'}, {'pos':5909, 'WT':'T', 'mut':'C'}],
                                  'plot_flank': 800},
               'cer3_red118_25kb':{'gDNA_size': 25391, 'LTR5_start': 5864, 'LTR5_end': 6287, 'LTR3_start':14543, 'LTR3_end': 14966,
                                  'F58H7_5_start':13285, 'F58H7_5_end':14107, 'oma1_start': 7700, 'oma1_end': 8082,
                                  'MM': [ {'pos':7723, 'WT':'C', 'mut':'G'},
                                         {'pos':7753, 'WT':'G', 'mut':'A'}, {'pos':7783, 'WT':'A', 'mut':'G'},
                                         {'pos':7813, 'WT':'T', 'mut':'C'}, {'pos':7843, 'WT':'G', 'mut':'A'},
                                         {'pos':7873, 'WT':'C', 'mut':'T'}, {'pos':7903, 'WT':'A', 'mut':'G'},
                                         {'pos':7916, 'WT':'T', 'mut':'C'}, {'pos':7943, 'WT':'G', 'mut':'A'},
                                         {'pos':7973, 'WT':'C', 'mut':'T'}, {'pos':8003, 'WT':'T', 'mut':'C'},
                                         {'pos':8033, 'WT':'T', 'mut':'C'}, {'pos':8063, 'WT':'T', 'mut':'C'}],
                                  'plot_flank': 800},
                'cer3_red122_25kb':{'gDNA_size': 25408, 'LTR5_start': 5864, 'LTR5_end': 6287, 'LTR3_start':14560, 'LTR3_end':14983 ,
                                  'F58H7_5_start': 13306, 'F58H7_5_end': 14124,  'oma1_start': 7700, 'oma1_end': 7700 + 399, 
                          'MM': [{'pos':7700+29, 'WT':'C', 'mut':'G'},{'pos':7700+59, 'WT':'T', 'mut':'A'},{'pos':7700+89, 'WT':'T', 'mut':'A'},{'pos':7700+119, 'WT':'A', 'mut':'T'},
                          {'pos':7700+149, 'WT':'A', 'mut':'T'},{'pos':7700+179, 'WT':'T', 'mut':'A'},{'pos':7700+209, 'WT':'A', 'mut':'T'},{'pos':7700+239, 'WT':'G', 'mut':'C'},
                          {'pos':7700+269, 'WT':'G', 'mut':'C'},{'pos':7700+299, 'WT':'A', 'mut':'T'},{'pos':7700+329, 'WT':'T', 'mut':'A'},{'pos':7700+359, 'WT':'T', 'mut':'A'},
                          {'pos':7700+389, 'WT':'T', 'mut':'A'}],'plot_flank': 800},
                'cer3_red123_25kb':{'gDNA_size': 25408, 'LTR5_start': 5864, 'LTR5_end': 6287, 'LTR3_start':14560, 'LTR3_end':14983 ,
                                  'F58H7_5_start': 13306, 'F58H7_5_end': 14124,  'oma1_start': 7700, 'oma1_end': 7700 + 399, 
                          'MM': [{'pos':7700+29, 'WT':'G', 'mut':'C'},{'pos':7700+59, 'WT':'A', 'mut':'T'},{'pos':7700+89, 'WT':'A', 'mut':'T'},{'pos':7700+119, 'WT':'A', 'mut':'T'},
                          {'pos':7700+149, 'WT':'C', 'mut':'G'},{'pos':7700+179, 'WT':'G', 'mut':'C'},{'pos':7700+209, 'WT':'G', 'mut':'C'},{'pos':7700+239, 'WT':'A', 'mut':'T'},
                          {'pos':7700+269, 'WT':'A', 'mut':'T'},{'pos':7700+299, 'WT':'C', 'mut':'G'},{'pos':7700+329, 'WT':'T', 'mut':'A'},{'pos':7700+359, 'WT':'G', 'mut':'C'},
                          {'pos':7700+389, 'WT':'G', 'mut':'C'}],'plot_flank': 800},
                'cer3_red124_25kb':{'gDNA_size': 25408, 'LTR5_start': 5864, 'LTR5_end': 6287, 'LTR3_start':14560, 'LTR3_end':14983 ,
                                  'F58H7_5_start': 13306, 'F58H7_5_end': 14124,  'oma1_start': 7700, 'oma1_end': 7700 + 399, 
                          'MM': [{'pos':7700+29, 'WT':'C', 'mut':'G'},{'pos':7700+59, 'WT':'C', 'mut':'G'},{'pos':7700+89, 'WT':'G', 'mut':'C'},{'pos':7700+119, 'WT':'T', 'mut':'A'},
                          {'pos':7700+149, 'WT':'T', 'mut':'A'},{'pos':7700+179, 'WT':'A', 'mut':'T'},{'pos':7700+209, 'WT':'C', 'mut':'G'},{'pos':7700+239, 'WT':'C', 'mut':'G'},
                          {'pos':7700+269, 'WT':'A', 'mut':'T'},{'pos':7700+299, 'WT':'C', 'mut':'G'},{'pos':7700+329, 'WT':'T', 'mut':'A'},{'pos':7700+359, 'WT':'T', 'mut':'A'},
                          {'pos':7700+389, 'WT':'T', 'mut':'A'}],'plot_flank': 800},
                'cer3_red125_25kb':{'gDNA_size': 25408, 'LTR5_start': 5864, 'LTR5_end': 6287, 'LTR3_start':14560, 'LTR3_end':14983 ,
                                  'F58H7_5_start': 13306, 'F58H7_5_end': 14124,  'oma1_start': 7700, 'oma1_end': 7700 + 399,
                          'MM': [{'pos':7700+29, 'WT':'T', 'mut':'A'},{'pos':7700+59, 'WT':'A', 'mut':'T'},{'pos':7700+89, 'WT':'G', 'mut':'C'},{'pos':7700+119, 'WT':'T', 'mut':'A'},
                          {'pos':7700+149, 'WT':'C', 'mut':'G'},{'pos':7700+179, 'WT':'G', 'mut':'C'},{'pos':7700+209, 'WT':'T', 'mut':'A'},{'pos':7700+239, 'WT':'T', 'mut':'A'},
                          {'pos':7700+269, 'WT':'A', 'mut':'T'},{'pos':7700+299, 'WT':'A', 'mut':'T'},{'pos':7700+329, 'WT':'A', 'mut':'T'},{'pos':7700+359, 'WT':'C', 'mut':'G'},
                          {'pos':7700+389, 'WT':'C', 'mut':'G'}],'plot_flank': 800},
                'cer3_red126_25kb':{'gDNA_size': 25407, 'LTR5_start': 5864, 'LTR5_end': 6287, 'LTR3_start':14559, 'LTR3_end':14982 ,
                                  'F58H7_5_start':13305, 'F58H7_5_end': 14123,  'oma1_start': 7700, 'oma1_end': 7700 + 398, 
                          'MM': [{'pos':7700+29, 'WT':'T', 'mut':'A'},{'pos':7700+59, 'WT':'A', 'mut':'T'},{'pos':7700+89, 'WT':'T', 'mut':'A'},{'pos':7700+119, 'WT':'T', 'mut':'A'},
                          {'pos':7700+149, 'WT':'T', 'mut':'A'},{'pos':7700+179, 'WT':'A', 'mut':'T'},{'pos':7700+209, 'WT':'G', 'mut':'C'},{'pos':7700+239, 'WT':'A', 'mut':'T'},
                          {'pos':7700+269, 'WT':'A', 'mut':'T'},{'pos':7700+299, 'WT':'C', 'mut':'G'},{'pos':7700+329, 'WT':'T', 'mut':'A'},{'pos':7700+359, 'WT':'T', 'mut':'A'},
                          {'pos':7700+389, 'WT':'C', 'mut':'G'}],'plot_flank': 800},
                'cer3_red127_25kb':{'gDNA_size': 25408, 'LTR5_start': 5864, 'LTR5_end': 6287, 'LTR3_start':14560, 'LTR3_end':14983 ,
                                  'F58H7_5_start': 13306, 'F58H7_5_end': 14124,  'oma1_start': 7700, 'oma1_end': 7700 + 399,
                          'MM': [{'pos':7700+29, 'WT':'T', 'mut':'A'},{'pos':7700+59, 'WT':'T', 'mut':'A'},{'pos':7700+89, 'WT':'A', 'mut':'T'},{'pos':7700+119, 'WT':'C', 'mut':'G'},
                          {'pos':7700+149, 'WT':'C', 'mut':'G'},{'pos':7700+179, 'WT':'A', 'mut':'T'},{'pos':7700+209, 'WT':'A', 'mut':'T'},{'pos':7700+239, 'WT':'C', 'mut':'G'},
                          {'pos':7700+269, 'WT':'C', 'mut':'G'},{'pos':7700+299, 'WT':'T', 'mut':'A'},{'pos':7700+329, 'WT':'T', 'mut':'A'},{'pos':7700+359, 'WT':'G', 'mut':'C'},
                          {'pos':7700+388, 'WT':'A', 'mut':'T'}],'plot_flank': 800},
                'cer3_red129_25kb':{'gDNA_size': 25407, 'LTR5_start': 5864, 'LTR5_end': 6287, 'LTR3_start':14559, 'LTR3_end':14982 ,
                                  'F58H7_5_start':13305, 'F58H7_5_end': 14123,  'oma1_start': 7700, 'oma1_end': 7700 + 398, 
                          'MM': [{'pos':7700+29, 'WT':'A', 'mut':'T'},{'pos':7700+59, 'WT':'T', 'mut':'A'},{'pos':7700+89, 'WT':'G', 'mut':'C'},{'pos':7700+119, 'WT':'T', 'mut':'A'},
                          {'pos':7700+149, 'WT':'A', 'mut':'T'},{'pos':7700+179, 'WT':'G', 'mut':'C'},{'pos':7700+209, 'WT':'T', 'mut':'A'},{'pos':7700+239, 'WT':'G', 'mut':'C'},
                          {'pos':7700+269, 'WT':'A', 'mut':'T'},{'pos':7700+299, 'WT':'G', 'mut':'C'},{'pos':7700+329, 'WT':'T', 'mut':'A'},{'pos':7700+359, 'WT':'A', 'mut':'T'},
                          {'pos':7700+389, 'WT':'C', 'mut':'G'}],'plot_flank': 800},
                'cer3_red130_25kb':{'gDNA_size': 25408, 'LTR5_start': 5864, 'LTR5_end': 6287, 'LTR3_start':14560, 'LTR3_end':14983 ,
                                  'F58H7_5_start': 13306, 'F58H7_5_end': 14124,  'oma1_start': 7700, 'oma1_end': 7700 + 399,
                          'MM': [{'pos':7700+29, 'WT':'G', 'mut':'C'},{'pos':7700+59, 'WT':'G', 'mut':'C'},{'pos':7700+89, 'WT':'C', 'mut':'G'},{'pos':7700+119, 'WT':'A', 'mut':'T'},
                          {'pos':7700+149, 'WT':'C', 'mut':'G'},{'pos':7700+179, 'WT':'C', 'mut':'G'},{'pos':7700+209, 'WT':'T', 'mut':'A'},{'pos':7700+239, 'WT':'T', 'mut':'A'},
                          {'pos':7700+269, 'WT':'G', 'mut':'C'},{'pos':7700+299, 'WT':'T', 'mut':'A'},{'pos':7700+329, 'WT':'A', 'mut':'T'},{'pos':7700+359, 'WT':'T', 'mut':'A'},
                          {'pos':7700+389, 'WT':'T', 'mut':'A'}],'plot_flank': 800},
                'cer3_red131_25kb':{'gDNA_size': 25402, 'LTR5_start': 5864, 'LTR5_end': 6287, 'LTR3_start':14554, 'LTR3_end': 14977,
                                  'F58H7_5_start':13300, 'F58H7_5_end': 14118,  'oma1_start': 7700, 'oma1_end': 7700 + 393, 
                          'MM': [{'pos':7700+29, 'WT':'T', 'mut':'A'},{'pos':7700+59, 'WT':'T', 'mut':'A'},{'pos':7700+89, 'WT':'A', 'mut':'T'},{'pos':7700+119, 'WT':'A', 'mut':'T'},
                          {'pos':7700+149, 'WT':'A', 'mut':'T'},{'pos':7700+179, 'WT':'C', 'mut':'G'},{'pos':7700+209, 'WT':'A', 'mut':'T'},{'pos':7700+239, 'WT':'T', 'mut':'A'},
                          {'pos':7700+269, 'WT':'T', 'mut':'A'},{'pos':7700+299, 'WT':'C', 'mut':'G'},{'pos':7700+329, 'WT':'A', 'mut':'T'},{'pos':7700+359, 'WT':'G', 'mut':'C'},
                          {'pos':7700+389, 'WT':'C', 'mut':'G'}],'plot_flank': 800},
                'cer3_red132a_25kb':{'gDNA_size': 25408, 'LTR5_start': 5864, 'LTR5_end': 6287, 'LTR3_start':14560, 'LTR3_end':14983 ,
                                  'F58H7_5_start': 13306, 'F58H7_5_end': 14124,  'oma1_start': 7700, 'oma1_end': 7700 + 399,
                          'MM': [{'pos':7700+29, 'WT':'C', 'mut':'G'},{'pos':7700+59, 'WT':'C', 'mut':'G'},{'pos':7700+89, 'WT':'C', 'mut':'G'},{'pos':7700+119, 'WT':'T', 'mut':'A'},
                          {'pos':7700+149, 'WT':'T', 'mut':'A'},{'pos':7700+179, 'WT':'T', 'mut':'A'},{'pos':7700+209, 'WT':'A', 'mut':'T'},{'pos':7924, 'WT':'C', 'mut':'T'},{'pos':7700+239, 'WT':'A', 'mut':'T'},
                          {'pos':7700+269, 'WT':'A', 'mut':'T'},{'pos':7700+299, 'WT':'A', 'mut':'T'},{'pos':7700+329, 'WT':'C', 'mut':'G'},{'pos':7700+359, 'WT':'C', 'mut':'G'},
                          {'pos':7700+389, 'WT':'A', 'mut':'T'}],'plot_flank': 800},
               
               'cer3_red147a_25kb':{'gDNA_size': 25400, 'LTR5_start': 5864, 'LTR5_end': 6287, 'LTR3_start':8091+6461, 'LTR3_end': 8091+6884,
                              'F58H7_5_start':8091+5207, 'F58H7_5_end':8091+6025, 'oma1_start': 7700, 'oma1_end': 8091,
                                  'MM':[],
                              'plot_flank': 800},
               'cer3_red146_25kb':{'gDNA_size': 25400, 'LTR5_start': 5864, 'LTR5_end': 6287, 'LTR3_start':8101+6461, 'LTR3_end': 8101+6884,
                              'F58H7_5_start':8101+5207, 'F58H7_5_end':8101+6025, 'oma1_start': 7700, 'oma1_end': 8101,
                                  'MM':[],
                              'plot_flank': 800},
               'cer3_red145_25kb':{'gDNA_size': 25400, 'LTR5_start': 5864, 'LTR5_end': 6287, 'LTR3_start':8103+6461, 'LTR3_end': 8103+6884,
                              'F58H7_5_start':8103+5207, 'F58H7_5_end':8103+6025, 'oma1_start': 7700, 'oma1_end': 8103,
                                  'MM':[],
                              'plot_flank': 800},
               'cer3_red153_25kb':{'gDNA_size': 25426, 'LTR5_start': 5864, 'LTR5_end': 6287, 'LTR3_start':14578, 'LTR3_end': 15001,
                              'F58H7_5_start':13324, 'F58H7_5_end':14142, 'oma1_start': 7700, 'oma1_end': 8117,
                                  'MM':[],
                              'plot_flank': 800},
               'cer3_red154_25kb':{'gDNA_size': 25426, 'LTR5_start': 5864, 'LTR5_end': 6287, 'LTR3_start':14551, 'LTR3_end': 14974,
                              'F58H7_5_start':13297, 'F58H7_5_end':14115, 'oma1_start': 7700, 'oma1_end': 8090,
                                  'MM':[],
                              'plot_flank': 800},
               'cer3_red152_25kb':{'gDNA_size': 25426, 'LTR5_start': 5864, 'LTR5_end': 6287, 'LTR3_start':14580, 'LTR3_end': 15003,
                              'F58H7_5_start':13326, 'F58H7_5_end':14144, 'oma1_start': 7700, 'oma1_end': 8119,
                                  'MM':[],
                              'plot_flank': 800},
               'cer3_red155_25kb':{'gDNA_size': 25426, 'LTR5_start': 5864, 'LTR5_end': 6287, 'LTR3_start':14551, 'LTR3_end': 14974,
                              'F58H7_5_start':13297, 'F58H7_5_end':14115, 'oma1_start': 7700, 'oma1_end': 8090,
                                  'MM':[],
                              'plot_flank': 800},
               'cer3_red166b_25kb':{'gDNA_size': 25426, 'LTR5_start': 5864, 'LTR5_end': 6287, 'LTR3_start':14550, 'LTR3_end': 14973,
                              'F58H7_5_start':13296, 'F58H7_5_end':14114, 'oma1_start': 7700, 'oma1_end': 8089,
                                  'MM':[],
                              'plot_flank': 800},
               'cer3_red167b_25kb':{'gDNA_size': 25426, 'LTR5_start': 5864, 'LTR5_end': 6287, 'LTR3_start':14337, 'LTR3_end': 14760,
                              'F58H7_5_start':13083, 'F58H7_5_end':13901, 'oma1_start': 7700, 'oma1_end': 7882,
                                  'MM':[],
                              'plot_flank': 800},
               
              }
    
    bowtie_prog=bowtie_folder+'bowtie ' + bowtie_ref + ' -5 '+str(barcode_size)+' -f -t -v 0 -a -p '+str(core_num)+' '+fa_file_name+' '+bowtie_output
    p=subprocess.run(bowtie_prog, shell=True, check=True,capture_output=True)
    #print(p.stderr)
    
    # a dictionary for cer3 reads, distinguishing oma-1 reads (definitively generator-derived), oma-1 reads (ambiguous), cer3-only reads, flanking reads
    cer3_reads_dict={'oma1_gen':[],'oma1_ambi':[],'cer3only':[],'flanking':[],'LTR':[]}
    infile_bowtie_map=open(bowtie_output, 'r')
    c=0
    for L1 in infile_bowtie_map:
        c+=1
        A1=L1.strip().split('\t')
        (t_name, t_start, t_size, t_end, t_strand, t_else) =(A1[0], int(A1[3])+1, len(A1[4]), int(A1[3])+len(A1[4]), A1[1],int(A1[6]))
        (t_oma1_gen, t_oma1_ambi,t_cer3only,t_flanking) = (False,False,False,False)
        # check if the read is in the flanking region:
        if t_end<cer3_dict[a]['LTR5_start'] or t_start>cer3_dict[a]['LTR3_end'] :
            cer3_reads_dict['flanking'].append({'start': t_start,'end': t_end,'strand':t_strand,'name': t_name,'color':'gray'})
        # check if the read is cer3 only:
        elif (t_end>=cer3_dict[a]['LTR5_start'] and t_end<cer3_dict[a]['oma1_start']) or (t_start>cer3_dict[a]['oma1_end'] and t_start<=cer3_dict[a]['LTR3_end']):
            cer3_reads_dict['cer3only'].append({'start': t_start,'end': t_end,'strand':t_strand,'name': t_name,'color':'black'})
            #check if the read is in LTR
            if (t_start>=cer3_dict[a]['LTR5_start'] and t_end<=cer3_dict[a]['LTR5_end']) or (t_start>=cer3_dict[a]['LTR3_start'] and t_end<=cer3_dict[a]['LTR3_end']):
                cer3_reads_dict['LTR'].append({'start': t_start,'end': t_end,'strand':t_strand,'name': t_name,'color':'black'})

        # check if the read is definitely derivded from oma1 sequence in the cer3
        # cover a MM position or cer3:oma-1/oma-1:cer3 junction


        else:
            for M in cer3_dict[a]['MM']:
                if M['pos'] in range(t_start,t_end+1):
                    if t_start>=cer3_dict[a]['oma1_start'] and t_end<=cer3_dict[a]['oma1_end']:
                        t_oma1_gen=True
                        cer3_reads_dict['oma1_gen'].append({'start': t_start,'end': t_end,'strand':t_strand,'name': t_name,'color':'red'})
                    break
            if t_oma1_gen==False:
                if ((cer3_dict[a]['oma1_start']-1) in range(t_start,t_end)) or ((cer3_dict[a]['oma1_end']+1) in range(t_start+1,t_end+1)): # label the junction reads
                    t_oma1_gen=True
                    cer3_reads_dict['oma1_gen'].append({'start': t_start,'end': t_end,'strand':t_strand,'name': t_name,'color':'grey'})
            if t_oma1_gen==False:
                cer3_reads_dict['oma1_ambi'].append({'start': t_start,'end': t_end,'strand':t_strand,'name': t_name,'color':'gray'})            
    infile_bowtie_map.close()
    rm_cmd='rm '+bowtie_output
    subprocess.run(rm_cmd, shell=True)

    # Calculate the number of siRNA in different categories.
    report_name='report_'+l+'_'+a+'.txt'
    report_out=open(report_name,'w')
    N_total_col=0 # total number of collapsed reads
    N_total_unCol=0 # total number of un-collapsed reads

    infile=open(fa_file_name,'r')
    for L1 in infile:
        if L1[0]=='>':
            N_total_col+=1
            N_total_unCol+=int(L1.strip().split('_')[-1])
    infile.close()
    report_out.write(l+'\n')
    report_out.write("Reference sequence for the alignment: %s" % a)
    report_out.write("\nTotal number of collapsed reads: %i" % N_total_col)
    report_out.write("\nTotal number of Uncollapsed reads: %i" % N_total_unCol)

    # align to whole ce6
    bowtie_ref=bowtie_folder+'indexes/ce6'
    N_ce6_col=0
    N_ce6_unCol=0
    bowtie_output=a+'_ce6_align.map'
    bowtie_prog=bowtie_folder+'bowtie ' + bowtie_ref + ' -5 '+str(barcode_size)+' -f -t -v 0 -a -p '+ str(core_num)+' '+fa_file_name+' '+bowtie_output
    p=subprocess.run(bowtie_prog, shell=True, check=True,capture_output=True)
    #print(p.stderr)
    
    infile_map=open(bowtie_output,'r')
    for L1 in infile_map:
        A1=L1.strip().split('\t')
        N_ce6_col+=1.0/(1+int(A1[6]))
        N_ce6_unCol+=int((A1[0].split('_'))[-1])/(1+int(A1[6]))
    report_out.write("\nTotal number of reads that perfectly match to ce6:")
    report_out.write("\nCollapsed reads: %i" % N_ce6_col)
    report_out.write("\nUncollapsed reads: %i" % N_ce6_unCol)
    infile_map.close()
    rm_cmd='rm '+bowtie_output
    subprocess.run(rm_cmd, shell=True)
    
    # Total number of cer3 siRNA,  no oma-1 sequence
    (N_antisense, N_sense)=(0,0)
    for T in cer3_reads_dict['cer3only']:
        if ((T['end']-T['start']+1) in range(sRNA_size_min, sRNA_size_max+1)) :
            if T['strand']=='-': N_antisense+=1
            if T['strand']=='+': N_sense+=1
    report_out.write("\n Total number of cer3 siRNA, antisense, including oma-1 junction: %i" % N_antisense)
    report_out.write("\n Total number of cer3 siRNA, sense, including oma-1 junction: %i" % N_sense)

    # Total number of oma-1 siRNAs with MM or cer3 junction
    (N_antisense, N_sense)=(0,0)
    for T in cer3_reads_dict['oma1_gen']:
        if ((T['end']-T['start']+1) in range(sRNA_size_min, sRNA_size_max+1)) :
            if T['strand']=='-': N_antisense+=1
            if T['strand']=='+': N_sense+=1

    report_out.write("\n Total number of generator siRNAs with MM or cer3 junction, antisense: %i" % N_antisense)
    report_out.write("\n Total number of generator siRNAs with MM or cer3 junction, sense: %i" % N_sense)

    # Total number of red20-matching oma-1 siRNAs without MM or cer3 junction
    (N_antisense, N_sense)=(0,0)
    for T in cer3_reads_dict['oma1_ambi']:
        if  ((T['end']-T['start']+1) in range(sRNA_size_min, sRNA_size_max+1)) :
            if T['strand']=='-': N_antisense+=1
            if T['strand']=='+': N_sense+=1
    report_out.write("\nTotal number of generator siRNAs without MM or cer3 junction, antisense: %i" % N_antisense)
    report_out.write("\nTotal number of generator siRNAs without MM or cer3 junction, sense: %i" % N_sense)
    
    
    # track dictionary: key: line (1, 2, etc), value: a list of tracks
    track_dict_sense = {1:[]}
    track_dict_antisense = {1:[]}
    for T_type in ['oma1_gen','oma1_ambi','cer3only','flanking']:
        for T1 in cer3_reads_dict[T_type]:
            if T1['strand'] == '+' and ((T1['end']-T1['start']+1) in range(sRNA_size_min, sRNA_size_max+1)) and (T1['start']>=cer3_dict[a]['LTR5_start']-cer3_dict[a]['plot_flank'] and T1['end']<=cer3_dict[a]['LTR3_end']+cer3_dict[a]['plot_flank']):
                overlap_flag='no'
                for track_line in range(1,len(track_dict_sense)+1): # check the latest track line
                    overlap_flag='no'
                    for T2 in track_dict_sense[track_line]:
                        # check if T1 and T2 overlap (extend start and end by ext_nt nt)
                        if T1['start'] in range(T2['start']-ext_nt, T2['end']+ext_nt) or T1['end'] in range(T2['start']-ext_nt, T2['end']+ext_nt):
                            overlap_flag='yes'
                            break
                        if T2['start'] in range(T1['start']-ext_nt, T1['end']+ext_nt) or T2['end'] in range(T1['start']-ext_nt, T1['end']+ext_nt) :
                            overlap_flag='yes'
                            break
                    if overlap_flag == 'no': # if no overlap with any of the existing tracks in the current line
                        track_dict_sense[track_line].append({'start': T1['start'],
                                                             'end': T1['end'],
                                                             'type':T_type,
                                                             'color': T1['color']})
                        break
                if overlap_flag == 'yes': # if it overlaps with a track in all existing line, add a new line, add the track to the new line
                    track_dict_sense[len(track_dict_sense)+1]=[{'start': T1['start'],
                                                             'end': T1['end'],
                                                             'type':T_type,
                                                             'color': T1['color']}]

            if T1['strand'] == '-' and ((T1['end']-T1['start']+1) in range(sRNA_size_min, sRNA_size_max+1)) and (T1['start']>=cer3_dict[a]['LTR5_start']-cer3_dict[a]['plot_flank'] and T1['end']<=cer3_dict[a]['LTR3_end']+cer3_dict[a]['plot_flank']):
                for track_line in range(1,len(track_dict_antisense)+1): # check the latest track line
                    overlap_flag='no'
                    for T2 in track_dict_antisense[track_line]:
                        overlap_flag='no'
                        # check if T1 and T2 overlap (extend start and end by ext_nt nt)
                        if T1['start'] in range(T2['start']-ext_nt, T2['end']+ext_nt) or T1['end'] in range(T2['start']-ext_nt, T2['end']+ext_nt):
                            overlap_flag='yes'
                            break
                        if T2['start'] in range(T1['start']-ext_nt, T1['end']+ext_nt) or T2['end'] in range(T1['start']-ext_nt, T1['end']+ext_nt) :
                            overlap_flag='yes'
                            break
                    if overlap_flag == 'no': # if no overlap with any of the existing tracks in the current line
                        track_dict_antisense[track_line].append({'start': T1['start'],
                                                             'end': T1['end'],
                                                             'type':T_type,
                                                             'color': T1['color']})
                        break
                if overlap_flag == 'yes': # if it overlaps with a track in all existing line, add a new line, add the track to the new line
                    track_dict_antisense[len(track_dict_antisense)+1]=[{'start': T1['start'],
                                                             'end': T1['end'],
                                                             'type':T_type,
                                                             'color': T1['color']}]


    # VSG plot
    VSG_py_tracks_name=l+'_'+a+'.py'
    VSG_plot_tracks_name=l+'_'+a+'.svg'
    outfile_VSG_py=open(VSG_py_tracks_name,'w')
    outfile_VSG_py.write('from VSG_Module import * \n')

    (LTR5_start,LTR5_end)=(cer3_dict[a]['LTR5_start'],cer3_dict[a]['LTR5_end'])  #<<<
    (LTR3_start,LTR3_end)=(cer3_dict[a]['LTR3_start'],cer3_dict[a]['LTR3_end']) # <<<
    (F58H7_5_start,F58H7_5_end)=(cer3_dict[a]['F58H7_5_start'],cer3_dict[a]['F58H7_5_end'])  # <<<
    (plot_start,plot_end)= (LTR5_start-cer3_dict[a]['plot_flank'], LTR3_end+cer3_dict[a]['plot_flank']) # start and end of the reference DNA used for plotting
    (oma1_start,oma1_end)=(cer3_dict[a]['oma1_start'],cer3_dict[a]['oma1_end'])

    #write experiment info on the plot: #skip this part for now because linux doesn't deal with vtext 
    #VSG_x=6000*L_per_nt
    #VSG_y=-1.0*(len(track_dict_sense)-10)*(track_width+2*line_gap)
    #outfile_VSG_py.write("vtext(text='"+exp+sample_info+"', yc="+str(VSG_y)+", x1="+str(VSG_x)+", font='courier 50 Bold', stroke=black)\n")
    #VSG_y=-1.0*(len(track_dict_sense)-8)*(track_width+2*line_gap)
    #outfile_VSG_py.write("vtext(text='"+a+"', yc="+str(VSG_y)+", x1="+str(VSG_x)+", font='courier 50 Bold', stroke=black)\n")
    #VSG_y=-1.0*(len(track_dict_sense)-6)*(track_width+2*line_gap)
    #outfile_VSG_py.write("vtext(text='total collapsed reads: "+str(N_total_col)+"', yc="+str(VSG_y)+", x1="+str(VSG_x)+", font='courier 50 Bold', stroke=black)\n")
    #VSG_y=-1.0*(len(track_dict_sense)-4)*(track_width+2*line_gap)
    #outfile_VSG_py.write("vtext(text='total collapsed reads that match to C. elegans genome: "+str(N_ce6_col)+"', yc="+str(VSG_y)+", x1="+str(VSG_x)+", font='courier 50 Bold',  stroke=black)\n")

    for track_line in range(len(track_dict_sense),0,-1):
        #print track_line
        VSG_y=-1.0*(len(track_dict_sense)-track_line)*(track_width+line_gap)
        for T in track_dict_sense[track_line]:
            (T_start, T_end)=(T['start'],T['end'])
            if T_start in range(plot_start,plot_end+1) and T_end in range(plot_start,plot_end+1) :
                outfile_VSG_py.write('vline(x1='+str(T_start*L_per_nt)+', x2='+str(T_end*L_per_nt)+', y1='+str(VSG_y)+', y2='+str(VSG_y)+', strokewidth='+str(track_width)+' ,stroke='+T['color']+')\n')

    # plot the gDNA b/t plot_start and plot_end
    gDNA_xc=(plot_end+plot_start)/2*L_per_nt
    gDNA_yc=-1.0*(len(track_dict_sense))*(track_width+line_gap) - 0.5*gDNA_track_width
    gDNA_width=(plot_end-plot_start+1)*L_per_nt
    outfile_VSG_py.write("vrect(xc="+str(gDNA_xc)+", yc="+str(gDNA_yc)+", width="+str(gDNA_width)+", height="+ str(gDNA_track_width) +", fill=grey,stroke=grey )\n")

    # plot the 5LTR_start to 3LTR_end 
    gDNA_xc=(LTR5_start+LTR3_end)/2*L_per_nt
    gDNA_yc=-1.0*(len(track_dict_sense))*(track_width+line_gap) - 0.5*gDNA_track_width
    gDNA_width=(LTR3_end-LTR5_start+1)*L_per_nt
    outfile_VSG_py.write("vrect(xc="+str(gDNA_xc)+", yc="+str(gDNA_yc)+", width="+str(gDNA_width)+", height="+ str(gDNA_track_width) +", fill=orange,stroke=orange )\n")

    # plot the 5 LTR
    gDNA_xc=(LTR5_start+LTR5_end)/2*L_per_nt
    gDNA_yc=-1.0*(len(track_dict_sense))*(track_width+line_gap) - 0.5*gDNA_track_width
    gDNA_width=(LTR5_end-LTR5_start+1)*L_per_nt
    outfile_VSG_py.write("vrect(xc="+str(gDNA_xc)+", yc="+str(gDNA_yc)+", width="+str(gDNA_width)+", height="+ str(gDNA_track_width) +", fill=blue,stroke=blue )\n")

    # plot oma1 insert
    gDNA_xc=(oma1_start+oma1_end)/2*L_per_nt
    gDNA_yc=-1.0*(len(track_dict_sense))*(track_width+line_gap) - 0.5*gDNA_track_width
    gDNA_width=(oma1_end-oma1_start+1)*L_per_nt
    outfile_VSG_py.write("vrect(xc="+str(gDNA_xc)+", yc="+str(gDNA_yc)+", width="+str(gDNA_width)+", height="+ str(gDNA_track_width) +", fill=purple,stroke=purple )\n")

    # plot the mismatches:
    for M in cer3_dict[a]['MM']:
        (L_x, L_y1, L_y2)=(M['pos']*L_per_nt,-1.0*(len(track_dict_sense))*(track_width+line_gap) - 1.0*gDNA_track_width,-1.0*(len(track_dict_sense))*(track_width+line_gap))
        outfile_VSG_py.write("vline(x1="+str(L_x)+", y1="+str(L_y1)+", x2="+str(L_x)+", y2="+ str(L_y2) +", strokewidth=2, stroke=white )\n")

    # plot the 3 LTR
    gDNA_xc=(LTR3_start+LTR3_end)/2*L_per_nt
    gDNA_yc=-1.0*(len(track_dict_sense))*(track_width+line_gap) - 0.5*gDNA_track_width
    gDNA_width=(LTR3_end-LTR3_start+1)*L_per_nt
    outfile_VSG_py.write("vrect(xc="+str(gDNA_xc)+", yc="+str(gDNA_yc)+", width="+str(gDNA_width)+", height="+ str(gDNA_track_width) +", fill=blue,stroke=blue )\n")

    # plot the F58H7.5
    gDNA_xc=(F58H7_5_start+F58H7_5_end)/2*L_per_nt
    gDNA_yc=-1.0*(len(track_dict_sense))*(track_width+line_gap) - 0.5*gDNA_track_width
    gDNA_width=(F58H7_5_end-F58H7_5_start+1)*L_per_nt
    outfile_VSG_py.write("vrect(xc="+str(gDNA_xc)+", yc="+str(gDNA_yc)+", width="+str(gDNA_width)+", height="+ str(gDNA_track_width) +", fill=red,stroke=red )\n")

    for track_line in range(1,len(track_dict_antisense)+1):
        #print track_line
        VSG_y=-1.0*(len(track_dict_sense)+track_line)*(track_width+line_gap) - 1.5*gDNA_track_width
        for T in track_dict_antisense[track_line]:
            (T_start, T_end)=(T['start'],T['end'])
            if T_start in range(plot_start,plot_end+1) and T_end in range(plot_start,plot_end+1) :
                outfile_VSG_py.write('vline(x1='+str(T_start*L_per_nt)+', x2='+str(T_end*L_per_nt)+', y1='+str(VSG_y)+', y2='+str(VSG_y)+', strokewidth='+str(track_width)+' ,stroke='+T['color']+')\n')


    outfile_VSG_py.write("vset(bg='white')\n")
    outfile_VSG_py.write("vdisplay('"+VSG_plot_tracks_name+"')")
    outfile_VSG_py.close()

    vsg_cmd='python '+VSG_py_tracks_name
    subprocess.run(vsg_cmd,shell=True)
    
    # calculate the pval using Wilcoxon-Mann-Whitney test, 11/30/2020
    # create bins for oma-1 insert
    cer3_allele=a
    bin_size=50
    bin_insert={}
    insert_counts=0 # number of generator 22Gs
    for k in range(cer3_dict[cer3_allele]['oma1_start'],cer3_dict[cer3_allele]['oma1_end'],bin_size):
        if k+bin_size<cer3_dict[cer3_allele]['oma1_end']:
            bin_insert[k]=0

    #create bins for flanking regions:
    bin_flank={}
    flank_size=400
    flank_left_counts=0 # number of antisense 22Gs in the left flank regions
    flank_right_counts=0 # number of antisense 22Gs in the right flank regions

    for k in range(cer3_dict[cer3_allele]['oma1_start']-flank_size,cer3_dict[cer3_allele]['oma1_start'],bin_size):
        if k+bin_size<=cer3_dict[cer3_allele]['oma1_start']:
            bin_flank[k]=0

    for k in range(cer3_dict[cer3_allele]['oma1_end']+1,cer3_dict[cer3_allele]['oma1_end']+1+flank_size,bin_size):
        if k+bin_size<=cer3_dict[cer3_allele]['oma1_end']+1+flank_size:
            bin_flank[k]=0

    # assign alignments to bins
    for T in cer3_reads_dict['oma1_gen']:
        if T['strand']=='-' and (T['end']-T['start']+1) in range(sRNA_size_min, sRNA_size_max+1) and (T['color']=='red'): # note that junction reads are excluded in this analysis
            insert_counts+=1
            T_bin=int((T['end']-cer3_dict[cer3_allele]['oma1_start'])/bin_size)*bin_size+cer3_dict[cer3_allele]['oma1_start']
            if T_bin in bin_insert:
                bin_insert[T_bin]+=1

    for T in cer3_reads_dict['oma1_ambi']:
        if T['strand']=='-' and (T['end']-T['start']+1) in range(sRNA_size_min, sRNA_size_max+1):
            insert_counts+=1
            T_bin=int((T['end']-cer3_dict[cer3_allele]['oma1_start'])/bin_size)*bin_size+cer3_dict[cer3_allele]['oma1_start']
            if T_bin in bin_insert:
                bin_insert[T_bin]+=1

    for T in cer3_reads_dict['cer3only']:
        if T['strand']=='-' and (T['end']-T['start']+1) in range(sRNA_size_min, sRNA_size_max+1) and T['start'] in range(cer3_dict[cer3_allele]['oma1_end']+1,(cer3_dict[cer3_allele]['oma1_end']+flank_size)):
            flank_right_counts+=1
            T_bin=int((T['end']-cer3_dict[cer3_allele]['oma1_end'])/bin_size)*bin_size+cer3_dict[cer3_allele]['oma1_end']+1
            if T_bin in bin_flank:
                bin_flank[T_bin]+=1

    for T in cer3_reads_dict['cer3only']:
        if T['strand']=='-' and (T['end']-T['start']+1) in range(sRNA_size_min, sRNA_size_max+1) and T['end'] in range(cer3_dict[cer3_allele]['oma1_start']-flank_size,cer3_dict[cer3_allele]['oma1_start']):
            flank_left_counts+=1
            T_bin=int((T['end']-(cer3_dict[cer3_allele]['oma1_start']-flank_size))/bin_size)*bin_size+(cer3_dict[cer3_allele]['oma1_start']-flank_size)
            if T_bin in bin_flank:
                bin_flank[T_bin]+=1
    (generator_counts,flank_counts)=([],[])
    for k in bin_insert:
        generator_counts.append(bin_insert[k])
    for k in bin_flank:
        flank_counts.append(bin_flank[k])
    
    report_out.write("\n\nWilcoxon-Mann-Whitney test. H0: generator 22G level = flanking 22G level. H1: generator 22G level < flanking 22G level")
    report_out.write("\nbin size: %i NT; generator bin number: %i; flanking bin number: %i" % (bin_size, len(bin_insert), len(bin_flank)) )
    report_out.write("n\counts in the generator bin: ")
    for w in generator_counts:
        report_out.write('%i,' % w)
    report_out.write("n\counts in the flank bin: ")
    for w in flank_counts:
        report_out.write('%i,' % w)
    report_out.write("\np-value: %f" % mannwhitneyu(generator_counts,flank_counts,alternative='less').pvalue)
    
    # report the degree of 22G suppression:
    report_out.write("\n\n22G-RNA counts in the left flanking region: %i (%.2f RPKM)" % (flank_left_counts, flank_left_counts*1000000/N_ce6_col/(flank_size/1000)))
    report_out.write("\n\n22G-RNA counts in the right flanking region: %i (%.2f RPKM)" % (flank_right_counts, flank_right_counts*1000000/N_ce6_col/(flank_size/1000)))
    report_out.write("\n\n22G-RNA counts in the generator region: %i (%.2f RPKM)" % (insert_counts, insert_counts*1000000/N_ce6_col/((cer3_dict[cer3_allele]['oma1_end']-cer3_dict[cer3_allele]['oma1_start'])/1000)))
    ss=repr('The degree of 22G suppression (generator 22G densitiy / flanking region 22G desity)')
    report_out.write("\n\n%s: %f" % (ss,(insert_counts/(cer3_dict[cer3_allele]['oma1_end']-cer3_dict[cer3_allele]['oma1_start']))/((flank_left_counts+flank_right_counts)/(2*flank_size))))
    report_out.close()
    
    #print(generator_counts,flank_counts)
    
def siRNA_track_plot_cer3_gag(a,l,m): 
# 10/27/2020 make a VSG plot for various Cer3 alleles, developed for siRNA suppression project
# 
# a: Cer3 allele name, l: library name, m: True or False for plotting MM and distinguishing ambiguous and definitive generator reads
# ins: True or False for highlighting the insertion
# assumes 4nt barcodes
    fa_extension='_3linker_removed_collapsed.fa'
    fa_folder=dir_path(l+fa_extension,fastq_folder)
    fa_file_name=fa_folder+l+fa_extension
    barcode_size=4

    # VSG options (set for figure making for illustrator):
    ext_nt = 4 # number (nt) of extension when considering if the two tracks overlap, aka, minimal horizontal gap
    L_per_nt=0.12 # number of pixel per nt
    line_gap=0.5 #
    gDNA_track_width = 1.5
    track_width=1.3 # track width for sRNA

    bowtie_ref=bowtie_folder+'indexes/'+a
    bowtie_output=l+'_'+a+'_perfect_align.map'

    (sRNA_size_min, sRNA_size_max)=(20,24) #inclusive for both 
    # cer3 alleles information, numbers are 1-based offset
    cer3_dict={'cer3_red20_25kb':{'gag_start':7000, 'oma1_start': 7700, 'oma1_end': 8118,'gag_end':8118+700,
                                  'MM': [{'pos':7709, 'WT':'T', 'mut':'C'}, {'pos':7739, 'WT':'C', 'mut':'G'},
                                         {'pos':7769, 'WT':'G', 'mut':'A'}, {'pos':7799, 'WT':'A', 'mut':'G'},
                                         {'pos':7829, 'WT':'T', 'mut':'C'}, {'pos':7859, 'WT':'G', 'mut':'A'},
                                         {'pos':7889, 'WT':'C', 'mut':'T'}, {'pos':7919, 'WT':'A', 'mut':'G'},
                                         {'pos':7952, 'WT':'A', 'mut':'C'}, {'pos':7979, 'WT':'G', 'mut':'A'},
                                         {'pos':8009, 'WT':'C', 'mut':'T'}, {'pos':8039, 'WT':'T', 'mut':'C'},
                                         {'pos':8069, 'WT':'T', 'mut':'C'}, {'pos':8099, 'WT':'T', 'mut':'C'}]},
               'cer3_red40_25kb':{'gag_start':7000, 'oma1_start': 7700, 'oma1_end': 8118,'gag_end':8118+700,
                                  'MM': [{'pos':7719, 'WT':'A', 'mut':'G'}, {'pos':7749, 'WT':'A', 'mut':'G'},
                                     {'pos':7779, 'WT':'A', 'mut':'G'}, {'pos':7809, 'WT':'G', 'mut':'A'},
                                     {'pos':7839, 'WT':'C', 'mut':'T'}, {'pos':7869, 'WT':'T', 'mut':'G'},
                                     {'pos':7899, 'WT':'T', 'mut':'C'}, {'pos':7929, 'WT':'G', 'mut':'A'},
                                     {'pos':7959, 'WT':'C', 'mut':'T'}, {'pos':7989, 'WT':'A', 'mut':'G'},
                                     {'pos':8019, 'WT':'T', 'mut':'C'}, {'pos':8049, 'WT':'C', 'mut':'T'},
                                     {'pos':8079, 'WT':'G', 'mut':'C'}, {'pos':8109, 'WT':'A', 'mut':'G'}]},
                'cer3_red47_25kb':{'gag_start':7000, 'oma1_start': 7700, 'oma1_end': 8164,'gag_end':8164+700,
                                  'MM': [{'pos':7729, 'WT':'A', 'mut':'G'}, {'pos':7759, 'WT':'T', 'mut':'C'},
                                         {'pos':7789, 'WT':'A', 'mut':'C'}, {'pos':7819, 'WT':'A', 'mut':'G'},
                                         {'pos':7849, 'WT':'T', 'mut':'C'}, {'pos':7879, 'WT':'A', 'mut':'G'},
                                         {'pos':7909, 'WT':'A', 'mut':'G'}, {'pos':7939, 'WT':'T', 'mut':'C'},
                                         {'pos':7969, 'WT':'T', 'mut':'C'}, {'pos':7999, 'WT':'A', 'mut':'G'},
                                         {'pos':8029, 'WT':'C', 'mut':'A'}, {'pos':8059, 'WT':'A', 'mut':'G'},
                                         {'pos':8089, 'WT':'A', 'mut':'G'}, {'pos':8119, 'WT':'C', 'mut':'T'}]},
               'cer3_red49_25kb':{'gag_start':7000, 'oma1_start': 7700, 'oma1_end': 8337,'gag_end':8337+700,
                                  'MM':[],
                              'plot_flank': 800},
               
               'cer3_red52_25kb':{'gag_start':7000, 'oma1_start': 7700, 'oma1_end': 8168,'gag_end':8168+700,
                                  'MM': [{'pos':7717, 'WT':'T', 'mut':'A'}, {'pos':7747, 'WT':'T', 'mut':'A'}, 
                                         {'pos':7777, 'WT':'C', 'mut':'G'}, {'pos':7807, 'WT':'G', 'mut':'C'},
                                         {'pos':7837, 'WT':'T', 'mut':'A'}, {'pos':7867, 'WT':'G', 'mut':'C'},
                                         {'pos':7897, 'WT':'G', 'mut':'C'}, {'pos':7927, 'WT':'A', 'mut':'T'},
                                         {'pos':7957, 'WT':'G', 'mut':'C'}, {'pos':7987, 'WT':'A', 'mut':'T'},
                                         {'pos':8017, 'WT':'G', 'mut':'C'}, {'pos':8047, 'WT':'C', 'mut':'G'},
                                         {'pos':8077, 'WT':'C', 'mut':'G'}, {'pos':8107, 'WT':'C', 'mut':'G'},
                                         {'pos':8137, 'WT':'A', 'mut':'T'}, {'pos':8167, 'WT':'G', 'mut':'C'}]},
               'cer3_red66_25kb':{'gag_start':7000, 'oma1_start': 7703, 'oma1_end': 8649,'gag_end':8649+700},
               'cer3_wt_25kb':{'gag_start':7000, 'oma1_start': 7463, 'oma1_end': 7929,'gag_end':8400, # oma1_start and oma1_end mark the positions of cer3 sequence inserted in oma-1::cer3(red57)
                                  'MM': [{'pos':7472, 'WT':'A', 'mut':'G'}, {'pos':7502, 'WT':'G', 'mut':'A'}, # WT: red57, mut: wt cer3
                                         {'pos':7532, 'WT':'T', 'mut':'C'}, {'pos':7562, 'WT':'A', 'mut':'G'},
                                         {'pos':7592, 'WT':'T', 'mut':'C'}, {'pos':7622, 'WT':'A', 'mut':'G'},
                                         {'pos':7652, 'WT':'T', 'mut':'C'}, {'pos':7682, 'WT':'A', 'mut':'G'},
                                         {'pos':7712, 'WT':'C', 'mut':'T'}, {'pos':7742, 'WT':'G', 'mut':'A'},
                                         {'pos':7772, 'WT':'A', 'mut':'G'}, {'pos':7802, 'WT':'G', 'mut':'A'},
                                         {'pos':7832, 'WT':'G', 'mut':'A'}, {'pos':7862, 'WT':'G', 'mut':'A'},
                                         {'pos':7892, 'WT':'A', 'mut':'G'}, {'pos':7922, 'WT':'A', 'mut':'G'}]},
               'cer8_red35_15kb':{'gag_start':5507-700, 'oma1_start': 5507, 'oma1_end': 5928,'gag_end':5928+700, # gag_start and end here do not mark gag, but rather plotting area
                              'MM': [{'pos':5519, 'WT':'T', 'mut':'C'}, {'pos':5549, 'WT':'C', 'mut':'G'},
                                     {'pos':5579, 'WT':'G', 'mut':'A'}, {'pos':5609, 'WT':'A', 'mut':'G'},
                                     {'pos':5639, 'WT':'T', 'mut':'C'}, {'pos':5669, 'WT':'G', 'mut':'A'},
                                     {'pos':5699, 'WT':'C', 'mut':'T'}, {'pos':5729, 'WT':'A', 'mut':'G'},
                                     {'pos':5762, 'WT':'A', 'mut':'C'}, {'pos':5789, 'WT':'G', 'mut':'A'},
                                     {'pos':5819, 'WT':'C', 'mut':'T'}, {'pos':5849, 'WT':'T', 'mut':'C'},
                                     {'pos':5879, 'WT':'T', 'mut':'C'}, {'pos':5909, 'WT':'T', 'mut':'C'}]},
               'cer3_red118_25kb':{'gag_start':7000, 'oma1_start': 7700, 'oma1_end': 8082,'gag_end':8082+700,
                                  'MM': [ {'pos':7723, 'WT':'C', 'mut':'G'},
                                         {'pos':7753, 'WT':'G', 'mut':'A'}, {'pos':7783, 'WT':'A', 'mut':'G'},
                                         {'pos':7813, 'WT':'T', 'mut':'C'}, {'pos':7843, 'WT':'G', 'mut':'A'},
                                         {'pos':7873, 'WT':'C', 'mut':'T'}, {'pos':7903, 'WT':'A', 'mut':'G'},
                                         {'pos':7916, 'WT':'T', 'mut':'C'}, {'pos':7943, 'WT':'G', 'mut':'A'},
                                         {'pos':7973, 'WT':'C', 'mut':'T'}, {'pos':8003, 'WT':'T', 'mut':'C'},
                                         {'pos':8033, 'WT':'T', 'mut':'C'}, {'pos':8063, 'WT':'T', 'mut':'C'}]},
                'cer3_red122_25kb':{'gag_start':7000, 'oma1_start': 7700, 'oma1_end': 7700 + 399,'gag_end':7700 + 399 +700, 
                          'MM': [{'pos':7700+29, 'WT':'C', 'mut':'G'},{'pos':7700+59, 'WT':'T', 'mut':'A'},{'pos':7700+89, 'WT':'T', 'mut':'A'},{'pos':7700+119, 'WT':'A', 'mut':'T'},
                          {'pos':7700+149, 'WT':'A', 'mut':'T'},{'pos':7700+179, 'WT':'T', 'mut':'A'},{'pos':7700+209, 'WT':'A', 'mut':'T'},{'pos':7700+239, 'WT':'G', 'mut':'C'},
                          {'pos':7700+269, 'WT':'G', 'mut':'C'},{'pos':7700+299, 'WT':'A', 'mut':'T'},{'pos':7700+329, 'WT':'T', 'mut':'A'},{'pos':7700+359, 'WT':'T', 'mut':'A'},
                          {'pos':7700+389, 'WT':'T', 'mut':'A'}]},
                'cer3_red123_25kb':{'gag_start':7000, 'oma1_start': 7700, 'oma1_end': 7700 + 399,'gag_end':7700 + 399 +700, 
                          'MM': [{'pos':7700+29, 'WT':'G', 'mut':'C'},{'pos':7700+59, 'WT':'A', 'mut':'T'},{'pos':7700+89, 'WT':'A', 'mut':'T'},{'pos':7700+119, 'WT':'A', 'mut':'T'},
                          {'pos':7700+149, 'WT':'C', 'mut':'G'},{'pos':7700+179, 'WT':'G', 'mut':'C'},{'pos':7700+209, 'WT':'G', 'mut':'C'},{'pos':7700+239, 'WT':'A', 'mut':'T'},
                          {'pos':7700+269, 'WT':'A', 'mut':'T'},{'pos':7700+299, 'WT':'C', 'mut':'G'},{'pos':7700+329, 'WT':'T', 'mut':'A'},{'pos':7700+359, 'WT':'G', 'mut':'C'},
                          {'pos':7700+389, 'WT':'G', 'mut':'C'}]},
                'cer3_red124_25kb':{'gag_start':7000, 'oma1_start': 7700, 'oma1_end': 7700 + 399,'gag_end':7700 + 399 +700, 
                          'MM': [{'pos':7700+29, 'WT':'C', 'mut':'G'},{'pos':7700+59, 'WT':'C', 'mut':'G'},{'pos':7700+89, 'WT':'G', 'mut':'C'},{'pos':7700+119, 'WT':'T', 'mut':'A'},
                          {'pos':7700+149, 'WT':'T', 'mut':'A'},{'pos':7700+179, 'WT':'A', 'mut':'T'},{'pos':7700+209, 'WT':'C', 'mut':'G'},{'pos':7700+239, 'WT':'C', 'mut':'G'},
                          {'pos':7700+269, 'WT':'A', 'mut':'T'},{'pos':7700+299, 'WT':'C', 'mut':'G'},{'pos':7700+329, 'WT':'T', 'mut':'A'},{'pos':7700+359, 'WT':'T', 'mut':'A'},
                          {'pos':7700+389, 'WT':'T', 'mut':'A'}]},
                'cer3_red125_25kb':{'gag_start':7000, 'oma1_start': 7700, 'oma1_end': 7700 + 399,'gag_end':7700 + 399 +700, 
                          'MM': [{'pos':7700+29, 'WT':'T', 'mut':'A'},{'pos':7700+59, 'WT':'A', 'mut':'T'},{'pos':7700+89, 'WT':'G', 'mut':'C'},{'pos':7700+119, 'WT':'T', 'mut':'A'},
                          {'pos':7700+149, 'WT':'C', 'mut':'G'},{'pos':7700+179, 'WT':'G', 'mut':'C'},{'pos':7700+209, 'WT':'T', 'mut':'A'},{'pos':7700+239, 'WT':'T', 'mut':'A'},
                          {'pos':7700+269, 'WT':'A', 'mut':'T'},{'pos':7700+299, 'WT':'A', 'mut':'T'},{'pos':7700+329, 'WT':'A', 'mut':'T'},{'pos':7700+359, 'WT':'C', 'mut':'G'},
                          {'pos':7700+389, 'WT':'C', 'mut':'G'}]},
                'cer3_red126_25kb':{'gag_start':7000, 'oma1_start': 7700, 'oma1_end': 7700 + 398,'gag_end':7700 + 398 +700, 
                          'MM': [{'pos':7700+29, 'WT':'T', 'mut':'A'},{'pos':7700+59, 'WT':'A', 'mut':'T'},{'pos':7700+89, 'WT':'T', 'mut':'A'},{'pos':7700+119, 'WT':'T', 'mut':'A'},
                          {'pos':7700+149, 'WT':'T', 'mut':'A'},{'pos':7700+179, 'WT':'A', 'mut':'T'},{'pos':7700+209, 'WT':'G', 'mut':'C'},{'pos':7700+239, 'WT':'A', 'mut':'T'},
                          {'pos':7700+269, 'WT':'A', 'mut':'T'},{'pos':7700+299, 'WT':'C', 'mut':'G'},{'pos':7700+329, 'WT':'T', 'mut':'A'},{'pos':7700+359, 'WT':'T', 'mut':'A'},
                          {'pos':7700+389, 'WT':'C', 'mut':'G'}]},
                'cer3_red127_25kb':{'gag_start':7000, 'oma1_start': 7700, 'oma1_end': 7700 + 399,'gag_end':7700 + 399 +700, 
                          'MM': [{'pos':7700+29, 'WT':'T', 'mut':'A'},{'pos':7700+59, 'WT':'T', 'mut':'A'},{'pos':7700+89, 'WT':'A', 'mut':'T'},{'pos':7700+119, 'WT':'C', 'mut':'G'},
                          {'pos':7700+149, 'WT':'C', 'mut':'G'},{'pos':7700+179, 'WT':'A', 'mut':'T'},{'pos':7700+209, 'WT':'A', 'mut':'T'},{'pos':7700+239, 'WT':'C', 'mut':'G'},
                          {'pos':7700+269, 'WT':'C', 'mut':'G'},{'pos':7700+299, 'WT':'T', 'mut':'A'},{'pos':7700+329, 'WT':'T', 'mut':'A'},{'pos':7700+359, 'WT':'G', 'mut':'C'},
                          {'pos':7700+388, 'WT':'A', 'mut':'T'}]},
                'cer3_red129_25kb':{'gag_start':7000, 'oma1_start': 7700, 'oma1_end': 7700 + 398,'gag_end':7700 + 398 +700, 
                          'MM': [{'pos':7700+29, 'WT':'A', 'mut':'T'},{'pos':7700+59, 'WT':'T', 'mut':'A'},{'pos':7700+89, 'WT':'G', 'mut':'C'},{'pos':7700+119, 'WT':'T', 'mut':'A'},
                          {'pos':7700+149, 'WT':'A', 'mut':'T'},{'pos':7700+179, 'WT':'G', 'mut':'C'},{'pos':7700+209, 'WT':'T', 'mut':'A'},{'pos':7700+239, 'WT':'G', 'mut':'C'},
                          {'pos':7700+269, 'WT':'A', 'mut':'T'},{'pos':7700+299, 'WT':'G', 'mut':'C'},{'pos':7700+329, 'WT':'T', 'mut':'A'},{'pos':7700+359, 'WT':'A', 'mut':'T'},
                          {'pos':7700+389, 'WT':'C', 'mut':'G'}]},
                'cer3_red130_25kb':{'gag_start':7000, 'oma1_start': 7700, 'oma1_end': 7700 + 399,'gag_end':7700 + 399 +700, 
                          'MM': [{'pos':7700+29, 'WT':'G', 'mut':'C'},{'pos':7700+59, 'WT':'G', 'mut':'C'},{'pos':7700+89, 'WT':'C', 'mut':'G'},{'pos':7700+119, 'WT':'A', 'mut':'T'},
                          {'pos':7700+149, 'WT':'C', 'mut':'G'},{'pos':7700+179, 'WT':'C', 'mut':'G'},{'pos':7700+209, 'WT':'T', 'mut':'A'},{'pos':7700+239, 'WT':'T', 'mut':'A'},
                          {'pos':7700+269, 'WT':'G', 'mut':'C'},{'pos':7700+299, 'WT':'T', 'mut':'A'},{'pos':7700+329, 'WT':'A', 'mut':'T'},{'pos':7700+359, 'WT':'T', 'mut':'A'},
                          {'pos':7700+389, 'WT':'T', 'mut':'A'}]},
                'cer3_red131_25kb':{'gag_start':7000, 'oma1_start': 7700, 'oma1_end': 7700 + 393,'gag_end':7700 + 393 +700, 
                          'MM': [{'pos':7700+29, 'WT':'T', 'mut':'A'},{'pos':7700+59, 'WT':'T', 'mut':'A'},{'pos':7700+89, 'WT':'A', 'mut':'T'},{'pos':7700+119, 'WT':'A', 'mut':'T'},
                          {'pos':7700+149, 'WT':'A', 'mut':'T'},{'pos':7700+179, 'WT':'C', 'mut':'G'},{'pos':7700+209, 'WT':'A', 'mut':'T'},{'pos':7700+239, 'WT':'T', 'mut':'A'},
                          {'pos':7700+269, 'WT':'T', 'mut':'A'},{'pos':7700+299, 'WT':'C', 'mut':'G'},{'pos':7700+329, 'WT':'A', 'mut':'T'},{'pos':7700+359, 'WT':'G', 'mut':'C'},
                          {'pos':7700+389, 'WT':'C', 'mut':'G'}]},
                'cer3_red132a_25kb':{'gag_start':7000, 'oma1_start': 7700, 'oma1_end': 7700 + 399,'gag_end':7700 + 399 +700, 
                          'MM': [{'pos':7700+29, 'WT':'C', 'mut':'G'},{'pos':7700+59, 'WT':'C', 'mut':'G'},{'pos':7700+89, 'WT':'C', 'mut':'G'},{'pos':7700+119, 'WT':'T', 'mut':'A'},
                          {'pos':7700+149, 'WT':'T', 'mut':'A'},{'pos':7700+179, 'WT':'T', 'mut':'A'},{'pos':7700+209, 'WT':'A', 'mut':'T'},{'pos':7924, 'WT':'C', 'mut':'T'},{'pos':7700+239, 'WT':'A', 'mut':'T'},
                          {'pos':7700+269, 'WT':'A', 'mut':'T'},{'pos':7700+299, 'WT':'A', 'mut':'T'},{'pos':7700+329, 'WT':'C', 'mut':'G'},{'pos':7700+359, 'WT':'C', 'mut':'G'},
                          {'pos':7700+389, 'WT':'A', 'mut':'T'}]},
                'cer3_red134_25kb':{'gag_start':7000, 'oma1_start': 7700, 'oma1_end': 7700 + 1223,'gag_end':7700 + 1223 +700, 
                          'MM': [{'pos':7700+29, 'WT':'G', 'mut':'C'},{'pos':7700+59, 'WT':'C', 'mut':'G'},{'pos':7700+89, 'WT':'A', 'mut':'T'},{'pos':7700+119, 'WT':'T', 'mut':'A'},
                          {'pos':7700+149, 'WT':'T', 'mut':'A'},{'pos':7700+179, 'WT':'A', 'mut':'T'},{'pos':7700+209, 'WT':'C', 'mut':'G'},{'pos':7700+239, 'WT':'G', 'mut':'C'},
                          {'pos':7700+269, 'WT':'T', 'mut':'A'},{'pos':7700+299, 'WT':'A', 'mut':'T'},{'pos':7700+329, 'WT':'A', 'mut':'T'},{'pos':7700+359, 'WT':'T', 'mut':'A'},
                          {'pos':7700+389, 'WT':'A', 'mut':'T'},{'pos':7700+419, 'WT':'G', 'mut':'C'},{'pos':7700+449, 'WT':'G', 'mut':'C'},{'pos':7700+479, 'WT':'T', 'mut':'A'},
                          {'pos':7700+509, 'WT':'A', 'mut':'T'},{'pos':7700+539, 'WT':'C', 'mut':'G'},{'pos':7700+569, 'WT':'T', 'mut':'A'},{'pos':7700+599, 'WT':'T', 'mut':'A'},
                          {'pos':7700+629, 'WT':'G', 'mut':'C'},{'pos':7700+659, 'WT':'A', 'mut':'T'},{'pos':7700+689, 'WT':'T', 'mut':'A'},{'pos':7700+719, 'WT':'A', 'mut':'T'},
                          {'pos':7700+749, 'WT':'A', 'mut':'T'},{'pos':7700+779, 'WT':'A', 'mut':'T'},{'pos':7700+809, 'WT':'C', 'mut':'G'},{'pos':7700+839, 'WT':'T', 'mut':'A'},
                          {'pos':7700+869, 'WT':'T', 'mut':'A'},{'pos':7700+899, 'WT':'A', 'mut':'T'},{'pos':7700+929, 'WT':'T', 'mut':'A'},{'pos':7700+959, 'WT':'G', 'mut':'C'},
                          {'pos':7700+989, 'WT':'T', 'mut':'A'},{'pos':7700+1019, 'WT':'T', 'mut':'A'},{'pos':7700+1049, 'WT':'T', 'mut':'A'},{'pos':7700+1079, 'WT':'T', 'mut':'A'},
                          {'pos':7700+1109, 'WT':'T', 'mut':'A'},{'pos':7700+1139, 'WT':'G', 'mut':'C'},{'pos':7700+1169, 'WT':'G', 'mut':'C'},{'pos':7700+1199, 'WT':'G', 'mut':'C'},
                         ]},
               'cer3_red147a_25kb':{'gag_start':7000, 'oma1_start': 7700, 'oma1_end': 8091,'gag_end':8091+700,
                                  'MM':[],
                              'plot_flank': 800},
               'cer3_red146_25kb':{'gag_start':7000, 'oma1_start': 7700, 'oma1_end': 8101,'gag_end':8101+700,
                                  'MM':[],
                              'plot_flank': 800},
               'cer3_red145_25kb':{'gag_start':7000, 'oma1_start': 7700, 'oma1_end': 8103,'gag_end':8103+700,
                                  'MM':[],
                              'plot_flank': 800},
               'cer3_red153_25kb':{ 'gag_start':7000,'oma1_start': 7700, 'oma1_end': 8117, 'gag_end':700+8117,
                                  'MM':[],
                              'plot_flank': 800},
               'cer3_red154_25kb':{'gag_start':7000, 'oma1_start': 7700, 'oma1_end': 8090,'gag_end':700+8090,
                                  'MM':[],
                              'plot_flank': 800},
               'cer3_red152_25kb':{'gag_start':7000, 'oma1_start': 7700, 'oma1_end': 8119,'gag_end':700+8119,
                                  'MM':[],
                              'plot_flank': 800},
               'cer3_red155_25kb':{'gag_start':7000, 'oma1_start': 7700, 'oma1_end': 8090,'gag_end':700+8090,
                                  'MM':[],
                              'plot_flank': 800},
               'cer3_red166b_25kb':{'gag_start':7000,'oma1_start': 7700, 'oma1_end': 8089,'gag_end':700+8089,
                                  'MM':[],
                              'plot_flank': 800},
               'cer3_red167b_25kb':{'gag_start':7000, 'oma1_start': 7700, 'oma1_end': 7882,'gag_end':700+7882,
                                  'MM':[],
                              'plot_flank': 800},

               }
    
    bowtie_prog=bowtie_folder+'bowtie ' + bowtie_ref + ' -5 '+str(barcode_size)+' -f -t -v 0 -a -p '+str(core_num)+' '+fa_file_name+' '+bowtie_output
    p=subprocess.run(bowtie_prog, shell=True, check=True,capture_output=True)
    #print(p.stderr)
    
    # a dictionary for cer3 reads, distinguishing oma-1 reads (definitively generator-derived), oma-1 reads (ambiguous), cer3-only reads, flanking reads
    cer3_reads_dict={'oma1_gen':[],'oma1_ambi':[],'cer3only':[]}
    infile_bowtie_map=open(bowtie_output, 'r')
    c=0
    for L1 in infile_bowtie_map:
        c+=1
        A1=L1.strip().split('\t')
        (t_name, t_start, t_size, t_end, t_strand, t_else) =(A1[0], int(A1[3])+1, len(A1[4]), int(A1[3])+len(A1[4]), A1[1],int(A1[6]))
        (t_oma1_gen, t_oma1_ambi,t_cer3only,t_flanking) = (False,False,False,False)
        # check if the read is in the gag region, but not covering the insertion:
        if (t_start>=cer3_dict[a]['gag_start'] and t_end<cer3_dict[a]['oma1_start']) or (t_start>cer3_dict[a]['oma1_end'] and t_end<=cer3_dict[a]['gag_end']):
            cer3_reads_dict['cer3only'].append({'start': t_start,'end': t_end,'strand':t_strand,'name': t_name,'color':'black'})
        # check if the read covers a MM position or cer3:oma-1/oma-1:cer3 junction
        if m: # plot MM and distinguish ambigous and definitive generator reads
            if (t_end>=cer3_dict[a]['oma1_start'] and t_start<=cer3_dict[a]['oma1_end']):
                if ((cer3_dict[a]['oma1_start']-1) in range(t_start,t_end)) or ((cer3_dict[a]['oma1_end']+1) in range(t_start+1,t_end+1)): # label the junction reads
                    t_oma1_gen=True
                    cer3_reads_dict['oma1_gen'].append({'start': t_start,'end': t_end,'strand':t_strand,'name': t_name,'color':'orange'})
                else:
                    for M in cer3_dict[a]['MM']:
                        if M['pos'] in range(t_start,t_end+1):
                            t_oma1_gen=True
                            cer3_reads_dict['oma1_gen'].append({'start': t_start,'end': t_end,'strand':t_strand,'name': t_name,'color':'red'})
                            break
                if t_oma1_gen==False:
                    cer3_reads_dict['oma1_ambi'].append({'start': t_start,'end': t_end,'strand':t_strand,'name': t_name,'color':'gray'})    
        else:
            if (t_start>=cer3_dict[a]['oma1_start'] and t_end<=cer3_dict[a]['oma1_end']):
                t_oma1_gen=True
                cer3_reads_dict['oma1_gen'].append({'start': t_start,'end': t_end,'strand':t_strand,'name': t_name,'color':'red'})
            if ((cer3_dict[a]['oma1_start']-1) in range(t_start,t_end)) or ((cer3_dict[a]['oma1_end']+1) in range(t_start+1,t_end+1)): # label the junction reads
                t_oma1_gen=True
                cer3_reads_dict['oma1_gen'].append({'start': t_start,'end': t_end,'strand':t_strand,'name': t_name,'color':'orange'})
  
    infile_bowtie_map.close()
    rm_cmd='rm '+bowtie_output
    subprocess.run(rm_cmd, shell=True)

    # Calculate the number of siRNA in different categories.
    report_name='report_'+l+'_'+a+'_2.txt'
    report_out=open(report_name,'w')
    N_total_col=0 # total number of collapsed reads
    N_total_unCol=0 # total number of un-collapsed reads

    infile=open(fa_file_name,'r')
    for L1 in infile:
        if L1[0]=='>':
            N_total_col+=1
            N_total_unCol+=int(L1.strip().split('_')[-1])
    infile.close()
    report_out.write(l+'\n')
    report_out.write("Reference sequence for the alignment: %s" % a)
    report_out.write("\nTotal number of collapsed reads: %i" % N_total_col)
    report_out.write("\nTotal number of Uncollapsed reads: %i" % N_total_unCol)

    # align to whole ce6
    bowtie_ref=bowtie_folder+'indexes/ce6'
    N_ce6_col=0
    N_ce6_unCol=0
    bowtie_output=a+'_ce6_align.map'
    bowtie_prog=bowtie_folder+'bowtie ' + bowtie_ref + ' -5 '+str(barcode_size)+' -f -t -v 0 -a -p '+ str(core_num)+' '+fa_file_name+' '+bowtie_output
    p=subprocess.run(bowtie_prog, shell=True, check=True,capture_output=True)
    #print(p.stderr)
    
    infile_map=open(bowtie_output,'r')
    for L1 in infile_map:
        A1=L1.strip().split('\t')
        N_ce6_col+=1.0/(1+int(A1[6]))
        N_ce6_unCol+=int((A1[0].split('_'))[-1])/(1+int(A1[6]))
    report_out.write("\nTotal number of reads that perfectly match to ce6:")
    report_out.write("\nCollapsed reads: %i" % N_ce6_col)
    report_out.write("\nUncollapsed reads: %i" % N_ce6_unCol)
    infile_map.close()
    rm_cmd='rm '+bowtie_output
    subprocess.run(rm_cmd, shell=True)
    

    # Total number of oma-1 siRNAs with MM or cer3 junction
    (N_antisense, N_sense)=(0,0)
    for T in cer3_reads_dict['oma1_gen']:
        if ((T['end']-T['start']+1) in range(sRNA_size_min, sRNA_size_max+1)) :
            if T['strand']=='-': N_antisense+=1
            if T['strand']=='+': N_sense+=1

    report_out.write("\n Total number of oma-1 siRNAs with MM or cer3 junction, antisense: %i" % N_antisense)
    report_out.write("\n Total number of oma-1 siRNAs with MM or cer3 junction, sense: %i" % N_sense)

    # Total number of red20-matching oma-1 siRNAs without MM or cer3 junction
    (N_antisense, N_sense)=(0,0)
    for T in cer3_reads_dict['oma1_ambi']:
        if  ((T['end']-T['start']+1) in range(sRNA_size_min, sRNA_size_max+1)) :
            if T['strand']=='-': N_antisense+=1
            if T['strand']=='+': N_sense+=1
    report_out.write("\nTotal number of red40-matching oma-1 siRNAs without MM or cer3 junction, antisense: %i" % N_antisense)
    report_out.write("\nTotal number of red40-matching oma-1 siRNAs without MM or cer3 junction, sense: %i" % N_sense)
    report_out.close()
    
    # track dictionary: key: line (1, 2, etc), value: a list of tracks
    track_dict_sense = {1:[]}
    track_dict_antisense = {1:[]}
    for T_type in ['oma1_gen','oma1_ambi','cer3only']:
        for T1 in cer3_reads_dict[T_type]:
            if T1['strand'] == '+' and ((T1['end']-T1['start']+1) in range(sRNA_size_min, sRNA_size_max+1)) and (T1['start']>=cer3_dict[a]['gag_start'] and T1['end']<=cer3_dict[a]['gag_end']):
                overlap_flag='no'
                for track_line in range(1,len(track_dict_sense)+1): # check the latest track line
                    overlap_flag='no'
                    for T2 in track_dict_sense[track_line]:
                        # check if T1 and T2 overlap (extend start and end by ext_nt nt)
                        if T1['start'] in range(T2['start']-ext_nt, T2['end']+ext_nt) or T1['end'] in range(T2['start']-ext_nt, T2['end']+ext_nt):
                            overlap_flag='yes'
                            break
                        if T2['start'] in range(T1['start']-ext_nt, T1['end']+ext_nt) or T2['end'] in range(T1['start']-ext_nt, T1['end']+ext_nt) :
                            overlap_flag='yes'
                            break
                    if overlap_flag == 'no': # if no overlap with any of the existing tracks in the current line
                        track_dict_sense[track_line].append({'start': T1['start'],
                                                             'end': T1['end'],
                                                             'type':T_type,
                                                             'color': T1['color']})
                        break
                if overlap_flag == 'yes': # if it overlaps with a track in all existing line, add a new line, add the track to the new line
                    track_dict_sense[len(track_dict_sense)+1]=[{'start': T1['start'],
                                                             'end': T1['end'],
                                                             'type':T_type,
                                                             'color': T1['color']}]

            if T1['strand'] == '-' and ((T1['end']-T1['start']+1) in range(sRNA_size_min, sRNA_size_max+1)) and (T1['start']>=cer3_dict[a]['gag_start'] and T1['end']<=cer3_dict[a]['gag_end']):
                for track_line in range(1,len(track_dict_antisense)+1): # check the latest track line
                    overlap_flag='no'
                    for T2 in track_dict_antisense[track_line]:
                        overlap_flag='no'
                        # check if T1 and T2 overlap (extend start and end by ext_nt nt)
                        if T1['start'] in range(T2['start']-ext_nt, T2['end']+ext_nt) or T1['end'] in range(T2['start']-ext_nt, T2['end']+ext_nt):
                            overlap_flag='yes'
                            break
                        if T2['start'] in range(T1['start']-ext_nt, T1['end']+ext_nt) or T2['end'] in range(T1['start']-ext_nt, T1['end']+ext_nt) :
                            overlap_flag='yes'
                            break
                    if overlap_flag == 'no': # if no overlap with any of the existing tracks in the current line
                        track_dict_antisense[track_line].append({'start': T1['start'],
                                                             'end': T1['end'],
                                                             'type':T_type,
                                                             'color': T1['color']})
                        break
                if overlap_flag == 'yes': # if it overlaps with a track in all existing line, add a new line, add the track to the new line
                    track_dict_antisense[len(track_dict_antisense)+1]=[{'start': T1['start'],
                                                             'end': T1['end'],
                                                             'type':T_type,
                                                             'color': T1['color']}]


    # VSG plot
    VSG_py_tracks_name=l+'_'+a+'_zoom.py'
    VSG_plot_tracks_name=l+'_'+a+'_zoom.svg'
    outfile_VSG_py=open(VSG_py_tracks_name,'w')
    outfile_VSG_py.write('from VSG_Module import * \n')

    (plot_start,plot_end)= (cer3_dict[a]['gag_start'], cer3_dict[a]['gag_end']) 
    (oma1_start,oma1_end)=(cer3_dict[a]['oma1_start'],cer3_dict[a]['oma1_end'])

    for track_line in range(len(track_dict_sense),0,-1):
        #print track_line
        VSG_y=-1.0*(len(track_dict_sense)-track_line)*(track_width+line_gap)
        for T in track_dict_sense[track_line]:
            (T_start, T_end)=(T['start'],T['end'])
            if T_start in range(plot_start,plot_end+1) and T_end in range(plot_start,plot_end+1) :
                outfile_VSG_py.write('vline(x1='+str(T_start*L_per_nt)+', x2='+str(T_end*L_per_nt)+', y1='+str(VSG_y)+', y2='+str(VSG_y)+', strokewidth='+str(track_width)+' ,stroke='+T['color']+')\n')

    # plot the gDNA b/t plot_start and plot_end
    gDNA_xc=(plot_end+plot_start)/2*L_per_nt
    gDNA_yc=-1.0*(len(track_dict_sense))*(track_width+line_gap) - 0.5*gDNA_track_width
    gDNA_width=(plot_end-plot_start+1)*L_per_nt
    outfile_VSG_py.write("vrect(xc="+str(gDNA_xc)+", yc="+str(gDNA_yc)+", width="+str(gDNA_width)+", height="+ str(gDNA_track_width) +", fill=grey,stroke=grey )\n")

    # plot oma1 insert
    gDNA_xc=(oma1_start+oma1_end)/2*L_per_nt
    gDNA_yc=-1.0*(len(track_dict_sense))*(track_width+line_gap) - 0.5*gDNA_track_width
    gDNA_width=(oma1_end-oma1_start+1)*L_per_nt
    outfile_VSG_py.write("vrect(xc="+str(gDNA_xc)+", yc="+str(gDNA_yc)+", width="+str(gDNA_width)+", height="+ str(gDNA_track_width) +", fill=purple,stroke=purple )\n")

    # plot the mismatches:
    if m:
        for M in cer3_dict[a]['MM']:
            (L_x, L_y1, L_y2)=(M['pos']*L_per_nt,-1.0*(len(track_dict_sense))*(track_width+line_gap) - 1.0*gDNA_track_width,-1.0*(len(track_dict_sense))*(track_width+line_gap))
            outfile_VSG_py.write("vline(x1="+str(L_x)+", y1="+str(L_y1)+", x2="+str(L_x)+", y2="+ str(L_y2) +", strokewidth=0.5, stroke=white )\n")

    for track_line in range(1,len(track_dict_antisense)+1):
        #print track_line
        VSG_y=-1.0*(len(track_dict_sense)+track_line)*(track_width+line_gap) - 1.5*gDNA_track_width
        for T in track_dict_antisense[track_line]:
            (T_start, T_end)=(T['start'],T['end'])
            if T_start in range(plot_start,plot_end+1) and T_end in range(plot_start,plot_end+1) :
                outfile_VSG_py.write('vline(x1='+str(T_start*L_per_nt)+', x2='+str(T_end*L_per_nt)+', y1='+str(VSG_y)+', y2='+str(VSG_y)+', strokewidth='+str(track_width)+' ,stroke='+T['color']+')\n')


    outfile_VSG_py.write("vset(bg='white')\n")
    outfile_VSG_py.write("vdisplay('"+VSG_plot_tracks_name+"')")
    outfile_VSG_py.close()

    vsg_cmd='python '+VSG_py_tracks_name
    subprocess.run(vsg_cmd,shell=True)
    
def siRNA_track_plot_gene_102920(a,l,m):
    # VSG options (set for figure making for illustrator):
    # a:  allele name, l: library name, m: True or False for plotting MM and distinguishing ambiguous and definitive generator reads

    ext_nt = 4 # number (nt) of extension when considering if the two tracks overlap, aka, minimal horizontal gap
    L_per_nt=0.12 # number of pixel per nt
    line_gap=0.5 #
    intron_track_width = 0.5 #
    exon_track_width = 1.7
    UTR_track_width = 1.3
    track_width=1.3 # track width for sRNA

    bowtie_ref=bowtie_folder+'indexes/'+a
    bowtie_output=l+'_'+a+'_perfect_align.map'

    (sRNA_size_min, sRNA_size_max)=(20,24) # inclusive for both

    fa_folder=dir_path(l+'_3linker_removed_collapsed.fa',fastq_folder)
    fa_file_name=fa_folder+l+'_3linker_removed_collapsed.fa'
    barcode_size=4

    # plot oma-1 tracks
    # plot area: chrIV:8,888,758-8,891181 (1-based and inclusive)
    allele_dict={'oma1_wt_102920':{'tx_start':201, 'tx_end': 2224, 'exon_starts': [311,397,611,884,1076,1351],'exon_ends':[348,564,841,1031,1307,1757], # these numbers are all 1-based and inclusive, the exons are all coding sequences
                                   'highlight_starts':[1129,1351],'highlight_ends':[1307,1590],'highlight_type':'exon',
                                  'MM': [{'pos':1138, 'WT':'C', 'mut':'T'}, {'pos':1168, 'WT':'G', 'mut':'C'}, # 'mut': nt used in WT oma-1, WT" nt used in Cer3:oma-1 ASG (red20)
                                         {'pos':1198, 'WT':'A', 'mut':'G'}, {'pos':1228, 'WT':'G', 'mut':'A'},
                                         {'pos':1258, 'WT':'C', 'mut':'T'}, {'pos':1288, 'WT':'A', 'mut':'G'},
                                         {'pos':1361, 'WT':'T', 'mut':'C'}, {'pos':1391, 'WT':'G', 'mut':'A'},
                                         {'pos':1421, 'WT':'C', 'mut':'A'}, {'pos':1451, 'WT':'A', 'mut':'G'},
                                         {'pos':1481, 'WT':'T', 'mut':'C'}, {'pos':1511, 'WT':'C', 'mut':'T'},
                                         {'pos':1541, 'WT':'C', 'mut':'T'}, {'pos':1571, 'WT':'C', 'mut':'T'}]},
                 'oma1_red57_102920':{'tx_start':201, 'tx_end': 2694, 'exon_starts': [311,397,611,884,1076,1351],'exon_ends':[348,564,841,1031,1307,1751], # these numbers are all 1-based and inclusive
                                      'highlight_starts':[1752],'highlight_ends':[2218], 'highlight_type':'UTR',
                                  'MM': [{'pos':1761, 'WT':'G', 'mut':'A'}, {'pos':1791, 'WT':'A', 'mut':'G'}, # WT: wt Cer3, mut: red57
                                         {'pos':1821, 'WT':'C', 'mut':'T'}, {'pos':1851, 'WT':'G', 'mut':'A'},
                                         {'pos':1881, 'WT':'C', 'mut':'T'}, {'pos':1911, 'WT':'G', 'mut':'A'},
                                         {'pos':1941, 'WT':'C', 'mut':'T'}, {'pos':1971, 'WT':'G', 'mut':'A'},
                                         {'pos':2001, 'WT':'T', 'mut':'C'}, {'pos':2031, 'WT':'A', 'mut':'G'},
                                         {'pos':2061, 'WT':'G', 'mut':'A'}, {'pos':2091, 'WT':'A', 'mut':'G'},
                                         {'pos':2121, 'WT':'A', 'mut':'G'}, {'pos':2151, 'WT':'A', 'mut':'G'},
                                         {'pos':2181, 'WT':'G', 'mut':'A'}, {'pos':2211, 'WT':'G', 'mut':'A'}]},
                 'oma1_red44_102920':{'tx_start':201, 'tx_end': 3171, 'exon_starts': [311,397,611,884,1076,1351,1841,2109,2255,2483],'exon_ends':[348,564,841,1031,1307,1796,2057,2203,2431,2704], # these numbers are all 1-based and inclusive
                                      'highlight_starts':[2034,2109,2255,2483],'highlight_ends':[2057,2203,2431,2655], 'highlight_type':'exon',
                                  'MM': [{'pos':2051, 'WT':'A', 'mut':'T'}, {'pos':2132, 'WT':'A', 'mut':'T'}, 
                                         {'pos':2162, 'WT':'G', 'mut':'C'}, {'pos':2192, 'WT':'C', 'mut':'G'},
                                         {'pos':2273, 'WT':'A', 'mut':'T'}, {'pos':2303, 'WT':'C', 'mut':'G'},
                                         {'pos':2333, 'WT':'C', 'mut':'G'}, {'pos':2363, 'WT':'T', 'mut':'A'},
                                         {'pos':2393, 'WT':'C', 'mut':'G'}, {'pos':2423, 'WT':'T', 'mut':'A'},
                                         {'pos':2504, 'WT':'C', 'mut':'G'}, {'pos':2534, 'WT':'G', 'mut':'C'},
                                         {'pos':2564, 'WT':'G', 'mut':'C'}, {'pos':2594, 'WT':'G', 'mut':'C'},
                                         {'pos':2624, 'WT':'T', 'mut':'A'}, {'pos':2654, 'WT':'C', 'mut':'G'}]},
                 'meg2_wt_021621a':{'tx_start':184, 'tx_end':3670 , 'exon_starts': [205,397,981,2644,3032,3217],'exon_ends':[343,699,2600,2875,3172,3241], # these numbers are all 1-based and inclusive  #chrX:13,922,301-13,926,235
                                   'highlight_starts':[657,981],'highlight_ends':[699,1337],'highlight_type':'exon',
                                  'MM': [{'pos':686, 'WT':'C', 'mut':'G'}, {'pos':997, 'WT':'C', 'mut':'G'}, # 'mut': nt used in WT meg-2, WT" nt used in Cer3:meg-2 (red130)
                                         {'pos':1027, 'WT':'G', 'mut':'C'}, {'pos':1057, 'WT':'T', 'mut':'A'},
                                        {'pos':1087, 'WT':'G', 'mut':'C'}, {'pos':1117, 'WT':'G', 'mut':'C'},
                                        {'pos':1147, 'WT':'A', 'mut':'T'}, {'pos':1177, 'WT':'A', 'mut':'T'},
                                        {'pos':1207, 'WT':'C', 'mut':'G'}, {'pos':1237, 'WT':'A', 'mut':'T'},
                                        {'pos':1267, 'WT':'T', 'mut':'A'}, {'pos':1297, 'WT':'A', 'mut':'T'},
                                        {'pos':1327, 'WT':'A', 'mut':'T'}]},
                 'mes1_wt_080921':{'tx_start':24, 'tx_end':4417 , 'exon_starts': [24,202,340,474,932,1318,1528,2504,2862,3092,3411,4291],'exon_ends':[143,285,439,717,1249,1480,1859,2809,3044,3332,4113,4417], # these numbers are all 1-based and inclusive  
                                   'highlight_starts':[2953,3092],'highlight_ends':[3044,3328],'highlight_type':'exon',
                                  'MM': [{'pos':2982, 'WT':'T', 'mut':'A'},{'pos':3012, 'WT':'A', 'mut':'T'},
                                        {'pos':3119, 'WT':'A', 'mut':'T'},{'pos':3149, 'WT':'T', 'mut':'A'},
                                         {'pos':3179, 'WT':'C', 'mut':'G'},{'pos':3209, 'WT':'A', 'mut':'T'},
                                         {'pos':3239, 'WT':'C', 'mut':'G'},{'pos':3269, 'WT':'T', 'mut':'A'},
                                         {'pos':3299, 'WT':'C', 'mut':'G'}]}, # WT; NT USED IN Cer3::mes1, mut: used in wt mes-1
    
                }

    bowtie_prog=bowtie_folder+'bowtie ' + bowtie_ref + ' -5 '+str(barcode_size)+' -f -t -v 0 -a -p '+str(core_num)+' '+fa_file_name+' '+bowtie_output
    p=subprocess.run(bowtie_prog, shell=True, check=True,capture_output=True)
    #print(p.stderr)

    # a dictionary for reads, distinguishing featured reads (definitively gene-derived), ambiguous, non-featured reads
    # if m==False, all reads go to the non-featured
    reads_dict={'feature_definitive':[],'feature_ambi':[],'non_feature':[]}
    infile_bowtie_map=open(bowtie_output, 'r')
    c=0
    for L1 in infile_bowtie_map:
        c+=1
        A1=L1.strip().split('\t')
        (t_name, t_start, t_size, t_end, t_strand, t_else) =(A1[0], int(A1[3])+1, len(A1[4]), int(A1[3])+len(A1[4]), A1[1],int(A1[6]))
        (t_feature_definitive, t_feature_ambi,non_feature) = (False,False,False)

        # check if the read covers a MM position or feature/non-feature junction
        if m: # plot MM and distinguish ambigous and definitive generator reads
            # check if the read is outside of the featured region
            if (t_start>=allele_dict[a]['tx_start'] and t_end<allele_dict[a]['highlight_starts'][0]) or (t_start>allele_dict[a]['highlight_ends'][-1] and t_end<=allele_dict[a]['tx_end']):
                reads_dict['non_feature'].append({'start': t_start,'end': t_end,'strand':t_strand,'name': t_name,'color':'black'})
                non_feature=True
            if (t_end>=allele_dict[a]['highlight_starts'][0] and t_start<=allele_dict[a]['highlight_ends'][-1]):
                if ((allele_dict[a]['highlight_starts'][0]-1) in range(t_start,t_end)) or ((allele_dict[a]['highlight_ends'][-1]+1) in range(t_start+1,t_end+1)): # label the junction reads
                    t_feature_definitive=True
                    reads_dict['feature_definitive'].append({'start': t_start,'end': t_end,'strand':t_strand,'name': t_name,'color':'orange'})
                else:
                    if 'MM' in allele_dict[a]:
                        for M in allele_dict[a]['MM']:
                            if M['pos'] in range(t_start,t_end+1):
                                t_feature_definitive=True
                                reads_dict['feature_definitive'].append({'start': t_start,'end': t_end,'strand':t_strand,'name': t_name,'color':'red'})
                                break
                if t_feature_definitive==False:
                    flag=True
                    for h_n in range(len(allele_dict[a]['highlight_starts'])):
                        if (t_start>=allele_dict[a]['highlight_starts'][h_n] and t_end<=allele_dict[a]['highlight_ends'][h_n]):
                            reads_dict['feature_ambi'].append({'start': t_start,'end': t_end,'strand':t_strand,'name': t_name,'color':'gray'}) 
                            flag=False
                            break
                    if flag:
                        reads_dict['non_feature'].append({'start': t_start,'end': t_end,'strand':t_strand,'name': t_name,'color':'black'}) 
        else:
            if (t_start>=allele_dict[a]['tx_start'] and t_end<=allele_dict[a]['tx_end']):
                reads_dict['non_feature'].append({'start': t_start,'end': t_end,'strand':t_strand,'name': t_name,'color':'black'})

    infile_bowtie_map.close()
    rm_cmd='rm '+bowtie_output
    subprocess.run(rm_cmd, shell=True)

    # track dictionary: key: line (1, 2, etc), value: a list of tracks
    track_dict_sense = {1:[]}
    track_dict_antisense = {1:[]}
    for T_type in ['feature_definitive','feature_ambi','non_feature']:
        for T1 in reads_dict[T_type]:
            if T1['strand'] == '+' and ((T1['end']-T1['start']+1) in range(sRNA_size_min, sRNA_size_max+1)) :
                overlap_flag='no'
                for track_line in range(1,len(track_dict_sense)+1): # check the latest track line
                    overlap_flag='no'
                    for T2 in track_dict_sense[track_line]:
                        # check if T1 and T2 overlap (extend start and end by ext_nt nt)
                        if T1['start'] in range(T2['start']-ext_nt, T2['end']+ext_nt) or T1['end'] in range(T2['start']-ext_nt, T2['end']+ext_nt):
                            overlap_flag='yes'
                            break
                        if T2['start'] in range(T1['start']-ext_nt, T1['end']+ext_nt) or T2['end'] in range(T1['start']-ext_nt, T1['end']+ext_nt) :
                            overlap_flag='yes'
                            break
                    if overlap_flag == 'no': # if no overlap with any of the existing tracks in the current line
                        track_dict_sense[track_line].append({'start': T1['start'],
                                                             'end': T1['end'],
                                                             'type':T_type,
                                                             'color': T1['color']})
                        break
                if overlap_flag == 'yes': # if it overlaps with a track in all existing line, add a new line, add the track to the new line
                    track_dict_sense[len(track_dict_sense)+1]=[{'start': T1['start'],
                                                             'end': T1['end'],
                                                             'type':T_type,
                                                             'color': T1['color']}]

            if T1['strand'] == '-' and ((T1['end']-T1['start']+1) in range(sRNA_size_min, sRNA_size_max+1)):
                for track_line in range(1,len(track_dict_antisense)+1): # check the latest track line
                    overlap_flag='no'
                    for T2 in track_dict_antisense[track_line]:
                        overlap_flag='no'
                        # check if T1 and T2 overlap (extend start and end by ext_nt nt)
                        if T1['start'] in range(T2['start']-ext_nt, T2['end']+ext_nt) or T1['end'] in range(T2['start']-ext_nt, T2['end']+ext_nt):
                            overlap_flag='yes'
                            break
                        if T2['start'] in range(T1['start']-ext_nt, T1['end']+ext_nt) or T2['end'] in range(T1['start']-ext_nt, T1['end']+ext_nt) :
                            overlap_flag='yes'
                            break
                    if overlap_flag == 'no': # if no overlap with any of the existing tracks in the current line
                        track_dict_antisense[track_line].append({'start': T1['start'],
                                                             'end': T1['end'],
                                                             'type':T_type,
                                                             'color': T1['color']})
                        break
                if overlap_flag == 'yes': # if it overlaps with a track in all existing line, add a new line, add the track to the new line
                    track_dict_antisense[len(track_dict_antisense)+1]=[{'start': T1['start'],
                                                             'end': T1['end'],
                                                             'type':T_type,
                                                             'color': T1['color']}]


    # VSG plot
    VSG_py_tracks_name=l+'_'+a+'.py'
    VSG_plot_tracks_name=l+'_'+a+'.svg'
    outfile_VSG_py=open(VSG_py_tracks_name,'w')
    outfile_VSG_py.write('from VSG_Module import * \n')

    (plot_start,plot_end)= (allele_dict[a]['tx_start'], allele_dict[a]['tx_end']) 

    for track_line in range(len(track_dict_sense),0,-1):
        #print track_line
        VSG_y=-1.0*(len(track_dict_sense)-track_line)*(track_width+line_gap)
        for T in track_dict_sense[track_line]:
            (T_start, T_end)=(T['start'],T['end'])
            if T_start in range(plot_start,plot_end+1) and T_end in range(plot_start,plot_end+1) :
                outfile_VSG_py.write('vline(x1='+str(T_start*L_per_nt)+', x2='+str(T_end*L_per_nt)+', y1='+str(VSG_y)+', y2='+str(VSG_y)+', strokewidth='+str(track_width)+' ,stroke='+T['color']+')\n')

    # plot the gDNA b/t plot_start and plot_end
    gDNA_xc=(plot_end+plot_start)/2*L_per_nt
    gDNA_yc=-1.0*(len(track_dict_sense))*(track_width+line_gap) - 0.5*exon_track_width
    gDNA_width=(plot_end-plot_start+1)*L_per_nt
    outfile_VSG_py.write("vrect(xc="+str(gDNA_xc)+", yc="+str(gDNA_yc)+", width="+str(gDNA_width)+", height="+ str(intron_track_width) +", fill=blue,stroke=blue )\n")
    # plot introns
    gDNA_xc=(allele_dict[a]['exon_starts'][0]+allele_dict[a]['exon_ends'][-1])/2*L_per_nt
    gDNA_yc=-1.0*(len(track_dict_sense))*(track_width+line_gap) - 0.5*exon_track_width
    gDNA_width=(allele_dict[a]['exon_ends'][-1]-allele_dict[a]['exon_starts'][0]+1)*L_per_nt
    outfile_VSG_py.write("vrect(xc="+str(gDNA_xc)+", yc="+str(gDNA_yc)+", width="+str(gDNA_width)+", height="+ str(intron_track_width) +", fill=gray,stroke=gray )\n")
    # plot exons:
    for e_n in range(len(allele_dict[a]['exon_starts'])):
        gDNA_xc=(allele_dict[a]['exon_starts'][e_n]+allele_dict[a]['exon_ends'][e_n])/2*L_per_nt
        gDNA_yc=-1.0*(len(track_dict_sense))*(track_width+line_gap) - 0.5*exon_track_width
        gDNA_width=(allele_dict[a]['exon_ends'][e_n]-allele_dict[a]['exon_starts'][e_n]+1)*L_per_nt
        outfile_VSG_py.write("vrect(xc="+str(gDNA_xc)+", yc="+str(gDNA_yc)+", width="+str(gDNA_width)+", height="+ str(exon_track_width) +", fill=blue,stroke=blue )\n")
    # plot highlight regions:
    if m:
        if allele_dict[a]['highlight_type']=='exon':
            feature_track_width=exon_track_width
        if allele_dict[a]['highlight_type']=='UTR':
            feature_track_width=intron_track_width
        for h_n in range(len(allele_dict[a]['highlight_starts'])):
            gDNA_xc=(allele_dict[a]['highlight_starts'][h_n]+allele_dict[a]['highlight_ends'][h_n])/2*L_per_nt
            gDNA_yc=-1.0*(len(track_dict_sense))*(track_width+line_gap) - 0.5*exon_track_width
            gDNA_width=(allele_dict[a]['highlight_ends'][h_n]-allele_dict[a]['highlight_starts'][h_n]+1)*L_per_nt
            outfile_VSG_py.write("vrect(xc="+str(gDNA_xc)+", yc="+str(gDNA_yc)+", width="+str(gDNA_width)+", height="+ str(feature_track_width) +", fill=purple,stroke=purple )\n")   
        if 'MM' in allele_dict[a]:
            for M in allele_dict[a]['MM']:
                (L_x, L_y1, L_y2)=(M['pos']*L_per_nt,-1.0*(len(track_dict_sense))*(track_width+line_gap) - 1.0*exon_track_width,-1.0*(len(track_dict_sense))*(track_width+line_gap))
                outfile_VSG_py.write("vline(x1="+str(L_x)+", y1="+str(L_y1)+", x2="+str(L_x)+", y2="+ str(L_y2) +", strokewidth=0.5, stroke=white )\n")

    for track_line in range(1,len(track_dict_antisense)+1):
        #print track_line
        VSG_y=-1.0*(len(track_dict_sense)+track_line)*(track_width+line_gap) - 1.5*exon_track_width
        for T in track_dict_antisense[track_line]:
            (T_start, T_end)=(T['start'],T['end'])
            if T_start in range(plot_start,plot_end+1) and T_end in range(plot_start,plot_end+1) :
                outfile_VSG_py.write('vline(x1='+str(T_start*L_per_nt)+', x2='+str(T_end*L_per_nt)+', y1='+str(VSG_y)+', y2='+str(VSG_y)+', strokewidth='+str(track_width)+' ,stroke='+T['color']+')\n')

    outfile_VSG_py.write("vset(bg='white')\n")
    outfile_VSG_py.write("vdisplay('"+VSG_plot_tracks_name+"')")
    outfile_VSG_py.close()

    vsg_cmd='python '+VSG_py_tracks_name
    subprocess.run(vsg_cmd,shell=True)  
    

###
def siRNA_track_plot_oma1_smg1_fusion(l): 
# 4/21/21 make a VSG plot for the oma-1::smg-1 transgene, 
# 
# l: library name,

# assumes 4nt barcodes
    a='oma1_smg1_fusion_mRNA'
    fa_extension='_3linker_removed_collapsed.fa'
    fa_folder=dir_path(l+fa_extension,fastq_folder)
    fa_file_name=fa_folder+l+fa_extension
    barcode_size=4

    # VSG options (set for figure making for illustrator):
    ext_nt = 4 # number (nt) of extension when considering if the two tracks overlap, aka, minimal horizontal gap
    L_per_nt=0.12 # number of pixel per nt
    line_gap=0.5 #
    gDNA_track_width = 1.5
    track_width=1.3 # track width for sRNA

    bowtie_ref=bowtie_folder+'indexes/'+a
    bowtie_output=l+'_'+a+'_perfect_align.map'

    (sRNA_size_min, sRNA_size_max)=(20,24) #inclusive for both 
    ref_dict={'oma1_smg1_fusion_mRNA':{'oma1_start': 189, 'oma1_end': 680, # start and end position of MM fragment
                                  'MM': [{'pos':191, 'WT':'A', 'mut':'G'}, {'pos':221, 'WT':'T', 'mut':'C'},
                                         {'pos':251, 'WT':'C', 'mut':'G'}, {'pos':281, 'WT':'T', 'mut':'A'},
                                         {'pos':311, 'WT':'T', 'mut':'C'}, {'pos':341, 'WT':'A', 'mut':'G'},
                                         {'pos':371, 'WT':'T', 'mut':'C'}, {'pos':401, 'WT':'T', 'mut':'C'},
                                         {'pos':431, 'WT':'A', 'mut':'G'}, {'pos':461, 'WT':'G', 'mut':'A'},
                                         {'pos':491, 'WT':'T', 'mut':'C'}, {'pos':521, 'WT':'G', 'mut':'A'},
                                         {'pos':551, 'WT':'A', 'mut':'G'}, {'pos':581, 'WT':'G', 'mut':'A'},
                                         {'pos':611, 'WT':'G', 'mut':'C'},
                                         {'pos':641, 'WT':'G', 'mut':'A'}, {'pos':671, 'WT':'A', 'mut':'G'}
                                        ]}}
    
    bowtie_prog=bowtie_folder+'bowtie ' + bowtie_ref + ' -5 '+str(barcode_size)+' -f -t -v 0 -a -p '+str(core_num)+' '+fa_file_name+' '+bowtie_output
    p=subprocess.run(bowtie_prog, shell=True, check=True,capture_output=True)
    #print(p.stderr)
    
    reads_dict={'oma1_MM':[],'ambi':[]}
    infile_bowtie_map=open(bowtie_output, 'r')
    c=0
    for L1 in infile_bowtie_map:
        c+=1
        A1=L1.strip().split('\t')
        (t_name, t_start, t_size, t_end, t_strand, t_else) =(A1[0], int(A1[3])+1, len(A1[4]), int(A1[3])+len(A1[4]), A1[1],int(A1[6]))

        
        t_MM=False
        if (t_end>=ref_dict[a]['oma1_start'] and t_start<=ref_dict[a]['oma1_end']):
            for M in ref_dict[a]['MM']:
                if M['pos'] in range(t_start,t_end+1):
                    reads_dict['oma1_MM'].append({'start': t_start,'end': t_end,'strand':t_strand,'name': t_name,'color':'red'})
                    t_MM=True
                    break
            if t_MM==False:
                reads_dict['ambi'].append({'start': t_start,'end': t_end,'strand':t_strand,'name': t_name,'color':'gray'})
        else: 
            if (t_end<=ref_dict[a]['oma1_start'] or t_start>=ref_dict[a]['oma1_end']):
                reads_dict['ambi'].append({'start': t_start,'end': t_end,'strand':t_strand,'name': t_name,'color':'gray'})
            else:
                #junction reads
                reads_dict['oma1_MM'].append({'start': t_start,'end': t_end,'strand':t_strand,'name': t_name,'color':'red'})
            

  
    infile_bowtie_map.close()
    rm_cmd='rm '+bowtie_output
    subprocess.run(rm_cmd, shell=True)

    # Calculate the number of siRNA in different categories.
    report_name='report_'+l+'_'+a+'_2.txt'
    report_out=open(report_name,'w')
    N_total_col=0 # total number of collapsed reads
    N_total_unCol=0 # total number of un-collapsed reads

    infile=open(fa_file_name,'r')
    for L1 in infile:
        if L1[0]=='>':
            N_total_col+=1
            N_total_unCol+=int(L1.strip().split('_')[-1])
    infile.close()
    report_out.write(l+'\n')
    report_out.write("Reference sequence for the alignment: %s" % a)
    report_out.write("\nTotal number of collapsed reads: %i" % N_total_col)
    report_out.write("\nTotal number of Uncollapsed reads: %i" % N_total_unCol)

    # align to whole ce6
    bowtie_ref=bowtie_folder+'indexes/ce6'
    N_ce6_col=0
    N_ce6_unCol=0
    bowtie_output=a+'_ce6_align.map'
    bowtie_prog=bowtie_folder+'bowtie ' + bowtie_ref + ' -5 '+str(barcode_size)+' -f -t -v 0 -a -p '+ str(core_num)+' '+fa_file_name+' '+bowtie_output
    p=subprocess.run(bowtie_prog, shell=True, check=True,capture_output=True)
    #print(p.stderr)
    
    infile_map=open(bowtie_output,'r')
    for L1 in infile_map:
        A1=L1.strip().split('\t')
        N_ce6_col+=1.0/(1+int(A1[6]))
        N_ce6_unCol+=int((A1[0].split('_'))[-1])/(1+int(A1[6]))
    report_out.write("\nTotal number of reads that perfectly match to ce6:")
    report_out.write("\nCollapsed reads: %i" % N_ce6_col)
    report_out.write("\nUncollapsed reads: %i" % N_ce6_unCol)
    infile_map.close()
    rm_cmd='rm '+bowtie_output
    subprocess.run(rm_cmd, shell=True)
    

    # Total number of juection or MM reads 
    (N_antisense, N_sense)=(0,0)
    for T in reads_dict['oma1_MM']:
        if ((T['end']-T['start']+1) in range(sRNA_size_min, sRNA_size_max+1)) :
            if T['strand']=='-': N_antisense+=1
            if T['strand']=='+': N_sense+=1

    report_out.write("\n Total number of oma-1 siRNAs with MM or  junction, antisense: %i" % N_antisense)
    report_out.write("\n Total number of oma-1 siRNAs with MM or  junction, sense: %i" % N_sense)


    report_out.close()
    
    # track dictionary: key: line (1, 2, etc), value: a list of tracks
    track_dict_sense = {1:[]}
    track_dict_antisense = {1:[]}
    for T_type in ['oma1_MM','ambi']:
        for T1 in reads_dict[T_type]:
            if T1['strand'] == '+' and ((T1['end']-T1['start']+1) in range(sRNA_size_min, sRNA_size_max+1)) :
                overlap_flag='no'
                for track_line in range(1,len(track_dict_sense)+1): # check the latest track line
                    overlap_flag='no'
                    for T2 in track_dict_sense[track_line]:
                        # check if T1 and T2 overlap (extend start and end by ext_nt nt)
                        if T1['start'] in range(T2['start']-ext_nt, T2['end']+ext_nt) or T1['end'] in range(T2['start']-ext_nt, T2['end']+ext_nt):
                            overlap_flag='yes'
                            break
                        if T2['start'] in range(T1['start']-ext_nt, T1['end']+ext_nt) or T2['end'] in range(T1['start']-ext_nt, T1['end']+ext_nt) :
                            overlap_flag='yes'
                            break
                    if overlap_flag == 'no': # if no overlap with any of the existing tracks in the current line
                        track_dict_sense[track_line].append({'start': T1['start'],
                                                             'end': T1['end'],
                                                             'type':T_type,
                                                             'color': T1['color']})
                        break
                if overlap_flag == 'yes': # if it overlaps with a track in all existing line, add a new line, add the track to the new line
                    track_dict_sense[len(track_dict_sense)+1]=[{'start': T1['start'],
                                                             'end': T1['end'],
                                                             'type':T_type,
                                                             'color': T1['color']}]

            if T1['strand'] == '-' and ((T1['end']-T1['start']+1) in range(sRNA_size_min, sRNA_size_max+1)):
                for track_line in range(1,len(track_dict_antisense)+1): # check the latest track line
                    overlap_flag='no'
                    for T2 in track_dict_antisense[track_line]:
                        overlap_flag='no'
                        # check if T1 and T2 overlap (extend start and end by ext_nt nt)
                        if T1['start'] in range(T2['start']-ext_nt, T2['end']+ext_nt) or T1['end'] in range(T2['start']-ext_nt, T2['end']+ext_nt):
                            overlap_flag='yes'
                            break
                        if T2['start'] in range(T1['start']-ext_nt, T1['end']+ext_nt) or T2['end'] in range(T1['start']-ext_nt, T1['end']+ext_nt) :
                            overlap_flag='yes'
                            break
                    if overlap_flag == 'no': # if no overlap with any of the existing tracks in the current line
                        track_dict_antisense[track_line].append({'start': T1['start'],
                                                             'end': T1['end'],
                                                             'type':T_type,
                                                             'color': T1['color']})
                        break
                if overlap_flag == 'yes': # if it overlaps with a track in all existing line, add a new line, add the track to the new line
                    track_dict_antisense[len(track_dict_antisense)+1]=[{'start': T1['start'],
                                                             'end': T1['end'],
                                                             'type':T_type,
                                                             'color': T1['color']}]


    # VSG plot
    VSG_py_tracks_name=l+'_'+a+'_zoom.py'
    VSG_plot_tracks_name=l+'_'+a+'_zoom.svg'
    outfile_VSG_py=open(VSG_py_tracks_name,'w')
    outfile_VSG_py.write('from VSG_Module import * \n')

    (plot_start,plot_end)= (1, 1843) 
    (oma1_start,oma1_end)=(ref_dict[a]['oma1_start'],ref_dict[a]['oma1_end'])

    for track_line in range(len(track_dict_sense),0,-1):
        #print track_line
        VSG_y=-1.0*(len(track_dict_sense)-track_line)*(track_width+line_gap)
        for T in track_dict_sense[track_line]:
            (T_start, T_end)=(T['start'],T['end'])
            if T_start in range(plot_start,plot_end+1) and T_end in range(plot_start,plot_end+1) :
                outfile_VSG_py.write('vline(x1='+str(T_start*L_per_nt)+', x2='+str(T_end*L_per_nt)+', y1='+str(VSG_y)+', y2='+str(VSG_y)+', strokewidth='+str(track_width)+' ,stroke='+T['color']+')\n')

    # plot the gDNA b/t plot_start and plot_end
    gDNA_xc=(plot_end+plot_start)/2*L_per_nt
    gDNA_yc=-1.0*(len(track_dict_sense))*(track_width+line_gap) - 0.5*gDNA_track_width
    gDNA_width=(plot_end-plot_start+1)*L_per_nt
    outfile_VSG_py.write("vrect(xc="+str(gDNA_xc)+", yc="+str(gDNA_yc)+", width="+str(gDNA_width)+", height="+ str(gDNA_track_width) +", fill=grey,stroke=grey )\n")

    # plot oma1 insert
    gDNA_xc=(oma1_start+oma1_end)/2*L_per_nt
    gDNA_yc=-1.0*(len(track_dict_sense))*(track_width+line_gap) - 0.5*gDNA_track_width
    gDNA_width=(oma1_end-oma1_start+1)*L_per_nt
    outfile_VSG_py.write("vrect(xc="+str(gDNA_xc)+", yc="+str(gDNA_yc)+", width="+str(gDNA_width)+", height="+ str(gDNA_track_width) +", fill=purple,stroke=purple )\n")

    # plot the mismatches:
    
    for M in ref_dict[a]['MM']:
        (L_x, L_y1, L_y2)=(M['pos']*L_per_nt,-1.0*(len(track_dict_sense))*(track_width+line_gap) - 1.0*gDNA_track_width,-1.0*(len(track_dict_sense))*(track_width+line_gap))
        outfile_VSG_py.write("vline(x1="+str(L_x)+", y1="+str(L_y1)+", x2="+str(L_x)+", y2="+ str(L_y2) +", strokewidth=0.5, stroke=white )\n")

    for track_line in range(1,len(track_dict_antisense)+1):
        #print track_line
        VSG_y=-1.0*(len(track_dict_sense)+track_line)*(track_width+line_gap) - 1.5*gDNA_track_width
        for T in track_dict_antisense[track_line]:
            (T_start, T_end)=(T['start'],T['end'])
            if T_start in range(plot_start,plot_end+1) and T_end in range(plot_start,plot_end+1) :
                outfile_VSG_py.write('vline(x1='+str(T_start*L_per_nt)+', x2='+str(T_end*L_per_nt)+', y1='+str(VSG_y)+', y2='+str(VSG_y)+', strokewidth='+str(track_width)+' ,stroke='+T['color']+')\n')


    outfile_VSG_py.write("vset(bg='white')\n")
    outfile_VSG_py.write("vdisplay('"+VSG_plot_tracks_name+"')")
    outfile_VSG_py.close()

    vsg_cmd='python '+VSG_py_tracks_name
    subprocess.run(vsg_cmd,shell=True)
    

###
def RNA_track_plot_oma1_smg1_fusion(l): 
# 4/21/21 make a VSG plot for the oma-1::smg-1 transgene, 
# 
# l: library name,

# assumes 4nt barcodes
    a='oma1_smg1_fusion_mRNA'
    fa_extension='_R1.fastq.gz'
    fa_folder=dir_path(l+fa_extension,fastq_folder)
    fa_file_name=fa_folder+l+fa_extension
    barcode_size=4

    # VSG options (set for figure making for illustrator):
    ext_nt = 4 # number (nt) of extension when considering if the two tracks overlap, aka, minimal horizontal gap
    L_per_nt=0.12 # number of pixel per nt
    line_gap=0.5 #
    gDNA_track_width = 1.5
    track_width=1.3 # track width for sRNA

    bowtie_ref=bowtie_folder+'indexes/'+a
    bowtie_output=l+'_'+a+'_perfect_align.map'

#    (sRNA_size_min, sRNA_size_max)=(20,24) #inclusive for both 
    ref_dict={'oma1_smg1_fusion_mRNA':{'oma1_start': 189, 'oma1_end': 680, # start and end position of MM fragment
                                  'MM': [{'pos':191, 'WT':'A', 'mut':'G'}, {'pos':221, 'WT':'T', 'mut':'C'},
                                         {'pos':251, 'WT':'C', 'mut':'G'}, {'pos':281, 'WT':'T', 'mut':'A'},
                                         {'pos':311, 'WT':'T', 'mut':'C'}, {'pos':341, 'WT':'A', 'mut':'G'},
                                         {'pos':371, 'WT':'T', 'mut':'C'}, {'pos':401, 'WT':'T', 'mut':'C'},
                                         {'pos':431, 'WT':'A', 'mut':'G'}, {'pos':461, 'WT':'G', 'mut':'A'},
                                         {'pos':491, 'WT':'T', 'mut':'C'}, {'pos':521, 'WT':'G', 'mut':'A'},
                                         {'pos':551, 'WT':'A', 'mut':'G'}, {'pos':581, 'WT':'G', 'mut':'A'},
                                         {'pos':611, 'WT':'G', 'mut':'C'},
                                         {'pos':641, 'WT':'G', 'mut':'A'}, {'pos':671, 'WT':'A', 'mut':'G'}
                                        ]}}
    
    bowtie_prog=bowtie_folder+'bowtie ' + bowtie_ref + ' -5 '+str(barcode_size)+' -q -t -v 0 -a -3 101 -p '+str(core_num)+' '+fa_file_name+' '+bowtie_output
    p=subprocess.run(bowtie_prog, shell=True, check=True,capture_output=True)
    #print(p.stderr)
    
    reads_dict={'oma1_MM':[],'ambi':[]}
    infile_bowtie_map=open(bowtie_output, 'r')
    c=0
    for L1 in infile_bowtie_map:
        c+=1
        A1=L1.strip().split('\t')
        (t_name, t_start, t_size, t_end, t_strand, t_else) =(A1[0], int(A1[3])+1, len(A1[4]), int(A1[3])+len(A1[4]), A1[1],int(A1[6]))

        
        t_MM=False
        if (t_end>=ref_dict[a]['oma1_start'] and t_start<=ref_dict[a]['oma1_end']):
            for M in ref_dict[a]['MM']:
                if M['pos'] in range(t_start,t_end+1):
                    reads_dict['oma1_MM'].append({'start': t_start,'end': t_end,'strand':t_strand,'name': t_name,'color':'red'})
                    t_MM=True
                    break
            if t_MM==False:
                reads_dict['ambi'].append({'start': t_start,'end': t_end,'strand':t_strand,'name': t_name,'color':'gray'})
        else: 
            if (t_end<=ref_dict[a]['oma1_start'] or t_start>=ref_dict[a]['oma1_end']):
                reads_dict['ambi'].append({'start': t_start,'end': t_end,'strand':t_strand,'name': t_name,'color':'gray'})
            else:
                #junction reads
                reads_dict['oma1_MM'].append({'start': t_start,'end': t_end,'strand':t_strand,'name': t_name,'color':'red'})
            

  
    infile_bowtie_map.close()
    rm_cmd='rm '+bowtie_output
    subprocess.run(rm_cmd, shell=True)

    # Calculate the number of siRNA in different categories.
    report_name='report_'+l+'_'+a+'_2.txt'
    report_out=open(report_name,'w')
    N_total_col=0 # total number of collapsed reads
    N_total_unCol=0 # total number of un-collapsed reads

    report_out.write(l+'\n')
    report_out.write("Reference sequence for the alignment: %s" % a)


    # align to whole ce6
    bowtie_ref=bowtie_folder+'indexes/ce6'
    N_ce6_col=0
    N_ce6_unCol=0
    bowtie_output=a+'_ce6_align.map'
    bowtie_prog=bowtie_folder+'bowtie ' + bowtie_ref + ' -5 '+str(barcode_size)+' -q -t -v 0 -a -3 101 -p '+ str(core_num)+' '+fa_file_name+' '+bowtie_output
    p=subprocess.run(bowtie_prog, shell=True, check=True,capture_output=True)
    #print(p.stderr)
    
    infile_map=open(bowtie_output,'r')
    for L1 in infile_map:
        A1=L1.strip().split('\t')
        N_ce6_col+=1.0/(1+int(A1[6]))
    report_out.write("\nTotal number of reads that perfectly match to ce6:")
    report_out.write("\n %i" % N_ce6_col)
    infile_map.close()
    rm_cmd='rm '+bowtie_output
    subprocess.run(rm_cmd, shell=True)
    

    # Total number of juection or MM reads 
    (N_antisense, N_sense)=(0,0)
    for T in reads_dict['oma1_MM']:
        if T['strand']=='-': N_antisense+=1
        if T['strand']=='+': N_sense+=1

    report_out.write("\n Total number of oma-1 siRNAs with MM or  junction, antisense: %i" % N_antisense)
    report_out.write("\n Total number of oma-1 siRNAs with MM or  junction, sense: %i" % N_sense)


    report_out.close()
    
    # track dictionary: key: line (1, 2, etc), value: a list of tracks
    track_dict_sense = {1:[]}
    track_dict_antisense = {1:[]}
    for T_type in ['oma1_MM','ambi']:
        for T1 in reads_dict[T_type]:
            if T1['strand'] == '+'  :
                overlap_flag='no'
                for track_line in range(1,len(track_dict_sense)+1): # check the latest track line
                    overlap_flag='no'
                    for T2 in track_dict_sense[track_line]:
                        # check if T1 and T2 overlap (extend start and end by ext_nt nt)
                        if T1['start'] in range(T2['start']-ext_nt, T2['end']+ext_nt) or T1['end'] in range(T2['start']-ext_nt, T2['end']+ext_nt):
                            overlap_flag='yes'
                            break
                        if T2['start'] in range(T1['start']-ext_nt, T1['end']+ext_nt) or T2['end'] in range(T1['start']-ext_nt, T1['end']+ext_nt) :
                            overlap_flag='yes'
                            break
                    if overlap_flag == 'no': # if no overlap with any of the existing tracks in the current line
                        track_dict_sense[track_line].append({'start': T1['start'],
                                                             'end': T1['end'],
                                                             'type':T_type,
                                                             'color': T1['color']})
                        break
                if overlap_flag == 'yes': # if it overlaps with a track in all existing line, add a new line, add the track to the new line
                    track_dict_sense[len(track_dict_sense)+1]=[{'start': T1['start'],
                                                             'end': T1['end'],
                                                             'type':T_type,
                                                             'color': T1['color']}]

            if T1['strand'] == '-' :
                for track_line in range(1,len(track_dict_antisense)+1): # check the latest track line
                    overlap_flag='no'
                    for T2 in track_dict_antisense[track_line]:
                        overlap_flag='no'
                        # check if T1 and T2 overlap (extend start and end by ext_nt nt)
                        if T1['start'] in range(T2['start']-ext_nt, T2['end']+ext_nt) or T1['end'] in range(T2['start']-ext_nt, T2['end']+ext_nt):
                            overlap_flag='yes'
                            break
                        if T2['start'] in range(T1['start']-ext_nt, T1['end']+ext_nt) or T2['end'] in range(T1['start']-ext_nt, T1['end']+ext_nt) :
                            overlap_flag='yes'
                            break
                    if overlap_flag == 'no': # if no overlap with any of the existing tracks in the current line
                        track_dict_antisense[track_line].append({'start': T1['start'],
                                                             'end': T1['end'],
                                                             'type':T_type,
                                                             'color': T1['color']})
                        break
                if overlap_flag == 'yes': # if it overlaps with a track in all existing line, add a new line, add the track to the new line
                    track_dict_antisense[len(track_dict_antisense)+1]=[{'start': T1['start'],
                                                             'end': T1['end'],
                                                             'type':T_type,
                                                             'color': T1['color']}]


    # VSG plot
    VSG_py_tracks_name=l+'_'+a+'_zoom.py'
    VSG_plot_tracks_name=l+'_'+a+'_zoom.svg'
    outfile_VSG_py=open(VSG_py_tracks_name,'w')
    outfile_VSG_py.write('from VSG_Module import * \n')

    (plot_start,plot_end)= (1, 1843) 
    (oma1_start,oma1_end)=(ref_dict[a]['oma1_start'],ref_dict[a]['oma1_end'])

    for track_line in range(len(track_dict_sense),0,-1):
        #print track_line
        VSG_y=-1.0*(len(track_dict_sense)-track_line)*(track_width+line_gap)
        for T in track_dict_sense[track_line]:
            (T_start, T_end)=(T['start'],T['end'])
            if T_start in range(plot_start,plot_end+1) and T_end in range(plot_start,plot_end+1) :
                outfile_VSG_py.write('vline(x1='+str(T_start*L_per_nt)+', x2='+str(T_end*L_per_nt)+', y1='+str(VSG_y)+', y2='+str(VSG_y)+', strokewidth='+str(track_width)+' ,stroke='+T['color']+')\n')

    # plot the gDNA b/t plot_start and plot_end
    gDNA_xc=(plot_end+plot_start)/2*L_per_nt
    gDNA_yc=-1.0*(len(track_dict_sense))*(track_width+line_gap) - 0.5*gDNA_track_width
    gDNA_width=(plot_end-plot_start+1)*L_per_nt
    outfile_VSG_py.write("vrect(xc="+str(gDNA_xc)+", yc="+str(gDNA_yc)+", width="+str(gDNA_width)+", height="+ str(gDNA_track_width) +", fill=grey,stroke=grey )\n")

    # plot oma1 insert
    gDNA_xc=(oma1_start+oma1_end)/2*L_per_nt
    gDNA_yc=-1.0*(len(track_dict_sense))*(track_width+line_gap) - 0.5*gDNA_track_width
    gDNA_width=(oma1_end-oma1_start+1)*L_per_nt
    outfile_VSG_py.write("vrect(xc="+str(gDNA_xc)+", yc="+str(gDNA_yc)+", width="+str(gDNA_width)+", height="+ str(gDNA_track_width) +", fill=purple,stroke=purple )\n")

    # plot the mismatches:
    
    for M in ref_dict[a]['MM']:
        (L_x, L_y1, L_y2)=(M['pos']*L_per_nt,-1.0*(len(track_dict_sense))*(track_width+line_gap) - 1.0*gDNA_track_width,-1.0*(len(track_dict_sense))*(track_width+line_gap))
        outfile_VSG_py.write("vline(x1="+str(L_x)+", y1="+str(L_y1)+", x2="+str(L_x)+", y2="+ str(L_y2) +", strokewidth=0.5, stroke=white )\n")

    for track_line in range(1,len(track_dict_antisense)+1):
        #print track_line
        VSG_y=-1.0*(len(track_dict_sense)+track_line)*(track_width+line_gap) - 1.5*gDNA_track_width
        for T in track_dict_antisense[track_line]:
            (T_start, T_end)=(T['start'],T['end'])
            if T_start in range(plot_start,plot_end+1) and T_end in range(plot_start,plot_end+1) :
                outfile_VSG_py.write('vline(x1='+str(T_start*L_per_nt)+', x2='+str(T_end*L_per_nt)+', y1='+str(VSG_y)+', y2='+str(VSG_y)+', strokewidth='+str(track_width)+' ,stroke='+T['color']+')\n')


    outfile_VSG_py.write("vset(bg='white')\n")
    outfile_VSG_py.write("vdisplay('"+VSG_plot_tracks_name+"')")
    outfile_VSG_py.close()

    vsg_cmd='python '+VSG_py_tracks_name
    subprocess.run(vsg_cmd,shell=True)

def bigWig_chip_pairend_ce10(lib,end_3_remove): # generate bigwig file using paired end chip-seq sequencing
# 7/31/21 
# lib, eg, "SG0621_lib1"
# end_3_remove: size of nt removed from 3'end
    (fastq1,fastq2)=('%s_R1.fastq.gz' % lib, '%s_R2.fastq.gz' % lib,)
    fq_file1=(dir_path(fastq1,fastq_folder))+fastq1
    fq_file2=(dir_path(fastq2,fastq_folder))+fastq2
    SAM_file=fq_align_bowtie2_SAM(fq_file1, fq_file2,lib,'ce10',0,end_3_remove)
    c1='%ssamtools view --threads 12 -b -S %s > %s.bam' % (samtoolsfolder,SAM_file,SAM_file)
    c2='%ssamtools sort --threads 12 %s.bam -o %s.sorted.bam' % (samtoolsfolder,SAM_file,SAM_file)
    c3='%ssamtools index -b %s.sorted.bam' % (samtoolsfolder,SAM_file)
    c4='%sbamCoverage -b %s.sorted.bam -o %s_ce10_rpm.bw --normalizeUsing CPM --binSize 10' % (deeptoolfolder,SAM_file,lib)
#    print(c1)
    p=subprocess.run(c1, shell=True, check=True,capture_output=True)
#    print(c1)
    print(p.stderr)
    p=subprocess.run(c2, shell=True, check=True,capture_output=True)
#    print(c2)
    print(p.stderr)
    p=subprocess.run(c3, shell=True, check=True,capture_output=True)
#    print(c3)
    print(p.stderr)
    p=subprocess.run(c4, shell=True, check=True,capture_output=True)
#    print(c4)
    print(p.stderr)

def bigWig_chip_singleEnd_ce10(lib,end_3_remove): # generate bigwig file using paired end chip-seq sequencing
# 7/31/21 
# lib, eg, "SG0621_lib1"
# end_3_remove: size of nt removed from 3'end
    #fastq1='%s_linker_removed_collapsed.fa' % lib
    fastq1='%s.fastq' % lib
    fq_file1=(dir_path(fastq1,fastq_folder))+fastq1
    SAM_file=fq_align_bowtie_SAM(fq_file1,'ce10',0,end_3_remove)

    c1='%ssamtools view --threads 12 -b -S %s > %s.bam' % (samtoolsfolder,SAM_file,SAM_file)
    c2='%ssamtools sort --threads 12 %s.bam -o %s.sorted.bam' % (samtoolsfolder,SAM_file,SAM_file)
    c3='%ssamtools index -b %s.sorted.bam' % (samtoolsfolder,SAM_file)
    c4='%sbamCoverage  -b %s.sorted.bam -o %s_ce10_rpkm.bw --extend 400 --normalizeUsing RPKM' % (deeptoolfolder,SAM_file,lib)
#    print(c1)
    p=subprocess.run(c1, shell=True, check=True,capture_output=True)
#    print(c1)
    print(p.stderr)
    p=subprocess.run(c2, shell=True, check=True,capture_output=True)
#    print(c2)
    print(p.stderr)
    p=subprocess.run(c3, shell=True, check=True,capture_output=True)
#    print(c3)
    print(p.stderr)
    p=subprocess.run(c4, shell=True, check=True,capture_output=True)
#    print(c4)
    print(p.stderr)
    
    
def bigWig_chip_pairend_hg19(lib,end_3_remove): # generate bigwig file using paired end chip-seq sequencing
# 7/31/21 
# lib, eg, "SG0621_lib1"
# end_3_remove: size of nt removed from 3'end
    (fastq1,fastq2)=('%s_R1.fastq.gz' % lib, '%s_R2.fastq.gz' % lib,)
    fq_file1=(dir_path(fastq1,fastq_folder))+fastq1
    fq_file2=(dir_path(fastq2,fastq_folder))+fastq2
    SAM_file=fq_align_bowtie2_SAM(fq_file1, fq_file2,lib,'hg19',0,end_3_remove)
    c1='%ssamtools view --threads 12 -b -S %s > %s.bam' % (samtoolsfolder,SAM_file,SAM_file)
    c2='%ssamtools sort --threads 12 %s.bam -o %s.sorted.bam' % (samtoolsfolder,SAM_file,SAM_file)
    c3='%ssamtools index -b %s.sorted.bam' % (samtoolsfolder,SAM_file)
    c4='%sbamCoverage -b %s.sorted.bam -o %s_ce10_rpkm.bw --normalizeUsing RPKM' % (deeptoolfolder,SAM_file,lib)
#    print(c1)
    p=subprocess.run(c1, shell=True, check=True,capture_output=True)
#    print(c1)
    print(p.stderr)
    p=subprocess.run(c2, shell=True, check=True,capture_output=True)
#    print(c2)
    print(p.stderr)
    p=subprocess.run(c3, shell=True, check=True,capture_output=True)
#    print(c3)
    print(p.stderr)
    p=subprocess.run(c4, shell=True, check=True,capture_output=True)
#    print(c4)
    print(p.stderr)
    
def sRNA_fastq_to_ucsc_bed_ce6(lib):
    # generate ucsc bed tracks for ucsc genome browser bed files, using ce6 
    # + and - strands reads are in different tracks
    # identical alignments are combined 
    # counts for each alignment are indicated
    # unique alignments are in black, nonuniq alignments are in grey
    ref='ce6'
    r5=4
    r3=0
    fq_extension='_3linker_removed.fastq.gz'
    fq_folder=dir_path(lib+fq_extension,fastq_folder)
    fp=fq_folder+lib+fq_extension
    bowtie_map=fa_fq_align_bowtie(fp,ref,r5,r3)
    FH=open(bowtie_map,'r')
    align_dict={'+':{},
                '-':{}}
    for ch in ['I','II','III','IV','V','X']:
        align_dict['+']['chr'+ch]={}
        align_dict['-']['chr'+ch]={}
    for L in FH:
        A=L.strip().split('\t')
        (strand,chrom,start,length,e)=(A[1],A[2],int(A[3]),len(A[4]),int(A[-1]))
        align_ID=A[1]+A[2]+A[3]+'_'+A[4]
        if align_ID in align_dict[strand][chrom]:
            align_dict[strand][chrom][align_ID]['count']+=1
        else:
            align_dict[strand][chrom][align_ID]={'start':start,
                                                 'end':start+length,
                                                 'count':1,
                                                 'unique':e
            }
    FH.close()
    outfile=open('%s_%s.bed' % (lib,ref),'w')
    outfile.write('browser position chrIV:8,888,958-8,890,981\n')
    outfile.write('track name="%s_plus_strand" description="%s_plus_strand" visibility=pack itemRgb="On"\n' % (lib,lib))
    for ch in ['I','II','III','IV','V','X']:
        (strand,chrom)=('+','chr'+ch,)
        for align_ID in align_dict[strand][chrom]:
            if align_dict[strand][chrom][align_ID]['unique']==0:
                rgb='0,0,0'
            else:
                rgb='255,0,0'
            outfile.write('%s\t%i\t%i\t%i\t%i\t%s\t%i\t%i\t%s\n' % (chrom,align_dict[strand][chrom][align_ID]['start'],align_dict[strand][chrom][align_ID]['end'],
                                                 align_dict[strand][chrom][align_ID]['count'],0,strand,
                                                 align_dict[strand][chrom][align_ID]['start'],align_dict[strand][chrom][align_ID]['end'],rgb))
    outfile.write('track name="%s_minus_strand" description="%s_minus_strand" visibility=pack itemRgb="On"\n' % (lib,lib))
    for ch in ['I','II','III','IV','V','X']:
        (strand,chrom)=('-','chr'+ch,)
        for align_ID in align_dict[strand][chrom]:
            if align_dict[strand][chrom][align_ID]['unique']==0:
                rgb='0,0,0'
            else:
                rgb='255,0,0'
            outfile.write('%s\t%i\t%i\t%i\t%i\t%s\t%i\t%i\t%s\n' % (chrom,align_dict[strand][chrom][align_ID]['start'],align_dict[strand][chrom][align_ID]['end'],
                                                 align_dict[strand][chrom][align_ID]['count'],0,strand,
                                                 align_dict[strand][chrom][align_ID]['start'],align_dict[strand][chrom][align_ID]['end'],rgb))
    outfile.close()
    scriptrun='gzip %s_%s.bed' % (lib,ref)
    p=subprocess.run(scriptrun, shell=True, check=True,capture_output=True)
    
def primer3_v1(seq,size_min,size_max):
    # use primer3 to design PCR primers 8/14/21
    # seq: input sequence, size_min, size_max: minimal and maximal sizes of the PCR product
    # return a dictionary: {"PCR_seq": sequence of PCR product, "left_primer": ####, "right_primer": ####}
    # this is a simple function only return the first pair of the PCR primers reported by primer3
    # sys.exit()  if there is no PCR primers can be found or Primer3 reports an error
    input_file=open('primer3_input_temp','w')
    
    input_file.write('SEQUENCE_ID=example\n')
    input_file.write('SEQUENCE_TEMPLATE=%s\n' % seq)
    input_file.write('PRIMER_TASK=generic\n')
    input_file.write('PRIMER_PICK_LEFT_PRIMER=1\n')
    input_file.write('PRIMER_PICK_INTERNAL_OLIGO=0\n')
    input_file.write('PRIMER_PICK_RIGHT_PRIMER=1\n')
    input_file.write('PRIMER_OPT_SIZE=20\n')
    input_file.write('PRIMER_MIN_SIZE=18\n')
    input_file.write('PRIMER_MAX_SIZE=22\n')
    input_file.write('PRIMER_PRODUCT_SIZE_RANGE=%i-%i\n' % (size_min,size_max))
    input_file.write('PRIMER_EXPLAIN_FLAG=1\n')
    input_file.write('=\n')
    input_file.close()
    primer3_cmd='/home/sam/Dropbox/bioinformatics/primer3/primer3/src/primer3_core primer3_input_temp > primer3_output_temp'
    p=subprocess.run(primer3_cmd, shell=True, check=True,capture_output=True)
    #print(p.stderr)
    out_file=open('primer3_output_temp','r')
    for L1 in out_file:
        A1=L1.strip().split('=')
        if L1.strip()=='PRIMER_PAIR_NUM_RETURNED=0' or A1[0]=='PRIMER_ERROR':
            sys.exit(L1)
        if A1[0]=='PRIMER_LEFT_0_SEQUENCE':
            p_l=A1[1]
        if A1[0]=='PRIMER_RIGHT_0_SEQUENCE':
            p_r=A1[1]
        if A1[0]+'='=='PRIMER_LEFT_0=':
            p_s=int(A1[1].split(',')[0]) # Primer3 uses 0 for first base
        if A1[0]+'='=='PRIMER_RIGHT_0=':
            p_e=int(A1[1].split(',')[0])+1 # Primer3 uses 0 for first base
    out_file.close()
    return ({'PCR_seq':seq[p_s:p_e],'left_primer':p_l,'right_primer':p_r})            

def get_cDNA_fasta(g):
    # get ce6 cDNA sequence using sanger gene name (not the three letter ones)
    fastafile=open(HOME+'/Dropbox/bioinformatics/Genome_ce6/ce6_cDNA_v1.fa','r')
    c=0
    for L1 in fastafile:
        c+=1
        if c%2==1:
            gn=L1.strip()
        else:
            if gn=='>'+g:
                print("find %s" % g)
                return L1.strip()
    fastafile.close()
    
def bigWig_RNAseq_singleEnd_ce10_v1(lib,end_3_remove): 
# generate bigwig file using single-ended RNA sequencing
# can take both fastq and fa files
# for old 3' linker ligation based, file name: SG##_lib##_3linker_removed.fastq.gz
#  08/21/21
# lib, eg, "SG0621_lib1"
# end_3_remove: size of nt removed from 3'end, use 0 for 3'linker removed fastq
# use 10 nt for bamCoverage binSize
#    (fastq1)=('%s_3linker_removed.fastq.gz' % lib)
    fastq1='%s_linker_removed_collapsed.fa' % lib
    fq_file1=(dir_path(fastq1,fastq_folder))+fastq1
    SAM_file=fq_align_bowtie_SAM(fq_file1,'ce10',4,end_3_remove)
    c1='%ssamtools view --threads 12 -b -S %s > %s.bam' % (samtoolsfolder,SAM_file,SAM_file)
    c2='%ssamtools sort --threads 12 %s.bam -o %s.sorted.bam' % (samtoolsfolder,SAM_file,SAM_file)
    c3='%ssamtools index -b %s.sorted.bam' % (samtoolsfolder,SAM_file)
    c4='%sbamCoverage -b %s.sorted.bam -o %s_ce6_rpkm_rev.bw --binSize 10 --normalizeUsing RPKM --filterRNAstrand forward' % (deeptoolfolder,SAM_file,lib) # note that bamcoverage assumes dUTP-based RNA-seq, so opposite
    c5='%sbamCoverage -b %s.sorted.bam -o %s_ce6_rpkm_fr.bw --binSize 10 --normalizeUsing RPKM --filterRNAstrand reverse' % (deeptoolfolder,SAM_file,lib)
#    print(c1)
    p=subprocess.run(c1, shell=True, check=True,capture_output=True)
#    print(c1)
    print(p.stderr)
    p=subprocess.run(c2, shell=True, check=True,capture_output=True)
#    print(c2)
    print(p.stderr)
    p=subprocess.run(c3, shell=True, check=True,capture_output=True)
#    print(c3)
    print(p.stderr)
    p=subprocess.run(c4, shell=True, check=True,capture_output=True)
#    print(c4)
    print(p.stderr)
    p=subprocess.run(c5, shell=True, check=True,capture_output=True)
#    print(c5)
    print(p.stderr)