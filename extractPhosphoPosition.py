#!/usr/bin/env python
# coding: utf-8

# In[117]:


# packages
import os
import sys #not sure if needed with click
import click
import getopt #not sure if needed
import shutil
import csv
import re
import pandas as pd
from Bio import SeqIO


# In[118]:


# classes and functions
class SeqContainer:
    general_mod = []
    rt = float
    accession = str
    accession_first = str
    ph_init_seq = str
    ph_trans_seq = str
    db_seq = str
    ph_pep_aa = []
    ph_pep_pos = []
    ph_prot_pos = []

# fills inital SeqContainer datastructure
# returns list of SeqContainer [0] and list of modification [1]
def buildSeqContainer(db_path, analysis_path):
    # read input 
    df = pd.read_csv(analysis_path)
    seq_col = df.sequence
    accession_col = df.accessions
    rt_col = df.rt_cf
    
    if(len(seq_col) != len(accession_col)):
        print("Error: sequence and accessions columns are not equal in lenght.")
        sys.exit()

    # find all modifictions apart from Phosphorylation ()
    l_mod = []
    for entry in seq_col:
        for m in re.findall('\(.+?\)', entry):
            l_mod.append(m)
    l_mod = list(set(l_mod))

    # remove 'Phospho' modification will be used later on.
    l_mod.remove('(Phospho)')

    # extract allways first accssion if multiple in string
    l_ac_first =[]
    for element in accession_col:
        ac_first = element.split(';',1)[0]
        l_ac_first.append(ac_first)    
    
    # list accession phosphosites
    d_ac_ph = {}
    l_ac_ph_rt = []
    count_seq = 0
    for entry in seq_col:
        if (entry != "UNIDENTIFIED_PEPTIDE"):
            l_ac_ph_rt.append([l_ac_first[count_seq],accession_col[count_seq],seq_col[count_seq],rt_col[count_seq]])
        count_seq += 1
            
    # dictionary accession sequence (fasta)
    d_id_seq = {}
    for seq_record in SeqIO.parse(db_path, "fasta"):
        d_id_seq[seq_record.id] = seq_record.seq
        
    # combine the list and dict in class SeqContainer 
    l_seqc = [] 
    l_seqc_mod = []
    for element in l_ac_ph_rt:
        for k_db, v_db in d_id_seq.iteritems():
            if (element[0] == k_db):
                seqc = SeqContainer()
                seqc.general_mod = l_mod
                accession_first = element[0]
                seqc.accession = element[1]
                seqc.ph_init_seq = element[2]
                seqc.rt = element[3]
                seqc.db_seq = v_db
                l_seqc.append(seqc)  
        
    return l_seqc
    
# extract phosphosite position and aminoacid in peptide and protein
# returns updated SeqContainer
def extractPhosphoPosition(SeqContainer):
    ph_str = SeqContainer.ph_init_seq
    if (ph_str == ""):
        print("Warning: Empty peptide squence was provided")
    ph_str = ph_str.replace('.','')
    # remove all modification apart from 'Phospho'
    for element in SeqContainer.general_mod:
        if element in ph_str:
            ph_str = ph_str.replace(element,'')
    ph_str_p = ph_str
    ph_str_a = ph_str_p.replace('(Phospho)','')
    SeqContainer.ph_trans_seq = ph_str_a

    l_site = []
    for m in re.finditer('(Phospho)', ph_str_p):
        l_site.append(m.start())

    # aminoacid = elem-2
    # position = elem-1
    l_pep_aa = []
    l_pep_pos = []
    for elem in l_site:
        l_pep_aa.append(ph_str_p[elem-2])
        l_pep_pos.append(elem-1)

    # correct for the lenght of "(Phospho)"
    for i in range(len(l_pep_pos)):
        l_pep_pos[i] = l_pep_pos[i] - (i*9)

    SeqContainer.ph_pep_aa = l_pep_aa
    SeqContainer.ph_pep_pos = l_pep_pos

    # if position was located find position in protein
    if (len(SeqContainer.ph_pep_aa) != 0 and len(SeqContainer.ph_pep_pos) != 0):
        peptide_position = []
        if (len(SeqContainer.db_seq) != 0):
        # look for position of peptide in db_sequence 
            for m in re.finditer(str(SeqContainer.ph_trans_seq), str(SeqContainer.db_seq)):
                peptide_position.append(m.start())
        else: 
            print("Database sequence was not found")

        if (len(peptide_position) != 0):
            l_prot_pos = []
            # add the posotion to the value in the peptide
            for element in SeqContainer.ph_pep_pos:
                l_prot_pos.append(element + peptide_position[0])
            SeqContainer.ph_prot_pos = l_prot_pos
        else:
            print("Was not able to localise peptide in database seqence")
    
    return SeqContainer

# build new dataframe
def buildDataFrameFromContainer(l_container):
    # get number of columns which are need in the csv
    max_len_position = 0
    for entry in l_container:
        len_position = len(entry.ph_pep_pos)
        if len_position > max_len_position:
            max_len_position = len_position

    # columns pep_position, prot_position
    # saved as AAPos S5 Y12 - S143 Y 150 
    # max_len_position * 2 + 2 (asccession, pep_seq)
    header_pep = []
    header_prot = []
    for i in range(max_len_position):
        header_pep.append('position_pep_'+str(i))
        header_prot.append('position_prot_'+str(i))

    header = [['rt_cf','accessions','sequence']]
    header.append(header_pep)
    header.append(header_prot)

    # flatten header list 
    header = [y for x in header for y in x]

    # make new dataframe
    df_container = pd.DataFrame(columns = header)

    # df[row,col]
    i=0
    for element in l_container:
        df_container.at[i,'rt_cf'] = element.rt    
        df_container.at[i,'accessions'] = element.accession
        df_container.at[i,'sequence'] = element.ph_init_seq
        if (len(element.ph_pep_aa) > 0):
            for j in range(len(element.ph_pep_aa)):
                df_container.at[i,'position_pep_'+str(j)] = str(str(element.ph_pep_aa[j])+str(element.ph_pep_pos[j]))
                df_container.at[i,'position_prot_'+str(j)] = str(str(element.ph_pep_aa[j])+str(element.ph_prot_pos[j]))
        i+=1
    
    # convert column rt_cf to numeric (float64)
    df_container['rt_cf'] = df_container.rt_cf.astype(float)

    return df_container


# In[119]:


def main():
    
    database = "/Volumes/elements/Ph_analysis/database/database/2018/swissprot_human_crap_decoy_20181017.fasta" 
    analysis = "/Volumes/elements/Ph_analysis/results/5min_HILIC/ConsensusMapNormalizer/ConsensusMapNormalizer.csv"
    output = "/Volumes/elements/Ph_analysis/results/5min_HILIC/ConsensusMapNormalizer/ConsensusMapNormalizer_output.csv"
    output_new = "/Volumes/elements/Ph_analysis/results/5min_HILIC/ConsensusMapNormalizer/ConsensusMapNormalizer_new_output.csv"
    
    # how to make command line tool with input 
    '''
    analysis = ''
    database = ''
    output = ''
    
    try:
      opts, args = getopt.getopt(argv,"hi:o:",["analysis=","database=","output="])
    except getopt.GetoptError:
      print 'extractPhosphoPosition.py -a <analysis> -d <database> -o <output>'
      sys.exit(2)
    for opt, arg in opts:
      if opt == '-h':
         print 'extractPhosphoPosition.py -a <analysiså> -d <database> -o <output>'
         sys.exit()
      elif opt in ("-a", "--analysis"):
         analysis = arg
      elif opt in ("-d", "--database"):
         database = arg
      elif opt in ("-o", "--output"):
         output = arg
    print('Input file is "', analysis)
    print('Database file is "', database)
    print('Output file is "', output)
    '''
    
    print("Starting the extraction")

    print("Building the datastructure")
    # build the inital SeqContainer datastructure
    l_container = buildSeqContainer(database, analysis)

    print("Extract the phosphosites")
    # extract phosphosite position and aminoacid in peptide and protein
    for entry in l_container:
        extractPhosphoPosition(entry)

    print("Construct dataframe")
    # build new dataframe from data in the container
    df_new = buildDataFrameFromContainer(l_container)
    df_new.sort_values('rt_cf')
    
    # analysis_path dataframe
    df_in = pd.read_csv(analysis_path)
    df_in.sort_values('rt_cf')
    
    print("Merging in progress")
    # merge the two dataframes by accessions and sequence -> generate final output
    df_out = pd.merge(df_in, df_new,  how='left', left_on=['rt_cf','accessions','sequence'], right_on = ['rt_cf','accessions','sequence'])
    df_out.sort_values('rt_cf')
    
    # write output
    df_new.to_csv(output_new, sep=',', encoding='utf-8')
    df_out.to_csv(output, sep=',', encoding='utf-8')
    
    print("Done")
    
    return 0


# In[120]:


if __name__ == "__main__":
    main()
    #main(sys:argv[1:])

