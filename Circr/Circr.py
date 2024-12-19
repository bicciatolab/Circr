import pandas as pd
pd.options.mode.chained_assignment = None
import numpy as np
import sys,re,os,time,argparse,collections,csv,pybedtools
import logging
from re import compile
from re import search
from sys import argv
from collections import defaultdict
from functools import partial
from operator import itemgetter
from multiprocessing import Pool
from itertools import islice
from collections import Counter

#@#@#@# DICTIONARY SECTION #@#@#@#

match = {'CG': '|',
		 'GC': '|',
		 'AU': '|',
		 'UA': '|',
		 'UG': ':',
		 'GU': ':'}

seeds_a = {' |||||| ': '7A1',
	  	   ':|||||| ': '7mer-A1',
		   '||||||| ': '8mer',
		   '||| ||| ': 'SM',
		   '|||:||| ': 'SM',
		   '||||||  ': 'Off-6mer'}

seeds_b = {' |||||| ': '6mer',
		   ':|||||| ': '6mer',
		   '||||||| ': '7mer-m8'}

DataTree = {'human': {
				'file': 'support_files/miRNA/hsa_mature.fa',
				'acronym': 'hsa',
				'hg19':{
						'AGO': 'support_files/human/hg19/hg19.AGO.bed',
						'INT': 'support_files/human/hg19/hg19.INT.bed',
						'GTF': 'support_files/human/hg19/hg19.ensGene.gtf',
		 				'TAXA': 9606,
		 				'Genome': 'support_files/human/hg19/hg19.fa',
		 				'rRNA_Coords': 'support_files/human/hg19/hg19.rRNA.bed'
						},
		 		'hg38':{
						'AGO': 'support_files/human/hg38/hg38.AGO.bed',
						'INT': 'support_files/human/hg38/hg38.INT.bed',
 						'GTF': 'support_files/human/hg38/hg38.ensGene.gtf',
 		 				'TAXA': 9606,
 		 				'Genome': 'support_files/human/hg38/hg38.fa',
 		 				'rRNA_Coords': 'support_files/human/hg38/hg38.rRNA.bed'
		 				}
		  		},
			'mouse': {
				'file' : 'support_files/miRNA/mmu_mature.fa',
				'acronym' : 'mmu',
				'mm9':{
						'AGO': 'support_files/mouse/mm9/mm9.AGO.bed',
						'INT': 'support_files/mouse/mm9/mm9.INT.bed',
						'GTF': 'support_files/mouse/mm9/mm9.ensGene.gtf',
		 				'TAXA': 10090,
		 				'Genome': 'support_files/mouse/mm9/mm9.fa',
		 				'rRNA_Coords': 'support_files/mouse/mm9/mm9.rRNA.bed'
					  },
		 		'mm10':{
						'AGO': 'support_files/mouse/mm10/mm10.AGO.bed',
						'INT': 'support_files/mouse/mm10/mm10.INT.bed',
 						'GTF': 'support_files/mouse/mm10/mm10.ensGene.gtf',
 		 				'TAXA': 10090,
 		 				'Genome': 'support_files/mouse/mm10/mm10.fa',
 		 				'rRNA_Coords': 'support_files/mouse/mm10/mm10.rRNA.bed'
		 			   }
			     },
			'worm': {
				'file': 'support_files/miRNA/cel_mature.fa',
				'acronym': 'cel',
				'ce10':{
						'AGO': 'support_files/worm/ce10/ce10.AGO.bed',
						'INT': 'support_files/worm/ce10/ce10.INT.bed',
						'GTF': 'support_files/worm/ce10/ce10.ensGene.gtf',
		 				'TAXA': 6239,
		 				'Genome': 'support_files/worm/ce10/ce10.fa',
		 				'rRNA_Coords': 'support_files/worm/ce10/ce10.rRNA.bed'
					   },
			    'ce11':{
						'AGO': 'support_files/worm/ce11/ce11.AGO.bed',
						'INT': 'support_files/worm/ce11/ce11.INT.bed',
 						'GTF': 'support_files/worm/ce11/ce11.ensGene.gtf',
 		 				'TAXA': 6239,
 		 				'Genome': 'support_files/worm/ce11/ce11.fa',
 		 				'rRNA_Coords': 'support_files/worm/ce11/ce11.rRNA.bed'
		 			   }
				},
			'fruitfly': {
					'file': 'support_files/miRNA/dme_mature.fa',
					'acronym': 'dme',
					'dm3':{
							'AGO': 'support_files/fruitfly/dm3/dm3.AGO.bed',
							'INT': 'support_files/fruitfly/dcm3/dm3.INT.bed',
							'GTF': 'support_files/fruitfly/dm3/dm3.ensGene.gtf',
		 					'TAXA': 7227,
		 					'Genome': 'support_files/fruitfly/dm3/dm3.fa',
		 					'rRNA_Coords': 'support_files/fruitfly/dm3/dm3.rRNA.bed'
						  },
		    		'dm6':{
							'AGO': 'support_files/fruitfly/dm6/dm6.AGO.bed',
							'INT': 'support_files/fruitfly/dm6/dm6.INT.bed',
 							'GTF': 'support_files/fruitfly/dm6/dm6.ensGene.gtf',
 		 					'TAXA': 7227,
 		 					'Genome': 'support_files/fruitfly/dm6/dm6.fa',
 		 					'rRNA_Coords': 'support_files/fruitfly/dm6/dm6.rRNA.bed'
		 				  }
		   	  		}
		}

organisms = {'human':['hg19', 'hg38'],
			 'fruitfly':['dm3', 'dm6'],
			 'worm':['ce10', 'ce11'],
			 'mouse': ['mm9', 'mm10']}

#@#@#@# DEFINITION SECTION #@#@#@#

def getCircrPath():
	c_path = os.path.realpath(__file__).replace("Circr.py","")
	return c_path

def hybrid_parser(lines, DATA, key):
	df = pd.DataFrame(columns=range(len(lines[9].strip()[10:-3])))
	df.loc[len(df)] = [x for x in lines[9].strip()[10:-3]]
	df.loc[len(df)] = [x for x in lines[10].strip('\n')[10:-3]]
	df.loc[len(df)] = [x for x in lines[11].strip('\n')[10:-3]]
	df.loc[len(df)] = [x for x in lines[12].strip()[10:-3]]
	m = []
	target = df.head(2).max().to_list()
	mirna = df.tail(2).max().to_list()
	tuples = tuple(zip(target, mirna))
	for i in tuples:
		if ''.join(i) in match.keys():
			m.append(match[''.join(i)])
		else:
			m.append(' ')
	ll = m.index('|')
	rl = -(m[::-1].index('|')+1)
	try:
		if (-(m[::-1].index(':'))+1) > rl:
			rl = -(m[::-1].index(':')+1)
	except ValueError:
		pass
	upper = ''.join([''.join(target[:ll]).lower(), ''.join(target[ll:rl]).replace(' ','-'), ''.join(target[rl:]).lower()])
	lower = ''.join([''.join(mirna[:ll]).lower(), ''.join(mirna[ll:rl]).replace(' ','-'), ''.join(mirna[rl:]).lower()])
	if (upper.startswith('a')) and (''.join(m[-8:]) in seeds_a):
		DATA[key].append(seeds_a[''.join(m[-8:])])
	elif (not upper.startswith('a')) and (''.join(m[-8:]) in seeds_b):
		DATA[key].append(seeds_a[''.join(m[-8:])])
	else:
		DATA[key].append('No Seed')

def run_process(obj):
    os.system('{}'.format(obj))

def examine(fileIn, DB, acronym):
    DATA = {}
    key = ''
    with open('mir_TMP_'+fileIn,'r') as Miranda:
        for line in Miranda:
            if "Performing" in line.strip().split(' '):
                key = ''.join(line.strip().split(' ')[2:])
                if key not in DATA.keys():
                    DATA[key] = []
            elif "Query:" in line.strip().split(' '):
                DATA[key].append(' '.join(line.strip().split(' ')[4:]))
            elif "|" in line.strip():
                DATA[key].append(line.strip('\n')[16:])
            elif "Ref:" in line.strip().split(' '):
                DATA[key].append(' '.join(line.strip().split(' ')[6:]))
            elif line.strip().split(' ')[0].startswith('>{0}'.format(acronym)):
                DATA[key].append(line.strip().split('\t')[2])
                DATA[key].append(line.strip().split('\t')[3])
                DATA[key].extend(line.strip().split('\t')[5].split(' '))
    DATA = {key: value for key, value in DATA.items() if value}
    for key in DATA.keys():
        if(DATA[key][2][-4] == 'a') and (DATA[key][1][-8:] in seeds_a):
            DATA[key].append(seeds_a[DATA[key][1][-8:]])
        elif(DATA[key][2][-4] != 'a') and (DATA[key][1][-8:] in seeds_b):
            DATA[key].append(seeds_b[DATA[key][1][-8:]])
        else:
            DATA[key].append('No Seed')
        DB.loc[len(DB)] = [key.split('vs')[0], key.split('vs')[1], DATA[key][-1], DATA[key][-3], DATA[key][-2], 'INT_M_'+str(len(DB))]
    os.system('rm -rf {0}'.format('mir_TMP_'+fileIn))
    return DB

def hybrid(fileIn):
	DF= pd.DataFrame(columns = ['miRNA Name','Circ Name', 'Seed Category', 'Start', 'End'])
	DATA = {}
	key = ''
	counter = 15
	with open('hyb_TMP_'+fileIn,'r') as file:
		while True:
			try:
				lines = list(islice(file, counter))
				if lines[0].startswith('target too long'):
					del lines[0]
					counter = 16
					key = lines[0].strip().split(': ')[1] + 'vs' + lines[2].strip().split(': ')[1]
					if key not in DATA.keys():
						DATA[key] = [int(lines[8].strip().split(' ')[-1]),int(lines[3].strip().split(': ')[1])+int(lines[8].strip().split(' ')[-1])]
						hybrid_parser(lines, DATA, key)
						DF.loc[len(DF)] = [key.split('vs')[1], key.split('vs')[0], DATA[key][-1], DATA[key][0], DATA[key][1]]
				elif lines[0].startswith('\n'):
					del lines[0]
					counter = 15
					key = lines[0].strip().split(': ')[1] + 'vs' + lines[2].strip().split(': ')[1]
					if key not in DATA.keys():
						DATA[key] = [int(lines[8].strip().split(' ')[-1]),int(lines[3].strip().split(': ')[1])+int(lines[8].strip().split(' ')[-1])]
						hybrid_parser(lines, DATA, key)
						DF.loc[len(DF)] = [key.split('vs')[1], key.split('vs')[0], DATA[key][-1], DATA[key][0], DATA[key][1]]
				elif lines[0].startswith('target:'):
					counter = 15
					key = lines[0].strip().split(': ')[1] + 'vs' + lines[2].strip().split(': ')[1]
					if key not in DATA.keys():
						DATA[key] = [int(lines[8].strip().split(' ')[-1]),int(lines[3].strip().split(': ')[1])+int(lines[8].strip().split(' ')[-1])]
						hybrid_parser(lines, DATA, key)
						DF.loc[len(DF)] = [key.split('vs')[1], key.split('vs')[0], DATA[key][-1], DATA[key][0], DATA[key][1]]
			except (IndexError, ValueError):
				break
	DF.to_csv('hyb_Parsed_'+fileIn, index=False)
	os.system('rm -rf {0}'.format('hyb_TMP_'+fileIn))

def TS_parser(fileIn, DB, acronym):
	df = pd.read_csv('TS_TMP_'+fileIn, sep='\t')
	df = df[['a_Gene_ID','miRNA_family_ID','UTR_start','UTR_end','Site_type']]
	df = df.rename(columns={'miRNA_family_ID':'miRNA Name','a_Gene_ID':'Circ Name','Site_type':'Seed Category','UTR_start':'Start','UTR_end':'End'})
	df['ID'] = np.arange(len(DB),(len(DB)+len(df)))
	df['ID'] = 'INT_TS_'+df['ID'].astype(str)
	df['miRNA Name'] = acronym + '-' + df['miRNA Name']
	os.system('rm -rf {0}'.format('TS_TMP_'+fileIn))
	return(pd.concat([DB,df]).reset_index(drop=True))

def RH_parser(fileIn, DB):
	df = pd.read_csv('hyb_Parsed_'+fileIn)
	df['ID'] = np.arange(len(DB),(len(DB)+len(df)))
	df['ID'] = 'INT_RH_'+df['ID'].astype(str)
	os.system('rm -rf {0}'.format('hyb_Parsed_'+fileIn))
	return(pd.concat([DB,df]).reset_index(drop=True))

def single_coord(strand, coord_final, DB, circ):
	chunk = coord_final.loc[coord_final['name'] == circ]
	circ_chunk = DB.loc[DB['Circ Name'] == circ]
	if strand == '-':
		circ_chunk['fac_Start'] = chunk.iloc[0]['end'] - circ_chunk['End'].astype(int)
		circ_chunk['fac_End'] = chunk.iloc[0]['end'] - circ_chunk['Start'].astype(int)
	else:
		circ_chunk['fac_Start'] = circ_chunk['Start'].astype(int) + chunk.iloc[0]['start']
		circ_chunk['fac_End'] = circ_chunk['End'].astype(int) + chunk.iloc[0]['start']
	circ_chunk['Start'] = circ_chunk['fac_Start']
	circ_chunk['End'] = circ_chunk['fac_End']
	circ_chunk = circ_chunk.drop(['fac_Start','fac_End'], axis = 1)
	circ_chunk['Chrom'] = chunk.iloc[0]['chrom']
	circ_chunk['Strand'] = chunk.iloc[0]['strand']
	return circ_chunk

def multiple_coord(strand, coord_final, DB, circ):
	chunk = coord_final.loc[coord_final['name'] == circ]
	chunk['score'] = pd.to_numeric(chunk['score'])
	if strand == '-':
		chunk = chunk.reindex(index=chunk.index[::-1])
		circ_chunk = DB.loc[DB['Circ Name'] == circ].reset_index(drop=True)
		IDX = circ_chunk.index[circ_chunk['Start'].astype(float) < chunk['score'].iloc[0]].to_list()
		df_c = circ_chunk.iloc[IDX, :]
		circ_chunk = circ_chunk.drop(circ_chunk.index[IDX]).reset_index(drop=True)
		df_c['fac_Start'] = chunk['end'].iloc[0] - df_c['End'].astype(int)
		df_c['fac_End'] = chunk['end'].iloc[0] - df_c['Start'].astype(int)
	else:
		circ_chunk = DB.loc[DB['Circ Name'] == circ].reset_index(drop=True)
		IDX = circ_chunk.index[circ_chunk['Start'].astype(float) < chunk['score'].iloc[0]].to_list()
		df_c = circ_chunk.iloc[IDX, :]
		circ_chunk = circ_chunk.drop(circ_chunk.index[IDX]).reset_index(drop=True)
		df_c['fac_Start'] = chunk['start'].iloc[0] + df_c['Start'].astype(int)
		df_c['fac_End'] = chunk['start'].iloc[0] + df_c['End'].astype(int)
	df_c['End'] = df_c['fac_End']
	df_c['Start'] = df_c['fac_Start']
	df_c = df_c.drop(['fac_Start', 'fac_End'], axis = 1)
	df_c['Chrom'] = chunk.iloc[0]['chrom']
	df_c['Strand'] = chunk.iloc[0]['strand']
	for i in range(1, len(chunk)):
		chunk['score'].iloc[i] = chunk['score'].iloc[i-1:i+1].sum()
		IDX = circ_chunk.index[circ_chunk['Start'].astype(float) < chunk['score'].iloc[i]].to_list()
		df = circ_chunk.iloc[IDX, :]
		circ_chunk = circ_chunk.drop(circ_chunk.index[IDX]).reset_index(drop=True)
		if strand == '-':
			df['fac_Start'] = chunk['end'].iloc[i] - (df['End'].astype(int) - chunk['score'].iloc[i-1])
			df['fac_End'] = chunk['end'].iloc[i] - (df['Start'].astype(int) - chunk['score'].iloc[i-1])
		else:
			df['fac_Start'] = chunk['start'].iloc[i] + (df['Start'].astype(int) - chunk['score'].iloc[i-1])
			df['fac_End'] = chunk['start'].iloc[i] + (df['End'].astype(int) - chunk['score'].iloc[i-1])
		df['Start'] = df['fac_Start']
		df['End'] = df['fac_End']
		df = df.drop(['fac_Start','fac_End'], axis = 1)
		df['Chrom'] = chunk.iloc[i]['chrom']
		df['Strand'] = chunk.iloc[i]['strand']
		df_c = pd.concat([df_c, df])
	return df_c

def add_interactions(db, interactions, column):
	if interactions.empty == False:
		for item in list(interactions[7].unique()):
			db.loc[db.ID == item, column] = 'Yes'
	return db

def calculate_exons(bedfile, gtf, rrna):
	logging.info("Loading input bed file: %s", str(bedfile))
	input_bed = pybedtools.BedTool(bedfile)
	bed_to_df = input_bed.to_dataframe()
	logging.info("Loading input GTF file: %s",str(gtf))
	gtf_file = pybedtools.BedTool(gtf)
	logging.info("Loading rRNA file: %s", str(rrna))
	rrna_file = pybedtools.BedTool(rrna)
	data = gtf_file.intersect(input_bed, wo=True, s=True)
	logging.info("Generating exons dataframe and polishing data")
	df_data = data.to_dataframe(names=['chr','gene_biotype','feature','start','end','6','strand','7','gene_info','c_chr','c_start','c_end','c_name','c_length','c_strand','overlap']).drop(['6','7','gene_biotype'], axis=1)
	exons = df_data.loc[df_data['feature'] == 'exon']
	exons = pd.concat([exons, exons['gene_info'].str.split(';', expand=True).rename(columns={0:'gene_id',1:'transcript_id',2:'exon_number',3:'gene_name',4:'gene_biotype',5:'transcript_name'})], axis=1).drop(['gene_info'], axis=1)
	exons['gene_id'] = exons['gene_id'].map(lambda x: x.lstrip(' gene_id ')).str.replace(r'"', '')
	exons['transcript_id'] = exons['transcript_id'].map(lambda x: x.lstrip(' transcript_id ')).str.replace(r'"', '')
	exons['exon_number'] = exons['exon_number'].map(lambda x: x.lstrip(' exon_number ')).str.replace(r'"', '')
	exons['gene_name'] = exons['gene_name'].map(lambda x: x.lstrip(' gene_name ')).str.replace(r'"', '')
	exons['gene_biotype'] = exons['gene_biotype'].map(lambda x: x.lstrip(' gene_biotype ')).str.replace(r'"', '')
	exons['transcript_name'] = exons['transcript_name'].map(lambda x: x.lstrip(' transcript_name ')).str.replace(r'"', '')
	exons['start'] = exons['start']-1
	logging.info("Returning dataframe")
	return exons, bed_to_df

def create_circular_DB(exons, bed_dataframe, rrna_file, genome):
	coordinates = pd.DataFrame()
	circulars = list(exons['c_name'].unique())
	logging.info("Calculating circulars coordinates")
	for c in circulars:
		trans = list(set(exons.loc[exons['c_name'] == c, 'transcript_id'].to_list()))
		for t in trans:
			slice = exons.loc[(exons['c_name'] == c)&(exons['transcript_id'] == t)].sort_values(by=['start'])
			if(slice['start'].iloc[0] == slice['c_start'].iloc[0]) and (slice['end'].iloc[-1] == slice['c_end'].iloc[-1]):
				coordinates = coordinates.append(slice[['chr','start','end','c_name','overlap','strand']])
	coordinates = coordinates.append(bed_dataframe.loc[bed_dataframe['name'].isin(list(set(list(bed_dataframe['name'].unique())).symmetric_difference(set(list(coordinates['c_name'].unique())))))].rename(columns={'chrom':'chr','name':'c_name','score':'overlap'})).drop_duplicates()
	coord_bed = pybedtools.BedTool.from_dataframe(coordinates)
	coordinates = coord_bed.intersect(rrna_file, wa=True, v=True)
	coord_final = coordinates.to_dataframe()
	coordinates.sequence(fi=genome,fo='output_fasta.fa', tab=True, name=True, s=True)
	fasted = pd.read_csv('output_fasta.fa', sep='\t', names=['circs','sequence'])
	fasted['circs'] = fasted['circs'].str.replace(r"\(.*\)","",regex=True)
	os.system('rm output_fasta.fa')
	fasted = pd.concat([fasted, fasted['circs'].str.split('::', expand=True).rename(columns={0:'circ_name',1:'position'})], axis=1).drop(['circs'],axis=1)
	circulars = list(fasted['circ_name'].unique())
	logging.info("Retrieving FASTA sequences")
	Data = {}
	for circ_name in circulars:
		Data['>'+circ_name] = ''.join(fasted.loc[fasted['circ_name']==circ_name]['sequence'].to_list())
	return Data, circulars, coord_final

def starting_from_coord(bedfile, rrna, genome):
	logging.info("Loading input bed file with coordinates: %s", str(bedfile))
	input_bed = pybedtools.BedTool(bedfile)
	bed_to_df = input_bed.to_dataframe()
	logging.info("Loading rRNA file: %s", str(rrna))
	rrna_file = pybedtools.BedTool(rrna)
	coordinates = input_bed.intersect(rrna_file, wa=True, v=True)
	coord_final = coordinates.to_dataframe()
	coordinates.sequence(fi=genome,fo='output_fasta.fa', tab=True, name=True, s=True)
	fasted = pd.read_csv('output_fasta.fa', sep='\t', names=['circs','sequence'])
	fasted['circs'] = fasted['circs'].str.replace(r"\(.*\)","",regex=True)
	os.system('rm output_fasta.fa')
	fasted = pd.concat([fasted, fasted['circs'].str.split('::', expand=True).rename(columns={0:'circ_name',1:'position'})], axis=1).drop(['circs'],axis=1)
	circulars = list(fasted['circ_name'].unique())
	logging.info("Retrieving FASTA sequences")
	Data = {}
	for circ_name in circulars:
		Data['>'+circ_name] = ''.join(fasted.loc[fasted['circ_name']==circ_name]['sequence'].to_list())
	return Data, circulars, coord_final

def getUTR(taxa):
	if taxa == 9606:
		return '3utr_human'
	elif taxa == 10090:
		return '3utr_human'
	elif taxa == 6239:
		return '3utr_worm'
	else:
		return '3utr_fly'

def third_party_processes(data_dict, taxa_number, cores, mature_file, acronym, miranda_sc, miranda_en, miranda_scale, miranda_go, miranda_ge, rnahybid_max):
	#Set the processes for the run of the three software
	logging.info("Setting the processes for the three software run up")
	DB = pd.DataFrame(columns = ['miRNA Name','Circ Name', 'Seed Category', 'Start', 'End', 'ID'])
	processes = []
	out_list = []
	file_path = getCircrPath()
	UTR = getUTR(taxa_number)
	for item in data_dict.keys():
		output = str(item)[1:]+'.txt'
		out_list.append(output)
		inp = str(item)[1:]+'_inputfile.fa'
		w = open(inp,'w')
		k = open('TS_'+inp,'w')
		#w = csv.writer(open(inp,'w'))
		#k = csv.writer(open('TS_'+inp,'w'))
		#w.writerow([item])
		w.write(item+'\n')
		#w.writerow([data_dict[item]])
		w.write(data_dict[item]+'\n')
		#k.writerow([item.lstrip('>')+'\t'+str(taxa_number)+'\t'+data_dict[item]])
		k.write(item.lstrip('>')+'\t'+str(taxa_number)+'\t'+data_dict[item]+'\n')
		w.close()
		k.close()
		processes.append('miranda {2} {0} -sc {3} -en {4} -scale {5} -go {6} -ge {7} -out {1}'.format(inp,'mir_TMP_'+output, mature_file, miranda_sc, miranda_en, miranda_scale, miranda_go, miranda_ge))						#<-- check input file for miranda
		processes.append('{0}/tools/targetscan_70.pl {0}/support_files/miRNA/miR_family_info.txt {1} {2}'.format(file_path,'TS_'+inp, 'TS_TMP_'+output))#<-- check support file for targetscan
		processes.append('RNAhybrid -m {4} -t {0} -q {2} -s {3} > {1}'.format(inp,'hyb_TMP_'+output, mature_file, UTR, rnahybid_max))	#<-- check support file for RNAhybrid

	#pool the analysis of the 3 software (RNAhybrid is the bottleneck) - only if there are more than one circular
	if len(out_list)==1:
		logging.info("Running miRNA:circRNA binding site predictions")
		for p in processes:
			run_process(p)
		logging.info("Parsing RNAhybrid's output files")
		hybrid(out_list[0])
		os.system('rm -rf *_inputfile.fa')
	else:
		logging.info("Running miRNA:circRNA binding site predictions")
		pool = Pool(processes=cores)
		pool.map(run_process, processes)
		pool.close()
		pool.join()
		#clean input and generate new pool for RNAhybrid
		os.system('rm -rf *_inputfile.fa')
		logging.info("Parsing RNAhybrid's output files")
		dpool = Pool(processes=cores)
		dpool.map(hybrid, out_list)
		dpool.close()
		dpool.join()
	#create a single dataframe with all predicted circulars from 3 software
	logging.info("Collecting results into a single location")
	for output in out_list:
		try:
			DB = examine(output, DB, acronym)
			DB = TS_parser(output, DB, acronym)
			DB = RH_parser(output, DB)
		except FileNotFoundError:
			logging.warning('File %s not found. Skipping', output)
			continue
	#return the generated database
	logging.info("Returning collected data")
	return DB

def compare_predicted(coordinates, circulars, database):
	tab = pd.DataFrame(columns = ['Chrom', 'Start', 'End', 'miRNA Name','Circ Name', 'Strand', 'Seed Category', 'ID'])
	for circ in circulars:
		if(coordinates.loc[coordinates['name'] == circ, 'strand'].max() == '-') and (len(coordinates.loc[coordinates['name'] == circ]) == 1):
			tab = pd.concat([tab, single_coord('-', coordinates, database, circ)])
		elif(coordinates.loc[coordinates['name'] == circ, 'strand'].max() == '-') and (len(coordinates.loc[coordinates['name'] == circ]) > 1):
			tab = pd.concat([tab, multiple_coord('-', coordinates, database, circ)])
		elif(coordinates.loc[coordinates['name'] == circ, 'strand'].max() == '+') and (len(coordinates.loc[coordinates['name'] == circ]) == 1):
			tab = pd.concat([tab, single_coord('+', coordinates, database, circ)])
		elif(coordinates.loc[coordinates['name'] == circ, 'strand'].max() == '+') and (len(coordinates.loc[coordinates['name'] == circ]) > 1):
			tab = pd.concat([tab, multiple_coord('+', coordinates, database, circ)])
	return tab

#Select miRanda's output when multiple softwares are able o find same interaction
#and report number of software that find that interaction
def filter_interactions(valInt, AGO, tab):
	validated_interactions = pybedtools.BedTool(valInt)
	AGO_interactions = pybedtools.BedTool(AGO)
	filtered = tab.merge(tab[tab.duplicated(subset=['miRNA Name', 'Circ Name'], keep=False)][['miRNA Name', 'Circ Name']], how='left', indicator=True)
	filtered = filtered[filtered['_merge'] == 'left_only'].drop(['_merge'],axis=1)
	filtered['Software Matched'] = 1
	duplicates = list(set(list(tab[tab.duplicated(subset=['miRNA Name', 'Circ Name'], keep=False)][['miRNA Name', 'Circ Name']].apply(tuple, axis=1))))
	to_add = []
	for i in duplicates:
		kkk = tab.loc[(tab['miRNA Name'] == i[0])&(tab['Circ Name'] == i[1])]
		kkk[['INT','Software','Number']] = kkk['ID'].str.split('_', expand=True)
		if (len(kkk['Software'].unique()) > 1) and ('M' in kkk['Software'].unique()):
			kkk['Software Matched'] = len(kkk['Software'].unique())
			to_add.append(kkk.loc[kkk['Software'] == 'M'].drop(['INT','Software','Number'], axis=1).values.tolist()[0])
		elif (len(kkk['Software'].unique()) == 1) and (len(kkk) > 1):
			kkk['Software Matched'] = len(kkk['Software'].unique())
			[to_add.append(x) for x in kkk.drop(['INT','Software','Number'], axis=1).values.tolist()]
	to_add = pd.DataFrame.from_records(to_add, columns=filtered.columns)
	filtered = pd.concat([filtered,to_add])
	filtered['Validated'] = 'No'
	filtered['AGO'] = 'No'
	filtered_bed = pybedtools.BedTool.from_dataframe(filtered)
	val_intersected = filtered_bed.intersect(validated_interactions, wa=True, wb=True, f=0.50).to_dataframe(names=list(range(0,len(filtered.columns) + len(validated_interactions.to_dataframe().columns))))
	AGO_intersected = filtered_bed.intersect(AGO_interactions, wa=True, wb=True, f=1).to_dataframe(names=list(range(0,len(filtered.columns) + len(AGO_interactions.to_dataframe().columns))))
	validated = val_intersected.loc[val_intersected[3] == val_intersected.iloc[:,-3]]
	filtered = add_interactions(filtered, validated, 'Validated')
	filtered = add_interactions(filtered, AGO_intersected, 'AGO')
	return filtered

#Add the notation from the CircBase database dump file
#from the supprot files
def add_CircBase_annotation(table, infile, genome_version):
	#get path of the CircBase ref file
	file_path = getCircrPath()
	CircBase = pd.read_csv(file_path+'/support_files/circBase_circRNA.txt', sep='\t', names=['Chrom', 'Start', 'End', 'CircName', 'Strand', 'Organism'])
	input_bed = pd.read_csv(infile, sep='\t', names=['Chrom', 'Start', 'End', 'CircName', 'Score', 'Strand'])
	#selecting only circulars of the same genome version investigated
	CircBase = CircBase.loc[CircBase['Organism'] == genome_version]
	#creating the reference column
	CircBase['ToBeMatched'] = CircBase['Chrom'] + CircBase['Start'].astype(str) + CircBase['End'].astype(str) + CircBase['Strand'].astype(str)
	# creating the reference dictionary from circBase
	CircBase_dict = dict(zip(CircBase.ToBeMatched, CircBase.CircName))
	match_dict = {}
	for ix, xg in input_bed.groupby(['Chrom', 'CircName', 'Strand']):
		xg = xg.sort_values(['Start'])
		start = xg['Start'].iloc[0]
		end = xg['End'].iloc[-1]
		# creating the string to be matched in circBase
		tbm = ix[0] + start.astype(str) + end.astype(str) + ix[2]
		# if there is a match in circBase
		if tbm in CircBase_dict:
			# add circBase ID
			match_dict[ix[1]] = CircBase_dict[tbm]
		else:
			# else, add 'Not Available'
			match_dict[ix[1]] = 'NA'

	# now update the output interaction file
	table['circBase ID'] = table['Circ Name']
	table['circBase ID'] = table['circBase ID'].map(match_dict)
	return table

############## END CODE #################
def main(bedfile, gtf, rrna, AGO, TAXA, outfile, Genome, cores, valid_interactions, mirna_mature, mirna_acronym, coord, genome_version,
		 miranda_sc, miranda_en, miranda_scale, miranda_go, miranda_ge, rnahybrid_max):  # valid_interactions
	begin = time.time()
	if coord:
		logging.info("Creating Circular Database and Retrieving Coordinates")
		circ_dict, circulars, coordinates = starting_from_coord(bedfile, rrna, Genome)
	else:
		logging.info("Retrieving exon information")
		data_exons, df_bed = calculate_exons(bedfile, gtf, rrna)
		logging.info("Creating Circular Database and Retrieving Coordinates")
		circ_dict, circulars, coordinates = create_circular_DB(data_exons, df_bed, rrna, Genome)
	logging.info("Performing analysis with RNAhybrid, miRanda and TargetScan")
	circular_database = third_party_processes(circ_dict, TAXA, cores, mirna_mature, mirna_acronym, miranda_sc, miranda_en, miranda_scale, miranda_go, miranda_ge, rnahybrid_max)
	logging.info("Merging overlapping seed predictions")
	table = compare_predicted(coordinates, circulars, circular_database)
	logging.info("Adding Annotation to Predicted Interactions")
	interactions = filter_interactions(valid_interactions, AGO, table)
	logging.info('Adding information from CircBase annotation')
	interactions = add_CircBase_annotation(interactions, bedfile, genome_version)
	logging.info('Printing results table in %s', outfile)
	interactions.to_csv(outfile, index=False)
	end = time.time() - begin
	logging.info('Analysis complete. Closing Circr')
	logging.info('Circr analysis performed in %s minutes', str(round(end / 60)))


if __name__ == "__main__":
	######## PARSER SECTION ###########
	parser = argparse.ArgumentParser(
	    prog='python3 Circr.py',
	    formatter_class=argparse.RawDescriptionHelpFormatter,
	    description="Circr help section")
	parser.add_argument("-i","--input", help="List of Circular RNA or their exon coordinates.")
	parser.add_argument("-c","--coord", action="store_true", help="Defines if the input file contains exon coordinates rather than Circular RNA coordinates")
	parser.add_argument("-s","--organism", type=str,help="Defines the selected organism for analysis. Available organisms are: human, mouse, worm, fruitfly. Default is human.", default="human")
	parser.add_argument("-v","--genome_version", type=str, help="Defines the genome version of the selected organism. Versions available are: human (hg19, hg38), mouse (mm9, mm10), worm (ce10, ce11), fruitfly (dm3, dm6). Default for human is hg38", default="none")
	parser.add_argument("--gtf", type=str,help="Alternative location for GTF file for the analysis. Default uses data for selected organism.", default="none")
	parser.add_argument("--genome", type=str,help="Alternative location for genome file for the analysis. Default uses data for selected organism.", default="none")
	parser.add_argument("--rRNA", type=str,help="Alternative location for ribosomal RNA file for the analysis. Default uses data for selected organism.", default="none")
	parser.add_argument("--miRNA", type=str,help="Alternative location for miRNA file for the analysis. Default uses data for selected organism.", default="none")
	parser.add_argument("--AGO", type=str, help="Alternative AGO peaks file. Default uses data for selected organism.", default="none")
	parser.add_argument("--validated_interactions", type=str, help="Alternative validated interactions file. Default uses data for selected organism.", default="none")
	parser.add_argument("--threads", type=int, help="Set the number of threads for multiprocess. Default is 8 cores", default=8)
	parser.add_argument("-o","--output", type=str,help="Defines output file name. Default is Circr_Analysis.csv", default="Circr_Analysis.csv")
	parser.add_argument("-Msc", type=float, help="Set score threshold, default 140. miRanda Parameter", default="140.0")
	parser.add_argument("-Men", type=float, help="Set energy threshold to -value Kcal/mol. Default 1. miRanda Parameter", default="1.0")
	parser.add_argument("-Mscale", type=float, help="Set scaling parameter. Default 4. miRanda Parameter", default="4.0")
	parser.add_argument("-Mgo", type=float, help="Set gap-open penalty to -value. Default 4 (-4). miRanda Parameter", default="-4.0")
	parser.add_argument("-Mge", type=float, help="Set gap-extend penalty to -value. Default 9 (-9). miRanda Parameter", default="-9.0")
	parser.add_argument("-RHmax", type=int, help="The  maximum  allowed  length of a target sequence. Default is 10000000. RNAhybrid Parameter", default="10000000")

	args = parser.parse_args()

	logging.basicConfig(filename='run.log', format='%(asctime)s - %(levelname)s:	%(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.INFO)

	### VARIABLE CHECKS ###
	if (args.genome_version == 'none') and (args.organism == 'human'):
		args.genome_version = 'hg38'

	if args.organism not in organisms.keys():
		logging.error("Organism not supported. Exiting...")
		sys.exit()

	if args.genome_version not in organisms[args.organism]:
		logging.error("Genome version not supported. Exiting...")
		sys.exit()

	### SETTING DEFAULTS ###

	if args.gtf == 'none':
		args.gtf = getCircrPath() + DataTree[args.organism][args.genome_version]['GTF']
	if args.genome == 'none':
		args.genome = getCircrPath() + DataTree[args.organism][args.genome_version]['Genome']
	if args.rRNA == 'none':
		args.rRNA = getCircrPath() + DataTree[args.organism][args.genome_version]['rRNA_Coords']
	if args.miRNA == 'none':
		args.miRNA = getCircrPath() + DataTree[args.organism]['file']
	if args.AGO == 'none':
		args.AGO = getCircrPath() + DataTree[args.organism][args.genome_version]['AGO']
	if args.validated_interactions == 'none':
		args.validated_interactions = getCircrPath() + DataTree[args.organism][args.genome_version]['INT']

	os.system('chmod +x {0}/tools/targetscan_70.pl'.format(getCircrPath()))

	logging.info('Running analysis on %s; genome set to %s and version %s', args.input, args.organism, args.genome_version)

	main(args.input,
		 args.gtf,
		 args.rRNA,
		 args.AGO,
		 DataTree[args.organism][args.genome_version]['TAXA'],
		 args.output,
		 args.genome,
		 args.threads,
		 args.validated_interactions,
		 args.miRNA, DataTree[args.organism]['acronym'],
		 args.coord,
		 args.genome_version,
		 args.Msc,
		 args.Men,
		 args.Mscale,
		 args.Mgo,
		 args.Mge,
		 args.RHmax)
