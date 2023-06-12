#! /home/sujieyi/anaconda3/envs/python3.9/bin/python3.9
import re
import click
import glob
import os
import subprocess
import shutil
import sys
import string
import numpy as np
import pandas as pd
import shutil
import nr_pie_stat

@click.command()
@click.option('--type', default='nucl', help='nucl or prot, default = nucl')
@click.option('--evalue', default=1e-5, help='aligment evalue, default = 1e-5')
@click.option('--name', default='Unigene', help='prefix name, default = Unigene')
@click.option('--outdir', type=str, default='anno_result', help='Output directory for annotation results, default=anno_result')
@click.option('--mode', type=click.Choice(['PLANT','ANIMAL','MICRO']), required=True, help='PLANT ANIMAL or MICRO')
@click.option('--db', type=str, default='all', help='all nr swissprot kog eggnog kegg go pfam card cazy phi vfdb rna dna mir, default=all')
@click.option('--dbpath', type=str, default = "/public/mid/rna/denovo_database_2022/New/", help='path of database, default=/public/mid/rna/denovo_database_2022/New/')
@click.option('--cpu', type=int, default=2,help='cpu, default=2')
@click.option('--seq', help='input sequence')


def run(outdir, mode, db, seq, type, evalue, name, dbpath, cpu):
    outdir = os.path.abspath(outdir)
    dbpath = os.path.abspath(dbpath)
    db = db.split(',')
    if 'all' in db or 'nr' in db or 'rna' in db or 'dna' in db or 'mir' in db:
        if not os.path.exists(outdir + "/NR/"):
            os.mkdir(os.path.join(outdir,'NR'))
        alignment_diamond(os.path.join(dbpath,'NR',mode+'.dmnd'), os.path.join(outdir,'NR',name+'.NR.blast.best.xls'), type, seq, cpu, evalue)
        analysize_diamond(os.path.join(outdir,'NR',name+'.NR.blast.best.xls'), os.path.join(dbpath,'NR',mode+'_anno.txt'), os.path.join(outdir,'NR',name+'.NR.blast.best.anno.xls'), 'Database_ID')
        nr_pie_stat.nr_pie_stat(os.path.join(outdir,'NR',name+".NR.blast.best.anno.xls"), os.path.join(outdir,'NR','NR.species.top10.xls'), 10, os.path.join(dbpath,'oeanno','nr_pie_stat.r'))

    if 'all' in db or 'kegg' in db or 'rna' in db or 'dna' in db or 'mir' in db:
        if not os.path.exists(os.path.join(outdir,'KEGG')):
            os.makedirs(os.path.join(outdir,'KEGG'))
        alignment_diamond(os.path.join(dbpath,'kegg',mode+'.dmnd'), os.path.join(outdir,'KEGG',name+'.KEGG.blast.best.xls'), type, seq, cpu, evalue)
        analysize_KEGG(os.path.join(outdir,'KEGG',name+'.KEGG.blast.best.xls'),
                       os.path.join(dbpath,'kegg','KEGG_anno.txt'), os.path.join(outdir,'KEGG',name+'.KEGG.blast.best.anno.xls'),
                       os.path.join(outdir,'KEGG','kegg.backgroud.xls'),os.path.join(outdir,'KEGG','anno-kegg.backgroud.xls'), mode, dbpath, outdir, name)
        best_anno = pd.read_csv(os.path.join(outdir,'KEGG',name+'.KEGG.blast.best.anno.xls'),sep='\t')
        best_anno = best_anno[best_anno['Pathway'].notnull()][['#GeneID','KO','Pathway','Pathway_definition']]
        pathways= []
        for idx,row in best_anno.iterrows():
            pathway = row['Pathway'].split(',')
            pathway_definition = row['Pathway_definition'].split('|')
            for i in range(len(pathway)):
                row['Pathway'] = pathway[i]
                row['Pathway_definition'] = pathway_definition[i]
                pathways.append(row[['Pathway_definition','Pathway','#GeneID','KO']])
        outfile = pd.DataFrame(pathways)
        GeneID = outfile.groupby(['Pathway_definition','Pathway'])['#GeneID'].apply(';'.join).reset_index()
        KO = outfile.groupby(['Pathway_definition','Pathway'])['KO'].apply(';'.join).reset_index()
        kegg_pathway = pd.merge(GeneID,KO,on=['Pathway_definition','Pathway'])[['Pathway_definition','Pathway','#GeneID','KO']]
        kegg_pathway = pd.concat([pd.concat([kegg_pathway.iloc[:,0:2],kegg_pathway['KO'].str.split(';').apply(lambda x: len(x))],axis=1),kegg_pathway.iloc[:,2:4]],
                  axis=1)
        kegg_pathway.columns = ['Pathway_definition', 'Pathway', 'Gene_number', 'Gene_id', 'KOs']
        keggpathway_three_levels = pd.read_csv(os.path.join(dbpath,'kegg','KEGGpathway_three_levels_v2_191022.xls'),sep='\t',header=None,names=['Pathway',0,1,2,3])
        keggpathway_three_levels['#class'] = keggpathway_three_levels[0]+'--'+keggpathway_three_levels[1]
        KEGG_pathway = pd.merge(keggpathway_three_levels[['#class','Pathway']],kegg_pathway,on='Pathway')[['#class','Pathway_definition','Pathway','Gene_number','Gene_id','KOs']]
        KEGG_pathway = KEGG_pathway.sort_values(by='#class')
        KEGG_pathway.to_csv(os.path.join(outdir,'KEGG',name+'.KEGG.pathway.xls'),sep='\t',index=False)
        totalnum = pd.read_csv(os.path.join(outdir,'KEGG','kegg.backgroud.xls'),sep='\t').shape[0]
        print(keggpathway_three_levels)
        print(kegg_pathway)
        print (kegg_pathway.columns)
        print(keggpathway_three_levels['#class'].str.split('--', expand=True))
        kegg_pathway[['level1','level2']] = keggpathway_three_levels['#class'].str.split('--',expand=True)
        #kegg_pathway.drop('#class', axis=1, inplace=True)
        kegg_pathway = kegg_pathway[['level1','level2','Pathway_definition','Pathway','Gene_number','Gene_id','KOs']]
        gene = []
        for idx,row in kegg_pathway[['level2', 'level1', 'Gene_id']].iterrows():
            gene_id = row['Gene_id'].split(';')
            for i in range(len(gene_id)):
                row['Gene_id'] = gene_id[i]
                gene.append(row)
        kegg_path = pd.DataFrame(gene)
        kegg_path = kegg_path.drop_duplicates()
        kegg_path = kegg_path.groupby(['level2', 'level1'])['Gene_id'].apply(';'.join).reset_index()
        kegg_path = pd.concat([kegg_path,kegg_path['Gene_id'].str.split(';').apply(lambda x:len(x))],axis=1)
        kegg_path.columns = ['level2', 'level1','Gene_id','Gene_number']
        kegg_path = kegg_path[['level2', 'level1','Gene_number','Gene_id']]
        kegg_path = kegg_path.sort_values(by='level1',ascending=True)
        kegg_path['percentage'] = kegg_pathway['Gene_number']/totalnum
        kegg_path = kegg_path[['level2', 'level1','Gene_number','percentage','Gene_id']]
        kegg_path.columns = ['Classification_level2', 'Classification_level1','Gene_number','percentage','Genes']
        kegg_path.to_csv(os.path.join(outdir,'KEGG',name+'.KEGG.classification.xls'),sep='\t',index=False)
        cmd = 'cd %s && /home/sujieyi/anaconda3/envs/r4.0/bin/Rscript %s -i %s -o %s' % (os.path.join(outdir,'KEGG'),os.path.join(dbpath,'oeanno','KObarplot.r'),name+'.KEGG.classification.xls','./')
        subprocess.call(cmd, shell=True)

    if 'all' in db or 'swissprot' in db or 'rna' in db or 'dna' in db or 'mir' in db:
        if not os.path.exists(os.path.join(outdir, "Swissprot")):
            os.makedirs(os.path.join(outdir, "Swissprot"))
            alignment_diamond(os.path.join(dbpath,'Swissprot',mode+'.dmnd'), os.path.join(outdir,'Swissprot',name+'.Swissprot.blast.best.xls'), type, seq, cpu, evalue)
            analysize_diamond(os.path.join(outdir,'Swissprot',name+'.Swissprot.blast.best.xls'), os.path.join(dbpath,'Swissprot','SWISSPROT_anno.txt'),
                              os.path.join(outdir,'Swissprot',name+'.Swissprot.blast.best.anno.xls'), 'Database_ID')
    if  ('all' in db or 'swissprot' in db or 'rna' in db or 'dna' in db or 'mir' in db) and ('all' in db or 'go' in db or 'rna' in db or 'dna' in db or 'mir' in db): 
        if not os.path.exists(os.path.join(outdir,'GO')):
            os.makedirs(os.path.join(outdir,'GO'))
        df = pd.read_table(os.path.join(outdir,'Swissprot',name+'.Swissprot.blast.best.anno.xls'),sep='\t')
        df = df[['#GeneID','Database_ID']]
        anno = pd.read_table(os.path.join(dbpath,'GO','GO_anno.txt'))
        out = pd.merge(df, anno, left_on='Database_ID',right_on='id',how='left')
        out = out[['#GeneID','GO_ID','GO_Term','GO_Category']]
        out = out[out['GO_ID'].notnull()]
        go_backgroud = out[['#GeneID','GO_ID','GO_Term']]
        out.to_csv(os.path.join(outdir,'GO',name+'.GO.anno.xls'),index=False,sep='\t')
        go_backgroud.to_csv(os.path.join(outdir,'GO','go.backgroud.xls'),index=False,sep='\t',header=0)
        GO_anno = out
        GO_ID = []
        for idx,row in GO_anno.iterrows():
            ID = row['GO_ID'].split(',')
            for i in range(len(ID)):
                row['GO_ID'] = ID[i]
                GO_ID.append(row[['#GeneID','GO_ID']])
        go = pd.DataFrame(GO_ID)
        goanno = pd.read_table(os.path.join(dbpath,'GO','GO_classification.xls'),sep='\t')
        outclass = pd.merge(go, goanno, left_on='GO_ID',right_on='GO_id', how='left')
        outclass = outclass[outclass['GO_classify2'].notnull()]
        outclass = outclass.drop('GO_ID',1)
        outclass.to_csv(os.path.join(outdir,'GO',name+'.gene2go.ancestor.xls'),index=False,sep='\t')
        totalnum = GO_anno.shape[0]
        outclass = outclass[['GO_classify1','GO_classify2_definition','#GeneID']].drop_duplicates()
        outclass = outclass.groupby(['GO_classify1','GO_classify2_definition'])['#GeneID'].apply(';'.join).reset_index()
        outclass = pd.concat([outclass[['GO_classify1','GO_classify2_definition','#GeneID']],outclass['#GeneID'].str.split(';').apply(lambda x:len(x))], axis=1)
        outclass.columns = ['#GO_classify1','GO_classify2','Gene','Number']
        outclass = outclass.sort_values(by='#GO_classify1')
        outclass = outclass[['#GO_classify1','GO_classify2','Number','Gene']]
        first_add = pd.DataFrame(['#Total_gene', '', totalnum, '']).T
        first_add.columns = outclass.columns
        outclass = pd.concat([first_add,outclass], axis=0 ,ignore_index=True)
        outclass.to_csv(os.path.join(outdir,'GO',name+'.GO.classification.stat.xls'),sep='\t',index=False)
        outclass_tmp = outclass[['#GO_classify1','GO_classify2','Number']]
        outclass_tmp.to_csv(os.path.join(outdir,'GO',name+'.GO.classification.stat.xls_temp'),sep='\t',index=False)
        cmd = 'cd %s && /home/sujieyi/anaconda3/envs/r4.0/bin/Rscript %s --infile %s --outpath %s --fileName %s && rm %s ' % (os.path.join(outdir,'GO'),os.path.join(dbpath,'oeanno','GOClassificationMap.r'),
                                                                                                                    os.path.join(outdir,'GO',name+'.GO.classification.stat.xls_temp'),os.path.join(outdir,'GO'),
                                                                                                                    name+'.GO.classification.stat',os.path.join(outdir,'GO',name+'.GO.classification.stat.xls_temp'))
        subprocess.call(cmd, shell=True)

    if 'all' in db or 'kog' in db or 'rna' in db or 'dna' in db:
        if not os.path.exists(os.path.join(outdir,'KOG')):
            os.makedirs(os.path.join(outdir,'KOG'))
        alignment_diamond(os.path.join(dbpath,'KOG','kyva.dmnd'),os.path.join(outdir,'KOG',name+'.KOG.blast.best.xls'), type, seq, cpu, evalue)
        analysize_diamond(os.path.join(outdir,'KOG',name+'.KOG.blast.best.xls'), os.path.join(dbpath,'KOG','KOG_anno.txt'),os.path.join(outdir,'KOG',name+'.KOG.blast.best.anno.xls'),'KOG id')
        KOG_anno = pd.read_csv(os.path.join(outdir,'KOG',name+'.KOG.blast.best.anno.xls'),sep='\t').iloc[:,[0,8,9,10,6]]
        KOG_anno.iloc[:,1] = KOG_anno.iloc[:,1].apply(lambda x:x.replace('[','',1).replace(']','',1).replace('','|').strip('|'))
        print (KOG_anno)
        KOG = []
        KOG_anno.iloc[:,2].fillna(' ',inplace=True) ##出现空的情况,应该是X对应的元素
        KOG_anno.iloc[:,3].fillna(' ', inplace=True)
        for idx,row in KOG_anno.iterrows():
            function_code = row[1].split('|')
            Functional_categories = row[2].split(',')
            Function_class_defination = row[3].split('|')
            for i in range(len(function_code)):
                row[1] = function_code[i]
                row[2] = Functional_categories[i]
                row[3] = Function_class_defination[i]
                KOG.append(row.iloc[[2,3,1,0,4]])
        KOG_anno = pd.DataFrame(KOG)
        KOG_anno = KOG_anno.drop_duplicates(KOG_anno.columns[0:4].tolist())  ##这里往下的列名都用columns代替，因为不知道名字是什么
        a = KOG_anno.groupby(KOG_anno.columns[0:3].tolist())[KOG_anno.columns[3]].apply(';'.join).reset_index()
        b = KOG_anno.groupby(KOG_anno.columns[0:3].tolist())[KOG_anno.columns[4]].apply(';'.join).reset_index()
        KOG_anno = pd.merge(a,b,on=KOG_anno.columns[0:3].tolist())
        KOG_anno = KOG_anno.sort_values(by=KOG_anno.columns[2])
        KOG_anno = pd.concat([KOG_anno,KOG_anno[KOG_anno.columns[4]].str.split(';').apply(lambda x:len(x))], axis=1)
        col = KOG_anno.columns.tolist()
        col[-1] = 'gene number'
        KOG_anno.columns = col
        KOG_anno = KOG_anno.iloc[:,[0,1,2,5,3,4]]
        KOG_anno.columns = ['#Functional categories','Function class definition','Function code','gene number','gene id','KOG id']
        if 'X' in KOG_anno['Function code'].tolist():
            KOG_anno.drop(KOG_anno[KOG_anno['Function code'] == 'X'].index, axis=0, inplace=True)
        KOG_anno.to_csv(os.path.join(outdir,'KOG',name+'.KOG.class.xls'),sep='\t',index=False)
        Data = {}
        Class = {"J":"Translation, ribosomal structure and biogenesis","A" :"RNA processing and modification", "K":"Transcription", "L":"Replication, recombination and repair",
                 "B":"Chromatin structure and dynamics", "D":"Cell cycle control, cell division, chromosome partitioning", "Y":"Nuclear structure", "V":"Defense mechanisms",
                 "T":"Signal transduction mechanisms", "M":"Cell wall/membrane/envelope biogenesis", "N":"Cell motility", "Z":"Cytoskeleton", "W":"Extracellular structures",
                 "U":"Intracellular trafficking, secretion, and vesicular transport", "O":"Posttranslational modification, protein turnover, chaperones",
                 "C":"Energy production and conversion", "G":"Carbohydrate transport and metabolism", "E":"Amino acid transport and metabolism", "F":"Nucleotide transport and metabolism",
                 "H":"Coenzyme transport and metabolism", "I":"Lipid transport and metabolism", "P":"Inorganic ion transport and metabolism", "Q":"Secondary metabolites biosynthesis, transport and catabolism",
                 "R":"General function prediction only", "S":"Function unknown"}
        with open(os.path.join(outdir, 'KOG', name + '.KOG.blast.best.anno.xls'),'r') as f, open(os.path.join(outdir,'KOG',name+'.KOG'+'.cluster.stat'),'w') as outfile:
            outfile.write('#ID'+'\t'+'Class_Name'+'\t'+'Numbers'+'\n')
            for line in f:
                if line.startswith('#') or line.strip() == '':
                    continue
                l = line.split('\t')
                l[8] = l[8].replace('[','').replace(']','')
                for i in l[8]:
                    Data[i] = Data.get(i,0)+1
            Data = dict(sorted(Data.items(), key=lambda x:x[1]))
            if 'X' in Data.keys():
                Data.pop('X')
            for i,j in Data.items():
                outfile.write(i+'\t'+Class[i]+'\t'+str(j)+'\n')
        cmd = 'cd %s && /home/sujieyi/anaconda3/envs/r4.0/bin/Rscript %s %s %s' % (os.path.join(outdir,'KOG'), os.path.join(dbpath,'oeanno','kogcog_anno_plot.r'),
                                                                                   os.path.join(outdir,'KOG',name+'.KOG'+'.cluster.stat'),os.path.join(os.path.join(outdir,'KOG',name+'KOG.cluster')))
        subprocess.call(cmd, shell=True)

    if 'cog' in db or 'mir' in db :
        if not os.path.exists(os.path.join(outdir,'COG')):
            os.makedirs(os.path.join(outdir,'COG'))
        alignment_diamond(os.path.join(dbpath,'cog','myva.dmnd'), os.path.join(outdir,'COG',name+'.COG.blast.best.xls'), type, seq, cpu, evalue)
        analysize_diamond(os.path.join(outdir,'COG',name+'.COG.blast.best.xls'), os.path.join(dbpath,'cog','COG_anno.txt'),
                          os.path.join(outdir,'COG',name+'.COG.blast.best.anno.xls'), 'COG id')
        cog_blast = pd.read_csv(os.path.join(outdir,'COG',name+'.COG.blast.best.anno.xls'),sep='\t').iloc[:,[0,8,9,10,6]]
        cog_blast.iloc[:,1] = cog_blast.iloc[:,1].apply(lambda x:x.replace('[','').replace(']','').replace('','|').strip('|'))
        cog = []
        for idx,row in cog_blast.iterrows():
            f = row[1].split('|')
            b = row[2].split(',')
            c = row[3].split('|')
            for i in range(len(f)):
                row[1] = f[i]
                row[2] = b[i]
                row[3] = c[i]
                cog.append(row.iloc[[2,3,1,0,4]])
        cog_blast = pd.DataFrame(cog)
        cog_blast = cog_blast.drop_duplicates(cog_blast.columns[0:4].tolist())
        a = cog_blast.groupby(cog_blast.columns[0:3].tolist())[cog_blast.columns[3]].apply(';'.join).reset_index()
        b = cog_blast.groupby(cog_blast.columns[0:3].tolist())[cog_blast.columns[4]].apply(';'.join).reset_index()
        cog_blast = pd.merge(a, b, on=cog_blast.columns[0:3].tolist())
        cog_blast = cog_blast.sort_values(by=cog_blast.columns[2])
        cog_blast = pd.concat([cog_blast,cog_blast[cog_blast.columns[4]].str.split(';').apply(lambda x:len(x))], axis=1)
        cog_blast = cog_blast.iloc[:,[0,1,2,5,3,4]]
        cog_blast.columns = ['#Functional categories', 'Function class definition', 'Function code', 'gene number', 'gene id', 'COG id']
        if 'X' in cog_blast['Function code'].tolist():
            cog_blast.drop(cog_blast[cog_blast['Function code'] == 'X'].index, axis=0, inplace=True)
        cog_blast.to_csv(os.path.join(outdir,'COG',name+'.COG.class.xls'), sep='\t', index=False)
        Data = {}
        Class = {"J": "Translation, ribosomal structure and biogenesis", "A": "RNA processing and modification",
                 "K": "Transcription", "L": "Replication, recombination and repair",
                 "B": "Chromatin structure and dynamics",
                 "D": "Cell cycle control, cell division, chromosome partitioning", "Y": "Nuclear structure",
                 "V": "Defense mechanisms",
                 "T": "Signal transduction mechanisms", "M": "Cell wall/membrane/envelope biogenesis",
                 "N": "Cell motility", "Z": "Cytoskeleton", "W": "Extracellular structures",
                 "U": "Intracellular trafficking, secretion, and vesicular transport",
                 "O": "Posttranslational modification, protein turnover, chaperones",
                 "C": "Energy production and conversion", "G": "Carbohydrate transport and metabolism",
                 "E": "Amino acid transport and metabolism", "F": "Nucleotide transport and metabolism",
                 "H": "Coenzyme transport and metabolism", "I": "Lipid transport and metabolism",
                 "P": "Inorganic ion transport and metabolism",
                 "Q": "Secondary metabolites biosynthesis, transport and catabolism",
                 "R": "General function prediction only", "S": "Function unknown"}
        with open(os.path.join(outdir, 'COG', name + '.COG.blast.best.anno.xls'), 'r') as f, open(os.path.join(outdir,'COG',name+'.COG'+'.cluster.stat'),'w') as outfile:
            outfile.write('#ID'+'\t'+'Class_Name'+'\t'+'Numbers'+'\n')
            for line in f:
                if line.startswith('#') or line.strip() == '':
                    continue
                l = line.split('\t')
                l[8] = l[8].replace('[','').replace(']','')
                for i in l[8]:
                    Data[i] = Data.get(i,0)+1
            Data = dict(sorted(Data.items(), key=lambda x:x[1]))
            if 'X' in Data.keys():
                Data.pop('X')
            for i,j in Data.items():
                outfile.write(i+'\t'+Class[i]+'\t'+str(j)+'\n')
        cmd = 'cd %s && /home/sujieyi/anaconda3/envs/r4.0/bin/Rscript %s %s %s' % (os.path.join(outdir,'COG'), os.path.join(dbpath,'oeanno','kogcog_anno_plot.r'),
                                                                                   os.path.join(outdir,'COG',name+'.COG'+'.cluster.stat'),os.path.join(os.path.join(outdir,'COG',name+'COG.cluster')))
        subprocess.call(cmd, shell=True)

    if 'all' in db or 'eggnog' in db or 'rna' in db or 'dna' in db or 'mir' in db:
        if not os.path.exists(os.path.join(outdir,'eggNOG')):
            os.makedirs(os.path.join(outdir,'eggNOG'))
        #alignment_diamond(os.path.join(dbpath,'eggNOG','eggNOG.dmnd'), os.path.join(outdir,'eggNOG',name+'.eggNOG.blast.best.xls'), type, seq, cpu, evalue)
        analysize_diamond(os.path.join(outdir,'eggNOG',name+'.eggNOG.blast.best.xls'), os.path.join(dbpath,'eggNOG','eggNOG_anno.txt'),
                          os.path.join(outdir,'eggNOG',name+'.eggNOG.blast.best.anno.xls_temp'), 'eggNOG_ID')
        egg_blast = pd.read_csv(os.path.join(outdir,'eggNOG',name+'.eggNOG.blast.best.anno.xls_temp'),sep='\t')
        dfegg = egg_blast.drop_duplicates(subset='#GeneID',keep='first',inplace=False)
        dfegg.to_csv(os.path.join(outdir,'eggNOG',name+'.eggNOG.blast.best.anno.xls'), sep='\t', index=False)
        dfegg = dfegg.iloc[:,[0,6]]
        dfegg.iloc[:,1] = dfegg.iloc[:,1].apply(lambda x:x.replace('','|').strip('|'))
        egg = []
        for idx,row in dfegg.iterrows():
            f = row[1].split('|')
            for i in range(len(f)):
                row[1] = f[i]
                egg.append(row.iloc[[1,0]])
        egg = pd.DataFrame(egg).drop_duplicates()
        egg = egg.groupby(egg.columns[0])[egg.columns[1]].apply(';'.join).reset_index().sort_values(by=egg.columns[0])
        egg = pd.concat([egg,egg[egg.columns[1]].str.split(';').apply(lambda x:len(x))], axis=1)
        egg = egg.iloc[:,[0,2]]
        egg.columns = ['eggNOG','number']
        eggNOG_summary = pd.merge(egg, pd.read_csv(os.path.join(dbpath,'eggNOG','fun.xls'),sep='\t'), on='eggNOG', how='right').fillna(0)
        print (eggNOG_summary)
        eggNOG_summary.to_csv(os.path.join(outdir,'eggNOG',name+'.eggNOG.summary.xls'),sep='\t',index=False)
        cmd = 'cd %s && /home/sujieyi/anaconda3/envs/r4.0/bin/Rscript %s %s %s ' % (os.path.join(outdir,'eggNOG'), os.path.join(dbpath,'oeanno','NOG.bar.R'),
                                                                                 os.path.join(outdir,'eggNOG',name+'.eggNOG.summary.xls'), os.path.join(outdir,'eggNOG',name+'.eggNOG'))
        subprocess.call(cmd, shell=True)
        os.remove(os.path.join(outdir,'eggNOG',name+'.eggNOG.blast.best.anno.xls_temp'))

    if 'all' in db or 'pfam' in db or 'rna' in db or 'dna' in db or 'mir' in db:
        if not os.path.exists(os.path.join(outdir,'Pfam')):
            os.makedirs(os.path.join(outdir,'Pfam'))
        alignment_hmmer(seq, os.path.join(dbpath,'pfam','Pfam-A.hmm'), os.path.join(outdir,'Pfam'), os.path.join(outdir,'Pfam',name+'.Pfam.align.txt'), type, cpu, evalue)
        shutil.rmtree(os.path.join(outdir,'Pfam.__checkpoints_longorfs'))
        with open(os.path.join(outdir, 'Pfam', name + '.Pfam.align.txt'), 'r') as file:
            lines = file.readlines()
            lines = [line.strip() for line in lines if not line.startswith('#')]
        output_lines = []
        for line in lines:
            line = line.split(None)
            out = line[2]+'\t'+line[1]
            for i in range(18,len(line)):
                out += ' '+line[i]
            output_lines.append(out.replace(' ','\t',1))
        result = []
        if type == 'nucl':
            # for line in output_lines:
            #     line = re.split('::|\t', line)
            #     result.append([line[1],line[4],line[5]])
            # result = pd.DataFrame(result)
            tmp = pd.DataFrame([x.split('\t') for x in output_lines],columns=['#Gene_ID', 'Pfam_IDs', 'Pfam_Description']).drop_duplicates(['#Gene_ID', 'Pfam_IDs'])
            result = pd.concat([tmp.groupby('#Gene_ID')['Pfam_IDs'].apply(','.join),tmp.groupby('#Gene_ID')['Pfam_Description'].apply('|'.join)], axis=1).reset_index()
        elif type == 'prot':
            for line in output_lines:
                line = line.split('\t')
                result.append([line[0],line[1],line[2]])
            result = pd.DataFrame(result)
        # result = result.drop_duplicates(result.columns[0:2].tolist())
        # a = result.groupby(result.columns[0])[result.columns[1]].apply(','.join).reset_index()
        # b = result.groupby(result.columns[0])[result.columns[2]].apply('|'.join).reset_index()
        # result = pd.merge(a,b,on=a.columns[0])
        # result = result.sort_values(by=result.columns[0])
        # result.columns = ['#Gene_ID','Pfam_IDs','Pfam_Description']
        result.to_csv(os.path.join(outdir,'Pfam',name+'.Pfam.align.anno.xls'), sep='\t', index=False)

    if 'all' in db or 'cazy' in db or 'mir' in db:
        if not os.path.exists(os.path.join(outdir,'CAZy')):
            os.makedirs(os.path.join(outdir,'CAZy'))
        alignment_hmmer(seq, os.path.join(dbpath,'CAZy','dbCAN-HMMdb-V8.txt'), os.path.join(outdir,'CAZy'), os.path.join(outdir,'CAZy',name+'.CAZy.align.txt'), type, cpu, evalue)
        shutil.rmtree(os.path.join(outdir, 'CAZy.__checkpoints_longorfs'))
        with open(os.path.join(outdir,'CAZy',name+'.CAZy.align.txt'),'r') as file:
            lines = file.readlines()
            lines = [line.strip() for line in lines if not line.startswith('#')]
        output_lines = []
        for line in lines:
            line = line.split(None)
            output_lines.append([line[2], line[7], line[0]])
        output_lines = pd.DataFrame(output_lines)
        result = output_lines.drop_duplicates(output_lines.columns[0])
        result[2] = result[2].apply(lambda x:x.replace('.hmm',''))
        if type == 'nucl':
            result = pd.concat([result.iloc[:,0:2],result.iloc[:,2].str.split('_', expand=True)], axis=1)
            result = result.iloc[:,[0,2,1]]
            result.columns = ['#GeneID','Family','Evalue']
            #result.to_csv(os.path.join(outdir,'CAZy',name+'.CAZy.tmp.xls'), sep='\t', index=False)
        elif type == 'prot':
            result = pd.concat([result.iloc[:,0:2],result.iloc[:,2].str.split('_',expand=True)], axis=1)
            result = result.iloc[:,[0,2,1]]
            result.columns = ['#Gene_ID','Family','Evalue']
            #result.to_csv(os.path.join(outdir,'CAZy',name+'.CAZy.tmp.xls'), sep='\t', index=False)
        anno = pd.read_csv(os.path.join(dbpath,'CAZy','CAZy_anno.txt'),sep='\t')
        out = pd.merge(result, anno, left_on='Family', right_on='id', how='left')
        out = out.drop('id',1)
        out.to_csv(os.path.join(outdir,'CAZy',name+'.CAZy.anno.xls'), sep='\t', index=False)
        out = out.iloc[:,[1,3,4,0]]  ###代码写到awk -F"\t" -v OFS="\t" \'{print $2,$4,$5,$1}\' %s这里
        out = out.groupby(out.columns.tolist()[0:3])[out.columns[3]].apply(','.join).reset_index()
        out = pd.concat([out,out[out.columns[3]].str.split(',').apply(lambda x:len(x))], axis=1)
        out = out.sort_values(by=out.columns[0])
        out = out.iloc[:,[0,1,2,4,3]]
        out.columns = ['Family','Class','Class_Definition','Genes_Count','Genes_List']
        out.to_csv(os.path.join(outdir,'CAZy',name+'.CAZy.family.xls'), sep='\t', index=False)
        out1 = pd.read_csv(os.path.join(outdir,'CAZy',name+'.CAZy.anno.xls'), sep='\t')
        out1 = out1.iloc[:,[3,4,0]]
        out1 = out1.groupby(out1.columns.tolist()[0:2])[out1.columns[2]].apply(','.join).reset_index()
        out1 = pd.concat([out1,out1[out1.columns[-1]].str.split(',').apply(lambda x:len(x))], axis=1)
        out1 = out1.sort_values(by=out1.columns[0])
        out1 = out1.iloc[:,[0,1,3,2]]
        out1.columns = ['Class','Class_Definition','Genes_Count','Genes_List']
        out1.to_csv(os.path.join(outdir,'CAZy',name+'.CAZy.class.xls'), sep='\t', index=False)
        cmd = 'cd %s && /home/sujieyi/anaconda3/envs/r4.0/bin/Rscript %s -i %s -o %s -n %s' % (os.path.join(outdir,'CAZy'), os.path.join(dbpath,'oeanno','cazy_anno.r'),
                                                                                    os.path.join(outdir,'CAZy',name+'.CAZy.class.xls'), os.path.join(outdir,'CAZy'), name)
        subprocess.call(cmd, shell=True)

    if 'all' in db or 'card' in db or 'mir' in db:
        if not os.path.exists(os.path.join(outdir,'CARD')):
            os.makedirs(os.path.join(outdir,'CARD'))
        alignment_diamond(os.path.join(dbpath,'CARD','CARD.dmnd'), os.path.join(outdir,'CARD',name+'.CARD.blast.best.xls'), type, seq, cpu, evalue)
        analysize_diamond(os.path.join(outdir,'CARD',name+'.CARD.blast.best.xls'), os.path.join(dbpath,'CARD','CARD_anno.txt'), os.path.join(outdir,'CARD',name+'.CARD.anno.xls'), 'Description')
        H = {}
        with open(os.path.join(outdir,'CARD',name+'.CARD.anno.xls'),'r') as file:
            for line in file:
                if not line.startswith('#') or line.strip() == '':
                    line = line.split('\t')
                    Name = line[-1].strip('\n')
                    H[Name] = H.get(Name,0)+1
        H = dict(sorted(H.items(), key=lambda x:x[1], reverse=True))
        limit_max = 10
        count = {}
        n = 0
        for i,j in H.items():
            if n < limit_max:
                count[i] = H[i]
                n += 1
            else:
                count['other'] = count.get('other',0)+H[i]
        print(count)
        with open(os.path.join(outdir,'CARD',name+'.statics.txt'),'w') as out:
            out.write('AROID'+'\t'+'Num'+'\n')
            for key in count.keys():
                out.write(key+'\t'+str(count[key])+'\n')
        cmd = 'cd %s && /home/sujieyi/anaconda3/envs/r4.0/bin/Rscript %s -i %s ' % (os.path.join(outdir,'CARD'), os.path.join(dbpath,'oeanno','card_pie_stat.r'),
                                                                              os.path.join(outdir,'CARD',name+'.statics.txt'))
        subprocess.call(cmd, shell=True)

    if 'all' in db or 'phi' in db or 'mir' in db:
        if not os.path.exists(os.path.join(outdir, 'PHI')):
            os.makedirs(os.path.join(outdir,'PHI'))
        alignment_diamond(os.path.join(dbpath,'PHI','phi_accessions.dmnd'), os.path.join(outdir,'PHI',name+'.PHI.blast.best.xls'), type, seq, cpu, evalue)
        analysize_diamond(os.path.join(outdir,'PHI',name+'.PHI.blast.best.xls'), os.path.join(dbpath,'PHI','PHI_anno.txt'), os.path.join(outdir,'PHI',name+'.PHI.anno.xls'), 'PHI_MolConn_ID')

    if 'all' in db or 'vfdb' in db or 'mir' in db:
        if not os.path.exists(os.path.join(outdir,'VFDB')):
            os.makedirs(os.path.join(outdir,'VFDB'))
        alignment_diamond(os.path.join(dbpath,'VFDB','VFDB_setB_pro.dmnd'), os.path.join(outdir,'VFDB',name+'.VFDB.blast.best.xls'), type, seq, cpu, evalue)
        analysize_diamond(os.path.join(outdir,'VFDB',name+'.VFDB.blast.best.xls'), os.path.join(dbpath,'VFDB','VFDB_anno.txt'), os.path.join(outdir,'VFDB',name+'.VFDB.anno.xls'), 'Database_ID')


def alignment_diamond(database_diamond, outfile, type, seq, cpu, evalue):
    if type == 'nucl':
        cmd ='/home/rna/softwares/diamond/diamond-2.0.13/bin/diamond blastx -q %s -d %s -p %s -k 1 --more-sensitive --outfmt 6 -e %s -o %s' % (seq, database_diamond, cpu, evalue, outfile)
        print (cmd)
    elif type == 'prot':
        cmd = '/home/rna/softwares/diamond/diamond-2.0.13/bin/diamond blastp -q %s -d %s -p %s -k 1 --more-sensitive --outfmt 6 -e %s -o %s' % (seq, database_diamond, cpu, evalue, outfile)
        print (cmd)
    else:
        print ('please check the seqence type!')
        sys.exit()
    subprocess.call(cmd, shell=True)

def alignment_hmmer(seq, database_hmmer, tmpdir, outfile, type, cpu, evalue):
    if type == 'nucl':
        shutil.copy(seq, tmpdir)
        (seq, filename) = os.path.split(os.path.abspath(seq))
        cmd = '/data/software/perl/perl-v5.24.1/bin/perl  /data/software/TransDecoder/v5.5.0/TransDecoder.LongOrfs -m 20 -t %s --output_dir %s' % (filename, tmpdir)
        cmd +='&& /home/fanyucai/software/hmmer/hmmer-v3.1b2/bin/hmmscan --cpu %s --acc --notextw -E %s --tblout %s %s %s' % (cpu, evalue, outfile, database_hmmer, os.path.join(tmpdir,'longest_orfs.pep'))
    elif type == 'prot':
        cmd = '/home/fanyucai/software/hmmer/hmmer-v3.1b2/bin/hmmscan --cpu %s --acc --notextw -E %s --tblout %s %s %s' % (cpu, evalue, outfile, database_hmmer, seq)
    else:
        print ('please check the seqence type!')
        sys.exit()
    subprocess.call(cmd, shell=True)

def analysize_diamond(inputfile, annofile, outfile, keyword):
    df = pd.read_table(inputfile, header=None, names=['#GeneID','Database_ID','Identity','alignment length','mismatches','gap openings','q. start','q. end','s. start','s. end','E_value','Score'])
    df = df[['#GeneID','Database_ID','E_value','Identity','Score']]
    anno = pd.read_table(annofile)
    out = pd.merge(df, anno, left_on='Database_ID',right_on='id',how='left')
    out = out.drop('id', 1)
    out = out[out[keyword].notnull()]
    out.to_csv(outfile, index=False,sep='\t')

def analysize_KEGG(inputfile, annofile, outfile, kegg_backfile, anno_keggfile, mode, dbpath, outdir, name):
    df = pd.read_table(inputfile, header=None, names=['#GeneID','Database_ID','Identity','alignment length','mismatches','gap openings','q. start','q. end','s. start','s. end','E_value','Score'])
    df = df[['#GeneID','Database_ID','E_value','Identity']]
    anno = pd.read_table(annofile)
    out = pd.merge(df, anno, left_on='Database_ID',right_on='id', how='left')
    out = out.drop('id', 1)
    out = out[out['KO'].notnull()]
    out.to_csv(outfile+"tmp", index=False, sep='\t')
    if mode in ['PLANT','ANIMAL']:
        out_file = out
        pathways = []
        for idx, row in out_file.iterrows():
            pathway_list = row['Pathway'].split(',')
            pw_definition = row['Pathway_definition'].split('|')
            for i in range(len(pathway_list)):
                new_row = row.copy()
                new_row['Pathway'] = pathway_list[i]
                new_row['Pathway_definition'] = pw_definition[i]
                pathways.append(new_row)
        out_file = pd.DataFrame(pathways)
        KEGG_three_levels = pd.read_csv(os.path.join(dbpath,'kegg','KEGGpathway_three_levels_v2_191022.xls'),sep='\t',names=['Pathway',0,1,2,'type'])
        tmp = pd.merge(out_file, KEGG_three_levels[['Pathway','type']], on='Pathway')
        tmp = tmp[tmp['type'].str.contains(mode.lower())].drop('type',1)
        Pathway =tmp.groupby(['#GeneID', 'Database_ID', 'E_value', 'Identity', 'KO', 'gene_name', 'description'])['Pathway'].apply(','.join).reset_index()
        Pathway_definition = tmp.groupby(['#GeneID', 'Database_ID', 'E_value', 'Identity', 'KO', 'gene_name', 'description'])['Pathway_definition'].apply('|'.join).reset_index()
        kegganno = pd.merge(Pathway,Pathway_definition,on=['#GeneID', 'Database_ID', 'E_value', 'Identity', 'KO', 'gene_name', 'description'])
        kegganno.to_csv(outfile, sep='\t', index=False)
    else:
        os.rename(outfile+"tmp", outfile)
        kegganno = out
    kegg_backgroud = kegganno[['#GeneID','Pathway','Pathway_definition']]
    kegg_backgroud = kegg_backgroud[kegg_backgroud['Pathway'].notnull()]
    kegg_backgroud.to_csv(kegg_backfile, index=False, sep='\t',header=0)
    anno_kegg_backgroud = kegganno[['#GeneID','KO','Pathway']]
    anno_kegg_backgroud = anno_kegg_backgroud[anno_kegg_backgroud['Pathway'].notnull()]
    anno_kegg_backgroud.to_csv(anno_keggfile, index=False, sep='\t', header=0)
    
if __name__ == '__main__':
    run()















