#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 13:12:35 2021

@author: ambuj
"""

import glob
import pandas as pd
import numpy as np
import os
import shutil
from Bio import SeqIO
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO
#cols = ['residue','pos']
class disorder_proteins:
    
    
    def __init__(self, cols=['residue','pos'],disorder_folder='DOTs/disorder_cont3_06April2021',
                 DOT_folder='DOTs/DOT_confirm08Mar2021',non_DOT_folder='DOTs/non_DOT08April2021',
                 issues_folder='DOTs/DOT_ISSUE08April2021', binding_residue_folder='DOTs/binding_residues08April2021_dist',
                 binding_DOT_folder='DOTs/DOT_BINDING10April2021_dist'):
        self.cols = cols
        self.disorder_folder = disorder_folder
        self.DOT_folder = DOT_folder
        self.non_DOT_folder = non_DOT_folder
        self.issues_folder = issues_folder
        self.binding_residue_folder = binding_residue_folder
        self.binding_DOT_folder = binding_DOT_folder
        return

    
    def f6_4_div (self,l1,l2):
        res = []
        for i in range(len(l1)):
            if l1[i] > l2[i]:
                val = float(l2[i]) / l1[i]
            else:
                val = float(l1[i]) / l2[i]
            res.append(val)
        return res
            
    def f6_3_diso_regions(self,key1,type1="complex"):
        dict1 = {'chain':'','DRs':0,'diso_str':''}
        dr = 0
        final_df = pd.DataFrame()
        if (type1 == "complex"):
            csv_file = pd.read_csv('complex_diso_05feb2021_csv/'+key1.split('|')[0]+'.csv')
        else:
            csv_file = pd.read_csv('free_diso_05feb2021_csv/'+key1.split('|')[0]+'.csv')
        if len(csv_file) ==0:
            return final_df
        chains = list(set(csv_file['chain']))
    #    print (chains)
        for ch in chains:
            dr = 0
            csv_file2 = csv_file[csv_file['chain'] == ch]
    #        print (len(csv_file2['pos']))
            diso_str = ''
            hold=-99999
            ctr=0
            for i in csv_file2['pos']:
    #            print ()
                
                if (i-hold) != 1:
                    dr += 1
                    diso_str += str(ctr)+'_'
                    ctr = 0
                ctr += 1
                hold=i
            diso_str += str(ctr)+'_'
            dict1['chain'] = ch
            dict1['DRs'] = dr
            dict1['diso_str'] = diso_str
            final_df = final_df.append(dict1,ignore_index=True)
        return
    #    return csv_file2
    #    print (csv_file)
        return final_df
# =============================================================================
#     def f6_2_calc_dot(self,comp_csv,free_csv,cmplx_pdb,cmplx_name,free_name):
#         DOT1 = pd.DataFrame()
#         DOT_issue = pd.DataFrame()
#         non_DOT1 = pd.DataFrame()
#     #    pos_error = []
#         atom_rec = [x for x in cmplx_pdb if x[0:4] == "ATOM"]
#     #    print (free_csv)
#     #    aa_list = [x[17:20] for x in atom_rec]
#     #    print (aa_list)
#     #    try:
#     #        pos_list = [int(x[22:28]) for x in atom_rec]
#     #    except ValueError:
#     #        pos_error.append(x[22:28])
#         for i in range(len(free_csv)):
#             if len(comp_csv) > 0:
#                 tmp1 = comp_csv[comp_csv['residue'] == free_csv['residue'][i]]
#                 tmp2 = tmp1[tmp1['pos'] == free_csv['pos'][i]]
#                 if (len(tmp2) > 0):
#     #                print (tmp2)
#                     non_DOT1 = non_DOT1.append(free_csv.iloc[i],ignore_index=True)
#                 else:
#                     try:
#     #                    tmp3 = [int(x[22:30].strip()) for x in atom_rec]
#     #                    print (free_csv.iloc[i])
#                         DOT1 = DOT1.append(free_csv.iloc[i],ignore_index=True)
#                         
#                     except ValueError:
#                         print (cmplx_pdb[0])
# #                        t = 0
#             else:
#                 DOT1 = DOT1.append(free_csv.iloc[i],ignore_index=True)
#     #            print ("no disorder in comp")
#     #    print (DOT1)
#         DOT_final = pd.DataFrame()
#         for j in range(len(DOT1)):
#             temp_dict = {'residue':'','pos':''}
#             res_sort = [x for x in  atom_rec if DOT1['residue'][j]==x[17:20]] 
#             b1 = [x for x in res_sort if str(int(DOT1['pos'][j]))==x[22:28].strip()]
#             if len(b1) > 0:
#                 temp_dict['residue'] = DOT1['residue'][j]
#                 temp_dict['pos'] = DOT1['pos'][j]
#                 DOT_final = DOT_final.append(temp_dict,ignore_index=True)
#             else:
#                 temp_dict['residue'] = DOT1['residue'][j]
#                 temp_dict['pos'] = DOT1['pos'][j]
#                 DOT_issue = DOT_issue.append(temp_dict,ignore_index=True)
#     #    print (DOT_final)
#         DOT_issue.to_csv(self.issues_folder+'/'+cmplx_name+'_'+free_name+'.txt')
#         DOT_final.to_csv('DOTs/DOT_confirm04Mar2021/'+cmplx_name+'_'+free_name+'.txt')
#         non_DOT1.to_csv('DOTs/non_DOT04Mar2021/'+cmplx_name+'_'+free_name+'.txt')
#         
#     #        print ([x for x in atom_rec if DOT1['residue'][j] in x] )
#     #    print (non_DOT1)
#         return DOT_final
# =============================================================================
    
    def f6_coverage2(self):
        complx_seq = SeqIO.index("blastdb/Proteins_complexes2.fasta", "fasta")
        free_seq = SeqIO.index("blastdb/free_protein_seq2.fasta", "fasta")
        
        for i in glob.glob('blastdb/final_result_04feb2021/*.out'):
        #    print (i)
            free_list = []
            df1 = pd.read_csv(i,delimiter='\t')
            cmplx = df1['query_acc'][0]
            cov2 = []
        #    diso_df1 = diso_regions(cmplx, "complex")
        #    diso_df1.to_csv('diso_summary_complex_06feb2021/'+cmplx.split('|')[0]+'.out')
#            diso_df1 = pd.read_csv('diso_summary_complex_06feb2021/'+cmplx.split('|')[0]+'.out')
#            complex_csv = pd.read_csv('complex_diso_05feb2021_csv/'+cmplx.split('|')[0].split('|')[0]+'.csv')
#            cmplx_pdb = open('complex_protein_str/'+cmplx.split('|')[0].split('|')[0].split('_')[0]+'.pdb').readlines()
            
            print (cmplx)
            len_cmplx = len(complx_seq[cmplx].seq)
        #    print (len_cmplx)
            
            for free in df1['sub_acc']:
                if free not in free_list:
                    len_free = len(free_seq[free].seq)
        #            print (len_free)
                    if len_cmplx > len_free:
                        div1 = float(len_free) / len_cmplx
                    else:
                        div1 = float(len_cmplx) / len_free
                    print (div1)
                    cov2.append(div1)
                    if div1 > 0.9:
                        free_list.append(free)
                    
        #            diso_df2 = diso_regions(free,"free")
        #            diso_df2.to_csv('diso_summary_free_06feb2021/'+free.split('|')[0]+'.out')
#                        diso_df2 = pd.read_csv('diso_summary_free_06feb2021/'+free.split('|')[0]+'.out')
#                        free_csv = pd.read_csv('free_diso_05feb2021_csv/'+free.split('|')[0].split('|')[0]+'.csv')
#                        DOT1 = self.f6_2_calc_dot(complex_csv,free_csv,cmplx_pdb,cmplx.split('|')[0],free.split('|')[0].split('|')[0])
        #        break
        #    break
            df1['coverage2'] = cov2
            df2 = df1[df1['coverage2'] > 0.9]
            if (len(df2) > 0):
                df2.to_csv('blastdb/final_result_04Mar2021/'+i.split('/')[-1], sep='\t')
        return
    def f7_2_calc_dot(comp_csv,free_csv,cmplx_pdb,cmplx_name,free_name):
        DOT1 = pd.DataFrame()
        DOT_issue = pd.DataFrame()
        non_DOT1 = pd.DataFrame()
        
    #    print (cmplx_pdb)
    #    pos_error = []
    #    atom_rec = [x for x in cmplx_pdb if x[0:4] == "ATOM"]
    #    print (free_csv)
    #    aa_list = [x[17:20] for x in atom_rec]
    #    print (aa_list)
    #    try:
    #        pos_list = [int(x[22:28]) for x in atom_rec]
    #    except ValueError:
    #        pos_error.append(x[22:28])
    #    print (free_csv)
        for i in range(len(free_csv)):
            if len(comp_csv) > 0:
    #            print (free_csv['residue'][i])
                tmp1 = comp_csv[comp_csv['residue'] == free_csv['residue'][i]]
                tmp2 = tmp1[tmp1['pos'] == free_csv['pos'][i]]
                if (len(tmp2) > 0):
    #                print (tmp2)
                    non_DOT1 = non_DOT1.append(free_csv.iloc[i],ignore_index=True)
                else:
                    try:
    #                    tmp3 = [int(x[22:30].strip()) for x in atom_rec]
    #                    print (free_csv.iloc[i])
                        DOT1 = DOT1.append(free_csv.iloc[i],ignore_index=True)
                        
                    except ValueError:
    #                    print (cmplx_pdb[0])
                        print ('error : '+str(i))
            else:
                DOT1 = DOT1.append(free_csv.iloc[i],ignore_index=True)
    #            print ("no disorder in comp")
    #    print (DOT1)
        DOT_final = pd.DataFrame()
        for j in range(len(DOT1)):
            temp_dict = {'residue':'','pos':''}
            try:
                resname = cmplx_pdb[int(DOT1['pos'][j])].get_resname()
    #            resname = cmplx_pdb[int(DOT1['pos'][j])].get_resname()
                
                if resname == DOT1['residue'][j]:
                    temp_dict['residue'] = DOT1['residue'][j]
                    temp_dict['pos'] = DOT1['pos'][j]
                    DOT_final = DOT_final.append(temp_dict,ignore_index=True)
                else:
                    temp_dict['residue'] = DOT1['residue'][j]
                    temp_dict['pos'] = DOT1['pos'][j]
                    DOT_issue = DOT_issue.append(temp_dict,ignore_index=True)
            except KeyError:
    #            print ('Y')
                temp_dict['residue'] = DOT1['residue'][j]
                temp_dict['pos'] = DOT1['pos'][j]
                DOT_issue = DOT_issue.append(temp_dict,ignore_index=True)
    #        res_sort = [x for x in  atom_rec if DOT1['residue'][j]==x[17:20]] 
    #        b1 = [x for x in res_sort if str(int(DOT1['pos'][j]))==x[22:28].strip()]
        print (cmplx_name+'_'+free_name)
    #    print (DOT_issue)
        DOT_issue.to_csv('DOTs/DOT_ISSUE04Mar2021/'+cmplx_name+'_'+free_name+'.txt')
        DOT_final.to_csv('DOTs/DOT_confirm04Mar2021/'+cmplx_name+'_'+free_name+'.txt')
        non_DOT1.to_csv('DOTs/non_DOT04Mar2021/'+cmplx_name+'_'+free_name+'.txt')
        
    #        print ([x for x in atom_rec if DOT1['residue'][j] in x] )
    #    print (non_DOT1)
        return DOT_final
    def f7_DOT_calc(self):
        complx_seq = SeqIO.index("blastdb/Proteins_complexes2.fasta", "fasta")
        free_seq = SeqIO.index("blastdb/free_protein_seq2.fasta", "fasta")
        
        for i in glob.glob('blastdb/final_result_04Mar2021/*.out'):
            print (i)
            free_list = []
            df1 = pd.read_csv(i,delimiter='\t')
            cmplx = df1['query_acc'][0]
            
            dict2 = {}
            diso_df1 = pd.read_csv('diso_summary_complex_06feb2021/'+cmplx.split('|')[0]+'.out')
            complex_csv = pd.read_csv('complex_diso_05feb2021_csv/'+cmplx.split('|')[0].split('|')[0]+'.csv')
            
            cmplx_ch = complx_seq[cmplx].description.split('|')[1].split()[1].split(',')[0]
            if len(complex_csv) > 0:
                complex_csv2 = complex_csv[complex_csv['chain'] == cmplx_ch]
                complex_csv2.reset_index(drop=True,inplace=True)
            else:
                complex_csv2 = pd.DataFrame()
            cmplx_pdb = open('complex_protein_str/'+cmplx.split('|')[0].split('|')[0].split('_')[0]+'.pdb').readlines()
            parser = PDBParser()
            structure = parser.get_structure('pdb', 'complex_protein_str/'+cmplx.split('|')[0].split('|')[0].split('_')[0]+'.pdb')
            cmplx_pdb = structure[0][cmplx_ch]
        #    print (len_cmplx)
            
            for free in df1['sub_acc']:
                if free not in free_list:
        
        #            print (free)
        
                    free_list.append(free)
                    
        #            diso_df2 = diso_regions(free,"free")
        #            diso_df2.to_csv('diso_summary_free_06feb2021/'+free.split('|')[0]+'.out')
                    diso_df2 = pd.read_csv('diso_summary_free_06feb2021/'+free.split('|')[0]+'.out')
                    free_csv = pd.read_csv('free_diso_05feb2021_csv/'+free.split('|')[0].split('|')[0]+'.csv')
                    free_ch = free_seq[free].description.split('|')[1].split()[1].split(',')[0]
                    free_csv2 = free_csv[free_csv['chain'] == free_ch]
                    free_csv2.reset_index(drop=True,inplace=True)
                    DOT1 = self.f7_2_calc_dot(complex_csv2,free_csv2,cmplx_pdb,cmplx.split('|')[0],free.split('|')[0].split('|')[0])
        return DOT1
    
    def f7_3_join_diso_DOT(self):
        
        for i in glob.glob('blastdb/final_result_04Mar2021/*.out'):
        #    fn1 = i.split('/')[-1].split('.')[0]
            all_diso = pd.read_csv(i,sep = '\t')
            for l1 in range(len(all_diso)):
                fn1 = all_diso['query_acc'][l1].split('|')[0]
                fn2 = all_diso['sub_acc'][l1].split('|')[0]
                fn3 = fn1+'_'+fn2
                
                DOT_confirm = pd.read_csv(self.DOT_folder+'/'+fn3+'.txt')
                non_DOT = pd.read_csv('DOTs/non_DOT04Mar2021/'+fn3+'.txt')
        #        free_csv = pd.read_csv('free_diso_05feb2021_csv/'+fn2+'.csv')
        #        issue_csv = pd.read_csv('DOTs/DOT_ISSUE04Mar2021/'+fn3+'.txt')
        #        df2=free_csv.append(issue_csv,ignore_index=True)
        #        index1 = set(df2.drop_duplicates(['pos','residue']).index)
        #        index2 = set(df2.drop_duplicates(['pos','residue'],keep='last').index)
        #        final_index = list(index1.intersection(index2))
        #        print(df2,final_index)
                df2=DOT_confirm.append(non_DOT,ignore_index=True)
        #        df3 = df2.loc[final_index]
                if (len(df2) > 0):
                    df2.sort_values(by=['pos'], inplace=True)
                    df3 = df2[self.cols]
                    df3.to_csv(self.disorder_folder+'/'+fn3)
                else:
                    df2.to_csv(self.disorder_folder+'/'+fn3)
        return
    def f8_2_get_regions(self,df1):
        temp_dict = {'no_of_regions':0,'lengths':[]}
        r1 = 1
        ctr = 1
#        print (df1)
        drlen = []
        if len(df1) == 1:
            r1 = 1
        else:
            for i in range(1,len(df1)):
                if (df1['pos'][i] - df1['pos'][i-1]) != 1:
                    r1 += 1
                    drlen.append (ctr)
                    ctr = 0
                ctr += 1
        drlen.append(ctr)
        temp_dict ['no_of_regions'] = r1
        temp_dict ['lengths'] = drlen
        return temp_dict
    
    def f8_combine_dataframes(self,folder, file1='DOTs/combined_dataframe.csv', file2='DOTs/regions_distribution.csv'):
        sumdf = pd.DataFrame()
        regions = pd.DataFrame()
        for i in glob.glob(folder+'/*'):
            df1 = pd.read_csv(i)
            
            sumdf = sumdf.append(df1,ignore_index=True)
            if (len(df1) >0):
                temp_dict1 = self.f8_2_get_regions (df1)
                regions = regions.append(temp_dict1,ignore_index = True)
        sumdf.to_csv(file1)
        regions.to_csv(file2)
        return
    
    def f9_2_molecule_type (self,seq1):
        rna = {'A': 0, 'G':0, 'C':0, 'U':0, 'I':0, 'R':0,'N':0,'D':0, 'Q':0, 'E':0, 'G':0, 'H':0, 'I':0, 'L':0, 'K':0, 'M':0,'F':0, 'P':0, 'S':0, 'T':0, 'W':0, 'Y':0, 'V':0,'X':0}
        for x in list(seq1):
            rna[x] += 1
        if rna['A']+rna['G']+rna['C']+rna['U']+rna['X'] > len(seq1)*0.7:
    #        print (all_seq[k].seq)
            mt = 'RNA'
        else:
            mt = 'Protein'
        return mt
    
    def f9_3_interaction_calc (self,cmplx_id, cmplx_ch, RNA_id, RNA_seqIO,dist):
        RNA_chs = RNA_seqIO [RNA_id].description.split('|')[1].split()[1].split(',')
        parser = PDBParser()
        structure = parser.get_structure('pdb', 'complex_protein_str/'+cmplx_id.split('_')[0]+'.pdb')
        cmplx_pdb = structure[0][cmplx_ch]
        contact_res = {'residue':'','pos':0}
        contact_df = pd.DataFrame()
        hold = []
        for chain1 in RNA_chs:
            rna_chain = structure[0][chain1]
            for res_cmplx in cmplx_pdb:
                for atom_cmplx in res_cmplx:
                    for res_rna in rna_chain:
                        for atom_rna in res_rna:
                            if (atom_cmplx-atom_rna) < dist:
                                if res_cmplx.get_resname()+'_'+str(atom_cmplx.get_full_id()[3][1]) not in hold:
                                    
                                    contact_res ['residue'] = res_cmplx.get_resname()
                                    contact_res ['pos'] = atom_cmplx.get_full_id()[3][1]
                                    contact_df = contact_df.append(contact_res,ignore_index=True)
                                    hold.append(contact_res ['residue']+'_'+str(contact_res ['pos']))
        
        return (contact_df)
    def f9_binding_residue(self, dist=3.5):
#        print ('Abc')
        complx_seqs = SeqIO.index("blastdb/Proteins_complexes2.fasta", "fasta")
        free_seqs = SeqIO.index("blastdb/free_protein_seq2.fasta", "fasta")
#        fasta_miss = []
        for i in glob.glob(self.non_DOT_folder+'/*.csv'):
#            print (i)
            df1 = pd.read_csv(i)
            cmplx_id = i.split('/')[-1].split('_')[0] + '_'+i.split('/')[-1].split('_')[1]
            free_id = i.split('/')[-1].split('_')[2] + '_'+i.split('/')[-1].split('_')[3].split('.')[0]
            if os.path.exists('DOTs/binding_residues08April2021_dist'+str(dist)+'A/'+cmplx_id+'.csv'):
                print (i)
                shutil.copyfile('DOTs/binding_residues08April2021_dist'+str(dist)+'A/'+cmplx_id+'.csv', 'Unique_complexes/binding_residues10April2021_dist'+str(dist)+'A/'+cmplx_id+'.csv')
            else:
                
                try:
                    cmplx_seq2 = complx_seqs[cmplx_id +'|Chain']
                except KeyError:
                    cmplx_seq2 = complx_seqs[cmplx_id +'|Chains']
                    
                try:
                    free_seq2 = free_seqs[free_id +'|Chain']
                except KeyError:
                    free_seq2 = free_seqs[free_id +'|Chains']
                    
                cmplx_ch = cmplx_seq2.description.split('|')[1].split()[1].split(',')[0]
                free_ch = free_seq2.description.split('|')[1].split()[1].split(',')[0]
            #    if os.path.exists('fasta_new2/'+cmplx_id.split('_')[0]+'.fasta'):
            #        a = 0
            #    else:
            #        if cmplx_id.split('_')[0] not in fasta_miss:
            #            fasta_miss.append (cmplx_id.split('_')[0])
                cmplx_name = cmplx_id.split('_')[0]
                for j in glob.glob('fasta_new2/'+cmplx_name+'*.fasta'):
            #        print (j)
                    if j.split('/')[-1].split('.')[0] != cmplx_id:
                        cmplx_seq3 = SeqIO.index(j, "fasta")
                        for k1 in cmplx_seq3.keys():
            #                print (k1)
                            mt = inst.f9_2_molecule_type (cmplx_seq3[k1].seq)
                            if mt == 'RNA':
                                contact_df = inst.f9_3_interaction_calc (cmplx_id,cmplx_ch, k1,cmplx_seq3,dist)
                                contact_df.to_csv('Unique_complexes/binding_residues10April2021_dist'+str(dist)+'A/'+cmplx_id+'.csv')
        return
    def f10_binding_DOT(self, dist=3.5):

#        cols = ['residue','pos']
        total = 0
        for i in glob.glob('DOTs/disorder_cont3_06April2021/*'):
            df1 = pd.read_csv(i)
#            df2 = pd.read_csv ('DOTs/DOT_confirm08April2021/'+i.split('/')[-1])
            if len(df1) > 0:
#                df3 = df1
        #        df3 = df3[cols]
        #        df3.append(df2,ignore_index=True)
                cmplx_id = i.split('/')[-1].split('_')[0] + '_'+i.split('/')[-1].split('_')[1]
                free_id = i.split('/')[-1].split('_')[2] + '_'+i.split('/')[-1].split('_')[3].split('.')[0]

#                cmplx_ch = cmplx_seq2.description.split('|')[1].split()[1].split(',')[0]
#                free_ch = free_seq2.description.split('|')[1].split()[1].split(',')[0]
                if os.path.exists(self.binding_residue_folder+str(dist)+'A/'+cmplx_id+'.csv'):
                    contact_df = pd.read_csv(self.binding_residue_folder+str(dist)+'A/'+cmplx_id+'.csv')
                else:
        #            print ('not exist: '+i)
                    contact_df = pd.DataFrame()
                if (len(contact_df) > 0):
                    dup1 = df1[df1.duplicated(['residue','pos'])]
                    if len(dup1) == 0:
                        df2 = df1.merge(contact_df, how='outer')
                        dup2 = df2[df2.duplicated(['residue','pos'])]
                        if len(dup2) >0:
        #                    print (df4)
                            total += len(dup2)
                            dup2[self.cols].to_csv(self.binding_DOT_folder+str(dist)+'A/'+cmplx_id+'_'+free_id+'.csv')
                    else:
                        print ('duplicate '+i)
#            break
        return

    def f11_segregate_cont3(self,source_folder='DOTs/disorder_confirm08April2021',target_folder='DOTs/disorder_cont3_06April2021'):
        
        n=3
        
        for i in glob.glob(source_folder+'/*'):
            result = False
            df1 = pd.read_csv(i)
            if len(df1)>2:
                df2 = df1[self.cols]
#                pos = [int(x[20:26]) for x in imp_remarks[1:]]
    
                for j in range (len(df2['pos'])-n+1):
#                    print (str(df2['pos'][n+j-1])+' '+str(df2['pos'][j]))
                    if df2['pos'][n+j-1]-df2['pos'][j] == n-1:
                        result = True
                
#                print (df2)
                        
            if result == True:
                shutil.copyfile(i,target_folder+'/'+i.split('/')[-1])
                print (i)
#        return result

    def f12_copy_files (self,source_folder,target_folder):
        for i in glob.glob(self.disorder_folder+'/*'):
            fn = i.split('/')[-1]
            shutil.copyfile(source_folder+'/'+fn+'.csv', target_folder+'/'+fn+'.csv')
        return
    def f13_extract_uniq_complex(self,source_folder,target_folder):
        cmplx_fn_list = []
        final_list = []
        for i in glob.glob(source_folder+'/*'):
            cmplx_fn = i.split('/')[-1].split('_')[0]
            
            if cmplx_fn not in cmplx_fn_list:
                cmplx_fn_list.append(cmplx_fn)
        for i in cmplx_fn_list:
            fns_byids = glob.glob(source_folder+'/'+i+'_*')
            if (len(fns_byids) == 1):
                final_list.append(fns_byids[0])
            else:
                temp = fns_byids[0]
                for j in fns_byids:
                    df1 = pd.read_csv(temp)
                    df2 = pd.read_csv(j)
                    if len(df2) >len(df1):
                        temp = j
                        
                final_list.append (temp)
#                break
        for f1 in final_list:
            shutil.copyfile(f1,target_folder+'/'+f1.split('/')[-1])
                
#        print (final_list)
        return final_list
    
inst = disorder_proteins()
inst2 = disorder_proteins(disorder_folder='Unique_complexes/disorder_cont3_10April2021',DOT_folder='Unique_complexes/DOT_confirm10Mar2021',non_DOT_folder='Unique_complexes/non_DOT10April2021',
                 issues_folder='Unique_complexes/DOT_ISSUE10April2021', binding_residue_folder='Unique_complexes/binding_residues10April2021_dist',
                 binding_DOT_folder='Unique_complexes/DOT_BINDING10April2021_dist')
inst.f8_combine_dataframes('DOTs/DOT_BINDING10April2021_dist6A',file1='DOTs/DOT_BINDING10April2021_df_6A.csv',file2='DOTs/DOT_BINDING10April2021_regions_6A.csv')
#inst2.f9_binding_residue(dist=6)
#inst2.f10_binding_DOT(dist=6)
#inst2.f8_combine_dataframes('Unique_complexes/binding_residues10April2021_dist3.5A',file1='Unique_complexes/binding_residues_df_3.5A.csv', file2='Unique_complexes/binding_region_3.5A.csv')
#inst.f9_binding_residue(dist=6)
#inst.f10_binding_DOT(dist=6)
#inst.f12_copy_files('DOTs/non_DOT08April2021', 'Unique_complexes/non_DOT10April2021', ref_folder='Unique_complexes/disorder_cont3_10April2021')
#uniq_cmplx = inst.f13_extract_uniq_complex()

            