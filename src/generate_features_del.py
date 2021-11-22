import vcf
import csv
import pyBigWig
import numpy as np
import os
import pysam
import time
import multiprocessing
import argparse
from sklearn import preprocessing
from functools import partial
import random

features_head = ['Location', 'ID', 'SV_LEN', 'CADD_RAW', 'CADD_PHRED', 'GDI_score', 'GDI_Phred', 'Lof_score', 'MCAP_score', 'DamagePredCount_mean', 'vest_score_mean',
                 'MVP_score_mean', 'MPC_score_mean', 'DANN_score_mean', 'GenoCanyon_score_mean', 'Integrated_fitCons_score_mean', 'GM12878_fitCons_score_mean',
                 'H1_hESC_fitCons_score_mean', 'HUVEC_fitCons_score_mean', 'LINSIGHT_mean', 'GERP_NR_mean', 'GERP_RS_mean', 'PhyloP100way_vertebrate_mean',
                 'PhyloP30way_mammalian_mean', 'PhyloP17way_primate_mean', 'PhastCons100way_vertebrate_mean', 'PhastCons30way_mammalian_mean', 'PhastCons17way_primate_mean',
                 'RVIS_ExAC_005_1', 'RVIS_ExAC_005_2', 'REVEL_score', 'Roadmap_SRS004212', 'E003_DNase', 'E003_H3K27ac', 'E003_H3K27me3', 'E003_H3K36me3',
                 'E003_H3K4me1_signal', 'E003_H3K4me3_signal', 'E003_H3K9me3_signal', 'E003_Methylation_signal', 'ENCFF225MAO_signal',
                 'EncodeUwRepliSeqGm12878WaveSignal_signal', 'phyloP100way_signal', '3utr_overlap', '5utr_overlap', 'CDS_overlap', 'Prom_overlap', 'SS_overlap',
                 'TAD_overlap', 'sensitive_overlap', 'conserved_overlap', 'Heterochromatin_overlap',
                 'All_disease_causing_High', 'All_disease_causing_Low', 'All_disease_causing_Medium', 'All_Mendelian_High', 'All_Mendelian_Low', 'All_Mendelian_Medium',
                 'Mendelian_AD_High', 'Mendelian_AD_Low', 'Mendelian_AD_Medium', 'Mendelian_AR_High', 'Mendelian_AR_Low', 'Mendelian_AR_Medium',
                 'all_PID_High', 'all_PID_Low', 'all_PID_Medium', 'PID_AD_High', 'PID_AD_Low', 'PID_AD_Medium', 'PID_AR_High', 'PID_AR_Low', 'PID_AR_Medium',
                 'All_cancer_High', 'All_cancer_Low', 'All_cancer_Medium', 'Cancer_recessive_High', 'Cancer_recessive_Low', 'Cancer_recessive_Medium',
                 'Cancer_dominant_High', 'Cancer_dominant_Low', 'Cancer_dominant_Medium',
                 'ExonicFunc.refGene_frameshift_deletion', 'ExonicFunc.refGene_nonframeshift_deletion', 'ExonicFunc.refGene_startloss',
                 'ExonicFunc.refGene_stopgain', 'ExonicFunc.refGene_stoploss', 'ExonicFunc.refGene_unknown',
                 'SIFT_pred_.', 'SIFT_pred_D', 'SIFT_pred_T', 'SIFT4G_pred_.', 'SIFT4G_pred_D', 'SIFT4G_pred_T',
                 'Polyphen2_HDIV_pred_.', 'Polyphen2_HDIV_pred_B', 'Polyphen2_HDIV_pred_D', 'Polyphen2_HDIV_pred_P',
                 'Polyphen2_HVAR_pred_.', 'Polyphen2_HVAR_pred_B', 'Polyphen2_HVAR_pred_D', 'Polyphen2_HVAR_pred_P']


def readVCF(VCF_file):
    vcf_reader = vcf.Reader(filename=VCF_file)
    vcfINFO = []

    for record in vcf_reader:
        temp = []
        try:
            region = "".join(str(i) for i in record.INFO['Func.refGene'])
            if region != 'exonic':
                continue

            chrom = record.CHROM
            pos = record.POS
            variantID = record.ID
            ref = record.REF

            alt = "".join(str(i) for i in record.ALT)
            if alt == '.':
                alt = ref[0]
            sv_len = abs(len(ref) - len(alt))
            temp.append(chrom)
            temp.append(pos)
            if variantID:
                temp.append(variantID)
            else:
                temp.append('None')
            temp.append(ref)
            temp.append(alt)
            temp.append(sv_len)

            gene_list = []
            gene = "".join(str(i) for i in record.INFO['Gene.refGene']).split('\\x3b')
            for g in gene:
                gene_list.append(g)

            if not gene_list:
                try:
                    GENEINFO = record.INFO['GENEINFO']
                    for i in range(len(GENEINFO)):
                        gene_list.append(GENEINFO[i].split(':')[0])
                except:
                    pass

            if not gene_list:
                try:
                    AAChange_list = record.INFO['AAChange.refGene']
                    for i in range(len(AAChange_list)):
                        gene_list.append(AAChange_list[i].split(':')[0])
                except:
                    pass
            temp.append(gene_list)
            ExonicFunc = "".join(str(i) for i in record.INFO['ExonicFunc.refGene'])
            temp.append(ExonicFunc)

            if len(ref) >= 2:
                vcfINFO.append(temp)
        except:
            pass
    return vcfINFO


def get_sv_score_file_list(chrom, sv_start, sv_end, step, ref_length, file_dir, which_score):
    bed_start = sv_start // step * step
    sv_file_list = []
    while bed_start <= sv_end:
        end = bed_start + step
        if end > ref_length:
            end = ref_length
        if not str(chrom).startswith('chr'):
            chrom = 'chr' + chrom
        file_name = os.path.join(file_dir, which_score + '_' + chrom + '-' + str(bed_start + 1) + '-' + str(end) + '.txt')
        sv_file_list.append(file_name)
        bed_start += step
    return sv_file_list


def cadd_score(chrom, pos, sv_len, ref_length, snpINFO, CADD_dir, step):
    cadd_files_list = get_sv_score_file_list(chrom, pos, pos + sv_len, step, ref_length, CADD_dir, 'cadd')
    scope4currentsv = []
    for cadd_file in cadd_files_list:
        with open(cadd_file, 'r') as read_cadd:
            isFinished = False
            while not isFinished:
                line = read_cadd.readline()
                if line:
                    list_text = line.split()
                    temp = []
                    if list_text[0] == str(chrom) and int(list_text[1]) >= pos and int(list_text[2]) < pos + sv_len + 1:
                        temp.append(int(list_text[1]))  # pos
                        temp.append(str(list_text[3]))  # ref
                        temp.append(str(list_text[4]))  # alt
                        temp.append(float(list_text[5]))  # raw
                        temp.append(float(list_text[6]))  # phred
                        scope4currentsv.append(temp)
                    elif list_text[0] == str(chrom) and int(list_text[2]) >= pos + sv_len + 1:
                        isFinished = True
                else:
                    isFinished = True
    count = 0
    raw = 0
    phred = 0
    cadd_raw = 0
    cadd_phred = 0
    curses = 0
    if scope4currentsv:
        for i in range(len(snpINFO)):
            for j in range(curses, len(scope4currentsv), 1):
                if scope4currentsv[j][0] > snpINFO[i][0]:
                    break
                elif snpINFO[i][0] == scope4currentsv[j][0] and snpINFO[i][1] == scope4currentsv[j][1] and snpINFO[i][2] == scope4currentsv[j][2]:
                    count += 1
                    raw += scope4currentsv[j][3]
                    phred += scope4currentsv[j][4]
                    curses = j
                    break
    if count == 0 and len(scope4currentsv) > 0:
        cadd_raw = np.mean([x[3] for x in scope4currentsv])
        cadd_phred = np.mean([x[4] for x in scope4currentsv])
    else:
        cadd_raw = round(raw / count, 6)
        cadd_phred = round(phred / count, 6)
    return [cadd_raw, cadd_phred]


def Lof_Score(gene_list, gene_dir):
    lof_file = os.path.join(gene_dir, 'Lof.txt')
    count = 0
    lof_score = 0
    with open(lof_file, 'r') as read_lof:
        while True:
            line = read_lof.readline()
            if line:
                list_text = line.split()
                if str(list_text[0]) in gene_list:
                    lof_score += float(list_text[1])
                    count += 1
                if count == len(gene_list):
                    break
            else:
                break
    if count:
        lof_score = lof_score / count
    return lof_score


def str2float(str):
    if str == '.':
        return 0.0
    else:
        return float(str)

def minDT(c1, c2):
    if c1 == 'D' or c2 =='D':
        return 'D'
    elif c1 == '.' and c2 == '.':
        return '.'
    else:
        return 'T'

def maxDPB(c1, c2):
    if c1 == 'D' or c2 == 'D':
        return 'D'
    elif c1 == 'P' or c2 == 'P':
        return 'P'
    elif c1 == '.' and c2 == '.':
        return '.'
    else:
        return 'B'

def dbnsfp(chrom, pos, sv_len, ref_length, snpINFO, dbnsfp_dir, step):
    scope4currentsv = []
    dbnsfp_files_list = get_sv_score_file_list(chrom, pos, pos + sv_len, step, ref_length, dbnsfp_dir, 'dbnsfp')
    for dbnsfp_file in dbnsfp_files_list:
        with open(dbnsfp_file, 'r') as read_dbnsfp:
            while True:
                line = read_dbnsfp.readline()
                if line:
                    list_text = line.split()
                    temp = []
                    if str(list_text[0]) == str(chrom) and int(list_text[1]) >= pos and int(
                            list_text[2]) < pos + sv_len + 1:
                        temp.append(int(list_text[1]))  # 0  pos
                        temp.append(list_text[3])  # 1  ref
                        temp.append(list_text[4])  # 2  alt

                        temp.append(list_text[6])  # 3  SIFT_pred
                        temp.append(list_text[7])  # 4  SIFT4G_pred
                        temp.append(list_text[8])  # 5  Polyphen2_HDIV_pred
                        temp.append(list_text[9])  # 6  Polyphen2_HVAR_pred

                        temp.append(str2float(list_text[5]))  # 7  DamagePredCount
                        temp.append(str2float(list_text[15]))  # 8  VEST4
                        temp.append(str2float(list_text[21]))  # 9  MVP
                        temp.append(str2float(list_text[22]))  # 10  MPC
                        temp.append(str2float(list_text[31]))  # 11  DANN_score
                        temp.append(str2float(list_text[38]))  # 12  GenoCanyon_score
                        temp.append(str2float(list_text[39]))  # 13  integrated_fitCons_score
                        temp.append(str2float(list_text[40]))  # 14  GM12878_fitCons_score
                        temp.append(str2float(list_text[41]))  # 15  H1-hESC_fitCons_score
                        temp.append(str2float(list_text[42]))  # 16 HUVEC_fitCons_score
                        temp.append(str2float(list_text[43]))  # 17 LINSIGHT
                        temp.append(str2float(list_text[44]))  # 18 GERP++_NR
                        temp.append(str2float(list_text[45]))  # 19 GERP++_RS
                        temp.append(str2float(list_text[46]))  # 20 phyloP100way_vertebrate
                        temp.append(str2float(list_text[47]))  # 21 phyloP30way_mammalian
                        temp.append(str2float(list_text[48]))  # 22 phyloP17way_primate
                        temp.append(str2float(list_text[49]))  # 23 phastCons100way_vertebrate
                        temp.append(str2float(list_text[50]))  # 24 phastCons30way_mammalian
                        temp.append(str2float(list_text[51]))  # 25 phastCons17way_primate

                        scope4currentsv.append(temp)
                    elif list_text[0] == str(chrom) and int(list_text[2]) >= pos + sv_len + 1:
                        break
                else:
                    break
    count = 0
    DamagePredCount = 0
    SIFT_pred = '.'
    SIFT4G_pred = '.'
    Polyphen2_HDIV_pred = '.'
    Polyphen2_HVAR_pred = '.'
    if scope4currentsv:
        SIFT_pred = scope4currentsv[0][3]
        SIFT4G_pred = scope4currentsv[0][4]
        Polyphen2_HDIV_pred = scope4currentsv[0][5]
        Polyphen2_HVAR_pred = scope4currentsv[0][6]
    vest_score = 0
    mvp_score = 0
    mpc_score = 0
    DANN_score = 0
    GenoCanyon_score = 0
    integrated_fitCons_score = 0
    GM12878_fitCons_score = 0
    H1_hESC_fitCons_score = 0
    HUVEC_fitCons_score = 0
    LINSIGHT = 0
    GERP_NR = 0
    GERP_RS = 0
    phyloP100way_vertebrate = 0
    phyloP30way_mammalian = 0
    phyloP17way_primate = 0
    phastCons100way_vertebrate = 0
    phastCons30way_mammalian = 0
    phastCons17way_primate = 0
    curses = 0
    if len(scope4currentsv):
        for i in range(len(snpINFO)):
            if snpINFO[i][1] != snpINFO[i][2]:
                for j in range(curses, len(scope4currentsv), 1):
                    if scope4currentsv[j][0] > snpINFO[i][0]:
                        break
                    if snpINFO[i][0] == scope4currentsv[j][0] and snpINFO[i][1] == scope4currentsv[j][1] and snpINFO[i][2] == scope4currentsv[j][2]:
                        count += 1
                        curses = j
                        if SIFT_pred != 'D':
                            SIFT_pred = minDT(SIFT_pred, scope4currentsv[j][3])
                        if SIFT4G_pred != 'D':
                            SIFT4G_pred = minDT(SIFT4G_pred, scope4currentsv[j][4])
                        if Polyphen2_HDIV_pred != 'P':
                            Polyphen2_HDIV_pred = maxDPB(Polyphen2_HDIV_pred, scope4currentsv[j][5])
                        if Polyphen2_HVAR_pred != 'P':
                            Polyphen2_HVAR_pred = maxDPB(Polyphen2_HVAR_pred, scope4currentsv[j][6])

                        DamagePredCount += scope4currentsv[j][7]
                        vest_score += scope4currentsv[j][8]
                        mvp_score += scope4currentsv[j][9]
                        mpc_score += scope4currentsv[j][10]
                        DANN_score += scope4currentsv[j][11]
                        GenoCanyon_score += scope4currentsv[j][12]
                        integrated_fitCons_score += scope4currentsv[j][13]
                        GM12878_fitCons_score += scope4currentsv[j][14]
                        H1_hESC_fitCons_score += scope4currentsv[j][15]
                        HUVEC_fitCons_score += scope4currentsv[j][16]
                        LINSIGHT += scope4currentsv[j][17]
                        GERP_NR += scope4currentsv[j][18]
                        GERP_RS += scope4currentsv[j][19]
                        phyloP100way_vertebrate += scope4currentsv[j][20]
                        phyloP30way_mammalian += scope4currentsv[j][21]
                        phyloP17way_primate += scope4currentsv[j][22]
                        phastCons100way_vertebrate += scope4currentsv[j][23]
                        phastCons30way_mammalian += scope4currentsv[j][24]
                        phastCons17way_primate += scope4currentsv[j][25]
    if count != 0:
        DamagePredCount_mean = round(DamagePredCount / count, 6)
        vest_score_mean = round(vest_score / count, 6)
        mvp_score_mean = round(mvp_score / count, 6)
        mpc_score_mean = round(mpc_score / count, 6)
        DANN_score_mean = round(DANN_score / count, 6)
        GenoCanyon_score_mean = round(GenoCanyon_score / count, 6)
        integrated_fitCons_score_mean = round(integrated_fitCons_score / count, 6)
        GM12878_fitCons_score_mean = round(GM12878_fitCons_score / count, 6)
        H1_hESC_fitCons_score_mean = round(H1_hESC_fitCons_score / count, 6)
        HUVEC_fitCons_score_mean = round(HUVEC_fitCons_score / count, 6)
        LINSIGHT_mean = round(LINSIGHT / count, 6)
        GERP_NR_mean = round(GERP_NR / count, 6)
        GERP_RS_mean = round(GERP_RS / count, 6)
        phyloP100way_vertebrate_mean = round(phyloP100way_vertebrate / count, 6)
        phyloP30way_mammalian_mean = round(phyloP30way_mammalian / count, 6)
        phyloP17way_primate_mean = round(phyloP17way_primate / count, 6)
        phastCons100way_vertebrate_mean = round(phastCons100way_vertebrate / count, 6)
        phastCons30way_mammalian_mean = round(phastCons30way_mammalian / count, 6)
        phastCons17way_primate_mean = round(phastCons17way_primate / count, 6)
        return [SIFT_pred, SIFT4G_pred, Polyphen2_HDIV_pred, Polyphen2_HVAR_pred, DamagePredCount_mean, vest_score_mean,
                mvp_score_mean, mpc_score_mean, DANN_score_mean,
                GenoCanyon_score_mean, integrated_fitCons_score_mean, GM12878_fitCons_score_mean,
                H1_hESC_fitCons_score_mean, HUVEC_fitCons_score_mean, LINSIGHT_mean, GERP_NR_mean, GERP_RS_mean,
                phyloP100way_vertebrate_mean, phyloP30way_mammalian_mean, phyloP17way_primate_mean,
                phastCons100way_vertebrate_mean, phastCons30way_mammalian_mean, phastCons17way_primate_mean]
    elif scope4currentsv:
        avg = [SIFT_pred, SIFT4G_pred, Polyphen2_HDIV_pred, Polyphen2_HVAR_pred]
        for column in range(7, len(scope4currentsv[0])):
            currentCol = [x[column] for x in scope4currentsv]
            avg.append(round(np.mean(currentCol), 6))
        return avg
    else:
        return ['T', 'T', 'B', 'B', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

def maxGDI(current, newGDI):
    if current == 'High' or newGDI == 'High':
        return 'High'
    elif current == 'Low' and newGDI == 'Low':
        return 'Low'
    else:
        return 'Medium'

def GDI(gene_list, gene_dir):
    gdi_file = os.path.join(gene_dir, 'GDI.txt')
    sum_gdi = 0
    sum_gdi_phred = 0
    all_disease = 'Low'
    all_Mendelian = 'Low'
    mendelian_AD = 'Low'
    mendelian_AR = 'Low'
    all_PID = 'Low'
    pid_AD = 'Low'
    pid_AR = 'Low'
    all_cancer = 'Low'
    cancer_recessive = 'Low'
    cancer_dominant = 'Low'
    count = 0
    with open(gdi_file, 'r') as read_GDI:
        while True:
            line = read_GDI.readline()
            if line:
                list_text = line.split()
                if str(list_text[0]) in gene_list:
                    sum_gdi += float(list_text[1])
                    sum_gdi_phred += float(list_text[2])
                    all_disease = maxGDI(all_disease, list_text[3])
                    all_Mendelian = maxGDI(all_Mendelian, list_text[4])
                    mendelian_AD = maxGDI(mendelian_AD, list_text[5])
                    mendelian_AR = maxGDI(mendelian_AR, list_text[6])
                    all_PID = maxGDI(all_PID, list_text[7])
                    pid_AD = maxGDI(pid_AD, list_text[8])
                    pid_AR = maxGDI(pid_AR, list_text[9])
                    all_cancer = maxGDI(all_cancer, list_text[10])
                    cancer_recessive = maxGDI(cancer_recessive, list_text[11])
                    cancer_dominant = maxGDI(cancer_dominant, list_text[12])
                    count += 1
                if count == len(gene_list):
                    break
            else:
                break
    if count:
        sum_gdi = sum_gdi/count
        sum_gdi_phred = sum_gdi_phred/count
    return [sum_gdi, sum_gdi_phred, all_disease, all_Mendelian, mendelian_AD, mendelian_AR, all_PID, pid_AD, pid_AR, all_cancer, cancer_recessive, cancer_dominant]




def mcap(chrom, pos, sv_len, ref_length, snpINFO, mcap_dir, step):
    scope4currentsv = []
    mcap_files_list = get_sv_score_file_list(chrom, pos, pos + sv_len, step, ref_length, mcap_dir, 'mcap')
    for mcap_file in mcap_files_list:
        with open(mcap_file, 'r') as read_mcap:
            while True:
                line = read_mcap.readline()
                if line:
                    list_text = line.split()
                    temp = []
                    if list_text[0] == str(chrom) and int(list_text[1]) >= pos and int(list_text[2]) < pos + sv_len + 1:
                        temp.append(int(list_text[1]))  # pos
                        temp.append(str(list_text[3]))  # ref
                        temp.append(str(list_text[4]))  # alt
                        temp.append(float(list_text[5]))  # mcap
                        scope4currentsv.append(temp)
                    elif list_text[0] == str(chrom) and int(list_text[2]) >= pos + sv_len + 1:
                        break
                else:
                    break
    count = 0
    mcap_score = 0
    mcap_score_mean = 0
    curses = 0
    if len(scope4currentsv) > 0:
        for i in range(len(snpINFO)):
            if snpINFO[i][1] != snpINFO[i][2]:
                for j in range(curses, len(scope4currentsv), 1):
                    if scope4currentsv[j][0] > snpINFO[i][0]:
                        break
                    if snpINFO[i][0] == scope4currentsv[j][0] and snpINFO[i][1] == scope4currentsv[j][1] and snpINFO[i][2] == scope4currentsv[j][2]:
                        count += 1
                        mcap_score += scope4currentsv[j][3]
                        curses = j
    if count != 0:
        mcap_score_mean = round(mcap_score / count, 6)
    else:
        tempmcapscoreList = [x[3] for x in scope4currentsv]
        if tempmcapscoreList:
            mcap_score_mean = np.mean(tempmcapscoreList)
        else:
            mcap_score_mean = 0
    return mcap_score_mean


def RVIS_ExAC_4KW(gene_list, gene_dir):
    rvis_file = os.path.join(gene_dir, 'RVIS.txt')
    count = 0
    rvis1 = 0
    rvis2 = 0
    with open(rvis_file, 'r') as read_RVIS:
        while True:
            line = read_RVIS.readline()
            if line:
                list_text = line.split()
                if list_text[0] in gene_list or list_text[1] in gene_list or list_text[2] in gene_list or list_text[3] in gene_list or \
                        list_text[4] in gene_list:
                    rvis1 += float(list_text[5])
                    rvis2 += float(list_text[6])
                    count += 1
                if count == len(gene_list):
                    break
            else:
                break
    if count:
        rvis1 = rvis1 / count
        rvis2 = rvis2 / count
    return [rvis1, rvis2]


def revel(chrom, pos, sv_len, ref_length, snpINFO, revel_dir, step):
    scope4currentsv = []
    revel_files_list = get_sv_score_file_list(chrom, pos, pos + sv_len, step, ref_length, revel_dir, 'revel')
    for revel_file in revel_files_list:
        with open(revel_file, 'r') as read_revel:
            while True:
                line = read_revel.readline()
                if line:
                    list_text = line.split()
                    temp = []
                    if list_text[0] == str(chrom) and int(list_text[1]) >= pos and int(list_text[2]) < pos + sv_len + 1:
                        temp.append(int(list_text[1]))  # pos
                        temp.append(str(list_text[3]))  # ref
                        temp.append(str(list_text[4]))  # alt
                        temp.append(float(list_text[5]))  # revel
                        scope4currentsv.append(temp)
                    elif list_text[0] == str(chrom) and int(list_text[2]) >= pos + sv_len + 1:
                        break
                else:
                    break
    count = 0
    revel_score = 0
    revel_score_mean = 0
    curses = 0
    if len(scope4currentsv) > 0:
        for i in range(len(snpINFO)):
            if snpINFO[i][1] != snpINFO[i][2]:
                for j in range(curses, len(scope4currentsv), 1):
                    if scope4currentsv[j][0] > snpINFO[i][0]:
                        break
                    if snpINFO[i][0] == scope4currentsv[j][0] and snpINFO[i][1] == scope4currentsv[j][1] and snpINFO[i][2] == scope4currentsv[j][2]:
                        count += 1
                        revel_score += scope4currentsv[j][-1]
                        curses = j
        if count != 0:
            revel_score_mean = round(revel_score / count, 6)
        else:
            temprevelscoreList = [x[-1] for x in scope4currentsv]
            if temprevelscoreList:
                revel_score_mean = np.mean(temprevelscoreList)
            else:
                revel_score_mean = 0
    return revel_score_mean


def featuresInbwFiles(bwFiles, chrom, start, end):
    bwValueList = []
    for bwfile in bwFiles:
        bw = pyBigWig.open(bwfile)
        if 'chr' not in chrom:
            chrom = 'chr' + chrom
        avg = []
        try:
            avg = bw.stats(chrom, start, end)
        except:
            avg = [0]
        bwValue = avg[0]
        if not bwValue:
            bwValue = 0
        bwValueList.append(bwValue)
        bw.close()
    return bwValueList


def featuresForOverlap(bedFiles, chrom, start, end):
    overlapValueList = []
    if 'chr' not in str(chrom):
        chrom = 'chr' + chrom
    if not os.path.exists('tmp'):
        os.makedirs('tmp')

    # Prevent multi-process conflicts
    randomA = random.randint(1, 1000)
    randomB = random.randint(1000, 2000)

    tmpFileName = 'tmp/' + chrom + '_' + str(start) + '_' + str(end) + str(randomA) + '_' + str(randomB)
    with open(tmpFileName, 'w+') as writeTmp:
        writeTmp.write(chrom + '\t' + str(start) + '\t' + str(end))
    for bedfile in bedFiles:
        intersectionFileName = tmpFileName + '_' + bedfile.split('/')[-1] + '_ints.tsv'
        bedtoolsCommand = 'bedtools intersect -a ' + tmpFileName + ' -b ' + bedfile + ' -wao > ' + intersectionFileName
        os.system(bedtoolsCommand)
        intersectDict = []
        with open(intersectionFileName, 'r') as intersect:
            d = {}
            d['filename'] = bedfile
            for line in intersect:
                lineList = line.split()
                overlapResult = lineList[-1]
                if (chrom, start, end) in d:
                    d[(chrom, start, end)] += int(overlapResult)
                else:
                    d[(chrom, start, end)] = int(overlapResult)
            avgResult = d[(chrom, start, end)] / (int(end) - int(start))
            overlapValueList.append(avgResult)
        deleteInterSectTmp = 'rm ' + intersectionFileName
        os.system(deleteInterSectTmp)
    deleteTmpBed = 'rm ' + tmpFileName
    os.system(deleteTmpBed)
    return overlapValueList

def onehot_GDI(gdi_field):
    gdi_onehot = [0, 0, 0]
    if gdi_field == 'High':
        gdi_onehot = [1, 0, 0]
    elif gdi_field == 'Low':
        gdi_onehot = [0, 1, 0]
    elif gdi_field == 'Medium':
        gdi_onehot = [0, 0, 1]
    return gdi_onehot


def onehot_sift(sift_value):
    SIFT_onehot = [0, 0, 0]
    if sift_value == '.':
        SIFT_onehot = [1, 0, 0]
    elif sift_value == 'D':
        SIFT_onehot = [0, 1, 0]
    elif sift_value == 'T':
        SIFT_onehot = [0, 0, 1]
    return SIFT_onehot


def onehot_polyphen(polyphen_value):
    polyphen = [0, 0, 0, 0]
    if polyphen_value == '.':
        polyphen = [1, 0, 0, 0]
    elif polyphen_value == 'B':
        polyphen = [0, 1, 0, 0]
    elif polyphen_value == 'D':
        polyphen = [0, 0, 1, 0]
    elif polyphen_value == 'P':
        polyphen = [0, 0, 0, 1]
    return polyphen

def onehot_exonicFunc(exonicFunc_refGene):
    ExonicFuncRefGene = [0, 0, 0, 0, 0, 1]
    if exonicFunc_refGene == 'frameshift_deletion':
        ExonicFuncRefGene = [1, 0, 0, 0, 0, 0]
    elif exonicFunc_refGene == 'nonframeshift_deletion':
        ExonicFuncRefGene = [0, 1, 0, 0, 0, 0]
    elif exonicFunc_refGene == 'startloss':
        ExonicFuncRefGene = [0, 0, 1, 0, 0, 0]
    elif exonicFunc_refGene == 'stopgain':
        ExonicFuncRefGene = [0, 0, 0, 1, 0, 0]
    elif exonicFunc_refGene == 'stoploss':
        ExonicFuncRefGene = [0, 0, 0, 0, 1, 0]
    return ExonicFuncRefGene

def extract_features(reference, cadd_dir, mcap_dir, dbnsfp_dir, revel_dir, gene_dir, bw_dir, overlap_dir, step, variationINFO):
    featuresValue4currentSV = []
    chrom = variationINFO[0]
    pos = variationINFO[1]
    variantID = variationINFO[2]
    ref = variationINFO[3]
    alt = variationINFO[4]
    sv_len = variationINFO[5]
    location = str(chrom) + ':' + str(pos) + '-' + str(int(pos) + sv_len)
    featuresValue4currentSV.append(location)
    featuresValue4currentSV.append(variantID)
    featuresValue4currentSV.append(sv_len)
    try:
        fasta = pysam.FastaFile(reference)

        tmpchrom = chrom
        if not str(chrom).startswith('chr'):
            tmpchrom = 'chr' + chrom

        snpALT = fasta.fetch(tmpchrom, pos + sv_len, pos + 2 * sv_len)
        snpALT = snpALT.upper()

        ref_length = fasta.get_reference_length(tmpchrom)

        snpINFO = []
        if len(ref) == len(snpALT) + 1:
            for i in range(len(snpALT)):
                snppos = pos + i + 1
                snpref = ref[i + 1]
                snpalt = snpALT[i]
                temp4snp = []
                if snpref != snpalt:
                    temp4snp.append(snppos)
                    temp4snp.append(snpref)
                    temp4snp.append(snpalt)
                    snpINFO.append(temp4snp)
            cadd_rawScore, cadd_phred = cadd_score(chrom, pos, sv_len, ref_length, snpINFO, cadd_dir, step)
            featuresValue4currentSV.append(cadd_rawScore)
            featuresValue4currentSV.append(cadd_phred)

            gene = variationINFO[6]
            GDI_lists = GDI(gene, gene_dir)

            GDI_score = GDI_lists[0]
            featuresValue4currentSV.append(GDI_score)
            GDI_Phred = GDI_lists[1]
            featuresValue4currentSV.append(GDI_Phred)

            All_disease_causing = onehot_GDI(GDI_lists[2])
            All_Mendelian = onehot_GDI(GDI_lists[3])
            Mendelian_AD = onehot_GDI(GDI_lists[4])
            Mendelian_AR = onehot_GDI(GDI_lists[5])
            all_PID = onehot_GDI(GDI_lists[6])
            PID_AD = onehot_GDI(GDI_lists[7])
            PID_AR = onehot_GDI(GDI_lists[8])
            All_cancer = onehot_GDI(GDI_lists[9])
            Cancer_recessive = onehot_GDI(GDI_lists[10])
            Cancer_dominant = onehot_GDI(GDI_lists[11])

            lof_score = Lof_Score(gene, gene_dir)
            featuresValue4currentSV.append(lof_score)

            exonicFunc_refGene = variationINFO[7]
            ExonicFuncRefGene = onehot_exonicFunc(exonicFunc_refGene)

            mcap_score = mcap(chrom, pos, sv_len, ref_length, snpINFO, mcap_dir, step)
            featuresValue4currentSV.append(round(mcap_score, 6))

            dbnsfp_scores = dbnsfp(chrom, pos, sv_len, ref_length, snpINFO, dbnsfp_dir, step)
            SIFT_pred = onehot_sift(dbnsfp_scores[0])
            SIFT4G_pred = onehot_sift(dbnsfp_scores[1])
            Polyphen2_HDIV = onehot_polyphen(dbnsfp_scores[2])
            Polyphen2_HVAR = onehot_polyphen(dbnsfp_scores[3])

            for i in range(4, len(dbnsfp_scores), 1):
                featuresValue4currentSV.append(dbnsfp_scores[i])

            RVIS_ExAC_scores = RVIS_ExAC_4KW(gene, gene_dir)
            for i in range(len(RVIS_ExAC_scores)):
                featuresValue4currentSV.append(RVIS_ExAC_scores[i])

            revel_score = revel(chrom, pos, sv_len, ref_length, snpINFO, revel_dir, step)
            featuresValue4currentSV.append(revel_score)

            bwFeatures_scores = featuresInbwFiles(bw_dir, chrom, pos, pos + sv_len)
            for i in range(len(bwFeatures_scores)):
                featuresValue4currentSV.append(round(bwFeatures_scores[i], 6))

            overlapFeatures = featuresForOverlap(overlap_dir, chrom, pos, pos + sv_len)
            for i in range(len(overlapFeatures)):
                featuresValue4currentSV.append(overlapFeatures[i])

            for i in range(len(All_disease_causing)):
                featuresValue4currentSV.append(All_disease_causing[i])

            for i in range(len(All_Mendelian)):
                featuresValue4currentSV.append(All_Mendelian[i])

            for i in range(len(Mendelian_AD)):
                featuresValue4currentSV.append(Mendelian_AD[i])

            for i in range(len(Mendelian_AR)):
                featuresValue4currentSV.append(Mendelian_AR[i])

            for i in range(len(all_PID)):
                featuresValue4currentSV.append(all_PID[i])

            for i in range(len(PID_AD)):
                featuresValue4currentSV.append(PID_AD[i])

            for i in range(len(PID_AR)):
                featuresValue4currentSV.append(PID_AR[i])

            for i in range(len(All_cancer)):
                featuresValue4currentSV.append(All_cancer[i])

            for i in range(len(Cancer_recessive)):
                featuresValue4currentSV.append(Cancer_recessive[i])

            for i in range(len(Cancer_dominant)):
                featuresValue4currentSV.append(Cancer_dominant[i])

            for i in range(len(ExonicFuncRefGene)):
                featuresValue4currentSV.append(ExonicFuncRefGene[i])

            for i in range(len(SIFT_pred)):
                featuresValue4currentSV.append(SIFT_pred[i])

            for i in range(len(SIFT4G_pred)):
                featuresValue4currentSV.append(SIFT4G_pred[i])

            for i in range(len(Polyphen2_HDIV)):
                featuresValue4currentSV.append(Polyphen2_HDIV[i])

            for i in range(len(Polyphen2_HVAR)):
                featuresValue4currentSV.append(Polyphen2_HVAR[i])
            print('chr' + str(chrom) + ', pos ' + str(pos) + ' deletion variation features calculation completed')
            return featuresValue4currentSV
    except:
        print('chr' + str(chrom) + '-' + str(pos) + ' failed.')


def output_result(result, output_file):
    locationAndId = [x[:3] for x in result if x]
    featureValues = [x[3:] for x in result if x]

    finalFeatures = []
    featureValues = preprocessing.normalize(np.array(featureValues), norm='l2')
    for i in range(len(locationAndId)):
        itemFeatures = []
        itemFeatures.append(locationAndId[i][0])
        itemFeatures.append(locationAndId[i][1])
        itemFeatures.append(locationAndId[i][2])
        for j in range(len(featureValues[i])):
            itemFeatures.append(featureValues[i][j])
        finalFeatures.append(itemFeatures)

    print('Writing feature values to ' + output_file + ' files...')
    f = open(output_file, 'w', encoding='utf-8')
    csv_writer = csv.writer(f)
    csv_writer.writerow(features_head)
    for line in finalFeatures:
        try:
            csv_writer.writerow(line)
        except:
            pass
    f.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate a feature matrix from deletion SV data')
    parser.add_argument('-i', '--input', required=True, help='A VCF file annotated by Annovar')
    parser.add_argument('-o', '--output', required=True, help='Output file in csv format')
    parser.add_argument('-R', '--reference', required=True, help='Reference genome file')
    parser.add_argument('-p', '--processes', type=int, default=8, help='Number of processes for parallel extraction of feature values')
    parser.add_argument('-g', '--gene', required=True, help='The directory for storing gene-based files, including GDI.txt, RVIS.txt and Lof.txt')
    parser.add_argument('-d', '--dbnsfp', required=True, help='The directory for storing dbnsfp features of each bed, the bed size is 1000000. The naming format of '
                                                              'a single file in the directory, such as chromosome 1, 1-1000000, is dbnsfp_chr1-1-1000000.txt')
    parser.add_argument('-m', '--mcap', required=True, help='The directory for storing M-CAP scores of each bed, the bed size is 1000000. The naming format of '
                                                            'a single file in the directory, such as chromosome 1, 1-1000000, is mcap_chr1-1-1000000.txt')
    parser.add_argument('-r', '--revel', required=True, help='The directory for storing REVEL scores of each  bed, the bed size is 1000000. The naming format of '
                                                             'a single file in the directory, such as chromosome 1, 1-1000000, is revel_chr1-1-1000000.txt')
    parser.add_argument('-b', '--bigwig', nargs='+', required=True, help='Feature files in BigWig format')
    parser.add_argument('-l', '--overlap', nargs='+', required=True, help='Gene coordinates for overlap features')
    parser.add_argument('-c', '--cadd', required=True, help='The directory for storing CADD scores of each  bed, the bed size is 1000000. The naming format of '
                                                            'a single file in the directory, such as chromosome 1, 1-1000000, is cadd_chr1-1-1000000.txt')
    print("Parsing arguments...")

    args = parser.parse_args()
    output_file = args.output
    reference = args.reference
    vcfFile = args.input
    num_process = args.processes
    cadd_dir = args.cadd
    mcap_dir = args.mcap
    dbnsfp_dir = args.dbnsfp
    revel_dir = args.revel
    gene_dir = args.gene
    bw_dir = args.bigwig
    overlap_dir = args.overlap

    time_start = time.time()
    print('Reading vcf file information...')
    vcfINFO = readVCF(vcfFile)
    print('The number of variants on exon is', len(vcfINFO))
    print('Calculating feature values from multiple files...')
    step = 1000000
    partial_func = partial(extract_features, reference, cadd_dir, mcap_dir, dbnsfp_dir, revel_dir, gene_dir, bw_dir, overlap_dir, step)
    pool = multiprocessing.Pool(processes=num_process)
    result = pool.map(partial_func, vcfINFO)
    pool.close()
    pool.join()

    output_result(result, output_file)

    print('All SV feature values calculations have been completed.')
    time_end = time.time()
    print('Total running time = ', round(time_end - time_start, 2), 'seconds')
