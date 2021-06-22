# source ./script/test_data_preparation.py from https://github.com/kanamekojima/rnnimp

from argparse import ArgumentParser
from io import TextIOWrapper
import json
import os
from posixpath import sep
import sys
import re
import os
from tqdm.notebook import tqdm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# from lib.utils import general as g
# from lib.config.config_class import vcf_config, mani_config
from ..utils import general as g
from ..config.config_class import vcf_config, mani_config, legend_config

def complement_base(base):
    if base == 'A':
        return 'T'
    if base == 'T':
        return 'A'
    if base == 'G':
        return 'C'
    if base == 'C':
        return 'G'
    return 'N'

def get_sequence(region, hg_refgenome):
    command = 'samtools faidx {:s} {:s}'.format(hg_refgenome, region)
    lines = g.system_with_stdout(command)
    sequence = ''.join(lines[1:])
    return sequence

def align(ref_seq, seq):
    for i in range(len(ref_seq) - len(seq)):
        mismatch_flag = False
        for j in range(len(seq)):
            if ref_seq[i + j] != seq[j]:
                mismatch_flag = True
                break
        if not mismatch_flag:
            return i
    return -1

def parse_manifest_from_vcfgenotype(vcf_file:str,chr_nums:list,hg_refgenome:str)->dict:
    marker_dict = {}
    with g.reading(vcf_file) as vcfp:
        header_dict = process_vcf_header(vcfp)
        for line in tqdm(vcfp, desc='create marker from manifest'):
            # đọc line và tách ra thành từng item
            chrom, position, snp, ref, alts, info_dict, genotypes = get_data_from_line(
                line, header_dict)
            if chrom not in chr_nums:
                continue
            if chrom not in marker_dict:
                marker_dict[chrom] = {}
            if position not in marker_dict[chrom]:
                marker_dict[chrom][position] = []
            marker_dict[chrom][position].append([ref]+alts)
    return marker_dict

def parse_manifest_from_manifest(manifest_file:str, chr_nums: list, hg_refgenome:str)->dict:
    marker_dict = {}
    with g.reading(manifest_file) as fp:
        for line in fp:
            if line.startswith(mani_config.start_header_end_line):
                break
        line = fp.readline()
        items = g.line_to_items(line, split=',')
        IlmnStrand_col = items.index(mani_config.IlmnStrand_col)
        Chr_col = items.index(mani_config.Chr_col)
        MapInfo_col = items.index(mani_config.MapInfo_col)
        SNP_col = items.index(mani_config.SNP_col)
        SourceStrand_col = items.index(mani_config.SourceStrand_col)
        SourceSeq_col = items.index(mani_config.SourceSeq_col)
        RefStrand_col = items.index(mani_config.RefStrand_col)
        for line in tqdm(fp, desc='create marker from manifest'):
            if line.startswith(mani_config.start_data_end_line):
                break
            items = g.line_to_items(line, mani_config.data_line_split_params)
            chrom = items[Chr_col]
            if chrom not in chr_nums:
                continue
            if chrom not in marker_dict:
                marker_dict[chrom] = {}

            position = items[MapInfo_col]
            alleles = items[SNP_col][1:-1].split(mani_config.alleles_split_params)
            if alleles[0] in mani_config.alleles_position_0:
                source_seq = items[SourceSeq_col]
                upstream_seq, a0, a1, downstream_seq = re.split(
                    mani_config.sourceSeq_split_params, source_seq)
                upstream_seq = upstream_seq.upper()
                downstream_seq = downstream_seq.upper()
                if a0 in mani_config.alleles_missing:
                    a0 = mani_config.alleles_character_replace_missing
                if a1 in mani_config.alleles_missing:
                    a1 = mani_config.alleles_character_replace_missing
                a0_seq = upstream_seq + a0 + downstream_seq
                a1_seq = upstream_seq + a1 + downstream_seq
                margin = 10
                region_start = int(position) - len(upstream_seq) - margin
                region_end = region_start + margin \
                    + max(len(a0_seq), len(a1_seq)) + margin
                region = 'chr{:s}:{:d}-{:d}'.format(
                    chrom, region_start, region_end)
                ref_seq = get_sequence(region, hg_refgenome)
                a0_align = align(ref_seq, a0_seq)
                a1_align = align(ref_seq, a1_seq)
                indel_position = max(a0_align, a1_align) + len(upstream_seq) \
                    + region_start - 1
                alleles[0] = upstream_seq[-1] + a0
                alleles[1] = upstream_seq[-1] + a1
                position = str(indel_position)
            else:
                if items[RefStrand_col] in mani_config.alleles_missing:
                    alleles[0] = complement_base(alleles[0])
                    alleles[1] = complement_base(alleles[1])
            if position not in marker_dict[chrom]:
                marker_dict[chrom][position] = []
            marker_dict[chrom][position].append(alleles)
    return marker_dict

def parse_manifest(manifest_file:str, chr_nums: list, hg_refgenome:str)->dict:
    isvcf = any([manifest_file.endswith(tail) for tail in vcf_config.vcf_tails])
    if isvcf:
        parse_manifest_from_vcfgenotype(manifest_file,chr_nums,hg_refgenome)
    else:
        parse_manifest_from_manifest(manifest_file,chr_nums,hg_refgenome)

def mapping_marker(chrom, position, ref, alts:list, marker_dict):
    flag = legend_config.unobserve
    if chrom in marker_dict:
        maker_position = marker_dict[chrom]
        if position in maker_position:
            for alleles in maker_position[position]:
                if ref in alleles and np.any([alt in alleles for alt in alts]):
                    flag = legend_config.observe
                    break
                # if alleles[0] == ref and alleles[1] in alts:
                #     flag = legend_config.observe
                #     break
                # if alleles[1] in alts and alleles[0] == ref:
                #     flag = legend_config.observe
                #     break
    return flag

def get_data_from_line(line, header_dict):
    items = g.line_to_items(line, vcf_config.vcf_data_split_params)
    # Item đầu tiên là chrom
    chrom = items[header_dict[vcf_config.chrom_col]]
    # Item thứ 2 là position
    position = items[header_dict[vcf_config.position_col]]
    # Item thứ 3 là snp
    snp = items[header_dict[vcf_config.snp_col]]
    # Item thứ 4 là REF
    ref = items[header_dict[vcf_config.ref_col]]
    # Item thứ  5 là ALT
    alts = items[header_dict[vcf_config.alt_col]].split(',')
    # get info data
    info_dict = g.list_expression_to_dict(
        items[header_dict[vcf_config.info_col]].split(';'))
    # Dữ liệu về genotypes từ index 9 đổ về sau
    genotypes = items[vcf_config.vcf_header_sample_startindex:]
    return chrom, position, snp, ref, alts, info_dict, genotypes

def is_missing_genotype(genotype):
    if genotype == '.':
        return True
    if genotype == './.':
        return True
    if genotype == '.|.':
        return True
    return False

def vcf2haplegend(vcf_file, keep_sample_list, marker_dict, output_prefix):
    g.mkdir(os.path.dirname(output_prefix))
    # get file path
    hap_file = legend_config.get_hap_file_name(output_prefix)
    legend_file = legend_config.get_legend_file_name(output_prefix)
    with g.reading(vcf_file) as fp:
        header_dict = process_vcf_header(fp)
        samples = list(header_dict[vcf_config.header_dict_sample_key].keys())
        # sort sample by index
        samples.sort(
            key=lambda sample: header_dict[vcf_config.header_dict_sample_key][sample])
        if keep_sample_list is None:
            keep_sample_id_list = list(range(len(samples)))
        else:
            keep_sample_id_list = [
                samples.index(sample) for sample in keep_sample_list
                if sample in samples
            ]
        sample_file = legend_config.get_sample_file_name(output_prefix)
        with open(sample_file, 'wt') as s_fp:
            s_fp.write('{}\n'.format(legend_config.sample_header_line))
            s_fp.write('{}\n'.format(legend_config.sample_first_line))
            for sample_id in keep_sample_id_list:
                sample_dict = legend_config.get_sample_dict(samples[sample_id], samples[sample_id], '0')
                sample_value_line = legend_config.get_sample_value_line(sample_dict)
                s_fp.write('{}\n'.format(sample_value_line))
        sample_size = len(keep_sample_id_list)
        
        with g.writing(hap_file) as h_fp, \
                g.writing(legend_file) as l_fp:

            l_fp.write('{}\n'.format(legend_config.legend_header_line))
            allele_id_list = [None] * (2 * sample_size)
            # danh sách ALT ID
            hap_record = [None] * (2 * sample_size)
            # Haplotype record
            for line in tqdm(fp, desc='vcf to haplegend'):
                # đọc line và tách ra thành từng item
                chrom, position, snp, ref, alts, info_dict, genotypes = get_data_from_line(
                    line, header_dict)
                # Loop qua danh sách sample và lấy dữ liệu allele_id của nó
                for i, sample_id in enumerate(keep_sample_id_list):
                    genotype = genotypes[sample_id]
                    if is_missing_genotype(genotype):
                        allele_id_list[2 * i] = None
                        allele_id_list[2 * i + 1] = None
                    else:
                        allele_id_pair = map(int, genotype.split(vcf_config.genotype_split_params))
                    for j, allele_id in enumerate(allele_id_pair):
                        allele_id_list[2 * i + j] = allele_id
                nb_alt = len(alts)
                afs = info_dict[vcf_config.info_key.af]
                flag = mapping_marker(chrom, position, ref, alts, marker_dict)
                for i, alt in enumerate(alts):
                    # duyệt qua danh sách alt
                    if alt == '.' or alt.startswith('<'):
                        continue
                    alt_allele_id = i + 1
                    # duyệt qua danh sách allele của các sample
                    for j, allele_id in enumerate(allele_id_list):
                        if allele_id is None:
                            # nếu allete_id ko tồn tại
                            hap_record[j] = '?'
                        elif allele_id == alt_allele_id:
                            # giống với alt thì là 1
                            hap_record[j] = '1'
                        else:
                            # khác alt thì là 0
                            hap_record[j] = '0'
                    h_fp.write(legend_config.hap_split_params.join(hap_record))
                    h_fp.write('\n')
                    legend_dict = legend_config.get_legend_dict(snp, chrom, position, ref, alt, afs[i], flag, nb_alt)
                    legend_value_line = legend_config.get_legend_value_line(legend_dict)
                    l_fp.write('{}\n'.format(legend_value_line))
    return hap_file, legend_file

def prepare_test_hap(hap_file, legend_file, output_prefix):
    array_marker_flag_list = []
    output_legend_file = legend_config.get_legend_file_name(output_prefix)
    with g.reading(legend_file) as fp, \
            g.writing(output_legend_file) as w_fp:
        line = fp.readline()
        items = g.line_to_items(line)
        try:
            array_marker_flag_col = items.index(
                legend_config.array_marker_flag)
        except ValueError:
            g.print_error()
            g.print_error('Error: Header {} not found in '.format(legend_config.array_marker_flag)
                          + legend_file)
            g.print_error()
            sys.exit(0)
        w_fp.write(line)
        for line in tqdm(fp, desc='prepare marker to legend'):
            items = g.line_to_items(line)
            array_marker_flag = items[array_marker_flag_col] == '1'
            array_marker_flag_list.append(array_marker_flag)
            if array_marker_flag:
                w_fp.write(line)
    output_hap_file = legend_config.get_hap_file_name(output_prefix)
    with g.reading(hap_file) as fp, \
            g.writing(output_hap_file) as w_fp:
        for i, line in enumerate(tqdm(fp, desc='prepare marker to hap')):
            if array_marker_flag_list[i]:
                w_fp.write(line)

# TODO move to genhelper
def process_vcf_header(rf: TextIOWrapper, wf: TextIOWrapper = None, write_callback=None):
    '''
    input:
        rf (TextIOWrapper): read file
        wf (TextIOWrapper): write file, default is None
        write_callback (str,func): define what you wanna write in file. Default is None:
            if func, it will get param are current line and previous line write_callback(line,old_line)
    '''
    header_dict = None
    sample_dict = None
    old_line = None
    for line in tqdm(rf, desc='process header vcf'):
        # reading region
        if line.startswith(vcf_config.vcf_header_line_startswith):
            items = g.line_to_items(line, vcf_config.vcf_data_split_params)
            header_dict = g.list_to_dict(
                items[:vcf_config.vcf_header_sample_startindex])
            sample_dict = g.list_to_dict(
                items[vcf_config.vcf_header_sample_startindex:], vcf_config.vcf_header_sample_startindex)
            header_dict[vcf_config.header_dict_sample_key] = sample_dict
            break
        # writing region
        if wf is None:
            continue
        w_line = line
        if type(write_callback) is type(lambda x: x):
            w_line = write_callback(line, old_line)
        elif type(write_callback) is str:
            w_line = write_callback
        wf.write(w_line)
        old_line = line
    return header_dict

def get_inter_bin_width(df:pd.DataFrame,inter_bin_width_percen:float):
    return int(df.shape[0]*inter_bin_width_percen)

def get_head_observe(old_df: pd.DataFrame,inter_bin_width_percen:float):
    return old_df.tail(get_inter_bin_width(old_df,inter_bin_width_percen))

def get_nb_bins_each_nb_variants(variant_labels):
    '''
    return:
        return dataframe contain data number bins with each number variants
    '''
    labels, nb_values = np.unique(variant_labels,return_counts=True)
    values, nb_bins = np.unique(nb_values,return_counts=True)
    return pd.DataFrame({'nb variants':values,'nb bin':nb_bins}), labels

# TODO move to general
def alert_yes_no(header,contain,question):
    # return False
    print('{}\n'.format(header))
    print(contain)
    print('{} [y/n]'.format(question))
    check = input()
    return check.startswith('n') or check.startswith('N')

def draw_hist_observe(df:pd.DataFrame,bin_col:str):
    df_observe = df[df[legend_config.array_marker_flag] == int(legend_config.observe)][bin_col].value_counts()
    nb_obser_col = 'nb_observe'
    df_observe = pd.DataFrame({'bin':df_observe.index,nb_obser_col:df_observe.values})
    print(df_observe[nb_obser_col].describe(percentiles=[0.1,0.25,0.5,0.75,0.9]))
    print('\n')

def make_bin_legend_data(legend_data:pd.DataFrame,number_bin,bin_observe = False):
    nb_data = len(legend_data)
    bin_col = 'bin'
    if bin_observe:
        bin_col_data = pd.DataFrame({bin_col:np.full(nb_data,np.NaN)})
        observes = legend_data[legend_config.array_marker_flag] == int(legend_config.observe)
        nb_observe = np.count_nonzero(observes)
        bin_col_data[bin_col].values[observes] = pd.qcut(np.arange(nb_observe),number_bin,labels=False)
        bin_col_data[bin_col].fillna(method='pad',inplace=True)
        bin_col_data[bin_col].fillna(method='bfill',inplace=True)
        bin_col_data[bin_col].astype(int)
        return bin_col_data[bin_col].values
    else:
        return pd.qcut(np.arange(nb_data),number_bin,labels=False)

def process_legend_to_region(legend_file:str,hap_file:str,number_bin:int,inter_bin_width_percen:float = 0,no_observation=False,output_folder=None,bin_observe=False)->str:
    #make region folder
    legend_region_folder = ''
    hap_region_folder = ''
    if output_folder is not None:
        os.makedirs(output_folder,exist_ok=True)
        if not os.path.isdir(output_folder):
            output_folder, _ = g.get_dir_and_base_name(output_folder)
        legend_region_folder = legend_config.make_region_dir(output_folder)
        hap_region_folder = legend_region_folder
    else:
        legend_region_folder = legend_config.make_region_dir(legend_file)
        hap_region_folder = legend_config.make_region_dir(hap_file)

    legend_data = pd.read_csv(legend_file,sep=legend_config.legend_split_params)
    bin_col = 'bin'
    legend_data[bin_col] =  make_bin_legend_data(legend_data,number_bin,bin_observe)
    draw_hist_observe(legend_data,bin_col)
    contain_alert, labels = get_nb_bins_each_nb_variants(legend_data[bin_col].values)
    header_alert = 'Info about number variant and number bin each variant'
    question_alert = 'Agree to div bin by this?'    
    while alert_yes_no(header_alert,contain_alert,question_alert):
        print('Chose new number bin')
        old_nb_bin = number_bin
        try:
            number_bin = int(input())
        except:
            number_bin = old_nb_bin
            print('number bin contain number character only!\n')
            continue
        legend_data[bin_col] = make_bin_legend_data(legend_data,number_bin,bin_observe)
        draw_hist_observe(legend_data,bin_col)
        contain_alert, labels = get_nb_bins_each_nb_variants(legend_data[bin_col].values)
    with g.reading(hap_file) as hap_file_stream:
        # label are sorted and orderly in data.
        # When we loop over label each line (exclude header) at legend file is map with line at hap file.
        max_label = max(labels)
        nb_character = len(str(max_label))
        legend_files = []
        legend_gtrue_files = []
        hap_files = []
        hap_gtrue_files = []
        next_header_legend = pd.DataFrame()
        next_header_hap = []
        for i, label in enumerate(tqdm(labels,desc='create region data from legend')):
            # get new data at this loop
            bin_data = legend_data[legend_data[bin_col] == label] # this bin data
            nb_data = len(bin_data)
            this_privious_index = int(nb_data*inter_bin_width_percen)
            this_next_index = nb_data - int(nb_data*inter_bin_width_percen)
            # file_path
            legend_region_file = legend_config.get_legend_region_file_name(legend_file,label,nb_character,output_folder=legend_region_folder)
            legend_region_gtrue_file = legend_config.get_legend_gtrue_file(legend_region_file)
            hap_region_file = legend_config.get_hap_region_file_name(hap_file,label,nb_character,output_folder=hap_region_folder)
            hap_region_gtrue_file = legend_config.get_hap_gtrue_file(hap_region_file)
            # write hap file
            previous_tail_hap = []            
            with g.writing(hap_region_file) as hap_region_file_stream, \
                g.writing(hap_region_gtrue_file) as hap_region_gtrue_file_stream:
                
                # writting concat data of 2 region
                if len(next_header_hap) > 0:
                    hap_region_file_stream.writelines(next_header_hap)
                    if no_observation == False:
                        hap_region_gtrue_file_stream.writelines(next_header_hap)
                    next_header_hap.clear()                
                
                for i in tqdm(range(nb_data),desc='process region {}'.format(label),leave=False):
                    line = hap_file_stream.readline()
                    # Nếu flag là observe thì ghi vô cả 2 file
                    flag = str(bin_data.iloc[i][legend_config.array_marker_flag])
                    if flag == legend_config.observe:
                        hap_region_file_stream.write(line)
                        if i >= this_next_index:
                            next_header_hap.append(line)
                        if i < this_privious_index:
                            previous_tail_hap.append(line)                            
                        if no_observation:
                            continue                    
                    hap_region_gtrue_file_stream.write(line)
            # write data to previous hap file
            if len(hap_files)>0:
                with g.appending(hap_files[-1]) as previous_hap_file:
                    previous_hap_file.writelines(previous_tail_hap)
            if len(hap_gtrue_files)>0 and no_observation == False:
                with g.appending(hap_gtrue_files[-1]) as previous_hap_gtrue_file:
                    previous_hap_gtrue_file.writelines(previous_tail_hap)
            previous_tail_hap.clear()
            # Write legend file
            #Previous legend
            previous_tail_legend = pd.DataFrame(bin_data.iloc[:this_privious_index])
            previous_tail_legend = previous_tail_legend[previous_tail_legend[legend_config.array_marker_flag] == int(legend_config.observe)]
            if len(previous_tail_legend) > 0:
                if len(legend_files) > 0:
                    previous_legend = pd.read_csv(legend_files[-1],sep=legend_config.legend_split_params)
                    previous_legend = pd.concat([previous_legend,previous_tail_legend])
                    legend_config.legend_dataframe_to_csv(legend_files[-1],previous_legend)
                    del previous_legend
                
                if len(legend_gtrue_files) > 0 and no_observation == False:
                    previous_legend_gtrue = pd.read_csv(legend_gtrue_files[-1],sep=legend_config.legend_split_params)
                    previous_legend_gtrue = pd.concat([previous_legend_gtrue,previous_tail_legend])
                    legend_config.legend_dataframe_to_csv(legend_gtrue_files[-1],previous_legend_gtrue)
                    del previous_legend_gtrue            
            del previous_tail_legend
            # make data for next region
            temp_next_header_legend = pd.DataFrame(bin_data.iloc[this_next_index:])
            temp_next_header_legend = temp_next_header_legend[temp_next_header_legend[legend_config.array_marker_flag]==int(legend_config.observe)]
            #This legend
            bin_data = pd.concat([next_header_legend,bin_data])
            next_header_legend = temp_next_header_legend
            del temp_next_header_legend
            # save input file
            this_legend = bin_data[bin_data[legend_config.array_marker_flag] == int(legend_config.observe)]
            legend_config.legend_dataframe_to_csv(legend_region_file,this_legend)
            del this_legend
            # save gtrue
            this_legend_gtrue = bin_data
            del bin_data
            if no_observation:
                this_legend_gtrue = this_legend_gtrue[this_legend_gtrue[legend_config.array_marker_flag] == int(legend_config.unobserve)]
            legend_config.legend_dataframe_to_csv(legend_region_gtrue_file,this_legend_gtrue)
            del this_legend_gtrue
            # add file
            hap_files.append(hap_region_file)
            hap_gtrue_files.append(hap_region_gtrue_file)
            legend_files.append(legend_region_file)
            legend_gtrue_files.append(legend_region_gtrue_file)
    return legend_region_folder

def find_fw_index(position,positions):
    bigger_indexs = np.where(positions >= position)[0]
    index = None
    if len(bigger_indexs) > 0:
        index = np.min(bigger_indexs)        
    else:
        index = len(positions)
    index = index -1
    if index == -1:
        index = None
    return index

def find_bw_index(position,positions):
    bigger_indexs = np.where(positions > position)[0]
    index = None
    if len(bigger_indexs) > 0:
        index = np.min(bigger_indexs)
    return index

def process_region_config(legend_folder:str,default_config_path:str):
    # get legend path, not gtrue path
    paths = []
    for file_name in os.listdir(legend_folder):
        path = os.path.join(legend_folder,file_name)
        if os.path.isfile(path) and path.endswith(legend_config.legend_tail) and not path.endswith(legend_config.gtrue_legend_tail):
            paths.append(path)
    assert np.all([os.path.isfile(legend_config.get_legend_gtrue_file(path)) for path in paths]), 'gtrue file must in same folder'
    paths.sort()
    for path in tqdm(paths,desc='make region config'):
        indata = pd.read_csv(path,sep=legend_config.legend_split_params)
        gtrue_path = legend_config.get_legend_gtrue_file(path)
        gtdata = pd.read_csv(gtrue_path,sep=legend_config.legend_split_params)        
        ids = indata[legend_config.position].values        
        fw_indexs = list(gtdata[legend_config.position].apply(lambda x: find_fw_index(x,ids)))
        bw_indexs = list(gtdata[legend_config.position].apply(lambda x: find_bw_index(x,ids)))
        num_inputs = len(indata)
        legend_config.to_region_config(num_inputs,fw_indexs,bw_indexs,path,default_config_path)
    
def process_genotyping_vcf(vcf_file, manifest_file, hg_refgenome, output_prefix, chroms):

    g.check_required_software('samtools')
    g.check_required_file(vcf_file)
    g.check_required_file(manifest_file)
    g.check_required_file(hg_refgenome)

    marker_dict = parse_manifest(manifest_file, chroms, hg_refgenome)
    with g.reading(vcf_file) as rf, \
            g.writing(output_prefix+'.vcf.gz') as wf:
        header_dict = process_vcf_header(rf, wf)
        for line in tqdm(rf, desc='mapping vcf with manifest file'):
            chrom, position, _, ref, alts, _, _ = get_data_from_line(
                line, header_dict)
            flag = mapping_marker(chrom, position, ref, alts, marker_dict)
            if flag == legend_config.observe:
                wf.write(line)
    print('\ngenotyping from vcf done!\n')

def process_data_to_legend(vcf_file, manifest_file, hg_refgenome, chroms, output_prefix, sample_list_file=None):

    # g.check_required_software('samtools')
    g.check_required_file(vcf_file)
    g.check_required_file(manifest_file)
    g.check_required_file(hg_refgenome)

    true_output_prefix = g.get_true_file_prefix(output_prefix)

    keep_sample_list = []
    try:
        with open(sample_list_file) as fp:
            for line in fp:
                keep_sample_list.append(line.rstrip())
    except:
        keep_sample_list = None
    marker_dict = parse_manifest(manifest_file, chroms, hg_refgenome)
    true_hap_file, true_legend_file = vcf2haplegend(vcf_file, keep_sample_list, marker_dict,true_output_prefix)
    # Dont need anymore bcuz vcf2haplegend was do that
    # set_marker_flags(
    #     true_legend_file,
    #     marker_dict,
    #     true_manifest_legend_file)
    prepare_test_hap(
        true_hap_file,
        true_legend_file,
        output_prefix)
    print('\nprepare data from vcf done!\n')

# if __name__ == '__main__':

#     parser = ArgumentParser(description='Prepare data for the model imputation ', add_help=True)

#     parser.add_argument('--mode',type=str,required=True,
#                         dest='mode',choices=['process','genotyping'],
#                         help='Chose what will you do.\n\tprocess: will call process data.\n\tgenotyping: will create genotyping data from vcf')
#     parser.add_argument('--vcf', type=str, required=True,
#                         dest='vcf_file', help='File vcf of data')
#     parser.add_argument('--omni', type=str, required=True,
#                         dest='manifest_file', help='Omni file data')
#     parser.add_argument('--hgref', type=str, required=True,
#                         dest='hg_refgenome', help='Genome Reference file')
#     parser.add_argument('--output_prefix', type=str, required=True,
#                         dest='output_prefix', help='Input file prefix')
#     parser.add_argument('--chroms',type=str,required=True,
#                         dest='chroms',help='chrom for extract information.\n File (with first line) or str by format "chromosome_id,chromosome_id,..."')
#     parser.add_argument('--test_sample', type=str, required=False, default=None,
#                         dest='test_sample_list_file', help='Sample will use for test. Each sample id on one line.')
#     args = parser.parse_args()

#     mode = args.mode
#     output_prefix = args.output_prefix
#     test_sample_list_file = args.test_sample_list_file
#     vcf_file = args.vcf_file
#     manifest_file = args.manifest_file
#     hg_refgenome = args.hg_refgenome
#     chroms_str = parser.chroms
#     chroms = None
#     if os.path.isfile(chroms_str):
#         with open(chroms_str,'r') as chroms_file:
#             line = chroms_file.readline()
#             chroms = line.strip().split(',')
#     else:
#         chroms = chroms_str.split(',')
#     if mode == 'process':
#         process_data_to_legend(vcf_file,
#                             manifest_file,
#                             hg_refgenome,
#                             chroms,
#                             output_prefix,
#                             test_sample_list_file = test_sample_list_file)
#     elif mode == 'genotyping':
#         genotyping_vcf(vcf_file, manifest_file, hg_refgenome, output_prefix, chroms)
