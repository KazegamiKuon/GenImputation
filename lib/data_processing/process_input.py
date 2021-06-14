# source ./script/test_data_preparation.py from https://github.com/kanamekojima/rnnimp

from argparse import ArgumentParser
from io import TextIOWrapper
import json
from lib import config
import os
import sys
import re
from tqdm import tqdm
import pandas as pd
import numpy as np
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


def parse_manifest(manifest_file, chr_nums: list, hg_refgenome):
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


def mapping_marker(chrom, position, ref, alts, marker_dict):
    flag = legend_config.unobserve
    if chrom in marker_dict:
        maker_position = marker_dict[chrom]
        if position in maker_position:
            for alleles in maker_position[position]:
                if alleles[0] == ref and alleles[1] in alts:
                    flag = legend_config.observe
                    break
                if alleles[1] in alts and alleles[0] == ref:
                    flag = legend_config.observe
                    break
    return flag


# def set_marker_flags(legend_file, marker_dict, output_file):
#     g.mkdir(os.path.dirname(output_file))
#     with g.reading(legend_file) as fp, \
#             g.writing(output_file) as w_fp:
#         w_fp.write(fp.readline().rstrip() + ' ' +
#                    legend_config.array_marker_flag+'\n')
#         for line in tqdm(fp, desc='set marker flags'):
#             items = g.line_to_items(line)
#             chrom = items[legend_config.chrom_index]
#             position = items[legend_config.position_index]
#             ref = items[legend_config.reference_index]
#             alts = [items[legend_config.alternate_allele_index]]
#             flag = mapping_marker(chrom, position, ref, alts, marker_dict)
#             w_fp.write('{:s} {:s}\n'.format(line.rstrip(), flag))


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

# def legend_to_region(legend_file,hap_file,number_bin,inter_bin_width_percen:float = 0,no_observation=False):
#     legend_data = pd.read_csv(legend_file,sep=legend_config.legend_split_params)    
#     nb_line = len(legend_data)
#     bin_col = 'bin'    
#     legend_data[bin_col] = pd.qcut(np.arange(nb_line),number_bin,labels=False)
#     labels, nb_values = np.unique(legend_data[bin_col],return_counts=True)
#     values,nb_bin = np.unique(nb_values,return_counts=True)
#     print('Info about number variant and number bin each variant\n')
#     print(pd.DataFrame({'nb variants':values,'nb bin':nb_bin}))
#     print('Agree to div bin by this? [y/n]')
#     check = input()
#     if check.startswith('n') or check.startswith('N'):
#         print('Stop process data.')
#         return
#     with g.reading(hap_file) as hap_file_stream:
#         # label are sorted and orderly in data.
#         # When we loop over label each line (exclude header) at legend file is map with line at hap file.
#         for i, label in enumerate(tqdm(labels,desc='create region data from legend')):
#             # get new data at this loop
#             bin_data = legend_data[legend_data[bin_col] == label] # this bin data
#             legend_region_file, legend_region_true_file = legend_config.get_hap_legend_region_file_name(legend_file,label)
#             hap_region_file, hap_region_true_file = legend_config.get_hap_legend_region_file_name(hap_file,label)
#             nb_bin_line = len(bin_data)
#             # Ghi ra file hap
#             with g.writing(hap_region_file) as hap_region_file_stream, \
#                 g.writing(hap_region_true_file) as hap_region_true_file_stream:
#                 for i in tqdm(range(nb_bin_line),desc='process region {:03d}'.format(label)):
#                     line = hap_file_stream.readline()
#                     # Nếu flag là observe thì ghi vô cả 2 file
#                     flag = str(temp.iloc[i][legend_config.array_marker_flag])
#                     if flag == legend_config.observe:
#                         hap_region_file_stream.write(line)
#                         if no_observation:
#                             continue
#                     hap_region_true_file_stream.write(line)
#             if no_observation:
#                 temp = temp[temp[legend_config.array_marker_flag] == int(legend_config.unobserve)]
#             temp.to_csv(legend_region_true_file,sep = legend_config.legend_split_params)

def genotyping_vcf(vcf_file, manifest_file, hg_refgenome, output_prefix, chroms):

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


def process_data_to_legend(vcf_file, manifest_file, hg_refgenome, chroms, output_prefix, test_sample_list_file=None):

    g.check_required_software('samtools')
    g.check_required_file(vcf_file)
    g.check_required_file(manifest_file)
    g.check_required_file(hg_refgenome)

    true_output_prefix = g.get_true_file_prefix(output_prefix)

    keep_sample_list = []
    try:
        with open(test_sample_list_file) as fp:
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
