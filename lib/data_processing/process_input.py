# source ./script/test_data_preparation.py from https://github.com/kanamekojima/rnnimp

from argparse import ArgumentParser
from io import FileIO, TextIOWrapper
import os
import sys
import re
from numpy import not_equal, split
from numpy.core.fromnumeric import sort
from tqdm import tqdm
from lib.config import config_class
from lib.utils import general as g
from lib.config.config_class import VCFConfigClass

vcf_config = VCFConfigClass()

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

def get_sequence(region, hg_fasta_file):
    command = 'samtools faidx {:s} {:s}'.format(hg_fasta_file, region)
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

def parse_manifest(manifest_file, chr_nums:list, hg_fasta_file):
    marker_dict = {}
    with g.reading(manifest_file) as fp:
        for line in fp:
            if line.startswith('[Assay]'):
                break
        line = fp.readline()
        items = g.line_to_items(line,split=',')
        IlmnStrand_col = items.index('IlmnStrand')
        Chr_col = items.index('Chr')
        MapInfo_col = items.index('MapInfo')
        SNP_col = items.index('SNP')
        SourceStrand_col = items.index('SourceStrand')
        SourceSeq_col = items.index('SourceSeq')
        RefStrand_col = items.index('RefStrand')
        for line in tqdm(fp,desc='create marker from manifest'):
            if line.startswith('[Controls]'):
                break
            items = g.line_to_items(line,',')
            chrom = items[Chr_col]            
            if chrom not in chr_nums:
                continue
            if chrom not in marker_dict:
                marker_dict[chrom] = {}
            
            position = items[MapInfo_col]
            alleles = items[SNP_col][1:-1].split('/')
            if alleles[0] == 'I' or alleles[0] == 'D':
                source_seq = items[SourceSeq_col]
                upstream_seq, a0, a1, downstream_seq = re.split(
                    '[\[\/\]]', source_seq)
                upstream_seq = upstream_seq.upper()
                downstream_seq = downstream_seq.upper()
                if a0 == '-':
                    a0 = ''
                if a1 == '-':
                    a1 = ''
                a0_seq = upstream_seq + a0 + downstream_seq
                a1_seq = upstream_seq + a1 + downstream_seq
                margin = 10
                region_start = int(position) - len(upstream_seq) - margin
                region_end = region_start + margin \
                             + max(len(a0_seq), len(a1_seq)) + margin
                region = 'chr{:s}:{:d}-{:d}'.format(
                    chrom, region_start, region_end)
                ref_seq = get_sequence(region, hg_fasta_file)
                a0_align = align(ref_seq, a0_seq)
                a1_align = align(ref_seq, a1_seq)
                indel_position = max(a0_align, a1_align) + len(upstream_seq) \
                                 + region_start - 1
                alleles[0] = upstream_seq[-1] + a0
                alleles[1] = upstream_seq[-1] + a1
                position = str(indel_position)
            else:
                if items[RefStrand_col] == '-':
                    alleles[0] = complement_base(alleles[0])
                    alleles[1] = complement_base(alleles[1])
            if position not in marker_dict[chrom]:
                marker_dict[chrom][position] = []
            marker_dict[chrom][position].append(alleles)
    return marker_dict

def mapping_marker(chrom,position,ref,alts,marker_dict):
    flag = '0'
    if chrom in marker_dict:
        maker_position = marker_dict[chrom]
        if position in maker_position:
            for alleles in maker_position[position]:
                if alleles[0] == ref and alleles[1] in alts:
                    flag = '1'
                    break
                if alleles[1] in alts and alleles[0] == ref:
                    flag = '1'
                    break
    return flag

def set_marker_flags(chrom,legend_file, marker_dict, output_file):
    g.mkdir(os.path.dirname(output_file))
    with g.reading(legend_file) as fp, \
         g.writing(output_file) as w_fp:
        w_fp.write(fp.readline().rstrip() + ' array_marker_flag\n')
        for line in tqdm(fp,desc='set marker flags'):
            items = g.line_to_items(line)
            position = items[1]
            ref = items[2]
            alts = [items[3]]
            flag = mapping_marker(chrom,position,ref,alts,marker_dict)
            w_fp.write('{:s} {:s}\n'.format(line.rstrip(), flag))

def get_snp_id(snp,chrom,position,ref,alt,nb_alt):
    snp_id = snp
    if snp_id == '.':
        snp_id = '{:s}:{:s}:{:s}:{:s}'.format(
            chrom, position, ref, alt)
    elif nb_alt >= 2:
        snp_id += ':{:s}:{:s}'.format(ref, alt)
    return snp_id

def get_snp_id_from_line(line,header_dict):
    chrom, position, snp, ref, alts, info_dict, genotypes = get_data_from_line(line,header_dict)
    nb_alt = len(alts)
    snp_ids = [get_snp_id(snp,chrom,position,ref,alt,nb_alt) for alt in alts]
    return snp_ids

def get_data_from_line(line,header_dict):
    items = g.line_to_items(line,vcf_config.vcf_data_split_params)
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
    info_dict = g.list_expression_to_dict(items[header_dict[vcf_config.info_col]].split(';'))
    # Dữ liệu về genotypes từ index 9 đổ về sau
    genotypes = items[vcf_config.vcf_header_sample_startindex:]
    return chrom, position, snp, ref, alts, info_dict, genotypes

def vcf2haplegend(vcf_file, keep_sample_list, output_prefix):
    def is_missing_genotype(genotype):
        if genotype == '.':
            return True
        if genotype == './.':
            return True
        if genotype == '.|.':
            return True
        return False

    g.mkdir(os.path.dirname(output_prefix))
    with g.reading(vcf_file) as fp:
        header_dict = process_vcf_header(fp)
        samples = list(header_dict[vcf_config.header_dict_sample_key].keys())
        # sort sample by index
        samples.sort(key=lambda sample: header_dict[vcf_config.header_dict_sample_key][sample])
        if keep_sample_list is None:
            keep_sample_id_list = list(range(len(samples)))
        else:
            keep_sample_id_list = [
                samples.index(sample) for sample in keep_sample_list
                if sample in samples
            ]
        with open(output_prefix + '.sample', 'wt') as s_fp:
            s_fp.write('ID_1 ID_2 missing\n')
            s_fp.write('0 0 0\n')
            for sample_id in keep_sample_id_list:
                s_fp.write(
                    '{0:s} {0:s} 0\n'.format(samples[sample_id]))
        sample_size = len(keep_sample_id_list)
        with g.writing(output_prefix + '.hap.gz') as h_fp, \
             g.writing(output_prefix + '.legend.gz') as l_fp:
            l_fp.write('id position ref alt af maf\n')
            allele_id_list = [None] * (2 * sample_size)
            # danh sách ALT ID
            hap_record = [None] * (2 * sample_size)
            # Haplotype record
            for line in tqdm(fp,desc='vcf to haplegend'):
                # đọc line và tách ra thành từng item
                chrom, position, snp, ref, alts, info_dict, genotypes = get_data_from_line(line,header_dict)
                # Loop qua danh sách sample và lấy dữ liệu allele_id của nó
                for i, sample_id in enumerate(keep_sample_id_list):
                    genotype = genotypes[sample_id]
                    if is_missing_genotype(genotype):
                        allele_id_list[2 * i] = None
                        allele_id_list[2 * i + 1] = None
                    else:
                        allele_id_pair = map(int, genotype.split('|'))
                    for j, allele_id in enumerate(allele_id_pair):
                        allele_id_list[2 * i + j] = allele_id
                nb_alt = len(alts)
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
                    h_fp.write(' '.join(hap_record))
                    h_fp.write('\n')
                    snp_id = get_snp_id(snp,chrom,position,ref,alt,nb_alt)
                    # convert to float AF and calculate MAF
                    af = float(info_dict['AF'][i])
                    maf = af if af <= 0.5 else 1 - af
                    # convert to str
                    af = str(af)
                    maf = str(maf)
                    # write line
                    l_fp.write('{:s} {:s} {:s} {:s} {:s} {:s}\n'.format(
                        snp_id, position, ref, alt, af, maf))

def prepare_test_hap(hap_file, legend_file, output_prefix):
    array_marker_flag_list = []
    with g.reading(legend_file) as fp, \
         g.writing(output_prefix + '.legend.gz') as w_fp:
        line = fp.readline()
        items = g.line_to_items(line)
        try:
            array_marker_flag_col = items.index('array_marker_flag')
        except ValueError:
            g.print_error()
            g.print_error('Error: Header "array_marker_flag" not found in '
                        + legend_file)
            g.print_error()
            sys.exit(0)
        w_fp.write(line)
        for line in tqdm(fp,desc='prepare marker to legend'):
            items = g.line_to_items(line)
            array_marker_flag = items[array_marker_flag_col] == '1'
            array_marker_flag_list.append(array_marker_flag)
            if array_marker_flag:
                w_fp.write(line)
    with g.reading(hap_file) as fp, \
         g.writing(output_prefix + '.hap.gz') as w_fp:
        for i, line in enumerate(tqdm(fp,desc='prepare marker to hap')):
            if array_marker_flag_list[i]:
                w_fp.write(line)

def process_vcf_header(rf:TextIOWrapper,wf:TextIOWrapper = None,write_callback = None):
    '''
    input:
        rf (TextIOWrapper): read file
        wf (TextIOWrapper): write file, default is None
        write_callback (str,func): define what you wanna write in file. Default is None:
            if func, it will get param are current line and previous line write_callback(line,old_line)
    '''
    header_dict=None
    sample_dict=None
    old_line = None
    for line in tqdm(rf,desc='process header vcf'):
        # reading region
        if line.startswith(vcf_config.vcf_header_line_startswith):
            items = g.line_to_items(line,vcf_config.vcf_data_split_params)
            header_dict = g.list_to_dict(items[:vcf_config.vcf_header_sample_startindex])
            sample_dict = g.list_to_dict(items[vcf_config.vcf_header_sample_startindex:],vcf_config.vcf_header_sample_startindex)
            header_dict[vcf_config.header_dict_sample_key] = sample_dict
            break
        # writing region
        if wf is None:
            continue
        w_line = line
        if type(write_callback) is type(lambda x: x):
            w_line = write_callback(line,old_line)
        elif type(write_callback) is str:
            w_line = write_callback
        wf.write(w_line)
        old_line = line
    return header_dict

def genotyping_vcf(vcf_file,manifest_file,hgref,output_prefix,chroms):
    marker_dict = parse_manifest(manifest_file,chroms,hgref)
    with g.reading(vcf_file) as rf, \
         g.writing(output_prefix+'.vcf.gz') as wf:
        header_dict = process_vcf_header(rf,wf)
        for line in tqdm(rf,desc='mapping vcf with manifest file'):
            chrom, position, _, ref, alts, _, _ = get_data_from_line(line,header_dict)
            flag = mapping_marker(chrom,position,ref,alts,marker_dict)
            if flag == '1':
                wf.write(line)

def process_data(vcf_file,manifest_file, hg_fasta_file, chroms,hap_prefix, true_hap_prefix, test_sample_list_file=None):
    g.check_required_software('samtools')
    g.check_required_file(vcf_file)
    g.check_required_file(manifest_file)
    g.check_required_file(hg_fasta_file)
    keep_sample_list = []
    try:
        with open(test_sample_list_file) as fp:
            for line in fp:
                keep_sample_list.append(line.rstrip())
    except:
        keep_sample_list = None
    vcf2haplegend(vcf_file, keep_sample_list, true_hap_prefix)
    marker_dict = parse_manifest(manifest_file, chroms, hg_fasta_file)
    set_marker_flags(
        true_hap_prefix + '.legend.gz', marker_dict,
        true_hap_prefix + '.omni.legend.gz')
    prepare_test_hap(
        true_hap_prefix + '.hap.gz', true_hap_prefix + '.omni.legend.gz',
        hap_prefix)

if __name__ == '__main__':
    
    parser = ArgumentParser(description='Prepare data for the model imputation ', add_help=True)
    
    parser.add_argument('--mode',type=str,required=True,
                        dest='mode',choices=['process','genotyping'],
                        help='Chose what will you do.\n\tprocess: will call process data.\n\tgenotyping: will create genotyping data from vcf')
    parser.add_argument('--vcf', type=str, required=True,
                        dest='vcf_file', help='File vcf of data')
    parser.add_argument('--omni', type=str, required=True,
                        dest='manifest_file', help='Omni file data')
    parser.add_argument('--hgref', type=str, required=True,
                        dest='hg_fasta_file', help='Genome Reference file')
    parser.add_argument('--output_prefix', type=str, required=True,
                        dest='output_prefix', help='Input file prefix')
    parser.add_argument('--chroms',type=str,required=True,
                        dest='chroms',help='chrom for extract information.\n File (with first line) or str by format "chromosome_id,chromosome_id,..."')
    parser.add_argument('--test_sample', type=str, required=False, default=None,
                        dest='test_sample_list_file', help='Sample will use for test. Each sample id on one line.')
    args = parser.parse_args()
    
    mode = args.mode
    output_prefix = args.output_prefix
    test_sample_list_file = args.test_sample_list_file
    vcf_file = args.vcf_file
    manifest_file = args.manifest_file
    hg_fasta_file = args.hg_fasta_file
    chroms_str = parser.chroms
    chroms = None
    if os.path.isfile(chroms_str):
        with open(chroms_str,'r') as chroms_file:
            line = chroms_file.readline()
            chroms = line.strip().split(',')
    else:
        chroms = chroms_str.split(',')
    if mode == 'process':
        process_data(vcf_file,
                    manifest_file,
                    hg_fasta_file,
                    chroms,
                    output_prefix,
                    output_prefix+'_true',
                    test_sample_list_file = test_sample_list_file)
    elif mode == 'genotyping':
        genotyping_vcf(vcf_file, manifest_file, hg_fasta_file, output_prefix, chroms)