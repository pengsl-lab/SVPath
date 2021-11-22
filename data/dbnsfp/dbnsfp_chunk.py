import argparse

import pysam
import multiprocessing
import os

def get_file_list(chromName, refStart, refEnd, step):
    file_list = []
    for bed in range(refStart, refEnd, step):
        start = bed
        end = bed + step - 1
        if end > refEnd:
            end = refEnd
        file_list.append([chromName, start, end])
    return file_list

def get_sv_dbnsfp_file_list(chrom, sv_start, sv_end, step, ref_length, dbnsfp_dir):
    bed_start = sv_start // step * step
    bed_end = sv_end // step * step + step
    if bed_end > ref_length:
        bed_end = ref_length
    sv_file_list = []
    while bed_start <= sv_end:
        end = bed_start + step
        if end > ref_length:
            end = ref_length
        file_name = os.path.join(dbnsfp_dir, 'dbnsfp_' + chrom + '-' + str(bed_start + 1) + '-' + str(end) + '.txt')
        sv_file_list.append(file_name)
        bed_start += step
    return sv_file_list

def get_chr_list(reference):
    chrom_list = [str(i) for i in range(1, 23, 1)]
    chrom_list.append('X')
    chrom_list.append('Y')

    refFasta = pysam.FastaFile(reference)
    chr_list = []
    for chrom in chrom_list:
        chr_list.append(['chr' + chrom, 1, refFasta.get_reference_length('chr' + chrom), 'dbnsfp_'+chrom+'.txt'])
    return chr_list

def para_chunking(chrom):
    chromName = chrom[0]
    chromRef_start = chrom[1]
    chromRef_end = chrom[2]
    chrom_file = chrom[3]
    step = 1000000

    chunk_list = get_file_list(chromName, chromRef_start, chromRef_end, step)
    with open(chrom_file, 'r') as dbnsfp_reader:
        flag = True
        tmpLine = ''
        for chunk in chunk_list:
            chromName = chunk[0]
            bed_start = chunk[1]
            bed_end = chunk[2]
            print(chunk)
            chunk_file_name = 'dbnsfp_' + chromName + '-' + str(bed_start) + '-' + str(bed_end) + '.txt'

            with open(chunk_file_name, 'w') as chunk_writer:
                while True:
                    if flag == True:
                        line = dbnsfp_reader.readline()
                        tmpLine = line
                    else:
                        line = tmpLine
                    if line:
                        line = line.split()
                        if 'chr' + line[0] == chromName and int(line[1]) >= bed_start and int(line[1]) <= bed_end:
                            for i in line:
                                chunk_writer.write(i + '\t')
                            chunk_writer.write('\n')
                            flag = True
                        else:
                            flag = False
                            break
                    else:
                        break

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Splitting dbnsfp file')
    parser.add_argument('-r', '--reference', required=True, help='Reference genome file')
    parser.add_argument('-p', '--processes', default=8, help='Number of processes')

    args = parser.parse_args()
    refer = args.reference
    num_process = args.processes
    chrom_list = get_chr_list(refer)
    pool = multiprocessing.Pool(processes=num_process)
    result = pool.map(para_chunking, chrom_list)
    pool.close()
    pool.join()