
import os
import re
import subprocess
import argparse

def pipeline_step_1(sample_name, R1_input, R2_input, unpaired_input, strandness, directory, Known_splice_site, GTF_file, ExonLength_file, sequence_type):
    # HISAT2Index = '/project/umw_chan_zhou/Data/GRCh38_hg38/HISAT2_index/genome'
    HISAT2Index = '/pi/chan.zhou-umw/Annotation/GRCh38_hg38/HISAT2_index/genome'
    # Known_splice_site = '/project/umw_chan_zhou/Data/GRCh38_hg38/gencode.v29.annotation.splice_site.txt'
    # GTF_file = '/project/umw_chan_zhou/Data/GRCh38_hg38/gencode.v29.annotation.gtf'
    # ExonLength_file = '/project/umw_chan_zhou/Data/GRCh38_hg38/gencode.v29.annotation.ExonLength'

    sambamba_view_path = directory + "/sambamba_view.sh"
    sambamba_sort_path = directory + "/sambamba_sort.sh"
    sed_1_path = directory + "/sed_1.sh"
    sed_2_path = directory + "/sed_2.sh"
    sed_3_path = directory + "/sed_3.sh"

    if strandness == "first":
        if sequence_type == "Paired":
            HISAT2_Strandness = "RF"
        elif sequence_type == "Single":
            HISAT2_Strandness = "R"
        HTSeq_Strandness = "reverse"
    elif strandness == "second":
        if sequence_type == "Paired":
            HISAT2_Strandness = "FR"
        elif sequence_type == "Single":
            HISAT2_Strandness = "F"
        HTSeq_Strandness = "yes"
    elif strandness == "unstrand":
        HTSeq_Strandness = "no"





    if not os.path.isdir('./dump'):
        os.system('mkdir dump')
    if not os.path.isdir('./log'):
        os.system('mkdir log')
    dump_child_path = "./dump/" + sample_name
    mkdir_dump_child = 'mkdir dump/' + sample_name
    if not os.path.isdir(dump_child_path):
        os.system(mkdir_dump_child)

    HISAT2_novel_splice_site_path = "dump/" + sample_name + "/HISAT2_novel_splice_site.txt"
    if sequence_type == "Paired":
        if strandness == "unstrand":
            run_hisat2_1 = subprocess.Popen(['hisat2', '-p', '10', '--dta', '-x', HISAT2Index,
                                             '-1', R1_input,'-2', R2_input,
                                             '--add-chrname', '--fr', '--known-splicesite-infile',
                                             Known_splice_site, '--novel-splicesite-outfile',HISAT2_novel_splice_site_path,
                                             '--novel-splicesite-infile', HISAT2_novel_splice_site_path, '--seed', '168', '--phred33',
                                             '--min-intronlen', '20', '--max-intronlen', '500000'], stdout=subprocess.PIPE,
                                            stderr=subprocess.PIPE)
        else:
            run_hisat2_1 = subprocess.Popen(['hisat2', '-p', '10', '--dta', '-x', HISAT2Index,
                                             '-1', R1_input,'-2', R2_input,
                                             '--add-chrname', '--rna-strandness', HISAT2_Strandness, '--fr', '--known-splicesite-infile',
                                             Known_splice_site, '--novel-splicesite-outfile',HISAT2_novel_splice_site_path,
                                             '--novel-splicesite-infile', HISAT2_novel_splice_site_path, '--seed', '168', '--phred33',
                                             '--min-intronlen', '20', '--max-intronlen', '500000'], stdout=subprocess.PIPE,
                                            stderr=subprocess.PIPE)

    elif sequence_type == "Single":
        if strandness == "unstrand":
            run_hisat2_1 = subprocess.Popen(['hisat2', '-p', '10', '--dta', '-x', HISAT2Index,
                                             '-U', unpaired_input,
                                             '--add-chrname', '--fr', '--known-splicesite-infile',
                                             Known_splice_site, '--novel-splicesite-outfile', HISAT2_novel_splice_site_path,
                                             '--novel-splicesite-infile', HISAT2_novel_splice_site_path, '--seed', '168',
                                             '--phred33',
                                             '--min-intronlen', '20', '--max-intronlen', '500000'], stdout=subprocess.PIPE,
                                            stderr=subprocess.PIPE)
        else:
            run_hisat2_1 = subprocess.Popen(['hisat2', '-p', '10', '--dta', '-x', HISAT2Index,
                                             '-U', unpaired_input,
                                             '--add-chrname', '--rna-strandness', HISAT2_Strandness, '--fr',
                                             '--known-splicesite-infile',
                                             Known_splice_site, '--novel-splicesite-outfile', HISAT2_novel_splice_site_path,
                                             '--novel-splicesite-infile', HISAT2_novel_splice_site_path, '--seed', '168',
                                             '--phred33',
                                             '--min-intronlen', '20', '--max-intronlen', '500000'], stdout=subprocess.PIPE,
                                            stderr=subprocess.PIPE)


    run_hisat2_2 = subprocess.Popen([sambamba_view_path, '-h -S -f bam -t 10 /dev/stdin'], stdin=run_hisat2_1.stdout,
                                    stdout=subprocess.PIPE, )

    run_hisat2_3_str = '--sort-by-name -t 10 --tmpdir dump/' + sample_name + ' -o dump/' + sample_name + '/' + sample_name + '.sortedbyname.bam /dev/stdin'

    run_hisat2_3 = subprocess.Popen([sambamba_sort_path, run_hisat2_3_str], stdin=run_hisat2_2.stdout, )

    run_hisat2_3.communicate()

    stderr_value = run_hisat2_1.communicate()[1]

    run_hisat2_3_HISAT2_output_path = 'log/' + sample_name + '.HISAT2'

    file(run_hisat2_3_HISAT2_output_path, 'w').write(stderr_value.encode('utf-8'))

    print("hisat2 done")

    run_sambamba_1_str = '-t 10 -h dump/' + sample_name + '/' + sample_name + '.sortedbyname.bam'
    run_sambamba_1 = subprocess.Popen([sambamba_view_path, run_sambamba_1_str],
                                      stdout=subprocess.PIPE, )

    run_sambamba_2 = subprocess.Popen([sed_1_path], stdin=run_sambamba_1.stdout, stdout=subprocess.PIPE, )

    run_sambamba_3 = subprocess.Popen([sed_2_path], stdin=run_sambamba_2.stdout, stdout=subprocess.PIPE, )

    run_sambamba_4 = subprocess.Popen([sed_3_path], stdin=run_sambamba_3.stdout, stdout=subprocess.PIPE, )

    run_sambamba_5_str = '-t 10 -f bam -S -o dump/' + sample_name + '/' + sample_name + '.sortedbyname.renamed.bam /dev/stdin'

    run_sambamba_5 = subprocess.Popen([sambamba_view_path, run_sambamba_5_str], stdin=run_sambamba_4.stdout)

    run_sambamba_5.communicate()

    # run_sambamba_3 = subprocess.Popen("echo 'SN:chrGL000225.1' | ./sed_3.sh", shell=True)

    run_sambamba_6_str = '-t 10 --tmpdir=dump/' + sample_name + ' -o dump/' + sample_name + '/' + sample_name + '.sortedbycoord.bam dump/' + sample_name + '/' + sample_name + '.sortedbyname.renamed.bam'

    run_sambamba_6 = subprocess.Popen([sambamba_sort_path, run_sambamba_6_str])

    run_sambamba_6.communicate()

    print("sambamba done")

    # os.system("rm -f dump/KT1/KT1.sortedbyname.bam dump/KT1/KT1.sortedbyname.renamed.bam")
    #
    # print("remove done")

    # Quantification

    run_htseq_count_str = 'dump/' + sample_name + '/' + sample_name + '.sortedbycoord.bam'
    run_htseq_count_str_2 = '--stranded=' + HTSeq_Strandness
    run_htseq_count = subprocess.Popen(
        ['htseq-count', '-f', 'bam', '-r', 'name', run_htseq_count_str_2, '-m', 'union', run_htseq_count_str,
         GTF_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    stdout_value, stderr_value = run_htseq_count.communicate()

    run_htseq_count_output_path = 'dump/' + sample_name + '/HTSEQ_count.txt'
    file(run_htseq_count_output_path, 'w').write(stdout_value.encode('utf-8'))

    run_htseq_count_err_path = 'log/' + sample_name + '.HTSEQ'
    file(run_htseq_count_err_path, 'w').write(stderr_value.encode('utf-8'))

    print("htseq-count done")

    # count='''head -1 log/KT1.HISAT2 | cut -d" " -f 1'''
    #
    # count_tem = os.popen(count).read()
    #
    # count_result = re.split('\n', count_tem)[0]

    count_result_path = 'log/' + sample_name + '.HISAT2'
    # with open(count_result_path, 'r') as lines:
    #     for line in lines:
    #         if 'reads; of these:' in line:
    #             count_result = line.split(' ')[0]
    #             break
    # with open(count_result_path, 'r') as fp:
    #     for i, line in enumerate(fp):
    #         if i == 3:
    #             count_result = line.split()[0]
    #             break
    if sequence_type == "Paired":
        with open(count_result_path, 'r') as lines:
            for line in lines:
                if 'aligned concordantly exactly 1 time' in line and "of these" not in line:
                    count_result = line.split()[0]
                    break
    elif sequence_type == "Single":
        with open(count_result_path, 'r') as lines:
            for line in lines:
                if 'aligned exactly 1 time' in line:
                    count_result = line.split()[0]
                    break

    print(line)
    print(count_result)

    #
    #
    # with open(count_result_path, 'r') as lines:
    #     for line in lines:
    #         if 'aligned exactly 1 time' in line:
    #             count_result = line.split()[0]
    #             break



    run_perl_str = 'dump/' + sample_name + '/HTSEQ_count.txt'
    rpkm_path = directory + '/src/rpkm.pl'
    run_perl = subprocess.Popen(['perl', rpkm_path, run_perl_str, ExonLength_file, count_result],
                                stdout=subprocess.PIPE, )

    stdout_value = run_perl.communicate()[0]

    run_perl_output_path = 'dump/' + sample_name + '/Gene_FPKM.txt'

    file(run_perl_output_path, 'w').write(stdout_value.encode('utf-8'))

    print("perl done")




    with open('done.txt', 'w') as rsh:
        rsh.write("")

    rm_str = "dump/" + str(sample_name) + "/*.bam*"
    tem_str = "rm -rf " + rm_str
    os.system(tem_str)

    print("All done")


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--sample_name", default=None, type=str, required=True, help="sample name")
    parser.add_argument("--R1_input", default=None, type=str, required=True, help="R1 input path")
    parser.add_argument("--R2_input", default=None, type=str, required=True, help="R2 input path")
    parser.add_argument("--unpaired_input", default=None, type=str, required=True, help="unpaired input path")
    parser.add_argument("--strandness", choices=["first", "second", "unstrand"], default=None, required=True,
                        help="strandness type")
    parser.add_argument("--directory", default=None, type=str, required=True, help="directory path")
    parser.add_argument("--Known_splice_site", default=None, type=str, required=True, help="Known_splice_site path")
    parser.add_argument("--GTF_file", default=None, type=str, required=True, help="GTF_file path")
    parser.add_argument("--ExonLength_file", default=None, type=str, required=True, help="ExonLength_file path")
    parser.add_argument("--sequence_type", default=None, type=str, required=True, help="sequence type")
    args = parser.parse_args()

    pipeline_step_1(sample_name=args.sample_name, R1_input=args.R1_input, R2_input=args.R2_input, unpaired_input=args.unpaired_input, strandness=args.strandness, directory=args.directory,
                    Known_splice_site=args.Known_splice_site, GTF_file=args.GTF_file, ExonLength_file=args.ExonLength_file, sequence_type=args.sequence_type)

if __name__ == '__main__':
    main()


