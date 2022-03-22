import os,sys,re
gff_window_file=sys.argv[1]
window_result_file=sys.argv[2]
gff_window={}
with open(gff_window_file) as fp:
    for line in fp:
        line=line.strip()
        if re.search('^#',line) or not line:
            continue
        line=line.split()
        gff_window['-'.join(line[0:3])]=line[3]
print("#Chrom\tStart\tEnd\tD\tFdm\tnumerator\tdenominator\tFDnumerator\tFDMdenominator\tAll_snp\tAvaible_snpType")
with open(window_result_file) as fp:
    for line in fp:
        line=line.strip()
        if re.search('^#',line) or not line:
            continue
        line=line.split()
        print('\t'.join(line+[gff_window['-'.join(line[0:3])]]))
