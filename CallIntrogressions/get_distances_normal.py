import os,sys,re,gzip
import argparse
import numpy as np

def read_vcf(file_):
    opener=gzip.open if file_.endswith(".gz") else open
    vcf_dict={}
    with opener(file_,'rt') as fp:
        for line in fp:
            line=line.strip()
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    acc_vec = line.split()[9:]
                continue
            fields=line.split()
            tags = fields[7].split(";")
            start = int(fields[1])

            # Get the support vector
            supp_vec = None
            for j in tags:
                if j.startswith("SUPP_VEC="):
                    supp_vec = [int(i) for i in j[9:]]
            if supp_vec is None:
                try:
                    raise ValueError("Missing 'SUPP_VEC' field")
                except ValueError:
                    continue
            vcf_dict.setdefault(fields[0],{}).setdefault(start,supp_vec)
    return acc_vec,vcf_dict            

parser = argparse.ArgumentParser(description='Get Jaccard similarity between one_group and others.')
parser.add_argument("vcf", metavar="<SVs.vcf>", type=str, help="SV vcf file with support vectors.")
parser.add_argument("species_file", metavar="<group.txt>", type=str, help="Second column is the phylogenetic group, first column is the accession")
parser.add_argument("species", metavar="<SP>", type=str, help="Group to compare to other groups")
parser.add_argument("fai", metavar="<reference.fasta.fai>", type=str, help="Fasta index file for the reference genome")
parser.add_argument("w", metavar="<100000>", type=int, default=1000000, help="Introgression window size.")
parser.add_argument("-m", metavar="5", type=int, default=5, help='minimum number of SVs needed to calculate Jaccard')
parser.add_argument("-o", metavar="<out_dir>", type=str, default=".", help='the out dir of saving results')

args = parser.parse_args()
vcf_file = args.vcf
fai_file = args.fai
species_file = args.species_file
target_species = args.species
window_size = args.w
min_den = args.m
out_dir= args.o

if min_den < 2:
    raise ValueError("-m must be at least 2")

# Associate the list of samples with a species
comp_species_dict={}
all_samples=[]
with open(species_file) as f:
    for line in f:
        line=line.strip()
        if not line:
            continue
        acc ,species = line.split()
        comp_species_dict.setdefault(species,[]).append(acc)
        all_samples.append(acc)
sll=comp_species_dict[target_species]

# Get the chromosome sizes
chr_lens = {}
with open(fai_file, "r") as f:
    for line in f:
        header, length, x, y, z = line.rstrip().split("\t")
        chr_lens[header] = int(length)

acc_vec,vcf_dict=read_vcf(vcf_file)
overs=set(all_samples)-set(acc_vec)
if len(overs) > 0:
    raise ValueError("Missing samples in vcf: %s"%(overs))

for chr_used in vcf_dict:
    poses=sorted(vcf_dict[chr_used].keys())
    n_windows = chr_lens[chr_used] // window_size
    for comp_species in comp_species_dict:
        if comp_species==target_species:
            continue
        distances = np.zeros((len(sll), n_windows))
        comp_max_accs = np.zeros((len(sll), n_windows), dtype=np.int32)
        current_window = 0
        supp_matrix = []
        for start in poses:
            widx = start // window_size
            supp_vec=vcf_dict[chr_used][start]
            # Build the support matrix for this window or start a new matrix for a new window
            if widx == current_window:
                supp_matrix.append(supp_vec)
            if widx > current_window or start==poses[-1]:
                # We have moved on to the next window. Save the distances for the finished window
                sm = np.asarray(supp_matrix)
                # For now, I will calculate the distances in a 'for' loop. Perhaps vectorize in the future
                # Iterate over the SLLs
                t_distances = [[] for i in range(len(sll))]
                if len(sm)==0:
                    pass
                else:
                    for i in range(len(sll)):
                        this_acc = sll[i]
                        supp_idx = acc_vec.index(this_acc)
                        this_vec = sm[:, supp_idx]
                    # Iterate over the comp species:
                        for comp_acc in comp_species_dict[comp_species]:
                            comp_supp_idx = acc_vec.index(comp_acc)
                            this_comp_vec = sm[:, comp_supp_idx]
                            if np.count_nonzero(this_comp_vec) >= min_den and np.count_nonzero(this_vec) >= min_den:
                                # Get the Jaccard distance
                                num = np.count_nonzero(np.logical_and(this_vec, this_comp_vec))  # Intersection
                                den = np.count_nonzero(np.logical_or(this_vec, this_comp_vec))  # Union
                                t_distances[i].append(num / den)
                            else:
                                t_distances[i].append(-1)
                    # Find which comp sample gave the max
                    t_distances_argmax = np.asarray([np.argmax(i) for i in t_distances])
                    comp_max_accs[:, current_window] = t_distances_argmax
                    # Get the max % shared SVs between a given SLL and each comp species sample
                    t_distances = np.asarray([np.max(i) for i in t_distances])
                    distances[:, current_window] = t_distances
                # Now that we have calculated the distances for the finished window, start the next one
                if widx == n_windows:
                    break
                current_window = widx
                supp_matrix = []
                supp_matrix.append(supp_vec)
        with open("%s/%s.%s_%s.max_map.txt"%(out_dir,chr_used,target_species,comp_species),'w') as f:
            f.write("Sample\t" + "\t".join( [str(i*window_size) for i in range(n_windows)]) + "\n")
            for i in range(len(sll)):
                f.write(sll[i] + "\t" + "\t".join( [comp_species_dict[comp_species][j] for j in list(comp_max_accs[i, :])] ) + "\n")
        with open("%s/%s.%s_%s.max_value.txt"%(out_dir,chr_used,target_species,comp_species),'w') as f:
            print("Sample\t" + "\t".join( [str(i*window_size) for i in range(n_windows)] ),file=f)
            for i in range(len(sll)):
                print(sll[i] + "\t" + "\t".join( [str(j) for j in list(distances[i, :])] ).replace("-1.0", "NA"),file=f) 
