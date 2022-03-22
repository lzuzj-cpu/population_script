import os,sys,re
import logging
import intron_window
gff_file=sys.argv[1]
try:
    reg=int(sys.argv[2]) #Regulation area
except:
    reg=5000

gene_dict={}
#{chr:[[gene_start,gene_end,gene_name]]}
cds_dict={}
#{gene_name:[[cds_start,cds_end]]}
with open(gff_file) as fp:
    for line in fp:
        line=line.strip()
        if not line or re.search("^#",line):
            continue
        line=line.split()
        if line[2]=="mRNA":
            gene_name=re.search("ID=([^;]*)",line[8]).group(1)
            gene_dict.setdefault(line[0],[]).append([int(line[3]),int(line[4]),gene_name])
        elif line[2]=="CDS":
            gene_name=re.search("Parent=([^;]*)",line[8]).group(1)
            cds_dict.setdefault(gene_name,[]).append([int(line[3]),int(line[4])])

def find_intron(contig):
    split_window=[]
    #[[chr,start,end,type]]
    chr_genes=gene_dict[contig]
    for gene in chr_genes:
        gene_cds=cds_dict[gene[2]]
        gene_introns=intron_window.intron_wins(gene,gene_cds)
        if type(gene_introns) != list:
            logging.warning('gene: %s is wrong, the cds and intron in this gene will not be writed',gene[2])
            continue
        for cds in gene_cds:
            split_window.append([contig,str(cds[0]),str(cds[1]),"gene=%s;type=CDS"%(gene[2])])
        if gene_introns:
            for intron in gene_introns:
                split_window.append([contig,str(intron[0]),str(intron[1]),"gene=%s;type=INTRON"%(gene[2])])
    return split_window

def find_Regulation_area(contig):
    split_window=[]
    chr_genes=gene_dict[contig]
    chr_genes.sort(key=lambda x:x[0])
    last_gene=[1]
    for gene in chr_genes:
        if last_gene==[1]:
            if gene[0]>1:
                area_start=1 if reg>=gene[0] else gene[0]-reg
                area_end=gene[0]-1
                split_window.append([contig,str(area_start),str(area_end),"gene=%s;type=REGULATION"%(gene[2])])
            last_gene=gene
        else:
            if gene[0]>last_gene[1]+1:
                if gene[0]-last_gene[1]-1<=reg:
                    area_start=last_gene[1]+1
                    area_end=gene[0]-1
                    split_window.append([contig,str(area_start),str(area_end),"gene=%s;type=REGULATION|gene=%s;type=REGULATION"%(last_gene[2],gene[2])])
                elif gene[0]-last_gene[1]-1>=2*reg:
                    area_start=last_gene[1]+1
                    area_end=area_start+reg-1
                    split_window.append([contig,str(area_start),str(area_end),"gene=%s;type=REGULATION"%(last_gene[2])])
                    area_end=gene[0]-1
                    area_start=area_end-reg+1
                    split_window.append([contig,str(area_start),str(area_end),"gene=%s;type=REGULATION"%(gene[2])])
                elif gene[0]-last_gene[1]-1>reg and gene[0]-last_gene[1]-1<2*reg:
                    over_len=2*reg-(gene[0]-last_gene[1]-1)
                    area_start=last_gene[1]+1
                    area_end=area_start+reg-over_len-1
                    split_window.append([contig,str(area_start),str(area_end),"gene=%s;type=REGULATION"%(last_gene[2])])
                    area_start=area_end+1
                    area_end=area_start+over_len-1
                    split_window.append([contig,str(area_start),str(area_end),"gene=%s;type=REGULATION|gene=%s;type=REGULATION"%(last_gene[2],gene[2])])
                    area_start=area_end+1
                    area_end=gene[0]-1
                    split_window.append([contig,str(area_start),str(area_end),"gene=%s;type=REGULATION"%(gene[2])])
            last_gene=gene
    area_start=last_gene[1]+1
    area_end=area_start+reg-1
    split_window.append([contig,str(area_start),str(area_end),"gene=%s;type=REGULATION"%(last_gene[2])])
    return split_window

all_split=[]
for contig in gene_dict:
    all_split=all_split+find_intron(contig)+find_Regulation_area(contig)
all_split.sort(key=lambda x:[x[0],int(x[1])])
print("#CHROM\tSTART\tEND\tTYPE")
for win in all_split:
    print("\t".join(win))
