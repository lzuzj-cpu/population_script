import os,sys,re
from decimal import *
from multiprocessing import Pool

window_gff_result=sys.argv[1]

def read_result(file_):
    with open(file_) as fp:   
        D_genes={}
        for line in fp:
            line=line.strip()
            if not line or re.search("^#",line):
                continue
            line=line.split()
            Data=[Decimal(x) for x in line[5:11]]
            for i in line[11].split("|"):
                gene_name=re.search('gene=([^;]*)',i).group(1)
                gene_type=re.search('type=(.*)$',i).group(1)
                D_genes.setdefault(gene_name,{}).setdefault(gene_type,[]).append(Data)
    return D_genes

def gene_D_Fdm(gene):
    Anumerator,Adenominator,AFDnumerator,AFDMdenominator=0,0,0,0
    AD,AFDM="NA","NA"
    AAll_snp,AAvaible_snp=0,0
    region={"CDS":{"D":"NA","FDM":"NA","All_snp":0,"Avaible_snp":0},"INTRON":{"D":"NA","FDM":"NA","All_snp":0,"Avaible_snp":0},"REGULATION":{"D":"NA","FDM":"NA","All_snp":0,"Avaible_snp":0}}
    for gene_type in D_genes[gene]:
        numerator,denominator,FDnumerator,FDMdenominator=0,0,0,0
        D,FDM="NA","NA"
        All_snp,Avaible_snp=0,0
        for data in D_genes[gene][gene_type]:
            numerator+=data[0]
            denominator+=data[1]
            FDnumerator+=data[2]
            FDMdenominator+=data[3]
            All_snp+=data[4]
            Avaible_snp+=data[5]
        Anumerator+=numerator
        Adenominator+=denominator
        AFDnumerator+=FDnumerator
        AFDMdenominator+=FDMdenominator
        AAll_snp+=All_snp
        AAvaible_snp+=Avaible_snp
        region[gene_type]['All_snp']=All_snp
        region[gene_type]['Avaible_snp']=Avaible_snp
        if denominator !=0:
            region[gene_type]['D']=numerator/denominator
        if FDMdenominator !=0:
            region[gene_type]['FDM']=FDnumerator/FDMdenominator
    if Adenominator !=0:
        AD=Anumerator/Adenominator
    if AFDMdenominator !=0:
        AFDM=AFDnumerator/AFDMdenominator
    return [gene,AD,AFDM,AAll_snp,AAvaible_snp,region["CDS"]["D"],region["CDS"]["FDM"],region["CDS"]["All_snp"],region["CDS"]["Avaible_snp"],region["INTRON"]["D"],region["INTRON"]["FDM"],region["INTRON"]["All_snp"],region["INTRON"]["Avaible_snp"],region["REGULATION"]["D"],region["REGULATION"]["FDM"],region["REGULATION"]["All_snp"],region["REGULATION"]["Avaible_snp"]]

D_genes=read_result(window_gff_result)
po1=Pool(20)
results=po1.map(gene_D_Fdm,list(D_genes.keys()))
po1.close()
print("#Gene\tGene_D\tGene_FDM\tGene_All_snp\tGene_Avaible_snp\tCDS_D\tCDS_FDM\tCDS_All_snp\tCDS_Avaible_snp\tINTRON_D\tINTRON_FDM\tINTRON_All_snp\tINTRON_Avaible_snp\tREGULATION_D\tREGULATION_FDM\tREGULATION_All_snp\tREGULATION_Avaible_snp")
for result in results:
    print('\t'.join([str(x) for x in result]))
