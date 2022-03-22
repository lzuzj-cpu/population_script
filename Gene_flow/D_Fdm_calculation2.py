import os,sys,re
import gzip
import argparse
from multiprocessing import Pool
import logging
#*------------------------
#use for calculate D and Fd
#zhangjin lzu 2022/2/25
#-------------------------
def read_windows(file_):
    windows=[]
    with open(file_) as fp:
        for line in fp:
            line=line.strip()
            if not line:
                continue
            line=line.split()
            windows.append([line[0],int(line[1]),int(line[2])])
    return windows

def read_group_list(file_):
    groups={}
    with open(file_) as fp:
        for line in fp:
            line=line.strip()
            if not line:
                continue
            line=line.split()
            groups[line[0]]=line[1]
    return groups

def read_group_order(file_):
    with open(file_) as fp:
        [P1,P2,P3,P4]=fp.readline().strip().split(',')
    return P1,P2,P3,P4

def read_vcf(file_):
    opener=gzip.open if file_.endswith(".gz") else open
    with opener(file_,"rt") as fp:
        for line in fp:
            line=line.strip()
            if re.search("^##",line):
                continue
            yield line 

def get_point_vcf(lines):
    vcf_point_dict={}
    for line in lines:
        line=line.split()
        vcf_point_dict.setdefault(line[0],{}).setdefault(int(line[1]),line[9:])
    return vcf_point_dict

def get_window_windows(window,step,vcf_point_dict):
    windows=[]
    for chr_ in vcf_point_dict:
        max_point=max(vcf_point_dict[chr_].keys())
        min_point=min(vcf_point_dict[chr_].keys())
        start=min_point
        end=start+window-1
        while end<=max_point:
            windows.append([chr_,start,end])
            start=start+step
            if start>=max_point:
                break
            end=end+step
        if start<max_point and end>max_point:
            windows.append([chr_,start,max_point])
    return windows

def get_single_windows(vcf_point_dict):
    windows=[]
    for chr_ in vcf_point_dict:
        for point in vcf_point_dict[chr_]:
            windows.append([chr_,point,point])
    return windows

def check_samples(groups,samples):
    no_in=[]
    for sample in groups:
        if sample not in samples:
            no_in.append(sample)
    if no_in:
        logging.warning("there are some samples not in vcf :%s",no_in)
        sys.exit()

def calculate_D_fdm(window_list):
    [Chr_,PosStart,PosEnd]=window_list
    numerator,denominator,FDnumerator,FDMdenominator=0,0,0,0
    all_snp=0
    avaible_snp=0
    D,FDM="NA","NA"
    if Chr_ not in vcf_point_dict:
        return [Chr_,PosStart,PosEnd,D,FDM,numerator,denominator,FDnumerator,FDMdenominator,all_snp,avaible_snp]
    for i in range(PosStart,PosEnd+1):
        if i not in vcf_point_dict[Chr_]:
            continue
        all_snp+=1
        ind=vcf_point_dict[Chr_][i]
        g={P1:{0:0,1:0},P2:{0:0,1:0},P3:{0:0,1:0},P4:{0:0,1:0}}
        RefAltSupport={P1:{0:{},1:{}},P2:{0:{},1:{}},P3:{0:{},1:{}},P4:{0:{},1:{}}}
        for j in range(0,len(ind)):
            name=samples[j]
            if name in groups:
                group=groups[name]
            else:
                continue
            types=re.split("[/|]",ind[j][0:3])
            for type_ in types:
                if type_=="0":
                    g[group][0]+=1
                    RefAltSupport[group][0][name]=1
                elif type_=="1":
                    g[group][1]+=1
                    RefAltSupport[group][1][name]=1
        failGroupSupport=0
        for group in RefAltSupport:
            #加 每个群体需要至少有一个个体：两个单倍型
            if len(RefAltSupport[group][1].keys())+len(RefAltSupport[group][0].keys())<2:
                failGroupSupport=1
                break
        if failGroupSupport==1:
            continue
        if g[P4][0]==g[P4][1]:
            continue
        DerivedAllele=1 if g[P4][0]>g[P4][1] else 0
        Pi1=g[P1][DerivedAllele]/(g[P1][0]+g[P1][1])
        Pi2=g[P2][DerivedAllele]/(g[P2][0]+g[P2][1])
        Pi3=g[P3][DerivedAllele]/(g[P3][0]+g[P3][1])
        Pi4=g[P4][DerivedAllele]/(g[P4][0]+g[P4][1])
        PiD=Pi2 if Pi2>Pi3 else Pi3
        PiDM=Pi1 if Pi1>Pi3 else Pi3
        Cabbai=(1-Pi1)*Pi2*Pi3*(1-Pi4)
        Cbabai=Pi1*(1-Pi2)*Pi3*(1-Pi4)
        if Pi3==0:
            continue
        numerator+=Cabbai-Cbabai
        denominator+=Cabbai+Cbabai
        FDnumerator+=Cabbai-Cbabai
        if Pi2 >= Pi1:
            FDMdenominator+=(1-Pi1)*PiD*PiD*(1-Pi4) - Pi1*(1-PiD)*PiD*(1-Pi4)
            if Cabbai+Cbabai !=0 and ((1-Pi1)*PiD*PiD*(1-Pi4) - Pi1*(1-PiD)*PiD*(1-Pi4))!=0:
                avaible_snp+=1
        else:
            FDMdenominator-=(1-PiDM)*Pi2*PiDM*(1-Pi4) - PiDM*(1-Pi2)*PiDM*(1-Pi4)
            if Cabbai+Cbabai !=0 and ((1-PiDM)*Pi2*PiDM*(1-Pi4) - PiDM*(1-Pi2)*PiDM*(1-Pi4))!=0:
                avaible_snp+=1
    if denominator!=0:
        D=numerator/denominator
    if FDMdenominator!=0:
        FDM=FDnumerator/FDMdenominator
    return [Chr_,PosStart,PosEnd,D,FDM,numerator,denominator,FDnumerator,FDMdenominator,all_snp,avaible_snp]

def write_result(result_,file_):
    with open(file_,'w') as fp:
        print("#Chrom\tStart\tEnd\tD\tFdm\tnumerator\tdenominator\tFDnumerator\tFDMdenominator\tAll_snp\tAvaible_snp",file=fp)
        for list_ in result_:
            print('\t'.join([str(x) for x in list_]),file=fp)

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,description='D_Fdm PROCESS')
    parser.add_argument('-o', '--out', dest='out_file', required=True, help='the file where result out')
    parser.add_argument('-v', '--vcf', dest='vcf_file', required=True, help='the vcf file you must be give,gzip or no')
    parser.add_argument('-j', '--theards', dest='theards', type=int, nargs='?', default=1, const=1, help='Please give the number of threads , defult = 1')
    parser.add_argument('-l', '--list', dest='group_list', required=True, type=str, metavar='group_list', help="the group list should be four group,every line is sample_name and group_name \nexample:\nHH1\tclassA\nHH2\tclassA\nHy1\tclassB")
    parser.add_argument('--order', dest='group_order', required=True, type=str, metavar='group_order', help="the four group order,just give only one row, split by\',\', order from one to four is like: P1,P2,P3,P4")
    parser.add_argument('-c', '--choice', dest='cal_type', nargs='?', default="Single", const="Single", choices=['Single', 'Window', 'Window_file'],help='choice the calculate type')
    parser.add_argument('-ws', '--window_step', dest='window_step', nargs='?', default="100000,100000", const="100000,100000",help='when calculate = Window,give the window and step, default=100000,100000')
    parser.add_argument('-wf', '--winfile', dest='winfile', nargs='?', default="", const="",help='when calculate =Window_file ,you mast give a file like that:\nchr1\t100\t10000\nchr2\t233\t32435')
    args = parser.parse_args()

    if args.cal_type=='Window':
        [window,step]=args.window_step.split(",")
        window=int(window)
        step=int(step)
    elif args.cal_type=='Window_file':
        windows=read_windows(args.winfile)
    global groups,P1,P2,P3,P4
    groups=read_group_list(args.group_list)
    P1,P2,P3,P4=read_group_order(args.group_order)
    lines=read_vcf(args.vcf_file)
    global samples
    samples=next(lines).split()[9:]
    check_samples(groups,samples)
    global vcf_point_dict
    vcf_point_dict=get_point_vcf(lines)
    if args.cal_type=='Window':
        windows=get_window_windows(window,step,vcf_point_dict)
    elif args.cal_type=='Single':
        windows=get_single_windows(vcf_point_dict)
    po1=Pool(args.theards)
    D_fdm_result=po1.map(calculate_D_fdm,windows)
    po1.close()
    D_fdm_result.sort(key=lambda x:(x[0],x[1]))
    write_result(D_fdm_result,args.out_file)

if __name__=="__main__":
    logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
    main()
