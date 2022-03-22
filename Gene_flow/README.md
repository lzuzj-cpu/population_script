# **Gene_flow**
---
## **Description**
---
`D_Fdm_calculation2.py` is a script for calculating `D` and `fdm`:
Parameter|descripted in 
---|:--
`D`|[Martin et al. 2015](https://doi.org/10.1093/molbev/msu269)
`fdm`|[Malinsky et al. 2015](https://doi.org/10.1126/science.aac9927)

## **Usage**
---
```
usage: D_Fdm_calculation2.py [-h] -o OUT_FILE -v VCF_FILE [-j [THEARDS]] -l
                             group_list --order group_order
                             [-c [{Single,Window,Window_file}]]
                             [-ws [WINDOW_STEP]] [-wf [WINFILE]]

D_Fdm PROCESS

optional arguments:
  -h, --help            show this help message and exit
  -o OUT_FILE, --out OUT_FILE
                        the file where result out
  -v VCF_FILE, --vcf VCF_FILE
                        the vcf file you must be give,gzip or no
  -j [THEARDS], --theards [THEARDS]
                        Please give the number of threads , defult = 1
  -l group_list, --list group_list
                        the group list should be four group,every line is sample_name and group_name
                        example:
                        HH1     classA
                        HH2     classA
                        Hy1     classB
  --order group_order   the four group order,just give only one row, split by',', order from one to fouris like: P1,P2,P3,P4
  -c [{Single,Window,Window_file}], --choice [{Single,Window,Window_file}]
                        choice the calculate type
  -ws [WINDOW_STEP], --window_step [WINDOW_STEP]
                        when calculate = Window,give the window and step, default=100000,100000
  -wf [WINFILE], --winfile [WINFILE]
                        when calculate =Window_file ,you mast give a file like that:
                        chr1    100     10000
                        chr2    233     32435
```

there are three cases the script can calculate

1. just calculate Gene flow at a single c locus：

`python3 D_Fdm_calculation2.py -o single.result -v VCF_FILE -l group_list --order group_order`

2. calculate Gene flow in window, the window start from the first locus locus in each chr :

`python3 D_Fdm_calculation2.py -o window.result -v VCF_FILE -l group_list --order group_order -c Window -ws 10000,2000`

3. you can also create a WINFILE to calculate Gene flow of destination area

`python3 D_Fdm_calculation2.py -o window.result -v VCF_FILE -l group_list --order group_order -c Window_file -wf WINFILE`

---
## **Other**
There are three scripts for calculating gene flow of genes
