# **CallIntrogressions**

---

## **Description**

---

This script is modified according to https://github.com/malonge/CallIntrogressions/blob/master/get_distances.py, a method used in ：

> Alonge, Michael, et al. "Major impacts of widespread structural variation on gene expression and crop improvement in tomato." Cell 182.1 (2020): 145-161.

## **Usage**

---

```python
usage: get_distances_normal.py [-h] [-m 5] [-o <out_dir>]
                              <SVs.vcf> <group.txt> <SP>
                               <reference.fasta.fai> <100000>


Get Jaccard similarity between one_group and others.

positional arguments:
  <SVs.vcf>             SV vcf file with support vectors.
  <group.txt>           Second column is the phylogenetic group, first column
                        is the accession
  <SP>                  Group to compare to other groups
  <reference.fasta.fai>
                        Fasta index file for the reference genome
  <100000>              Introgression window size.

optional arguments:
  -h, --help            show this help message and exit
  -m 5                  minimum number of SVs needed to calculate Jaccard
  -o <out_dir>          the out dir of saving results

```
