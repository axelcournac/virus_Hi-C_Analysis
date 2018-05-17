# Codes and functions for the analysis of virus genome contacts with human chromosomes
### Table of contents

* [Dependencies](https://github.com/axelcournac/virus_HiC_ANALYSIS/edit/master/README.md#dependencies)
* [Raw data extraction and alignment](https://github.com/axelcournac/virus_HiC_ANALYSIS/edit/master/README.md#raw-data-extraction-and-alignment)
* [Building of the contacts map](https://github.com/axelcournac/virus_HiC_ANALYSIS/edit/master/README.md#building-of-the-contacts-map)


### Dependencies

Scripts and codes can be run on OS X and other Unix-based systems. It basically requires to have Python installed on your machine which is commonly installed on Unix-based systems. 
For windows, you can have a look to https://www.python.org/downloads/windows/. Then, a few python modules are necessary for diverses operations on arrays and vizualisation. 

#### Python (>=2.7)
* Numpy
* Matplotlib (>=1.0)
* Scipy
* Biopython

#### External programs

* `Bowtie2 ` / [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* `R` / [R](https://cran.r-project.org/)



## Raw data extraction and alignment
#### Data extraction
Data can be dowloaded on Short Read Archive server at the following address **http://www.ncbi.nlm.nih.gov/sra**.

A SRA executable called fastq-dump from SRA can be used to extract and split both mates of a library (to use it, you can go with your terminal to the directory containg the executables files by using the bash command cd).Then the program can be used like this:  /fastq-dump library_identification --split-3 -O /path_to_a_directory
 
```bash
./fastq-dump SRR639031 --split-3 -O /run/media/axel/human_data/
```

#### Alignment

For the alignement step, we will use the sofware Bowtie2 and an iterative procedure like the one of [DADE] (https://github.com/scovit/dade). 

We process the pairs of reads so that every read has a mapping quality superior to 30. 
Here, some lines that can be used to do this task:

```bash
#  Keeping only the columns of the sam file that contain necessary information:
awk '{print $1,$3,$4,$2,$5;}' p1.sam > p1.sam.0
awk '{print $1,$3,$4,$2,$5;}' p2.sam > p2.sam.0

# Sort according to the read identification to have both mates in the same order
# if sort does not have -V option try -d
sort -V -k1 p1.sam.0 > p1.sam.0.sorted
sort -V -k1 p2.sam.0 > p2.sam.0.sorted

# Pairing of both mates in a single file
paste p1.sam.0.sorted p2.sam.0.sorted > p1_p2_merged

# Removal of intermediar files
rm p1.sam.0.sorted
rm p2.sam.0.sorted

# Filtering of paires of reads that both have a Mapping Quality above 30
awk '{if($1 eq $6 && $5>= 30 && $10 >= 30) print $2,$3,$4,$7,$8,$9}'  p1_p2_merged  > output_alignment_idpt.dat

# Removal of intermediar file
rm p1_p2_merged
```
At this stage, you should have a file containing these information:
```
chr1 104180 16 chr1 104057 0
chr1 3570510 16 chr1 3570450 0
chr1 4255981 0 chr1 4256104 16
chr1 159457 16 chr1 159370 0
chr1 4113710 16 chr1 4113584 0
chr1 4259849 16 chr1 4259818 0
chr1 3753874 0 chr1 3754001 16
chr1 2856270 16 chr1 2856124 0
chr1 4134782 16 chr1 4134678 0
```

chr1 corresponds here to the chromosome of human genome. We used hg19 reference genome. We name this file output_alignment_idpt.dat.ind3.

We then assigned eahc read to its corresponding restriction fragment as described previously in [https://github.com/axelcournac/3C_tutorial](https://github.com/axelcournac/3C_tutorial). 


## Building of the contacts map
To build the contact map and and filtered the non informative events, we use the python code 3Cevents_MATRICE.py [`3Cevents_MATRICE.py`](python_codes/3Cevents_MATRICE.py):
```bash
python Matrice_Creator_ANNOTATIONS10_4_several.py 1 100000 HBVayw HBX_mutant_100kb 0 /media/axel/RSG53/CONCATENATED_PIERRICK_DATA_ITERATIV/output_alignment_idpt_BC78_ACTT_XMutant_fused20runs_cap_BC78_capBC78.indices.filtered.pcr5
```
The first argument corresponds to the path of the file named output_alignment_idpt.dat.ind3 containing the informations of the mapped pairs of reads. 
The second argument is the size of the bin, here 5000 bp. We generally use this resolution for the whole study which is a good compromise between resolution and signa robustness.
The third argument is the name of the prefixe for the file of the contacts maps. 


![alt tag](https://github.com/axelcournac/virus_Hi-C_Analysis/blob/master/pictures/chr3_CM_all_runs_HBVayw_400kb_HBVayw_0.png)



