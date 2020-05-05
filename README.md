[![CircleCI](https://circleci.com/gh/jbelyeu/unfazed/tree/master.svg?style=svg)](https://circleci.com/gh/jbelyeu/unfazed/tree/master)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/unfazed/README.html)

# Unfazed: phasing tool for de novo SNVs and SVs
Unfazed identifies the parent of origin for de novo variants, accepting input from either a vcf file or bed file of variant information. Unfazed works for point mutations (SNVs and INDELs) as well as larger structural mutations.

## How it works
### Extended read-backed phasing
Unfazed identifies 'informative sites' 5kb upstream or downstream from a de novo variant, using a VCF/BCF of SNVs for the trio (the child and both parents). These informative sites are variants inherited from the parents that allow identification of the origin of the read (maternal or paternal). 

Informative sites must be HET in the child and discernibly different in parents, specifically HOM_REF|HOM_ALT, HET|HOM_ALT, or HET|HOM_REF. These patterns allow identification of the parent of origin for the allele found at that site in each read spanning the region.

Extended read-backed phasing adds sensitivity by chaining reads together using mutually overlapped heterozygous sites, i.e. if two reads have the same allele for a given het site, those reads must come from the same parent. This allows the search distance from a given de novo variant to extend farther than possible with a single read pair.

### Allele-balance CNV phasing
Copy-number variants allow an additional method of phasing, wherein potential heterozygous sites are found inside the copy-altered region. 
* In a deletion, the allele of the de novo CNV's origin parent should disappear and all sites should be HOM_REF for the other parent's allele. 
* In a duplication, the allele balance of the de novo CNV's origin parent should be about double in proportion to the allele from the other parent.

## How to use it 
Unfazed is available for install from conda. Requires at least Python 3.5.

`conda install unfazed `

<details>
  <summary>Unfazed options:</summary>
  
  ```

usage: unfazed [-h] [-v] -d DNMS -s SITES -p PED [-b BAM_DIR]
               [--bam-pairs [BAM_PAIRS [BAM_PAIRS ...]]] [-t THREADS]
               [-o {vcf,bed}] [--include-ambiguous] [--verbose]
               [--outfile OUTFILE]

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Installed version (0.1.5)
  -d DNMS, --dnms DNMS  valid VCF OR BED file of the DNMs of interest> If BED,
                        must contain chrom, start, end, kid_id, var_type
                        columns (default: None)
  -s SITES, --sites SITES
                        sorted/bgzipped/indexed VCF/BCF file of SNVs to
                        identify informative sites. Must contain each kid and
                        both parents (default: None)
  -p PED, --ped PED     ped file including the kid and both parent IDs
                        (default: None)
  -b BAM_DIR, --bam-dir BAM_DIR
                        directory where bam/cram files (named {sample_id}.bam
                        or {sample_id}.cram) are stored for offspring. If not
                        included, --bam-pairs must be set (default: None)
  --bam-pairs [BAM_PAIRS [BAM_PAIRS ...]]
                        space-delimited list of pairs in the format
                        {sample_id}:{bam_path} where {sample_id} matches an
                        offspring id from the dnm file. Can be used with
                        --bam-dir arg, must be used in its absence (default:
                        None)
  -t THREADS, --threads THREADS
                        number of threads to use (default: 2)
  -o {vcf,bed}, --output-type {vcf,bed}
                        choose output type. If --dnms is not a VCF/BCF, output
                        must be to BED format. Defaults to match --dnms input
                        file (default: None)
  --include-ambiguous   include ambiguous phasing results (default: False)
  --verbose             print verbose output including sites and reads used
                        for phasing. Only applies to BED output (default:
                        False)
  --outfile OUTFILE     name for output file. Defaults to stdout (default:
                        /dev/stdout)
```
</details>

### A simple use case is:

```
unfazed\
  -d mydenovos.vcf.gz\
  -s sites.vcf.gz\
  -p myped.ped\
  --bam-pairs a_sample:a_sample.bam
```
This will print an annotated vcf file of phased variants.

### Unfazed will also accept a bed file as input:

```
unfazed\
  -d mydenovos.bed\
  -s sites.vcf.gz\
  -p myped.ped\
  --bam-pairs a_sample:a_sample.bam
```

This will print a bed file of phased variants. The input bed file must have the following tab-separated columns: chrom, start, end, kid_id, var_type, where var_type is SNV, INDEL, POINT, DEL, DUP, INV, INS, MEI, or BND.

## Interpreting unfazed output
The output options for unfazed are either an annotated version of the input VCF file or a BED file.

### VCF annotations
Unfazed adds three tags to the FORMAT field of the VCF.
* UOP: origin parent, which may be paternal:0, maternal:1, or missing:-1
* UOPS: support for the UOP call (count of informative sites)
* UET: evidence type(s) for the UOP call, which may be 0:readbacked, 1:allele-balance (for CNVs only), 2:both, 3:ambiguous-A
readbacked, 4:ambiguous-allele-balance, 5:ambiguous-both, -1:missing.

VCF output is only possible when given `--dnms` is a VCF file.


<details>
  <summary>VCF lines before/after annotation with unfazed:</summary>
  
#### Before:
```
22	16256430	.	A	G	453874	PASS	AC=3;AF=0.495;AN=6;BaseQRankSum=-0.217;ClippingRankSum=-0.212;DP=82746;ExcessHet=2.14748e+09;FS=0;InbreedingCoeff=-0.9859;MLEAC=888;MLEAF=0.497;MQ=21.85;MQ0=0;MQRankSum=-0.138;NEGATIVE_TRAIN_SITE;QD=11.03;ReadPosRankSum=0.072;SOR=8.989;VQSLOD=-3.497;culprit=InbreedingCoeff;set=Intersection;het_var=NA12878	GT:AD:DP:GQ:PL	0/1:17,12:29:99:267,0,418	0/1:11,9:20:99:206,0,266	0/1:17,15:32:99:347,0,410
22	16256512	.	T	C	834865	PASS	AC=3;AF=0.494;AN=6;BaseQRankSum=2.47;ClippingRankSum=-0.031;DP=147650;ExcessHet=2.14748e+09;FS=0.554;InbreedingCoeff=-0.9823;MLEAC=887;MLEAF=0.496;MQ=25.78;MQ0=0;MQRankSum=-0.328;QD=11.42;ReadPosRankSum=0.223;SOR=0.774;VQSLOD=-3.334;culprit=InbreedingCoeff;set=Intersection;het_var=NA12878	GT:AD:DP:GQ:PL	0/1:10,14:24:99:343,0,224	0/1:13,11:24:99:254,0,312	0/1:26,20:46:99:478,0,660
22	16266920	.	C	CAG	584855	PASS	AC=3;AF=0.363;AN=6;BaseQRankSum=0.697;ClippingRankSum=-0.113;DP=83342;ExcessHet=2.14748e+09;FS=1.952;InbreedingCoeff=-0.571;MLEAC=650;MLEAF=0.364;MQ=16.98;MQ0=0;MQRankSum=-0.005;QD=8.93;ReadPosRankSum=-0.21;SOR=0.445;VQSLOD=2.41;culprit=FS;set=variant;het_var=NA12878	GT:AD:DP:GQ:PL	0/1:84,51:135:99:1778,0,3267	0/1:108,66:174:99:2334,0,4264	0/1:65,65:130:99:2452,0,2492
22	16277852	.	C	T	274318	VQSRTrancheSNP99.90to100.00	AC=3;AF=0.284;AN=6;BaseQRankSum=4.91;ClippingRankSum=-0.03;DP=157044;ExcessHet=2.14748e+09;FS=0.667;InbreedingCoeff=-0.4018;MLEAC=509;MLEAF=0.285;MQ=17.03;MQ0=0;MQRankSum=-1.706;NEGATIVE_TRAIN_SITE;QD=5.23;ReadPosRankSum=2.73;SOR=0.603;VQSLOD=-17.88;culprit=QD;set=Intersection;het_var=NA12878	GT:AD:DP:GQ:PGT:PID:PL	0/1:84,40:124:99:.:.:764,0,1872	0/1:106,43:149:99:.:.:772,0,2527	0/1:154,58:212:99:.:.:1076,0,3816
22	18571008	.	G	A	557243	PASS	AC=3;AF=0.667;AN=6;BaseQRankSum=1.17;ClippingRankSum=0;DP=47734;ExcessHet=4.1068;FS=0.612;InbreedingCoeff=-0.015;MLEAC=1184;MLEAF=0.676;MQ=54.25;MQ0=0;MQRankSum=0.077;POSITIVE_TRAIN_SITE;QD=25.94;ReadPosRankSum=0.72;SOR=0.618;VQSLOD=1.91;culprit=FS;set=Intersection;het_var=NA12878	GT:AD:DP:GQ:PGT:PID:PL	0/1:6,13:19:99:.:.:409,0,140	0/1:8,11:19:99:.:.:354,0,232	0/1:7,22:29:99:.:.:747,0,163
22	18844942	.	C	T	678344	VQSRTrancheSNP99.90to100.00	AC=3;AF=0.539;AN=6;BaseQRankSum=3.17;ClippingRankSum=0.044;DP=85970;ExcessHet=2.14748e+09;FS=18.41;InbreedingCoeff=-0.7775;MLEAC=895;MLEAF=0.54;MQ=34.81;MQ0=0;MQRankSum=-1.376;NEGATIVE_TRAIN_SITE;QD=16.04;ReadPosRankSum=0.025;SOR=1.79;VQSLOD=-14.27;culprit=InbreedingCoeff;set=Intersection;het_var=NA12878	GT:AD:DP:GQ:PL	0/1:55,26:81:99:530,0,1495	0/1:77,11:88:52:52,0,2094	0/1:35,46:81:99:1202,0,836
22	21088146	.	C	G	633192	PASS	AC=3;AF=0.449;AN=6;BaseQRankSum=0.768;ClippingRankSum=-0.065;DP=92196;ExcessHet=13.309;FS=1.195;InbreedingCoeff=-0.0559;MLEAC=804;MLEAF=0.45;MQ=37.08;MQ0=0;MQRankSum=0.358;NEGATIVE_TRAIN_SITE;POSITIVE_TRAIN_SITE;QD=19.29;ReadPosRankSum=0.43;SOR=0.591;VQSLOD=1.74;culprit=InbreedingCoeff;set=Intersection;het_var=NA12878	GT:AD:DP:GQ:PL	0/1:34,38:72:99:991,0,812	0/1:29,26:55:99:758,0,707	0/1:46,57:103:99:1655,0,1073
22	21141300	.	T	C	636960	PASS	AC=3;AF=0.449;AN=6;BaseQRankSum=0.716;ClippingRankSum=-0.094;DP=82210;ExcessHet=12.0763;FS=0;InbreedingCoeff=-0.0533;MLEAC=804;MLEAF=0.45;MQ=39.81;MQ0=0;MQRankSum=0.034;NEGATIVE_TRAIN_SITE;POSITIVE_TRAIN_SITE;QD=19.83;ReadPosRankSum=0.383;SOR=0.674;VQSLOD=6.33;culprit=FS;set=Intersection;het_var=NA12878	GT:AD:DP:GQ:PGT:PID:PL	0/1:45,58:103:99:.:.:1562,0,1256	0/1:49,62:111:99:.:.:1754,0,1331	0/1:69,49:118:99:.:.:1459,0,2013
22	30857373	.	A	C	939637	PASS	AC=3;AF=0.63;AN=6;BaseQRankSum=0.718;ClippingRankSum=0.209;DP=87540;ExcessHet=3.2255;FS=0.567;InbreedingCoeff=-0.0028;MLEAC=1126;MLEAF=0.63;MQ=38.77;MQ0=0;MQRankSum=0.027;NEGATIVE_TRAIN_SITE;POSITIVE_TRAIN_SITE;QD=24.49;ReadPosRankSum=0.463;SOR=0.629;VQSLOD=2.19;culprit=FS;set=Intersection;het_var=NA12878	GT:AD:DP:GQ:PL	0/1:47,47:94:99:1478,0,1355	0/1:48,38:86:99:1166,0,1416	0/1:72,73:145:99:2241,0,1931
22	30857448	.	A	G	815602	PASS	AC=3;AF=0.606;AN=6;BaseQRankSum=0.673;ClippingRankSum=-0.03;DP=82240;ExcessHet=5.3075;FS=0;InbreedingCoeff=-0.0184;MLEAC=1083;MLEAF=0.606;MQ=41.17;MQ0=0;MQRankSum=0.044;NEGATIVE_TRAIN_SITE;POSITIVE_TRAIN_SITE;QD=23.27;ReadPosRankSum=0.243;SOR=0.69;VQSLOD=6.57;culprit=FS;set=Intersection;het_var=NA12878	GT:AD:DP:GQ:PL	0/1:32,31:63:99:861,0,932	0/1:33,40:73:99:1058,0,876	0/1:60,66:126:99:1904,0,1635
```
#### After:
```
22	16256430	.	A	G	453874	PASS	AC=3;AF=0.495;AN=6;BaseQRankSum=-0.217;ClippingRankSum=-0.212;DP=82746;ExcessHet=2.14748e+09;FS=0;InbreedingCoeff=-0.9859;MLEAC=888;MLEAF=0.497;MQ=21.85;MQ0=0;MQRankSum=-0.138;NEGATIVE_TRAIN_SITE;QD=11.03;ReadPosRankSum=0.072;SOR=8.989;VQSLOD=-3.497;culprit=InbreedingCoeff;set=Intersection;het_var=NA12878	GT:AD:DP:GQ:PL:UOP:UOPS:UET	0/1:17,12:29:99:267,0,418:-1:-1:-1	0/1:11,9:20:99:206,0,266:-1:-1:-1	0/1:17,15:32:99:347,0,410:-1:-1:-1
22	16256512	.	T	C	834865	PASS	AC=3;AF=0.494;AN=6;BaseQRankSum=2.47;ClippingRankSum=-0.031;DP=147650;ExcessHet=2.14748e+09;FS=0.554;InbreedingCoeff=-0.9823;MLEAC=887;MLEAF=0.496;MQ=25.78;MQ0=0;MQRankSum=-0.328;QD=11.42;ReadPosRankSum=0.223;SOR=0.774;VQSLOD=-3.334;culprit=InbreedingCoeff;set=Intersection;het_var=NA12878	GT:AD:DP:GQ:PL:UOP:UOPS:UET	0/1:10,14:24:99:343,0,224:-1:-1:-1	0/1:13,11:24:99:254,0,312:-1:-1:-1	0/1:26,20:46:99:478,0,660:-1:-1:-1
22	16266920	.	C	CAG	584855	PASS	AC=3;AF=0.363;AN=6;BaseQRankSum=0.697;ClippingRankSum=-0.113;DP=83342;ExcessHet=2.14748e+09;FS=1.952;InbreedingCoeff=-0.571;MLEAC=650;MLEAF=0.364;MQ=16.98;MQ0=0;MQRankSum=-0.005;QD=8.93;ReadPosRankSum=-0.21;SOR=0.445;VQSLOD=2.41;culprit=FS;set=variant;het_var=NA12878	GT:AD:DP:GQ:PL:UOP:UOPS:UET	0/1:84,51:135:99:1778,0,3267:-1:-1:-1	0/1:108,66:174:99:2334,0,4264:-1:-1:-1	0/1:65,65:130:99:2452,0,2492:-1:-1:-1
22	16277852	.	C	T	274318	VQSRTrancheSNP99.90to100.00	AC=3;AF=0.284;AN=6;BaseQRankSum=4.91;ClippingRankSum=-0.03;DP=157044;ExcessHet=2.14748e+09;FS=0.667;InbreedingCoeff=-0.4018;MLEAC=509;MLEAF=0.285;MQ=17.03;MQ0=0;MQRankSum=-1.706;NEGATIVE_TRAIN_SITE;QD=5.23;ReadPosRankSum=2.73;SOR=0.603;VQSLOD=-17.88;culprit=QD;set=Intersection;het_var=NA12878	GT:AD:DP:GQ:PGT:PID:PL:UOP:UOPS:UET	0/1:84,40:124:99:.:.:764,0,1872:-1:-1:-1	0/1:106,43:149:99:.:.:772,0,2527:-1:-1:-1	0/1:154,58:212:99:.:.:1076,0,3816:-1:-1:-1
22	18571008	.	G	A	557243	PASS	AC=3;AF=0.667;AN=6;BaseQRankSum=1.17;ClippingRankSum=0;DP=47734;ExcessHet=4.1068;FS=0.612;InbreedingCoeff=-0.015;MLEAC=1184;MLEAF=0.676;MQ=54.25;MQ0=0;MQRankSum=0.077;POSITIVE_TRAIN_SITE;QD=25.94;ReadPosRankSum=0.72;SOR=0.618;VQSLOD=1.91;culprit=FS;set=Intersection;het_var=NA12878	GT:AD:DP:GQ:PGT:PID:PL:UOP:UOPS:UET	0/1:6,13:19:99:.:.:409,0,140:-1:-1:-1	0/1:8,11:19:99:.:.:354,0,232:-1:-1:-1	0/1:7,22:29:99:.:.:747,0,163:-1:-1:-1
22	18844942	.	C	T	678344	VQSRTrancheSNP99.90to100.00	AC=3;AF=0.539;AN=6;BaseQRankSum=3.17;ClippingRankSum=0.044;DP=85970;ExcessHet=2.14748e+09;FS=18.41;InbreedingCoeff=-0.7775;MLEAC=895;MLEAF=0.54;MQ=34.81;MQ0=0;MQRankSum=-1.376;NEGATIVE_TRAIN_SITE;QD=16.04;ReadPosRankSum=0.025;SOR=1.79;VQSLOD=-14.27;culprit=InbreedingCoeff;set=Intersection;het_var=NA12878	GT:AD:DP:GQ:PL:UOP:UOPS:UET	0/1:55,26:81:99:530,0,1495:1:1:0	0/1:77,11:88:52:52,0,2094:-1:-1:-1	0/1:35,46:81:99:1202,0,836:-1:-1:-1
22	21088146	.	C	G	633192	PASS	AC=3;AF=0.449;AN=6;BaseQRankSum=0.768;ClippingRankSum=-0.065;DP=92196;ExcessHet=13.309;FS=1.195;InbreedingCoeff=-0.0559;MLEAC=804;MLEAF=0.45;MQ=37.08;MQ0=0;MQRankSum=0.358;NEGATIVE_TRAIN_SITE;POSITIVE_TRAIN_SITE;QD=19.29;ReadPosRankSum=0.43;SOR=0.591;VQSLOD=1.74;culprit=InbreedingCoeff;set=Intersection;het_var=NA12878	GT:AD:DP:GQ:PL:UOP:UOPS:UET	0/1:34,38:72:99:991,0,812:1:1:0	0/1:29,26:55:99:758,0,707:-1:-1:-1	0/1:46,57:103:99:1655,0,1073:-1:-1:-1
22	21141300	.	T	C	636960	PASS	AC=3;AF=0.449;AN=6;BaseQRankSum=0.716;ClippingRankSum=-0.094;DP=82210;ExcessHet=12.0763;FS=0;InbreedingCoeff=-0.0533;MLEAC=804;MLEAF=0.45;MQ=39.81;MQ0=0;MQRankSum=0.034;NEGATIVE_TRAIN_SITE;POSITIVE_TRAIN_SITE;QD=19.83;ReadPosRankSum=0.383;SOR=0.674;VQSLOD=6.33;culprit=FS;set=Intersection;het_var=NA12878	GT:AD:DP:GQ:PGT:PID:PL:UOP:UOPS:UET	0/1:45,58:103:99:.:.:1562,0,1256:1:2:0	0/1:49,62:111:99:.:.:1754,0,1331:-1:-1:-1	0/1:69,49:118:99:.:.:1459,0,2013:-1:-1:-1
22	30857373	.	A	C	939637	PASS	AC=3;AF=0.63;AN=6;BaseQRankSum=0.718;ClippingRankSum=0.209;DP=87540;ExcessHet=3.2255;FS=0.567;InbreedingCoeff=-0.0028;MLEAC=1126;MLEAF=0.63;MQ=38.77;MQ0=0;MQRankSum=0.027;NEGATIVE_TRAIN_SITE;POSITIVE_TRAIN_SITE;QD=24.49;ReadPosRankSum=0.463;SOR=0.629;VQSLOD=2.19;culprit=FS;set=Intersection;het_var=NA12878	GT:AD:DP:GQ:PL:UOP:UOPS:UET	0/1:47,47:94:99:1478,0,1355:0:2:0	0/1:48,38:86:99:1166,0,1416:-1:-1:-1	0/1:72,73:145:99:2241,0,1931:-1:-1:-1
22	30857448	.	A	G	815602	PASS	AC=3;AF=0.606;AN=6;BaseQRankSum=0.673;ClippingRankSum=-0.03;DP=82240;ExcessHet=5.3075;FS=0;InbreedingCoeff=-0.0184;MLEAC=1083;MLEAF=0.606;MQ=41.17;MQ0=0;MQRankSum=0.044;NEGATIVE_TRAIN_SITE;POSITIVE_TRAIN_SITE;QD=23.27;ReadPosRankSum=0.243;SOR=0.69;VQSLOD=6.57;culprit=FS;set=Intersection;het_var=NA12878	GT:AD:DP:GQ:PL:UOP:UOPS:UET	0/1:32,31:63:99:861,0,932:0:2:0	0/1:33,40:73:99:1058,0,876:-1:-1:-1	0/1:60,66:126:99:1904,0,1635:-1:-1:-1
```
UOP, UOPS, and UET tags are added to each record, with -1 for variants/samples that were not phased.
</details>

### BED output
Either VCF or BED input can produce a BED as output, with the following columns: 
```
chrom	start	end	vartype	kid	origin_parent	other_parent	evidence_count	evidence_types
```

Additional information can be included by using the `--verbose` option, which will append the following columns

```
chrom	start	end	vartype	kid	origin_parent	other_parent	evidence_count	evidence_types	origin_parent_sites	origin_parent_reads	other_parent_sites	other_parent_reads
```

Evidence counts and types in BED output match those in VCF output.


<details>
  <summary>BED input/ouput entries:</summary>
  
#### Input
```
#chrom	start	end	kid_id	var_type
22	16256429	16256430	NA12878	POINT
22	16256511	16256512	NA12878	POINT
22	16266919	16266920	NA12878	POINT
22	16277851	16277852	NA12878	POINT
22	18571007	18571008	NA12878	POINT
22	18844941	18844942	NA12878	POINT
22	21088145	21088146	NA12878	POINT
22	21141299	21141300	NA12878	POINT
22	30857372	30857373	NA12878	POINT
22	30857447	30857448	NA12878	POINT

```
#### Output (no entries for variants that cannot be phased)
```
#chrom	start	end	vartype	kid	origin_parent	other_parent	evidence_count	evidence_types
22	18844941	18844942	POINT	NA12878	NA12892	NA12891	1	READBACKED
22	21088145	21088146	POINT	NA12878	NA12892	NA12891	1	READBACKED
22	21141299	21141300	POINT	NA12878	NA12892	NA12891	2	READBACKED
22	30857372	30857373	POINT	NA12878	NA12891	NA12892	2	READBACKED
22	30857447	30857448	POINT	NA12878	NA12891	NA12892	2	READBACKED

```
</details>

**Ambiguous results** derive from inconsistent phasing (different parent of origin indicated by different informative sites or reads. These may indicate sequencing errors or mosaic events and will *not* be reported unless the `--include-ambiguous` argument is included.

## Performance
Many variants lack informative sites and are therefore can't be phased. Unfazed also makes no attempt to phase multiallelic sites (which should be very rare among de novo calls). Generally about 30% of de novo SNV/INDEL variants are phaseable via unfazed, and about 50% of CNV/SV variants. 

The runtime of unfazed is highly dependent on the size of the sites VCF, as well as the number of variants. A multithreaded approach is used to improved performance; however, as the performance is bound by file IO, more than 2 threads yield dimenishing returns (and can even cause a slowdown due to race conditions). Running with 2 threads (default option) is therefore recommended. (Expert note: running with 1 thread can often produce more informative error messages in the case of a silent failure)
