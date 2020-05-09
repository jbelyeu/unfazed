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
Unfazed also applies an additional phasing technique to copy-number variants (CNVs), by using the allele balance of heterozygous sites are found inside the copy-altered region. 
* In a deletion, the allele of the de novo CNV's origin parent disappears and all sites should be HOM_REF for the other parent's allele. 
* In a duplication, the allele balance of the de novo CNV's origin parent should be about double in proportion to the allele from the other parent.

## How to use it 
Unfazed is available for install from conda. Requires at least Python 3.5.

`conda install -c bioconda unfazed `

<details>
  <summary>Unfazed options:</summary>
  
  ```

UNFAZED v0.2.3
usage: unfazed [-h] [-v] -d DNMS -s SITES -p PED [-b BAM_DIR]
               [--bam-pairs [BAM_PAIRS [BAM_PAIRS ...]]] [-t THREADS]
               [-o {vcf,bed}] [--include-ambiguous] [--verbose]
               [--outfile OUTFILE] [-r REFERENCE]

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Installed version (0.2.3)
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
  -r REFERENCE, --reference REFERENCE
                        reference fasta file (required for crams) (default:
                        None)
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

## Unfazed input and output
Unfazed will accept as input either a valid VCF file of de novo variants or a BED file with fields described below. Output can be either an annotated VCF or a BED file.

### VCF annotations
Unfazed adds two tags to the FORMAT field of the VCF.
* UOPS: support for the UOP call (count of informative sites)
* UET: evidence type(s) for the UOP call, which may be 0:readbacked, 1:allele-balance (for CNVs only), 2:both, 3:ambiguous-A
readbacked, 4:ambiguous-allele-balance, 5:ambiguous-both, -1:missing.

#### Important notes:
* Unfazed replaces `/` with `|` in the GT tag for phased variants to indicate phasing, following the phase order paternal|maternal.

* VCF output is only possible when `--dnms` is a VCF file.


#### VCF lines before annotation with unfazed:

```
22	16256430	.	A	G	453874	PASS	AC=3;AF=0.495;AN=6	GT:AD:DP:GQ:PL	0/1:17,12:29:99:267,0,418	0/1:11,9:20:99:206,0,266	0/1:17,15:32:99:347,0,410
22	16256512	.	T	C	834865	PASS	AC=3;AF=0.494;AN=6  GT:AD:DP:GQ:PL	0/1:10,14:24:99:343,0,224	0/1:13,11:24:99:254,0,312	0/1:26,20:46:99:478,0,660
22	30857373	.	A	C	939637	PASS	AC=3;AF=0.63;AN=6	GT:AD:DP:GQ:PL	0/1:47,47:94:99:1478,0,1355	0/1:48,38:86:99:1166,0,1416	0/1:72,73:145:99:2241,0,1931
22	30857448	.	A	G	815602	PASS	AC=3;AF=0.606;AN=6	GT:AD:DP:GQ:PL	0/1:32,31:63:99:861,0,932	0/1:33,40:73:99:1058,0,876	0/1:60,66:126:99:1904,0,1635
```
#### After:
<pre>
22	16256430	.	A	G	453874	PASS	AC=3;AF=0.495;AN=6	GT:AD:DP:GQ:PL:UOPS:UET	0/1:17,12:29:99:267,0,418:-1:-1	0/1:11,9:20:99:206,0,266:-1:-1	0/1:17,15:32:99:347,0,410:-1:-1
22	16256512	.	T	C	834865	PASS	AC=3;AF=0.494;AN=6	GT:AD:DP:GQ:PL:UOPS:UET	0/1:10,14:24:99:343,0,224:-1:-1	0/1:13,11:24:99:254,0,312:-1:-1	0/1:26,20:46:99:478,0,660:-1:-1
22	30857373	.	A	C	939637	PASS	AC=3;AF=0.63;AN=6	GT:AD:DP:GQ:PL:UOPS:UET	<b>1|0:47,47:94:99:1478,0,1355:2:0</b>	0/1:48,38:86:99:1166,0,1416:-1:-1	0/1:72,73:145:99:2241,0,1931:-1:-1
22	30857448	.	A	G	815602	PASS	AC=3;AF=0.606;AN=6	GT:AD:DP:GQ:PL:UOPS:UET	<b>1|0:32,31:63:99:861,0,932:2:0</b>	0/1:33,40:73:99:1058,0,876:-1:-1	0/1:60,66:126:99:1904,0,1635:-1:-1
</pre>
The GT alleles are separated with `|` for phased variants (where FORMAT fields are __bold__). The added UOPS and UET tags are added to each record, with -1 for variants/samples that were not phased.
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


#### BED input
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
#### BED output (no entries for variants that cannot be phased)
```
#chrom	start	end	vartype	kid	origin_parent	other_parent	evidence_count	evidence_types
22	18844941	18844942	POINT	NA12878	NA12892	NA12891	1	READBACKED
22	21088145	21088146	POINT	NA12878	NA12892	NA12891	1	READBACKED
22	21141299	21141300	POINT	NA12878	NA12892	NA12891	2	READBACKED
22	30857372	30857373	POINT	NA12878	NA12891	NA12892	2	READBACKED
22	30857447	30857448	POINT	NA12878	NA12891	NA12892	2	READBACKED

```

**Ambiguous results** derive from inconsistent phasing (different parent of origin indicated by different informative sites or reads. These may indicate sequencing errors or mosaic events and will *not* be reported unless the `--include-ambiguous` argument is included.

## Performance
Many variants lack informative sites and are therefore can't be phased. Unfazed also makes no attempt to phase multiallelic sites (which should be very rare among de novo calls). Generally a little under 30% of de novo SNVs are phaseable via unfazed, and about 50% of CNV/SV variants.

The runtime of unfazed is highly dependent on the size of the sites VCF, as well as the number of variants. A multithreaded approach is used to improved performance; however, as the performance is bound by file IO, more than 2 threads yield dimenishing returns (and can even cause a slowdown due to race conditions). Running with 2 threads (default option) is therefore recommended. (Expert note: running with 1 thread can often produce more informative error messages in the case of a silent failure)
