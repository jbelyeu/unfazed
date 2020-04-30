[![CircleCI](https://circleci.com/gh/jbelyeu/unfazed/tree/master.svg?style=svg)](https://circleci.com/gh/jbelyeu/unfazed/tree/master)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/unfazed/README.html)

# unfazed: phasing tool for both SNVs and SVs
Unfazed identifies the parent of origin for de novo variants, accepting input from either a vcf file or bed file of variant information. Unfazed works for point mutations (SNVs and INDELs) as well as larger structural mutations.

## How it works
### Extended read-backed phasing
Unfazed identifies 'informative sites' 5kb upstream or downstream from a de novo variant, using a VCF/BCF of SNVs for the trio (the child and both parents). These informative sites are variants inherited from the parents that allow identification of the origin of the read (maternal or paternal). 

Informative sites must be heterozygous in the child and discernibly different in parents, specifically HOM_REF|HOM_ALT, HET|HOM_ALT, or HET|HOM_REF. These patterns allow identification of the parent of origin for the allele found at that site in each read spanning the region.

Extended read-backed phasing adds sensitivity by chaining reads together using mutually overlapped heterozygous sites, i.e. if two reads have the same allele for a given het site, those reads must come from the same parent. This allows the search distance from a given de novo variant to extend farther than possible with a single read pair.

### Allele-balance CNV phasing
Copy-number variants allow an additional method of phasing, wherein potential heterozygous sites are found inside the copy-altered region. 
* In a deletion, the allele of the de novo CNV's origin parent should disappear and all sites should be HOM_REF for the other parent's allele. 
* In a duplication, the allele balance of the de novo CNV's origin parent should be about double in proportion to the allele from the other parent.

## How to use it 
Unfazed is available for install from conda. You must have or install at least Python 3.5.

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

