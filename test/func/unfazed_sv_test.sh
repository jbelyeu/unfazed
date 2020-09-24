#!/bin/bash

test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
. ssshtest

STOP_ON_FAIL=1
data_path="test/data/"
func_path="test/func/"

bam=$data_path"NA12878.bam"
sites_vcf=$data_path"trio_snvs_chr22.vcf.gz"
sv_hets_vcf=$data_path"trio_hets_svs_chr22.vcf.gz"
sv_hets_bed=$data_path"trio_hets_svs_chr22.bed"
ped=$data_path"trio.ped"
missing_kid_ped=$data_path"trio_missing_kid.ped"
missing_dad_ped=$data_path"trio_missing_dad.ped"


printf "\n\nSV phasing tests" 
echo "##########################################################################"

run phase_sv_vcf_to_bed \
    unfazed \
        -d $sv_hets_vcf \
        -s $sites_vcf \
        -p $ped \
        --verbose 
if [ $phase_sv_vcf_to_bed ]; then
    assert_exit_code 0
fi

