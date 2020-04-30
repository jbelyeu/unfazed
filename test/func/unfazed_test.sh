#!/bin/bash

test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
. ssshtest

STOP_ON_FAIL=0
data_path="test/data/"
func_path="test/func/"

bam=$data_path"NA12878.bam"
sites_vcf=$data_path"trio_snvs_chr22.vcf.gz"
snv_hets_vcf=$data_path"trio_hets_snvs_chr22.vcf.gz"
sv_hets_vcf=$data_path"trio_hets_svs_chr22.vcf.gz"
snv_hets_bed=$data_path"trio_hets_snvs_chr22.bed"
sv_hets_bed=$data_path"trio_hets_svs_chr22.bed"
ped=$data_path"trio.ped"
missing_kid_ped=$data_path"trio_missing_kid.ped"
missing_dad_ped=$data_path"trio_missing_dad.ped"


echo "SNV phasing tests"
echo "##########################################################################"

run phase_snv_vcf_to_bed \
    unfazed \
        -d $snv_hets_vcf \
        -s $sites_vcf \
        -p $ped \
        -o bed\
        --bam-pairs "NA12878":$bam
if [ $phase_snv_vcf_to_bed ]; then
    assert_exit_code 0
    assert_in_stdout "#chrom	start	end	vartype	kid	origin_parent	other_parent	evidence_count	evidence_types"
    assert_in_stdout "22	21141299	21141300	POINT	NA12878	NA12892	NA12891	8	READBACKED"
    assert_in_stdout "22	30864691	30864692	POINT	NA12878	NA12891	NA12892	8	READBACKED"
    assert_in_stdout "22	36556963	36556964	POINT	NA12878	NA12892	NA12891	17	READBACKED"
    assert_in_stdout "22	41609689	41609690	POINT	NA12878	NA12892	NA12891	26	READBACKED"
    assert_in_stdout "22	41652845	41652846	POINT	NA12878	NA12892	NA12891	8	READBACKED"
    assert_in_stdout "22	42072911	42072912	POINT	NA12878	NA12891	NA12892	8	READBACKED"
    assert_in_stdout "22	42523635	42523636	POINT	NA12878	NA12891	NA12892	8	READBACKED"
fi

run phase_snv_bed_to_bed \
    unfazed \
        -d $snv_hets_bed \
        -s $sites_vcf \
        -p $ped \
        -o bed\
        --bam-pairs "NA12878":$bam
if [ $phase_snv_bed_to_bed ]; then
    assert_exit_code 0
    assert_in_stdout "#chrom	start	end	vartype	kid	origin_parent	other_parent	evidence_count	evidence_types"
    assert_in_stdout "22	21141299	21141300	POINT	NA12878	NA12892	NA12891	8	READBACKED"
    assert_in_stdout "22	30864691	30864692	POINT	NA12878	NA12891	NA12892	8	READBACKED"
    assert_in_stdout "22	36556963	36556964	POINT	NA12878	NA12892	NA12891	17	READBACKED"
    assert_in_stdout "22	41609689	41609690	POINT	NA12878	NA12892	NA12891	26	READBACKED"
    assert_in_stdout "22	41652845	41652846	POINT	NA12878	NA12892	NA12891	8	READBACKED"
    assert_in_stdout "22	42072911	42072912	POINT	NA12878	NA12891	NA12892	8	READBACKED"
    assert_in_stdout "22	42523635	42523636	POINT	NA12878	NA12891	NA12892	8	READBACKED"
fi

run phase_snv_vcf_to_vcf \
    unfazed \
        -d $snv_hets_vcf \
        -s $sites_vcf \
        -p $ped \
        --bam-pairs "NA12878":$bam
if [ $phase_snv_vcf_to_vcf ]; then
    assert_exit_code 0
    assert_in_stdout '##FORMAT=<ID=UOP,Number=1,Type=Float,Description="Unfazed-identified origin parent. Paternal:`0`, maternal:`1`, missing:`-1`">'
    assert_in_stdout '##FORMAT=<ID=UOPS,Number=1,Type=Float,Description="Count of pieces of evidence supporing the unfazed-identified origin parent or `-1` if missing">'
    assert_in_stdout '##FORMAT=<ID=UET,Number=1,Type=Float,Description="Unfazed evidence type: `0` (readbacked), `1` (non-readbacked, for CNVs only), `2` (both), `3` (ambiguous readbacked), `4` (ambiguous non-readbacked), `5` (ambiguous both) or `-1` (missing)">'
    assert_in_stdout '22	16256512	.	T	C	834865	PASS	AC=3;AF=0.494;AN=6;BaseQRankSum=2.47;ClippingRankSum=-0.031;DP=147650;ExcessHet=2.14748e+09;FS=0.554;InbreedingCoeff=-0.9823;MLEAC=887;MLEAF=0.496;MQ=25.78;MQ0=0;MQRankSum=-0.328;QD=11.42;ReadPosRankSum=0.223;SOR=0.774;VQSLOD=-3.334;culprit=InbreedingCoeff;set=Intersection;het_var=NA12878	GT:AD:DP:GQ:PL:UOP:UOPS:UET	0/1:10,14:24:99:343,0,224:-1:-1:-1	0/1:13,11:24:99:254,0,312:-1:-1:-1	0/1:26,20:46:99:478,0,660:-1:-1:-1'
    assert_in_stdout '22	16266920	.	C	CAG	584855	PASS	AC=3;AF=0.363;AN=6;BaseQRankSum=0.697;ClippingRankSum=-0.113;DP=83342;ExcessHet=2.14748e+09;FS=1.952;InbreedingCoeff=-0.571;MLEAC=650;MLEAF=0.364;MQ=16.98;MQ0=0;MQRankSum=-0.005;QD=8.93;ReadPosRankSum=-0.21;SOR=0.445;VQSLOD=2.41;culprit=FS;set=variant;het_var=NA12878	GT:AD:DP:GQ:PL:UOP:UOPS:UET	0/1:84,51:135:99:1778,0,3267:-1:-1:-1	0/1:108,66:174:99:2334,0,4264:-1:-1:-1	0/1:65,65:130:99:2452,0,2492:-1:-1:-1'
    assert_in_stdout '22	16277852	.	C	T	274318	VQSRTrancheSNP99.90to100.00	AC=3;AF=0.284;AN=6;BaseQRankSum=4.91;ClippingRankSum=-0.03;DP=157044;ExcessHet=2.14748e+09;FS=0.667;InbreedingCoeff=-0.4018;MLEAC=509;MLEAF=0.285;MQ=17.03;MQ0=0;MQRankSum=-1.706;NEGATIVE_TRAIN_SITE;QD=5.23;ReadPosRankSum=2.73;SOR=0.603;VQSLOD=-17.88;culprit=QD;set=Intersection;het_var=NA12878	GT:AD:DP:GQ:PGT:PID:PL:UOP:UOPS:UET	0/1:84,40:124:99:.:.:764,0,1872:-1:-1:-1	0/1:106,43:149:99:.:.:772,0,2527:-1:-1:-1	0/1:154,58:212:99:.:.:1076,0,3816:-1:-1:-1'
    assert_in_stdout '22	18571008	.	G	A	557243	PASS	AC=3;AF=0.667;AN=6;BaseQRankSum=1.17;ClippingRankSum=0;DP=47734;ExcessHet=4.1068;FS=0.612;InbreedingCoeff=-0.015;MLEAC=1184;MLEAF=0.676;MQ=54.25;MQ0=0;MQRankSum=0.077;POSITIVE_TRAIN_SITE;QD=25.94;ReadPosRankSum=0.72;SOR=0.618;VQSLOD=1.91;culprit=FS;set=Intersection;het_var=NA12878	GT:AD:DP:GQ:PGT:PID:PL:UOP:UOPS:UET	0/1:6,13:19:99:.:.:409,0,140:-1:-1:-1	0/1:8,11:19:99:.:.:354,0,232:-1:-1:-1	0/1:7,22:29:99:.:.:747,0,163:-1:-1:-1'
    assert_in_stdout '22	18844942	.	C	T	678344	VQSRTrancheSNP99.90to100.00	AC=3;AF=0.539;AN=6;BaseQRankSum=3.17;ClippingRankSum=0.044;DP=85970;ExcessHet=2.14748e+09;FS=18.41;InbreedingCoeff=-0.7775;MLEAC=895;MLEAF=0.54;MQ=34.81;MQ0=0;MQRankSum=-1.376;NEGATIVE_TRAIN_SITE;QD=16.04;ReadPosRankSum=0.025;SOR=1.79;VQSLOD=-14.27;culprit=InbreedingCoeff;set=Intersection;het_var=NA12878	GT:AD:DP:GQ:PL:UOP:UOPS:UET	0/1:55,26:81:99:530,0,1495:-1:-1:-1	0/1:77,11:88:52:52,0,2094:-1:-1:-1	0/1:35,46:81:99:1202,0,836:-1:-1:-1'
    assert_in_stdout '22	18846088	.	A	G	200380	VQSRTrancheSNP99.90to100.00	AC=3;AF=0.346;AN=6;BaseQRankSum=0.541;ClippingRankSum=0.069;DP=64536;ExcessHet=2.14748e+09;FS=32.721;InbreedingCoeff=-0.5462;MLEAC=632;MLEAF=0.353;MQ=34.18;MQ0=0;MQRankSum=0.893;QD=4.27;ReadPosRankSum=0.535;SOR=2.939;VQSLOD=-37.68;culprit=FS;set=variant2;het_var=NA12878	GT:AD:DP:GQ:PGT:PID:PL:UOP:UOPS:UET	0/1:24,7:31:99:.:.:125,0,666:-1:-1:-1	0/1:38,16:54:99:.:.:346,0,1046:-1:-1:-1	0/1:168,24:192:99:.:.:187,0,4851:-1:-1:-1'
    assert_in_stdout '22	19163577	.	C	CAT	778472	PASS	AC=3;AF=0.596;AN=6;BaseQRankSum=0.679;ClippingRankSum=0.073;DP=61526;ExcessHet=40.4149;FS=0;InbreedingCoeff=-0.1275;MLEAC=1067;MLEAF=0.597;MQ=48.28;MQ0=0;MQRankSum=0.112;QD=29.23;ReadPosRankSum=-0.063;SOR=0.709;VQSLOD=5.5;culprit=FS;set=Intersection;het_var=NA12878	GT:AD:DP:GQ:PGT:PID:PL:UOP:UOPS:UET	0/1:19,18:37:99:.:.:697,0,832:-1:-1:-1	0/1:10,13:23:99:.:.:504,0,428:-1:-1:-1	0/1:26,32:58:99:.:.:1242,0,1098:-1:-1:-1'
    assert_in_stdout '22	19950235	.	C	T	607700	PASS	AC=3;AF=0.498;AN=6;BaseQRankSum=0.441;ClippingRankSum=-0.148;DP=70408;ExcessHet=4.7114;FS=0.547;InbreedingCoeff=-0.0148;MLEAC=889;MLEAF=0.498;MQ=45.52;MQ0=0;MQRankSum=0.058;POSITIVE_TRAIN_SITE;QD=21.89;ReadPosRankSum=0.49;SOR=0.705;VQSLOD=4.74;culprit=FS;set=Intersection;het_var=NA12878	GT:AD:DP:GQ:PGT:PID:PL:UOP:UOPS:UET	0/1:13,17:30:99:.:.:517,0,329:-1:-1:-1	0/1:27,11:38:99:.:.:273,0,803:-1:-1:-1	0/1:38,36:74:99:.:.:969,0,1043:-1:-1:-1'
    #below this line phased, above did not
    assert_in_stdout '22	21141300	.	T	C	636960	PASS	AC=3;AF=0.449;AN=6;BaseQRankSum=0.716;ClippingRankSum=-0.094;DP=82210;ExcessHet=12.0763;FS=0;InbreedingCoeff=-0.0533;MLEAC=804;MLEAF=0.45;MQ=39.81;MQ0=0;MQRankSum=0.034;NEGATIVE_TRAIN_SITE;POSITIVE_TRAIN_SITE;QD=19.83;ReadPosRankSum=0.383;SOR=0.674;VQSLOD=6.33;culprit=FS;set=Intersection;het_var=NA12878	GT:AD:DP:GQ:PGT:PID:PL:UOP:UOPS:UET	0/1:45,58:103:99:.:.:1562,0,1256:1:8:0	0/1:49,62:111:99:.:.:1754,0,1331:-1:-1:-1	0/1:69,49:118:99:.:.:1459,0,2013:-1:-1:-1'
    assert_in_stdout '22	30864692	.	G	A	849895	PASS	AC=3;AF=0.612;AN=6;BaseQRankSum=-0.183;ClippingRankSum=0.168;DP=83652;ExcessHet=3.2909;FS=0.591;InbreedingCoeff=-0.0039;MLEAC=1094;MLEAF=0.612;MQ=40.06;MQ0=0;MQRankSum=0.113;NEGATIVE_TRAIN_SITE;POSITIVE_TRAIN_SITE;QD=23.8;ReadPosRankSum=0.465;SOR=0.584;VQSLOD=2.16;culprit=FS;set=Intersection;het_var=NA12878	GT:AD:DP:GQ:PL:UOP:UOPS:UET	0/1:41,30:71:99:785,0,1319:0:8:0	0/1:59,47:106:99:1253,0,1769:-1:-1:-1	0/1:60,57:117:99:1644,0,1948:-1:-1:-1'
    assert_in_stdout '22	36556964	.	C	T	317075	PASS	AC=3;AF=0.324;AN=6;BaseQRankSum=-0.787;ClippingRankSum=0.027;DP=63460;ExcessHet=0;FS=0.554;InbreedingCoeff=0.1498;MLEAC=581;MLEAF=0.325;MQ=43.33;MQ0=0;MQRankSum=-0.016;POSITIVE_TRAIN_SITE;QD=18.91;ReadPosRankSum=0.429;SOR=0.65;VQSLOD=2.58;culprit=InbreedingCoeff;set=Intersection;het_var=NA12878	GT:AD:DP:GQ:PGT:PID:PL:UOP:UOPS:UET	0/1:7,16:23:99:.:.:323,0,169:1:17:0	0/1:4,11:15:99:.:.:315,0,99:-1:-1:-1	0/1:65,56:121:99:.:.:1335,0,1984:-1:-1:-1'
    assert_in_stdout '22	41609690	.	C	T	280642	PASS	AC=3;AF=0.349;AN=6;BaseQRankSum=-0.391;ClippingRankSum=0.067;DP=48808;ExcessHet=0.2302;FS=0;InbreedingCoeff=0.0433;MLEAC=632;MLEAF=0.354;MQ=52.97;MQ0=0;MQRankSum=0.026;POSITIVE_TRAIN_SITE;QD=19.6;ReadPosRankSum=0.329;SOR=0.729;VQSLOD=6.46;culprit=FS;set=Intersection;het_var=NA12878	GT:AD:DP:GQ:PGT:PID:PL:UOP:UOPS:UET	0/1:12,9:21:99:.:.:231,0,378:1:26:0	0/1:9,9:18:99:.:.:264,0,252:-1:-1:-1	0/1:23,37:60:99:.:.:1070,0,666:-1:-1:-1'
    assert_in_stdout '22	41652846	.	G	T	359359	PASS	AC=3;AF=0.404;AN=6;BaseQRankSum=-0.083;ClippingRankSum=-0.014;DP=56194;ExcessHet=2.4615;FS=0;InbreedingCoeff=0.0034;MLEAC=724;MLEAF=0.405;MQ=51.41;MQ0=0;MQRankSum=0.204;POSITIVE_TRAIN_SITE;QD=20.1;ReadPosRankSum=0.066;SOR=0.67;VQSLOD=5.62;culprit=FS;set=Intersection;het_var=NA12878	GT:AD:DP:GQ:PGT:PID:PL:UOP:UOPS:UET	0/1:12,11:23:99:.:.:322,0,359:1:8:0	0/1:13,12:25:99:.:.:361,0,372:-1:-1:-1	0/1:14,22:36:99:.:.:604,0,417:-1:-1:-1'
    assert_in_stdout '22	42072912	.	C	A	2.4957e+06	PASS	AC=3;AF=0.777;AN=6;BaseQRankSum=0.929;ClippingRankSum=0.052;DP=171198;ExcessHet=0.6194;FS=0;InbreedingCoeff=0.0368;MLEAC=1390;MLEAF=0.777;MQ=26.98;MQ0=0;MQRankSum=-2.388;NEGATIVE_TRAIN_SITE;POSITIVE_TRAIN_SITE;QD=30.04;ReadPosRankSum=0.746;SOR=0.698;VQSLOD=0.983;culprit=MQRankSum;set=Intersection;het_var=NA12878	GT:AD:DP:GQ:PGT:PID:PL:UOP:UOPS:UET	0/1:65,70:135:99:.:.:2030,0,1951:0:8:0	0/1:69,60:129:99:.:.:1697,0,2070:-1:-1:-1	0/1:190,183:373:99:.:.:5271,0,6469:-1:-1:-1'
    assert_in_stdout '22	42523636	.	C	A	45951.7	VQSRTrancheSNP99.90to100.00	AC=3;AF=0.138;AN=6;BaseQRankSum=0.846;ClippingRankSum=-0.118;DP=35891;ExcessHet=147.547;FS=63.584;InbreedingCoeff=-0.1976;MLEAC=270;MLEAF=0.151;MQ=35.41;MQ0=0;MQRankSum=-4.198;QD=3.68;ReadPosRankSum=0.059;SOR=9.899;VQSLOD=-99.4;culprit=FS;set=variant2;het_var=NA12878	GT:AD:DP:GQ:PGT:PID:PL:UOP:UOPS:UET	0/1:42,10:52:99:.:.:185,0,1261:0:8:0	0/1:44,12:56:99:.:.:238,0,1396:-1:-1:-1	0/1:29,5:34:67:.:.:67,0,1001:-1:-1:-1'
fi

run missing_kid_from_ped \
    unfazed \
        -d $snv_hets_vcf \
        -s $sites_vcf \
        -p $missing_kid_ped \
        -o bed\
        --bam-pairs "NA12878":$bam
if [ $missing_kid_from_ped ]; then
    assert_in_stderr "NA12878 missing from pedigree file, will be skipped"
    assert_in_stderr "No phaseable variants"
fi

run missing_dad_from_ped \
    unfazed \
        -d $snv_hets_vcf \
        -s $sites_vcf \
        -p $missing_dad_ped \
        -o bed\
        --bam-pairs "NA12878":$bam
if [ $missing_dad_from_ped ]; then
    assert_in_stderr "Parent of sample NA12878 missing from pedigree file, will be skipped"
    assert_in_stderr "No phaseable variants"
fi


run invalid_output_to_vcf \
    unfazed \
        -d $snv_hets_bed \
        -s $sites_vcf \
        -p $ped \
        -o vcf\
        --bam-pairs "NA12878":$bam
if [ $invalid_output_to_vcf ]; then
    assert_in_stderr "Invalid option: --output-type is vcf, but input is not a vcf type. Rerun with \`--output-type bed\` or input dnms as one of the following: vcf, vcf.gz, bcf"
fi

run invalid_bam \
    unfazed \
        -d $snv_hets_bed \
        -s $sites_vcf \
        -p $ped \
        -o vcf\
        --bam-pairs "NA12878":"bob"
if [ $invalid_bam ]; then
    assert_in_stderr "invalid filename bob"
fi

printf "\n\nSV phasing tests" 
echo "##########################################################################"

run phase_sv_vcf_to_bed \
    unfazed \
        -d $sv_hets_vcf \
        -s $sites_vcf \
        -p $ped \
        -o bed\
        --bam-pairs "NA12878":$bam
if [ $phase_sv_vcf_to_bed ]; then
    assert_exit_code 0
    assert_in_stdout "#chrom	start	end	vartype	kid	origin_parent	other_parent	evidence_count	evidence_types"
    assert_in_stdout "22	42523630	42537254	DEL	NA12878	NA12891	NA12892	35	READBACKED"
fi

run phase_sv_bed_to_bed \
    unfazed \
        -d $sv_hets_bed \
        -s $sites_vcf \
        -p $ped \
        -o bed\
        --bam-pairs "NA12878":$bam
if [ $phase_sv_vcf_to_bed ]; then
    assert_exit_code 0
    assert_in_stdout "#chrom	start	end	vartype	kid	origin_parent	other_parent	evidence_count	evidence_types"
    assert_in_stdout "22	42523630	42537254	DEL	NA12878	NA12891	NA12892	35	READBACKED"
fi

run phase_sv_vcf_to_vcf \
    unfazed \
        -d $sv_hets_vcf \
        -s $sites_vcf \
        -p $ped \
        --bam-pairs "NA12878":$bam 
if [ $phase_sv_vcf_to_vcf ]; then
    assert_exit_code 0
    assert_in_stdout '##FORMAT=<ID=UOP,Number=1,Type=Float,Description="Unfazed-identified origin parent. Paternal:`0`, maternal:`1`, missing:`-1`">'
    assert_in_stdout '##FORMAT=<ID=UOPS,Number=1,Type=Float,Description="Count of pieces of evidence supporing the unfazed-identified origin parent or `-1` if missing">'
    assert_in_stdout '##FORMAT=<ID=UET,Number=1,Type=Float,Description="Unfazed evidence type: `0` (readbacked), `1` (non-readbacked, for CNVs only), `2` (both), `3` (ambiguous readbacked), `4` (ambiguous non-readbacked), `5` (ambiguous both) or `-1` (missing)">'
    assert_in_stdout '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA12878	NA12891	NA12892'
    assert_in_stdout '22	42523631	10839	N	<DEL>	1055.41	.	SVTYPE=DEL;SVLEN=-13623;END=42537254;STRANDS=+-:2863;IMPRECISE;CIPOS=-184,454;CIEND=-474,112;CIPOS95=-78,78;CIEND95=-70,70;SU=2863;PE=2863;SR=0;'
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOP:UOPS:UET	0/1:173:287.68:-31,-2,-20:45:30:14:30:14:0:0:0:30:14:0.32:0:35:0	0/1:149:479.61:-49,-1,-16:51:29:21:29:21:0:0:0:29:21:0.42:-1:-1:-1	0/0:200:0:-0,-26,-88:89:89:0:88:0:56:0:0:32:0:0:-1:-1:-1'
    assert_in_stdout '22	42934847	10843	N	<DEL>	797.65	.	SVTYPE=DEL;SVLEN=-1071;END=42935918;STRANDS=+-:1279;IMPRECISE;CIPOS=-78,429;CIEND=-486,113;CIPOS95=-78,80;CIEND95=-65,65;SU=1279;PE=1279;SR=0;'
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOP:UOPS:UET	0/1:146:221.74:-24,-2,-17:39:26:12:25:11:13:0:2:12:9:0.31:-1:-1:-1	0/0:165:0:-0,-17,-55:56:56:0:55:0:31:0:0:24:0:0:-1:-1:-1	0/1:62:62.74:-12,-5,-27:41:33:7:33:6:19:0:1:14:5:0.15:-1:-1:-1'
fi
