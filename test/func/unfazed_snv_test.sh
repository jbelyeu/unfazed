#!/bin/bash

test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
. ssshtest

STOP_ON_FAIL=1
data_path="test/data/"
func_path="test/func/"

bam=$data_path"NA12878.bam"
sites_vcf=$data_path"trio_snvs_chr22.vcf.gz"
snv_hets_vcf=$data_path"trio_hets_snvs_chr22.vcf.gz"
snv_hets_bed=$data_path"trio_hets_snvs_chr22.bed"
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
        --verbose \
        --bam-pairs "NA12878":$bam
if [ $phase_snv_vcf_to_bed ]; then
    assert_exit_code 0
    assert_in_stdout "#chrom	start	end	vartype	kid	origin_parent	other_parent	evidence_count	evidence_types"
    assert_in_stdout '22	18844941	18844942	POINT	NA12878	NA12892	NA12891	1	READBACKED	18844298'
    assert_in_stdout '22	21088145	21088146	POINT	NA12878	NA12892	NA12891	1	READBACKED	21087760'
    assert_in_stdout '22	21141299	21141300	POINT	NA12878	NA12892	NA12891	2	READBACKED	21141433,21141733'
    assert_in_stdout '22	30857372	30857373	POINT	NA12878	NA12891	NA12892	2	READBACKED	30856412,30856856'
    assert_in_stdout '22	30857447	30857448	POINT	NA12878	NA12891	NA12892	2	READBACKED	30856412,30856856'
    assert_in_stdout '22	30862399	30862400	POINT	NA12878	NA12891	NA12892	4	READBACKED	30861914,30861930,30864609,30864791'
    assert_in_stdout '22	30864691	30864692	POINT	NA12878	NA12891	NA12892	4	READBACKED	30861914,30861930,30864609,30864791'
    assert_in_stdout '22	36556963	36556964	POINT	NA12878	NA12892	NA12891	2	READBACKED	36556697,36556822'
    assert_in_stdout '22	41566594	41566595	POINT	NA12878	NA12892	NA12891	4	READBACKED	41566822,41567163,41567534,41567721'
    assert_in_stdout '22	41609689	41609690	POINT	NA12878	NA12892	NA12891	19	READBACKED	41608461,41608525,41608628,41608669,41608947,41609429,41609442,41609858,41610023,41610241,41610390,41610405,41610448,41610666,41610667,41611062,41611239,41611609,41611618'
    assert_in_stdout '22	41613187	41613188	POINT	NA12878	NA12892	NA12891	3	READBACKED	41612540,41612542,41612964'
    assert_in_stdout '22	41613302	41613303	POINT	NA12878	NA12892	NA12891	3	READBACKED	41612540,41612542,41612964'
    assert_in_stdout '22	41652845	41652846	POINT	NA12878	NA12892	NA12891	2	READBACKED	41652344,41652732'
    assert_in_stdout '22	42072911	42072912	POINT	NA12878	NA12891	NA12892	1	READBACKED	42073133'
    assert_in_stdout '22	42523635	42523636	POINT	NA12878	NA12891	NA12892	3	READBACKED	42523210,42523408,42523527'
    assert_in_stdout '22	50617725	50617726	POINT	NA12878	NA12891	NA12892	1	READBACKED	50617982'
fi

run phase_snv_bed_to_bed \
    unfazed \
        --verbose \
        -d $snv_hets_bed \
        -s $sites_vcf \
        -p $ped \
        -o bed\
        --bam-pairs "NA12878":$bam
if [ $phase_snv_bed_to_bed ]; then
    assert_exit_code 0
    assert_in_stdout "#chrom	start	end	vartype	kid	origin_parent	other_parent	evidence_count	evidence_types"
    assert_in_stdout '22	18844941	18844942	POINT	NA12878	NA12892	NA12891	1	READBACKED	18844298'
    assert_in_stdout '22	21088145	21088146	POINT	NA12878	NA12892	NA12891	1	READBACKED	21087760'
    assert_in_stdout '22	21141299	21141300	POINT	NA12878	NA12892	NA12891	2	READBACKED	21141433,21141733'
    assert_in_stdout '22	30857372	30857373	POINT	NA12878	NA12891	NA12892	2	READBACKED	30856412,30856856'
    assert_in_stdout '22	30857447	30857448	POINT	NA12878	NA12891	NA12892	2	READBACKED	30856412,30856856'
    assert_in_stdout '22	30862399	30862400	POINT	NA12878	NA12891	NA12892	4	READBACKED	30861914,30861930,30864609,30864791'
    assert_in_stdout '22	30864691	30864692	POINT	NA12878	NA12891	NA12892	4	READBACKED	30861914,30861930,30864609,30864791'
    assert_in_stdout '22	36556963	36556964	POINT	NA12878	NA12892	NA12891	2	READBACKED	36556697,36556822'
    assert_in_stdout '22	41566594	41566595	POINT	NA12878	NA12892	NA12891	4	READBACKED	41566822,41567163,41567534,41567721'
    assert_in_stdout '22	41609689	41609690	POINT	NA12878	NA12892	NA12891	19	READBACKED	41608461,41608525,41608628,41608669,41608947,41609429,41609442,41609858,41610023,41610241,41610390,41610405,41610448,41610666,41610667,41611062,41611239,41611609,41611618'
    assert_in_stdout '22	41613187	41613188	POINT	NA12878	NA12892	NA12891	3	READBACKED	41612540,41612542,41612964'
    assert_in_stdout '22	41613302	41613303	POINT	NA12878	NA12892	NA12891	3	READBACKED	41612540,41612542,41612964'
    assert_in_stdout '22	41652845	41652846	POINT	NA12878	NA12892	NA12891	2	READBACKED	41652344,41652732'
    assert_in_stdout '22	42072911	42072912	POINT	NA12878	NA12891	NA12892	1	READBACKED	42073133'
    assert_in_stdout '22	42523635	42523636	POINT	NA12878	NA12891	NA12892	3	READBACKED	42523210,42523408,42523527'
    assert_in_stdout '22	50617725	50617726	POINT	NA12878	NA12891	NA12892	1	READBACKED	50617982'
fi

#split the next two tests to decrease the amount of test to stdout
rm -f out.tmp
unfazed \
    -d $snv_hets_vcf \
    -s $sites_vcf \
    -p $ped \
    --outfile out.tmp\
    --bam-pairs "NA12878":$bam

if [ "$(uname)" == "Darwin" ]; then
    grep_header=(grep -E 'FORMAT=<ID=(UOP|UOPS|UET)' out.tmp )
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    grep_header=(grep -E 'FORMAT=<ID=(UOP|UOPS|UET)' out.tmp )
fi

run phase_snv_vcf_to_vcf_header \
    "${grep_header[@]}"
if [ $phase_snv_vcf_to_vcf_header ]; then
    assert_exit_code 0
    assert_in_stdout '##FORMAT=<ID=UOPS,Number=1,Type=Float,Description="Count of pieces of evidence supporting the unfazed-identified origin parent or `-1` if missing">'
    assert_in_stdout '##FORMAT=<ID=UET,Number=1,Type=Float,Description="Unfazed evidence type: `0` (readbacked), `1` (allele-balance, for CNVs only), `2` (both), `3` (ambiguous readbacked), `4` (ambiguous allele-balance), `5` (ambiguous both) or `-1` (missing)">'
fi

run phase_snv_vcf_to_vcf_body \
    grep -v "#" out.tmp | grep "UOPS" | cut -f 9-12
if [ $phase_snv_vcf_to_vcf_body ]; then
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:200:650.44:-67,-2,-25:75:45:30:44:29:27:11:4:17:13:0.4:-1:-1	0/1:165:239.71:-26,-2,-19:43:29:13:28:12:17:3:2:11:6:0.3:-1:-1	0/1:200:452.56:-47,-2,-23:61:39:22:38:21:24:9:1:14:10:0.36:-1:-1'	#16256430
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0|1:100:137.83:-15,-2,-12:25:17:7:17:7:17:5:1:0:0:0.29:9:0	0/0:81:0:-0,-8,-27:27:27:0:27:0:27:0:0:0:0:0:-1:-1	0/1:74:74.77:-12,-4,-23:36:29:6:29:6:29:5:0:0:0:0.17:-1:-1'	#16256512
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0|1:71:305.78:-31,-1,-8:29:15:13:15:13:15:11:1:0:0:0.46:9:0	0/0:78:0:-0,-8,-26:26:26:0:26:0:26:0:0:0:0:0:-1:-1	0/1:123:257.75:-27,-1,-14:35:22:12:22:12:22:10:1:0:0:0.35:-1:-1'	#16266920
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:66:152.86:-16,-1,-8:20:12:7:12:7:12:5:1:0:0:0.37:-1:-1	0/0:75:0:-0,-8,-25:25:25:0:25:0:25:0:0:0:0:0:-1:-1	0/1:98:98.81:-13,-3,-16:28:21:6:21:6:21:5:0:0:0:0.22:-1:-1'	#16277852
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:124:386.68:-40,-1,-14:42:24:17:24:17:24:17:0:0:0:0.41:-1:-1	0/1:99:293.75:-30,-1,-11:33:19:13:19:13:19:13:0:0:0:0.41:-1:-1	0/1:113:113.78:-14,-3,-19:33:25:7:25:7:25:7:0:0:0:0.22:-1:-1'	#18571008
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:124:386.68:-40,-1,-14:42:24:17:24:17:24:17:0:0:0:0.41:-1:-1	0/1:99:293.75:-30,-1,-11:33:19:13:19:13:19:13:0:0:0:0.41:-1:-1	0/1:113:113.78:-14,-3,-19:33:25:7:25:7:25:7:0:0:0:0.22:-1:-1'	#18844942
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0|1:200:281.6:-32,-4,-28:59:42:16:41:15:25:7:0:16:7:0.27:3:0	0/0:144:0:-0,-14,-48:49:49:0:48:0:26:0:0:22:0:0:-1:-1	1/1:117:1358.95:-138,-14,-2:49:1:47:0:46:0:26:1:0:18:1:-1:-1'	#18846088
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:116:182.79:-20,-2,-13:30:20:9:20:9:0:0:0:20:9:0.31:-1:-1	0/1:89:89.86:-11,-2,-11:21:15:5:15:5:0:0:0:15:5:0.25:-1:-1	0/1:59:59.8:-10,-4,-20:31:25:5:25:5:0:0:0:25:5:0.17:-1:-1'	#19163577
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:147:350.67:-37,-2,-16:44:27:16:27:16:0:0:0:27:16:0.37:-1:-1	0/1:141:275.72:-29,-2,-16:39:25:13:25:13:0:0:0:25:13:0.34:-1:-1	0/1:139:302.71:-32,-2,-15:40:25:14:25:14:0:0:0:25:14:0.36:-1:-1'	#19950235
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:147:350.67:-37,-2,-16:44:27:16:27:16:0:0:0:27:16:0.37:-1:-1	0/1:134:278.72:-29,-2,-15:38:24:13:24:13:0:0:0:24:13:0.35:-1:-1	0/1:132:305.71:-32,-1,-15:39:24:14:24:14:0:0:0:24:14:0.37:-1:-1'	#19951207
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:125:125.74:-16,-4,-23:39:30:8:30:8:0:0:0:30:8:0.21:-1:-1	0/1:128:128.74:-16,-4,-22:38:29:8:29:8:0:0:0:29:8:0.22:-1:-1	0/1:158:158.74:-19,-3,-20:38:28:9:28:9:0:0:0:28:9:0.24:-1:-1'	#19951271
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:89:89.54:-18,-9,-49:72:61:11:60:10:32:0:2:28:8:0.14:-1:-1	0/1:74:74.57:-17,-9,-46:67:56:10:56:9:33:0:1:23:8:0.14:-1:-1	0/1:158:164.75:-19,-3,-19:36:26:9:26:9:0:0:0:26:9:0.26:-1:-1'	#20937601
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:89:116.85:-13,-2,-11:22:15:6:15:6:15:3:2:0:0:0.29:-1:-1	1/1:53:620.39:-63,-6,-1:21:0:21:0:21:0:15:5:0:0:1:-1:-1	0/1:29:29.87:-6,-3,-14:21:17:3:17:3:17:2:0:0:0:0.15:-1:-1'	#20939123
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0|1:200:737.42:-75,-1,-23:77:43:33:42:32:25:14:2:17:15:0.43:11:0	0/0:200:0:-0,-29,-96:97:97:0:96:0:61:0:0:35:0:0:-1:-1	0/1:200:452.49:-48,-3,-30:72:48:23:47:22:30:7:4:17:10:0.32:-1:-1'	#21073122
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:86:86.79:-12,-3,-19:34:26:7:25:6:15:3:0:10:2:0.19:-1:-1	0/0:200:0:-0,-25,-83:84:84:0:83:0:49:0:0:34:0:0:-1:-1	1/1:97:1122.61:-114,-11,-2:39:0:39:0:38:0:14:5:0:18:1:-1:-1'	#21088146
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0|1:135:200.76:-22,-2,-15:34:23:10:23:10:23:7:3:0:0:0.3:30:0	0/0:108:0:-0,-11,-36:36:36:0:36:0:36:0:0:0:0:0:-1:-1	0/1:70:176.84:-19,-1,-8:21:13:8:13:8:13:3:4:0:0:0.38:-1:-1'	#21141300
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:41:41.76:-9,-5,-26:37:31:5:31:5:31:5:0:0:0:0.14:-1:-1	0/0:153:0:-0,-15,-51:52:52:0:51:0:26:0:0:25:0:0:-1:-1	0/0:75:0:-8,-15,-63:78:71:6:70:5:44:4:0:26:1:0.067:-1:-1'	#21167787
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:98:98.75:-14,-4,-23:38:30:7:30:7:0:0:0:30:7:0.19:-1:-1	0/1:44:44.77:-10,-5,-25:36:30:5:30:5:0:0:0:30:5:0.14:-1:-1	0/0:200:0:-0,-23,-76:78:77:0:76:0:47:0:0:29:0:0:-1:-1'	#22664216
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:47:47.84:-8,-3,-16:25:20:4:20:4:0:0:0:20:4:0.17:-1:-1	0/1:73:356.75:-37,-1,-8:32:16:15:16:15:0:0:0:16:15:0.48:-1:-1	0/1:135:200.76:-22,-2,-15:34:23:10:23:10:0:0:0:23:10:0.3:-1:-1'	#22868493
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:38:38.38:-6,-2,-6:82:72:9:71:8:40:0:4:31:3:0.1:-1:-1	0/0:69:0:-0,-7,-13:76:76:0:75:0:46:0:0:29:0:0:-1:-1	0/1:32:86.72:-10,-1,-5:90:76:13:75:12:45:0:4:30:8:0.14:-1:-1'	#22868497
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0|1:176:176.71:-21,-3,-22:44:32:11:31:10:17:1:1:14:7:0.24:3:0	0/1:154:347.67:-36,-2,-17:47:29:17:28:16:19:7:1:9:8:0.36:-1:-1	1/1:66:768.1:-78,-8,-1:27:0:27:0:26:0:7:2:0:16:1:-1:-1'	#22868773
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0|1:158:242.72:-26,-2,-18:42:28:13:27:12:14:3:1:13:7:0.31:3:0	0/1:127:566.59:-58,-1,-14:54:28:25:27:24:18:14:1:9:8:0.47:-1:-1	1/1:89:1033.98:-105,-11,-2:36:0:36:0:35:0:15:3:0:16:1:-1:-1'	#22868776
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0|1:143:326.69:-34,-2,-16:44:27:16:26:15:10:0:2:16:13:0.37:2:0	0/0:24:0.02:-5,-8,-31:39:35:3:35:3:35:0:2:0:0:0.079:-1:-1	1/1:58:679.48:-69,-7,-1:25:0:24:0:23:0:0:3:0:19:1:-1:-1'	#22869123
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:182:182.66:-23,-4,-28:50:38:11:38:11:38:9:1:0:0:0.22:-1:-1	0/0:180:0:-1,-19,-67:73:70:2:69:1:43:0:0:26:1:0.014:-1:-1	0/0:66:0:-3,-10,-37:43:40:2:40:2:40:1:0:0:0:0.048:-1:-1'	#22869209
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:122:122.8:-15,-3,-16:30:22:7:22:7:0:0:0:22:7:0.24:-1:-1	0/1:35:35.82:-8,-4,-20:29:24:4:24:4:0:0:0:24:4:0.14:-1:-1	0/1:98:98.81:-13,-3,-16:28:21:6:21:6:0:0:0:21:6:0.22:-1:-1'	#22869218
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:119:119.73:-16,-4,-24:41:32:8:32:8:32:7:0:0:0:0.2:-1:-1	0/1:59:59.67:-13,-7,-35:51:43:7:43:7:43:5:1:0:0:0.14:-1:-1	0/1:200:221.67:-25,-3,-24:47:34:12:34:12:34:9:2:0:0:0.26:-1:-1'	#22869538
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:200:515.49:-54,-2,-27:71:45:25:44:24:27:13:1:17:9:0.35:-1:-1	0/1:200:425.44:-47,-4,-38:81:57:23:56:22:37:10:2:19:9:0.28:-1:-1	0/1:159:656.52:-67,-1,-17:64:34:29:33:28:20:16:2:13:9:0.46:-1:-1'	#22869545
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	1|0:90:245.79:-26,-1,-10:29:17:11:17:11:17:7:3:0:0:0.39:2:0	1/1:56:649.93:-66,-7,-1:22:0:22:0:22:0:19:2:0:0:1:-1:-1	0/0:135:0:-0,-14,-45:45:45:0:45:0:45:0:0:0:0:0:-1:-1'	#22869548
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:29:29.81:-8,-5,-22:31:26:4:26:4:0:0:0:26:4:0.13:-1:-1	0/0:36:0:-2,-5,-20:23:21:1:21:1:0:0:0:21:1:0.045:-1:-1	0/1:86:86.79:-12,-3,-19:32:25:6:25:6:0:0:0:25:6:0.19:-1:-1'	#22869649
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:200:215.6:-27,-5,-33:61:46:14:45:13:26:5:0:19:8:0.22:-1:-1	0/1:200:215.53:-28,-6,-41:71:55:15:54:14:32:7:0:22:7:0.21:-1:-1	0/0:200:0:-0,-25,-82:83:83:0:82:0:51:0:0:31:0:0:-1:-1'	#22869742
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:200:356.37:-43,-7,-51:94:71:22:70:21:40:11:3:30:6:0.23:-1:-1	0/1:71:71.38:-21,-14,-70:99:85:13:84:12:51:5:3:33:3:0.12:-1:-1	0/0:200:0:-0,-22,-73:74:74:0:73:0:45:0:0:28:0:0:-1:-1'	#23438191
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:41:41.57:-14,-10,-49:69:59:9:58:8:34:1:0:24:7:0.12:-1:-1	0/1:131:131.56:-21,-8,-43:69:56:12:55:11:33:2:0:22:9:0.17:-1:-1	0/0:200:0:-0,-23,-77:78:78:0:77:0:49:0:0:28:0:0:-1:-1'	#23503121
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:133:434.65:-45,-1,-14:46:26:19:26:19:26:16:2:0:0:0.42:-1:-1	0/1:116:260.75:-27,-1,-13:34:21:12:21:12:21:9:2:0:0:0.36:-1:-1	0/1:118:311.72:-32,-1,-13:37:22:14:22:14:22:12:1:0:0:0.39:-1:-1'	#23503170
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:200:491.51:-51,-2,-27:68:44:24:43:23:27:7:2:16:13:0.35:-1:-1	0/0:192:0:-0,-19,-64:65:65:0:64:0:43:0:0:21:0:0:-1:-1	0/1:200:344.6:-37,-3,-25:58:39:18:38:17:22:8:2:16:6:0.31:-1:-1'	#23503205
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:29:29.81:-8,-5,-22:31:26:4:26:4:0:0:0:26:4:0.13:-1:-1	0/0:36:0:-2,-5,-20:23:21:1:21:1:0:0:0:21:1:0.045:-1:-1	0/1:86:86.79:-12,-3,-19:32:25:6:25:6:0:0:0:25:6:0.19:-1:-1'	#23503756
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:93:218.8:-23,-1,-11:28:17:10:17:10:0:0:0:17:10:0.37:-1:-1	0/1:133:149.79:-17,-2,-16:31:22:8:22:8:0:0:0:22:8:0.27:-1:-1	0/1:104:104.83:-13,-2,-14:26:19:6:19:6:0:0:0:19:6:0.24:-1:-1'	#23708474
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:200:212.53:-28,-7,-41:72:56:15:55:14:34:8:0:21:6:0.2:-1:-1	0/1:116:116.59:-19,-7,-41:64:52:11:51:10:32:4:0:19:6:0.16:-1:-1	0/1:182:182.6:-24,-6,-36:62:48:13:47:12:27:4:0:20:8:0.2:-1:-1'	#23712520
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0|1:170:392.63:-41,-2,-19:51:32:19:31:18:17:7:3:14:7:0.37:2:0	0/1:142:482.62:-49,-1,-15:51:29:22:28:21:15:12:0:13:8:0.43:-1:-1	1/1:127:1477.12:-150,-15,-2:53:1:51:0:50:0:19:5:0:25:1:-1:-1'	#23712647
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	1|0:156:191.74:-22,-2,-18:37:26:10:26:10:26:6:3:0:0:0.28:1:0	1/1:71:827.19:-84,-8,-1:28:0:28:0:28:0:22:5:0:0:1:-1:-1	0/0:153:0:-0,-15,-51:51:51:0:51:0:51:0:0:0:0:0:-1:-1'	#24579049
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0|1:200:470.53:-49,-2,-25:65:42:23:41:22:25:11:1:16:9:0.35:2:0	0/1:200:326.63:-35,-2,-23:54:36:17:35:16:23:5:1:12:10:0.31:-1:-1	1/1:125:1447.58:-147,-15,-2:50:0:50:0:49:0:18:6:0:24:1:-1:-1'	#24656854
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:74:200.82:-21,-1,-8:24:14:9:14:9:14:4:4:0:0:0.39:-1:-1	0/1:72:227.81:-24,-1,-8:25:14:10:14:10:14:6:3:0:0:0.42:-1:-1	0/1:105:161.81:-18,-2,-12:27:18:8:18:8:18:2:5:0:0:0.31:-1:-1'	#26709766
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0|1:111:236.77:-25,-1,-13:32:20:11:20:11:0:0:0:20:11:0.35:5:0	0/0:200:0:-0,-23,-77:78:78:0:77:0:49:0:0:28:0:0:-1:-1	0/1:56:56.79:-10,-4,-21:32:26:5:26:5:0:0:0:26:5:0.16:-1:-1'	#30184950
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0|1:167:167.63:-22,-5,-33:57:44:12:43:11:24:3:1:19:6:0.2:2:0	0/0:159:0:-0,-16,-53:55:54:0:53:0:30:0:0:23:0:0:-1:-1	1/1:61:709.02:-72,-7,-1:25:0:25:0:24:0:5:3:0:16:1:-1:-1'	#30202774
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0|1:104:104.63:-17,-7,-37:57:47:10:46:9:27:1:1:19:6:0.16:2:0	0/0:132:0:-0,-13,-44:45:44:0:44:0:24:0:0:20:0:0:-1:-1	1/1:53:620.39:-63,-6,-1:22:0:22:0:21:0:2:3:0:16:1:-1:-1'	#30384572
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:200:248.6:-29,-4,-31:60:44:15:43:14:25:0:2:18:12:0.25:-1:-1	0/0:150:0:-0,-15,-50:53:51:1:50:0:29:0:0:21:0:0:-1:-1	0/1:25:299.84:-31,-1,-3:23:9:13:8:12:7:0:4:1:8:0.6:-1:-1'	#30700607
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:68:68.82:-10,-3,-17:28:22:5:22:5:0:0:0:22:5:0.19:-1:-1	0/1:59:59.86:-8,-2,-12:21:16:4:16:4:0:0:0:16:4:0.2:-1:-1	0/0:189:0:-0,-19,-63:64:64:0:63:0:38:0:0:25:0:0:-1:-1'	#30733111
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:89:89.67:-15,-6,-34:53:43:9:42:8:23:0:1:19:7:0.16:-1:-1	0/1:59:59.8:-10,-4,-20:33:26:6:25:5:15:0:1:10:4:0.17:-1:-1	0/0:183:0:-0,-18,-61:62:62:0:61:0:34:0:0:27:0:0:-1:-1'	#30762140
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	1|0:141:275.72:-29,-2,-16:39:25:13:25:13:0:0:0:25:13:0.34:2:0	0/0:66:0:-2,-8,-30:33:31:1:31:1:0:0:0:31:1:0.031:-1:-1	0/1:140:140.77:-17,-3,-18:34:25:8:25:8:0:0:0:25:8:0.24:-1:-1'	#30765502
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:173:287.68:-31,-2,-20:45:30:14:30:14:0:0:0:30:14:0.32:-1:-1	0/1:149:479.61:-49,-1,-16:51:29:21:29:21:0:0:0:29:21:0.42:-1:-1	0/0:200:0:-0,-26,-88:89:89:0:88:0:56:0:0:32:0:0:-1:-1'	#30767729
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:146:221.74:-24,-2,-17:39:26:12:25:11:13:0:2:12:9:0.31:-1:-1	0/0:165:0:-0,-17,-55:56:56:0:55:0:31:0:0:24:0:0:-1:-1	0/1:62:62.74:-12,-5,-27:41:33:7:33:6:19:0:1:14:5:0.15:-1:-1'	#30771458
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0|1:200:248.66:-28,-3,-23:49:35:14:34:13:21:2:1:13:9:0.28:3:2	0/0:198:0:-0,-20,-66:67:67:0:66:0:40:0:0:26:0:0:-1:-1	1/1:64:817.19:-83,-7,-1:32:2:29:1:28:1:3:5:0:19:0.97:-1:-1'	#30776095
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0|1:200:248.66:-28,-3,-23:49:35:14:34:13:21:2:1:13:9:0.28:3:2	0/0:200:0:-0,-20,-67:68:68:0:67:0:41:0:0:26:0:0:-1:-1	1/1:64:817.19:-83,-7,-1:32:2:29:1:28:1:3:5:0:19:0.97:-1:-1'	#30857373
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:200:239.65:-27,-4,-26:52:37:14:37:13:23:2:1:14:9:0.26:-1:-1	0/0:195:0:-0,-20,-65:66:66:0:65:0:40:0:0:25:0:0:-1:-1	1/1:56:649.93:-66,-7,-1:23:0:23:0:22:0:1:1:0:19:1:-1:-1'	#30857448
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:60:206.84:-22,-1,-7:22:12:9:12:9:0:0:0:12:9:0.43:-1:-1	0/0:192:0:-0,-19,-64:65:65:0:64:0:39:0:0:25:0:0:-1:-1	1/1:48:561.31:-57,-6,-1:19:0:19:0:19:0:0:0:0:19:1:-1:-1'	#30862400
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:95:269.77:-28,-1,-11:32:18:13:18:12:10:0:3:8:9:0.4:-1:-1	0/1:50:50.78:-10,-5,-23:36:29:6:28:5:16:0:0:12:5:0.15:-1:-1	1/1:51:590.85:-60,-6,-1:21:0:21:0:20:0:0:4:0:16:1:-1:-1'	#30864692
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:106:290.75:-30,-1,-12:36:21:14:20:13:12:0:4:8:9:0.39:-1:-1	0/1:50:50.72:-11,-6,-30:46:38:7:37:6:23:0:1:14:5:0.14:-1:-1	0/1:25:584.7:-60,-2,-4:38:13:24:12:23:10:0:7:2:16:0.66:-1:-1'	#31838085
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:71:71.89:-9,-2,-9:17:12:4:12:4:0:0:0:12:4:0.25:-1:-1	0/0:81:0:-0,-8,-27:27:27:0:27:0:14:0:0:13:0:0:-1:-1	0/1:33:140.9:-15,-1,-4:14:7:6:7:6:0:0:0:7:6:0.46:-1:-1'	#32756703
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0|1:200:404.59:-43,-2,-22:58:37:20:36:19:21:3:2:15:13:0.35:66:0	0/0:189:0:-0,-19,-63:64:64:0:63:0:39:0:0:24:0:0:-1:-1	0/1:200:251.61:-29,-4,-30:59:43:15:42:14:26:5:3:16:6:0.25:-1:-1'	#32756744
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:38:38.89:-6,-2,-11:18:14:3:14:3:14:1:1:0:0:0.18:-1:-1	0/1:62:257.81:-27,-1,-7:25:13:11:13:11:13:7:3:0:0:0.46:-1:-1	0/1:77:173.83:-19,-1,-9:23:14:8:14:8:14:5:2:0:0:0.36:-1:-1'	#35463162
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:53:53.85:-8,-3,-14:22:18:4:18:4:0:0:0:18:4:0.18:-1:-1	0/1:135:200.76:-22,-2,-15:34:23:10:23:10:0:0:0:23:10:0.3:-1:-1	0/0:150:0:-0,-15,-50:50:50:0:50:0:29:0:0:21:0:0:-1:-1'	#35463179
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:26:26.87:-6,-3,-15:22:18:3:18:3:0:0:0:18:3:0.14:-1:-1	0/0:39:0:-2,-6,-21:24:22:1:22:1:0:0:0:22:1:0.043:-1:-1	0/1:21:197.91:-21,-1,-3:15:6:8:6:8:0:0:0:6:8:0.57:-1:-1'	#35742925
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:48:263.82:-27,-1,-6:23:11:11:11:11:0:0:0:11:11:0.5:-1:-1	0/1:75:122.86:-14,-1,-9:20:13:6:13:6:0:0:0:13:6:0.32:-1:-1	0/1:70:176.84:-19,-1,-8:21:13:8:13:8:0:0:0:13:8:0.38:-1:-1'	#35743124
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	1|0:200:560.45:-58,-2,-29:76:48:27:47:26:28:14:1:19:10:0.36:3:0	1/1:127:1477.12:-150,-15,-2:51:0:51:0:50:0:21:8:0:20:1:-1:-1	0/1:161:422.63:-44,-1,-18:52:31:20:30:19:19:9:2:11:7:0.39:-1:-1'	#36556964
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:200:587.44:-61,-2,-28:77:48:28:47:27:24:17:2:23:8:0.36:-1:-1	0/1:200:755.39:-77,-1,-24:81:46:34:45:33:23:22:1:22:10:0.42:-1:-1	0/1:200:791.34:-81,-2,-28:89:52:36:51:35:27:21:4:24:9:0.41:-1:-1'	#36710183
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:200:407.53:-44,-3,-29:67:45:21:44:20:24:8:3:20:8:0.31:-1:-1	0/1:200:230.63:-27,-4,-29:56:41:14:40:13:23:7:1:17:4:0.25:-1:-1	0/1:200:308.53:-35,-5,-35:70:51:18:50:17:33:6:4:17:6:0.25:-1:-1'	#37578214
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:100:137.83:-15,-2,-12:25:17:7:17:7:17:2:5:0:0:0.29:-1:-1	0/0:90:0:-0,-9,-30:30:30:0:30:0:30:0:0:0:0:0:-1:-1	0/1:50:50.85:-8,-3,-15:24:19:4:19:4:19:3:1:0:0:0.17:-1:-1'	#37578579
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:79:146.84:-16,-1,-9:22:14:7:14:7:14:1:5:0:0:0.33:-1:-1	0/1:52:158.87:-17,-1,-6:18:10:7:10:7:10:2:4:0:0:0.41:-1:-1	0/0:45:0:-0,-5,-15:15:15:0:15:0:15:0:0:0:0:0:-1:-1'	#37603021
    assert_in_stdout 'GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:UOPS:UET	0/1:20:20.89:-6,-4,-17:24:20:3:20:3:20:1:1:0:0:0.13:-1:-1	0/0:51:0:-0,-5,-17:18:17:0:17:0:17:0:0:0:0:0:-1:-1	0/0:21:0.03:-2,-4,-15:18:16:1:16:1:16:0:0:0:0:0.059:-1:-1'	#37603051
fi

run phase_snv_vcf_to_bed_ambig \
    unfazed \
        -d $snv_hets_bed \
        -s $sites_vcf \
        -p $ped \
        --include-ambiguous \
        --bam-pairs "NA12878":$bam
if [ $phase_snv_vcf_to_bed_ambig ]; then
    assert_in_stdout '#chrom	start	end	vartype	kid	origin_parent	other_parent	evidence_count	evidence_types'
    assert_in_stdout '22	18844941	18844942	POINT	NA12878	NA12892	NA12891	1	READBACKED'
    assert_in_stdout '22	21088145	21088146	POINT	NA12878	NA12892	NA12891	1	READBACKED'
    assert_in_stdout '22	21141299	21141300	POINT	NA12878	NA12892	NA12891	2	READBACKED'
    assert_in_stdout '22	24656853	24656854	POINT	NA12878	NA12891|NA12892	None	86	AMBIGUOUS_READBACKED'
    assert_in_stdout '22	30857372	30857373	POINT	NA12878	NA12891	NA12892	2	READBACKED'
    assert_in_stdout '22	30857447	30857448	POINT	NA12878	NA12891	NA12892	2	READBACKED'
    assert_in_stdout '22	30862399	30862400	POINT	NA12878	NA12891	NA12892	4	READBACKED'
    assert_in_stdout '22	30864691	30864692	POINT	NA12878	NA12891	NA12892	4	READBACKED'
    assert_in_stdout '22	32756743	32756744	POINT	NA12878	NA12891|NA12892	None	17	AMBIGUOUS_READBACKED'
    assert_in_stdout '22	36556963	36556964	POINT	NA12878	NA12892	NA12891	2	READBACKED'
    assert_in_stdout '22	39387557	39387558	POINT	NA12878	NA12891|NA12892	None	17	AMBIGUOUS_READBACKED'
    assert_in_stdout '22	41566594	41566595	POINT	NA12878	NA12892	NA12891	4	READBACKED'
    assert_in_stdout '22	41609689	41609690	POINT	NA12878	NA12892	NA12891	19	READBACKED'
    assert_in_stdout '22	41613187	41613188	POINT	NA12878	NA12892	NA12891	3	READBACKED'
    assert_in_stdout '22	41613302	41613303	POINT	NA12878	NA12892	NA12891	3	READBACKED'
    assert_in_stdout '22	41652845	41652846	POINT	NA12878	NA12892	NA12891	2	READBACKED'
    assert_in_stdout '22	42072911	42072912	POINT	NA12878	NA12891	NA12892	1	READBACKED'
    assert_in_stdout '22	42523635	42523636	POINT	NA12878	NA12891	NA12892	3	READBACKED'
    assert_in_stdout '22	50617725	50617726	POINT	NA12878	NA12891	NA12892	1	READBACKED'
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
