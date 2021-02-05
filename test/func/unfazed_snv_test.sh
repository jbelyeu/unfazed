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
        --quiet \
        -o bed\
        --verbose \
        --bam-pairs "NA12878":$bam
if [ $phase_snv_vcf_to_bed ]; then
    assert_exit_code 0
    assert_in_stdout '#chrom	start	end	vartype	kid	origin_parent	other_parent	evidence_count	evidence_types	origin_parent_sites	origin_parent_reads	other_parent_sites	other_parent_reads'
    assert_in_stdout '22	18844941	18844942	POINT	NA12878	NA12892	NA12891	1	READBACKED	18844298'
    assert_in_stdout '22	21088145	21088146	POINT	NA12878	NA12892	NA12891	1	READBACKED	21087760'
    assert_in_stdout '22	21141299	21141300	POINT	NA12878	NA12892	NA12891	1	READBACKED	21141433'
    assert_in_stdout '22	30857372	30857373	POINT	NA12878	NA12891	NA12892	1	READBACKED	30856856'
    assert_in_stdout '22	30857447	30857448	POINT	NA12878	NA12891	NA12892	1	READBACKED	30856856'
    assert_in_stdout '22	30862399	30862400	POINT	NA12878	NA12891	NA12892	2	READBACKED	30861914,30861930'
    assert_in_stdout '22	30864691	30864692	POINT	NA12878	NA12891	NA12892	2	READBACKED	30864609,30864791'
    assert_in_stdout '22	36556963	36556964	POINT	NA12878	NA12892	NA12891	2	READBACKED	36556697,36556822'
    assert_in_stdout '22	41566594	41566595	POINT	NA12878	NA12892	NA12891	1	READBACKED	41566822'
    assert_in_stdout '22	41609689	41609690	POINT	NA12878	NA12892	NA12891	3	READBACKED	41609429,41609442,41609858'
    assert_in_stdout '22	41613187	41613188	POINT	NA12878	NA12892	NA12891	3	READBACKED	41612540,41612542,41612964'
    assert_in_stdout '22	41613302	41613303	POINT	NA12878	NA12892	NA12891	3	READBACKED	41612540,41612542,41612964'
    assert_in_stdout '22	41652845	41652846	POINT	NA12878	NA12892	NA12891	2	READBACKED	41652344,41652732'
    assert_in_stdout '22	42072911	42072912	POINT	NA12878	NA12891	NA12892	1	READBACKED	42073133'
    assert_in_stdout '22	50617725	50617726	POINT	NA12878	NA12891	NA12892	1	READBACKED	50617982'
fi


run phase_snv_vcf_to_bed_multiread \
    unfazed \
        -d $snv_hets_vcf \
        -s $sites_vcf \
        --quiet \
        -p $ped \
        -o bed\
        --verbose \
        --multiread-proc-min 1 \
        --bam-pairs "NA12878":$bam
if [ $phase_snv_vcf_to_bed_multiread ]; then
    assert_exit_code 0
    assert_in_stdout '#chrom	start	end	vartype	kid	origin_parent	other_parent	evidence_count	evidence_types'
    assert_in_stdout '22	18844941	18844942	POINT	NA12878	NA12892	NA12891	1	READBACKED	18844298'
    assert_in_stdout '22	21088145	21088146	POINT	NA12878	NA12892	NA12891	1	READBACKED	21087760'
    assert_in_stdout '22	21141299	21141300	POINT	NA12878	NA12892	NA12891	1	READBACKED	21141433'
    assert_in_stdout '22	30857372	30857373	POINT	NA12878	NA12891	NA12892	1	READBACKED	30856856'
    assert_in_stdout '22	30857447	30857448	POINT	NA12878	NA12891	NA12892	1	READBACKED	30856856'
    assert_in_stdout '22	30862399	30862400	POINT	NA12878	NA12891	NA12892	2	READBACKED	30861914,30861930'
    assert_in_stdout '22	30864691	30864692	POINT	NA12878	NA12891	NA12892	2	READBACKED	30864609,30864791'
    assert_in_stdout '22	36556963	36556964	POINT	NA12878	NA12892	NA12891	2	READBACKED	36556697,36556822'
    assert_in_stdout '22	41566594	41566595	POINT	NA12878	NA12892	NA12891	1	READBACKED	41566822'
    assert_in_stdout '22	41609689	41609690	POINT	NA12878	NA12892	NA12891	3	READBACKED	41609429,41609442,41609858'
    assert_in_stdout '22	41613187	41613188	POINT	NA12878	NA12892	NA12891	3	READBACKED	41612540,41612542,41612964'
    assert_in_stdout '22	41613302	41613303	POINT	NA12878	NA12892	NA12891	3	READBACKED	41612540,41612542,41612964'
    assert_in_stdout '22	41652845	41652846	POINT	NA12878	NA12892	NA12891	2	READBACKED	41652344,41652732'
    assert_in_stdout '22	42072911	42072912	POINT	NA12878	NA12891	NA12892	1	READBACKED	42073133'
    assert_in_stdout '22	50617725	50617726	POINT	NA12878	NA12891	NA12892	1	READBACKED	50617982'
fi

run phase_snv_vcf_to_bed_multiread1 \
    unfazed \
        -d $snv_hets_vcf \
        -s $sites_vcf \
        --quiet \
        -p $ped \
        -o bed\
        -t 1\
        --verbose \
        --multiread-proc-min 1 \
        --bam-pairs "NA12878":$bam
if [ $phase_snv_vcf_to_bed_multiread1 ]; then
    assert_exit_code 0
    assert_in_stdout '#chrom	start	end	vartype	kid	origin_parent	other_parent	evidence_count	evidence_types'
    assert_in_stdout '22	18844941	18844942	POINT	NA12878	NA12892	NA12891	1	READBACKED	18844298'
    assert_in_stdout '22	21088145	21088146	POINT	NA12878	NA12892	NA12891	1	READBACKED	21087760'
    assert_in_stdout '22	21141299	21141300	POINT	NA12878	NA12892	NA12891	1	READBACKED	21141433'
    assert_in_stdout '22	30857372	30857373	POINT	NA12878	NA12891	NA12892	1	READBACKED	30856856'
    assert_in_stdout '22	30857447	30857448	POINT	NA12878	NA12891	NA12892	1	READBACKED	30856856'
    assert_in_stdout '22	30862399	30862400	POINT	NA12878	NA12891	NA12892	2	READBACKED	30861914,30861930'
    assert_in_stdout '22	30864691	30864692	POINT	NA12878	NA12891	NA12892	2	READBACKED	30864609,30864791'
    assert_in_stdout '22	36556963	36556964	POINT	NA12878	NA12892	NA12891	2	READBACKED	36556697,36556822'
    assert_in_stdout '22	41566594	41566595	POINT	NA12878	NA12892	NA12891	1	READBACKED	41566822'
    assert_in_stdout '22	41609689	41609690	POINT	NA12878	NA12892	NA12891	3	READBACKED	41609429,41609442,41609858'
    assert_in_stdout '22	41613187	41613188	POINT	NA12878	NA12892	NA12891	3	READBACKED	41612540,41612542,41612964'
    assert_in_stdout '22	41613302	41613303	POINT	NA12878	NA12892	NA12891	3	READBACKED	41612540,41612542,41612964'
    assert_in_stdout '22	41652845	41652846	POINT	NA12878	NA12892	NA12891	2	READBACKED	41652344,41652732'
    assert_in_stdout '22	42072911	42072912	POINT	NA12878	NA12891	NA12892	1	READBACKED	42073133'
    assert_in_stdout '22	50617725	50617726	POINT	NA12878	NA12891	NA12892	1	READBACKED	50617982'
fi


run phase_snv_bed_to_bed \
    unfazed \
        --verbose \
        -d $snv_hets_bed \
        -s $sites_vcf \
        --quiet \
        -p $ped \
        -o bed\
        --bam-pairs "NA12878":$bam
if [ $phase_snv_bed_to_bed ]; then
    assert_exit_code 0
    assert_in_stdout '#chrom	start	end	vartype	kid	origin_parent	other_parent	evidence_count	evidence_types'
    assert_in_stdout '22	18844941	18844942	POINT	NA12878	NA12892	NA12891	1	READBACKED	18844298'
    assert_in_stdout '22	21088145	21088146	POINT	NA12878	NA12892	NA12891	1	READBACKED	21087760'
    assert_in_stdout '22	21141299	21141300	POINT	NA12878	NA12892	NA12891	1	READBACKED	21141433'
    assert_in_stdout '22	30857372	30857373	POINT	NA12878	NA12891	NA12892	1	READBACKED	30856856'
    assert_in_stdout '22	30857447	30857448	POINT	NA12878	NA12891	NA12892	1	READBACKED	30856856'
    assert_in_stdout '22	30862399	30862400	POINT	NA12878	NA12891	NA12892	2	READBACKED	30861914,30861930'
    assert_in_stdout '22	30864691	30864692	POINT	NA12878	NA12891	NA12892	2	READBACKED	30864609,30864791'
    assert_in_stdout '22	36556963	36556964	POINT	NA12878	NA12892	NA12891	2	READBACKED	36556697,36556822'
    assert_in_stdout '22	41566594	41566595	POINT	NA12878	NA12892	NA12891	1	READBACKED	41566822'
    assert_in_stdout '22	41609689	41609690	POINT	NA12878	NA12892	NA12891	3	READBACKED	41609429,41609442,41609858'
    assert_in_stdout '22	41613187	41613188	POINT	NA12878	NA12892	NA12891	3	READBACKED	41612540,41612542,41612964'
    assert_in_stdout '22	41613302	41613303	POINT	NA12878	NA12892	NA12891	3	READBACKED	41612540,41612542,41612964'
    assert_in_stdout '22	41652845	41652846	POINT	NA12878	NA12892	NA12891	2	READBACKED	41652344,41652732'
    assert_in_stdout '22	42072911	42072912	POINT	NA12878	NA12891	NA12892	1	READBACKED	42073133'
    assert_in_stdout '22	50617725	50617726	POINT	NA12878	NA12891	NA12892	1	READBACKED	50617982'
fi

#split the next two tests to decrease the amount of test to stdout
rm -f out.tmp
unfazed \
    -d $snv_hets_vcf \
    -s $sites_vcf \
    --quiet \
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
    assert_in_stdout '##FORMAT=<ID=UET,Number=1,Type=Float,Description="Unfazed evidence type: `0` (readbacked), `1` (allele-balance, for CNVs only), `2` (both), `3` (ambiguous readbacked), `4` (ambiguous allele-balance), `5` (ambiguous both), `6` (auto-phased sex-chromosome variant in male), or `-1` (missing)">'
fi

run phase_snv_vcf_to_vcf_body \
    grep -v "#" out.tmp | grep "UOPS" | cut -f 9-12
if [ $phase_snv_vcf_to_vcf_body ]; then
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:6,13:19:99:.:.:409,0,140:-1:-1	0/1:8,11:19:99:.:.:354,0,232:-1:-1	0/1:7,22:29:99:.:.:747,0,163:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PL:UOPS:UET	0|1:55,26:81:99:530,0,1495:1:0	0/1:77,11:88:52:52,0,2094:-1:-1	0/1:35,46:81:99:1202,0,836:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:13,17:30:99:.:.:517,0,329:-1:-1	0/1:27,11:38:99:.:.:273,0,803:-1:-1	0/1:38,36:74:99:.:.:969,0,1043:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PL:UOPS:UET	0/1:16,8:24:99:208,0,474:-1:-1	0/1:12,17:29:99:543,0,366:-1:-1	0/1:21,13:34:99:393,0,718:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:14,17:31:99:.:.:436,0,315:-1:-1	0/1:14,10:24:99:.:.:267,0,406:-1:-1	0/1:13,13:26:99:.:.:374,0,423:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:90,65:155:99:.:.:1655,0,2354:-1:-1	0/1:77,84:161:99:.:.:2396,0,2053:-1:-1	0/1:36,41:77:99:.:.:1173,0,1097:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:18,30:48:99:.:.:918,0,422:-1:-1	0/1:24,19:43:99:.:.:523,0,620:-1:-1	0/1:31,37:68:99:.:.:1038,0,765:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:28,22:50:99:.:.:641,0,854:-1:-1	0/1:18,23:41:99:.:.:607,0,511:-1:-1	0/1:63,62:125:99:.:.:1703,0,2087:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PL:UOPS:UET	0|1:34,38:72:99:991,0,812:1:0	0/1:29,26:55:99:758,0,707:-1:-1	0/1:46,57:103:99:1655,0,1073:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0|1:45,58:103:99:.:.:1562,0,1256:1:0	0/1:49,62:111:99:.:.:1754,0,1331:-1:-1	0/1:69,49:118:99:.:.:1459,0,2013:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PL:UOPS:UET	0/1:31,21:52:99:506,0,965:-1:-1	0/1:48,36:84:99:1009,0,1455:-1:-1	0/1:56,64:120:99:1750,0,1828:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:15,72:87:99:.:.:3091,0,707:-1:-1	0/1:85,82:167:99:.:.:3190,0,5131:-1:-1	0/1:43,65:108:99:.:.:2660,0,2846:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:10,96:106:99:.:.:4137,0,346:-1:-1	0/1:112,72:184:99:.:.:2805,0,7213:-1:-1	0/1:85,102:187:99:.:.:3955,0,5468:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:12,74:86:99:.:.:2152,0,189:-1:-1	0/1:68,82:150:99:.:.:2239,0,1981:-1:-1	0/1:100,135:235:99:.:.:3607,0,2895:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:15,55:70:99:.:.:2404,0,713:-1:-1	0/1:84,93:177:99:.:.:3777,0,5325:-1:-1	0/1:79,89:168:99:.:.:3699,0,4902:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:18,45:63:99:.:.:2184,0,1822:-1:-1	0/1:87,61:148:99:.:.:2497,0,8299:-1:-1	0/1:77,84:161:99:.:.:3482,0,7651:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:18,63:81:99:.:.:1951,0,331:-1:-1	0/1:94,69:163:99:.:.:1773,0,2623:-1:-1	0/1:142,120:262:99:.:.:3317,0,3649:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PL:UOPS:UET	0/1:15,82:97:99:2576,0,290:-1:-1	0/1:74,79:153:99:2296,0,2237:-1:-1	0/1:157,163:320:99:4804,0,5206:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:74,90:164:99:.:.:2680,0,1842:-1:-1	0/1:88,91:179:99:.:.:2612,0,2316:-1:-1	0/1:29,29:58:99:.:.:809,0,773:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:44,48:92:99:.:.:1420,0,1149:-1:-1	0/1:40,45:85:99:.:.:1410,0,1093:-1:-1	0/1:40,29:69:99:.:.:831,0,1067:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PL:UOPS:UET	0/1:48,30:78:99:853,0,1418:-1:-1	0/1:49,36:85:99:1037,0,1327:-1:-1	0/1:43,42:85:99:1286,0,1258:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PL:UOPS:UET	0/1:43,20:63:99:423,0,1323:-1:-1	0/1:47,33:80:99:778,0,1460:-1:-1	0/1:41,39:80:99:981,0,1274:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PL:UOPS:UET	0/1:8,9:17:99:264,0,220:-1:-1	0/1:13,7:20:99:179,0,367:-1:-1	0/1:7,13:20:99:406,0,171:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PL:UOPS:UET	0/1:7,12:19:99:344,0,186:-1:-1	0/1:7,14:21:99:401,0,141:-1:-1	0/1:32,29:61:99:758,0,962:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:30,27:57:99:.:.:806,0,828:-1:-1	0/1:19,23:42:99:.:.:682,0,537:-1:-1	0/1:56,64:120:99:.:.:2166,0,1590:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PL:UOPS:UET	0/1:38,36:74:99:1097,0,1140:-1:-1	0/1:25,22:47:99:665,0,734:-1:-1	0/1:69,64:133:99:2028,0,2042:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:43,46:89:99:.:.:1158,0,1281:-1:-1	0/1:61,44:105:99:.:.:1097,0,1871:-1:-1	0/1:99,115:214:99:.:.:2794,0,2806:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:32,27:59:99:.:.:789,0,819:-1:-1	0/1:26,34:60:99:.:.:1026,0,631:-1:-1	0/1:94,69:163:99:.:.:2112,0,2525:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PL:UOPS:UET	0/1:58,51:109:99:1471,0,1706:-1:-1	0/1:40,37:77:99:1079,0,1153:-1:-1	0/1:191,145:336:99:4091,0,6567:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:20,22:42:99:.:.:591,0,447:-1:-1	0/1:21,15:36:99:.:.:425,0,560:-1:-1	0/1:93,104:197:99:.:.:3111,0,2409:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:49,34:83:99:.:.:1010,0,1552:-1:-1	0/1:37,23:60:99:.:.:596,0,1149:-1:-1	0/1:69,73:142:99:.:.:2053,0,2209:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:140,101:241:99:.:.:2626,0,3429:-1:-1	0/1:106,88:194:99:.:.:2352,0,2883:-1:-1	0/1:26,36:62:99:.:.:1106,0,682:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:16,24:40:99:.:.:752,0,421:-1:-1	0/1:27,23:50:99:.:.:699,0,811:-1:-1	0/1:15,18:33:99:.:.:525,0,348:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PL:UOPS:UET	0/1:23,14:37:99:346,0,674:-1:-1	0/1:31,22:53:99:573,0,889:-1:-1	0/1:60,52:112:99:1493,0,1618:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:42,28:70:99:.:.:785,0,1167:-1:-1	0/1:39,36:75:99:.:.:1074,0,1085:-1:-1	0/1:45,46:91:99:.:.:1511,0,1267:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:44,55:99:99:.:.:1621,0,1312:-1:-1	0/1:63,53:116:99:.:.:1440,0,1966:-1:-1	0/1:7,8:15:99:.:.:224,0,212:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PL:UOPS:UET	1|0:47,47:94:99:1478,0,1355:1:0	0/1:48,38:86:99:1166,0,1416:-1:-1	0/1:72,73:145:99:2241,0,1931:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PL:UOPS:UET	1|0:32,31:63:99:861,0,932:1:0	0/1:33,40:73:99:1058,0,876:-1:-1	0/1:60,66:126:99:1904,0,1635:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	1|0:39,22:61:99:.:.:513,0,1059:2:0	0/1:37,32:69:99:.:.:853,0,1105:-1:-1	0/1:99,110:209:99:.:.:3046,0,2841:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PL:UOPS:UET	1|0:41,30:71:99:785,0,1319:2:0	0/1:59,47:106:99:1253,0,1769:-1:-1	0/1:60,57:117:99:1644,0,1948:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PL:UOPS:UET	0/1:61,71:132:99:2062,0,1808:-1:-1	0/1:76,57:133:99:1574,0,2368:-1:-1	0/1:101,102:203:99:2862,0,3398:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:23,11:34:99:.:.:301,0,671:-1:-1	0/1:24,16:40:99:.:.:623,0,1369:-1:-1	0/1:60,45:105:99:.:.:1272,0,1773:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:22,9:31:99:.:.:156,0,653:-1:-1	0/1:21,17:38:99:.:.:623,0,1369:-1:-1	0/1:56,44:100:99:.:.:1190,0,1561:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:73,46:119:99:.:.:1079,0,1898:-1:-1	0/1:71,61:132:99:.:.:1652,0,1861:-1:-1	0/1:27,28:55:99:.:.:866,0,663:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PL:UOPS:UET	0/1:34,37:71:99:870,0,997:-1:-1	0/1:29,23:52:99:583,0,869:-1:-1	0/1:21,18:39:99:533,0,633:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0|1:7,16:23:99:.:.:323,0,169:2:0	0/1:4,11:15:99:.:.:315,0,99:-1:-1	0/1:65,56:121:99:.:.:1335,0,1984:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:33,14:47:99:.:.:351,0,842:-1:-1	0/1:34,38:72:99:.:.:1131,0,925:-1:-1	0/1:112,116:228:99:.:.:3424,0,2962:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PL:UOPS:UET	0/1:16,12:28:99:278,0,445:-1:-1	0/1:9,10:19:99:266,0,185:-1:-1	0/1:21,19:40:99:552,0,439:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:27,42:69:99:.:.:1230,0,668:-1:-1	0/1:27,35:62:99:.:.:1055,0,674:-1:-1	0/1:18,13:31:99:.:.:337,0,462:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PL:UOPS:UET	0/1:30,39:69:99:1140,0,718:-1:-1	0/1:37,46:83:99:1324,0,910:-1:-1	0/1:25,31:56:99:878,0,712:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PL:UOPS:UET	0/1:23,30:53:99:810,0,586:-1:-1	0/1:37,35:72:99:851,0,924:-1:-1	0/1:17,17:34:99:415,0,386:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:28,34:62:99:.:.:1027,0,759:-1:-1	0/1:35,42:77:99:.:.:1239,0,993:-1:-1	0/1:13,21:34:99:.:.:666,0,328:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:30,22:52:99:.:.:603,0,762:-1:-1	0/1:37,20:57:99:.:.:516,0,1040:-1:-1	0/1:28,26:54:99:.:.:742,0,698:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:26,24:50:99:.:.:643,0,631:-1:-1	0/1:20,15:35:99:.:.:370,0,551:-1:-1	0/1:21,19:40:99:.:.:561,0,650:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PL:UOPS:UET	0/1:24,29:53:99:824,0,658:-1:-1	0/1:34,30:64:99:918,0,995:-1:-1	0/1:20,18:38:99:529,0,667:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:37,44:81:99:.:.:1178,0,974:-1:-1	0/1:49,44:93:99:.:.:1165,0,1298:-1:-1	0/1:28,16:44:99:.:.:457,0,797:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:112,116:228:99:.:.:3094,0,2721:14:3	0/1:161,103:264:99:.:.:2647,0,4177:-1:-1	0/1:150,137:287:99:.:.:3611,0,3774:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PL:UOPS:UET	0/1:92,107:199:99:3106,0,2566:-1:-1	0/1:100,145:245:99:4149,0,2619:-1:-1	0/1:54,48:102:99:1519,0,1564:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PL:UOPS:UET	0|1:24,26:50:99:804,0,631:1:0	0/1:47,35:82:99:972,0,1276:-1:-1	0/1:87,82:169:99:2328,0,2791:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:20,22:42:99:.:.:656,0,494:-1:-1	0/1:16,31:47:99:.:.:940,0,339:-1:-1	0/1:46,40:86:99:.:.:1160,0,1142:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0|1:12,9:21:99:.:.:231,0,378:3:0	0/1:9,9:18:99:.:.:264,0,252:-1:-1	0/1:23,37:60:99:.:.:1070,0,666:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0|1:39,39:78:99:.:.:1145,0,955:3:0	0/1:39,38:77:99:.:.:1077,0,963:-1:-1	0/1:57,73:130:99:.:.:2175,0,1372:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0|1:7,10:17:99:.:.:310,0,224:3:0	0/1:10,8:18:99:.:.:219,0,325:-1:-1	0/1:20,17:37:99:.:.:435,0,678:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PL:UOPS:UET	0/1:27,29:56:99:850,0,705:-1:-1	0/1:25,36:61:99:1025,0,607:-1:-1	0/1:107,109:216:99:3038,0,2904:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0|1:12,11:23:99:.:.:322,0,359:2:0	0/1:13,12:25:99:.:.:361,0,372:-1:-1	0/1:14,22:36:99:.:.:604,0,417:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PL:UOPS:UET	0/1:29,30:59:99:880,0,796:-1:-1	0/1:37,39:76:99:1157,0,1071:-1:-1	0/1:47,47:94:99:1404,0,1419:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:56,58:114:99:.:.:1714,0,1486:-1:-1	0/1:45,59:104:99:.:.:1735,0,1120:-1:-1	0/1:24,41:65:99:.:.:1163,0,546:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:17,14:31:99:.:.:432,0,427:-1:-1	0/1:26,15:41:99:.:.:421,0,673:-1:-1	0/1:64,68:132:99:.:.:1839,0,1601:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:32,32:64:99:.:.:883,0,795:-1:-1	0/1:34,35:69:99:.:.:977,0,942:-1:-1	0/1:23,40:63:99:.:.:1161,0,564:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PL:UOPS:UET	0/1:8,13:21:99:395,0,218:-1:-1	0/1:18,12:30:99:270,0,456:-1:-1	0/1:8,8:16:99:250,0,237:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PL:UOPS:UET	0/1:59,46:105:99:1311,0,1923:-1:-1	0/1:68,62:130:99:1754,0,2140:-1:-1	0/1:74,52:126:99:1462,0,2636:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	1|0:65,70:135:99:.:.:2030,0,1951:1:0	0/1:69,60:129:99:.:.:1697,0,2070:-1:-1	0/1:190,183:373:99:.:.:5271,0,6469:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:32,47:79:99:.:.:1344,0,907:-1:-1	0/1:48,30:78:99:.:.:805,0,1412:-1:-1	0/1:91,79:170:99:.:.:2335,0,3055:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:81,86:167:99:.:.:2487,0,2260:-1:-1	0/1:70,64:134:99:.:.:1810,0,1969:-1:-1	0/1:26,24:50:99:.:.:730,0,823:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:22,27:49:99:.:.:826,0,561:-1:-1	0/1:20,33:53:99:.:.:957,0,516:-1:-1	0/1:52,26:78:99:.:.:667,0,1437:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:126,136:262:99:.:.:4092,0,3799:-1:-1	0/1:121,115:236:99:.:.:3238,0,3763:-1:-1	0/1:151,147:298:99:.:.:4085,0,4959:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PL:UOPS:UET	1|0:4,8:12:88:226,0,88:1:0	0/1:7,8:15:99:241,0,190:-1:-1	0/1:10,12:22:99:332,0,217:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PL:UOPS:UET	0/1:19,17:36:99:502,0,514:-1:-1	0/1:10,13:23:99:401,0,225:-1:-1	0/1:17,22:39:99:601,0,426:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PL:UOPS:UET	0/1:15,12:27:99:346,0,401:-1:-1	0/1:14,12:26:99:317,0,394:-1:-1	0/1:9,11:20:99:306,0,248:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:34,49:83:99:.:.:1336,0,937:-1:-1	0/1:41,38:79:99:.:.:1007,0,1222:-1:-1	0/1:22,26:48:99:.:.:623,0,650:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:93,75:168:99:.:.:1981,0,2722:-1:-1	0/1:80,65:145:99:.:.:1698,0,2271:-1:-1	0/1:72,70:142:99:.:.:1890,0,2104:-1:-1'
    assert_in_stdout 'GT:AD:DP:GQ:PGT:PID:PL:UOPS:UET	0/1:18,21:39:99:.:.:537,0,436:-1:-1	0/1:24,21:45:99:.:.:524,0,674:-1:-1	0/1:10,8:18:99:.:.:190,0,237:-1:-1'
fi

run phase_snv_vcf_to_bed_ambig \
    unfazed \
        -d $snv_hets_bed \
        -s $sites_vcf \
        --quiet \
        -p $ped \
        --include-ambiguous \
        --bam-pairs "NA12878":$bam
if [ $phase_snv_vcf_to_bed_ambig ]; then
    assert_exit_code 0
    assert_in_stdout '#chrom	start	end	vartype	kid	origin_parent	other_parent	evidence_count	evidence_types'
    assert_in_stdout '22	18844941	18844942	POINT	NA12878	NA12892	NA12891	1	READBACKED'
    assert_in_stdout '22	21088145	21088146	POINT	NA12878	NA12892	NA12891	1	READBACKED'
    assert_in_stdout '22	21141299	21141300	POINT	NA12878	NA12892	NA12891	1	READBACKED'
    assert_in_stdout '22	30857372	30857373	POINT	NA12878	NA12891	NA12892	1	READBACKED'
    assert_in_stdout '22	30857447	30857448	POINT	NA12878	NA12891	NA12892	1	READBACKED'
    assert_in_stdout '22	30862399	30862400	POINT	NA12878	NA12891	NA12892	2	READBACKED'
    assert_in_stdout '22	30864691	30864692	POINT	NA12878	NA12891	NA12892	2	READBACKED'
    assert_in_stdout '22	36556963	36556964	POINT	NA12878	NA12892	NA12891	2	READBACKED'
    assert_in_stdout '22	39387557	39387558	POINT	NA12878	NA12891|NA12892	None	10	AMBIGUOUS_READBACKED'
    assert_in_stdout '22	41566594	41566595	POINT	NA12878	NA12892	NA12891	1	READBACKED'
    assert_in_stdout '22	41609689	41609690	POINT	NA12878	NA12892	NA12891	3	READBACKED'
    assert_in_stdout '22	41613187	41613188	POINT	NA12878	NA12892	NA12891	3	READBACKED'
    assert_in_stdout '22	41613302	41613303	POINT	NA12878	NA12892	NA12891	3	READBACKED'
    assert_in_stdout '22	41652845	41652846	POINT	NA12878	NA12892	NA12891	2	READBACKED'
    assert_in_stdout '22	42072911	42072912	POINT	NA12878	NA12891	NA12892	1	READBACKED'
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
