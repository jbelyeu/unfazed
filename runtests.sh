#!/bin/bash
set -e
echo "running functional tests:"
bash test/func/unfazed_snv_test.sh
bash test/func/unfazed_sv_test.sh
echo "finished functional tests"
