#!/bin/sh

ls production_Run16dAu200GeV/run_*/CNT/CNT_* | awk -F "_" '{print $7}' | awk -F "-" '{print $2}' | sort -u

