#!/bin/csh

setenv ODBCINI /opt/phenix/etc/odbc.ini.mirror

source ~/.cshrc

root -b -q analyze_theo.C\($1\)


