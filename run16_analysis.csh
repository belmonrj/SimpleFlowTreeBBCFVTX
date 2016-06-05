#!/bin/csh

setenv ODBCINI /opt/phenix/etc/odbc.ini.mirror

source ~belmonrj/.cshrc

echo runnumber is $1 and segment number is $2

root -b -q analyze_theo.C\($1,$2\)


