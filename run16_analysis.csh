#!/bin/csh

setenv ODBCINI /opt/phenix/etc/odbc.ini.mirror

#set x = $PWD

#cd ~;
#echo $PWD
#source .cshrc
#cd $x;
source ~/.cshrc

root -b -q analyze_theo.C\($1\) >& $2

exit;

