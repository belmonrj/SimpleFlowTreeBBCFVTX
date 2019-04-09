#!/bin/csh
setenv HOME /phenix/u/$LOGNAME
source /etc/csh.login
foreach i (/etc/profile.d/*.csh)
   source $i
end
source $HOME/.login
source /opt/phenix/bin/phenix_setup.csh

echo now running with argument $1

./cumulants.exe $1

