#!/usr/local/bin/perl


#for ( $i = 360; $i < 363; $i++ )
$execfile = "runbyrun_fromlogfiles.dat";
open(PI,">$execfile") || die "Could not open temp file\n";

for ( $i = 0; $i < 363; $i++ )
{

    $logfile = "log/job_".$i.".out";
    $runn = `grep runNumber $logfile | awk \'\{print \$3\}\'`;
    chomp($runn);
    $cratio = `grep Ratio $logfile | awk \'\{print \$NF\}\' | sed -e \'s/(//g\' -e \'s/)//g\'`;
    chomp($cratio);
    $vratio = `grep outer $logfile | awk \'\{print \$NF\}\' | sed -e \'s/(//g\' -e \'s/)//g\'`;
    chomp($vratio);
    $eratio = `grep \'number of\' $logfile | awk \'\{print \$NF\}\' | sed -e \'s/(//g\' -e \'s/)//g\'`;
    chomp($eratio);
    print "log file is $logfile and run number is $runn and cluster ratio is $cratio and vertex ratio is $vratio and event ratio passing cluster cut is $eratio\n";
    print PI "$runn\t$cratio\t$vratio\t$eratio\n";

}

close(PI);
