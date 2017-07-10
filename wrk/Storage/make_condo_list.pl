#!/usr/local/bin/perl



$execfile = "all.list";
open(PI,"<$execfile") || die "Could not open temp file\n";

$nruns_200 = 0;
$nruns_62  = 0;
$nruns_39  = 0;
$nruns_20  = 0;

while(<PI>)
{
    $runn = $_;
    chomp($runn);

    if ( $runn >= 454744 && $runn <= 455639 ) { $nruns_200 = $nruns_200 + 1; }
    if ( $runn >= 455792 && $runn <= 456283 ) { $nruns_62  = $nruns_62  + 1; }
    if ( $runn >= 456652 && $runn <= 457298 ) { $nruns_20  = $nruns_20  + 1; }
    if ( $runn >= 457634 && $runn <= 458167 ) { $nruns_39  = $nruns_39  + 1; }
    print "Arguments = $runn\n";
    print "Queue 1\n";
    print "\n";

}

close(PI);

print "Number of runs in 200 GeV is $nruns_200\n";
print "Number of runs in 62.4 GeV is $nruns_62\n";
print "Number of runs in 39 GeV is $nruns_39\n";
print "Number of runs in 19.6 GeV is $nruns_20\n";

print "All done!\n";
