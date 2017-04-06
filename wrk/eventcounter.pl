#!/usr/local/bin/perl



$execfile = "indexed_runlist.dat";
open(PI,"<$execfile") || die "Could not open temp file\n";

$nevents_200 = 0;
$nevents_62 = 0;
$nevents_39 = 0;
$nevents_20 = 0;

while(<PI>)
{
    $line = $_;
    chomp($line);
    #print "the line is $line\n";
    @tokens = split(" ",$line);
    $index = $tokens[0];
    $runn = $tokens[1];
    #print "the index is $index and the run is $runn\n";
    $nevents_run = `grep -m 1 \'Processed\' log/job_$index.out \| awk -F \'\[\ \/\]\' \'\{print \$2\}\'`;
    chomp($nevents_run);
    print "the index is $index and the run is $runn and the number of events is $nevents_run\n";

    if ( not defined $nevents_run ) { $nevents_run = 0; }

    if ( $runn >= 454744 && $runn <= 455639 ) { $nevents_200 = $nevents_200 + $nevents_run; }
    if ( $runn >= 455792 && $runn <= 456283 ) { $nevents_62 =  $nevents_62  + $nevents_run; }
    if ( $runn >= 456652 && $runn <= 457298 ) { $nevents_20 =  $nevents_20  + $nevents_run; }
    if ( $runn >= 457634 && $runn <= 458167 ) { $nevents_39 =  $nevents_39  + $nevents_run; }

}

close(PI);

print "Number of events in 200 GeV is $nevents_200\n";
print "Number of events in 62.4 GeV is $nevents_62\n";
print "Number of events in 39 GeV is $nevents_39\n";
print "Number of events in 19.6 GeV is $nevents_20\n";

print "All done!\n";
