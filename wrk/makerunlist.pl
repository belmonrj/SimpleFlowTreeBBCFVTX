#!/usr/local/bin/perl

print "Advisory: generating the list of runs may take several minutes, please be patient\n";

print "Now making run list \n";
@runs = `ls input/tree_* | awk -F "_" '{print \$2}' | sort -u | awk '{print substr(\$1,5);}'`;

print "Now printing to screen \n";
$nsegments_total = 0;


$execfile = "indexed_runlist.dat";
open(PI,">$execfile") || die "Could not open temp file\n";

$index = 0;
foreach ( @runs )
{
    $runn = $_;
    chomp($runn);
    print "run number $index is $runn\n";

    print PI "$index $runn";
    print PI "\n";

    $index++;

}

close(PI);

print "All done!\n";
