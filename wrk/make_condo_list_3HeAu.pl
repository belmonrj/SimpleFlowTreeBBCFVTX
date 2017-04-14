#!/usr/local/bin/perl



$execfile = "heau_runlist.dat";
open(PI,"<$execfile") || die "Could not open temp file\n";



while(<PI>)
{
    $runn = $_;
    chomp($runn);

    print "Arguments = $runn\n";
    print "Queue 1\n";
    print "\n";

}

close(PI);


print "All done!\n";
