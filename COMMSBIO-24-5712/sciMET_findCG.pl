#!/usr/bin/perl
# find CG sites; used for Fig. 3j filtering in COMMSBIO-24-5712

use strict;
use warnings;

my $die = "

sciMET_findCG.pl [genome fasta].gz > [CG sites bed file]

";

if (!defined $ARGV[0]) { die $die };

my $IN;
if ($ARGV[0] =~ /gz$/) {
    open $IN, "zcat $ARGV[0] |" or die "Can't open pipe to zcat: $!";
} else {
    open $IN, "$ARGV[0]" or die "Can't open file: $!";
}

my ($pos, $chr, $Nmer);

while (my $l = <$IN>) {
    chomp $l;
    if ($l =~ /^\>/) {
        $pos = 0;
        $chr = $l;
        $chr =~ s/^\>//;
        $chr =~ s/\s.+$//;
        $Nmer = "NN";
    } else {
        my @S = split(//, $l);
        for (my $i = 0; $i < @S; $i++) {
            $pos++;
            $Nmer =~ s/^.//;
            $Nmer .= uc($S[$i]);

            # Check for forward Nmer context
            if ($Nmer eq "CG") {
                my $mBase = $pos - 1;
                print "$chr\t$mBase\t$Nmer\t\+\n";

                my $mBase = $pos;
                print "$chr\t$mBase\t$Nmer\t\-\n";
            }
        }
    }
}

close $IN;
exit;
