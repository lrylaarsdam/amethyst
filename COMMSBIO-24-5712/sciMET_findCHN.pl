#!/usr/bin/perl
# find CHN sites

use strict;
use warnings;

my $die = "

sciMET_findCHN.pl [genome fasta].gz > [CHN sites bed file]

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
        $Nmer = "NNN";
    } else {
        my @S = split(//, $l);
        for (my $i = 0; $i < @S; $i++) {
            $pos++;
            $Nmer =~ s/^.//;
            $Nmer .= uc($S[$i]);

            # Check for forward Nmer context
            if ($Nmer =~ /^C[^G]/) {
                my $mBase = $pos - 2;
                print "$chr\t$mBase\t$Nmer\t\+\n";
            }

            # Check for reverse Nmer context and print reverse complement
            if ($Nmer =~ /.[^C]G$/) {
                my $revcomp = reverse $Nmer;
                $revcomp =~ tr/ATGC/TACG/;
                my $mBase = $pos;
                print "$chr\t$mBase\t$revcomp\t\-\n";
            }
        }
    }
}

close $IN;
exit;
