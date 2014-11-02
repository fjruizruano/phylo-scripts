#!/usr/bin/perl -w

=head1 NAME

raxml_nexusPartConvert.pl

=head1 SYNOPSIS

perl raxml_nexusPartConvert.pl -m aln.nex -r DNA > raxml_partition_file

Options:

 --help        Show brief help and exit
 --matrixfile  Your input alignment in Nexus format
 --raxstring   DNA or AA model to be used

=head1 DESCRIPTION

Given a Nexus matrix, generate a RAxML partition guide file

=head1 AUTHOR

Apurva Narechania
anarechania *a|t* amnh.org

=head1 COPYRIGHT

Copyright (c) 2008-2009 American Museum of Natural History

This library is free software;  
you can redistribute it and/or modify 
it under the same terms as Perl itself.

=cut

# ----------------------------------------------------

#####TODO:
#####     

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;

my ($help, $matrixfile, $raxstring);
GetOptions(
    'h|help'          => \$help,
    'm|matrixfile=s'  => \$matrixfile,
    'r|raxstring=s'   => \$raxstring,
	   ) or pod2usage;

pod2usage if $help;

#####MAIN#####

# parse some basic stuff out of the Nexus file
my $nexus = parse_nexus ($matrixfile);

# write out the raxml format for partitions
	foreach my $cs (sort keys %{$nexus->{'charset'}}){
    	print "$raxstring, $cs=$nexus->{'charset'}->{$cs}\n";
	}

#####SUBS#####

sub parse_nexus{
    my $alignfile = shift;

    open (NEX, "$alignfile");
    my $charset = {};
    my $nexus   = {};
    while (my $line = <NEX>){
        chomp $line;

        # take only first instances of all of these things                                           
        # header only                                                                              
        ($nexus->{'nchar'}    = $1) if (($line =~m/nchar\s*=\s*(\d+)/i) and (!$nexus->{'nchar'}));
        ($nexus->{'ntax'}     = $1) if (($line =~m/ntax\s*=\s*(\d+)/i) and (!$nexus->{'ntax'}));
        ($nexus->{'datatype'} = $1) if (($line =~m/datatype\s*=\s*(\w+)/i) and (!$nexus->{'datatype'}));
        ($nexus->{'missing'}  = $1) if (($line =~m/missing\s*=\s*(.{1})/i) and (!$nexus->{'missing'}));
        ($nexus->{'gap'}      = $1) if (($line =~m/gap\s*=\s*(.{1})/i) and (!$nexus->{'gap'}));

        if ($line =~m/outgroup/i){
            $line =~s/outgroup//ig;
            $line =~s/\s+//g;
            $line =~s/\;//g;

            # any instances of more than one outgroup???? <====FIX                              
            $nexus->{'outgroup'} = $line;
        }

        if ($line =~m/charset/i){
            $line =~s/charset//ig;
            $line =~s/\s+//g;
            $line =~s/\;//g;

            my ($partition, $coords) =
		split (/\=/, $line);
            $charset->{$partition} = $coords;
        }
        $nexus->{'charset'} = $charset;
    }
    close (NEX);

    return ($nexus);
}
