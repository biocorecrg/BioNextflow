#! /usr/bin/env perl

=head1 NAME

=head1 SYNOPSIS

  perl makePeakAnnoReport.pl [-ext "extension to be removed"] [-o output.txt ][-h help] 

=head1 DESCRIPTION

  This script make a table suitable for a bar graph in multiQC

Typical usage is as follows:

  % perl makePeakAnnoReport.pl -ext "_peaks.xls.anno.stats" [-o peak_anno_stats_mqc.txt]

=head2 Options

The following options are accepted:

 --ext=<file extension> 		Specify file exetnsion

 --o=<output file> 				Specify file output file

 --help                   		This documentation.


=head1 AUTHOR

Luca Cozzuto <luca.cozzuto@crg.es> 

=cut
use warnings;
use strict;
use Data::Dumper;
use File::Basename;
use Pod::Usage;
use Getopt::Long;

my $USAGE = "perl makePeakAnnoReport.pl [-ext extension] [-o output] [-h help]";

my ($extension,$output,$show_help);

&GetOptions(    	
			'ext=s'			=> \$extension,
			'output|o=s'    => \$output,
			'help|h'        => \$show_help
			)
  or pod2usage(-verbose=>2);
pod2usage(-verbose=>2) if $show_help;


my @files = glob "'*${extension}*'";


my $header = "# id: peakanno
# plot_type: bargraph
# section_name: Peak annotation statistics
# description: Number of peaks associated with the indicated genomic feature
Sample";

my $fh;
my %results; 
my %features;

foreach my $file (@files) {
	open($fh, "<", $file)
 	or die "Can't open < $file: $!";
 	
	while (my $row = <$fh>) {
		chomp $row;
		if ($row !~ "Features") {
			my @fields = split("\t", $row);
			my $number = $fields[1];
			my $feature = $fields[0];
			$file =~ s/$extension//g;
			$results{$file}{$feature} = $number;
			$features{$feature} = $feature;
		 }
	}
	close($fh);
}

my $body;

my @features = sort {$b cmp $a} keys (%features);
foreach my $feature (@features) {
	$header .= "\t".$feature;
}
$header .= "\n";

foreach my $sampleName (keys (%results)) {
	$header .= $sampleName;
	foreach my $feature (@features) {
		my $featVal = $results{$sampleName}{$feature};
		$header .=  "\t".$featVal;
	}	
	$header .= "\n";

}

# open output file
open(my $fw, ">", $output);
print $fw $header;
close($fw);





