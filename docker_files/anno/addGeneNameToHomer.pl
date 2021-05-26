#! /usr/bin/env perl

=head1 NAME

=head1 SYNOPSIS

  perl addGeneNameToHomer.pl [-g gtf file] [-o output][-f feature] [-a homer gtf_annotation] [-h help]

=head1 DESCRIPTION

  This script reads a gtf_annotation file produced by home and attach the corresponding gene name
  
Typical usage is as follows:

  % perl addGeneNameToHomer.pl -g genome.gtf -a anno.txt 
  
=head2 Options

The following options are accepted:

 --g=<genome gtf file>   				(Mandatory).

 --a=<gtf_annotation file>  		 		(Mandatory).

 --o=<output gtf_annotation file>  		 	(Default: output.anno).

 --f=<feature>  		 				(Default: transcript).

 --fid=<featureid>					(Default: the feature specified with -f + "_id").

 --help                   				This documentation.


=head1 AUTHOR

Luca Cozzuto <luca.cozzuto@crg.es> 

=cut
use warnings;
use strict;
use Data::Dumper;
use File::Basename;
use Pod::Usage;
use Getopt::Long;
use Bio::Tools::GFF;


my $USAGE = "perl extractSeqs.pl [-g genome file] [-f feature] [-a gtf_annotation] [-h help]";

my ($input,$gtf_annotation,$output, $chosen_feature, $feature_id, $show_help);

&GetOptions(    	
			'annotation|a=s'		=> \$input,
			'output|o=s'			=> \$output,
			'gtf_annotation|g=s'	=> \$gtf_annotation,
			'feature|f=s'		=> \$chosen_feature,
			'featureid|fid=s'	=> \$feature_id,
			'help|h'        	=> \$show_help
			)
  or pod2usage(-verbose=>2);
pod2usage(-verbose=>2) if $show_help;

if (!$input) { die ("Please specify input file")}
if (!$gtf_annotation) { die ("Please specify a gtf file with annotation")}
if (!$output) { $output="output.anno"}
if (!$chosen_feature) { $chosen_feature="transcript"}
if (!$feature_id) { $feature_id=$chosen_feature."_id"}


# open peak gtf_annotation input file (made by homer)
open(my $fa, "<", $input)
        or die "Can't open < $input: $!";

# open genomics gtf_annotation file (GTF)
my $fl;
if ($gtf_annotation =~ /\.gz$/i) {
        open($fl, "gunzip -c $gtf_annotation |") or die "gunzip $gtf_annotation: $!";
} else {
	open($fl, "<", $gtf_annotation)
                or die "Can't open < $gtf_annotation: $!";
}

 # open output file
open(my $fw, ">", $output);

my %names; 

my $feature;
while (my $gtf_row = <$fl>) {
    chomp $gtf_row;
    if ($gtf_row !~ "#") {
	    my @fields = split("\t", $gtf_row);
	    my $desc = $fields[8];
	    $feature = $fields[2];
        if ($feature eq "$chosen_feature") {
	    	my $params = getParams($desc);
			my $id = $params->{$feature_id};
			$names{$id} = $params->{'gene_name'};
	    }
   	 }
}
close($fl);


my $genename;
 while (my $row = <$fa>) {
    chomp $row;
	$row =~ s/\t+$//;
    my  @fields = split("\t", $row);
    if ($fields[2] eq "Start" ) {
		splice @fields, 11;
		$row = join ("\t", @fields);
    	$row .="\tGene Name";
    }
    else {
    	my $geneid = $fields[10];
    	$genename = "NA";
    	if ($names{$geneid}) {
    		$genename = $names{$geneid};
    	} 
    	$row .="\t$genename";
    }
 	print ($fw $row."\n");
}
close($fw);

sub  getParams {
        my ($string) = @_;
        my %fulldesc;
        $string =~ s/"//g;
        my @descs = split(";",$string);
        foreach my $desc (@descs) {
                my @vals = split(" ", $desc);
                $fulldesc{$vals[0]} = $vals[1];
        }
        return \%fulldesc;
}

