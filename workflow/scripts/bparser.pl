#!/usr/local/bin/perl
use strict;
use warnings;
use Bio::SearchIO;

# Usage information
die "Usage: $0 <BLAST-report-file> <number-of-top-hits> <output-file>\n", if (@ARGV != 3);

my ($infile,$numHits,$outfile) = @ARGV;
print "Parsing the BLAST result ...";
my $in = Bio::SearchIO->new(-format => 'blast', -file => $infile);
open (OUT,">$outfile") or die "Cannot open $outfile: $!";

# print the header info for tab-deliminated columns
print OUT "query_name\tquery_length\taccession_number\tlength\tdescription\tE value\tbit score\tframe\tquery_start\t";
print OUT "query_end\thit_start\thit_end\tidentical\t% Query Aligned\n";

# extraction of information for each result recursively
while ( my $result = $in->next_result ) {
	# the name of the query sequence
   	print OUT $result->query_name . "\t";

        # the length of the query sequence
    	print OUT $result->query_length;

        # output "no hits found" if there is no hits
    	if ( $result->num_hits == 0 ) {
		print OUT "\tNo hits found\n";
    	} else {
		my $count = 0;

                # process each hit recursively
		while (my $hit = $result->next_hit) {
			if($count != 0){
				 # the name of the query sequence
	        		print OUT $result->query_name. "\t";

        			# the length of the query sequence
        			print OUT $result->query_length;
			}
			#print OUT "\t" if ($count = 0);
                        # get the accession numbers of the hits
			print OUT "\t" . $hit->accession . "\t";
                        # get the lengths of the hit sequences
                        print OUT $hit->length . "\t";
                        # get the description of the hit sequences
			print OUT $hit->description . "\t";
                        # get the E value of the hit
			print OUT $hit->significance . "\t";
                        #get the bit score of the hit
			print OUT $hit->bits . "\t";

                        my $hspcount = 0;

			my $tmp = $result->query_name."\t".$result->query_length."\t".$hit->accession."\t".$hit->length."\t".$hit->description."\t".$hit->significance."\t".$hit->bits."\t";
                        # process the top HSP for the top hit
			while (my $hsp = $hit->next_hsp) {
                        	print OUT $tmp, if ($hspcount > 0);
                        	# get the frame of the query sequence
				print OUT $hsp->query->frame . "\t";
                                # get the start and the end of the query sequence in the alignment
				print OUT $hsp->start('query') . "\t" . $hsp->end('query'). "\t";
                                # get the start and the end of the hit sequence in the alignment
				print OUT $hsp->start('hit') . "\t" . $hsp->end('hit') . "\t";
                                # get the identity value
				printf OUT "%.1f" , ($hsp->frac_identical * 100);
				print OUT "%\t";
				printf OUT "%.1f" , ((($hsp->end('query') - $hsp->start('query'))/$result->query_length) * 100); 
		       		print OUT "%\n";
                                $hspcount++;
				last if($hspcount == 1);
                        }
			$count++;

                        # flow control for the number of hits needed
			last if ($count == $numHits);
		}
    	}
}
close OUT;
print " DONE!!! -- http://www.bioinformatics-made-simple.com--\n";


