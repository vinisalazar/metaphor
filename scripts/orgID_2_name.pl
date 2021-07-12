#!/usr/bin/perl
######################################3
##
## converts output from diamond to organism counts
##
## makes a massive hash table and then picks out the necessary taxonomic lineage
##
##
##
##########################################
use 5.010;
use strict;
use warnings;
use Data::Dumper;
use XML::LibXML;
use File::Glob;
use File::Basename;
use Cwd qw(getcwd);

#my $orgID2name = shift;
my $tax_rank = shift;
my $full_lineage_dump = shift;
my $outputDir = shift;
my $domFiles = shift;

#chdir($outputDir);
#my @diamondTabOutputs = split(';', $domFiles);
my @diamondTabOutputs = glob("$outputDir/*.domTable.txt");

my $array_size=@diamondTabOutputs;
if($array_size==0) {print Dumper(@diamondTabOutputs); die print "Can't find diamond tab output: $array_size $outputDir\n"};

#global hash variables
my %orgHash;
my %countHash;
my @allSpeciesArray;
my $sampleTitle="Name\tKegg Code\tLineage\t";
my %lineageHash;
my %lineageMap;

#fille out hashes for use
&create_lineage_hash($full_lineage_dump);
&create_lineage_map($tax_rank);

#cuts diamond blast output to get the hits
foreach my $diamond (@diamondTabOutputs){
	my @sampleDiamondOutput = `cut -f 2 $diamond | cut -d ':' -f 1 | sort | uniq -c`;

	#get base name
	my $base=basename($diamond);
	my @splitFile = split('\.', $base);
	my $sampleName = $splitFile[0]."_".$splitFile[1];
	$sampleTitle = $sampleTitle.$sampleName."\t";


	foreach my $sample (@sampleDiamondOutput){

		my $sampleResults = &cleanOutput($sample);

		my @splitCutResults = split(' ', $sampleResults);
		my $speciesID = $splitCutResults[1];
		my $speciesCount = $splitCutResults[0];
		#my $speciesName = &returnName($speciesID);

		#push all species into array to make a unique array containing ALL species present
		push (@allSpeciesArray, $speciesID);		

		#print Dumper(@allSpeciesArray);
		#$sampleResults =~ s/^\s//g;
		#does this just mask empty names? Find out why it isn't returning a result
		if(defined($speciesID)){
			$countHash{$sampleName}{$speciesID}=$speciesCount;
			#print "$speciesID\t$splitCutResults[0]\n";
		}
	}
}

#proves data is there
#print Dumper(%countHash);

@allSpeciesArray = grep defined, @allSpeciesArray;
my @uniqSpeciesArray = &uniq(@allSpeciesArray);

#print TEst
#print Dumper(%countHash);

######print out table for all results
#go through files -> list of unique species
#add to title
print $sampleTitle."\n";
#for each sample file
my $keyFlag=0;
#go through each species
foreach my $species (@uniqSpeciesArray){
		my $lineageName;
		my $lineageString;
	foreach my $diamond (@diamondTabOutputs){
		#return base
		my $base=basename($diamond);
	        my @splitFile = split('\.', $base);
		
		### this is inefficient and doesn't work with files which have no underscore
		#must change
       	        my $baseName = $splitFile[0]."_".$splitFile[1];
	
		my $value = $countHash{$baseName}{$species};
		if(!defined($value)){$value='0';}
		
		#print "$baseName\t$species\t$value\n";
		my $linOut = &return_lineage($species);
		#split lineage output
		#why are some lineage out empty?
		if(defined($linOut)){
			my @lineageSplit = split(/_{3}/, $linOut);

			$lineageName = &cleanOutput($lineageSplit[0]);
                	$lineageString = &cleanOutput($lineageSplit[1]);

	                #print $keyFlag."\n";
	                if($keyFlag==0){
				#stylised version
        	                #print $lineageName."\t".$species."\t".$lineageString."\t".$value."\t";
				#biom version
				print $lineageName."\t".$species."\t".$lineageString."\t$value\t";
	                }else{
        	                print $value."\t";
	                }
        	        $keyFlag++;
		}else{
			#print $value."\t";
		}

		#Initially the code above was here
		#Find out why Some $linOut is empty
		# my @lineageSplit = split(/_{3}/, $linOut);

                 #       my $lineageName = &cleanOutput($lineageSplit[0]);
                #        my $lineageString = &cleanOutput($lineageSplit[1]);

                        #print $keyFlag."\n";
                 #       if($keyFlag==0){
                 #               print $lineageName."\t".$species."\t".$lineageString."\t".$value."\t";
                 #       }else{
                 #               print $value."\t";
                 #       }
                 #       $keyFlag++;

		 #end comment
		 #cut

	} 
	print "\n";
	$keyFlag=0;
}


################################################################
#Fucntions
#


#return uniq array
sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}

#return lineage
sub return_lineage{
	my $ID = shift;

	my $lineageOutput = $lineageMap{$ID};
	return $lineageOutput;
}

#creates hash from ncbi dump
sub create_lineage_hash{
	my $ncbi_input = shift;

	open my $FH, "<:encoding(utf8)", $ncbi_input or die "$ncbi_input: $!";
        while (my $line = <$FH>) {
                chomp($line);
                my @splitOrgFile = split(/\|/, $line);
                #print Dumper(@splitOrgFile);

                #clean split output
                my $ncbiKey = &cleanOutput($splitOrgFile[0]);
                my $root = &cleanOutput($splitOrgFile[1]);
                my $lineage_value = &cleanOutput($splitOrgFile[2]);

                $lineageHash{$ncbiKey}=$lineage_value;
	}
}

#merges tax anomic rank and ncbi full lineage table
sub create_lineage_map{
	my $tax_rank_input =shift;

	open my $TR, "<:encoding(utf8)", $tax_rank_input or die "$tax_rank_input: $!";
        while (my $line = <$TR>) {
                chomp($line);
                my @splitRankFile = split(/\t/, $line);
                #print Dumper(@splitRankFile);

                #clean split output
                my $keggID = &cleanOutput($splitRankFile[0]);
                my $ncbi_key = &cleanOutput($splitRankFile[1]);
                my $name = &cleanOutput($splitRankFile[5]);
                my $lineage = $lineageHash{$ncbi_key};
		my $lineage_name="$name"."___"."$lineage";

		$lineageMap{$keggID}=$lineage_name;
	}

}


sub returnName {
	my $id = shift;
	
	my $name = $orgHash{$id};
	return $name;
}


## simple cleaning function for removing empty spaces
sub cleanOutput {

	
	my $input = shift;

	#perfform clean if input is defined
	#find out why it isn't...
	if(defined($input)){

		#remove spaces at beginning
		$input =~ s/^\s+//;
		#remove spaces at end
		$input =~ s/\s$//;

		#remove new line at end
		$input =~ s/\R\z//;
		return $input;
	}
}



#print Dumper(%orgHash);
