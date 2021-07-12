#!/usr/bin/perl

use 5.010;
use strict;
use warnings;
use Data::Dumper;
use XML::LibXML;
use File::Glob;
use File::Basename;

#inputs files for various reasons
#may change to parameters in a bit
my $outputDir = shift;
my $outputFileName = shift;
my $analysisType = shift;
my $koFormatFile = shift;
my $speciesFile = shift;
my $xmlFiles = shift;

print "$xmlFiles\n";

my @xmlArray = split(';', $xmlFiles); 

my @geneArray; 
my @sampleArray;
my @geneKeyArray;

#hashes to be used as references
#Significantly inceases speed
my %resultsHash;
my %keggHash;
my %speciesHash;
my %keggToFunction;


#create hash for tables for quick usage
&create_species_hash;

#main function
&main;

#create kegg hash and ko list
&create_kegg_hash;

#use gene count table to create functional level counts
#uses the dome tables
&function_count;

#gets the dom tables from the output folder and creates function levels
#half the function is repeated. Object orient
#Update: quite a bit different to originally thought
sub function_count{
	 #glob all output files and create complete table entries from all files
	#print Dumper(%speciesHash);
        my @outputFiles = glob("$outputDir/*.domTable.txt");

	if(!@outputFiles){die "no dom tables in output directory\n";}	

	my @sampleNameArray;

	#create mega level hash
	my %megaLevelHash;

        foreach my $outFile (@outputFiles){
	my @level1Array;
	my @level2Array;
	my @level3Array;
	my @level4Array;

                my @geneList = `cut -f 2 $outFile | cut -d ' '  -f $analysisType`;
                #print Dumper(@geneList);
                foreach my $gene (@geneList){
                        chomp($gene);
                        if($gene ne "Kegg ID"){
                                if(defined($gene)){
					my @levelArray=%speciesHash{$gene};
					if(defined($levelArray[0])){
						#push(@geneArray, $levelArray[1][1]);
						my $kid = $levelArray[1][1];		
						#print Dumper(%keggHash);
						#Check to see if it's getting all the counts properly!
						if(defined($kid)){
							my @levOne = %keggHash{$kid};
							if(defined($levOne[1])){
								my @countArray = $levOne[1];
								#print Dumper (@countArray);
								
								#this dumper is important! It is showing that all hits are being printed!
								#print Dumper( @{$levOne[1]} );
								my $arrayCount=0;
								foreach my $test( @{$levOne[1]} ){
									#print $test."\n";
									$arrayCount++;
								}
								#print $arrayCount."\n";
								#divide Array count by 5 to get the number of iterations of the megeahash creation
								#Of course coudl be a more elegant way of doing this.
								my $iterations = $arrayCount/5;
								my $count;
								for ($count=0;$count<$arrayCount;$count+=5){
									#print $levOne[1][$count+1]."\n";
									#print $levOne[1][$count+2]."\n";
									#print $levOne[1][$count+3]."\n";
									#print "$levOne[1][$count+1]\t $levOne[1][$count+2]\t$levOne[1][$count+3]";

									#cycle through megaHashes
									my $levelOneString = $levOne[1][1];
                                                               		my $levelTwoString = $levOne[1][2];
                                                                	my $levelThreeString = $levOne[1][3];
									my $levelFourString;
                                                                   	
									if(defined($levelOneString)){
                                                                       	 	$levelOneString =~ s/^\s+|\s+$//g;
                                                                        	push(@level1Array, $levelOneString);
                                                                  	 }
                                                                   	if(defined($levelTwoString)){
                                                                     		$levelTwoString =~ s/^\s+|\s+$//g;
                                                                        	push(@level2Array, $levelTwoString);
                                                                   	 }
                                                                   	if(defined($levelThreeString)){
                                                                      		$levelThreeString =~ s/^\s+|\s+$//g;
                                                                        	push(@level3Array, $levelThreeString);
                                                                   	}
									if(defined($levelOneString) && defined($levelTwoString) && defined($levelThreeString)){
			
										#split last column into kegg and function
										my @functionSplit = split('  ', $levOne[1][$count+4]);
										#print Dumper(@functionSplit);
										#my @functionSplitOnDoubleSpace = split('  ', $functionSplit[1]);
	
										my $function = $functionSplit[2];										
										$function =~ s/^\s+|\s+$//g;

										$levelFourString = "$levOne[1][$count+1]\t $levOne[1][$count+2]\t$levOne[1][$count+3]\t$function";
										$levelFourString =~ s/^\s+|\s+$//g;
										push(@level4Array, $levelFourString);
									}
								}
								#print $arrayCount."\n";
								#print "-----\n";
								$arrayCount='0';
							}
						}
					}
				}
                        }
                }
		#data dump outputs
		#create sample name

		my ($base, $dir, $ext) = basename($outFile);

		#my @sampleName=split('\/', $outFile);
		my @sampleNameBase=split('\.', $base);
		my $sampleNamePrefix=$sampleNameBase[0];
		#print Dumper(@sampleName);
		push(@sampleNameArray, $sampleNamePrefix);

		#create megahash for table	
		my %level1Counts;
        	$level1Counts{$_}++ for @level1Array;

		my %level2Counts;
                $level2Counts{$_}++ for @level2Array;

		my %level3Counts;
                $level3Counts{$_}++ for @level3Array;
	
		my %level4Counts;
		$level4Counts{$_}++ for @level4Array;

		#store hash counts in mega hash	
		$megaLevelHash{$sampleNamePrefix}{"level1"}=\%level1Counts;
		$megaLevelHash{$sampleNamePrefix}{"level2"}=\%level2Counts;
		$megaLevelHash{$sampleNamePrefix}{"level3"}=\%level3Counts;
		$megaLevelHash{$sampleNamePrefix}{"level4"}=\%level4Counts;
        }
		#pass hash and key array to print level table function
		#print Dumper(%megaLevelHash);
		&print_level_table(\%megaLevelHash, \@sampleNameArray);
		&print_all_table(\%megaLevelHash, \@sampleNameArray);
}

#print table with all the levels for a gene on one line
sub print_all_table{
	 my $megaLevelHash_sub = shift;  my $sampleNameArray_sub = shift; my $hashLevel; my %tempHash; my @keyArray; my @sArray; my @roleArray; my %temp;

	 #set necessary data structures
	 my @levelArray=('level1', 'level2', 'level3', 'level4');
         @sArray=&ret_sample_array($sampleNameArray_sub);
	 my %temp_all_level_table_hash=%{$megaLevelHash_sub};
	 #print Dumper(%temp_all_level_table_hash);

         my $file = "Functional.table.counts.txt";
         open(my $fileHandle, '>', $outputDir."/".$file) or die "Could not open file '$file' $!";

	#foreach my $hashLevel (@levelArray){
	my @roleArray_all_table = &ret_role_array($sampleNameArray_sub, $megaLevelHash_sub, "level4");
	my @uniqArray=&uniq(@roleArray_all_table);

                my $title = "Sample\t";
                foreach my $sub (@{$sampleNameArray_sub}){
                        $title=$title.$sub."\t";
                }
                print $fileHandle $title."\n";

			#print Dumper(@uniqArray);
                        #need different order to print table
                        foreach my $r (@uniqArray){
                                #print "$r:\t";
                                print $fileHandle $r."\t";
                                foreach my $s (@sArray){
					#go through each hashlevel and print the brite levels first
					my $levelString="";
					foreach my $hashLevelInner (@levelArray){
                                        	my ($levelKey, $levelValue) = $temp_all_level_table_hash{$s}{$hashLevelInner}{$r};
						if(defined($levelKey)){$levelString = $levelString."\t".$r;}
					}
						my $value= $temp_all_level_table_hash{$s}{"level4"}{$r};
                                        	if(!defined($value)){$value ='0';}
						####print $levelString."\t".$value."\n";
                                        	print $fileHandle "$value\t";
                                }
                                print $fileHandle "\n";
                        }
                        print $fileHandle "\n";

                        #reset arrays for next leve     l
                        undef(@uniqArray);
                        undef(@sArray);
                        undef(%temp_all_level_table_hash);
                        undef(@roleArray);
	#}
}

#return all the roles in one array
sub ret_role_array{
	my $role_temp_array = shift;
	my $megaLevelHash_sub = shift;
	my $level = shift;
	my @roleArray;
	my @sArray;

	     foreach my $samp (@{$role_temp_array}){
                   push(@sArray, $samp);
                   my %temp = %{$megaLevelHash_sub};
                   my @temp2 = $temp{$samp}{$level};

                   for my $href ( @temp2 ) {
                          for my $role ( keys %$href ) {
                                 push(@roleArray, $role);
                          }
                  }
             }

	return @roleArray;
}

#function for returning sArray
#sArray is empty below
sub ret_sample_array{
	 my $role_temp_array = shift;
	 my @sArray;

	foreach my $samp (@{$role_temp_array}){
		 push(@sArray, $samp);
	}
	return @sArray;
}


#prints out level tables
sub print_level_table{	
	my $megaLevelHash_sub = shift; 	my $sampleNameArray_sub = shift; my %tempHash; my @keyArray; my @sArray; my @roleArray; my %temp;
		
	%temp=%{$megaLevelHash_sub};
	#print Dumper(%temp);
	#create level Array, levels will always stay the same
	my @levelArray=('level1', 'level2', 'level3');
	@sArray=&ret_sample_array($sampleNameArray_sub);

	print Dumper(@levelArray);

		foreach my $hashLevel (@levelArray){
		#open files for printing 
		#samplename => level => key #go trhough samples
                #need to perform two loops to get keys from hash
                #may be a more elegant way, will look into it later
		#print $hashLevel."\n";
                @roleArray=&ret_role_array($sampleNameArray_sub, $megaLevelHash_sub, $hashLevel);
                #remove duplicates from array

                #print Dumper (@roleArray);
                my @uniqArray=&uniq(@roleArray);
		#print Dumper(@uniqArray);
		
		print "$hashLevel\n";
		my $file = "$hashLevel.counts.txt";
		open(my $fileHandle, '>', $outputDir."/".$file) or die "Could not open file '$file' $!";	

		my $title = "Sample\t";
                foreach my $sub (@{$sampleNameArray_sub}){
                        $title=$title.$sub."\t";
                }
                print $fileHandle $title."\n";
	
			#need different order to print table
			foreach my $r (@uniqArray){
				#print "$r:\t";
				print $fileHandle $r."\t";
				foreach my $s (@sArray){
					my $value= $temp{$s}{$hashLevel}{$r};
					if(!defined($value)){$value ='0';}
					print $fileHandle "$value\t";
				}
				print $fileHandle "\n";
			}
			print $fileHandle "\n";

			#reset arrays for next leve	l
			undef(@uniqArray);
			
	        }
			undef(@sArray);
			undef(%temp);
			undef(@roleArray);
}

sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}


sub create_species_hash{
	open my $speciesFile, "<:encoding(utf8)", $speciesFile or die "$speciesFile: $!";
        while (my $line = <$speciesFile>) {
                chomp $line;
                my @splitSpeciesFile = split('\t', $line);
                push(@{ $speciesHash{$splitSpeciesFile[0]} }, @splitSpeciesFile);
        }
}

sub create_kegg_hash{

	open my $koFormatFile, "<:encoding(utf8)", $koFormatFile or die "$koFormatFile: $!";
        while (my $line = <$koFormatFile>) {
                chomp $line;
		my @splitKoFile = split('\t|\|', $line);
		push(@{ $keggHash{$splitKoFile[0]} }, @splitKoFile);
        }

}

sub main{

# find all xml files in directory

	if(!@xmlArray){ die "There are no xml files currently in the folder for processing\n"; }

	
	foreach my $xmlFile (@xmlArray){

	my $dom = XML::LibXML->load_xml(location => $xmlFile);
	my $xmlBase = basename($xmlFile);
	open(my $fh, ">", "$outputDir/$xmlBase.domTable.txt") or die "Can't open $xmlFile.domTable.txt for writing: debug - main fucntion\n";
	print $fh "Query ID\tKegg ID\tDescription of Hit\n";

		#takes xml output from diamond and parses into simple Query - ID - description file
		foreach my $hit ($dom->findnodes('/BlastOutput/BlastOutput_iterations/Iteration')) {
			if($hit->findvalue('./Iteration_query-ID') =~ m/Query/){
				print $fh $hit->findvalue('./Iteration_query-ID')."\t";
				ID:
				foreach my $id ($hit->findnodes('./Iteration_hits/Hit/Hit_id')) { if(defined($id->textContent)) { print $fh $id->textContent; print "\t"; last ID; } else {print "\n"; }}
	
				DEF:
				foreach my $def ($hit->findnodes('./Iteration_hits/Hit/Hit_def')){ if(defined($def->textContent)) { print $fh $def->textContent; last DEF; } else { print "\n"; }}
			
				print $fh "\n";	
			}else{
				print $fh "\n";
			} # end else
	
		}#end foreach
	close($fh);
	}

	#glob all output files and create complete table entries from all files
	my @outputFiles = glob("$outputDir/*domTable.txt");

	foreach my $outFile (@outputFiles){
		my @geneList = `cut -f 2 $outFile | cut -d ' '  -f $analysisType`;
		#print Dumper(@geneList);
		foreach my $gene (@geneList){
			chomp($gene);
			if($gene ne "Kegg ID"){
				push(@geneArray, $gene);
			}
		}
	}
	
	#remove duplicates from the array
	my @tmpArray=uniq(@geneArray);
	@geneArray=@tmpArray;
	
	#simplest way to do this is to take each dom file and sort 
	#and uniq and check to see if the gene is in there with 
	#the gene copy nomber. Place in a hash and extract using key
	foreach my $outFile2 (@outputFiles) {
		my @geneCopyNumber = `cut -f 2 $outFile2 | cut -d ' ' -f $analysisType | sort | uniq -c`;
			#take base of filename
			my($name,$path,$suffix)=fileparse($outFile2);
			my @fileNameSplit=split(/\./, $name);
                        my $basename=$fileNameSplit[0];

			#push names and keys into arrays for later
			push(@sampleArray, $basename);
			push(@geneKeyArray, $name);
			
			#chomp newline from end of string
			#remove spaces from the start of string
			#split into two parts and place into hash
			foreach my $gcn (@geneCopyNumber){
				chomp $gcn;
				$gcn =~ s/^\s+//;

				my @splitGCN=split(' ', $gcn);
				#hash by file name then gene and gene count number
				if(!defined($splitGCN[1])){$splitGCN[1]="Empty";}
				$resultsHash{$basename}{$splitGCN[1]}=$splitGCN[0];
			}
	}
		#print Dumper(%resultsHash);

	#Print out the table!!!!!!!
	my $geneCountTable=$outputFileName;
	
 	open(my $geneCountTableFH, ">", $outputDir."/".$geneCountTable) or die "CAn't open $geneCountTable for writing\n";

	 #headers
	 print $geneCountTableFH "Gene/Sample\t";
	
	 foreach my $samps (@sampleArray){
		print $geneCountTableFH $samps."\t";
	 }

	 #quicikly sort array according to alphabetic order
	 my @tmpSortArray = sort { lc($a) cmp lc($b) } @geneArray;
	 @geneArray = @tmpSortArray;

	 print $geneCountTableFH "\n";
	 	foreach my $geneKey (@geneArray){
			if($geneKey ne "" && $geneKey ne "Kegg"){
				print $geneCountTableFH "$geneKey\t";
				foreach my $sampleKey (@sampleArray){
					my $value=$resultsHash{$sampleKey}{$geneKey};
					if (defined($value)){print $geneCountTableFH "$value\t"}else{print $geneCountTableFH "0\t";}	
				}
				print $geneCountTableFH "\n";
			}
	 	} 	
} # end main
