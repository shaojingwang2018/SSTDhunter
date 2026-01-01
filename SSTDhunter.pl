#/usr/bin/perl -w
use strict;
use warnings;

if(!defined($ARGV[0])){
die " ___ ___ ___ ___ ___     ___ ___ ___ ___ ___
|___|___|___|___|___|^-^|___|___|___|___|___|

SSTDhunter: 
	A tool for investigating Steroid-17,20-desmolase (SSTD) encoding genes in metagenomic data.
	
Version: 
	v1.0

Usage:
	SSTDhunter.pl [prefix].sam\n
	
Output: 
	[prefix]_SSTDhunter_abundance.txt
	[prefix]_SSTDhunter_reads.txt
	

Recommended commands for Bowtie2
	\$ bowtie2 -p 30 -x SSTDhunterDB -U [prefix].fq.gz --reorder -S [prefix].sam\n
";
}

# Define the necessary variables
my(@mid, $i, $j, $k);
my(@table, %hash, %hashl);
my $total = 0;
my $prefix;

#Define prefix
@mid = split(/[\\\/]/, $ARGV[0]);
$prefix = $mid[$#mid];
$prefix =~ s/.sam$//g;

#Extract information from the sam file 
open(HEHE,"<".$ARGV[0]) or die $!;
open(OUT,">".$prefix."_SSTDhunter_reads.txt") or die $!;
while(<HEHE>){
	$_ =~ s/[\r\n]+//g;
	if($_ =~ /^@..\t/){
		# nt length of reference genes
		if($_ =~ /^\@SQ\t.+/){
			@mid = split(/\t/, $_);
			if($mid[1] =~ /^SN:.+/ && $mid[2] =~ /^LN:.+/){
				$mid[1] =~ s/^.+_//;
				$mid[1] =~ s/\d+$//;
				$k = $mid[1];
				
				$mid[2] =~ s/^LN://;
				if(!defined($hashl{$k})){
					$hashl{$k} = $mid[2];
				}
				elsif($hashl{$k} < $mid[2]){
					$hashl{$k} = $mid[2];
				}
			}
			else{
				die "\n".'[Error].#2.Please report bugs to shaojingwang@tmu.edu.cn. thanks!'."\n\n";
			}
		}
		# reference information
		# do nothing 
	}
	else{
		# Total data accumulation
		$total ++;
		
		# Calculate row by row 
		@mid = split(/\t/, $_);
		
		# Identifying matching reads
		if($mid[5] =~ /M/){
			# Verify the target genes
			if($mid[2] =~ /des/){
				$mid[2] =~ s/^.+_//;
				$mid[2] =~ s/\d+$//;
				$k = $mid[2];
				
				# Mapped data accumulation for each gene
				if(!defined($hash{$k})){
					push @table, [$k];
					$hash{$k} = 1;
				}
				else{
					$hash{$k} ++;
				}
			}
			else{
				die "\n[Error].#1.Wrong database used, please check the SSTDhunter database used during bowtie2 alignment.\n\n";
			}
			
			# Generage temporary file of target sequences in fasta format
			print OUT ">".$mid[0]." ".$mid[2]." ".$mid[5]."\n";
			print OUT $mid[9]."\n";
		}
	}
}
close OUT;
close HEHE;

@table = sort{$a->[0] cmp $b->[0]} @table;

for $i(0 .. $#table){
	#Assignment of gene length in kilobases
	$table[$i][1] = $hashl{$table[$i][0]}/1000;
	
	#Assignment of reads mapped to the gene
	$table[$i][2] = $hash{$table[$i][0]};
	
	#Total in millions
	$table[$i][3] = $total/1000000;
	
	#Calculation of RPKM 
	$table[$i][4] = $table[$i][2]/($table[$i][1]*$table[$i][3]);
}

open(OUT,">".$prefix."_SSTDhunter_abundance.txt") or die $!;
print OUT "gene_name\tRPKM\n";
for $i(0 .. $#table){
	print OUT $table[$i][0]."\t";
	print OUT $table[$i][4]."\n";
}
close OUT;




