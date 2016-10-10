## Attention, nouvelle version a partir des alignements sur NM_000546.5 pour bwa, gmap et blasr (commenter en fonction)
## Final table format : 
	# IDMM_PacBio.tsv : [Run Sample Fragment Read RefPosition ReadPosition Error_type Size Alteration]  
	# Coverage_PacBio.tsv : [Sample Fragment Coverage]  
# Just reads whith a MAPQ => 3
use warnings;
use strict;
use List::Util qw(sum);

## For bwa alignments
# my $dir_sam = "/home/jrudewicz/Bureau/TP53_MICADov2/alignment/bwa/sam/" ;
# my $dir_analyses = "/home/jrudewicz/Bureau/PACBiofiguresThese/AliStat/bwa" ;
# my $blasr = "False";

## For gmap alignments
my $dir_sam = "/home/jrudewicz/Bureau/TP53_MICADov2/alignment/gmap/sam/" ;
my $dir_analyses = "/home/jrudewicz/Bureau/PACBiofiguresThese/AliStat/gmap" ;
my $blasr = "False";

## For blasr alignments
# my $dir_sam = "/home/jrudewicz/Bureau/TP53_MICADov2/alignment/blasr/sam/" ;
# my $dir_analyses = "/home/jrudewicz/Bureau/PACBiofiguresThese/AliStat/blasr" ;
# my $blasr = "True";

my(

$fragment,
$sample,
$comp_reads,
$run,
$biopsy,
$read,
$nb,
$seq,
$qual,
$condition_SoftClip,
$nbSeqNt,
$position,
$SoftClip,
$readName,
@split_readName,
@cig,
@SeqNt,
@Quality,
@samfile,
@ligne,
@split_MD,
@ASCII,
%ASCII

);

@ASCII = split("","!\"#\$\%&'()*+,-./0123456789:;<=>?\@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~") ;

for (my $i=0; $i<=$#ASCII; $i++)
{

	$ASCII{$ASCII[$i]} = $i ;

}


@samfile = `ls $dir_sam`;
open ( ID , ">$dir_analyses/IDMM_PB.tsv") ;
open ( COVERAGE , ">$dir_analyses/Coverage_PB.tsv") ;
print ID "Run\tSample\tBiopsy\tFragment\tRead\tPosition\tReadPosition\tError_type\tSize\tAlteration\tQuality\n" ;
print COVERAGE "Sample\tBiopsy\tFragment\tCoverage\n" ;

foreach my $sam (@samfile)
{

	chomp($sam);
	print "$sam\n";
	$sam =~ m/([NC])_(.+)_(\d+)/ ;
	$fragment = $1 ;
	$sample = $2 ;
	$biopsy = $3 ;
	$comp_reads = 0 ;

	open ( SAM , "$dir_sam/$sam") ;

	while (<SAM>)
	{		
					
		if ($_ !~ m/^@/ and $_ ne "")
		{
											
			@ligne = split ( "\t" , $_ ) ;	

			next if $ligne[4] < 3 ; # Not sup 30 because GMAP seems to have bug for MAPQ [...]
			#next if $ligne[2] ne "17" ; 

			$comp_reads ++ ;

			if ($ligne[0] =~ m/(.+)\/(.+)\/ccs/)
			{
			
				$run = $1 ;
				$read = $2 ;

			}
			elsif ($ligne[0] =~ m/(.+)\/(.+)/ )
			{
				$readName = $1 ;
				if ($blasr eq "True") 
				{
					@split_readName = split("/", $readName) ;
					$run = $split_readName[0] ;
					$read = $split_readName[1] ;
				}
				else
				{
					$run = $1 ;
					$read = $2 ;
				}
			}
			else # Au cas où
			{

				print "Erreur pour le nom des read et run\n" ;

			}

			#print "SAMPLE $sample FRAGMENT $fragment READ $read\n";

			$condition_SoftClip = 0 ;
			$nbSeqNt = 0 ; # position on the read
			$position = $ligne[3]-1 ; # position on GRCh37
			$SoftClip = 0 ;
			@cig = split("", $ligne[5]) ;
			@SeqNt = split("", $ligne[9]) ;
			@Quality = split("", $ligne[10]) ;
			my @compo_cigard ;
			my @SeqNt_wHoutI ;
			my @SeqQual_wHoutI ;

			foreach  (my $Q = 0 ; $Q <= $#Quality ; $Q ++ )
			{

				$Quality[$Q] = $ASCII{$Quality[$Q]} ;

			}

			for (my $case_CIGARD = 0 ; $case_CIGARD <= $#cig ; $case_CIGARD ++ )
			{
			
				$nb = "" ;
											
				if ($cig[$case_CIGARD] !~ m/(\d)/)
				{
					
					next if $cig[$case_CIGARD] eq "H" ; 

					for (my $case_before = $case_CIGARD-1 ; $case_before >= 0 ; $case_before -- )
					{
														
						last if ($cig[$case_before] !~ /\d/) ;
							
						$nb = "$cig[$case_before]$nb" ;
														
					}

					
					if ($cig[$case_CIGARD] eq "S" and $condition_SoftClip == 0)
					{

						$nbSeqNt += $nb ; 
						$SoftClip = $nb ; 
						$condition_SoftClip = 1 ;

					}
					elsif ($cig[$case_CIGARD] eq "N")
					{

						$position += $nb ; 
						push (@compo_cigard,$nb,"N");

					}
					elsif ($cig[$case_CIGARD] eq "I")
					{

						$seq = join ("", @SeqNt[($nbSeqNt)..($nbSeqNt+$nb-1)]) ;
						
						$qual = mean(@Quality[($nbSeqNt)..($nbSeqNt+$nb-1)]) ;
						print ID "$run\t$sample\t$biopsy\t$fragment\t$read\t$position\t$nbSeqNt\t$cig[$case_CIGARD]\t$nb\t$seq\t$qual\n" ;
						$nbSeqNt += $nb ; 

					}
					else
					{

						push (@SeqNt_wHoutI,@SeqNt[($nbSeqNt)..($nbSeqNt+$nb-1)]) if (($cig[$case_CIGARD] eq "M") or ($cig[$case_CIGARD] eq"=") or ($cig[$case_CIGARD] eq "X"));
						push (@SeqQual_wHoutI,@Quality[($nbSeqNt)..($nbSeqNt+$nb-1)]) if (($cig[$case_CIGARD] eq "M") or ($cig[$case_CIGARD] eq"=") or ($cig[$case_CIGARD] eq "X"));
						$position += $nb ; 
						$nbSeqNt += $nb if (($cig[$case_CIGARD] eq "M") or ($cig[$case_CIGARD] eq"=") or ($cig[$case_CIGARD] eq "X")); 
						push (@compo_cigard,$nb,"M") if (($cig[$case_CIGARD] eq "M") or ($cig[$case_CIGARD] eq"=") or ($cig[$case_CIGARD] eq "X"));
						push (@compo_cigard,$nb,"D") if $cig[$case_CIGARD] eq "D"; 


					}
				}
			}

			($ligne[12] = $ligne[18]) if $blasr eq "True";
			$ligne[12] =~ m/MD:Z:(.+)/;
			$ligne[12] = $1 ;

 			if ($ligne[12] =~ m/[ATCG\^]/)
 			{

				@split_MD = split ("", $ligne[12]) ;
				$nbSeqNt = 0 ; 
				my $compt = 0;
				$position = $ligne[3] ;  

				for (my $case_MD = 0 ; $case_MD <= $#split_MD ; $case_MD ++ )
				{

					$nb = "" ;


					if ($split_MD[$case_MD] !~ m/(\d)/)
					{
					
						for (my $case_before = $case_MD-1 ; $case_before >= 0 ; $case_before -- )
						{
														
							last if ($split_MD[$case_before] !~ /\d/) ;
							
							$nb = "$split_MD[$case_before]$nb" ;

						}

						$nb = 0 if $nb eq "" ;
						my $nb2 = $nb ;
						$nbSeqNt += $nb ;

						for (my $case_compo = $compt ; $case_compo <= $#compo_cigard ; $case_compo += 2)
						{

							next if $compo_cigard[$case_compo+1] eq "D" ; 

							if ($compo_cigard[$case_compo+1] eq "N")
							{

								$position += $compo_cigard[$case_compo] ;

							}
							elsif ($compo_cigard[$case_compo] < ($nb2)) 
							{

								$position += $compo_cigard[$case_compo] ;
								$nb2 -= $compo_cigard[$case_compo] ;

							}
							else
							{

								$position += $nb2 ;
								$compo_cigard[$case_compo] -= $nb2 ; 
								$compt = $case_compo ;
								($position += $compo_cigard[$case_compo+2] and $compt = ($case_compo + 4) ) if ($compo_cigard[$case_compo] == 0 and $compo_cigard[$case_compo+3] eq "N") ;
								last ;
							

							}
						}

						if ($split_MD[$case_MD] eq "^")
						{

							$nb = 0 ;
							$seq = "" ;

							while ($split_MD[$case_MD+1] =~ m/[ATCG]/)
						 	{

								$seq .= $split_MD[$case_MD+1] ;							
								$case_MD ++ ;
								$nb ++ ;
																						
							}		
							
							$qual = mean(@Quality[($nbSeqNt-2)..($nbSeqNt+1)]) ;
							print ID "$run\t$sample\t$biopsy\t$fragment\t$read\t$position\t$nbSeqNt\tD\t$nb\t$seq\t$qual\n" ;
							$position += $nb ;

						}
						elsif ($split_MD[$case_MD] =~ m/[ATCG]/)
						{

							($position += $compo_cigard[$compt+4] and $compt = ($compt + 6) ) if ($compo_cigard[$compt] == 0 and $compo_cigard[$compt+3] eq "D" and $compo_cigard[$compt+5] eq "N") ;
							$nbSeqNt ++; 
							print ID "$run\t$sample\t$biopsy\t$fragment\t$read\t$position\t$nbSeqNt\tM\t1\t$SeqNt_wHoutI[$nbSeqNt-1]\t$SeqQual_wHoutI[$nbSeqNt-1]\n" ;
							$compo_cigard[$compt] -= 1 ;							
							$position ++ ;

						}												
					} 
				}
			}
		}		
	}									
	close SAM ;

	print COVERAGE "$sample\t$biopsy\t$fragment\t$comp_reads\n" ;

}
close ID ;

sub mean {
    return sum(@_)/@_;
}
