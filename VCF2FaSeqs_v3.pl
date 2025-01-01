#!/usr/bin/env perl
#This takes a vcf file (plain vcf or bcf) and generate sequence or pseudosequence
#to represent the genotype of the sames at 1 seq/sample for phylogeny or
#sequence similarity search, i.e. wgs-based genotyping
#This version tries to save memory by combining 4 individual hashes (AF,GF,REFALT,POS) into a multi-key hashes %VAR and removes variants not passing the filters
#Author: Ping Liang, 2023-4-23
#This version differs from v2 in using $SEQ{$pos,$sID} vs $SEQ{$pos}{$sID} for different memory usage, and process the variants in batches

use strict;
use File::Basename;
use Getopt::Std;

my (%opt);
getopts("h:v:b:f:r:c:p:s:a:o:F:a:l:P:I:",\%opt);
if ($opt{h}||!$opt{v}){
    print qq(
    This script is for running gt for polymorphic ref-ME, requiring samtools and blat in PATH
    Usage: $0 [options] 
    options: 
    -v: input file in vcf format, required
    -c: string, chrID, can be used to split variants by chromosome, optional
    -p: string, position range, e.g., 1-100, optional
    -P: string, file name for position list to generate only sequences for the listed position, ignoring the other filters
        this is used to generate a genotyping sequence for a new sample for the intended database, optional
    -s: int, sequence conversion option
       1 for SNPs as default 
       2 for non-SNPs
       3 hybrid, for all variants: SNP as in #1 and others as in #2
          1-3: 2bp/loci, 4-6 1bp/loci 
       4 by genotype, for all variants: 0/0->C,0/1->A,1/1->T, ./.->N
       5 by genotype, for SNPs: 0/0->C,0/1->A,1/1->T, ./.->N
       6 by genotype, for non-SNPs: 0/0->C,0/1->A,1/1->T, ./.->N
    
    -a: int, minimal allele frequency to retain, default: any, optional
    -o: string, pre-fix for output, default, [vcfFile].fa, optional
    -b: int, batch size for variant process, default is 10,000
    -f: string, PASS, filtering string for quality column, optional
    -l: 1 default to output a variant pos file [outfile.bed] in bed: chr\tpos pos\tref\talt; use 0 to disable or a file name;  incompatible with -P
    -F: 1, require all samples have genotype calls at the position to be retained, optional, default 0.5
    -I: "[._-]", keep the first part of the sample ID before the first specified character to remove unnessary part string section, default: none
    -r: int, 1 for simple progress, 10000 for also every 10000 lines, optional
    -h: this help menu
    );
    #to-do: allow the use of user input variant list, and calculate the corresponding sequence
    exit;
}

if ($opt{P} && $opt{l}){die "-p and -l are incompatible, use only one of them.\n"}
if (! -s $opt{v}){die "vcf file $opt{v} not exist or empty.\n"}
my $bcfview=$opt{v};
my $st=1; #sequence converting option, default is for SNPs only
if ($opt{s}){$st=$opt{s}}
if ($opt{c}){$bcfview .=" $opt{c}"}
if ($opt{p}){$bcfview .=":$opt{p}"}
open (VCF, "bcftools view $bcfview |") or die "$!.\n";
my $minAF="0";
my $minGF="0.5"; #default minimal ratio of genotyped samples among all
if ($opt{F}){$minGF=$opt{F}}
if ($opt{a}){$minAF=$opt{a}}
my $batch=10000;
if ($opt{b}){$batch=$opt{b}}
my ($outfile)=basename($opt{v})=~/^(.+)\./;
$outfile ="seq_${st}_$minAF.$outfile";
if ($opt{o}){$outfile=$opt{o}}
open(OUT,">$outfile.fa") or die "$!.\n";
#print STDERR "outfile: $outfile\n";

my $printloc=1;
if ($opt{l}){$printloc=$opt{l}}
my $posfile="$outfile.bed"; #by default, a corresponding bed file for variant positions is generated
if ($printloc ne "0" && !$opt{P}){
    if($printloc ne "1"){$posfile=$opt{l}} #use user speciifed file name for the position output file
    open(POS,">$posfile") or die "$!.\n";
}
#my (@SMP,%SEQ,$ln,$sn,%REFALT,%AF,%POS,%seen,@P,%GF);
my (@SMP,%SEQ,$sn,%VAR,%seen,@P,%FA);
#process the last header line to get the list of sample IDs
#if ($opt{r}){print STDERR "processing vcf header of $opt{v}.....\n"}

if ($opt{P} && ! -s $opt{P}){die "$opt{P} not found.\n"}
if ($opt{P}){#process pre-defined positoin list
    open(IN, "<$opt{P}") or die "$!\n";
    while(<IN>){
	if (!$_){next}
	my ($c,$p,$p,$ST,$tmp,$tmp,$tmp)=split /\t/;
	my $P="$c\t$p";
	if ($VAR{$P}){next}
	else{$VAR{$P}={C=>$c,P=>$p}}
	push @P, $P;
    }
    close IN;
}
my $header;
while(<VCF>){
    if ($_=~/^#CHROM/){
	$header=$_;
	last
    }
}

#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG00096.final   HG00097.final   HG00099.final   HG00100.final
chomp $header;
my @f=split /\t/, $header;
@SMP=@f[9 .. $#f];
$sn=scalar @SMP; 
foreach my $n (0 .. $#SMP){#remove any extra string in the sample ID
    my $tmp=$SMP[$n]; #($tmp)=basename($tmp)=~/^(.+?)\.{0,1}/;
    if ($opt{I} && $tmp=~/$opt{I}/){($tmp)=$tmp=~/^(.+?)[\.]/}; #keep only the part of sample ID before the first .
    $SMP[$n]=$tmp;
}
 
#processing data line
#if ($opt{r}){print STDERR "processing vcf file: $opt{v} with $sn sample(s).....\n"}
my ($ln,$gv,$snp,$snp,$ind)=(0,0,0,0); #total number of variant lines and variants kept from -a and -F
while (<VCF>){#starting from the 1st data line
    chomp;
    $ln++;
    my @f=split /\t/;
    #NC_012007.3     1019    .       C       T       79      .       DP=63;VDB=0.134913;SGB=-0.688148;RPB=0.856285;MQB=0.216462;MQSB=0.984879;BQB=0.884203;MQ0F=0.031746;ICB=1;HOB=0.5;AC=1;AN=2;DP4=26,18,5,10;MQ=32        GT:PL   0/1:113,0,255
    #NC_012007.3     1061    .       gtt     gTtt    7.99236 .       INDEL;IDV=7;IMF=0.12963;DP=54;VDB=0.00163602;SGB=-0.651104;MQSB=0.233365;MQ0F=0.0555556;ICB=1;HOB=0.5;AC=1;AN=2;DP4=18,28,6,2;MQ=32     GT:PL   0/1:42,0,214

    if ($opt{f} && $f[5]!~/$opt{f}/i){next} #skipping variants not passing the filter
    if ($opt{r} && $opt{r}>=100 && $ln % $opt{r}  == 0){print STDERR "processing variant line: $gv/$ln\r"}
    if ($f[1] !~/^\d+$/){next} #remove irregular lines
    my ($c,$p,$ref,$alt,$qual,$inf)=@f[0,1,3,4];
    if ($alt eq "."){next} #skip entries with only the reference allele (alt being the same as ref)
    my @ALT=split /\,/,$alt; #all alt alleles correspond to 1,2,3,...
    my $indel; #check to see if any alt allel represents indel, if so, the entire site is treated as an indel site
    if (length($ref)>1 || $inf=~/^INDEL/){$indel=1}#ref allele longer than alter allele
    else{
	foreach my $al (@ALT){
	    if (length($al)>1){
		$indel=1; #print "$c\t$p\t$ref\t$alt\n";
		last
	    }
	}
    }
    if (!$opt{P}){
	if ($st=~/[15]/ && ($ref!~/^[ACGT]$/i || $indel)){next} #exclude non-SNPs, even for a specified position list?
	if ($st=~/[26]/ && $ref=~/^[ACGT]$/i && !$indel){next} #exclude SNPs, unless an pre-defined variant list is specified
    }
    my $pos="$c\t$p";
    if ($opt{P} && !$VAR{$pos}->{C}){next}; #when a list of positions specified, ignore all other positions
    if (!$opt{P}){$VAR{$pos}->{C}=$c; $VAR{$pos}->{P}=$p} #save position info to %VAR
    my @GT=@f[9 .. $#f]; #grab the sample genotype columns
    my ($af,$ac,$ta,$gs)=(0,0,0,0); #allele frequency, count, total count, genotyped sample cnt
    foreach my $n (0 .. $#GT){#@GT correspond to @SMP
	my $sID=$SMP[$n];
	my ($a1,$a2)=$GT[$n]=~/^(.)\/(.)/; #numberic allele ID; igore any additional alleles in the same sample
	if ($a1 ne "."){
	    $ta++;
	    $ac++ if ($a1=~/^[1-9]$/)#all alternate alleles
	}
	if ($a2 ne "."){
	    $ta++;
	    $ac++ if ($a2=~/^[1-9]$/)#all alternatealleles
	}
	if ($a1 ne "." && $a2 ne "."){$gs++} #count the number of genotypable samples

	#converting gt to sequence
	#option 1: true sequence, SNPs only by keep ref and alternate sequence as they are
	#option 2: pseudo-sequence, for SV only, treat ref allele as C and alternative allele as A
	#option 3: hybrid, for all: for SNPs following option 1, for all others, following option 2
	#option 4: pseudo-sequnece for all: 0/0 as C, 0/1 as G, 1/1 as A, "./." as N
	#option 5: pseudo-sequnece for SNPs: 0/0 as C, 0/1 as G, 1/1 as A, "./." as N
	#option 6: pseudo-sequnece for non-SNPs: 0/0 as C, 0/1 as G, 1/1 as A, "./." as N
	my ($s1,$s2,$s)=("N","N","N");
	if ($st=~/1/){#collecting SNPs only, use the reported real sequences
	    #if ($ref=~/^[AGCT]$/i && !$indel){#limit to SNPs
		if ($a1=~/0/){$s1=$ref} #ref allele, take ref seq/SNP
		elsif($a1=~/[1-9]/){$s1=$ALT[$a1-1]} #assign the corresponding alt sequence
		if ($a2=~/0/){$s2=$ref}
		elsif($a2=~/[1-9]/){$s2=$ALT[$a2-1]}
	    #}
	}elsif($st=~/2/){#pseudo-sequence, 2bp/loci: C for 0, A for 1, N for no call
	    #if ($indel){#limit to indels
		if ($a1=~/0/){$s1="C"}
		elsif($a1=~/1/){$s1="A"}
		elsif($a1=~/2/){$s1="T"}
		elsif($a1=~/[3-9]/){$s1="G"}
	    
		if ($a2=~/0/){$s2="C"}
		elsif($a2=~/1/){$s2="A"}
		elsif($a2=~/2/){$s2="T"} #multi-allele locus
		elsif($a2=~/[3-9]/){$s2="G"} #multi-allele locus
	    #}
	}elsif($st=~/3/){#hybrid, use real seq for SNPs & pseudo-seq for non-SNPs, C for 0, A for 1, N for no call
	    if ($ref=~/^[AGCT]$/i && !$indel){#SNPs
		if ($a1=~/0/){$s1=$ref}
		elsif($a1=~/[1-9]/){$s1=$ALT[$a1-1]}
		if ($a2=~/0/){$s2=$ref}
		elsif($a2=~/[1-9]/){$s2=$ALT[$a2-1]}
	    }elsif ($indel){#non-SNPs, pseudo-seq in 2bp/loci
		if ($a1=~/0/){$s1="C"}
		elsif($a1=~/1/){$s1="A"}
		elsif($a1=~/2/){$s1="T"}
		elsif($a1=~/[3-9]/){$s1="G"}

		if ($a2=~/0/){$s2="C"}
		elsif($a2=~/1/){$s2="A"}
		elsif($a2=~/2/){$s2="T"} #multi-allele locus
		elsif($a2=~/[3-9]/){$s2="G"} #multi-allele locus
	    }
	}elsif($st=~/[456]/){#pseudo-sequence, 1bp/loci: C for 0/0, A for 0/1, T for 1/1, and N for ./. 
	    if ($GT[$n]=~/^0\/0/){$s="C"} #homozygous ref allele
	    elsif($GT[$n]=~/^0\/[1-9]/){$s="A"} #heterozygous
	    elsif($GT[$n]=~/^1\/1/){$s="T"} #homozgyous alt allele
	    elsif($GT[$n]=~/^[1-9]\/[2-9]/){$s="G"} #anything beyond 1/1 for multi-allele locus
	    #missing minor other genotypes, e.g.: 2/2, 1/3, etc
	}
	my $seq=$s1.$s2; #diploid
	if ($st=~/[456]/){$seq=$s} #1 bp per locus/smp
	$SEQ{$sID}{$pos}=$seq;
	if ($opt{r} && $opt{r}==2){print "$sID\t$c:$p\t$ref/$alt\t$GT[$n]\t$a1/$a2\t$s1/$s2\t$seq\n"}
	if ($opt{r} && $opt{r}==3){print "$pos\r"}
	
    }#sample loop
    $VAR{$pos}->{A}=0;
    if ($ac && $ta){$VAR{$pos}->{A}=sprintf("%0.3f", ($ac/$ta)*100)}
    if ($opt{a} && $VAR{$pos}->{A} < $minAF && !$opt{P}){del($pos);delete($VAR{$pos});next}
    $VAR{$pos}->{G}=sprintf("%0.5f",($gs/$sn));
    if ($VAR{$pos}->{G} < $minGF && !$opt{P}){del($pos); delete($VAR{$pos}); next}
    #if ($VAR{$pos}){$VAR{$pos}->{P}=$p}
    if ($opt{r} && $opt{r}==10){if ($ln>=$opt{r}){print "$pos\n$SMP[0]\t$SEQ{$pos}{$SMP[0]}\n$SMP[9]\t$SEQ{$pos}{$SMP[9]}\n"; last}} #for testing purpose, process only a small number of loci
    if ($printloc && !$opt{P}){#generate location list
	$VAR{$pos}->{GT}="$ref/$alt";
	my $ST=1; #default variant type to report in the position file as 2-bp true sequence for a diploid locus
	if ($indel){
	    $ST=2; #2-bp pseudo-seq for a locus
	    if($st=~/[456]/){$ST=4} #single base pseudo sequence for a diploid locus
	}
	$VAR{$pos}->{ST}=$ST;
	if ($seen{$pos}){next}else{$seen{$pos}=1} #avoid redundant entries of loci
	if ($VAR{$pos}->{GT} && $VAR{$pos}->{P} && $VAR{$pos}->{G}){push @P, $pos} #generate a position array in the original order
	$gv++; #record the number of qualified variants
	if ($indel){$ind++}else{$snp++}
	#print "$pos\t$VAR{$pos}->{GT}\t$VAR{$pos}->{A}\t$gs/$sn:$VAR{$pos}->{G}\n";
    }#store the info for printing later
    if ($ln >= $batch && $ln % $batch == 0){#generate the sequence by batches to save memory by deleting elements in %SEQ and save the sequence to %FA
	print_data();
	if ($printloc && !$opt{P}){print_loc()}
    }
}#loop for variants/vcf data lines

close VCF;
undef %seen;
print_data(); #print the last batch of data
if ($printloc ne "0" && !$opt{P}){
    print_loc();
}
#print the sequences
foreach my $s (@SMP){
    print OUT ">$s\n$FA{$s}\n" if ($FA{$s})
}

if ($opt{r}){print STDERR "processed $opt{v} with a total of $sn samples and collected $gv variants (snp:$snp; indel:$ind) from a total of $ln for st:$st.\n"}
close OUT; close POS if ($printloc);
exit 0;

sub print_data{ #print sequences batch by batch
    foreach my $s (@SMP){
	foreach my $p (@P){
	    if (!$VAR{$p}){next}
	    #if ($opt{a} && $VAR{$p}->{A} < $minAF && !$opt{P}){next} #skip position with af below the cutoff
	    #if ($VAR{$p}->{G} < $minGF && !$opt{P}){next} #skip position with missing genotype
	    my $b="NN";
	    if ($st>3){$b="N"} #single base verse 2-bp per locus as the default value
	    if ($SEQ{$s}{$p}){$b=$SEQ{$s}{$p}}#update
	    $FA{$s} .=$b;
	    if ($opt{r} && $opt{r}==11){print "$s\t$p\t$b\n"}
	    delete($SEQ{$s}{$p});
	}
	#delete any remaining variants not passing the filter
	foreach my $p (keys %{$SEQ{$s}}){
	    delete($SEQ{$s}{$p});
	}
    }
}

sub del{#delete all %SEQ for a pos
    my $pos=shift;
    foreach my $s (@SMP){delete($SEQ{$pos}{$s})}
}

#print position list
sub print_loc{
    foreach my $P (@P){
	#if ($VAR{$p}){next}
	#if ($opt{a} && $VAR{$p}->{A} < $minAF && !$opt{P}){next} #skip position with af below the cutoff
	#if ($VAR{$p}->{G} < $minGF && !$opt{P}){next} #skip position with missing genotype
	if ($VAR{$P}->{P} && $VAR{$P}->{GT} && $VAR{$P}->{G}){print POS "$P\t$VAR{$P}->{P}\t$VAR{$P}->{GT}\t$VAR{$P}->{ST}\t$VAR{$P}->{A}\t$VAR{$P}->{G}\n"}
	delete($VAR{$P});
    }
    @P=();#reset @P
}
