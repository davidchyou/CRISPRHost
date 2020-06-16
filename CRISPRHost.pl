use Cwd 'getcwd';
use Cwd 'abs_path';
use Bio::SeqIO;

my $cd_path = abs_path($0); 
$cd_path = get_path($cd_path);

my $path = "";
my $outpath = "";
my $kingdom = "AB";
my $dbpath = "$cd_path/DB/CRISPRBankSpacers_4_95_2555_100_all_refseq_nr.fa";
my $user_dbpath = "";
my $fasta_path = "";
my $e_cutoff = 1e-8;
my $max_target_seqs = 500;
my $dbsize = 1e4;
my $b_cutoff = 50;
my $pr_cutoff = 0.7;
my $min_nocshftsbg = 1;
my $non_crisprbank = 0;
my $mask_arrays = 0;

my $ind = 0;
foreach(@ARGV) {
	if (@ARGV[$ind] eq '-in') {
		$path = @ARGV[$ind + 1];
		#if (! (-e $path)) {
		#	die "cannot open file: " . $path;
		#}
	}
	
	if (@ARGV[$ind] eq '-out') {
		$outpath = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-kingdom') {
		$kingdom = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-db') {
		$dbpath = @ARGV[$ind + 1];
		$user_dbpath = @ARGV[$ind + 1];
		#if (! (-e $dbpath)) {
		#	die "cannot find DB: " . $dbpath;
		#}
	}
	
	if (@ARGV[$ind] eq '-dbsize') {
		$dbsize = int(@ARGV[$ind + 1]);
	}
	
	if (@ARGV[$ind] eq '-e_cutoff') {
		$e_cutoff = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-b_cutoff') {
		$b_cutoff = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-bp_cutoff') {
		$pr_cutoff = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-min_nocsh') {
		$min_nocshftsbg = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-fasta_path') {
		$fasta_path = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-non_crisprbank') {
		$non_crisprbank = 1;
	}
	
	if (@ARGV[$ind] eq '-max_target_seqs') {
		$max_target_seqs = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-mask_arrays') {
		$mask_arrays = 1;
	}
	
	if (@ARGV[$ind] eq '-h') {
		print "\n";
		system("cat help.txt");
		print "\n";
		exit;
	}
	
	if (@ARGV[$ind] eq '-help') {
		print "\n";
		system("cat help.txt");
		print "\n";
		exit;
	}
	
	$ind++;
}

if (! (-d "$cd_path/DB")) {
	system("unzip $cd_path/DB.zip -d $cd_path");
}

if ($kingdom eq "A") {
	$dbpath = "$cd_path/DB/CRISPRBankSpacers_4_95_2555_100_archaea_refseq_nr.fa";
} elsif ($kingdom eq "B") {
	$dbpath = "$cd_path/DB/CRISPRBankSpacers_4_95_2555_100_bacteria_refseq_nr.fa";
} else {
	$dbpath = "$cd_path/DB/CRISPRBankSpacers_4_95_2555_100_all_refseq_nr.fa";
}

if ((length($user_dbpath) > 0) && (-e $user_dbpath)) {
	$dbpath = $user_dbpath;
}

if (-d $outpath) {
	system("rm -rf $outpath");
}
mkdir($outpath);

if ((length($fasta_path)) > 0 && (-e $fasta_path)) {
	system("cp $fasta_path $outpath/input_spacers.fa");
	my $cmd_fadb = "makeblastdb -in $outpath/input_spacers.fa -out $outpath/input_spacers.fa -dbtype nucl -parse_seqids";
	$flag = system($cmd_fadb);
	$dbpath = "$outpath/input_spacers.fa";
}

my $path_2 = "$outpath/input_genomes.fa";
my $path_3 = "$outpath/input_genomes_masked.fa";
my $path_mask = "$outpath/c_arrays.gff";

my $in = Bio::SeqIO->new('-file' => "$path",  '-format' => 'Fasta');
open(FNA, ">>$path_2");
while (my $seq = $in->next_seq()) {
	my $seqid = $seq->id();
	my $seq_str = $seq->seq();
	
	print FNA ">$seqid\n";
	print FNA "$seq_str\n";
}
close(FNA);

my $minced_path=`which minced >&1 2>&1`; chomp $minced_path; $minced_path=~s/\r//g;
if(not defined $minced_path or $minced_path eq "" or $minced_path=~/:/) {
	if ($mask_arrays > 0) {
		print "CRISPR-array predictor \'minced\' not found, skip the array-masking step.\n";
	}
	$mask_arrays = 0;
}

if ($mask_arrays > 0) {
	my $flag = system("minced -gff $path_2 $path_mask >&1 2>&1");
	$flag = system("bedtools maskfasta -fi $path_2 -bed $path_mask -fo $path_3 >&1 2>&1");
}
my $flag = system("makeblastdb -in $path_2 -out $path_2 -dbtype nucl -parse_seqids");

my $qfile = $path_2;
if (($mask_arrays > 0) && (-e $path_3)) {
	$qfile = $path_3;
}

$flag = system("blastn -query $qfile -db $dbpath -word_size 7 -out $outpath/out.txt -gapopen 10 -max_target_seqs $max_target_seqs " . 
		       "-gapextend 2 -penalty -1 -reward 1 -task blastn-short -lcase_masking -outfmt " .
		       "\'6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq nident slen qlen sstrand\' " .
		       "-dbsize $dbsize -evalue $e_cutoff -num_threads 24 >/dev/null 2>&1");
		          
$flag = system("perl $cd_path/ProcessCRISPRHostOutput.pl $outpath/out.txt $outpath/out.txt.csv $path_2 $dbpath $b_cutoff $pr_cutoff $min_nocshftsbg $non_crisprbank $outpath/spacer_matched.fna");

open(RAWDAT, "$outpath/out.txt.csv");
my @contents = <RAWDAT>;
close(RAWDAT);

my $head = @contents[0];
$head =~ s/\n//g;

open(TEMP, ">>$outpath/out.txt.csv.temp");

for (my $j = 1; $j < scalar(@contents); $j++) {
	my $line = @contents[$j];
	$line =~ s/\n//g;
	print TEMP $line . "\n";

}
close(TEMP);

open(SDAT, ">>$outpath/out.txt.csv.temp2");
print SDAT $head . "\n";
 $flag = system("cat $outpath/out.txt.csv.temp" . 
	            " \| sort -t $\',\' -k17nr,17 -k10nr,10 >> " . 
	            "$outpath/out.txt.csv.temp2");
unlink("$outpath/out.txt.csv.temp");
unlink("$outpath/out.txt.csv");
$flag = system("mv $outpath/out.txt.csv.temp2 $outpath/out.txt.csv");

my $rcmd = "source(\"$cd_path/CSVOutput.r\")\; " .
           "flag <- summary_csv\(\"$outpath/out.txt.csv\"\, \"$outpath/summ.csv\"\)\;" .
           "flag <- reorder_csv\(\"$outpath/out.txt.csv\"\, \"$outpath/full.csv\"\)\;"; 
system("Rscript -e '$rcmd' 2>&1");

open(SUMM, "$outpath/summ.csv");
@contents = <SUMM>;
close(SUMM);

open(SUMMOUT, ">>$outpath/summary.csv");
for (my $k = 0; $k < scalar(@contents); $k++) {
	my $line = @contents[$k];
	$line =~ s/\n//g;
	$line =~ s/\"//g;
	print SUMMOUT $line . "\n";

}
close(SUMMOUT);
unlink("$outpath/summ.csv");

open(FULL, "$outpath/full.csv");
@contents = <FULL>;
close(FULL);

open(FULLOUT, ">>$outpath/full_results.csv");
for (my $k = 0; $k < scalar(@contents); $k++) {
	my $line = @contents[$k];
	$line =~ s/\n//g;
	$line =~ s/\"//g;
	print FULLOUT $line . "\n";

}
close(FULLOUT);
unlink("$outpath/full.csv");

sub get_path() {
	my $dir=shift(@_);
	my @arr_p1=split('\/',$cd_path);
	pop(@arr_p1);
	$dir=join("\/",@arr_p1);
		
	return $dir;
}