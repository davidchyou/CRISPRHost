use Cwd 'getcwd';
use Cwd 'abs_path';

my $cd_path = abs_path($0); 
$cd_path = get_path($cd_path);

my $file_in = @ARGV[0];
my $file_out = @ARGV[1];
my $db = @ARGV[2];
my $fa = @ARGV[3];
my $b_cutoff = @ARGV[4];
my $pr_cutoff = @ARGV[5];
my $min_nocshftsbg = @ARGV[6];
my $non_crisprbank = @ARGV[7];
my $matched_spacer_fasta =  @ARGV[8];

open(INFILE, $file_in);
my @contents = <INFILE>;
close(INFILE);

my $header = "VIRAL_GENOME,SPACER_ID,ARRAY_SCORE,KNOWN_SPECIES,".
             "SPACER_ALIGN_START,SPACER_ALIGN_END,TARGET_ALIGN_START,TARGET_ALIGN_END,E_VALUE,BITSCORE,SPACER_LENGTH,".
             "N_SPACER_NT_MISSED_5P,N_SPACER_NT_MISSED_3P,BLAST_PIDENT,N_NT_IN_HYBRID,N_NT_MATCHED_IN_HYBRID,".
             "PROP_NT_MATCHED_IN_HYBRID,SEQ_SPACER_ALIGNED,SEQ_TARGET_ALIGNED,TARGET_FLANK_5P,TARGET_FLANK_3P,ORIGINAL_SPACER,SPACER_ORI_ALIGN,SPACER_ASSEMBLY,N_HOST_SPACER";
             
open(OUTFILE, ">>$file_out");
print OUTFILE "$header\n";

my %elements = ();
my %spacer_matches = ();

my @all_lines = ();

for (my $i = 0; $i < scalar(@contents); $i++) {
	my $line = @contents[$i];
	$line =~ s/\n//g;
	
	my $strout = processLine($line);
	my @toks = split(",", $strout);
	my $nn = scalar(@toks);
	my $target_id = @toks[0];
	my $cd_spacer_id = @toks[1];
	my $bitscore = @toks[9];
	my $prop_nt_matched = @toks[16];
	my $spacer_matched = @toks[21];
	my $gcf = @toks[$nn - 1];
	
	if ($bitscore < $b_cutoff) {
		next;
	}
	
	if ($prop_nt_matched < $pr_cutoff) {
		next;
	}
	
	if (exists $elements{$target_id}{$gcf}{$cd_spacer_id}) {
		$elements{$target_id}{$gcf}{$cd_spacer_id}++;
	} else {
		$elements{$target_id}{$gcf}{$cd_spacer_id} = 1;
	}
	
	my $viral_spacer_head = ">" . $target_id . "|" . $cd_spacer_id;
	if (! (exists $spacer_matches{$viral_spacer_head})) {
		$spacer_matches{$viral_spacer_head} = $spacer_matched;
	}
	
	push(@all_lines, $strout);
}

for (my $i = 0; $i < scalar(@all_lines); $i++) {
	my $line = @all_lines[$i];
	$line =~ s/\n//g;
	
	my @toks = split(",", $line);
	my $nn = scalar(@toks);
	my $target_id = @toks[0];
	my $cd_spacer_id = @toks[1];
	my $gcf = @toks[$nn - 1];
	
	my $count = 0;
	foreach my $key (keys %{$elements{$target_id}{$gcf}}) {
		$count++;
	}
	
	if ($count < $min_nocshftsbg) {
		next;
	}
	
	print OUTFILE "$line,$count\n";
}

close(OUTFILE);

open(FA, ">>$matched_spacer_fasta");
foreach my $key (keys %spacer_matches) {
	print FA $key . "\n";
	print FA $spacer_matches{$key} . "\n";
}
close(FA);

sub processLine {
	my ($line) = @_;
	
	my @toks = split(/[\t]/, $line);
	
	my $spacer_id = @toks[1];
	my @spacer_id_toks = split(/\|/, $spacer_id);
	
	my $cd_spacer_src = "";#@spacer_id_toks[1];
	my $cd_spacer_id = "";#@spacer_id_toks[2];
	my $array_score = "";#@spacer_id_toks[3];
	my $gcf = "";
	my $known_species = "";
	
	if (($non_crisprbank + 0.0) > 0) {
		$cd_spacer_src = $spacer_id;
		$cd_spacer_id = $spacer_id;
		$array_score = "";
		$gcf = $spacer_id;
		$known_species = $spacer_id;
		
	} else {
		$cd_spacer_src = @spacer_id_toks[1];
		$cd_spacer_id = @spacer_id_toks[2];
		$array_score = @spacer_id_toks[3];
		for (my $i = 5; $i < scalar(@spacer_id_toks); $i++) {
			my $sp_name = @spacer_id_toks[$i];
			$sp_name =~ s/\_/ /g;
		
			if ($i == scalar(@spacer_id_toks) - 1) {
				$known_species .= "$sp_name";
			} else {
				$known_species .= "$sp_name;";
			}
		}
	
		my @src_toks = split("#", $cd_spacer_src);
		$gcf = @src_toks[1];
	
		if (length($known_species) == 0) {
			my $organism_details = @spacer_id_toks[1];
			$known_species = organism_short_name($organism_details);
		}
	}
	
	my $target_id = @toks[0];
	my $blast_pident = @toks[2];
	my $query_start = @toks[8];
	my $query_end = @toks[9];
	
	if ($query_start > $query_end) {
		my $temp = $query_start;
		$query_start = $query_end;
		$query_end = $temp;
	}
	
	my $target_start = @toks[6];
	my $target_end = @toks[7];
	
	if ($target_start > $target_end) {
		my $temp = $target_start;
		$target_start = $target_end;
		$target_end = $temp;
	}
	
	my $e_value = @toks[10];
	my $bitscore = @toks[11];
	
	my $query_seq_aligned = @toks[12];
	my $target_seq_aligned = @toks[13];
	my $n_match = @toks[14];
	my $target_length = @toks[15];
	my $query_length = @toks[16];
	my $spacer_ori_align = @toks[17];
	
	my $n_spacer_nt_missed_5p =  $query_start - 1;
	my $n_spacer_nt_missed_3p =  abs($target_length - $query_end);
	
	my $n_nt_align = num_bases_crt($query_seq_aligned, $target_length);
	my $n_nt_matched = num_bases_matched_crt($n_match);
	my $prop_nt_matched = $n_nt_matched / $n_nt_align;
	
	my $flank5p = extract_seq($db, $target_id, $target_start - 20, $target_start - 1, $query_length);
	my $flank3p = extract_seq($db, $target_id, $target_end + 1, $target_end + 20, $query_length);
	my $orig_spacer = extract_spacer($fa, $cd_spacer_id);
	
	my $out = "$target_id,$cd_spacer_id,$array_score,$known_species,".
	           "$query_start,$query_end,$target_start,$target_end,$e_value,$bitscore,$target_length,".
	           "$n_spacer_nt_missed_5p,$n_spacer_nt_missed_3p,$blast_pident,$n_nt_align,$n_nt_matched,".
	           "$prop_nt_matched,$target_seq_aligned,$query_seq_aligned,$flank5p,$flank3p,$orig_spacer,$spacer_ori_align,$gcf";
    
    return $out;
}

sub ngaps {
	my ($str) = @_;
	my $count = 0;
	my @arr = split('', $str);
	
	for (my $i = 0; $i < scalar(@arr); $i++) {
		my $char = @arr[$i];
		if ($char eq "-") {
			$count++;
		}
	}
	return $count;
}

sub num_bases_crt {
	my ($sseq_aln, $spacer_len) = @_;
	
	my $ngaps = ngaps($sseq_aln);
	my $nbases = length($sseq_aln) - ngaps;
	my $total = $spacer_len + $nbases;
	return $total;
}

sub num_bases_matched_crt {
	my ($n_matches) = @_;
	my $total = 2 * $n_matches;
	return $total;
}

sub extract_seq {
	my ($db_path, $contig, $start, $end, $length) = @_;
	
	my $x = $start;
	my $y = $end;
	
	if ($x < 1) {
		$x = 1;
	}
	
	if ($y < 1) {
		$y = 1;
	}
	
	if ($x > $length) {
		$x = $length;
	}
	
	if ($y > $length) {
		$y = $length;
	}
	
	if ($x > $y) {
		my $z = $x;
		$x = $y;
		$y = $z;
	}
	
	my $range = $x . "-" . $y;
	
	my $cmd_bdb = "$cd_path/blastdbcmd -db $db_path -entry \'$contig\' -range $range";
	
	my @rs = `$cmd_bdb`;
	
	if (scalar(@rs) < 2) {
		return "";
	}
	
	if (index(@rs[0], "Error") >= 0) {
		return "";
	}	
	
	if (index(@rs[0], ">") < 0) {
		return "";
	}
	
	my $seq = "";
	for (my $j = 1; $j < scalar(@rs); $j++) {
		my $seqstr = @rs[$j];
		$seqstr =~ s/\n//g;
		$seq .= $seqstr;
	}
	
	return $seq;
}

sub extract_spacer {
	my ($fa_path, $spacer_id) = @_;
	
	my $cmd = "grep $spacer_id $fa_path -A 1";
	
	my @rs = `$cmd`;
	
	if (scalar(@rs) < 2) {
		return "";
	}
	
	my $seq = @rs[1];
	$seq =~ s/\n//g;
	
	return $seq;
}

sub organism_short_name {
	my ($genome_details) = @_;
	
	my @gtoks_1 = split("#", $genome_details);
	my $gcf = "NA";
	if (scalar(@gtoks_1) > 1) {
		$gcf = @gtoks_1[1];
	}
	
	my @gtoks_2 = split(/[_]/, @gtoks_1[0]);
	my $organism = "Unknown";
	
	if (scalar(@gtoks_2) >= 1) {
		$organism = @gtoks_2[0];
	} 
	
	if ($organism eq "Candidatus") {
		my $str = "";
		
		if (scalar(@gtoks_2) == 2) {
			$str = @gtoks_2[1];
		}
		
		if (scalar(@gtoks_2) > 2) {
			$str = @gtoks_2[1] . "_" . @gtoks_2[2];
		}
		$organism .= "_" . "$str";
	} else {
		if (scalar(@gtoks_2) >= 2) {
			$organism .= "_" . "@gtoks_2[1]";
		}
	}
	
	if ((scalar(@gtoks_2) > 2) && ((@gtoks_2[1] eq "sp.") || (@gtoks_2[1] eq "Sp."))) {
		$organism .= "_" . "@gtoks_2[2]";
	}
	return $organism;
}

sub get_path() {
	my $dir=shift(@_);
	my @arr_p1=split('\/',$cd_path);
	pop(@arr_p1);
	$dir=join("\/",@arr_p1);
		
	return $dir;
}