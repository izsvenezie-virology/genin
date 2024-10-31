#! /usr/bin/perl

#USAGE: perl genotype_predictor.pl --fasta FILE --output DIRECTORY --Ncpu NUMBER --config FILE --NNparameters FILE --version NUMBER
#MUST-HAVE PARAMETERS:
#--fasta FILE: fasta file with unknown sequence to whom predict genotype;
#	N.B.: different segments belonging to the same virus must share the same exactly string for sample name; segments must be indicated with one of the following; PA, PB1, PB2, MP, NA, NS, NP; sample name and segment must be separated by '$config{'SEPARATOR'}' symbol;
#--out OUTPUT: directory for output files;
#OPTIONAL PARAMETERS:
#--Ncpu NUMBER: maximum number of cpus to be used for parallelization;
#DEBUGGING PARAMETERS:
#--config FILE: configuration file to be used; if not present, most recent version will be used;
#--NNparameters FILE: neural net parameter file to be used; if not present, most recent version will be used;
#--version NUMBER: version number to be used for the prediction;
#       N.B.2: parameters "--config", "--NNparameters" and "--version", if present, must be always together;

#cartella contenente i file di configurazione e dei parametri
my $path = $0;
$path =~ s/\/genin\.pl$//;
my $dir_data = $path.'/data';
my $dir_src = $path.'/src';
my (@content_list);
opendir (my $opened_dir, $dir_data) || print STDERR "WARNING: cannot open directory \"$dir_data\" containing similarity tables, seq2ver files and saved neural nets.\n";
while (my $line = readdir $opened_dir) {push (@content_list, $line);}
closedir $opened_dir;
opendir (my $opened_dir, $dir_src) || print STDERR "WARNING: cannot open directory \"$dir_src\" containing configuration file and neural neet parameters.\n";
while (my $line = readdir $opened_dir) {push (@content_list, $line);}
closedir $opened_dir;

#raccolgo le variabili indicate
my $fasta_file = '';
my $output_directory = '';
my $N_cpu = 1;
my $N_cpu_flag = 0;
my $config_file = '';
my $neural_net_parameters_file = '';
my $version = '';
my $index = 0;
while ($index <= $#ARGV)
{
	if ($ARGV[$index] eq '--fasta')
	{
		$index++;
		$fasta_file = $ARGV[$index];
		print STDOUT "--fasta\t".$fasta_file."\n";
	}
	elsif ($ARGV[$index] eq '--output')
	{
		$index++;
		$output_directory = $ARGV[$index];
		print STDOUT "--output\t".$output_directory."\n";
	}
	elsif ($ARGV[$index] eq '--Ncpu')
	{
		$index++;
		$N_cpu = $ARGV[$index];
		print STDOUT "--Ncpu\t".$N_cpu."\n";
		$N_cpu_flag = 1;
	}
	elsif ($ARGV[$index] eq '--config')
	{
		$index++;
		$config_file = $ARGV[$index];
		print STDOUT "--config\t".$config_file."\n";
	}
	elsif ($ARGV[$index] eq '--NNparameters')
	{
		$index++;
		$neural_net_parameters_file = $ARGV[$index];
		print STDOUT "--NNparameters\t".$neural_net_parameters_file."\n";
	}
	elsif ($ARGV[$index] eq '--version')
	{
		$index++;
		$version = $ARGV[$index];
		print STDOUT "--version\t".$version."\n";
	}
	else {print STDERR "WARNING: unknown parameter \"$ARGV[$index]\".\n";}
	$index++;
}

#se non indicate, trovo le versioni più recenti dei file di configurazione e coi parametri e del dataset
if ($config_file eq '')
{
	$config_file = &find_last_version ($dir_src, \@content_list, 'config_v', '.csv');
	print STDOUT "--config (automatically found)\t".$config_file."\n";
}
if ($neural_net_parameters_file eq '')
{
	$neural_net_parameters_file = &find_last_version ($dir_src, \@content_list, 'neural_net_parameters_v', '.csv');
	print STDOUT "--NNparameters (automatically found)\t".$neural_net_parameters_file."\n";
}
my $dataset_directory = '';
if ($version eq '')
{
	$dataset_directory = &find_last_version ($dir_data, \@content_list, 'version', '');
	$dataset_directory =~ /^$dir_data\/version([0-9])+$/;
	$version = $1;
	print STDOUT "--version (automatically found)\t".$version."\n";
	print STDOUT "\tdataset used (automatically found)\t".$dataset_directory."\n";
}
else {$dataset_directory = $dir_data.'/version'.$version;}

#se non settato, indico che ho usato 1 cpu
if ($N_cpu_flag == 0) {print STDOUT "used 1 cpu\n";}

#se necessario, creo la directory di output
&crea_dir ($output_directory);

#memorizzo i parametri della rete neurale da utilizzare
my (%param);
open (I, "<$neural_net_parameters_file") or die "ERROR: parameter \"--NNparameters\" is wrong, please check it, script stopped..";
my $header = <I>;
chomp $header;
my @header = split (/;/, $header);
while (my $line = <I>)
{
	chomp $line;
	$line =~ s/\/versionX\//\/version$version\//g;
	my @line = split (/;/, $line);
	for (my $i = 1; $i <= $#line; $i++) {$param{$header[$i]}{$line[0]} = $line[$i];
}
}
close I;
if (! -e $path.'/'.$param{'seq2ver_file'}{'PA'}) {die "ERROR: something is wrong with parameter \"--NNparameters\", please check it, script stopped.";}

#memorizzo le variabili generali
my (%config);
open (I, "<$config_file") or die "ERROR: parameter \"--config\" is wrong, please check it, script stopped..";
while (my $line = <I>)
{
	chomp $line;
	my @line = split ("\t", $line);
	$config{$line[0]} = $line[1];
}

#memorizzo le sequenze di ogni campione del dataset
my (%dataset_sequence, %N_dataset_sequences);
foreach $seg (sort(keys(%{$param{'seq2ver_file'}})))
{
	$N_dataset_sequences{$seg} = 0;
	my $file = $path.'/'.$param{'seq2ver_file'}{$seg};
	open (I, "<$file") or print STDERR "WARNING: file \"$file\" does not exists.\n";
	while (my $line = <I>)
	{
		chomp $line;
		my @line = split ("\t", $line);
		$dataset_sequence{$seg}{$line[0]} = $line[1];
		$N_dataset_sequences{$seg}++;
	}
	close I;
}

#memorizzo le sequenze dei campioni dei quale predire il genotipo
my @segments = split (';', $config{'SEGMENTS'});
my (%sample_sequence, %sample2modified_name, %N_sample_sequences);
my ($flag, $temp_sample_name, $temp_seg);
open (I, "<$fasta_file") or die "ERROR: file \"$file\" does not exists, script stopped.";
while (my $line = <I>)
{
	chomp $line;
	if ($line eq '') {next;}
	if ($line =~ /^>/)
	{
		$flag = 0;
		for (my $i = 0; $i <= $#segments; $i++)
		{
			if ($line =~ /$config{'SEPARATOR'}$segments[$i]$/)
			{
				$flag = 1;
				$temp_seg = $segments[$i];
				last;
			}
		}
		if ($flag == 1)
		{
			$line =~ /^>(.+)$config{'SEPARATOR'}$temp_seg$/;
			$temp_sample_name = $1;
			if (exists $sample_sequence{$temp_sample_name}{$temp_seg})
			{
				$flag = 0;
				print STDERR "WARNING: $line is present more than 1 time in fasta file \"$fasta_file\", it will be ignored.\n";
			}
			else
			{
				$sample2modified_name{$temp_sample_name} = &not_allowed_char2underscore ($temp_sample_name);
				$sample_sequence{$temp_sample_name}{$temp_seg} = "";
				if (exists $N_sample_sequences{$temp_seg}) {$N_sample_sequences{$temp_seg}++;}
				else {$N_sample_sequences{$temp_seg} = 1;}
			}
		}
		else {print STDERR "WARNING: $line in \"$fasta_file\" file is not properly formatted, it will be ignored.\n";}
	}
	else
	{
		if ($flag == 1)
		{
			$sample_sequence{$temp_sample_name}{$temp_seg}.= &maiuscolo ($line);
			$sample_sequence{$temp_sample_name}{$temp_seg} =~ s/N//g;
		}
	}
}

#calcolo le similarità di sequenza tra quelle in input e quelle nel dataset
my $N_total_jobs = 0;
foreach my $key (sort(keys(%N_sample_sequences))) {$N_total_jobs+= $N_sample_sequences{$key} * $N_dataset_sequences{$key};}
my $scratch_dir = $output_directory.'/computation';
&crea_dir ($scratch_dir);
my $chunk_size = 1 + sprintf("%.0f", $N_total_jobs / $N_cpu);
my $chunk_index = 1;
my $command_file = $scratch_dir.'/chunk_N'.$chunk_index.'.sh';
open (C, ">$command_file");
my $N_seq = 1;
foreach my $sample (sort(keys(%sample_sequence)))
{
	foreach my $seg (sort(keys(%{$sample_sequence{$sample}})))
	{
		foreach my $dataset_sample (sort(keys(%{$dataset_sequence{$seg}})))
		{
			print C 'echo -ne "'.$seg.'\t'.$sample.'\t'.$dataset_sample.'\t"'."\n";
			print C 'echo -e ">'.$sample.'\n'.$sample_sequence{$sample}{$seg}.'\n>'.$dataset_sample.'\n'.$dataset_sequence{$seg}{$dataset_sample}.'" | mafft --thread 1 --auto /dev/stdin 2>/dev/null | perl '.$path.'/'.$config{'SIMILARITY_SCRIPT'}.' 2>/dev/null'."\n";
			$N_seq++;
			if ($N_seq > $chunk_size)
			{
				print C 'echo "FINISHED"';
				close C;
				my $comando = 'nohup bash '.$command_file.' >'.$scratch_dir.'/chunk_N'.$chunk_index.'.txt &';
				system ($comando);
				$chunk_index++;
				$command_file = $scratch_dir.'/chunk_N'.$chunk_index.'.sh';
				open (C, ">$command_file");
				$N_seq = 1;
			}
		}
	}
}
if ($N_seq == 1)
{
	close C;
	$chunk_index--;
}
else
{
	print C 'echo "FINISHED"';
	close C;
	my $comando = 'nohup bash '.$command_file.' >'.$scratch_dir.'/chunk_N'.$chunk_index.'.txt &';
	system ($comando);
}

#controllo che tutti i valori di similarità necessari siano stati calcolati
while (1)
{
	my $count = 0;
	for (my $i = 1; $i <= $chunk_index; $i++)
	{
		my $file = $scratch_dir.'/chunk_N'.$i.'.txt';
		open (I, "<$file");
		while (my $line = <I>)
		{
			chomp $line;
			if ($line eq 'FINISHED') {$count++;}
		}
		close I;
	}
	if ($count == $chunk_index) {last;}
}

#raccolgo i valori di similarità calcolati
my (%simil_values);
for (my $i = 1; $i <= $chunk_index; $i++)
{
	my $file = $scratch_dir.'/chunk_N'.$i.'.txt';
	open (I, "<$file") or die "ERROR: cannot open file \"$file\".";
	while (my $line = <I>)
	{
		chomp $line;
		my @line = split ("\t", $line);
		if ($line ne 'FINISHED') {$simil_values{$line[1]}{$line[0]}{$line[2]} = $line[3];}
	}
	close I;
}

#definisco il sottotipo di NA
my @subtypes = split (';', $config{'NA_subtype'});
my (%sample2NAsub, %sample2NAval);
foreach my $sample (sort(keys(%sample_sequence)))
{
	$sample2NAsub{$sample} = 'unknown';
	$sample2NAval{$sample} = 0;
	if (exists $sample_sequence{$sample}{'NA'})
	{
		my $dataset_sample_with_max_value;
		foreach my $dataset_sample (sort(keys(%{$simil_values{$sample}{'NA'}})))
		{
			if ($simil_values{$sample}{'NA'}{$dataset_sample} > $sample2NAval{$sample})
			{
				$sample2NAval{$sample} = $simil_values{$sample}{'NA'}{$dataset_sample};
				$dataset_sample_with_max_value = $dataset_sample;
			}
		}
		if ($sample2NAval{$sample} >= $param{'ProbabilityThreshold'}{'NA'})
		{
			for (my $i = 0; $i <= $#subtypes; $i++)
			{
				if (exists $dataset_sequence{$subtypes[$i]}{$dataset_sample_with_max_value})
				{
					$sample2NAsub{$sample} = $subtypes[$i];
					last;
				}
			}
		}
		else {print STDERR "WARNING: it was not possible to find proper NA subtype for sample $sample, its genotype will be predicted ignoring NA subtype.\n";}
	}
}

#eseguo la predizione per i singoli segmenti di tutti i campioni salvando e stampando le probabilità
my (%P, %ver);
foreach my $sample (sort(keys(%sample_sequence)))
{
	&crea_dir ($output_directory.'/'.$sample2modified_name{$sample});
	foreach my $seg (sort(keys(%{$sample_sequence{$sample}})))
	{
		my $outdir_sample_seg_prediction = $output_directory.'/'.$sample2modified_name{$sample}.'/'.$seg.'_prediction';
		&crea_dir ($outdir_sample_seg_prediction);
		if ( ($seg eq 'NA') && ($sample2NAsub{$sample} ne 'unknown') && ($param{'OnlyVer'}{$sample2NAsub{$sample}} > 0) )
		{
			$ver{$sample}{'NA'} = $param{'OnlyVer'}{$sample2NAsub{$sample}};
			$P{$sample}{'NA'} = 1;
		}
		elsif ( ($seg ne 'NA') || ( ($seg eq 'NA') && ($sample2NAsub{$sample} ne 'unknown') && ($param{'OnlyVer'}{$sample2NAsub{$sample}} == 0) ) )
		{
			my $real_seg = $seg;
			if ($seg eq 'NA') {$real_seg = $sample2NAsub{$sample};}
			my $input_file = $outdir_sample_seg_prediction.'/input_data_nnet.txt';
			open (O, ">$input_file");
			print O $sample;
			foreach my $dataset_sample (sort(keys(%{$dataset_sequence{$real_seg}}))) {print O "\t".$simil_values{$sample}{$seg}{$dataset_sample};}
			print O "\n";
			close O;
			my $output_file = $outdir_sample_seg_prediction.'/prediction_nnet.txt';
			my $model = $path.'/'.$param{'NeuralNet'}{$real_seg}.'_model.Rdata';
			my $comando = 'cat '.$input_file.' | Rscript '.$path.'/'.$config{'NEURAL_NET_PREDICTOR'}.' '.$model.' '.$output_file.' >'.$output_file.'.log 2>'.$output_file.'.err';
			system ($comando);
			$output_file = &correggi_formato_predizione ($output_file, $path.'/'.$config{'NEURAL_NET_CORRECTOR'}, $model);
			my ($chosen_version, $best_P) = &take_best_probability ($output_file);
			$ver{$sample}{$seg} = $chosen_version;
			$P{$sample}{$seg} = $best_P;
		}
	}
}

#memorizzo la composizione dei genotipi
my (%composition);
$config{'GENOTYPE_COMPOSITION'} =~ s/\/versionX\//\/version$version\//g;
my $file = $path.'/'.$config{'GENOTYPE_COMPOSITION'};
open (I, "<$file");
my $header = <I>;
chomp $header;
my @header = split ("\t", $header);
while (my $line = <I>)
{
	chomp $line;
	my @line = split ("\t", $line);
	for (my $i = 1; $i <= $#line; $i++) {push (@{$composition{$header[$i]}{$line[$i]}}, $line[0]);}
}
close I;

#calcolo il genotipo di appartenenza e stampo i risultati
my $sep = ";";
my $results_file = $output_directory.'/result.txt';
open (O, ">$results_file");
print O "Sample";
for (my $i = 0; $i <= $#header; $i++) {print O $sep.$header[$i].$sep.$header[$i].' Probability';}
print O $sep."Note\n";
foreach my $sample (sort(keys(%sample_sequence)))
{
	print O $sample.$sep;
	my (%genotype_count, %reliability, %lacking);
	my $cumul_P = 1;
	my $flag_rel = 0, $flag_lack = 7;
	for (my $i = 0; $i <= $#segments; $i++)
	{
		my $seg = $segments[$i]; #variabile helper
		$reliability{$seg} = 1;
		$lacking{$seg} = 1;
		if (exists $ver{$sample}{$seg})
		{
			my $version = $ver{$sample}{$seg}; #variabile helper
			my $P = $P{$sample}{$seg}; #variabile helper
			for (my $j = 0; $j <= $#{$composition{$seg}{$version}}; $j++)
			{
				my $gen = ${$composition{$seg}{$version}}[$j];
				if (exists $genotype_count{$gen}) {$genotype_count{$gen}++;}
				else {$genotype_count{$gen} = 1;}
			}
			$lacking{$seg} = 0;
			$flag_lack--;
			if ($P < $param{'ProbabilityThreshold'}{$seg})
			{
				$reliability{$seg} = 0;
				$flag_rel++;
			}
			$cumul_P = $cumul_P * $P;
		}
	}
	my $best_score = 0;
	foreach my $genotype (sort(keys(%genotype_count)))
	{
		if ($genotype_count{$genotype} > $best_score) {$best_score = $genotype_count{$genotype};}
	}
	my (@chosen_genotpe);
	foreach my $genotype (sort(keys(%genotype_count)))
	{
		if ($genotype_count{$genotype} == $best_score) {push (@chosen_genotpe, $genotype);}
	}	
	print O $chosen_genotpe[0];
	for (my $i = 1; $i <= $#chosen_genotpe; $i++) {print O ",".$chosen_genotpe[$i];}
	if ( ($flag_rel > 0) || ($best_score < 7) ) {print O "*";}
	print O $sep.$cumul_P;
	for (my $i = 1; $i <= $#header; $i++)
	{		
		if (exists $ver{$sample}{$header[$i]})
		{
			my $piece = '';
			if ($reliability{$header[$i]} == 0) {$piece = '*';}
			print O $sep.$ver{$sample}{$header[$i]}.$piece.$sep.$P{$sample}{$header[$i]};
		}
		else {print O $sep.'NA'.$sep.'NA';}
	}
	if ($flag_rel > 0)
	{
		print O $sep."it was not possible to make a reliable prediction for segment";
		if ($flag_rel > 1) {print O "s";}
		my $piece = " ";
		foreach my $seg_not_reliable (sort(keys(%reliability)))
		{
			if ($reliability{$seg_not_reliable} == 0) {$piece.= $seg_not_reliable.",";}
		}
		chop $piece;
		print O $piece;
	}
	if ($flag_lack > 0)
	{
		if ($flag_rel > 0) {print O " - ";}
		else {print O $sep;}
		print O "segment";
		if ($flag_lack > 1) {print O "s";}
		my $piece = " ";
		foreach my $seg_lacking (sort(keys(%lacking)))
		{
			if ($lacking{$seg_lacking} == 1) {$piece.= $seg_lacking.",";}
		}
		chop $piece;
		print O $piece." ";
		if ($flag_lack > 1) {print O "were";}
		else {print O "was";}
		print O " not present in input fasta file";
	}
	if ( ($flag_rel == 0) && ($flag_lack == 0) && ($best_score < 7) ) {print O "potential novel genotype due to original segment version composition";}
	print O "\n";
}
close O;

#trova la versione più recente dell'elemento indicato
sub find_last_version
{
	#variabili iniziali
	my ($dir, $list, $pattern_left, $pattern_right);
	($dir, $list, $pattern_left, $pattern_right) = @_;

	#scorro la lista e seleziono solo gli elementi di interesse
	my (%temp);
	for (my $i = 0; $i <= $#{$list}; $i++)
	{
		if ($list->[$i] =~ /^$pattern_left[0-9]+$pattern_right$/)
		{
			$list->[$i] =~ /^$pattern_left([0-9]+)$pattern_right$/;
			my $version = $1;
			$temp{$version} = $version;
		}
	}

	#trovo la versione più recente
	my $file = $dir.'/'.$pattern_left;
	foreach my $version (sort { $temp{$b} <=> $temp{$a} } keys %temp)
	{
		$file.= $version;
		last;
	}
	$file.= $pattern_right;
	return $file;
}

#trasforma tutti i caratteri non-alfanumerici ad eccezione dell'underscore di una stringa in un underscore
sub not_allowed_char2underscore
{
	#variabili iniziali
	my ($stringa);
	($stringa) = @_;

	#converte tutti i caratteri non-alfanumerici in underscore
	my $stringa_m = "";
	my @stringa = split (//, $stringa);
	for (my $i = 0; $i <= $#stringa; $i++)
	{
		if ( ($stringa[$i] =~ /[a-z]/) || ($stringa[$i] =~ /[A-Z]/) || ($stringa[$i] =~ /[0-9]/) || ($stringa[$i] =~ /_/) ){$stringa_m.= $stringa[$i];}
		else {$stringa_m.= '_';}
	}
	return ($stringa_m);
}

#accerta che tutti i nucleotidi siano scritti in maiuscolo
sub maiuscolo
{
	#variabili iniziali
	my ($seq);
	($seq) = @_;

	#converte gli eventuali nucleotidi scritti in minuscolo in maiuscolo
	my $seq_m = "";
	my @seq = split (//, $seq);
	for (my $i = 0; $i <= $#seq; $i++)
	{
		if ( ($seq[$i] eq 'A') || ($seq[$i] eq 'a') ) {$seq_m.= 'A'; next;}
		if ( ($seq[$i] eq 'C') || ($seq[$i] eq 'c') ) {$seq_m.= 'C'; next;}
		if ( ($seq[$i] eq 'G') || ($seq[$i] eq 'g') ) {$seq_m.= 'G'; next;}
		if ( ($seq[$i] eq 'T') || ($seq[$i] eq 't') ) {$seq_m.= 'T'; next;}
		if ( ($seq[$i] eq 'N') || ($seq[$i] eq 'n') ) {$seq_m.= 'N'; next;}
		if ( ($seq[$i] eq 'R') || ($seq[$i] eq 'r') ) {$seq_m.= 'R'; next;}
		if ( ($seq[$i] eq 'Y') || ($seq[$i] eq 'y') ) {$seq_m.= 'Y'; next;}
		if ( ($seq[$i] eq 'S') || ($seq[$i] eq 's') ) {$seq_m.= 'S'; next;}
		if ( ($seq[$i] eq 'W') || ($seq[$i] eq 'w') ) {$seq_m.= 'W'; next;}
		if ( ($seq[$i] eq 'K') || ($seq[$i] eq 'k') ) {$seq_m.= 'K'; next;}
		if ( ($seq[$i] eq 'M') || ($seq[$i] eq 'm') ) {$seq_m.= 'M'; next;}
		if ( ($seq[$i] eq 'B') || ($seq[$i] eq 'b') ) {$seq_m.= 'B'; next;}
		if ( ($seq[$i] eq 'D') || ($seq[$i] eq 'd') ) {$seq_m.= 'D'; next;}
		if ( ($seq[$i] eq 'H') || ($seq[$i] eq 'h') ) {$seq_m.= 'H'; next;}
		if ( ($seq[$i] eq 'V') || ($seq[$i] eq 'v') ) {$seq_m.= 'V'; next;}
		if ( ($seq[$i] eq '.') || ($seq[$i] eq '-') ) {next;}
		print STDERR "WARNING: found a non-IUPAC nucleotide ($seq[$i]) in the fasta file.\n";
	}
	return ($seq_m);
}

#crea in modo controllato una directory
sub crea_dir
{
	#variabili iniziali
	my ($dir);
	($dir) = @_;

	#creazione directory
	if (-d $dir) {print STDERR "WARNING: directory \"$dir\" already exists.\n"}
	else
	{
		my $comando = 'mkdir '.$dir;
		system($comando);
	}
}

#corregge il formato delle predizioni con 2 sole alternative
sub correggi_formato_predizione
{
	#variabili iniziali
	my ($file, $script, $model);
	($file, $script, $model) = @_;

	#apro il file per controllarne il formato
	open (I, "<$file");
	$first_line = <I>;
	chomp $first_line;
	my (@name, @value);
	my $flag = 0;
	if ($first_line eq 'V1')
	{
		#se il file è formattato male, lo memorizzo per modificarlo successivamente
		$flag = 1;
		while (my $line = <I>)
		{
			$line =~ /^(.+)\t(.+)\n$/;
			my $name = $1;
			my $value = $2;
			push (@name, $name);
			push (@value, $value);
		}
	}
	close I;

	#in caso di file formattato male, trovo a cosa si riferiscono i valori salvati e ricreo il file con la corretta formattazione
	my $new_file = $file;
	if ($flag == 1)
	{
		#trovo a cosa si riferiscono i valori salvati
		$comando = 'Rscript '.$script.' '.$model;
		my $output = `$comando`;
		my @output = split("\n", $output);
		chomp $output[1];
		my @level = split(" ", $output[1]);
		chomp $output[3];
		my @class = split(" ", $output[3]);
		chomp $output[5];
		my @prob = split(" ", $output[5]);
		my %level;
		for (my $i = 0; $i <= $#level; $i++) {$level{$level[$i]} = 0;}
		for (my $i = 0; $i <= $#class; $i++)
		{
			if ($prob[$i] > 0.5) {$level{$class[$i]}++;}
		}

		#ricreo il file con la corretta formattazione
		$stampa = '';
		foreach my $key (reverse sort { $level{$a} <=> $level{$b} } keys %level) {$stampa.= $key."\t";}
		chop $stampa;
		$stampa.= "\n";
		for (my $i = 0; $i <= $#name; $i++)
		{
			my $other_values = 1 - $value[$i];
			$stampa.= $name[$i]."\t".$value[$i]."\t".$other_values."\n";
		}
		$new_file =~ s/\.txt$/_corrected.txt/;
		open (O, ">$new_file");
		print O $stampa;
		close O;
	}
	return ($new_file);
}

#trova il valore di probabilità più elevato
sub take_best_probability
{
	#variabili iniziali
	my ($file);
	($file) = @_;

	open (I, "<$file");
	my $versions = <I>;
	my $prob = <I>;
	close I;
	chomp $versions;
	chomp $prob;
	my @version = split ("\t", $versions);
	my @prob = split ("\t", $prob);
	my %temp;
	for (my $i = 1; $i <= $#prob; $i++) {$temp{$version[$i - 1]} = $prob[$i];}
	my ($chosen_version, $best_P);
	foreach my $ver (sort { $temp{$b} <=> $temp{$a} } keys %temp)
	{
		$chosen_version = $ver;
		$best_P = $temp{$ver};
		last;
	}
	return ($chosen_version, $best_P);
}
