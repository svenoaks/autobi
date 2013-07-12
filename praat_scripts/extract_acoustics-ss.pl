#!/usr/bin/perl

# Run without arguments for a help screen.

use strict;

my $UNDEF = 'undef';

# praat script for extracting standard acoustics.
my $ExStAc_PRAAT_SCRIPT = "/proj/speech/tools/speechlab/praat_scripts/extractStandardAcoustics2.praat";

# path to praat (linux only!)
#my $PRAAT = "/proj/speech/tools/praat/linux/praat";
my $PRAAT = "/Applications/Praat.app/Contents/MacOS/Praat";

# wav file (full path), minimum pitch and maximum pitch
my $wav_file = shift;
my $min_pitch = shift;
my $max_pitch = shift;

my $TRIM_PERCENTAGE = shift;

# sanity checks
die "Usage: extract_acoustics.pl WavFile MinPitch MaxPitch Trim% From_1 To_1 From_2 To_2...\n".q^
 This script receives several parameters:
   WavFile   full path of the wav file to process
   MinPitch  speaker's estimated minimum pitch (eg, 50 for male, 75 for female)
   MaxPitch  speaker's estimated maximum pitch (eg, 300 for male, 500 for female)
   Trim%     percentage of the data to trim on each side (eg, 5 will 
             trim the upper 5% and the lower 5% --or 10% in total).
   From_1    start point of interval 1
   To_1      end point of interval 1
   From_2    start point of interval 2
   To_2      end point of interval 2
   ...
 and computes a set of standard acoustic features over the given
 intervals only, ignoring the rest of the wav file.
 The intervals must not be empty (i.e. From_i < To_i), and they must
 be sorted incrementally (i.e. To_i <= From_i+1).
 There may be an arbitrary number of intervals, but at least one.

 Example:
 ./extract_acoustics.pl /proj/speech/projects/games/data/session_01/s01.objects.0.A.wav 50 300 5 69.41 70.66 75.47 76.03
^
unless ($wav_file and $min_pitch and $max_pitch and $TRIM_PERCENTAGE ne '');
die "File not found: '$wav_file'\n" unless -e $wav_file;
die "File path must be absolute: '$wav_file'\n" unless substr($wav_file, 0, 1) eq "/";

# (from_1, to_1, from_2, to_2, from_3, to_3, ...)
my @points = @ARGV;

my @from;
my @to;
for (my $i=0; $i<@points; $i++) {
	if ($i % 2 == 0) {
		push @from, $points[$i];
	}
	else {
		push @to, $points[$i];
	}
}
my $num_intervals = @to+0;


# more sanity checks
for my $i (0..$num_intervals-1) {
	die "Error: invalid interval\n" 
	unless $from[$i] < $to[$i];
}
for my $i (0..$num_intervals-2) {
	die "Error: intervals must be sorted\n" 
	unless $to[$i] <= $from[$i+1];
}

# - - - - - - - - - - - - - - - - - - - - - - - - -

# Now prepare the praat script to chunk all the intervals
# and concatenate them together.

# stem for temporary files
my $tmp = "/tmp/acoustics-".time;
	
open FILE, ">$tmp.praat";

print FILE "Open long sound file... $wav_file\n";
print FILE "Rename... main\n";

for my $i (0..$num_intervals-1) {
	print FILE "select LongSound main\n";
	print FILE "Extract part... $from[$i] $to[$i] yes\n";
	print FILE "Rename... part_$i\n";
}

print FILE "select Sound part_0\n";
for my $i (1..$num_intervals-1) {
	print FILE "plus Sound part_$i\n";
}

print FILE qq^
Concatenate
Rename... sound

select Sound sound
endTime = Get end time
To Pitch... 0 $min_pitch $max_pitch
vcd2tot_frames = undefined
vcd_frames = Count voiced frames
tot_frames = Get number of frames
if tot_frames > 0
	vcd2tot_frames = vcd_frames / tot_frames
endif
Down to PitchTier
Write to text file... $tmp.PitchTier

select Sound sound
To Intensity... $min_pitch 0 yes
Down to IntensityTier
Write to text file... $tmp.IntensityTier

printline 'vcd2tot_frames:3'
^;
close FILE;


# Run the script we just wrote, to create the PraatTiers of all
# the chunks concatenated, and get the voiced2total frames ratio
my $vcd2tot_frames = `$PRAAT $tmp.praat`;
chomp $vcd2tot_frames;


# Read the praat tiers
my %PitchTier     = read_PraatTier("$tmp.PitchTier");
my %IntensityTier = read_PraatTier("$tmp.IntensityTier");

# Remove the top and bottom $TRIM_PERCENTAGE of the data
# from the Pitch and Intensity tiers.
trim_PraatTier(\%PitchTier,     $TRIM_PERCENTAGE);
trim_PraatTier(\%IntensityTier, $TRIM_PERCENTAGE);

my @pch_mmms = max_min_mean_stdev(values %PitchTier);
my @int_mmms = max_min_mean_stdev(values %IntensityTier);

print "duration:--undefined--\n";
print "F0_MIN:".	nf($pch_mmms[1])."\n";
print "F0_MAX:".	nf($pch_mmms[0])."\n";
print "F0_MEAN:".	nf($pch_mmms[2])."\n";
print "F0_MEDIAN:--undefined--\n";
print "F0_STDV:".	nf($pch_mmms[3])."\n";
print "F0_MAS:--undefined--\n";
print "ENG_MAX:".	nf($int_mmms[0])."\n";
print "ENG_MIN:".	nf($int_mmms[1])."\n";
print "ENG_MEAN:".	nf($int_mmms[2])."\n";
print "ENG_STDV:".	nf($int_mmms[3])."\n";
print "VCD2TOT_FRAMES:".nf($vcd2tot_frames)."\n";

# Remove the temporary files.
# `rm $tmp.*`;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# number format
sub nf {
	my $n = shift;
	if ($n =~ m/undef/) {
		return '--undefined--';
	}
	else {
		return sprintf("%.3f", $n);
	}
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# returns a hash: time-->value read from a PraatTier file
sub read_PraatTier {
	# full path to the file
	my $filename = shift;

	# read the times
	my @times = `cat $filename | grep number | sed "s/.* = //"`;
	chomp @times;
    
    print STDERR "@times";

	# read the values
	my @values = `cat $filename | grep value | sed "s/.* = //"`;
	chomp @values;
	
	my %res;

	if (@values != @times) {
		print STDERR "Error reading '$filename'. Mismatch in the number of 'time' and 'value' fields.\n";
	}
	else {
		for (my $i=0; $i<@times; $i++) {
			$res{$times[$i]} = $values[$i];
		}
	}
	
	return %res;
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Discards a percentage of the data (most likely outliers) from a PraatTier hash.
# Call: trim_PraatTier(\@PraatTier, 5);
sub trim_PraatTier {
	# Reference to the PraatTier hash to trim. The hash has this form: time->value
	my $data = shift;

	# Percentage of the data to trim on each side before doing the calculations.
	# Example: if $alpha==5, then the upper 5% and the lower 5% (10% in total) 
	# are discarded, and then the max, min, mean, stdev are computed.
	my $alpha = shift;

	# - - -

	my @times  = keys %$data;

	my $n = @times+0;

	# sort the times according to the @values
	@times = sort { $$data{$a} <=> $$data{$b} } @times;

	# calculate the number of data points to discard...
	my $discard = int($n * $alpha / 100);

	# find out the times of the values that should be discarded
	my @discarded_times;
	for (1..$discard) {
		push @discarded_times, (shift @times);
		push @discarded_times, (pop @times);
	}

	# remove the discarded entries 	from the hash.
	for my $time (@discarded_times) {
		delete($$data{$time});
	}
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Returns the (max, min, mean, stdev) of the given array.
# Call: max_min_mean_stdev(@data);
sub max_min_mean_stdev {
	# array with the data
	my @data = @_;

	my $n = @data+0;

	my ($max, $min, $mean, $stdev);

	if ($n>0) {
		# calculate the max, min and mean
		$max = $data[0];
		$min = $data[0];
		my $sum = 0;

		for my $x (@data) {
			$x+=0;

			if ($x > $max) {$max = $x}
			if ($x < $min) {$min = $x}

			$sum += $x;
		}
		$mean = $sum/$n;
		
		if ($n>1) {
			# calculate the stdev
			$sum = 0;
			for my $x (@data) {
				my $delta = $x - $mean;
				$sum += $delta * $delta;
			}
			$stdev = sqrt($sum/($n-1));
		}
		else {
			$stdev = $UNDEF;
		}
	}
	else {
		($max, $min, $mean, $stdev) = ($UNDEF, $UNDEF, $UNDEF, $UNDEF);
	}
	
	return ($max, $min, $mean, $stdev);
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
