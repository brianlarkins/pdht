#!/usr/bin/perl
use strict 'vars';

sub partition {
	return -1 unless (scalar (@_) > 0);

	# Find the array index at which the deviation is greatest.
	my $max_deviation_index = -1;
	my $max_deviation = -1;

	my $i;
	for ($i = 0; $i < scalar(@_) - 1; $i++) {
		my $current_deviation = $_[$i + 1] - $_[$i];
		# print "Current deviation:  $current_deviation\n";

		if ($current_deviation > $max_deviation) {
			$max_deviation_index = $i;
			$max_deviation = $current_deviation;
		}
	}

	return $max_deviation_index;
}


sub do_partition {
	# Sort the input, so that we can compute deviations between
	# consecutive elements.
	@_ = sort {$a <=> $b} @_;

	my $p = partition(@_);
	my $i;
	for ($i = 0; $i <= $p; $i++) {
		print "$_[$i] ";
	}
	print "\n";
	for (; $i < scalar (@_); $i++)  {
		print "$_[$i] ";
	}
	print "\n";
}


my @data = (1, 2, 5, 6, 7);

do_partition(@data);
print "-------------------\n";
do_partition( (-3, 50, -39, 1, 2, 3) );
