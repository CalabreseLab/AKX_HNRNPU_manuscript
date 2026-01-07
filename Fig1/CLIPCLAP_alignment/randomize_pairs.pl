#!/usr/bin/env perl
use strict;
use warnings;


#------------------------------------------------------------
# set_orientation(flag, is_read1, r1_rev, r2_rev)
#
# Input:
#   flag     - original FLAG value for this read
#   is_read1 - 1 if this read is R1, 0 if it's R2
#   r1_rev   - 1 if R1 should be reverse-complement strand (i.e. '-' / FLAG 0x10)
#   r2_rev   - 1 if R2 should be reverse-complement strand
#
# Behavior:
#   We ONLY touch bit 0x10 (this read reverse) and 0x20 (mate reverse).
#   All other bits (paired, proper pair, read1/read2, etc.) are preserved.
#
#   For R1:
#       - 0x10 says if R1 itself is reverse.
#       - 0x20 says if R2 (mate) is reverse.
#   For R2:
#       - 0x10 says if R2 itself is reverse.
#       - 0x20 says if R1 (mate) is reverse.
#
# Output:
#   New FLAG integer with consistent strand bits for this read.
#------------------------------------------------------------
sub set_orientation {
    my ($flag, $is_read1, $r1_rev, $r2_rev) = @_;

    # Clear the strand-related bits:
    #  0x10 = this read is reverse
    #  0x20 = mate is reverse
    # Using bitwise AND with the inverse (~) to unset them.
    $flag &= ~(0x10 | 0x20);

    if ($is_read1) {
        # This read is R1
        $flag |= 0x10 if $r1_rev;   # set "this read reverse" if R1 is reverse
        $flag |= 0x20 if $r2_rev;   # set "mate reverse" if R2 is reverse
    } else {
        # This read is R2
        $flag |= 0x10 if $r2_rev;   # set "this read reverse" if R2 is reverse
        $flag |= 0x20 if $r1_rev;   # set "mate reverse" if R1 is reverse
    }

    return $flag;
}

#------------------------------------------------------------
# %pairs will temporarily store reads keyed by QNAME.
# Each key accumulates its two reads (R1 + R2), then we process and flush.
#------------------------------------------------------------
my %pairs;

# Counters for how many pairs end up with R1 on + vs -.
# We count per PAIR (based on R1 orientation), not per read.
my $r1_pos = 0;
my $r1_neg = 0;
my $orphans = 0;

#------------------------------------------------------------
# Main loop: read SAM from STDIN line by line
#------------------------------------------------------------
while (my $line = <STDIN>) {

    # Header lines (start with '@') are copied directly and not modified.
    if ($line =~ /^\@/) {
        print $line;
        next;
    }

    # For alignment lines: strip trailing newline and split into fields.
    chomp $line;
    my @f = split(/\t/, $line);

    # By SAM spec: field 0 = QNAME, field 1 = FLAG
    my ($qname, $flag) = @f[0,1];

    # Store this alignment under its QNAME.
    # We store a reference to @f so we can modify FLAG later.
    push @{ $pairs{$qname} }, \@f;

    # Once we have exactly 2 records for this QNAME, we assume:
    #   - they are the two mates (R1 and R2),
    #   - input is guaranteed properly paired + primary (your assumption).
    if (@{ $pairs{$qname} } == 2) {

        # Get the two reads for this pair.
        my ($a, $b) = @{ $pairs{$qname} };

        # Identify which is R1 and which is R2 via FLAG bits:
        #  0x40 -> this read is first in pair (R1)
        #  0x80 -> this read is second in pair (R2)
        my ($r1, $r2);
        if ($a->[1] & 0x40) {
            # read a is R1
            ($r1, $r2) = ($a, $b);
        } else {
            # otherwise b must be R1 (given the guarantees)
            ($r1, $r2) = ($b, $a);
        }

        #--------------------------------------------------------
        # Random assignment of strand for the pair:
        #
        # We flip a coin (0 or 1):
        #   flip = 0: R1 on +,  R2 on -
        #   flip = 1: R1 on -,  R2 on +
        #
        # This keeps them opposite, like a typical proper FR/RF layout,
        # but randomizes which strand R1 ends up on.
        #--------------------------------------------------------
        my $flip   = int(rand(2));
        my $r1_rev = $flip ? 1 : 0;       # 1 means R1 reverse ('-'), 0 means forward ('+')
        my $r2_rev = $r1_rev ? 0 : 1;     # mate is always opposite

        # Count number of pairs by final R1 orientation:
        if ($r1_rev) {
            $r1_neg++;
        } else {
            $r1_pos++;
        }

        #--------------------------------------------------------
        # Update FLAGs for R1 and R2
        #
        # We do NOT touch pairing bits, only 0x10 / 0x20 via set_orientation.
        #   - For R1: is_read1 = 1
        #   - For R2: is_read1 = 0
        #--------------------------------------------------------
        my $new_flag_r1 = set_orientation($r1->[1], 1, $r1_rev, $r2_rev);
        my $new_flag_r2 = set_orientation($r2->[1], 0, $r1_rev, $r2_rev);

        # Write new FLAG values back into the field arrays.
        $r1->[1] = $new_flag_r1;
        $r2->[1] = $new_flag_r2;

        #--------------------------------------------------------
        # Output:
        # We must preserve the original order of lines in the file.
        # $a and $b are "first seen" vs "second seen".
        # So we print them in that order, but with updated FLAGs.
        #--------------------------------------------------------
        print join("\t", @$a), "\n";
        print join("\t", @$b), "\n";

        # Done with this QNAME, free memory.
        delete $pairs{$qname};
    }
}

#------------------------------------------------------------
# Safety: if anything is unexpectedly left in %pairs (e.g., odd count),
# we just print those reads unchanged. This should not happen under
# your "properly paired" guarantee, but it's a harmless fallback.
#------------------------------------------------------------

if (scalar keys %pairs) {
    $orphans = scalar keys %pairs;
    warn "# WARNING: $orphans read(s) left unpaired — check input integrity!\n";

    # Output them unchanged just in case
    for my $recs (values %pairs) {
        for my $r (@$recs) {
            print join("\t", @$r), "\n";
        }
    }
}
#------------------------------------------------------------
# Print summary stats to STDERR so as not to corrupt SAM output.
# These tell you, for this input file, how many pairs ended up
# with R1 on + vs - after randomization.
#------------------------------------------------------------
warn "# R1 assigned to + strand: $r1_pos\n";
warn "# R1 assigned to - strand: $r1_neg\n";
