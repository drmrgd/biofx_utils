#!/usr/bin/env perl
# Read in either an Ion Torrent VCF or an Illumina VCF file and output a Ts/TV
# ratio, and a deamination score.  The deamination score will be the total C>T
# (or G>A) changes divded by the total changes.  
#
# 10/4/2019 - D Sims
################################################################################
use strict;
use warnings;
use autodie;

use File::Basename;
use Getopt::Long;
use Data::Dump;
use Parallel::ForkManager;
use Sort::Versions;

my $version = '2.0.111919';

# Get total number of available cores to determine how many to use by default.
# chomp(my $num_procs = sprintf("%d", qx(nproc --all)/2));

# NOTE: Disable auto parallel processing since the Mojave seems to choke on
# any attempts to parallelize things.  Probably a "security" patch for Mac OS,
# which would not be a problem on HPC systems.
my $num_procs = 0;

my $quiet = 0;
my $qflag;
($quiet == 0) ? ($qflag = "Off") : ($qflag = "On");

my $reference = "$ENV{HOME}/Dropbox/reference/hg19/hg19.fasta";
my $contexts = __make_context_hash();

my $program = basename($0);
my $description = <<"EOT";
Program to read in an Ion Torrent or Illumina VCF file and compute the transition
/ transversion ratio (Ts/Tv) and a deamination score.
EOT

my $usage = <<"EOT";
USAGE: $program [options] --source <VCF>

Options
    -s, --source     Source of the data. Allowed values are either 'ion' or 
                     'illumina'.
    -n, --numprocs   Number of parallel processes to run if running more than
                     one VCF through at a time.  Entering '0' will disable
                     parallel processing. DEFAULT: $num_procs.
    -q, --quiet      Do not print extra messages. Only print final result table.
                     DEFAULT: $qflag.
    -r, --reference  Reference genome to use for getting sequence context.
                     DEFAULT: $reference.
    -o, --outfile    Write data to file instead of just stdout.
    -v, --version    Print the version information and exit.
    -h, --help       Print this help text and exit.
EOT

my $help;
my $ver_info;
my $source;
my $outfile;

GetOptions(
    "source|s=s"    => \$source,
    "outfile|o=s"   => \$outfile,
    "numprocs|n=i"  => \$num_procs,
    "quiet|q"       => \$quiet,
    "reference|r=s" => \$reference,
    "version|v"     => \$ver_info,
    "help|h"        => \$help,
) or die "$usage";

sub help {
    print "$program - v$version\n\n$usage\n";
    exit;
}

sub version {
    print "$program - v$version\n";
    exit;
}

help() if $help;
version() if $ver_info;

die "ERROR: you must provide a VCF source!\n" unless $source;
unless (grep { $source eq $_ } qw(ion illumina)) {
    die "ERROR: valid sources are 'ion' or 'illumina'.\n";
}

unless (qx(which samtools)) {
    die "ERROR: You must have samtools installed and in your path!\n";
}

my $outfh;
if ($outfile) {
    print "Writing results to $outfile.\n";
    open($outfh, ">", $outfile)
} else {
    $outfh = \*STDOUT;
}

my @vcfs = @ARGV;
die "ERROR: You must load at least one VCF file!\n" unless @vcfs;

################################################################################
#                      End Prog Setup and CLI Arg Parsing                     #
################################################################################
my $base_data;
($num_procs > 0) 
    ? ($base_data = proc_vcfs_parallel(\@vcfs, \$source)) 
    : ($base_data = proc_vcfs(\@vcfs, \$source));

# Print it all out.
my @sbs_order = qw(C>A C>G C>T T>A T>C T>G C>T_at_CpG C>T_at_other);
my $samp_col_width = get_width([keys %$base_data]);
my $format_string = "%-${samp_col_width}s %-9s %-7s %-7s %-6s %-6s %-6s %-6s " .
    "%-6s %-6s %-12s %-13s %-20s\n";
printf {$outfh} $format_string, qw(Sample Tot_Vars Indels Ts/Tv), @sbs_order,
    "Deamination_Score";

for my $sample (sort{ versioncmp($a, $b) } keys %$base_data) {
    local $SIG{__WARN__} = sub {
        my $msg = shift;
        die "ERROR processing $sample:\n$msg\n";
    };
    local $SIG{__DIE__} = sub {
        my $msg = shift;
        die "ERROR: Can not process data from sample $sample: \n\t$msg\n";
    };

    # Print the Summary Table.
    printf {$outfh} $format_string, 
        $sample, 
        $base_data->{$sample}{'tot_vars'},
        $base_data->{$sample}{'num_indels'}, 
        $base_data->{$sample}{'titv'},
        @{$base_data->{$sample}{'sbs_6'}}{@sbs_order}, 
        $base_data->{$sample}{'deam_score'};
}

# Print a SBS-96 table for histogram plotting (through R script). Need to
# make sure it's in order first.
# TODO:  This works to print out the data, but skip for now until I can figure
# out what to do with it (probably write an R Script).
=cut
my @sorted_mut_classes;
for my $c (@sbs_order[0..5]) {
    push(@sorted_mut_classes, sort grep { $_ =~ /\[$c\]/ } keys %$contexts);
}

print join(",", "SBS", keys %$base_data), "\n";
for my $sbs (@sorted_mut_classes) {
    my @out_data = $sbs;
    for my $sample (sort { versioncmp($a, $b) } keys %$base_data) {
        push(@out_data, $base_data->{$sample}{'sbs_96'}{$sbs});
    }
    print join(",", @out_data), "\n";
}
=cut

sub proc_vcfs {
    my ($vcfs, $source) = @_;

    my %base_data;
    my @keys = qw(tot_vars num_indels titv deam_score sbs_6 sbs_96);

    for my $vcf (@$vcfs) {
        my ($sample_name, $data) = parse_vcf(\$vcf, $source);
        my ($tot_vars, $indels, $titv, $deam_score, $sbs_6, 
            $sbs_96) = get_var_metrics($data, $sample_name);

        @{$base_data{$$sample_name}}{@keys} = @{[$$tot_vars, $$indels, $$titv, 
            $$deam_score, $sbs_6, $sbs_96]};
    }
    return \%base_data;
}

sub proc_vcfs_parallel {
    my $vcfs = shift;
    my %base_data;

    my $pm = Parallel::ForkManager->new($num_procs);
    $pm->run_on_finish(
        sub {
            my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_ref) = @_;
            my @wanted = qw(tot_vars num_indels ti_count tv_count deam_events);

            my $vcf = $data_ref->{'input_file'};
            my $name = $data_ref->{'id'};
            $name //= basename($vcf);
            @{$base_data{$name}}{@wanted} = @{$data_ref}{@wanted};
        }
    );

    for my $vcf (@$vcfs) {
        $pm->start and next;
        my ($sample_name, $data) = parse_vcf(\$vcf, \$source);
        my ($tot_vars, $indels, $titv, $deam_score, $sbs_6, 
            $sbs_96) = get_var_metrics($data, $sample_name);

        $pm->finish(0,
            {
                input_file  => $vcf,
                id          => $$sample_name,
                tot_vars    => $$tot_vars,
                num_indels  => $$indels,
                titv        => $$titv,
                deam_score  => $$deam_score,
                sbs_6       => $sbs_6,
                sbs_96      => $sbs_96,
            }
        );
    }
    $pm->wait_all_children;
    return \%base_data;
}

sub parse_vcf {
    my ($vcf, $source) = @_;
    print "Processing $$vcf...\n" unless $quiet;
    ($$source eq 'ion')
        ? return run_ion_parser($vcf)
        : return run_illumina_parser($vcf);
}

sub run_illumina_parser {
    my $vcf = shift;
    my @vcf_data;

    # In Illumina files, there is paired tumor normal for some of these, but not
    # always.  Also, not always easy to figure out which is which from sample
    # field in VCF.  Rajesh recommneds using the file name for this instead.
    # my $sample_name = __get_sample_name($vcf);
    my $sample_name = (split(/\./, basename($$vcf)))[0];

    my $tumor_index;
    open(my $fh, "<", $$vcf);
    while (<$fh>) {
        # Since we sometimes have paired tumor normal, need to figure out which
        # column to grab up. Use $sample_name.
        next unless /^#CHROM/ || /^chr/;
        chomp(my @fields = split(/\t/));
        if (/^#CHROM/) {
            ($tumor_index) = grep { $fields[$_] eq $sample_name } 0..$#fields;
            next;
        }

        # Filter field informative in MuTect2 data; Keep only "PASS"
        next unless $fields[6] eq 'PASS';

        # TODO: need to fix this to work more precisely. For now just skip the
        # multi-allele entries.
        next if grep { $_ =~ /,/ } @fields[3,4];

        # Get VAF for filtering.
        my ($vaf, $ref_cov, $alt_cov) = __get_vaf(\@fields, $tumor_index);
        next if $vaf < 1;

        push(@vcf_data, ["$fields[0]:$fields[1]", "$fields[3]>$fields[4]"]);
    }
    return \$sample_name, \@vcf_data;
}

sub run_ion_parser {
    my $vcf = shift;
    my @vcf_data;

    my $sample_name = __get_sample_name($vcf);

    # Use VCF extractor here, despite the overhead, to hand complex variant
    # cases in the VCF.  Else, we can just get by parsing the REF / ALT cols.
    my $cmd = "vcfExtractor.pl -CNn $$vcf";
    open(my $stream, "-|", $cmd);

    while (<$stream>) {
        next unless /^chr/;
        local $SIG{__WARN__} = sub {
            my $msg = shift;
            warn "Issue with the following entry:\n";
            print $_;
        };
        chomp(my @data = split(/,/));

        # Filter out sub 1% variants.
        next if $data[3] < 1;
        push(@vcf_data, [$data[0], "$data[1]>$data[2]"]);
    }
    return \$sample_name, \@vcf_data;
}

sub get_var_metrics {
    my ($vcf_data, $sample_name) = @_;

    my ($tot_vars, $indels) = (0,0);

    my %sbs_96 = %$contexts;
    for my $data (@$vcf_data) {
        
        my ($coord, $change) = @$data;
        $tot_vars++;

        # Count, but don't process, indels.
        do { $indels++; next} unless length($change) == 3;

        my ($chr, $pos) = split(/:/, $coord);
        my $query = sprintf("%s:%s-%s", 
            ($chr =~ /^chr/) ? $chr : "chr$chr",
            $pos-1,
            $pos+1);

        my $context;
        if ($change =~ /^[CT]/) {
            # Change is one of the 6 types; don't need to get a reverse comp.
            $context = get_context($query);
        } else {
            # Need to get the reverse comp to match one of the 6 sigs.
            $context = rev_comp(get_context($query));
            my $comp = rev_comp($change);
            $change = join('>', reverse(split('', $comp)));
        }
        my ($first, $second, $third) = split('', $context);
        my $c_string = "${first}[$change]${third}";

        if (! exists $sbs_96{$c_string}) {
            print "Can't find '$c_string' in the 96 patterns hash($coord;$change ",
                "=> $context -> $c_string)\n";
        } else {
            $sbs_96{$c_string}++;
        }
    }

    unless ($quiet) {
        my $count = 0;
        $count += $sbs_96{$_} for keys %$contexts;
        print "\tFound $tot_vars variants in $$sample_name:\n";
        print "\t  SNVs:          " .  ($tot_vars - $indels) . "\n";
        print "\t  Indels:        $indels\n";
        print "\tComputing SBS-6 and variant motifs...";
    }

    my ($tstv, $var_motifs, $transitions, $transversions) = get_motifs(\%sbs_96);
    my $deam_score = sprintf("%0.2f", 
        $var_motifs->{'C>T_at_other'} / $var_motifs->{'C>T_at_CpG'});

    unless ($quiet) {
        print "Done!\n";
        print "\t  Transitions:   $transitions\n";
        print "\t  Transversions: $transversions\n";
        print "\t  Ts/Tv:         $tstv\n";
        print "\t  Deam Score:    $deam_score\n";
    }
    return \$tot_vars, \$indels, \$tstv, \$deam_score, $var_motifs, \%sbs_96;
}

sub get_motifs {
    my $var_contexts = shift;

    my ($transitions, $transversions) = (0,0);
    my %motifs = (
        'C>A'           => 0,
        'C>T'           => 0,
        'C>G'           => 0,
        'T>A'           => 0,
        'T>C'           => 0,
        'T>G'           => 0,
        'C>T_at_CpG'    => 0,
        'C>T_at_other'  => 0,
    );

    for my $context (keys %$var_contexts) {
        my ($first_base, $change, $last_base) =
            $context =~ /([ACTG])\[([ACTG]>[ACTG])\]([ACTG])/;
        $motifs{$change} += $var_contexts->{$context};
        if ($change eq 'C>T') {
            $transitions += $var_contexts->{$context};
            ($context =~ /G$/)
                ? ($motifs{"C>T_at_CpG"} += $var_contexts->{$context})
                : ($motifs{'C>T_at_other'} += $var_contexts->{$context});
        }
        elsif ($change eq 'T>C') {
            $transitions += $var_contexts->{$context};
        } else {
            $transversions += $var_contexts->{$context};
        }
    }
    my $tstv = sprintf("%0.3f", $transitions / $transversions);
    return $tstv, \%motifs, $transitions, $transversions;
}

sub get_context {
    my $query = shift;
    open(my $stream, "-|", "samtools faidx $reference $query");
    return(map{chomp; $_}<$stream>)[1];
}

sub rev_comp {
    my $str = shift;
    my %comp = ('A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A');
    my $new_str = '';
    my @bases = grep { /[ACTG]/} split('', $str);
    $new_str .= $comp{$_} for @bases;
    return reverse($new_str);
}

sub __get_sample_name {
    my $vcf = shift;
    open(my $fh, "<", $$vcf);
    while (<$fh>) {
        next unless /^#CHROM/;
        chomp(my $name = (split(/\t/))[-1]);
        return $name;
    }
    close $fh; #shouldn't get here and need this, but just in case....
}

sub __get_vaf {
    my ($vcf_elems, $info_index) = @_;
    $info_index //= 9;

    my %format;
    @format{split(/:/, $vcf_elems->[8])} = split(/:/, $vcf_elems->[$info_index]);

    my ($ref_reads, $alt_reads);

    # If we have raw Ion VCF, then need to get different fields than rest.
    if (exists $format{'FAO'}) {
        $ref_reads = $format{'FRO'};
        $alt_reads = $format{'FAO'};
    } else {
        ($ref_reads, $alt_reads) = split(/,/, $format{'AD'});
    }
    return (sprintf("%0.2f", ($alt_reads / ($alt_reads + $ref_reads)) * 100),
         $ref_reads, $alt_reads);
}

sub __make_context_hash {
    my %contexts;
    my @bases = qw(A C T G);
    my @mut_classes = qw([C>A] [C>G] [C>T] [T>A] [T>C] [T>G]);
    for my $first_base (@bases) {
        for my $class (@mut_classes) {
            for my $third_base (@bases) {
                $contexts{"${first_base}${class}${third_base}"} = 0;
            }
        }
    }
    return \%contexts;
}

sub get_width {
    my $strs = shift;
    my $longest = 0;
    for (@$strs) {
        $longest = length($_) if length($_) > $longest;
    }
    return $longest + 2;
}
