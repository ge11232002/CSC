############################################################################
# Primerfinder.pm - perl module for handeling a set of composite alignments
#
# Copyright
############################################################################

package AT::Tools::PrimerFinder;

use strict;
use Carp;
use Bio::LocatableSeq;
use AT::Tools::Run::Primer3;

=head1 NAME

Primerfinder - module for designing PCR primers and sequencing primers.

=head1 SYNOPSIS

    use Primerfinder;
    my $primerfinder = Primerfinder -> new(-seq         => Bio::Seq object,
				       -p_min_tm    => $p_min_tm,
				       -p_max_tm    => $p_max_tm,
				       -p_opt_tm    => $p_opt_tm,
				       -p_gc        => $p_gc,
				       -p_optl      => $pcr_optl,
				       -s_min_tm    => $s_min_tm,
				       -s_max_tm    => $s_max_tm,
				       -s_opt_tm    => $s_opt_tm,
				       -s_gc        => $s_gc,
				       -s_optl      => $s_optl,
				       -fp          => $fp,
				       -xp          => $xp,
				       -xs          => $xs,
				       -ir          => $ir,
				       -pcrl        => $pcrl,
				       -sl          => $sl,
                                       -non_ol      => $non_ol,
				       -d           => $d,
                                       -s_nret      => $s_nret,
				      );

    $primerfinder -> get_seq();
    $primerfinder -> get_p_min_tm();
    $primerfinder -> get_p_max_tm();
    $primerfinder -> get_p_opt_tm();
    $primerfinder -> get_p_gc();
    $primerfinder -> get_p_optl();
    $primerfinder -> get_s_min_tm();
    $primerfinder -> get_s_max_tm();
    $primerfinder -> get_s_opt_tm();
    $primerfinder -> get_s_gc();
    $primerfinder -> get_s_optl();
    $primerfinder -> get_fp();
    $primerfinder -> get_xp();
    $primerfinder -> get_xs();
    $primerfinder -> get_pcrl();
    $primerfinder -> get_sl();
    $primerfinder -> get_d();
    $primerfinder -> get_pcr_nret();
    $primerfinder -> get_s_nret();
    $primerfinder -> pcr_primers(-target => $target);
    $primerfinder -> seq_primers(-target => $target, -pp => $pp);


=head1 DESCRIPTION

This method creates a new Primerfinder object, and sets all parameters for the primer design.

=head1 METHODS DESCRIPTION

=cut

############################################################################

=head2 B<new()>

I<Title:> new()

I<Usage:>

$primerfinder = Primerfinder -> new(-seq         => Bio::Seq object,
				    -p_min_tm    => $p_min_tm,
				    -p_max_tm    => $p_max_tm,
				    -p_opt_tm    => $p_opt_tm,
				    -p_gc        => $p_gc,
                                    -p_optl      => $p_optl,
				    -s_min_tm    => $s_min_tm,
				    -s_max_tm    => $s_max_tm,
				    -s_opt_tm    => $s_opt_tm,
				    -s_gc        => $s_gc,
                                    -s_optl      => $s_optl,
				    -fp          => $fp,
				    -xp          => $xp,
				    -xs          => $xs,
                                    -ir          => $ir,
				    -pcrl        => $pcrl,
				    -sl          => $sl,
                                    -d           => $d,
                                    -pcr_nret    => $pcr_nret,
				    -s_nret      => $s_nret,
				     );

I<Input:>

          -seq        - The B<Bio::Seq> object in which primers should
                        be picked.
	  -p_min_tm   - Minimal tm value for pcr primer
	  -p_max_tm   - Maximal tm value for pcr primer
	  -p_opt_tm   - Optimal tm value for pcr primer
	  -p_gc       - length of 3'gc clamp for pcr primer
          -p_optl     - Optimal length of pcr primer
	  -s_min_tm   - Minimal tm value for pcr primer
	  -s_max_tm   - Maximal tm value for pcr primer
	  -s_opt_tm   - Optimal tm value for pcr primer
	  -s_gc       - Lenth of 3'gc clamp for sequencing primer
          -s_optl     - Optimal length of sequencing primer
	  -fp         - Accepted number of 3 prime false priming nucleotides
                        in primer (i.e. when a primer is searched against
			all possible "words" in the search sequence)
	  -xp         - Extra sequence outside of the specified target region
                        that should be included in the region to be amplified
                        using the PCR primers.
          -xs         - Extra sequence outside of the specified target
                        region that should be included in the region for
                        which Sequncing primers are designed (Since it is
                        difficult to get high quality sequence direclty after
                        the primre end)
          -ir         - Number of base pairs outside of the target + extra
                        region that should be included in the sequence used
                        for primer design. This will be used to set the
                        "Included region" in Primer 3. Default value is 2500.
          -d          - Minimum distance between first position of a PCR
                        primer and the first position of the corresponding
                        sequencing primer
	  -pcrl       - Range of pcr product length, format: \"lower
                        boundary-upper boundary\"
	  -sl         - Range of seq product lenth (distance between a pair
                        of sequence primers) , format: \"lower boundary-upper
                        boundary\"
	  -pcr_nret   - Numbers of PCR primers returned by primer 3 for
                        each region
	  -s_nret     - Numbers of sequencing primers returned by primer 3
                        for each region

I<Output:>

Returns a new B<PrimerFinder> object.

I<Description:>

=cut

sub new {
    my ( $class, %args ) = @_;

    #New prmerfactory...

    #Default values for Primersearch run:
    my $fp = 10;
    my $xp = 100;
    my $xs = 50;
    my $ir = 2500;
    my $d  = 15;

#These parameters have default values in Primer 3 which will be used if the parameters are not set here:
    my (
        $seq,    $p_min_tm, $p_max_tm, $p_opt_tm, $p_gc,
        $p_optl, $s_min_tm, $s_max_tm, $s_opt_tm, $s_gc,
        $s_optl, $pcrl,     $sl,       $pcr_nret, $s_nret
    );
    if ( defined $args{-seq} ) {
        $seq = $args{-seq};
        delete $args{-seq};
    }
    if ( defined $args{-p_min_tm} ) {
        $p_min_tm = $args{-p_min_tm};
        delete $args{-p_min_tm};
    }
    if ( defined $args{-p_max_tm} ) {
        $p_max_tm = $args{-p_max_tm};

        delete $args{-p_max_tm};
    }
    if ( defined $args{-p_opt_tm} ) {
        $p_opt_tm = $args{-p_opt_tm};

        delete $args{-p_opt_tm};
    }
    if ( defined $args{-p_gc} ) {
        $p_gc = $args{-p_gc};
        delete $args{-p_gc};
    }
    if ( defined $args{-p_optl} ) {
        $p_optl = $args{-p_optl};
        delete $args{-p_optl};
    }
    if ( defined $args{-s_min_tm} ) {
        $s_min_tm = $args{-s_min_tm};
        delete $args{-s_min_tm};
    }
    if ( defined $args{-s_max_tm} ) {
        $s_max_tm = $args{-s_max_tm};
        delete $args{-s_max_tm};
    }
    if ( defined $args{-s_opt_tm} ) {
        $s_opt_tm = $args{-s_opt_tm};
        delete $args{-s_opt_tm};
    }
    if ( defined $args{-s_gc} ) {
        $s_gc = $args{-s_gc};
        delete $args{-s_gc};
    }
    if ( defined $args{-s_optl} ) {
        $s_optl = $args{-s_optl};
        delete $args{-s_optl};
    }
    if ( defined $args{-xp} ) {
        $xp = $args{-xp};
        delete $args{-xp};
    }
    if ( defined $args{-xs} ) {
        $xs = $args{-xs};
        delete $args{-xs};
    }

    if ( defined $args{-pcrl} ) {
        $pcrl = $args{-pcrl};
        delete $args{-pcrl};
    }
    if ( defined $args{-sl} ) {
        $sl = $args{-sl};
        delete $args{-sl};
    }
    if ( defined $args{-sl} ) {
        $sl = $args{-sl};
        delete $args{-sl};
    }
    if ( defined $args{-pcr_nret} ) {
        $pcr_nret = $args{-pcr_nret};
        delete $args{-pcr_nret};
    }
    if ( defined $args{-s_nret} ) {
        $s_nret = $args{-s_nret};
        delete $args{-s_nret};
    }
    return bless {
        _seq      => $seq,
        _p_min_tm => $p_min_tm,
        _p_max_tm => $p_max_tm,
        _p_opt_tm => $p_opt_tm,
        _p_gc     => $p_gc,
        _p_optl   => $p_optl,
        _s_min_tm => $s_min_tm,
        _s_max_tm => $s_max_tm,
        _s_opt_tm => $s_opt_tm,
        _s_gc     => $s_gc,
        _s_optl   => $s_optl,
        _fp       => $fp,
        _xp       => $xp,
        _xs       => $xs,
        _pcrl     => $pcrl,
        _sl       => $sl,
        _ir       => $ir,
        _d        => $d,
        _pcr_nret => $pcr_nret,
        _s_nret   => $s_nret,
      },
      ref $class || $class;
}

=head2 B<get_seq()>

I<Title:> get_seq()

I<Usage:> $primerfinder -> get_seq();

I<Input:> No input

I<Output:> Returns the -seq attribute

I<Description:>

=cut

sub get_seq { $_[0]->{_seq} }

=head2 B<get_p_min_tm()>

I<Title:> get_p_min_tm()

I<Usage:> $primerfinder -> get_p_min_tm();

I<Input:> No input

I<Output:> Returns the -p_min_tm attribute

I<Description:>

=cut

sub get_p_min_tm { $_[0]->{_p_min_tm} }

=head2 B<get_p_max_tm()>

I<Title:> get_p_max_tm()

I<Usage:> $primerfinder -> get_p_max_tm();

I<Input:> No input

I<Output:> Returns the -p_max_tm attribute

I<Description:>

=cut

sub get_p_max_tm { $_[0]->{_p_max_tm} }

=head2 B<get_p_opt_tm()>

I<Title:> get_p_opt_tm()

I<Usage:> $primerfinder -> get_p_opt_tm();

I<Input:> No input

I<Output:> Returns the -p_opt_tm attribute

I<Description:>

=cut

sub get_p_opt_tm { $_[0]->{_p_opt_tm} }

=head2 B<get_p_gc()>

I<Title:> get_p_gc()

I<Usage:> $primerfinder -> get_p_gc();

I<Input:> No input

I<Output:> Returns the -p_gc attribute

I<Description:>

=cut

sub get_p_gc { $_[0]->{_p_gc} }

=head2 B<get_p_optl()>

I<Title:> get_p_optl()

I<Usage:> $primerfinder -> get_p_optl();

I<Input:> No input

I<Output:> Returns the -p_optl attribute

I<Description:>

=cut

sub get_p_optl { $_[0]->{_p_optl} }

=head2 B<get_s_min_tm()>

I<Title:> get_s_min_tm()

I<Usage:> $primerfinder -> get_s_min_tm();

I<Input:> No input

I<Output:> Returns the -s_min_tm attribute

I<Description:>

=cut

sub get_s_min_tm { $_[0]->{_s_min_tm} }

=head2 B<get_s_max_tm()>

I<Title:> get_s_max_tm()

I<Usage:> $primerfinder -> get_s_max_tm();

I<Input:> No input

I<Output:> Returns the -s_max_tm attribute

I<Description:>

=cut

sub get_s_max_tm { $_[0]->{_s_max_tm} }

=head2 B<get_s_opt_tm()>

I<Title:> get_s_opt_tm()

I<Usage:> $primerfinder -> get_s_opt_tm();

I<Input:> No input

I<Output:> Returns the -s_opt_tm attribute

I<Description:>

=cut

sub get_s_opt_tm { $_[0]->{_s_opt_tm} }

=head2 B<get_s_gc()>

I<Title:> get_s_gc()

I<Usage:> $primerfinder -> get_s_gc();

I<Input:> No input

I<Output:> Returns the -s_gc attribute

I<Description:>

=cut

sub get_s_gc { $_[0]->{_s_gc} }

=head2 B<get_s_optl()>

I<Title:> get_s_optl()

I<Usage:> $primerfinder -> get_s_optl();

I<Input:> No input

I<Output:> Returns the -s_optl attribute

I<Description:>

=cut

sub get_s_optl { $_[0]->{_s_optl} }

=head2 B<get_fp()>

I<Title:> get_fp()

I<Usage:> $primerfinder -> get_fp();

I<Input:> No input

I<Output:> Returns the -fp attribute

I<Description:>

=cut

sub get_fp { $_[0]->{_fp} }

=head2 B<get_xp()>

I<Title:> get_xp()

I<Usage:> $primerfinder -> get_xp();

I<Input:> No input

I<Output:> Returns the -xp attribute

I<Description:>

=cut

sub get_xp { $_[0]->{_xp} }

=head2 B<get_xs()>

I<Title:> get_xs()

I<Usage:> $primerfinder -> get_xs();

I<Input:> No input

I<Output:> Returns the -xs attribute

I<Description:>

=cut

sub get_xs { $_[0]->{_xs} }

=head2 B<get_pcrl()>

I<Title:> get_pcrl()

I<Usage:> $primerfinder -> get_pcrl();

I<Input:> No input

I<Output:> Returns the -pcrl attribute

I<Description:>

=cut

sub get_pcrl { $_[0]->{_pcrl} }

=head2 B<get_sl()>

I<Title:> get_sl()

I<Usage:> $primerfinder -> get_sl();

I<Input:> No input

I<Output:> Returns the -sl attribute

I<Description:>

=cut

sub get_sl { $_[0]->{_sl} }

=head2 B<get_d()>

I<Title:> get_d()

I<Usage:> $primerfinder -> get_d();

I<Input:> No input

I<Output:> Returns the -d attribute

I<Description:>

=cut

sub get_d { $_[0]->{_d} }

=head2 B<get_pcr_nret()>

I<Title:> get_pcr_ntet()

I<Usage:> $primerfinder -> get_pcr_nret();

I<Input:> No input

I<Output:> Returns the -pcr_nret attribute

I<Description:>

=cut

sub get_pcr_nret { $_[0]->{_pcr_nret} }

=head2 B<get_s_nret()>

I<Title:> get_s_ntet()

I<Usage:> $primerfinder -> get_s_nret();

I<Input:> No input

I<Output:> Returns the -s_nret attribute

I<Description:>

=cut

sub get_s_nret { $_[0]->{_s_nret} }

=head2 B<pcr_primers()>

I<Title:> pcr_primers()

I<Usage:> $primerfinder -> pcr_primers(-target => $target);

I<Input:>

    $target - reference to an array holding the start and end positions
              in the B<Bio::Seq> objecct of the target sequence.

I<Output:>

Returns a list of B<Primerpair> objects

I<Description:>

Uses Primer3Factory to find primers for the target region.First a new
B<Primer3Factory> object is created, then all Primer3 parameters are set
to the values relevant for PRC primer design found in the B<Primerfinder>
object (parameters are set through the B<Primer3Factory>), finally
appropriate PCR primers are designed.

=cut

sub pcr_primers {
    my ( $self, %args ) = @_;
    my $p3f = AT::Tools::Run::Primer3->new();
    my $target;
    my @valid_pp;
    if ( defined $args{-target} ) {
        $target = $args{-target};
        delete $args{-target};
        my $pcr_target =
          ( $target->[0] - $self->get_xp() ) . ","
          . ( ( $target->[1] + $self->get_xp() ) -
              ( $target->[0] - $self->get_xp() ) );
        $p3f->set_param( TARGET => $pcr_target );
    }

    if ( defined $args{-exclude_regions})  {
        my $regstring = join " ", map {$_->[0].",".($_->[1]-$_->[0]+1)} @{$args{-exclude_regions}};
        $p3f->set_param(EXCLUDED_REGION => $regstring);
        print STDERR "REGSTRING: '$regstring'\n";
    }

    # Check if target+extra lies within the sequence...
    if ( ( $target->[0] - $self->get_xp() ) < 0 ) {
        carp join( "",
            "Target region starts outside of provided sequence,",
            "impossible to design primer for this target!\n" );
        croak "Target region error!\n";
    }
    elsif ( ( $target->[1] + $self->get_xp() ) >
        ( $self->get_seq()->length() - 1 ) )
    {
        carp join( "",
            "Target region ends outside of provided sequence,",
            " impossible to design primer for this target!\n" );
        croak "Target region error!\n";
    }

    #If target is OK:
    else {

        #setting all Primer3Factory inputs...
        if ( $self->get_p_max_tm() ) {
            $p3f->set_param( PRIMER_MAX_TM => $self->get_p_max_tm() );
        }
        if ( $self->get_p_min_tm() ) {
            $p3f->set_param( PRIMER_MIN_TM => $self->get_p_min_tm() );
        }
        if ( $self->get_p_opt_tm() ) {
            $p3f->set_param( PRIMER_OPT_TM => $self->get_p_opt_tm() );
        }
        if ( $self->get_p_gc() ) {
            $p3f->set_param( PRIMER_GC_CLAMP => $self->get_p_gc() );
        }
        if ( $self->get_p_optl() ) {
            $p3f->set_param( PRIMER_OPT_SIZE => $self->get_p_optl() );
        }
        if ( $self->get_pcrl() ) {
            $p3f->set_param( PRIMER_PRODUCT_SIZE_RANGE => $self->get_pcrl() );
        }

        #Sequence in which primers should be picked (Included region):
        my ( $end, $start, $included_region );
        my $ir = $self->{_ir};

        #Does the start of the included region lie within the sequence?
        if ( ( $target->[0] - $self->get_xp() ) > ( $ir - 1 ) ) {
            $start = $target->[0] - $self->get_xp() - $ir;
        }
        else {
            $start = 0;
        }

        #Does the end of the included region lie within the sequence?
        if ( ( $self->get_seq()->length() - 1 ) >
            ( $target->[1] + $self->get_xp() + $ir ) )
        {
            $end = $target->[1] + $self->get_xp() + $ir;
        }
        else {
            $end = $self->get_seq()->length() - 1;
        }
        $included_region = $start . "," . ( $end - $start );
        $p3f->set_param( INCLUDED_REGION => $included_region );

        #Number of primerpairs returned
        if ( $self->get_pcr_nret() ) {
            $p3f->set_param( PRIMER_NUM_RETURN => $self->get_pcr_nret() );
        }

        #starts searching for primers:
        my $found_primerpair = 0;

        #run primer search:
        my $seq = $self->get_seq()->seq();
        $p3f->run( -sequence => $seq );
        #print STDERR "ERRORS: ".$p3f->resultstone->PRIMER_ERROR."\n";
        #Looping through all possible primers:
        foreach my $pp ( $p3f->get_primer_pairs ) {

            #Analysing the primer for false priming:
            if ( $pp->LEFT_SEQUENCE && $pp->RIGHT_SEQUENCE ) {
                my @poslength_right     = split /,/, $pp->RIGHT;
                my $startpos_right      = $poslength_right[0];
                my @poslength_left      = split /,/, $pp->LEFT;
                my $startpos_left       = $poslength_left[0];
                my $false_priming_right = 0;
                my $false_priming_left  = 0;
                $false_priming_right = search_false_primingsites(
                    -seq      => $self->get_seq->seq(),
                    -prim     => $pp->RIGHT_SEQUENCE,
                    -fp       => $self->get_fp(),
                    -start    => $startpos_right,
                    -is_right => 1
                );
                $false_priming_left = search_false_primingsites(
                    -seq      => $self->get_seq()->seq(),
                    -prim     => $pp->LEFT_SEQUENCE,
                    -fp       => $self->get_fp(),
                    -start    => $startpos_left,
                    -is_right => 0
                );

                if ( $false_priming_right == 0 && $false_priming_left == 0 ) {
                    push @valid_pp, $pp;
                    $found_primerpair = 1;
                }

                #else{
                # print "false priming primer\n";
                #}
            }
        }
        if ( $found_primerpair == 0 ) {
            carp "Could not find valid PCR primerpair for this region\n";
        }
    }
    return (@valid_pp);
}

=head2 B<seq_primers()>

I<Title:> seq_primers()

I<Usage:> $primerfinder -> seq_primers(-target => $target, -pp => $pp);

I<Input:>

    $target - reference to an array holding the start and end positions
              in the B<Bio::Seq> object  of the target sequence.

    $pp - B<Primerpair> object for the PCR reaction, sequencing primers
          should be designed within the genomic region spanned by the
          PCR primers.


I<Output:>

Returns a list of B<Primerpair> objects to be used as sequencing primers.

I<Description:>

Uses Primer3Factory to find primers for the target region.First a new
B<Primer3Factory> object is created, then all Primer3 parameters are set
to the values relevant for sequencing primer design found in the
B<Primerfinder> object (parameters are set through the B<Primer3Factory>),
finally appropriate PCR primers are designed.

=cut

sub seq_primers {
    my ( $self, %args ) = @_;
    my $p3f = AT::Tools::Run::Primer3->new();
    my ( $target, $pp );
    my @valid_pp;
    if ( defined $args{-target} ) {
        $target = $args{-target};
        delete $args{-target};
        my $seq_target =
          ( $target->[0] - $self->get_xs() ) . ","
          . ( ( $target->[1] + $self->get_xs() ) -
              ( $target->[0] - $self->get_xs() ) );
        $p3f->set_param( TARGET => $seq_target );
    }
    else {
        die "No target provided";
    }
    if ( defined $args{-exclude_regions})  {
        my $regstring = join " ", map {$_->[0].",".($_->[1]-$_->[0]+1)} @{$args{-exclude_regions}};
        $p3f->set_param(EXCLUDED_REGION => $regstring);
    }
    if ( defined $args{-pp} ) {
        $pp = $args{-pp};
        delete $args{-pp};
        my @poslength_right = split /,/, $pp->RIGHT;
        my $startpos_right  = $poslength_right[0];
        my @poslength_left  = split /,/, $pp->LEFT;
        my $startpos_left   = $poslength_left[0];

        my $included_region = (
            ( $startpos_left + $self->get_d() ) . ","
              . (
                ( $startpos_right - $self->get_d() ) -
                  ( $startpos_left + $self->get_d() )
              )
        );
        $p3f->set_param( INCLUDED_REGION => $included_region );
    }
    else {
        die "No pcr primer pair provided";
    }

    #setting all Primer3Factory inputs...
    if ( $self->get_s_max_tm() ) {
        $p3f->set_param( PRIMER_MAX_TM => $self->get_s_max_tm() );
    }
    if ( $self->get_s_min_tm() ) {
        $p3f->set_param( PRIMER_MIN_TM => $self->get_s_min_tm() );
    }
    if ( $self->get_s_opt_tm() ) {
        $p3f->set_param( PRIMER_OPT_TM => $self->get_s_opt_tm() );
    }
    if ( $self->get_s_gc() ) {
        $p3f->set_param( PRIMER_GC_CLAMP => $self->get_s_gc() );
    }
    if ( $self->get_s_optl() ) {
        $p3f->set_param( PRIMER_OPT_SIZE => $self->get_s_optl() );
    }
    if ( $self->get_sl() ) {
        $p3f->set_param( PRIMER_PRODUCT_SIZE_RANGE => $self->get_sl() );
    }
    if ( $self->get_s_nret() ) {
        $p3f->set_param( PRIMER_NUM_RETURN => $self->get_s_nret() );
    }

    #run primer search:
    my $seq = $self->get_seq()->seq();
    $p3f->run( -sequence => $seq );
    if ( $p3f->get_primer_pairs ) {
        @valid_pp = $p3f->get_primer_pairs;
        return (@valid_pp);
    }
    else {
        carp "Unable to find valid sequencing primer pairs for this region\n";
        return (@valid_pp);
    }
}

=head2 B<search_false_primingsites()>

I<Title:> search_false_primingsites()

I<Usage:> Used internally by B<Primerfinder::pcr_primers()>.

I<Output:> Returns "1" if the primer is false priming, "0" otherwise.

I<Description:> Checks if a PCR primer falsely matches another part of the
sequence.

=cut

sub search_false_primingsites {

    # fp = accepted number of nucleotides identical in the
    #      3'end of the primer in a false priming site.
    # $false_priming indicates if the primer's 3'end is identical to any
    #      part of the sequence besides the primer (Which could lead to
    #      false priming...).

    my $is_false_priming = 0;
    my %args             = @_;
    my $seq              = $args{-seq};
    my $primer           = $args{-prim};
    my $start            = $args{-start};
    my $fp               = $args{-fp};
    my $is_right         = $args{-is_right};

    if ($is_right) {
        my $primerObj = Bio::LocatableSeq->new(
            -seq   => $primer,
            -id    => $primer,
            -strand => 0, # to prevent warnings
            -start => 1,
            -end   => ( length($primer) - 1 )
        );
        $primerObj = $primerObj->revcom;
        $primer    = $primerObj->seq;
        foreach my $i ( 0 .. ( length($seq) - length($primer) ) ) {
            if (   ( substr $primer, 0, $fp ) eq ( substr $seq, $i, $fp )
                && ( $i + length($primer) - 1 ) != $start )
            {
                $is_false_priming = 1;
            }
        }
    }
    else {
        foreach my $i ( 0 .. ( length($seq) - length($primer) ) ) {
            if ( ( substr $primer, ( length($primer) - $fp ), $fp ) eq
                ( substr $seq, ( $i + length($primer) - $fp ), $fp )
                && $i != $start )
            {
                $is_false_priming = 1;
            }
        }
    }
    return $is_false_priming;
}

=head2 B<calculate_tm()>

I<Title:> calculate_tm()

I<Usage:> Used internally by B<Primerfinder::pcr_primers()>.

I<Output:> A real number.

I<Description:> Calculates the Tm value of a primer using the old
2(A+T)+4(C+G) method. (Perhaps good if you are used to calculate it
like this...)

=cut

sub calculate_tm {
    my %args   = @_;
    my $primer = $args{-primer};
    my $Tm     = 0;
    for my $pos ( 0 .. length($primer) ) {
        my $base_at_pos = substr $primer, $pos, 1;
        if (   $base_at_pos eq "A"
            || $base_at_pos eq "T"
            || $base_at_pos eq "a"
            || $base_at_pos eq "t" )
        {
            $Tm = $Tm + 2;
        }
        if (   $base_at_pos eq "G"
            || $base_at_pos eq "C"
            || $base_at_pos eq "g"
            || $base_at_pos eq "c" )
        {
            $Tm = $Tm + 4;
        }
    }
    return $Tm;
}

1;
