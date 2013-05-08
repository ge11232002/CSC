############################################################################
# Masker.pm - perl module for handeling a set of composite alignments
#
############################################################################

package AT::Tools::Run::Masker;

use strict;
use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long;
use File::Temp qw/ tempfile tempdir /;


=head1 NAME

Masker - module for masking sequences for known repeats

=head1 SYNOPSIS

    use Masker;
    my $masker =AT::Tools::Run:: Masker -> new(-seq         => Bio::Seq object,
			       );

    $masker -> get_in_seq();
    $masker -> get_masked_seq();
    $masker -> run();


=head1 DESCRIPTION

This method maskes the input sequence and returnes a reference to the masked sequenceobject.

=head1 METHODS DESCRIPTION

=cut

############################################################################


=head2 B<new()>

I<Title:> new()

I<Usage:> 

$masker = Masker -> new(-seq => Bio::Seq object);

I<Input:> $seq  - The B<Bio::Seq> object that should be masked.
	 
I<Output:> Returns a B<Masker> object. 

I<Description:>

=cut

sub new {
  my($class, %args) = @_;
  my $seq;
  if(defined $args{-seq}){
    $seq = $args{-seq};
  }
  return bless { _seq		=> $seq,
      		 _tempdir	=> tempdir("MaskerXXXXXXXX", TMPDIR => 1,
					    CLEANUP => 1)
	       }, ref $class || $class;
}

=head2 B<get_seq()>

I<Title:> get_seq()

I<Usage:> $masker -> get_seq();

I<Input:> No input

I<Output:> Returns a reference to the input sequence object

I<Description:> 

=cut

sub get_masked_seq        {$_[0] -> {_maskedseq} }

=head2 B<get_masked_seq()>

I<Title:> get_masked_seq()

I<Usage:> $masker -> get_masked_seq();

I<Input:> No input

I<Output:> Returns a refernce to the masked sequence object

I<Description:> 

=cut

sub get_seq        {$_[0] -> {_seq} }
=head2 B<get_tempfile()>

I<Title:> get_tempfile()

I<Usage:> $masker -> get_tempfile();

I<Input:> No input

I<Output:> Returns the temporary filename for the sequence object to be masked.

I<Description:> 

=cut

sub get_tempfile     {$_[0] -> {_tempfile} }

=head2 B<run()>

I<Title:> run()

I<Usage:> $masker -> run();

I<Input:> No input

I<Output:> Returns a reference to the masked B<Bio::Seq> object.

=cut 

sub run {
  my($self, %args)=@_;
  $self->print_seq_to_file(-seq => $self->get_seq());
  my $tempfile=$self->get_tempfile();
  print "Tempfile: ".$tempfile."\n";
  system ("repeatmasker $tempfile");
  my $maskedseqIO  = Bio::SeqIO -> new(-file => $tempfile.".masked" , '-format' => 'fasta');
  my $maskedseqObj = $maskedseqIO -> next_seq();
  $self->{_maskedseq} = $maskedseqObj;
  $self -> delete_tempfiles();
  return $maskedseqObj;
}

=head2 B<print_seq_to_file()>

I<Title:> print_seq_to_file()

I<Usage:> used internally by B<Run::Masker::run()>

I<Output:> Nothing returned

I<Description:> Prints the B<Bio::Seq> object to be masked to a temporary file that can be read by repeatmasker in B<Run::Blastz::run>, and thereafter deleted.

=cut

sub print_seq_to_file{
  my ($self, %args)=@_;
  my $tempfile       = tempnam($self->{_tempdir}, "seq");
  my $output_stream  = Bio::SeqIO->new(-file => ">".$tempfile, -format=>'Fasta');
  $output_stream    ->write_seq( $self->get_seq());
  $self->{_tempfile} = $tempfile;
}


=head2 B<delete_tempfiles()>

I<Description> Removes the temporary file holding the result from Blastz.

=cut


sub delete_tempfiles {
  my ($self,%args) = @_;
  my $tempfile= $self->get_tempfile();
  my $all_temp_files=$tempfile.".*";
  print "\nAll temp files: ".$all_temp_files."\n";
  system ("rm -f $tempfile");
  system ("rm -f $all_temp_files");
}


1;
