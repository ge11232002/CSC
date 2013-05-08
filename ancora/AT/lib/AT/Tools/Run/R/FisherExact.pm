#
# Copyright Par Engstrom and Boris Lenhard
#
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::Tools::Run::R::FisherExact - wrapper for the R method fisher.test

=head1 SYNOPSIS

 use AT::Tools::Run::R::FisherExact;

 my $data = [[2,5,2,1], [160,740,40,60]];
 my $p_values = AT::Tools::Run::R::FisherExact->new->test($data);
 foreach my $p (@$p_values) { print "$p\n"; }

=head1 APPENDIX

Runs R to perform a set of Fisher Exact tests on 2x2 contingency tables.
The p-values for the tests are returned.
It is advisable to use the module for batch runs since starting R takes a couple
of seconds. Once the program is started, the tests are performed very quickly.

The input data for this module's test() method should be one or more
2x2 contingency tables represented as 4-element arrays.
The elements are (in order):
 successes in first sample
 successes in second sample
 failures in first sample
 failures in second sample

If comparing a sample to the population from which it is drawn, the contingency
table can be calculated as:
 x
 pN-x
 n-x
 N-pN+x-n

Where n is the sample size, N the population size, x the number of successes in
the sample and p the proportion of successes in the population.

By default, R calculates p-values for 2-sided tests (over- and underrepresentation).
Use the options method to change this.

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::Tools::Run::R::FisherExact;
use strict;
use vars '@ISA';
use Carp;
use AT::Root;
use File::Temp qw(tempfile);
use IPC::Open2;
use IO::Handle;

@ISA = qw(AT::Root);

=head2 new

 Title     : new
 Usage     : my $fisher = AT::Tools::Run::R::FisherExact->new;
 Function  : Constructor
 Returns   : FisherExact object
 Args      : options	Optional. Array of options to be passed to
			the fisher.test method in R.
			E.g.: options => ["alternative='l'", "or=1.5"]
			Do help(fisher.test) in R for details.

=cut

sub new  {
    my ($caller, %args) = @_;
    my $self = bless {
	options => ($args{'options'} || [])
    }, ref $caller || $caller;
    return $self;
}


=head2 test

 Title     : test
 Usage     : my $p_value_ary = $fisher->test(
		[ [2,5,2,1],
 	          [160,740,40,60]
	        ]);
	     my $p_value_hash = $fisher->test(
		{ a_test => [2,5,2,1],
		  another_test => [160,740,40,60]
		});
 Function  : Do the test!
 Returns   : If the argument was an array ref:
	       Ref to array of p-values for the tests. The first p-value
	       corresponds to the first contingency table in the input
	       data and so on.
	     If the argument was a hash ref:
	       Ref to hash where the keys those in the input hash
	       and each value is a p-value.
 Args      : Ref to an array of arrays or to a hash of arrays.
	     Each subarray should cotain the elements of a 2x2
	     contingency table to be tested.
	     See module description above.

=cut

sub test {
    my ($self, $data) = @_;
    my $result;
    
    if(ref($data) eq 'ARRAY') {
	$result = $self->_do_test($data);
    }
    elsif(ref($data) eq 'HASH') {
	my @keys = keys %$data;
	my @data_ary;
	foreach my $key (@keys) {
	    push @data_ary, $data->{$key};
	}
	my $result_aryref = $self->_do_test(\@data_ary);
	my %result_hash;
    	for my $i (0..@keys-1) {
	    $result_hash{$keys[$i]} = $result_aryref->[$i];
	}
	$result = \%result_hash;
    }
    else {
	croak "Argument should be reference to array or hash";
    }

    return $result;
}


sub _do_test {
    my ($self, $data) = @_;
    
    # Write data to file that R can input
    my $data_fn = $self->_data_to_tmpfile($data);

    # Start R
    my($R_out, $R_in) = (IO::Handle->new, IO::Handle->new);
    my $R_pid = open2($R_out, $R_in, "R --slave --vanilla");

    # Run R commands
    print $R_in
	"fish <- function(data) fisher.test(",
	join(', ', ('matrix(data,nrow=2)', @{$self->options})),
	")\$p.value;\n";
    print $R_in "A <- read.table('$data_fn');\n";
    print $R_in "B <- sapply(data.frame(t(A)),fish);\n";
    print $R_in "write.table(B,row.names=FALSE,col.names=FALSE,quote=FALSE);\n";
    $R_in->close();

    # Read output (p-values) from R
    my @p_values;
    while(my $p = <$R_out>) {
	chomp $p;
	push @p_values, $p;
    }
    $R_out->close();
    waitpid($R_pid, 0);
    pop @p_values if($p_values[-1] eq '');
    warn("Unexpected number of p-values returned from R")
	unless(scalar(@p_values) eq scalar(@$data));
    
    return \@p_values;
}


sub _data_to_tmpfile {
    my ($self, $data) = @_;
    my ($tmp_fh, $tmp_fn) = tempfile(UNLINK => 1);
    foreach my $row (@$data) {
	print $tmp_fh join("\t", @$row), "\n";
    }
    close($tmp_fh);
    return $tmp_fn;
}


# Docs for AUTOLOADed methods follow...

=head2 options

 Title     : options
 Usage     : $fisher->options(["alternative='l'", "or=1.5"];
 Function  : Get/set options array
 Returns   : Array ref.
 Args      : Array of options to be passed to the fisher.test method in R.

=cut


1;

