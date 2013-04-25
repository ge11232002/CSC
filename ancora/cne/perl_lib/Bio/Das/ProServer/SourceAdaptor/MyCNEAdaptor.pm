package Bio::Das::ProServer::SourceAdaptor::MyCNEAdaptor;

use strict;
use warnings;
use base qw(Bio::Das::ProServer::SourceAdaptor::MyCNEAdaptorBase);

use constant DEFAULT_TYPE => 'HCNE';

sub init {
  my $self = shift;
  $self->{'capabilities'} = {'features'      => '1.0',
			     'stylesheet'    => '1.0'};

  my $dsn = $self->dsn;
  my $config = $self->config;
  my $transport = $self->transport;
  my $db = $transport->adaptor;

  # Parse dsn to get assembly ids and cutoffs
  my ($asm1_id, $asm2_id, $len, $identity) = $dsn =~ /^([A-z0-9]+)_cne_([A-z0-9]+)_len(\d+)_id(\d+)$/;
  defined($asm1_id) or die "Unable to parse dsn $dsn";
  $identity = $identity / 10; # Change to percent

  # Get more info about the two assemblies
  my $asm1_info = $db->get_assembly_info($asm1_id) or die "Unknown assembly $asm1_id";
  my $asm2_info = $db->get_assembly_info($asm2_id) or die "Unknown assembly $asm2_id";

  # Set information about DAS source unless specified in config file
  unless($config->{'title'}) {
      $self->{'title'} = $self->_make_short_organism_name($asm2_info).' HCNEs';
  }
  unless($config->{'description'}) {
      $self->{'description'} =
	  "Non-exonic sequences conserved in the ".$self->_make_long_organism_name($asm2_info).
	  " genome. Minimum $identity% identity over $len alignment columns.";
  }
  $self->_set_coordinates($asm1_info);
  $self->_set_mapmaster($asm1_id);

  # Set private parameters to be used when fetching features
  $self->{'cne_method_label'} = "Scan for non-exonic sequences with at least $identity% identity over $len alignment columns";
  $self->{'cne_db_table'} = $config->{'dbtable'} || $db->get_cne_table_name(assembly1 => $asm1_id, assembly2 => $asm2_id,
									    min_identity => $identity/100, min_length => $len);
  $self->{'cne_assembly'} = $asm1_id;

  $transport->disconnect(); # must do this because ProServer does not clean up after dsn and sources requests
}


sub das_stylesheet {
    my $self = shift;
    my $config = $self->config;

    if($config->{'stylesheet'} or $config->{'stylesheetfile'}) {
	return $self->SUPER::das_stylesheet();
    }

    my $color = $config->{'color'} || 'black';
    return qq(<?xml version="1.0" standalone="yes"?>
<!DOCTYPE DASSTYLE SYSTEM "http://www.biodas.org/dtd/dasstyle.dtd">
<DASSTYLE>
<STYLESHEET version="1.0">
  <CATEGORY id="default">
    <TYPE id="default">
        <GLYPH>
           <BOX>
	      <HEIGHT>6</HEIGHT>
              <BGCOLOR>$color</BGCOLOR>
              <FGCOLOR>$color</FGCOLOR>
	      <BUMP>yes</BUMP>
           </BOX>
        </GLYPH>
     </TYPE>
  </CATEGORY>
</STYLESHEET>
</DASSTYLE>
);
}


# The entry_points and types calls do not seem to be used by Ensembl.
# I have therefore disabled them.
# The two corresponding methods build_entry_points and build_types below
# are just scaffolds in case we want to implement this functionality.

sub build_entry_points
{
    my ($self) = @_;

    my @ep = ({ segment => '1', length => 247249719, subparts => 'no'},
	      { segment => 'chr1', length => 247249719, subparts => 'no'});
    return @ep;
}


sub build_types
{
    my ($self) = @_;

    my $dsn = $self->{'dsn'};
    my $type = $self->config->{'type'} || DEFAULT_TYPE;
    my $table = $self->config->{'dbtable'} || $dsn;
    my $method = $self->config->{'method'} || $table;
    my $type_text = $self->config->{'type_text'} || 'Highly conserved noncoding element';

    my @types = ( { type => $type,
		    method => $method,
		    category => 'similarity',
		    description => $type_text
		  }
	);
    return @types;
}

# Build features - this we support

sub build_features {
  my ($self, $opts) = @_;

  # Get dsn
  my $dsn = $self->dsn;

  # Get private parameters set by init()
  my $table = $self->{'cne_db_table'};
  my $asm = $self->{'cne_assembly'};
  my $method_label = $self->{'cne_method_label'};

  # Get parameters optionally set in config file
  my $config = $self->config;
  my $type = $config->{'type'} || DEFAULT_TYPE;
  my $type_text = $config->{'type_text'} || 'Highly conserved noncoding element';
  my $method = $config->{'method'} || $table;
  my $details_url = defined($config->{'details_script'}) ? $config->{'details_script'}."?type=cne&table=$table&asm=$asm&id=" : '';
  
  # Get location to query
  my $chr = $opts->{'segment'} || die "no segment";
  my $start = $opts->{'start'} || 1;
  my $end = $opts->{'end'} || 5e8;
  my $my_chr = $chr =~ /^chr/ ? $chr : "chr$chr";

  # Run query to obtain CNE objects
  my $transport = $self->transport;
  my $cnes = $transport->adaptor->get_cnes_in_region(chr => $my_chr,
							   start => $start,
							   end => $end,
							   table_name => $table,
							   assembly => $asm);

  # Make features of the expected format from the CNEs
  my @features;
  for my $cne (@$cnes) {
    push @features, {
		     id     => $cne->id,
                     label  => $cne->chr2.':'.commas($cne->start2),
                     type   => $type,
		     typetxt => $type_text,
		     typecategory => 'similarity',
                     method => $method,
		     method_label => $method_label,
                     start  => $cne->start1,
                     end    => $cne->end1,
		     target_id => $cne->chr2,
		     target_start => $cne->start2,
		     target_stop => $cne->end2,
		     link   => $details_url.$cne->id,
		     linktxt => "Details about this HCNE",
		     #ori    => '+',
		     #note   => "some_note",
                    };
  }

  $transport->disconnect(); # this may not be necessary but just in case...

  # Return the features
  return @features;
}


sub commas
{
    my ($pos) = @_;
    my @triplets;
    while(length($pos)) {
	unshift @triplets, substr($pos, -3, 3, '');
    }
    return join(',',@triplets);
}


1;
