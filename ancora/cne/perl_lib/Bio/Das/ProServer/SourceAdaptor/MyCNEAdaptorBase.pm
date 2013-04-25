package Bio::Das::ProServer::SourceAdaptor::MyCNEAdaptorBase;

use strict;
use warnings;
use Bio::Das::ProServer::SourceAdaptor::Transport::MyCNEDatabase;
use base qw(Bio::Das::ProServer::SourceAdaptor);

use constant MAPMASTER_PREFIX => 'http://genome.ucsc.edu/cgi-bin/das/';

sub _set_coordinates {
    my ($self, $asm_info) = @_;
    unless($self->config->{'coordinates'}) {
	my $coord_system = $asm_info->{'das_coord'} or
	    die "No coordinate system defined for assembly ".$asm_info->assembly_id;
	my $default_loc = $asm_info->{'default_ensembl_loc'} or 
	    die "No default Ensembl location defined for assembly ".$asm_info->assembly_id;
	$default_loc =~ s/-/,/;
	$self->{'coordinates'} = { $coord_system => $default_loc };
    }
}


sub _set_mapmaster {
    my ($self, $asm_id) = @_;
    unless($self->config->{'mapmaster'}) {
	$self->{'mapmaster'} =  MAPMASTER_PREFIX.$asm_id;
    }
}

sub _make_short_organism_name {
    my ($self, $asm_info) = @_;
    $asm_info->{organism_latin} =~ /^(.)\S+\s+(...)/;
    return $1.$2;
}

sub _make_long_organism_name {
    my ($self, $asm_info) = @_;
    if($asm_info->{organism_common}) {
	return $asm_info->{organism_common}.' ('.$asm_info->{organism_latin}.')';
    }
    else {
	return $asm_info->{organism_latin};
    }
}

1;
