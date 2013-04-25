package Bio::Das::ProServer::SourceAdaptor::Transport::MyCNEDatabase;

use strict;
use CNE::DB;
use base qw(Bio::Das::ProServer::SourceAdaptor::Transport::generic);

sub adaptor {
  my $self = shift;
  unless($self->{'_adaptor'}) {
    my $host     = $self->config->{'dbhost'} || "localhost";
    my $port     = $self->config->{'dbport'} || "3306";
    my $dbname   = $self->config->{'dbname'} || "cne";
    my $username = $self->config->{'dbuser'} || "nobody";
    my $password = $self->config->{'dbpass'};

    $self->{'_adaptor'} ||= CNE::DB->connect(-dbhost => $host,
					     -dbport => $port,
					     -dbuser => $username,
					     -dbname => $dbname,
					     -dbpass => $password);
    
    # cache connection?
  }

  return $self->{'_adaptor'};
}

sub disconnect {
    my $self = shift;
    undef $self->{'_adaptor'};
}

1;
