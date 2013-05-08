package Bio::Graphics::Browser::Plugin::CNEPlotInstance4;
use strict;
use vars '@ISA';
#use lib '/home/engstrom/DEVEL/CNE/lib';
use Bio::Graphics::Browser::Plugin::CNEPlot;
@ISA = qw(Bio::Graphics::Browser::Plugin::CNEPlot);


sub static_plugin_setting {
    my ($self, $name) = @_;
    return $self->browser_config->plugin_setting($name);
}
