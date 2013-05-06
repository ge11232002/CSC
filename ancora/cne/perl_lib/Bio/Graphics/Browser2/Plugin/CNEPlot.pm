package Bio::Graphics::Browser2::Plugin::CNEPlot;
use strict;
use Bio::Graphics::Browser2::Plugin;
use Bio::Graphics::Feature;
use Benchmark;
use POSIX qw(strtod);
use CGI qw(:standard *table);
use vars '$VERSION','@ISA';
use CNE::DB;
use AT::GFX::SeqColorIndex;
use AT::Tools::RangeHandler;
use Data::Dumper;

$VERSION = '0.1'; # This supposed to be the version number for the library, I think...

@ISA = qw(Bio::Graphics::Browser2::Plugin);

# Package-specific version number. Used to check whether a users settings for this plugin are
# from the same version of the plugin. If not, plugin settings will be set to the default for
# this version. This is to avoid strange effects when we add new functionality.
my $MY_VERSION = 0.91;

my $ABS_SCORE = 0;
my $DB_NAME = "cne";
my $DB_HOST = "localhost";
my $DB_USER = "nobody";
my $METHODS_URL = "/methods.html";

my $WIN_SIZE_MIN = 10;
my $WIN_SIZE_MAX = 1000;
my $Y_MIN = 0.1;
my $Y_MAX = 100;
my $NR_GRAPHS_MAX = 5;
my $DEFAULT_NR_CNE_SETS = 3;
my $DEFAULT_MIN_CNE_LENGTH = 50;
my $DEFAULT_CNE_DISPLAY_THRESHOLD = 1e7;
my $DEFAULT_MAX_INTENSITY = 140;
my %DEFAULT_COLORS = ( 3 => [[200,200,0], [220,130,0], [255,0,0]],
		       4 => [[200,200,0], [200,150,0], [220,100,0], [255,0,0]]);

# called by gbrowse to return name of plugin for popup menu
sub name { 
    shift->asm2_name . " HCNEs";
}
# called by gbrowse to return description of plugin
sub verb { 'Configure' }

# called by gbrowse to return description of plugin
sub description {
    my $self = shift;
    p("This track shows ".$self->name().
      ". To configure this track, select it from the Reports & Analysis menu or click on its density plot.");
}

# called by gbrowse to return type of plugin
sub type { 'annotator' }

sub init {
    my $self = shift;
    warn "I am in the init() CNE plogin now!!!";
    $self->{_nr_cne_sets} = $self->static_plugin_setting('nr_cne_sets') || $DEFAULT_NR_CNE_SETS;
    #$self->{_asm1} = $self->browser_config->source();
    $self->{_asm1} =  $self->browser_config->name();
    $self->{_asm2} = $self->static_plugin_setting('asm_id') || die "No second assembly configured for plugin";
    $self->{_asm2_name} = $self->static_plugin_setting('asm_name') || $self->{_asm2};
}

# called by gbrowse to configure default settings for plugin
sub config_defaults {
  my $self = shift;

  warn "I am in CNE config_defaults now!!!";
  my $color_cnes = $self->static_plugin_setting('color_cnes');
  $color_cnes = 'dark' unless defined $color_cnes;

  my $bump_cnes = $self->static_plugin_setting('bump_cnes');
  $bump_cnes = 1 unless defined $bump_cnes;

  my %def = (window_size => $self->static_plugin_setting('window_size') || 50,
	     max_score => $self->static_plugin_setting('max_score') || 20,
	     nr_graphs => $self->static_plugin_setting('nr_graphs') || 1,
	     rank_un_low => $self->static_plugin_setting('rank_un_low') || 0,
	     color_cnes => $color_cnes,
	     bump_cnes => $bump_cnes
      );

  for my $i (1..$self->nr_cne_sets) {
      $def{"min_cne_id$i"} = $self->static_plugin_setting("min_cne_id$i") || '';
      $def{"min_cne_len$i"} = $self->static_plugin_setting("min_cne_len$i") || $DEFAULT_MIN_CNE_LENGTH;
      $def{"show_cnes$i"} = $self->static_plugin_setting("show_cnes$i") || 0;
      $def{"show_graph$i"} = $self->static_plugin_setting("show_graph$i") || 0;
  }

  $def{'version'} = $MY_VERSION;

  return \%def;
}


sub nr_cne_sets { shift->{_nr_cne_sets} }
sub asm1 { shift->{_asm1} }
sub asm2 { shift->{_asm2} }
sub asm2_name { shift->{_asm2_name} }

sub _set_config_to_default {
    my $self = shift;
    $self->configuration($self->config_defaults);
}

sub _get_colors
{
    my ($self, $set_nrs) = @_;

    my $color_string = $self->static_plugin_setting('colors');
    my $nr_cne_sets = $self->nr_cne_sets;
    my @color_list;
    if($color_string) {
	@color_list =  split /\s+/, $color_string;
	die "Configuration error for plugin ".$self->name.": number of intensities does not match number of CNE sets."
	    unless(@color_list == $nr_cne_sets);
	@color_list = @color_list[map { $_ - 1} @$set_nrs] if($set_nrs);
	for my $i (0..@color_list-1) {
	    my ($r,$g,$b) = split ',', $color_list[$i];
	    die "error parsing color string $color_list[$i]; should have format red,green,blue; e.g. 0,0,255" unless(defined $b);
	    $color_list[$i] = [$r,$g,$b];
	}
    }
    elsif(my $def_colors = $DEFAULT_COLORS{$nr_cne_sets}) {
	@color_list = $set_nrs ? map { $def_colors->[$_ - 1] } @$set_nrs : @$def_colors;
    }
    else {
	my $max_int = $DEFAULT_MAX_INTENSITY;
	my $int_step = $nr_cne_sets > 1 ? int($max_int / ($nr_cne_sets-1)) : 0;
	@color_list = map { [($max_int - ($_-1) * $int_step) x 3] } ($set_nrs ? @$set_nrs : (1..$nr_cne_sets));
    }
    return \@color_list;
}


sub _parse_integer_input
{
    my ($self, $n, $min, $max, $default) = @_;
    return $default unless(defined($n) and $n =~ /\d/);
    $n = int($n);
    return $min if($n < $min);
    return $max if($n > $max);
    return $n;
}


sub _parse_decimal_input
{
    my ($self, $n, $min, $max, $default) = @_;
    return $default unless(defined $n);
    $n =~ s/^\s+//;
    $n =~ s/\s+$//;
    return $default if($n eq '');
    my $unparsed;
    $! = 0;
    ($n, $unparsed) = strtod($n);
    return $default if($unparsed != 0 || $!);
    return $min if($n < $min);
    return $max if($n > $max);
    return $n;
} 

# called by gbrowse to reconfigure plugin settings based on CGI parameters
sub reconfigure {
    my $self = shift;
    my $config = $self->configuration;
    my $defaults = $self->config_defaults;

    # Update HCNE selection
    for my $i (1..$self->nr_cne_sets) {
	$config->{"min_cne_len$i"} = $self->_parse_integer_input($self->config_param("min_cne_len$i"),
								 1, 999, $defaults->{"min_cne_len$i"});
	foreach my $name ("min_cne_id$i", "show_graph$i", "show_cnes$i") {
	    $config->{$name} = $self->config_param($name) || 0;
	}
	#$config->{"min_cne_id$i"} = $self->config_param("min_cne_id$i") || defaults->{"min_cne_id$i"};
	#$config->{"show_graph$i"} = $self->config_param("show_graph$i") || 0;
	#$config->{"show_cnes$i"} = $self->config_param("show_cnes$i") || 0;
    }
    
    # Update window size
    $config->{'window_size'} = $self->_parse_integer_input($self->config_param('window_size'),
							   $WIN_SIZE_MIN, $WIN_SIZE_MAX, 
							   $defaults->{'window_size'}); 
    
    # Update y-axis height
    $config->{'max_score'} = $self->_parse_decimal_input($self->config_param('max_score'),
							 $Y_MIN, $Y_MAX, 
							 $defaults->{'max_score'}); 
    
    # Update number of graphs to show
    $config->{'nr_graphs'} = $self->_parse_integer_input($self->config_param('nr_graphs'),
							 1, $NR_GRAPHS_MAX,
							 $defaults->{'nr_graphs'});
    $config->{'rank_un_low'} = $self->config_param('rank_un_low') || 0;
    
    # Update settings for individual CNEs
    $config->{color_cnes} = $self->config_param('color_cnes') || 0;
    $config->{bump_cnes} = $self->config_param('bump_cnes') ? 1 : 0;

    # Update version
    $config->{version} = $defaults->{version};
}

# called by gbrowse to create a <form> fragment for changing settings
sub configure_form {
    my $self = shift;

    $self->_set_config_to_default() unless(($self->configuration->{version} || 0) eq $MY_VERSION);
    my $current_config = $self->configuration;

    # Connect to DB
    my $db = CNE::DB->connect(dbhost => $DB_HOST,
			      dbname => $DB_NAME,
			      dbuser => $DB_USER)
	or die "could not connect to db $DB_NAME @ $DB_HOST";
    ## debug
    warn "I am in the CNE plugin:configure_form now!!!";
    $Data::Dumper::Indent = 1;
    $Data::Dumper::Sortkeys = 1;
    open (FH, '>/mnt/biggley/home/gtan/debug/debug_configure_form.txt');
    print FH Dumper($db);
    close FH;


    my $nr_cne_sets = $self->nr_cne_sets();
    my $asm1 = $self->asm1();
    my $asm2 = $self->asm2();
    my $cne_display_threshold = $self->static_plugin_setting('cne_display_threshold') || $DEFAULT_CNE_DISPLAY_THRESHOLD;

    my $popup_script =
	"if (! window.focus)return true;
         newWindow = window.open(this.href, 'cne_color_legend', 'width=400,height=240,resizable=yes,scrollbars=yes');
         newWindow.focus();
         return false;";
   
    # Get data for HCNE selection form
    my $table_names = $db->get_cne_table_names_for_assemblies($asm1, $asm2);
    my %table_infos = map { $_ => $db->get_cne_table_info($_) } @$table_names;
    $table_names =  
	[ sort { $table_infos{$a}{min_length} <=> $table_infos{$b}{min_length}
		 or $table_infos{$a}{min_identity} <=> $table_infos{$b}{min_identity} } @$table_names ];
    my (@id_cutoff_codes, %id_cutoff_labels);
    my $seen_100pc;
    foreach my $t (@$table_names) {
	my $table_info = $table_infos{$t};
	next unless($table_info->{version} == 1); # reserve tables with version nr != 1 for non-ancora use
	if($table_info->{min_identity} == 1) {    # if there are several perfect match sets, show that with smallest window size only
	    next if($seen_100pc);
	    $seen_100pc = 1;
	}
	my $id_cutoff_code = $table_info->{min_identity}.'/'.$table_info->{min_length};
	push @id_cutoff_codes, $id_cutoff_code;
	$id_cutoff_labels{$id_cutoff_code} = (int($table_info->{min_identity}*100)).
	    '% over '.$table_info->{min_length}.' columns';
    }
    my $density_colors = $self->_get_colors();

    # Create HTML elements for HCNE selection
    my @cne_selectors;
    my @row_color_strings = ('rgb(200,200,200)', 'rgb(220,220,220)');
    my $j = 0;
    for my $i (1..$nr_cne_sets) {
	my $density_color_str = 'rgb('.join(',',@{shift @$density_colors}).')';
	push @cne_selectors,
	Tr({-align => 'center', -style => 'background-color:'.$row_color_strings[$j] },
	   td($i).
	   td(popup_menu(-name => $self->config_name("min_cne_id$i"),
			 -values => \@id_cutoff_codes,
			 -labels => \%id_cutoff_labels,
			 -default => $current_config->{"min_cne_id$i"})).
	   td(textfield(-name => $self->config_name("min_cne_len$i"),
			-default => $current_config->{"min_cne_len$i"} || $DEFAULT_MIN_CNE_LENGTH,
			-size => 3, -maxlength => 3,).' bp').
	   td('&nbsp;&nbsp;').
	   td({-style => 'background-color:'.$density_color_str},
	      checkbox(-name => $self->config_name("show_graph$i"),
		       -checked => $current_config->{"show_graph$i"},
		       -label => '')).
	   td('&nbsp;&nbsp;').
	   td(checkbox(-name => $self->config_name("show_cnes$i"),
		       -checked => $current_config->{"show_cnes$i"},
		       -label => ''))
	    );
	$j = 1-$j;
    }
    
    # Create HTML string for form
    my %table_head_attr = (-align => 'center', -style => 'padding-left:10px; padding-right:10px');
    my $form = 
	hr().
	h3("HCNE selection")."\n".
	p("Up to $nr_cne_sets HCNE sets from this pairwise comparison may be shown simultaneously.")."\n".
	start_table({-style => 'width:auto', -cellspacing => '0px'}).
	Tr(th(\%table_head_attr, ["Set", "Minimum identity", "Minimum size"]).
	   th({%table_head_attr, -colspan => 3}, "Show density plot").
	   th(\%table_head_attr, "Show HCNEs"))."\n".
	join("\n",  @cne_selectors).
	end_table()."\n".
	p('The procedure that scans for HCNEs measures the fraction of base identities over a fixed number of alignment columns. See '.
	  a({-href => $METHODS_URL}, 'Methods').' for details.')."\n".
	hr().
	h3("Density plot display").
	b("Size of sliding window: ").
	textfield(-name => $self->config_name('window_size'),
		  -default => $current_config->{window_size},
		  -size => 4,
		  -maxlength => 4)
	." kb (range: $WIN_SIZE_MIN to $WIN_SIZE_MAX)\n".br().
	b("Height of y-axis: ").
	textfield(-name => $self->config_name('max_score'),
		  -default => $current_config->{max_score},
	          -size => 4, -maxlength => 4)
	." percent (range: $Y_MIN to $Y_MAX)\n".br().
	b("Separate densities based on chromosome in other genome: "). # No / into 1 graphs, into 2 graphs...
	popup_menu(-name => $self->config_name('nr_graphs'),
		   -values => [1..$NR_GRAPHS_MAX],
		   -labels => { 1 => "no", map { $_ => "into at most $_ plots" } (2..$NR_GRAPHS_MAX) },
		   -default => $current_config->{nr_graphs})
	."\n".br().
	p({style => 'margin-left:+30px; margin-top:5px'},
	  "If the option to separate densities is chosen, sequences from the other genome will be selected based on amount of HCNE sequence in the displayed region and half a sliding window to either side.".br().
	  checkbox(-name => $self->config_name('rank_un_low'),
		   -checked => $current_config->{rank_un_low},
		   -label => 'When selecting sequences, give lowest priority to alternative haplotypes and poorly assembled sequences (e.g. chrUn, chrN_random).'))."\n".
		   hr().
	h3("HCNE display").
	b("Color: ").
	popup_menu(-name => $self->config_name('color_cnes'),
		   -values => [0, 'dark', 'UCSC'],
		   -default => $current_config->{color_cnes},
		   -labels => { 0 => "all black",
				'dark' => "by chromosome in other genome (default colors)",
				'UCSC' => "by chromosome in other genome (UCSC colors)" })
	."\n".
	a({-href=>"../../cne_color_legend", -onClick=>$popup_script}, "Color legend").br()."\n".
	b("Layout: ").
	radio_group(-name => $self->config_name('bump_cnes'),
		    -values => [0,1],
		    -default => $current_config->{bump_cnes},
		    -labels => { 0 => "compact (all on one row)", 1 => "expanded (bump to multiple rows unless too many)" })
	."\n".
	p("Note: for performance reasons, individual HCNEs are not shown on genomic segments above ".
	$cne_display_threshold/1e6." Mb in size.").
	hr().
	"\n";

  # Return form
  return $form;
}

# called by gbrowse to annotate the DNA, returning features
sub annotate {
    my ($self, $segment) = @_;

    $self->_set_config_to_default() unless(($self->configuration->{version} || 0) eq $MY_VERSION);
    warn "I am in the annotate CNE now";
    $Data::Dumper::Indent = 1;
    $Data::Dumper::Sortkeys = 1;
    open (FH, '>/mnt/biggley/home/gtan/debug/debug_CNE_annotate_self.txt');
    print FH Dumper($self);
    close FH;
    open (FH, '>/mnt/biggley/home/gtan/debug/debug_CNE_annotate_segment.txt');
    print FH Dumper($segment);
    close FH;

    # Create a FeatureFile object to store the result in
    my $out_feature_list = Bio::Graphics::FeatureFile->new;
    $out_feature_list->mtime(time()); # To prevent gbrowse_img from retrieving a cached image
    
    return $out_feature_list unless($segment);
    warn "In annotate, still running";
    # Connect to DB
    my $db = CNE::DB->connect(dbhost => $DB_HOST,
			      dbname => $DB_NAME,
			      dbuser => $DB_USER)
	or die "could not connect to db $DB_NAME @ $DB_HOST";

    warn "In annotate,  the db is   ", Dumper($db);
    warn "In annotate, the out_feature_list is   ", Dumper($out_feature_list);
    my ($chr, $start, $end) = $self->_get_display_bounds($segment);
    warn "In annotate, the chr, start, end  ", $chr, " ", $start, "  ", $end;
    # Make the density features
    $self->_make_density_features($out_feature_list, $db, $chr, $start, $end);
    warn "In annotate, Am I running??";
    open (FH, '>/mnt/biggley/home/gtan/debug/debug_CNE_annotate_self_after_density_features.txt');
    print FH Dumper($self);
    close FH;
    open (FH, '>/mnt/biggley/home/gtan/debug/debug_CNE_annotate_out_feature_list_after_density_features.txt');
    print FH Dumper($out_feature_list);
    close FH;
    # Make the CNE features
    $self->_make_cne_features($out_feature_list, $db, $chr, $start, $end);
    open (FH, '>/mnt/biggley/home/gtan/debug/debug_CNE_annotate_self_after_CNE_features.txt');
    print FH Dumper($self);
    close FH;
    open (FH, '>/mnt/biggley/home/gtan/debug/debug_CNE_annotate_out_feature_list_after_CNE_features.txt');
    print FH Dumper($out_feature_list);
    close FH;
    # For debug use
    my $msg_string = "";
    # For debug usage
    if($msg_string) {
	$out_feature_list->add_type('message' => {glyph => 'text_in_box',
						  key => 0,
						  text => $msg_string});
	my $message = Bio::Graphics::Feature->new(-start => $segment->start,
						  -end => $segment->end);
	$out_feature_list->add_feature($message, 'message');
    }
    
    return $out_feature_list;
}


sub _get_display_bounds
{
    my ($self, $segment) = @_;

    my $id = $segment->seq_id;
    my $start = $segment->start;
    my $end = $segment->end;

    # Make sure start and end are not out of bounds
    # This can happen e.g. if user zooms out a lot
    my $whole_segment = $self->database->segment($id);
    my $min_start = $whole_segment->start;
    my $max_end = $whole_segment->end;
    $start = $min_start if($start < $min_start);
    $end = $max_end if($end > $max_end);

    return($id, $start, $end);
}


sub _get_sets_to_show
{
    warn "I am in _get_sets_to_show now!!!";

    my ($self, $prefix) = @_;
    my @set_nrs;
    for my $i (1..$self->nr_cne_sets()) {
	push @set_nrs, $i if($self->configuration->{$prefix.$i});
    }   
    return \@set_nrs;
}


sub _get_table_names
{
    warn "I am in _get_table_names now!!!";
    my ($self, $db, $set_nrs) = @_;
    my $asm1 = $self->asm1();
    my $asm2 = $self->asm2();
    my @tables;
    warn "In _get_table_names, the asm1 is  ", $asm1;
    warn "In _get_table_names, the asm2 is  ", $asm2;
    foreach my $i (@$set_nrs) {
        warn "In _get_table_names, doing ", $i;
        warn "The min_cne_id is", $self->configuration->{"min_cne_id$i"};
	my $id_cutoff_code = $self->configuration->{"min_cne_id$i"} or next;
	my ($min_id, $min_len) = split '/', $id_cutoff_code;
    warn "The min_id is ", $min_id;
    warn "The min_len is ", $min_len;
	my $table = $db->get_cne_table_name(assembly1 => $asm1, 
					    assembly2 => $asm2,
					    min_identity => $min_id,
					    min_length => $min_len,
	                                    version => 1);
    warn "In _get_table_names, check the table", Dumper($table);
	push @tables, $table if($table and $db->cne_table_exists($table));
    }
    return \@tables;
}


sub _get_min_lengths
{
    my ($self, $set_nrs) = @_;
    my @min_lengths;
    foreach my $i (@$set_nrs) {
	push @min_lengths, $self->configuration->{"min_cne_len$i"} || $DEFAULT_MIN_CNE_LENGTH;
    }
    return \@min_lengths;
}


sub _make_density_features
{
    warn "I am in _make_density_features now!!!";
    my ($self, $out_feature_list, $db, $chr, $start, $end) = @_;
    
    # For benchmarking
    my $t0;
    
    # Get length of the displayed segment 
    my $length = $end - $start + 1;

    # Get parameters for this plugin
    my $sets_to_show = $self->_get_sets_to_show('show_graph');
    warn "In _make_density_features, sets_to_show is ", @$sets_to_show;
    my $tables = $self->_get_table_names($db, $sets_to_show);
    warn "In _make_density_features, see tables is", Dumper($tables);
    my $min_lengths = $self->_get_min_lengths($sets_to_show);
    my $window_size = $self->configuration->{window_size}*1000;  # convert window size value from kb to bp
    my $max_score = $self->configuration->{max_score};
    my $nr_graphs = $self->configuration->{nr_graphs};
    warn "In _make_density_features, nr_graphs is ", $nr_graphs;
    my $rank_un_low = $self->configuration->{rank_un_low};
    warn "In _make_density_features, am I running? The rank_un_low ", $rank_un_low;

    my $label = $self->asm2_name . " HCNE density";
    warn "In _make_density_features, am I running? The label ", $label;
    warn "In _make_density_features, check the browse_config";
    open (FH, '>/mnt/biggley/home/gtan/debug/debug_CNE_annotate_self_browser_config.txt');
    print FH Dumper($self->browser_config);
    close FH;
    #my $pixel_width = $self->browser_config->width;
    my $pixel_width = 800;
    warn "In _make_density_features, am I running? The pixel_width  ", $pixel_width;
    my $assembly = $self->browser_config->name;
    warn "In _make_density_features, am I running? The assembly ", $assembly;
    #my $assembly = $self->browser_config->source;
    my ($my_last_name) = ref($self) =~ /(\w+)$/;

    # Return if no nothing to do
    return if($nr_graphs == 0 or @$sets_to_show == 0);
    # Compute the step size
    $t0 = Benchmark->new;
    my $step_size;
    if($length <= $pixel_width) {
	# The lower limit for the step size is 1 base
	$step_size = 1;
    } 
    else {
	# If we are not not zoomed in to base resolution, 
	# we want at least 1 step / pixel and at least 10 steps / window.
	$step_size = int($length / $pixel_width);
	$step_size = int($window_size/10) if($step_size > $window_size/10);
	# Adjust the step size so that window_size is evenly divisible by it.
	while($window_size % $step_size) {  # there must be a smarter algorithm for this...
	    $step_size--;
	}
    }
    warn "In _make_density_features, am I running? the step_size  ", $step_size;
    my $win_nr_steps = $window_size / $step_size;
#    print STDERR "Calculated step size ".timediff_str($t0,Benchmark->new)."\n";


    # Compute how large region we need to look at around the displayed region
    my $context_start = $start - int((($win_nr_steps-1)*$step_size)/2+0.5);
    $context_start  = 1 if($context_start < 1);
    warn "In _make_density_features, the context_start is ", $context_start;
    my $context_end = $end + int((($win_nr_steps-1)*$step_size)/2+$step_size+0.5);
    warn "In _make_density_features, the context_end is ", $context_end;
   
    # Define a feature type with drawing parameters for the curve
    $out_feature_list->add_type('density_curve' => {glyph => 'fast_xyplot',
						    graph_type => 'line',
						    key   => $label,
						    bgcolor => 'blue',
						    height => 50,
						    min_score => 0,
						    max_score => $max_score,					 
						    scale => 'left',
						    clip => 1
						});
    

    # Get features that overlap the region and convert them to a set of non-overlapping ranges
    $t0 = Benchmark->new;
    my %range_lists_by_graphs;
    foreach my $table (@$tables) {
	my $min_length = shift @$min_lengths;
	my @get_cne_args = (table_name => $table,
			    assembly => $assembly,
			    chr => $chr,
			    start => $context_start,
			    end => $context_end,
			    min_length => $min_length);
	if($nr_graphs > 1) {
	    my $ranges_by_chr2 = $db->get_cne_ranges_in_region_partitioned_by_other_chr(@get_cne_args);
	    while(my ($chr2, $ranges) = each %$ranges_by_chr2) {
		push @{$range_lists_by_graphs{$chr2}}, collapse_intervals($ranges); # Collapse overlapping ranges
	    }
	}
	else {
	    my $ranges = $db->get_cne_ranges_in_region(@get_cne_args);
	    push @{$range_lists_by_graphs{'all'}}, collapse_intervals($ranges);  # Collapse overlapping ranges
	}
    }
#    print STDERR "Got features ".timediff_str($t0,Benchmark->new)."\n";

    organize_range_lists(\%range_lists_by_graphs, $nr_graphs, $rank_un_low);

    my $chr2_list = sort_chr([keys %range_lists_by_graphs]);
    foreach my $chr2 (@$chr2_list) {
	my $range_lists = $range_lists_by_graphs{$chr2};

	# Calculate the window scores
	$t0 = Benchmark->new;
	my ($first_score_start, @win_scores_list);
	foreach my $ranges (@$range_lists) {
	    my $win_scores;
	    ($first_score_start, $win_scores) = calc_window_scores($start, $end, $ranges, $win_nr_steps, $step_size);
	    push @win_scores_list, $win_scores;
	}
#	print STDERR "Calculated window scores ".timediff_str($t0,Benchmark->new)."\n";
	
	# Create a Bio::Graphics::Feature object for the curve and add it to the FeatureFile object
	my %attributes;
	$attributes{'score_start'} = $first_score_start;
	$attributes{'score_span'} = $step_size;
	$attributes{'score_list'} = \@win_scores_list;
	$attributes{'color_list'} = $self->_get_colors($sets_to_show);
	my $density_curve = Bio::Graphics::Feature->new(-seq_id => $chr,
							-start=> $start,
							-end  => $end,
							-display_name => $chr2 eq 'all' ? undef : "Density $chr2",
							-type=> 'CneDensity',
							-strand => '+',
							-url => '?plugin_action=Go&plugin='.$my_last_name,
							-attributes => \%attributes);
	$out_feature_list->add_feature($density_curve,'density_curve');
    }

}


sub _make_cne_features
{
    warn "I am in _make_cne_features now!!!";
    my ($self, $out_feature_list, $db, $chr, $start, $end) = @_;
    
    # For benchmarking
    my $t0;
    
    # Get length of the displayed segment ; return if too long
    my $length = $end - $start + 1;
    my $cne_display_threshold = $self->static_plugin_setting('cne_display_threshold') || $DEFAULT_CNE_DISPLAY_THRESHOLD;
    return if($length > $cne_display_threshold);   

    # Get settings
    my $sets_to_show = $self->_get_sets_to_show('show_cnes');
    my $tables = $self->_get_table_names($db, $sets_to_show);
    my $min_lengths = $self->_get_min_lengths($sets_to_show);
    my $labels_string = $self->static_plugin_setting('cne_labels');   
    my @labels = $labels_string ? split /\s+/, $labels_string : ();
    my $color_cnes = $self->configuration->{color_cnes};
    my $bump_cnes = $self->configuration->{bump_cnes};
#    my $assembly = $self->browser_config->source;
    my $assembly = $self->browser_config->name;
    # Prepare to color CNEs
    my ($color_sub, $seq_color_index);
    if($color_cnes) {
	$seq_color_index = AT::GFX::SeqColorIndex->new(scheme => $color_cnes);
	if($seq_color_index) {
	    $color_sub = sub { shift->attributes('color') };
	}
	else {
	    $color_cnes = 0;
	    $color_sub = 'black';
	}
    }
    else {
	$color_sub = 'black'
    }

    # Get features that overlap the displayed region from database
    $t0 = Benchmark->new;
    my $n = 0;
    foreach my $i (0..$#$tables) {
	my $table = $tables->[$i];
	my $min_length = $min_lengths->[$i];
	my $track_id = "cne$i";
	my $table_info = $db->get_cne_table_info($table);
	my $label = $self->asm2_name.' HCNEs (min. identity: '.($table_info->{min_identity}*100).'% over '.
	    $table_info->{min_length}.' columns, min. size: '.$min_length.' bp)';
	my $cnes = $db->get_cnes_in_region(table_name => $table,
					   assembly => $assembly,
					   chr => $chr,
					   start => $start,
					   end => $end,
					   min_length => $min_length);
	$out_feature_list->add_type($track_id => {glyph => 'generic',
						  key   => $label,
						  bgcolor => $color_sub,
						  fgcolor => $color_sub,
						  $bump_cnes ? () : (bump => 0, label => 0)
				    });

	foreach my $cne (@$cnes) {
	    my %attributes;
	    if($color_cnes) {
		$attributes{color} = sprintf("#%02x%02x%02x", $seq_color_index->get_seq_color($cne->chr2));
	    }
	    my $url = "../../cne_details?type=cne&table=$table&asm=$assembly&id=".$cne->id;
	    my $cne_name = $cne->chr2.':'.$cne->start2;
	    my $feature = Bio::Graphics::Feature->new(-seq_id => $chr,
						      -start=> $cne->start1,
						      -end  => $cne->end1,
						      -type => 'HCNE',
						      -strand => $cne->strand,
						      -attributes => \%attributes,
						      -url    => $url,
						      -name  => $cne_name
		);
	    $out_feature_list->add_feature($feature, $track_id);
	    $n++;
	}
	unless(@$cnes) {
	    # If there are no HCNEs, add a dummy element so the emtpy track is drawn
	    my $feature = Bio::Graphics::Feature->new(-seq_id => $chr,
						      -start=> 0,
						      -end  => 0,
						      -type => 'dummy',
						      -strand => '+'
		);
	    $out_feature_list->add_feature($feature, $track_id);
	}
    }
#    print STDERR "Made $n individual CNE features".timediff_str($t0,Benchmark->new)."\n";
}


sub organize_range_lists
{
    my ($indexed_range_lists, $nr_graphs, $rank_un_low) = @_;
    my @all_keys = keys %$indexed_range_lists;
    my $nr_keys = scalar @all_keys;

    # Check for easy cases
    return if($nr_keys <= $nr_graphs);
    if($nr_graphs < 1) { 
	$indexed_range_lists = {};
	return;
    }

    # Split into two key lists if we should rank un-assembled chrom low
    my (@key_list1, @key_list2);
    if($rank_un_low) {
	foreach my $key (@all_keys) {
	    if($key =~ /chrU|chrNA|random|hap|Unknown/) {
		push @key_list2, $key;
	    }
	    else {
		push @key_list1, $key;
	    }
	}
    }
    else {
	@key_list1 = @all_keys;
    }

    # Sort the sets of range lists by how much CNE sequence they span
    my %spans;
    while(my ($key, $range_lists) = each %$indexed_range_lists) {
	$spans{$key} = AT::Tools::RangeHandler->sum($range_lists->[0]);
	# ^ We may want to intersect with the displayed region here...
    }

    my @sorted_keys = (sort({ $spans{$b} <=> $spans{$a} } @key_list1),
		       sort({ $spans{$b} <=> $spans{$a} } @key_list2));
    
    # Merge ranges for the last graph
    my @merged_range_lists;
    foreach my $key (@sorted_keys[$nr_graphs-1..$nr_keys-1]) {
	my $range_lists = $indexed_range_lists->{$key};
	delete $indexed_range_lists->{$key};
	for my $i (0..@$range_lists-1) {
	    push @{$merged_range_lists[$i]}, @{$range_lists->[$i]};
	}
    }
    for my $i (0..@merged_range_lists-1) {
	$merged_range_lists[$i] = collapse_intervals($merged_range_lists[$i]);
    }

    $indexed_range_lists->{'other'} = \@merged_range_lists;
}


sub collapse_intervals
{
    my ($exons_ref) = @_;

    return [] unless($exons_ref and @$exons_ref);

    my @hsps = sort {$a->[0] <=> $b->[0]} @$exons_ref;

    my @iv;
    my ($start, $end) = @{$hsps[0]};
    for my $i (1..@hsps-1) {
	my ($my_start, $my_end) = @{$hsps[$i]};
	if($my_start > $end) {
	    push @iv, [$start,$end];
	    ($start, $end) = ($my_start, $my_end); 
	}
	else {
	    $end = $my_end if ($end < $my_end);
	}
    }
    push @iv, [$start,$end];
    return \@iv;
}


# CODE TO COMPUTE RANDOM CURVE
#  for(my $pos = $abs_start; $pos <= $end; $pos += $step_size) {
#      my $value = int(rand(1)*101);
#      my $feature = Bio::Graphics::Feature->new(-start=>$pos,
#						-end  =>$pos+$step_size-1,
#						-type => 'density_value',
#						-score => $value,
#						-strand => '+',
#						);
#      $density_curve->add_segment($feature);
#  }


# 1. Get features in [start-win_size/2, end+win_size/2]
# 2. If there is a core region included in all blocks that are to be drawn?
#      2.1. Compute scores up to start of core region (incrementally)
#            store these in one array
#      2.1. Compute score for core region, store as a scalar value
#      2.2. Compute scores for blocks to be shown, by subtracting a pre-score and adding a post-score, incrementally
#      Memory usage is about 2*nr of blocks
# 3. else
#      3.1. Apply the current algorithm.
#      Memory usage is < 3*nr of blocks 


sub calc_window_scores
{
    my ($display_start, $display_end, $ivs, $win_nr_blocks, $blk_size) = @_;
    # start and end are start of region to be shown. need to calc context start
 
    my $display_size = $display_end - $display_start + 1;
    my $win_size = $win_nr_blocks * $blk_size;
    my $max_score = $blk_size * $win_nr_blocks;

    my $offset = int((($win_nr_blocks-1)*$blk_size)/2+0.5);
    my $compute_start = $display_start - $offset;
    my $compute_end = $display_end + $offset;
    my $compute_size = $compute_end - $compute_start;

    my $first_score_start; 
    my @scores;

    if($win_size > 2*$display_size) {
	# If win_size is much larger than display_size we can do some optimization.
	# by separately computing scores for a core region common to all displayed windows
	# This algorithm should work also in cases where the display region does no contain an even nr of blocks

	# Compute number of blocks n in display region. This is the same as the number of displayed window scores.
	my $nr_blocks = int($display_size / $blk_size) + ($display_size % $blk_size ? 1 : 0); 
	 # ^ add an extra block in case the region does not contain an even number of blocks

	# Compute n-1 blocks from compute_start,
	# keep track of individual block scores in @pre_blk_scores and their sum in $win_score
	my ($blk_start, $blk_end) = ($compute_start, $compute_start + $blk_size-1);
	my @pre_blk_scores = (0) x ($nr_blocks-1);
	my $win_score = 0;
	my $j = 0;   # $j = interval index
	for my $i (0..($nr_blocks-2)) {
	    my $blk_score = calc_block_score($blk_start, $blk_end, \$j, $ivs);
	    $win_score += $blk_score;
	    $pre_blk_scores[$i] = $blk_score;
	    $blk_start += $blk_size;
	    $blk_end += $blk_size;
	}

	# Compute the score for the first displayed window. This is the score for the win_nr_blocks first blocks.
	# We compute it as the score for the first n-1 blocks (already computed, in $win_score)
	# plus the score for the next win_nr_blocks - (n-1) blocks (the "core" region shared by all displayed windows)
	$blk_end += $blk_size * ($win_nr_blocks - $nr_blocks);
	$win_score += calc_block_score($blk_start, $blk_end, \$j, $ivs);
	push @scores, $ABS_SCORE ? $win_score : 100 * $win_score / $max_score;
	$blk_start += $blk_size * ($win_nr_blocks - $nr_blocks);
	$first_score_start = $blk_start - $offset;
	#die "first_score_start $first_score_start < display_start $display_start" unless($first_score_start >= $display_start);

	# Compute the scores for the remaining (n-1) displayed windows
	# Window 2 has score of window 1 - score of block 1 + score of the first block after the core region
	# Window 3 has score of window 2 - score of block 2 + score of the second block after the core region
	# etc

	# Debug checks
	#die "error: blk_start $blk_start != display_start $display_start + offset $offset"
	    #if($blk_start != $display_start + $offset);
#	print STDERR "blk_start $blk_start, display_start $display_start, offset $offset, blk_size $blk_size\n";
	die "error: blk_start $blk_start + $blk_size - 1 != $blk_end"
	    if($blk_start + $blk_size -1 != $blk_end);

	$blk_start += $blk_size;
	$blk_end += $blk_size;
	my $prev_win_score = $win_score;
	for my $i (0..($nr_blocks-2)) {
	    my $blk_score = calc_block_score($blk_start, $blk_end, \$j, $ivs);
	    $win_score = $prev_win_score + $blk_score - $pre_blk_scores[$i];
	    push @scores, $ABS_SCORE ? $win_score : 100 * $win_score / $max_score;
	    $prev_win_score = $win_score;
	    $blk_start += $blk_size;
	    $blk_end += $blk_size;
	}
    }
    else {

	# Compute number of blocks in the computed region.
	# This is the same as the number of windows scores we will compute.
	my $nr_blocks = int($compute_size / $blk_size) + ($compute_size % $blk_size ? 1 : 0); 
	 # ^ add an extra block in case the region does not contain an even number of blocks

	# Allocate at space for at least win_nr_blocks+1 blocks. This is needed to make sure
	# we don't go a full circle around the @blk_scores array when we find dead blocks below
	my @blk_scores = (0) x ($nr_blocks > $win_nr_blocks ? $nr_blocks : $win_nr_blocks + 1);

	my ($blk_start, $blk_end) = ($compute_start, $compute_start + $blk_size-1);
	my $prev_win_score = 0;
	my $j = 0;   # $j = interval index

	# Compute scores for all windows in the computed region.
	# The first windows (those outside the display) region will have incomplete scores
	# and will not be kept.
	for my $i (1..$nr_blocks) {  # $i = block index
	    my $blk_score = calc_block_score($blk_start, $blk_end, \$j, $ivs);
	    my $dead_block = $i - $win_nr_blocks;
	     # ^ the size of @blk_scores must be larger than $win_nr_blocks for this to work; see above
	    my $win_score = $prev_win_score + $blk_score - $blk_scores[$dead_block];
	    if($blk_start >= $display_start + $offset) {   
		$first_score_start = $blk_start - $offset unless(defined $first_score_start);
		push @scores, $ABS_SCORE ? $win_score : 100 * $win_score / $max_score;
	    }
	    $blk_scores[$i] = $blk_score;
	    $prev_win_score = $win_score;
	    $blk_start += $blk_size;
	    $blk_end += $blk_size;
	}
    }

    return ($first_score_start, \@scores);
}


sub calc_block_score
{
    my ($blk_start, $blk_end, $j_ref, $ivs) = @_;
    my $blk_score = 0;
    my $j = $$j_ref;
    for(;$j < @$ivs; $j++)
    {
	my ($iv_start, $iv_end) = @{$ivs->[$j]};
	last if($iv_start > $blk_end);
	if($iv_end >= $blk_start) {
	    my $ol_start = ($iv_start > $blk_start) ? $iv_start : $blk_start;
	    my $ol_end = ($iv_end < $blk_end) ? $iv_end : $blk_end;
	    $blk_score += $ol_end - $ol_start + 1;
	    last if($iv_end > $blk_end);
	}
    }
    $$j_ref = $j;
    return $blk_score;
}


sub sort_chr
{
    my ($names) = @_;

    return $names if(@$names < 2);

    my %name2num;

    # make conversion index (name -> number)
    foreach my $name (@$names) {
	my $nr;
	if($name =~ /^(?:chr|group)(.+)/i) {
	    $nr = $1;
	    if($nr =~ /^(\d+)/) { 
		$nr = $1;
		$nr += .5 if($name =~ /random/);
	    }
	    elsif($nr =~ /^U/) { $nr = 104; }
	    elsif($nr =~ /^X/) { $nr = 101; }
	    elsif($nr =~ /^Y/) { $nr = 102; }
	    elsif($nr =~ /^M/) { $nr = 103; }
	    elsif($nr =~ /^V/) { $nr = 5; }
	    elsif($nr =~ /^IV/) { $nr = 4; }
	    elsif($nr =~ /^III/) { $nr = 3; }
	    elsif($nr =~ /^II/) { $nr = 2; }
	    elsif($nr =~ /^I/) { $nr = 1; }
	    else{ $nr = 110; }
	}
	elsif($name =~ /^(?:scaffold|contig|supercont|super)_?(\d+)/i or $name =~ /^(\d+)/) {
	    $nr = $1;
	}
	elsif($name eq 'other') {
	    $nr = 1e6;
	}
	else {
	    $nr = 1e6-1;
	}
	$name2num{$name} = $nr;
    }

    # sort the names by their numbers
    my @sorted_names = sort { $name2num{$a} <=> $name2num{$b} } @$names;

    return \@sorted_names;
}


sub timediff_str
{
    my ($b1, $b2) = @_;
    my $t = timestr(timediff($b2, $b1));
    return "[$t sec]";
}


# Do

# investigate:
# can we have on-the-fly-drawn density tracks in the overview window? Could not find a way to do this.
# can gbrowse show entire chromosomes in the detail window? YES.
# note: we can buffer whole-chrom density tracks that are reused often if computing them takes time

# now:
# - try to rewrite xyplot
#     get the resolution (pixels/bp)
#     use jk's algorithm to find the value for each pixel
#     draw a filled polygon
#     try to draw multiple polygons on top of each other
# For fast drawing:
#  In this method:
#   - create just one Bio::Graphics::Feature object
#   - give it an attribute that holds the @win_data (start, end, score);
#     alternatively, give: start, block_size, and array of scores
#  In the glyph:
#   - obtain pixel boundary coords from the single feature object
#   - create and fill an array of y values, one per pixel
#   - draw the array, e.g. using the polygon function
# If further optimization is needed, we can compute the max score for each pixel directly in the score computation function
#  (probably we won't gain very much from this and it seems tricky to get the resolution info)

# To get bounds in pixels from panel object:
#       $left = $panel->left
#       $right = $panel->right
#       $top   = $panel->top
#       $bottom = $panel->bottom
# but how do we get a panel object? maybe it does not exist when plugin is called...


# test:
# adjustable window and step size: what should be set automatically?
#  test: automatic step size and adjustable window size, try 100-200 steps
# import hcne data (need to keep query chromosome) 
# split the plot up into several tracks, one for each query chromosome 


  # Then think about:
  # how we want the look-and-feel to be. several tracks, or one configurable track?
  # how can configure button be reached? through about button?
  # can we place it in a different section than "Analysis"?
  # do we need to optimize xyplot?
  # shall we write a custom xyplot to draw overlayed densities?

  # note: there is a glyph type "image" - we could write a C extension for the drawing if needed
  # maybe the UCSC code could be used :) 


# Main drawing code from jk:

#x1 = ((wi->start-seqStart) + (dataOffset * usingDataSpan)) * pixelsPerBase;
#x2 = x1 + (usingDataSpan * pixelsPerBase);
#for (i = x1; i <= x2; ++i)
#{
#    int xCoord = preDrawZero + i;
#    if ((xCoord >= 0) && (xCoord < preDrawSize))
#    {
#	double dataValue =
#	    BIN_TO_VALUE(datum,wi->lowerLimit,wi->dataRange);
#	
#	++preDraw[xCoord].count;
#	if (dataValue > preDraw[xCoord].max)
#	    preDraw[xCoord].max = dataValue;
#	if (dataValue < preDraw[xCoord].min)
#	    preDraw[xCoord].min = dataValue;
#	preDraw[xCoord].sumData += dataValue;
#	preDraw[xCoord].sumSquares += dataValue * dataValue;
#    }
#}
#}


1;




