# AT::WWW::ComparativeGenomeBrowser module
#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::WWW::ComparativeGenomeBrowser - comparative genome browser web-interface

=head1 SYNOPSIS

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::WWW::ComparativeGenomeBrowser;

use strict;
use vars '@ISA';
use CGI;
use File::Temp;
use AT::Root;
use AT::GFX::AlignedLoci;
use AT::GFX::Locus;
use AT::GFX::LocusConfig;

@ISA = qw(AT::Root);

# defaults for user-configurable options
my $DEF_IMG_WIDTH = 1000;
my $DEF_CONN_HEIGHT = 50;
my $DEF_MAX_SEC_SPACE = 1.5;
my $DEF_MIN_CONTEXT = 0;
my @DEF_LOC = ('chr1', '100000', '200000');
my $DEF_LOC_STR = 'chr1:100,000-200,000';

# non-user-configurable options
my $MIN_EXON_SIZE = 0;
my $OUTPUT_SVG = 0;

sub new
{
    my ($caller, %args) = @_;
    my $config = $args{config} || die "No config arg";
    my $fn = $config->{drawing_config_file} || die "No drawing_config_file parameter in config";
    $config->{drawing_config} = AT::GFX::LocusConfig->new(file => $fn);
    my $self = bless { config => $config }, ref $caller || $caller;
    return $self;   
}


sub html_dir
{
    shift->{config}->{html_dir} || die "No html_dir configured\n";
}


sub tmp_dir
{
    shift->{config}->{tmp_dir} || die "No tmp_dir configured\n";
}


sub tmp_dir_url
{
    shift->{config}->{tmp_dir_url} || die "No tmp_dir_url configured\n";
}


sub drawing_config
{
    shift->{config}->{drawing_config} || die "No drawing_config configured\n";
}


sub handle_request
{
    my ($self, $q) = @_;

    # html header
    print $q->header();
    print $q->start_html(-title => "Comparative Genome Browser",
			 -script => $self->_javascript_header());
    print $q->h2("Comparative Genome Browser");
    print $q->start_form(-name => 'form', -method=>'GET');

    # Get assemblies
    my ($asm1, $asm2, $asm3) = $self->_get_assemblies($q);

    if($asm1) {
	# If script run with some parameters
	# Get locations
	my ($location_result, $loc1, $loc2, $loc3) = $self->_get_locations($q, $asm1, $asm2, $asm3);
	if($location_result eq 'ok') {	
	    # If locations are unambiguous: draw
	    my $drawer = $self->_draw_loci($q, $asm1, $asm2, $asm3, $loc1, $loc2, $loc3);
	    $self->_print_upper_controls($q, $asm1, $asm2, $asm3);
	    if($drawer) {
		my ($img_fh, $img_fn) = $self->_new_tmp_file();
		print $img_fh $drawer->gd_image->png;
		close $img_fh;
		print $q->hr;
		print $q->img({-src => $img_fn});
	    }
	}
	else {
	    # Else print the controls
	    $self->_print_upper_controls($q, $asm1, $asm2, $asm3);
	    if($location_result eq 'search') {
		print $q->hr;
		# Output search results if any (if none there was an error, and message should already have been printed)
	    }
	}
    }
    else {
	# Script was run without parameters
	# exectute javascript to set defaults?
	($asm1) = values %{$self->drawing_config->assemblies};
	$q->param("asm1",$asm1->id);
	$q->param("loc1active","on");
	$q->param("loc1",$DEF_LOC_STR);
	$self->_print_upper_controls($q, $asm1, $asm2, $asm3);
    }
    
    # Output more controls
    print $q->hr;
    print $q->a({-name => 'options'});
    $self->_print_lower_controls($q, $asm1, $asm2, $asm3);

    # Html footer
    print $q->end_form;
    print $q->end_html;

    $self->_clean_tmp_files();
}


sub _draw_loci
{
    my ($self, $q, $asm1, $asm2, $asm3, $loc1, $loc2, $loc3) = @_;

    # Get drawing parameters from the form
    my @drawer_args = (width => $q->param('imgwidth'),
		       connector_height => $q->param('connheight'),
		       min_exon_size => $MIN_EXON_SIZE,
		       max_nr_secondary_panels => $q->param('secpanels'),
		       max_space_factor => $q->param('maxsecspace'),
		       min_context => $q->param('mincontext'),
		       svg => $OUTPUT_SVG);

    # Generate a drawing object based on inputs
    my $drawer;
    if($loc1 and $loc3 and $asm2) {
       # If button 3 and all assemblies selected: show sandwich in
	my @asm_params = ($asm1->get_drawing_parameters,
			  $asm2->get_drawing_parameters,
			  $asm3->get_drawing_parameters);
	if($asm_params[0]->{alignment_db}->alignment_exists($asm2->id)) {
	    if($asm_params[2]->{alignment_db}->alignment_exists($asm2->id)) {
		$drawer = AT::GFX::AlignedLoci->new(assembly_param => \@asm_params, @drawer_args);	
		$drawer->draw_locus([[0, @$loc1], [2, @$loc3]]);
	    }
	    else {
		print $q->p("No alignment from ".($asm3->id)." to ".($asm2->id)."!");
	    }
	}
	else {
	    print $q->p("No alignment from ".($asm1->id)." to ".($asm2->id)."!");
	}
    }
    elsif($loc2) {
       # elsif button 2: show reversed pair or sandwich out 
	if($asm3) {
	    # show sandwich out
	    my @asm_params = ($asm1->get_drawing_parameters,
			      $asm2->get_drawing_parameters,
			      $asm3->get_drawing_parameters);
	    if($asm_params[1]->{alignment_db}->alignment_exists($asm1->id)) {
		if($asm_params[1]->{alignment_db}->alignment_exists($asm3->id)) {
		    $drawer = AT::GFX::AlignedLoci->new(assembly_param => \@asm_params, @drawer_args);	
		    $drawer->draw_locus([[1, @$loc2]]);
		}
		else {
		    print $q->p("No alignment from ".($asm2->id)." to ".($asm3->id)."!");
		}
	    }
	    else {
		print $q->p("No alignment from ".($asm2->id)." to ".($asm1->id)."!");
	    }
	}
	else {
	    # show reversed pair
	    my @asm_params = ($asm1->get_drawing_parameters,
			      $asm2->get_drawing_parameters);
	    if($asm_params[1]->{alignment_db}->alignment_exists($asm1->id)) {
		$drawer = AT::GFX::AlignedLoci->new(assembly_param => \@asm_params, @drawer_args);
		$drawer->draw_locus([[1, @$loc2]]);
	    }
	    else {
		print $q->p("No alignment from ".($asm2->id)." to ".($asm1->id)."!");
	    }
	}
    }
    elsif($loc1) {
	# else: show pair or single
	if($asm2) {
	    # show pair
	    my @asm_params = ($asm1->get_drawing_parameters,
			      $asm2->get_drawing_parameters);
	    if($asm_params[0]->{alignment_db}->alignment_exists($asm2->id)) {
		$drawer = AT::GFX::AlignedLoci->new(assembly_param => \@asm_params, @drawer_args);	
		$drawer->draw_locus(@$loc1);
	    }
	    else {
		print $q->p("No alignment from ".($asm1->id)." to ".($asm2->id)."!");
	    }
	}
	else {
	    # show single
	    my $asm_params = $asm1->get_drawing_parameters;
	    $drawer = AT::GFX::Locus->new(%$asm_params, @drawer_args);
	    $drawer->draw_locus(@$loc1);
	}
    }

    if($drawer and $drawer->errstr) {
	print $q->p($drawer->errstr);
	$drawer = undef;
    }

    return $drawer;
}


sub _get_assemblies
{
    my ($self, $q) = @_;

    my @assemblies;

    for my $i (1..3) {
	my $asm_id = $q->param("db$i");
	my @sel_tracks = $q->param("tracks$i");
	my $asm;
	if($asm_id and $asm_id ne 'none') {
	    $asm = $self->drawing_config->assemblies->{$asm_id} or die "Assembly $asm_id not configured";
	    $asm->enabled_tracks({ map {$_ => 1} @sel_tracks });
	}
	push @assemblies, $asm;
    }

    return @assemblies;
}


sub _print_upper_controls
{
    my ($self, $q, @assemblies) = @_;

    my @no_assembly;

    print $q->start_table();
    print $q->Tr({-align => 'left'},
		 $q->th(['Assembly','', 'Location','']),
		 $q->th({-colspan => 2}, 'Move'),
		 $q->th({-colspan => 2}, 'Zoom in'),
		 $q->th({-colspan => 1}, 'Zoom out')
		 );

    for my $i (1..3) {
	my $asm = $assemblies[$i-1];
	my @asm_active = (disabled => 1) unless ($asm);
	my @loc_active = (disabled => 1) unless ($q->param("loc${i}active"));
	print $q->Tr($q->td([
			     $q->popup_menu(-name => "db$i",
					    -values => [@no_assembly, keys %{$self->drawing_config->assemblies}],
					    -onChange => "javascript:populateData($i, this.options[selectedIndex].text)"),
			     $q->checkbox(-name => "loc${i}active",
					  -label => '',
					  -onClick => "javascript:locActivityChange($i)",
					  @asm_active),
			     $q->textfield(-name => "loc$i",
					   -size => 30,
					   -maxlength => 50,
					   -onKeyDown => "javascript:refreshIfEnterPressed(event)",
			   #-onSubmit => "return confirm('Are you sure you want to view this page?');",
					   @loc_active),
			     $q->pre(" "),
			     $q->submit(-name => "move${i}L3", -value => '<<<', @loc_active).
			     $q->submit(-name => "move${i}L2", -value => '<<', @loc_active).
			     $q->submit(-name => "move${i}L1", -value => '<', @loc_active).
			     $q->submit(-name => "move${i}R1", -value => '>', @loc_active).
			     $q->submit(-name => "move${i}R2", -value => '>>', @loc_active).
			     $q->submit(-name => "move${i}R3", -value => '>>>', @loc_active),
			     $q->pre(" "),
			     $q->submit("zoom${i}I1",'1.5x', @loc_active).
			     $q->submit("zoom${i}I2",'3x', @loc_active).
			     $q->submit("zoom${i}I3",'10x', @loc_active),
			     $q->pre(" "),
			     $q->submit("zoom${i}O1",'1.5x', @loc_active).
			     $q->submit("zoom${i}O2",'3x', @loc_active).
			     $q->submit("zoom${i}O3",'10x', @loc_active),
			     ])
		     );
	@no_assembly = ('none');
    }

    print $q->end_table;

    print $q->start_table;
    print $q->Tr($q->td([$q->submit('Refresh'),
			 $q->pre(" "),
			 'Nonsyntenic',
			 $q->popup_menu('secpanels',[-1,1..9],2,
					{-1,'Show all',1,'Hide all', 2, 'Max 1', 3,'Max 2', 4,'Max 3',
					 5,'Max 4', 6,'Max 5',7,'Max 6',8,'Max 7',9,'Max 8'}),
			 $q->pre(" "),
			 'Primary panel width',
			 $q->textfield(-name => 'imgwidth',
				       -default => $DEF_IMG_WIDTH,
				       -size => 4,
				       -maxlength => 4,
				       -onKeyDown => "javascript:refreshIfEnterPressed(event)"),
			 $q->pre(" "),
			 'Connector height',
			 $q->textfield(-name => 'connheight',
				       -default => $DEF_CONN_HEIGHT,
				       -size => 3,
				       -maxlength => 3,
				       -onKeyDown => "javascript:refreshIfEnterPressed(event)"),
			 $q->pre(" "),
			 'Max space factor',
			 $q->textfield(-name => 'maxsecspace',
				       -default => $DEF_MAX_SEC_SPACE,
				       -size => 4,
				       -maxlength => 4,
				       -onKeyDown => "javascript:refreshIfEnterPressed(event)"),
			 $q->pre(" "),
			 'Min context',
			 $q->textfield(-name => 'mincontext',
				       -default => $DEF_MIN_CONTEXT,
				       -size => 7,
				       -maxlength => 7,
				       -onKeyDown => "javascript:refreshIfEnterPressed(event)"),
			 #$q->pre(" "),
			 #$q->a({-href => ($self->html_dir).'/help.html'}, "Help")
			 ])
		 );
    print $q->end_table;
}


sub _print_lower_controls
{
    my ($self, $q, @assemblies) = @_;

    my @track_lists;
    for my $i (1..3) {
	my (@all_tracks, @sel_tracks);
	my $asm = $assemblies[$i-1];
	if($asm) {
	    @all_tracks = map { $_->{name} } @{$asm->tracks};	
	    @sel_tracks = keys %{$asm->enabled_tracks};
	}
	else {
	    @all_tracks = ('No assembly selected');
	}
	push @track_lists, $q->scrolling_list("tracks$i",\@all_tracks,\@sel_tracks,15,'true');
    }

    print $q->start_table();
    print $q->Tr($q->th({-align => 'left', -colspan => 2}, ['Track selection']));
    print $q->Tr($q->td({-align => 'left'}, ['Assembly 1&nbsp;&nbsp;', 'Assembly 2&nbsp;&nbsp;', 'Assembly 3']));
    print $q->Tr($q->td(\@track_lists));
    print $q->end_table;   

}


#sub _get_button_and_value
#{
#    my ($self, $q) = @_;
#    my @buttons = ("Refresh", "zoom in", "zoom out", "move");
#    foreach my $button (@buttons) {
#	my $value = $q->param($button);
#	return ($button, $value) if($value);
#    }
#    return ("","");
#}


sub _get_locations
{
    my ($self, $q, @assemblies) = @_;
    my @locations = (undef,undef,undef);

    for my $i (1..3) {
	my $asm = $assemblies[$i-1];
	if($asm and $q->param("loc${i}active")) {
	    my $loc_str = $q->param("loc$i");
	    if($loc_str) {
		my ($chr, $start, $end) = $self->_id_string_to_region($loc_str);
		if($chr) {
		    ($start,$end) = $self->_calc_range($q,$i,$start,$end);
		    if($start) {
			$locations[$i-1] = [$chr, $start, $end];
			$q->param("loc$i", $self->_loc_str($chr,$start,$end));
		    }
		    else {
			return ('error');
		    }
		}
		else {
		    print "Invalid location $loc_str";
		    return ('error');
		}
	    }
	    else {
		$locations[$i-1] = \@DEF_LOC;
		$q->param("loc$i", $DEF_LOC_STR);
	    }
	}
    }

    return ('ok', @locations);
}   


sub _loc_str
{
    my ($self, $chr, $start, $end) = @_;
    my (@start_triplets, @end_triplets);
    while(length($start)) {
	unshift @start_triplets, substr($start, -3, 3, '');
    }
    while(length($end)) {
	unshift @end_triplets, substr($end, -3, 3, '');
    }
    $start = join(',',@start_triplets);
    $end = join(',',@end_triplets);
    return "$chr:$start-$end";
}


sub _id_string_to_region
{
    my ($self, $str) = @_;
    my ($chr,$pos) = split /:/, $str;
    my ($start, $end) = split '-', $pos if($pos);
    return () unless($chr and $start and $end);
    $start =~ s/,//g;
    $end =~ s/,//g;
    return ($chr, $start, $end);
}


sub _calc_range
{
    my ($self, $q, $i, $start, $end) = @_;
    
    my $max_size = 5000000;  

    my $cmd = "";
    my $arg;
    if($q->param("move${i}L1")) {
	$cmd = 'move'; $arg = -0.1;
    }
    elsif($q->param("move${i}L2")) {
	$cmd = 'move'; $arg = -0.5;
    }
    elsif($q->param("move${i}L3")) {
	$cmd = 'move'; $arg = -0.95;
    }
    elsif($q->param("move${i}R1")) {
	$cmd = 'move'; $arg = 0.1;
    }
    elsif($q->param("move${i}R2")) {
	$cmd = 'move'; $arg = 0.5;
    }
    elsif($q->param("move${i}R3")) {
	$cmd = 'move'; $arg = 0.95;
    }
    elsif($q->param("zoom${i}O1")) {
	$cmd = 'zoom'; $arg = 1.5;
    }
    elsif($q->param("zoom${i}O2")) {
	$cmd = 'zoom'; $arg = 3;
    }
    elsif($q->param("zoom${i}O3")) {
	$cmd = 'zoom'; $arg = 10;
    }
    elsif($q->param("zoom${i}I1")) {
	$cmd = 'zoom'; $arg = 1/1.5;
    }
    elsif($q->param("zoom${i}I2")) {
	$cmd = 'zoom'; $arg = 1/3;
    }
    elsif($q->param("zoom${i}I3")) {
	$cmd = 'zoom'; $arg = 1/10;
    }

    if($cmd eq 'move') {
	my $d = $arg * ($end-$start+1);
	$start = int($start + $d + .5);
	$end = int($end + $d + .5);
    }
    elsif($cmd eq 'zoom') {
	my $old_size = $end-$start+1;
	my $new_size = $old_size * $arg;
        my $d = ($new_size - $old_size) / 2;
	$start = int($start - $d + .5);
	$end = int($end + $d + .5);
    }

    $start = 1 if($start < 1);
    $end = 1 if($end < 1);

    if($end < $start) {
	print "Invalid coordinates: start greater than end!";
	return;
    }

    if($end >= $start + $max_size) {
	print "NOTE: region shrunk to maximum allowed size (",$max_size/1000000," Mb)";
	$end = $start + $max_size - 1;
    }
    
    return ($start,$end);
}


sub _new_tmp_file
{
    my $self = shift;
    my $dir = $self->tmp_dir;
    my ($fh, $fn) = File::Temp::tempfile(DIR => $dir,
					 SUFFIX => ".gfdb.png",
					 UNLINK => 0);
    #my ($html_path) = $fn =~ /(\/tmp\/\w+.gfdb.png$)/;
    my ($html_path) = $fn =~ /(\/\w+.gfdb.png$)/;
    $html_path = ($self->tmp_dir_url).$html_path;

    return ($fh, $html_path);
}


sub _clean_tmp_files
{
    my ($self) = @_;
    my @fn = glob(($self->tmp_dir)."/*gfdb.png");
    my $now = time;
    foreach my $fn (@fn) {
	my ($mtime) = (stat($fn))[9];
	my $age = $now - $mtime;
	unlink($fn) if($age > 15*60);  # remove if older than 15 minutes
    }
}


sub _javascript_header
{
    my ($self) = @_;
    my $javascript = "\tvar asm = new Object();\n";
    foreach my $asm (values %{$self->drawing_config->assemblies}) {
	my $asm_id = $asm->id;
	$javascript = $javascript."\tasm.$asm_id = new Object();\n";
	my @track_names = map { $_->{name} } @{$asm->tracks};
	my @track_status = map { $asm->enabled_tracks->{$_} ? 'true' : 'false' } @track_names;
	$javascript = $javascript."\tasm.$asm_id.tracks = [".join(',', map { "'$_'" } @track_names)."];\n";
	$javascript = $javascript."\tasm.$asm_id.tracksOn = [".join(',', @track_status)."];\n";
    }

    $javascript .= <<END_OF_JAVASCRIPT;
    asm.none = new Object();
    asm.none.tracks = ['No assembly selected'];
    asm.none.tracksOn = [0];
    function populateData( i, asmName ) { 
	select = window.document.form["tracks"+i];
	// Clear the old list of options
	select.options.length = 0; 
	// Create new options
        tracks = asm[asmName].tracks;
        tracksOn = asm[asmName].tracksOn;
	for( j = 0; j < tracks.length; j++ ) { 
	    select.options[j] = new Option( tracks[j], tracks[j], tracksOn[j] ); 
	} 
	checkbox = window.document.form["loc"+i+"active"];
	if(asmName == "none") {
	    checkbox.disabled = true;
	    if(checkbox.checked) {
		checkbox.checked = false;
		locActivityChange(i);
	    }
	}
	else {
	    checkbox.disabled = false;
	}
    } 
    function disableLocControls( i ) {
	window.document.form["loc"+i+"active"].checked = false;
	window.document.form["loc"+i].disabled = true;
	for( j = 1; j <= 3; j++) {
	    window.document.form["move"+i+"L"+j].disabled = true;
	    window.document.form["move"+i+"R"+j].disabled = true;
	    window.document.form["zoom"+i+"O"+j].disabled = true;
	    window.document.form["zoom"+i+"I"+j].disabled = true;
	}
    }
    function enableLocControls( i ) {
	window.document.form["loc"+i+"active"].checked = true;
	window.document.form["loc"+i].disabled = false;
	for( j = 1; j <= 3; j++) {
	    window.document.form["move"+i+"L"+j].disabled = false;
	    window.document.form["move"+i+"R"+j].disabled = false;
	    window.document.form["zoom"+i+"O"+j].disabled = false;
	    window.document.form["zoom"+i+"I"+j].disabled = false;
	}
    }
    function locActivityChange( i ) {
	if(i == 1) {
	    if(window.document.form.loc1active.checked) {	
		enableLocControls(1);
		disableLocControls(2);
	    }
	    else {
		disableLocControls(1);
		enableLocControls(2);
		disableLocControls(3);
	    }
	}
	else if(i == 2) {
	    if(window.document.form.loc2active.checked) {
		disableLocControls(1);
		enableLocControls(2);
		disableLocControls(3);
	    }
	    else {
		enableLocControls(1);
		disableLocControls(2);
	    }
	}
	else if(i == 3) {
	    if(window.document.form.loc3active.checked) {
		enableLocControls(1);
		disableLocControls(2);
		enableLocControls(3);
	    }
	    else {
		disableLocControls(3);
	    }
	}
    }
    function refreshIfEnterPressed(event) {
	if((event.which && event.which == 13) || (event.keyCode && event.keyCode == 13)) {
	    document.form.Refresh.click();
	    return false;
	} 
	else return true;
    }
END_OF_JAVASCRIPT

    return $javascript;
}



1;
