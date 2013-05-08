# AT::WWW::GenomeFeatureDB module
#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::WWW::GenomeFeatureDB - web-interface to genome feature database

=head1 SYNOPSIS

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::WWW::GenomeFeatureDB;

use strict;
use vars '@ISA';
use CGI;
use File::Temp;
use AT::Root;
use AT::DB::GenomeFeature;
use AT::DB::GenomeMapping;
use AT::DB::TranscriptSeq;
use AT::DB::GenomeAlignment;
use AT::DB::GenomeAssemblyNibs;
use AT::GFX::AlignedLociCL;
use AT::GFX::LocusCL;
use GeneLynx::MySQLdb;

@ISA = qw(AT::Root);

# Feature databases (should be in config file)
#my $FT_DBS = ['HS_NAT_01', 'HS_NAT_02', 'HS_NAT_03', 'HS_NAT_09', 'F3_NAT_04', 'F3_NAT_06', 'F3_NAT_07', 'F3_NAT_09'];
my $FT_DBS = ['HS_NAT_03','HS_NAT_09','F3_NAT_07','F3_NAT_09'];
my $FT_DB_LABELS = {#'HS_NAT_01' => "Human dataset hg17 v0.1 (HS_NAT_01)",
		    #'HS_NAT_02' => "Human dataset hg17 v0.2 (HS_NAT_02)",
		    'HS_NAT_03' => "Human dataset hg17 v0.4 (HS_NAT_03)",
		    'HS_NAT_09' => "Human dataset hg17 v0.5 (HS_NAT_09)",
		    #'F3_NAT_04' => "Mouse dataset mm5 v0.2 (F3_NAT_04)",
		    #'F3_NAT_06' => "Mouse dataset mm5 v0.3 (F3_NAT_06)",
		    'F3_NAT_07' => "Mouse dataset mm5 v0.4 (F3_NAT_07)",
		    'F3_NAT_09' => "Mouse dataset mm5 v0.5 (F3_NAT_09)" };
my $DEF_FT_DB1 = 'F3_NAT_09';
my $DEF_FT_DB2 = 'HS_NAT_09';
my %FT_DB_HOSTS  = {#'HS_NAT_01' => "nautilus",
		    #'HS_NAT_02' => "nautilus",
		    'HS_NAT_03' => "angband",
		    'HS_NAT_09' => "angband",
		    #'F3_NAT_04' => "nautilus",
		    #'F3_NAT_06' => "nautilus",
		    'F3_NAT_07' => "angband",
		    'F3_NAT_09' => "angband"};


# Tracks user can choose to display
my $MAIN_TRACKS = ['Chains','TUs','TFs','mRNAs','ESTs'];
my $MAIN_TRACK_LABELS = { };            # Displayed track names (optional)
my $DEF_MAIN_TRACKS = ['Chains','TUs','mRNAs'];  # Tracks displayed by default
my $XTRA_TRACKS = ['Assembly', 'Gaps', 'CpG islands'];
my $XTRA_TRACK_LABELS = { }; 
my $DEF_XTRA_TRACKS = [];               # Tracks displayed by default
my $XTRA_TRACK_COMMENT = { 'Assembly' => 'brown = finished; cyan = whole genome shotgun; gray = other (e.g. draft)'
			 };

sub new
{
    my ($caller, %args) = @_;
    my $config = $args{config} || die "No config arg";
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


sub nib_dir
{
    my ($self, $asm_name) = @_;
    my $base = $self->{config}->{nib_dir} || die "No nib_dir configured\n";
    return "$base/$asm_name/nib";
}


sub handle_request
{
    my ($self, $q, $dbuser, $dbpass) = @_;
    # dbuser and dbpass should be taken from $self->config

    # html header
    print $q->header,
        $q->start_html("Comparative Genome Browser"),
        $q->h1("Comparative Genome Browser");

    # get params
    my $dbname1 = $q->param('db1');
    my $dbname2 = $q->param('db2');
    my $loc = $q->param('loc');
    my $tracks1 = [$q->param('mtracks1'),$q->param('xtracks1')];
    my $tracks2 = [$q->param('mtracks2'),$q->param('xtracks2')];
    my ($button, $button_value) = _get_button_and_value($q);
    
    my ($chr, $start, $end, $drawer);
    if($dbname1 and $loc) {
        # prepare to draw image: connect to db and construct drawing object
        if($dbname2 eq '(none)') {
            my $dbs = $self->_db_connect1($dbname1, $dbuser, $dbpass);
            $drawer = AT::GFX::LocusCL->new(%$dbs,
    					 organism => ($dbname1 =~ /^F3/) ? 'mouse' : 'human',
    					 tracks => $tracks1);
            ($chr, $start, $end) = $dbs->{feature_db}->id_string_to_region($loc);
	    print $dbs->{feature_db}->errstr unless($chr);
        }
        else {
            my $dbs = $self->_db_connect2($dbname1, $dbname2, $dbuser, $dbpass);
            $drawer = AT::GFX::AlignedLociCL->new(%$dbs,
    					 organism1 => ($dbname1 =~ /^F3/) ? 'mouse' : 'human',
    					 organism2 => ($dbname2 =~ /^F3/) ? 'mouse' : 'human',
    					 tracks1 => $tracks1,
    					 tracks2 => $tracks2);
            ($chr, $start, $end) = $dbs->{feature_db1}->id_string_to_region($loc);
	    print $dbs->{feature_db1}->errstr unless($chr);
        }
        if($chr) {
	    ($start,$end) = _calc_range($start,$end,$button,$button_value);
	    if($start) {
		$q->param('loc', "$chr:$start-$end");
	    }
	    else {
		$drawer = undef;
	    }
	}
    }
 
    # input form
    print $q->start_form(-method=>'GET');
    print $q->start_table();
    print $q->Tr({-align => 'left'},$q->th(['Location','Database']),
                                    $q->th({-colspan => scalar(@$MAIN_TRACKS)},
                                           'Tracks ('.($q->a({-href=>'#extra_tracks'},'more')).')')
                                    
                );
    print $q->Tr($q->td([$q->textfield('loc','chr1:100,000-200,000',30,50),
			 $q->popup_menu('db1',
			    		 $FT_DBS,
					 $DEF_FT_DB1,
					 $FT_DB_LABELS),
			 $q->checkbox_group('mtracks1', $MAIN_TRACKS, $DEF_MAIN_TRACKS,
                                            labels => $MAIN_TRACK_LABELS)
			]));
    print $q->Tr($q->td({-align =>'right'}, 'versus'),
		 $q->td([$q->popup_menu('db2',
			    		 ['(none)',@$FT_DBS],
					 $DEF_FT_DB2,
					 $FT_DB_LABELS),
		         $q->checkbox_group('mtracks2',$MAIN_TRACKS, $DEF_MAIN_TRACKS,
                                            labels => $MAIN_TRACK_LABELS)
			 ])
		 );
    print $q->end_table;
    print $q->start_table;
    print $q->Tr($q->td([$q->submit('show'),
			 $q->pre(" "),
			 "move",
			 $q->submit('move','<<<').
			 $q->submit('move','<<').
			 $q->submit('move','<').
			 $q->submit('move','>').
			 $q->submit('move','>>').
			 $q->submit('move','>>>'),
			 $q->pre(" "),
			 "zoom in",
			 $q->submit('zoom in','1.5x').
			 $q->submit('zoom in','3x').
			 $q->submit('zoom in','10x'),
                         $q->pre(" "),
			 "zoom out",
			 $q->submit('zoom out','1.5x').
			 $q->submit('zoom out','3x').
			 $q->submit('zoom out','10x'),
			 $q->pre(" "),
			 $q->a({-href => '/f3/gfdb/help.html'}, "Help")
			])
		 );
    print $q->end_table;
    print $q->hr;
    
    # draw image and output   
    if($drawer) {
        print $q->h5("View of $loc");
	if($drawer->draw_locus($chr,$start,$end)) {
	    my ($img_fh, $img_fn) = $self->new_tmp_file();
	    print $img_fh $drawer->gd_image->png;
	    close $img_fh;
	    print $q->img({-src => $img_fn});
	}
	else {
	    print $q->p($drawer->errstr);
	}
    }

    # more tracks
    print $q->hr;
    print $q->a({-name => 'extra_tracks'});
    print $q->h4('Extra tracks');
    print $q->start_table();
    print $q->Tr($q->th({-align => 'left'}, ['Track', 'Primary panel&nbsp;&nbsp;', 'Secondary panels', '']));
    my %xtrack_blank_labels = map { $_ => '' } @$XTRA_TRACKS;
    my @xtracks1_checkbox_group = $q->checkbox_group('xtracks1',$XTRA_TRACKS,$DEF_XTRA_TRACKS,
                                                     labels => \%xtrack_blank_labels);
    my @xtracks2_checkbox_group = $q->checkbox_group('xtracks2',$XTRA_TRACKS,$DEF_XTRA_TRACKS,
                                                     labels => \%xtrack_blank_labels);
    foreach my $track_name (@$XTRA_TRACKS) {       
        my $track_label = $XTRA_TRACK_LABELS->{$track_name} || $track_name;
	my $track_comment = $XTRA_TRACK_COMMENT->{$track_name} || '';
        print $q->Tr($q->td($track_label),
                     $q->td({-align => 'center'},
                            [shift @xtracks1_checkbox_group,
                             shift @xtracks2_checkbox_group
                            ]),
		     $q->td($track_comment)
                     );
    }
    print $q->end_table;   

    # html footer
    print $q->end_form;
    print $q->end_html;

    $self->clean_tmp_files();
}


sub _get_button_and_value
{
    my $q = shift;
    my @buttons = ("show", "zoom in", "zoom out", "move");
    foreach my $button (@buttons) {
	my $value = $q->param($button);
	return ($button, $value) if($value);
    }
    return ("","");
}


sub _calc_range
{
    my ($start, $end, $cmd, $cmd_arg) = @_;
    
    my $max_size = 5000000;
    
    if($cmd eq 'move') {
	my $speed;
	my $arg_length = length($cmd_arg);
	if($arg_length == 1) { $speed = 0.1 }
	elsif($arg_length == 2) { $speed = 0.5 }
	elsif($arg_length == 3) { $speed = 0.95 }
	else { $speed = 0; } 
	my $d = $speed * ($end-$start+1);
	$d = -$d if(substr($cmd_arg,0,1) eq '<');
	$start = int($start + $d + .5);
	$end = int($end + $d + .5);
    }
    elsif($cmd eq 'zoom in' or $cmd eq 'zoom out') {
	unless(chop $cmd_arg eq 'x') { die "Zoom factor does not end in 'x'" }
	my $old_size = $end-$start+1;
	my $new_size = ($cmd eq 'zoom in') ?  $old_size / $cmd_arg : $old_size * $cmd_arg;
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


sub _db_connect1
{
    my ($self, $dbname, $DBUSER, $DBPASS) = @_;

    my $dbhost = $FT_DB_HOSTS{$dbname};

    # This should be taken from $self->config; as should username and password
    my ($seq_db_name, $map_db_name, $gl_db_name) =
	( $dbname =~ /^F3/) ?
	('F3_TRANSCRIPT_SEQ', 'AT_MM_MAY04', 'GLMOUSE_1_99') :
	('TRANSCRIPT_SEQ', 'AT_HS_MAY04', 'GLHUMAN_1_99') ;

    my $seqdb = AT::DB::TranscriptSeq->connect(-dbname => $seq_db_name,
					       -dbhost => 'nautilus.cgb.ki.se',
					       -dbuser => $DBUSER,
					       -dbpass => $DBPASS);
    my $ftdb = AT::DB::GenomeFeature->connect(-dbname => $dbname,
					  -dbhost => $dbhost,
					  -dbuser => $DBUSER,
					  -dbpass => $DBPASS);
    my $mapdb = AT::DB::GenomeMapping->connect(-dbname => $map_db_name,
					  -dbhost => 'nautilus.cgb.ki.se',
					  -dbuser => $DBUSER,
					  -dbpass => $DBPASS);
    my $gldb= GeneLynx::MySQLdb->connect(-dbname => $gl_db_name,
				       -dbhost => 'nautilus.cgb.ki.se',
				       -dbuser => $DBUSER,
				       -dbpass => $DBPASS);

    my %dbs =  ( trseq_db => $seqdb,
		 feature_db => $ftdb,
                 mapping_db => $mapdb,
		 genelynx_db => $gldb );
    return \%dbs;
}


sub _db_connect2
{
    my ($self, $dbname1, $dbname2, $DBUSER, $DBPASS) = @_;

    my $dbhost1 = $FT_DB_HOSTS{$dbname1};
    my $dbhost2 = $FT_DB_HOSTS{$dbname2};

    # This should be taken from $self->config; as should username and password
    my ($seq_db_name1, $seq_db_name2, $gl_db_name, $map_db_name1, $map_db_name2, $asm_db_name1, $asm_db_name2);
    if($dbname1 =~ /^F3/) {
	unless($dbname2 =~ /^HS/) { print "Cannot compare $dbname1 and $dbname2"; die; }
	($seq_db_name1, $seq_db_name2, $gl_db_name, $map_db_name1, $map_db_name2, $asm_db_name1, $asm_db_name2) =
	('F3_TRANSCRIPT_SEQ', 'TRANSCRIPT_SEQ', 'GLMOUSE_1_99', 'AT_MM_MAY04', 'AT_HS_MAY04', 'MM_MAY04', 'HS_MAY04');
    }
    else {
	unless($dbname2 =~ /^F3/) { print "Cannnot compare $dbname1 and $dbname2"; die; }
	($seq_db_name1, $seq_db_name2, $gl_db_name, $map_db_name1, $map_db_name2, $asm_db_name1, $asm_db_name2) =
	('TRANSCRIPT_SEQ', 'F3_TRANSCRIPT_SEQ', 'GLHUMAN_1_99', 'AT_HS_MAY04', 'AT_MM_MAY04', 'HS_MAY04', 'MM_MAY04') ;
    }
    
    my $ftdb1 = AT::DB::GenomeFeature->connect(-dbname => $dbname1,
					  -dbhost => $dbhost1,
					  -dbuser => $DBUSER,
					  -dbpass => $DBPASS);
    my $ftdb2 = AT::DB::GenomeFeature->connect(-dbname => $dbname2,
					  -dbhost => $dbhost2,
					  -dbuser => $DBUSER,
					  -dbpass => $DBPASS);

    my $mapdb1 = AT::DB::GenomeMapping->connect(-dbname => $map_db_name1,
					  -dbhost => 'nautilus.cgb.ki.se',
					  -dbuser => "engstrom",
					  -dbpass => "cog117");

    my $mapdb2 = AT::DB::GenomeMapping->connect(-dbname => $map_db_name2,
					  -dbhost => 'nautilus.cgb.ki.se',
					  -dbuser => "engstrom",
					  -dbpass => "cog117");

    my $alndb = AT::DB::GenomeAlignment->connect(-dbname => $map_db_name1,
					  -dbhost => 'nautilus.cgb.ki.se',
					  -dbuser => "engstrom",
					  -dbpass => "cog117");

    my $seqdb1 = AT::DB::TranscriptSeq->connect(-dbname => $seq_db_name1,
					       -dbhost => 'nautilus.cgb.ki.se',
					       -dbuser => $DBUSER,
					       -dbpass => $DBPASS);
    my $seqdb2 = AT::DB::TranscriptSeq->connect(-dbname => $seq_db_name2,
					       -dbhost => 'nautilus.cgb.ki.se',
					       -dbuser => $DBUSER,
					       -dbpass => $DBPASS);


    my $asmdb1 = AT::DB::GenomeAssemblyNibs->new(assembly_name => $asm_db_name1,
	dir => $self->nib_dir($asm_db_name1));

    my $asmdb2 = AT::DB::GenomeAssemblyNibs->new(assembly_name => $asm_db_name2,
	dir => $self->nib_dir($asm_db_name2));

    my %dbs =  ( feature_db1 => $ftdb1,
		 feature_db2 => $ftdb2,
                 mapping_db1 => $mapdb1,
                 mapping_db2 => $mapdb2,
		 alignment_db1 => $alndb,
		 assembly_db1 => $asmdb1,
		 assembly_db2 => $asmdb2,
		 trseq_db1 => $seqdb1,
		 trseq_db2 => $seqdb2,
		 );
    return \%dbs;
}


sub new_tmp_file
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


sub clean_tmp_files
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



# Quick-docs:
#
# my @param_keys = $q->param     get parameter names
# my @values = $q->param($key);  get param value(s)
# $q->delete(@keys)              delete parameters
# $q->delete_all();              delete all parameters
# 
# $q->self_url();                url of script
#
# $q->hr;                        horiz row


1;
