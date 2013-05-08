# AT::WWW::GenomeFeatureDB module
#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::WWW::ChainViewer - web-interface to chain images

=head1 SYNOPSIS

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::WWW::ChainViewer;

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
use AT::GFX::AlignedLoci;
use AT::GFX::Locus;
use GeneLynx::MySQLdb;

@ISA = qw(AT::Root);

# Feature databases (should be in config file)
my $FT_DB1 = 'HS_NAT_09';
my $FT_DB2 = 'F3_NAT_09';
my $DBHOST = 'localhost';
my $UCSC_DB1 = 'hg17';
my $UCSC_DB2 = 'mm5';
my $ORGANISM1 = 'human';
my $ORGANISM2 = 'mouse';

my $TABLE_DESC = "Chain images show locations in the genome of mRNA (cDNA) mappings and EST mappings, 
TUs inferred from those mappings and the span of chains formed by the TUs. Each of these types of 
features is show in a separate track. Locations of CpG islands, highlighting potential promoter
regions, are also indicated in a separate track. Mappings are labeled with DDBJ/EMBL/GenBank
 transcript sequence accession numbers, or with RIKEN clone ids for FANTOM3 sequences not assigned
 accession numbers. RIKEN EST mappings labeled with RIKEN clone ids are suffixed by F, R, or F+R, indicating EST read direction (forward, reverse or both)."; 

sub new
{
    my ($caller, %args) = @_;
    my $config = $args{config} || die "No config arg";
    my $self = bless { config => $config }, ref $caller || $caller;
    return $self;   
}


sub chain_img_dir
{
    shift->{config}->{chain_img_dir} || die "No chain_img_dir configured\n";
}



sub handle_request
{
    my ($self, $q, $dbuser, $dbpass) = @_;
    # dbuser and dbpass should be taken from $self->config

    # html header
    print $q->header,
        $q->start_html("Chains in Human and Mouse Genomes"),
        $q->h3("Chains in Human and Mouse Genomes");

    # get params
    my $org = $q->param('org');
    my $chr = $q->param('chr');
    my $order = $q->param('order');
    
    # input form
    print $q->start_form(-method=>'GET');
    print $q->start_table();
    print $q->Tr({-align => 'left'},$q->th(['Genome','Chr','Order','']));
    print $q->Tr($q->td([$q->popup_menu('org', ['Human','Mouse'],'Human'),
			 $q->popup_menu('chr', ['All',1..22,'X','Y']),
			 $q->popup_menu('order', ['Location','Size','ID'],'Location'),
			 $q->submit('List')
			 ]));
    print $q->end_table;
    
    # draw image and output   
    if($org) {
	print $q->hr;
	$org = lc $org;
	my $subhead = $chr eq 'All' ?
	    "All $org chains" : "Chains on $org Chromosome $chr";
        print $q->h4($subhead);
	my ($dbname, $ucsc_dbname) = $org eq $ORGANISM1 ? ($FT_DB1, $UCSC_DB1) : ($FT_DB2, $UCSC_DB2);
	$self->list_chains($q, $dbname, $dbuser, $dbpass, $org, $ucsc_dbname, $chr, $order);
    }

    # more tracks
    print $q->hr;

    # html footer
    print $q->end_form;
    print $q->end_html;

}


sub list_chains
{
    my ($self, $q, $dbname, $dbuser, $dbpass, $org, $ucsc_db, $chr_req, $order) = @_;

    my $db = AT::DB::GenomeFeature->connect(-dbname => $dbname,
					    -dbhost => $DBHOST,
					    -dbuser => $dbuser,
					    -dbpass => $dbpass);

    my $query = "select distinct chain_id, chr, start, end, size from chain_aspair_bdpair ";
    if($chr_req ne 'All') {
	$query .= " where chr = 'chr$chr_req' ";
    }
    my $sth = $db->dbh->prepare($query);
    $sth->execute();

    my @chains;
    while(my ($id,$chr,$start,$end,$size) = $sth->fetchrow_array) {
	$chr =~ s/^chr//;
	push @chains, { id => $id,
			chr => $chr,
			start => $start,
			end => $end,
			size => $size };
    }

    if($order eq 'Size') {
	@chains = sort { $b->{size} <=> $a->{size} or order_by_loc() } @chains;
    }
    elsif($order eq 'Location') {
	@chains = sort order_by_loc @chains;
    }
    elsif($order eq 'ID') {
	@chains = sort { $a->{id} <=> $b->{id} } @chains;
    }

    print $q->start_table();
    print $q->Tr({-align => 'center', bgcolor=> 'lightblue'},
		 $q->th(['ID','Chr','Start','End','Size (TUs)']),
		 $q->th({-colspan => 2}, 'Links')
		 );
    my $chain_img_dir = $self->chain_img_dir;

    my $flip = 0;

    foreach my $c (@chains) {
	my $id = $c->{id};
	my $chr = $c->{chr};
	my $start = $c->{start};
	my $end = $c->{end};
	my $size = $c->{size};
	my $ucsc_loc = 'chr'.$chr.':'.($start-5000).'-'.($end+5000);
	my $img_link = $q->a({-href=>"$chain_img_dir/$org/chain_$id.png"},'Image');
	my $ucsc_link = $q->a({-href=>" http://genome.ucsc.edu/cgi-bin/hgTracks?org=$org&db=$ucsc_db&position=$ucsc_loc"},
			      'UCSC');
	my $bg = $flip ? 'lightgray' : 'lightgray';
	print $q->Tr({-align => 'right', -bgcolor => $bg},
		     $q->td([$id,
			     $chr,
			     $start,
			     $end,
			     $size,
			     $img_link,
			     $ucsc_link
			     ]));
	$flip = 1-$flip;
    }
    print $q->end_table();

    print $q->p($TABLE_DESC);

}


{
    my %chr_order;
    foreach my $chr (1..22) {
	$chr_order{$chr} = $chr;
    }
    $chr_order{'X'} = 23;
    $chr_order{'Y'} = 24;
 
    sub order_by_loc
    {
	return ( $chr_order{$a->{chr}} <=> $chr_order{$b->{chr}} or $a->{start} <=> $b->{start} );
    }

}

1;
