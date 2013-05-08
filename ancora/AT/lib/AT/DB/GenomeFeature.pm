# AT::DB::GenomeFeature module
#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::DB::GenomeFeature - interface to MySQL relational database of genome
features


=head1 SYNOPSIS

=over 4

=item * creating a database object by connecting to an existing AT-type
database

my $db = AT::DB::GenomeMapping->connect(
    -dbname => "AT",
    -dbhost => "myhostn.mydomain",
    -dbuser => "myusername",
    -dbpass => "mypassword");

=back

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::DB::GenomeFeature;

use strict;
use vars '@ISA';
use Bio::Seq;
use Bio::LocatableSeq;
use Carp;
use AT::DB::Binner;
use AT::DB::MySQLdb;
use AT::DB::ResultStream;
use AT::FT::Gene::Dynamic;
use AT::FT::PEP;
use AT::FT::GFMapping;
use AT::Prediction::AntisensePair;

@ISA = qw(AT::DB::MySQLdb AT::DB::Binner);


sub _init
{
    my ($self, %args) = @_;
    $self->{assembly_db} = $args{assembly_db};
    $self->{mapping_db} = $args{mapping_db};
     $self->{get_primary_mappings} = $args{get_primary_mappings} || 0;
}


sub errstr
{
    return shift->{_errstr} || '';
}

#sub _assembly_db
#{
#    my ($self) = @_;
#    return $self->{_auxiliary_dbs}[0] || croak "no assembly db";
#}

#sub _mapping_db
#{
#    my ($self) = @_;
#    return $self->{_auxiliary_dbs}[1] || croak "no mapping db";
#}


# Antisense pair methods

sub store_antisense_pair
{
    my ($self, $pair) = @_;

    #croak "No pair_set_nr argument " unless(defined $pair_set_nr);

    my $gene1 = $pair->member1;
    my $gene2 = $pair->member2;

    if($pair->pair_id) {
	warn "Will not update database entry for pair ", $pair->pair_id, "\n";
	return $pair->pair_id;
    }
    #my $pair_id = $pair->pair_id or croak("No pair_id for AntisensePair ",$gene1->loc_str,"/",gene2->loc_str);

    # store genes if not done already
    $self->store_gene($gene1) unless($gene1->gene_id);
    $self->store_gene($gene2) unless($gene2->gene_id);

    # get table id
    #my $table_id = ($gene1->assembly).'_ps'.$pair_set_nr;

    my $query1 = "INSERT INTO aspair (gene_id1, gene_id2,
				      category, splicing, rel_ori,
				      exon_ol_sizes, exon_ol_starts,
                                      total_exon_ol, total_repeat_exon_ol,
				      rep_acc1, rep_acc2,
				      mRNA_accs1, mRNA_accs2, EST_accs1, EST_accs2,
                                      notes)
				     VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)";
    my $sth1 = $self->dbh->prepare_cached($query1);
    my ($g1_as_rep, $g2_as_rep) = $pair->representative_ids;
    my ($ol_mRNA_ids1, $ol_EST_ids1) = $pair->overlapper_ids(1);
    my ($ol_mRNA_ids2, $ol_EST_ids2) = $pair->overlapper_ids(2);
    my $repeat_ol = $pair->total_repeat_exon_ol;
    $repeat_ol = '\N' unless defined($repeat_ol);

    $sth1->execute($gene1->gene_id,
		  $gene2->gene_id,
		  $pair->category,
		  $pair->SU_category,
		  $pair->relative_orientation,
		  join(",", $pair->exon_overlaps),
		  join(",", $pair->exon_overlap_starts),
		  $pair->exon_overlap_sum,
		  $repeat_ol,
		  $g1_as_rep,
		  $g2_as_rep,
		  join(',',@$ol_mRNA_ids1),
		  join(',',@$ol_mRNA_ids2),
		  join(',',@$ol_EST_ids1),
		  join(',',@$ol_EST_ids2),
		  join('; ',@{$pair->notes}));
    $sth1->finish;

    my $pair_id = $self->_last_insert_id;
    $pair->pair_id($pair_id);

    if($pair->category < 3) {
	my $query2 = "INSERT INTO aspair_exon_ol_support (aspair_id, seq_types, splicing,
                                                          nr_mRNAs1, nr_mRNAs2,
                                                          nr_ESTs1, nr_ESTs2)
                                                          VALUES(?,?,?,?,?,?,?)";
	my $sth2 = $self->dbh->prepare_cached($query2);
	$sth2->execute($pair_id, $pair->exon_overlap_support);
	$sth2->finish;
    }

    return $pair_id;
}


sub get_antisense_pair_by_id
{
    my ($self, $pair_id) = @_;

    # Get two geneids from aspair table
    #my $table_id = "${assembly}_ps${pair_set_nr}";
    my $query = "SELECT gene_id1, gene_id2 FROM aspair WHERE aspair_id = ?";
    my $sth = $self->dbh->prepare_cached($query);
    $sth->execute($pair_id);
    my ($gene_id1, $gene_id2) = $sth->fetchrow_array();
    $sth->finish();
    return unless($gene_id2);

    # get parameters for this antisense pair set
#    my ($gene_set_nr, $stringency, $score_method) =
#	$self->_get_aspair_set_params($assembly, $pair_set_nr);
die "get_antisense_pair_by_id needs updating"; my($stringency, $score_method);

    # retrieve the two genes
    my $gene1 = $self->get_gene_by_id($gene_id1, $stringency, $score_method);
    my $gene2 = $self->get_gene_by_id($gene_id2, $stringency, $score_method);
    unless($gene1 and $gene2) {
	croak "Missing gene(s) for antisense pair $pair_id";
    }

    # return aspair object
    return AT::Prediction::AntisensePair->new(pair_id => $pair_id,
					      member1 => $gene1,
					      member2 => $gene2);
}


sub _get_aspair_set_params
{
    my ($self, $assembly, $pair_set_nr) = @_;
    my $key = "${assembly}_ps${pair_set_nr}";
    my $params_ref = $self->{_aspair_set_params}{$key};
    unless($params_ref) {
	my $query = "SELECT gene_set, stringency, score_method FROM aspair_param WHERE assembly = ? AND aspair_set = ?";
	my $sth = $self->dbh->prepare_cached($query);
	$sth->execute($assembly, $pair_set_nr);
	my @params = $sth->fetchrow_array();
	$sth->finish();
	unless(@params) {
	    croak "No params for aspair_set $key";
        }
	$self->{_aspar_set_params}{$key} = $params_ref = \@params;
    }
    return @$params_ref;
}


# Gene methods


sub store_gene  {
    my ($self, $gene) = @_;
    my ($query, $sth);

    # Check arguments
    unless ($gene->isa("AT::FT::Gene::Dynamic")) { 
        croak "Wrong type of gene object - not dynamic gene";
    }
    #$gene_set_nr = 1 unless(defined $gene_set_nr);

    if($gene->gene_id) {
	warn "Will not update database entry for gene ", $gene->gene_id, "\n";
	return $gene->gene_id;
    }

    print STDERR "storing gene ", $gene->loc_str, "\n";

    my $offset = $gene->start - 1;
    my @exons = $gene->strand == 1 ? $gene->exon_list : reverse $gene->exon_list;
    my $exon_string = join(";", map { ($_->abs_start - $offset).",".($_->abs_end - $offset) } @exons);

    # get table id
    #my $table_id = ($gene->assembly).'_gs'.$gene_set_nr;

    # Insert data into gene table
    $gene->gene_id(0);   # Make sure auto_increment is used to set id
    $query = "INSERT INTO gene (chr, min_start, max_end, start, end, strand, rep_acc, nr_ests, score, is_spliced, exons) VALUES(?,?,?,?,?,?,?,?,?,?,?)";
    $sth = $self->dbh->prepare_cached($query);
    $sth->execute($gene->chr,
		  $gene->min_start,
		  $gene->max_end,
		  $gene->start,
		  $gene->end,
		  $gene->strand,
		  $gene->representative_mRNA_acc || '',
		  $gene->nr_EST_gfmappings,
		  $gene->gfmapping_score_sum,
		  $gene->is_spliced,
		  $exon_string);
    $sth->finish;
    my $gene_id = $self->_last_insert_id;
    $gene->gene_id($gene_id);

#    # Insert data into pep table
#    my $gene_offset = $gene->min_start - 1;
#    foreach my $pep ($gene->pep_list) {
#	$pep->pep_id(0);
#        $query = "INSERT INTO pep_${table_id} (gene_id, rel_start, rel_end, score1, score2, open_score, close_score) VALUES(?,?,?,?,?,?,?)";
#        $sth = $self->dbh->prepare_cached($query);
#        $sth->execute($gene_id,
#		      $pep->start - $gene_offset,
#		      $pep->end - $gene_offset,
#		      $pep->score1,
#		      $pep->score2,
#		      $pep->open_score,
#		      $pep->close_score);
#        my $pep_id = $self->_last_insert_id;
#        $pep->pep_id($pep_id);
#    }

    # Store gfmappings
    my @gfmap_ids;
    foreach my $gfmap ($gene->gfmapping_list) {
	push @gfmap_ids, $self->store_gfmapping($gfmap);
    }

    # Insert data into gene_gfmap table
    $query = "INSERT INTO gene_gfmap (gene_id, gfmap_id) VALUES(?,?)";
    $sth = $self->dbh->prepare_cached($query);
    foreach my $gfmap_id (@gfmap_ids) {
        $sth->execute($gene_id, $gfmap_id);
    }
    $sth->finish();
    
    return $gene_id;
}


sub get_gene_by_id
{
    my ($self, $gene_id, $stringency, $score_method, %opt) = @_;
    my ($query, $sth);

    #my $table_id = "${assembly}_gs${gene_set_nr}";

    # get gene data
    $query = "SELECT chr,min_start,max_end,strand FROM gene WHERE gene_id = ?";
    $sth = $self->dbh->prepare_cached($query);
    $sth->execute($gene_id);
    my $gene_row = $sth->fetchrow_arrayref();
    $sth->finish();
    return unless($gene_row);
    my ($chr, $gene_start, $gene_end, $strand) = @$gene_row;

    # get gfmapping data
    my $gfmappings = $self->get_gfmappings_by_gene_id($gene_id)
	unless($opt{no_gfmappings});

    # get pep data and create pep objects
    $query = "SELECT * FROM pep WHERE gene_id = ?";
    $sth = $self->dbh->prepare_cached($query);
    $sth->execute($gene_id);
    my @peps;
    while(my $pep_row = $sth->fetchrow_arrayref()) {
	push @peps, AT::FT::PEP->new(#pep_id => $pep_row->[0],
				     start => $gene_start + $pep_row->[2] - 1,
				     end => $gene_start + $pep_row->[3] - 1,
				     score1 => $pep_row->[4],
				     score2 => $pep_row->[5],
				     open_score => $pep_row->[6],
				     close_score => $pep_row->[7]);
    }
    $sth->finish();

# store splices in table splice_${assembly}_gs${gene_set_nr}
# fields: left_pep_id, right_pep_id, score
# ie very easy!
# the tricky part is how to represent pep graph in memory with splices still there
    
    # get genome seq
    my $seq = $self->assembly_db->get_genome_seq(chr => $chr,
				     	         start => $gene_start,
						 end => $gene_end,
						 strand => $strand);
    
    # return gene object
    return AT::FT::Gene::Dynamic->new
	(seq => $seq,
	 peps => \@peps,
	 strand => $strand,
	 stringency => $stringency,
	 score_method => $score_method,
	 gfmapping_list => $gfmappings);    

}


# GFMapping methods




sub store_gfmapping
{
    my ($self, $m) = @_;
    my ($query, $sth);

    if($m->gfmapping_id) {
	#warn "Will not update database entry for GFMapping ", $m->gfmapping_id, "\n";
	return $m->gfmapping_id;
    }

    my $assembly = $m->assembly;
    my $table = $self->{gfmapping_table};

    my $hsp_string;
    my $offset = $m->start - 1;
    foreach my $hsp ($m->HSP_list) {
	my $gap = $hsp->right_gap;
	my $gap_type;
	if($gap) {
	    if($gap->is_broken) { $gap_type = "B"; }
	    elsif($gap->is_intron) { $gap_type = "I"; }
	    else { $gap_type = "G"; }
	}
	$hsp_string .= join(',',($hsp->start - $offset, $hsp->end - $offset, $gap ? ($gap_type,$gap->jnc_str) : ())).';';
    }
    
    $query = "INSERT INTO gfmap (cluster_id, chr, start, end, strand, trust_start, trust_end, query_bp_past_left_border,
				 query_bp_past_right_border, starts_near_T_stretch, ends_near_A_stretch, internally_primed,
				 max_qInsert_size, best_for_query, has_3p_EST, has_5p_EST, nr_introns, score1, score2, hsps)
				VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)";
    $sth = $self->dbh->prepare_cached($query);
    $sth->execute($m->cluster_id,
		  $m->chr,
		  $m->start,
		  $m->end,
		  $m->strand,
		  $m->trust_start,
		  $m->trust_end,
		  $m->query_bp_past_left_border,
		  $m->query_bp_past_right_border,
		  $m->starts_near_T_stretch,
		  $m->ends_near_A_stretch,
		  $m->internally_primed,
		  $m->max_qInsert,
		  $m->best_for_query,
		  $m->has_3p_EST,
		  $m->has_5p_EST,
		  $m->nr_introns,
		  $m->score1,
		  $m->score2,
		  $hsp_string);
    $sth->finish;
    my $gfmap_id = $self->_last_insert_id;
    $m->gfmapping_id($gfmap_id);

    $query = "INSERT INTO gfmap_primap (gfmap_id, primap_id, acc, version, type, reversed) VALUES (?,?,?,?,?,?)";
    $sth = $self->dbh->prepare_cached($query);
    foreach my $pri_map ($m->primary_mapping_list) {
	my $reversed = $m->strand == $pri_map->strand_numeric ? 0 : 1;
        $sth->execute($gfmap_id, $pri_map->mapping_id, $pri_map->qName, $pri_map->mRNAInfo->version, $pri_map->qType, $reversed);
    }
    $sth->finish;
    
    return $gfmap_id;
}


sub get_gfmappings_by_gene_id
{
    my ($self, $assembly, $gene_set_nr, $gene_id) = @_;
    my ($query, $sth);

    #my $table_id = "${assembly}_gs${gene_set_nr}";
    $query = "SELECT b.* FROM gene_gfmap a, gfmap_$assembly b WHERE a.gene_id = ? and a.gfmap_id = b.gfmap_id";
    $sth = $self->dbh->prepare_cached($query);
    $sth->execute($gene_id);
    my @gfmappings;
    while(my $row = $sth->fetchrow_hashref()) {
	push @gfmappings, $self->_create_gfmapping_from_row($assembly, $row);
    }
    $sth->finish();

    return \@gfmappings;
}


sub get_gfmapping_by_id
{
    my ($self, $id) = @_;
    my ($query, $sth);

    $query = "SELECT * FROM gfmap a WHERE a.gfmap_id = ?";
    $sth = $self->dbh->prepare_cached($query);
    $sth->execute($id);
    my $row = $sth->fetchrow_hashref();
    return unless($row);
    my $gfmapping = $self->_create_gfmapping_from_row($row);
    $sth->finish();

    return $gfmapping;
}


sub get_gfmapping_by_primary_mapping_id
{
    my ($self, $id) = @_;
    my ($query, $sth);

    $query = "SELECT b.* FROM gfmap_primap a, gfmap b WHERE a.primap_id = ? and a.gfmap_id = b.gfmap_id";
    $sth = $self->dbh->prepare_cached($query);
    $sth->execute($id);
    my $row = $sth->fetchrow_hashref();
    return unless($row);
    my $gfmapping = $self->_create_gfmapping_from_row($row);
    $sth->finish();

    return $gfmapping;
}


sub get_one_gfmapping_by_acc
{
    my ($self, $id) = @_;
    my ($query, $sth);

    $query = "SELECT b.* FROM gfmap_primap a, gfmap b WHERE a.acc = ? and a.gfmap_id = b.gfmap_id order by rand() limit 1";
    $sth = $self->dbh->prepare_cached($query);
    $sth->execute($id);
    my $row = $sth->fetchrow_hashref();
    return unless($row);
    my $gfmapping = $self->_create_gfmapping_from_row($row);
    $sth->finish();

    return $gfmapping;
}


sub get_gfmappings_by_cluster_id
{
    my ($self, $id, $table) = @_;
    my ($query, $sth);
    $table = "gfmap" unless ($table);
    $query = "SELECT * FROM $table a WHERE a.cluster_id = ?";
    $sth = $self->dbh->prepare_cached($query);
    $sth->execute($id);
    my @gfmappings;
    while(my $row = $sth->fetchrow_hashref()) {
        push @gfmappings, $self->_create_gfmapping_from_row($row);
    }
    $sth->finish();

    return \@gfmappings;
}



sub get_gfmappings_in_region_as_cluster_stream
{
    my ($self, $chr, $start, $end, $table) = @_;

    $start = 1 unless ($start);
    $end = $self->assembly_db->get_chr_size($chr) unless($end);

    my $cluster_query = "SELECT cluster_id FROM gfmap_cluster
			 WHERE chr = ? AND end >= ? AND start <= ?
			 ORDER BY chr, start";
    my $cluster_sth = $self->dbh->prepare($cluster_query);
    $cluster_sth->execute($chr, $start, $end);

    my $next_method = sub {
	my ($cluster_id) = $cluster_sth->fetchrow_array();
	return unless($cluster_id);
	$self->get_gfmappings_by_cluster_id($cluster_id, $table);
    };

    return AT::DB::ResultStream->new(next_method => $next_method);
}


sub get_gfmappings_in_region
{
    my ($self, $chr, $start, $end) = @_;
    unless($chr and $start and $end) {
	croak "get_gfmappings_in_region: chr, start or end argument missing";
    }

    my $query = "SELECT m.* FROM gfmap_cluster c, gfmap m
		 WHERE c.chr = ? AND c.end >= ? AND c.start <= ?
		   AND c.cluster_id = m.cluster_id
		   AND m.end >= ? AND m.start <= ?";
   my $sth = $self->dbh->prepare($query);
   $sth->execute($chr, $start, $end, $start, $end);

    my @gfmappings;
    while(my $row = $sth->fetchrow_hashref()) {
        push @gfmappings, $self->_create_gfmapping_from_row($row);
    }
    $sth->finish();

    return \@gfmappings;
}


sub count_gfmappings_in_region
{
    my ($self, $chr, $start, $end, %opt) = @_;

    my $query = 
"SELECT count(gfmap_id)
 FROM gfmap_cluster a, gfmap b
 WHERE a.chr = ? AND a.end >= ? AND a.start <= ?
   AND a.cluster_id = b.cluster_id
   AND b.end >= ? AND b.start <= ?";

    my @bind = ($chr, $start, $end, $start, $end);

    if(my $strand = $opt{strand}) {
	$query .= " AND strand = ?";
	push @bind, $strand;
    }

    if(my $min_score1 = $opt{min_score1}) {
	$query .= " AND score1 >= ?";
	push @bind, $min_score1;
    }

    if($opt{exclude_if_downstream_A_stretch}) {
	$query .= " AND ends_near_A_stretch = 0 AND starts_near_T_stretch = 0";
    }

    my $sth = $self->dbh->prepare_cached($query);
    $sth->execute(@bind);
    my ($count) = $sth->fetchrow_array();
    $sth->finish();
    
    return $count;
}


sub _create_gfmapping_from_row
{
    my ($self, $row) = @_;

    my $gfmapping;

    if($self->{get_primary_mappings}) {
        # Get primary mappings    
        my $mapdb = $self->mapping_db();
        my $query = "SELECT primap_id from gfmap_primap WHERE gfmap_id = ?";
        my $sth = $self->dbh->prepare_cached($query);
        $sth->execute($row->{gfmap_id});
        my @primary_maps;
        while(my ($primap_id) = $sth->fetchrow_array()) {
	    push @primary_maps, $mapdb->get_mapping_by_id($primap_id);
        }
	$sth->finish();
        # Create gfmapping object
        $gfmapping = AT::FT::GFMapping->new(%$row,
					   primary_mappings => \@primary_maps);
    }
    else {
	# Get accs of primary mappings
        my $query = "SELECT acc, type from gfmap_primap a WHERE gfmap_id = ?";
        my $sth = $self->dbh->prepare_cached($query);
        $sth->execute($row->{gfmap_id});
        my (@mRNA_accs, @EST_accs);
        while(my ($acc, $qType) = $sth->fetchrow_array()) {
	    if($qType eq 'mRNA') { push @mRNA_accs, $acc; }
	    elsif($qType eq 'EST') { push @EST_accs, $acc; }
        }
	$sth->finish();
        # Create gfmapping object
        $gfmapping = AT::FT::GFMapping->new(%$row,
					       mRNA_accs => \@mRNA_accs,
					       EST_accs => \@EST_accs);


    }

    # Add hsps and gaps
    my $offset = $gfmapping->{start} - 1;
    my $gap;
    foreach my $hsp_string (split /;/, $row->{hsps}) {
	my ($start, $end, $type, $jnc_str) = split (/,/,$hsp_string);
	my $hsp = AT::FT::GFMappingHSP->new(start => $start+$offset,
					    end => $end+$offset);
	$gfmapping->add_right_HSP($hsp, $gap);
	if($type) {
	    if($type eq 'I') {
	        $gap = AT::FT::GFMappingGap->new(is_broken => 0, is_intron => 1,
					     jnc_str => $jnc_str);
	    }
	    elsif($type eq 'B')  {
		$gap = AT::FT::GFMappingGap->new(is_broken => 1, is_intron => 0,
					     jnc_str => $jnc_str)
	    }
	    else {
	        $gap = AT::FT::GFMappingGap->new(is_broken => 0, is_intron => 0,
					     jnc_str => $jnc_str);
	    }
	}
    }

    return $gfmapping;
}


sub id_string_to_region
{
    my ($self, $id) = @_;

    my ($chr, $start, $end);

    if($id =~ /^(chr.+):([\d,]+)-([\d,]+)/) {
        ($chr, $start, $end) = ($1, $2, $3);
        $start =~ s/,//g;
        $end =~ s/,//g;   
    }
    elsif($id =~ /^gene:(\d+)/) {
        my $gene_id = $1;
        ($chr, $start, $end) = $self->dbh->selectrow_array(
    "select chr, start, end from gene g where gene_id = $gene_id");
    }
    elsif($id =~ /^chain:(\d+)/) {
        my $chain_id = $1;
        ($chr, $start, $end) = $self->dbh->selectrow_array(
    "select g.chr, min(g.start), max(g.end)
     from chain_aspair_bdpair c, gene g
     where chain_id = $chain_id and c.gene_id = g.gene_id group by chain_id");
    }
    elsif($id =~ /^aspair:(\d+)/) {
        my $aspair_id = $1;
        ($chr, $start, $end) = $self->dbh->selectrow_array(
    "select chr, min(start), max(end)
     from aspair p, gene g
     where aspair_id = $aspair_id and (p.gene_id1 = g.gene_id or p.gene_id2 = g.gene_id) group by aspair_id");
    }
    elsif($id =~ /^bdpair:(\d+)/) {
        my $bdpair_id = $1;
        ($chr, $start, $end) = $self->dbh->selectrow_array(
    "select chr, min(start), max(end)
     from bdpair p, gene g
     where bdpair_id = $bdpair_id and (p.gene_id1 = g.gene_id or p.gene_id2 = g.gene_id) group by bdpair_id");
    }
    else {
	my $sth = $self->dbh->prepare("select chr, start, end from gfmap natural join gfmap_primap where acc = ?");
	$sth->execute($id);
	my @locs;
	while(my @loc = $sth->fetchrow_array) {
	    push @locs, \@loc;
	}
	if(@locs == 1) {
	    ($chr, $start, $end) = @{$locs[0]};
	}
	elsif(@locs == 0) {
            $self->{_errstr} = "Could not find $id in database";
	    return;
	}
	else {
	    $self->{_errstr} = "$id has multiple locations: ".join(", ", map { ($_->[0]).':'.($_->[1]).'-'.($_->[2]) } @locs);
	    return;
	}
    }

    return defined($end) ? ($chr, $start, $end) : undef;
}


##################
# DB build methods
##################

sub index_database
{
    my $self = shift;
    
    warn "Indexing aspair table\n";
    $self->_do_or_die("alter table aspair add key (gene_id1)");
    $self->_do_or_die("alter table aspair add key (gene_id2)");
    warn "  done\n";

    warn "Indexing gene tables\n";
    $self->_do_or_die("alter table gene add unique (chr,min_start,max_end,strand)");
    $self->_do_or_die("alter table gene add unique chr_2 (chr,start,end,strand)");
    $self->_do_or_die("alter table gene_gfmap add key(gene_id)");
    $self->_do_or_die("alter table gene_gfmap add key(gfmap_id)");
    warn "  done\n";

    warn "Indexing gfmapping tables\n";
    $self->_do_or_die("alter table gfmap add key chr_start (chr, start)");
    $self->_do_or_die("alter table gfmap add key chr_end (chr, end)");
    $self->_do_or_die("alter table gfmap_primap add key (gfmap_id)");
    $self->_do_or_die("alter table gfmap_primap add unique (primap_id)");
    $self->_do_or_die("alter table gfmap_primap add key (acc)");
    warn "  done\n";
}


sub create_gfmapping_clusters
{
    my ($self) = @_;
    my $table = "gfmap_cluster";
    die "Table $table already exists" if($self->_table_exists($table));

    my $create_query =
"CREATE TABLE $table (
    cluster_id mediumint unsigned not null,
    chr char(5) not null default '',
    start int unsigned not null default 0,
    end int unsigned not null default 0,
    primary key (cluster_id)
)";    
    $self->dbh->do($create_query);

    my $chr_list = $self->dbh->selectcol_arrayref("SELECT DISTINCT chr from gfmap");
    my $select_sth = $self->dbh->prepare("SELECT gfmap_id, start, end from gfmap WHERE chr = ? ORDER BY start");
    my $insert_sth = $self->dbh->prepare("INSERT $table VALUES (?,?,?,?)");
    my $update_sth = $self->dbh->prepare("UPDATE gfmap SET cluster_id = ? WHERE gfmap_id = ?");
    my $cluster_id = 1;
    foreach my $chr (@$chr_list) {
	$select_sth->execute($chr);
	my ($map_id, $cl_start, $cl_end) = $select_sth->fetchrow_array();
	$update_sth->execute($cluster_id, $map_id);
	while(my ($map_id, $map_start, $map_end) = $select_sth->fetchrow_array()) {
	    if($map_start <= $cl_end) {
		$cl_end = $map_end if($cl_end < $map_end);
	    }
	    else {
		$insert_sth->execute($cluster_id, $chr, $cl_start, $cl_end);
		$cluster_id++;
		($cl_start, $cl_end) = ($map_start, $map_end);
	    }
	    $update_sth->execute($cluster_id, $map_id);
	}
	$insert_sth->execute($cluster_id, $chr, $cl_start, $cl_end);
        $cluster_id++;
    }

    $insert_sth->finish();
    $update_sth->finish();

    $self->_do_or_die("ALTER TABLE $table ADD INDEX chr_start (chr, start)");
    $self->_do_or_die("ALTER TABLE $table ADD INDEX chr_end (chr, end)");
    $self->_do_or_die("ALTER TABLE gfmap ADD INDEX (cluster_id)");
}


1;


