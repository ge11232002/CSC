#!/usr/bin/env perl -w

use AT::DB::GenomeMapping;
use AT::DB::GenomeAssembly;
use AT::CMP::CrossMatchFactory;
use DBSNP::MySQLdb;
use strict;
use TFBS::DB::JASPAR2;

#use Test;
#plan(tests => 3);



my $jaspardb = TFBS::DB::JASPAR2->connect("dbi:mysql:JASPAR2:forkhead.cgb.ki.se", "krivan", "wk3003");

my $matrixset = $jaspardb->get_MatrixSet(-sysgroups=>["vertebrate"],
					 -min_ic => 11);

my $mousemapdb = AT::DB::GenomeMapping->connect(-dbhost => "sql.cgb.ki.se",
						-dbname => "MAP_GLMOUSE_0_2",
						-dbuser => "borisl",
						-dbpass => "***");
my $humanmapdb = AT::DB::GenomeMapping->connect(-dbhost => "sql.cgb.ki.se",
						-dbname => "MAP_GLYNX_1_1",
						-dbuser => "borisl",
						-dbpass => "***");
my $mousegenomedb = AT::DB::GenomeAssembly->connect(-dbhost => "sql.cgb.ki.se",
						    -dbname => "MM_FEB02",
						    -dbuser => "borisl",
						    -dbpass => "***");
my $humangenomedb = AT::DB::GenomeAssembly->connect(-dbhost => "sql.cgb.ki.se",
						    -dbname => "HS_APR02",
						    -dbuser => "borisl",
						    -dbpass => "***");
my $dbsnp = DBSNP::MySQLdb->connect("NONmRNA", 
				    "iga.cgb.ki.se", 
				    "borisl", 
				    "***");


while (my $line = <>) {
    print STDERR $line;
    chomp $line;
    my ($hglid, $mglid) = split /\s+/, $line;    

    my $msth = $mousemapdb->dbh->prepare
	(qq! 
	 SELECT qName from MAPPING where qName like "$mglid\\_\\_%" 
	 order by matches desc limit 1
	 !);
    $msth->execute();
    my ($macc) = $msth->fetchrow_array;
    #print "macc: $macc\n";
    next unless $macc;
    
    my $hsth = $humanmapdb->dbh->prepare
	(qq! 
	 SELECT qName from MAPPING where qName like "$hglid\\_\\_%" 
	 order by matches desc limit 1
	 !);
    
    $hsth->execute();
    my ($hacc) = $hsth->fetchrow_array;
    #print "hacc: $hacc\n";
    next unless $hacc;
    
    
    my $cmf =  AT::CMP::CrossMatchFactory->new 
	(assembly1_db => $humangenomedb,
	 mapping1_db  => $humanmapdb,
	 assembly2_db => $mousegenomedb,
	 mapping2_db  => $mousemapdb,
	 min_coverage => 0.9,
	 min_identity => 0.96,
	 rel_start => "-3000",
	 rel_end   => "e2"
	 );
    
    $cmf->run(acc1=>$hacc,
	      acc2=>$macc,
	      extend=>1) or next;
    
    my $cm = $cmf->crossmatch;
    
    
#     my @hexons = $cm->alignment1->query_exon_list;
#     my $i=1;
#     foreach my $e (@hexons) {
# 	printf("  QE %s: %d-%d [%d]\n", $i++, $e->start, $e->end, $e->length);
#     }
    
#     my @htexons = $cm->alignment1->target_exon_list;
#     $i=1;
#     foreach my $e (@htexons) {
# 	printf("  TE %s: %d-%d [%d]\n", $i++, $e->start, $e->end, $e->length);
#     }
    
#     my @hintrons = $cm->alignment1->intron_list;
#     $i=1;
#     foreach my $e (@hintrons) {
# 	printf("  IN %s: %d-%d [%d]\n", $i++, $e->start, $e->end, $e->length);
#     }
    
    # print "-----------------------\n";
   

    my ($promostart, $promoend);
    if ($cm->alignment1->target_exon(1)->strand ==1)  {
	$promostart = 
	    ($cm->alignment1->target_exon(1)->each_tag_value("abs_start"))[0]
	    -2000;
	$promoend = 
	    ($cm->alignment1->target_exon(1)->each_tag_value("abs_start"))[0];
    }
    else  {
	$promostart = 
	    ($cm->alignment1->target_exon(1)->each_tag_value("abs_end"))[0];

	$promoend = 
	    ($cm->alignment1->target_exon(1)->each_tag_value("abs_end"))[0]
	    +2000;
    }
    print STDERR "PROMOTER: $promostart-$promoend "
	.$cm->alignment1->target_exon(1)->strand."\n";

    my %snpset = rs_set_in_region($humangenomedb->dbh,
				   $dbsnp,
				   $cm->mapping1->tName,
				   $promostart, $promoend);
    if ($cm->alignment1->nr_introns)  {
	my ($intron1start, $intron1end) = 
	    (($cm->alignment1->intron(1)->each_tag_value("abs_start"))[0],
	     ($cm->alignment1->intron(1)->each_tag_value("abs_end"))[0]); 
	print STDERR "INTRON: $intron1start-$intron1end\n";
	%snpset = (%snpset, rs_set_in_region($humangenomedb->dbh,
					     $dbsnp,
					     $cm->mapping1->tName,
					     $intron1start, $intron1end));

    }
    
    unless (%snpset) { #print "No SNPs.\n\n"; 
		       next; }
    foreach my $rs (keys %snpset)  {
	my $output;
	my @sequences = snp_to_sequences($cm->genomic1, %{$snpset{$rs}});
	my @sitesets;
	my $it = $matrixset->Iterator;
	while (my $pwm = $it->next) {
	    foreach my $seq (@sequences) {
		push @sitesets,  $pwm->search_seq(-seqstring=>$seq, -threshold=>"50%");
	    }
	    $output = compare_sitesets(@sitesets);
	}
	my $rel_tss;

	if( $cm->alignment1->target_exon(1)->strand == 1) {
	    $rel_tss = $snpset{$rs}->{'start'} - $promoend; 
	}
	else {
	    $rel_tss = $promostart - $snpset{$rs}->{'start'}; 
	    
	}
	$rel_tss++ if $rel_tss>=0;

	print STDERR "REL. TSS $rel_tss\n";
	print STDERR "SNP START ".$snpset{$rs}->{'start'}."\n";
	if ($output) {
	    print ">$hglid rs$rs "
		.$rel_tss
		." (". scalar(@sequences)." variants)\n";
	    
	    print $output;
	    print "//\n";

	}

    }
	
				   

    #print "-----------------------\n";


#     $i=1;
#     my @mexons = $cm->alignment2->query_exon_list;
#     foreach my $e (@mexons) {
# 	printf("    %s: %d-%d [%d]\n", $i++, $e->start, $e->end, $e->length);
#     }
    
#     my @mtexons = $cm->alignment2->target_exon_list;
#     $i=1;
#     foreach my $e (@mtexons) {
# 	printf("  TE %s: %d-%d [%d]\n", $i++, $e->start, $e->end, $e->length);
#     }
    
#     my @mintrons = $cm->alignment2->intron_list;
#     $i=1;
#     foreach my $e (@mintrons) {
# 	printf("  IN %s: %d-%d [%d]\n", $i++, $e->start, $e->end, $e->length);
#     }
    

}



exit(0);



#---------------------------------------------------------

sub rs_set_in_region  {
    my ($dbh, $dbsnp, $chr, $start, $end) = @_;
    my $sth = $dbh->prepare(q!SELECT distinct name, chromStart, chromEnd 
			    from snpNih
			    where chrom= ? 
			    and chromStart >= ?
			    and chromEnd <= ?!);
    $sth->execute($chr, $start, $end);
    my %snp;
    while (my ($rs, $snpstart, $snpend) = $sth->fetchrow_array)  {
	#print "--: $rs $chr$snpstart-$snpend\n";
	my $obj = $dbsnp->get_SNP_by_rs($rs);
	next unless $obj;
	$snp{$rs} = {chr => $chr,
		     start => $snpstart,
		     end   => $snpend,
		     obj => $obj
		     };
    }
    return %snp;
}

  
   
sub snp_to_sequences  {

    my ($seqobj, %snphash) = @_;
    my @sequences;
    my $seqbefore= 
	$seqobj->subseq
	($seqobj->column_from_residue_number($snphash{start}-24),
	 $seqobj->column_from_residue_number($snphash{start}-1));
    my $seqafter = 
	$seqobj->subseq
	($seqobj->column_from_residue_number($snphash{end}+1),
	 $seqobj->column_from_residue_number($snphash{end}+24));
    foreach my $allele ($snphash{obj}->allele_list)  {
	push @sequences, $seqbefore.$allele.$seqafter;
    }
    return @sequences;

}
		   

sub compare_sitesets  {
    my @sitesets = @_;
    my %pos_hash;
    my $output = "";
    foreach my $ss (@sitesets) {
	my $it = $ss->Iterator(-sort_by => "start") ;
	while (my $site = $it->next) {
	    if ($site->start<=25 and $site->end >=25)  {
		push @{$pos_hash{$site->start.".".$site->strand.".".$site->pattern->ID}}, $site;
	    }
	}
	
    }
    foreach my $pos (keys %pos_hash) {
	my @scores = sort {$a<=>$b} (map {$_->score} @{$pos_hash{$pos}});
	my @rel_scores = sort {$a<=>$b} (map {$_->rel_score} @{$pos_hash{$pos}});
	if ($rel_scores[-1] > 0.7 and (@scores < @sitesets or $scores[-1]-$scores[0]>=3))  {
	    my @sites = @{$pos_hash{$pos}};
	    next unless @sites;
	    $output .= join("\t", $pos, $sites[0]->pattern->name,  map {$_->score, $_->seq->seq} @sites)."\n";
	}
    }
    return $output;

}

