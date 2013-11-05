package CONSITE::Page::SelectRelativeGenomicSeq;

use strict;
use warnings;

use ConSNPWeb::Page::PageI;
use CONSNP::BioGraphics::TranscriptPanel;

use vars qw'@ISA';
@ISA = qw'ConSNPWeb::Page::PageI';


sub output  {
    my ($self) = @_;
    my $q = $self->query;
    my $env = $self->environment;
    my $template = $self->template;

    my ($SP1, $GLID1, $SP2, $GLID2) = ($q->param("sp1"),
                                       $q->param("glid1"),
                                       $q->param("sp2"),
                                       $q->param("glid2"));

    my %glinfo = ( human => { short => "HS",
                              url  => sub {"http://human.genelynx.org/cgi-bin/record?glid=".$_[0]}
                           },
                   mouse => { short => "MM",
                              url  => sub {"http://mouse.genelynx.org/cgi-bin/record?glid=".$_[0]}
                           },
                   rat => { short => "RN",
                              url  => sub {"http://rat.genelynx.org/cgi-bin/record?glid=".$_[0]}
                           }
                 );




    my @rows;

    my $cdna_res = GeneLynx::Resource->new("cdna_sequences");
    my $ucsc     = GeneLynx::Resource->new("ucsc_golden_path");

    my @sp = ({sp => $SP1, glid => $GLID1, color =>"#34AACD"}, {sp => $SP2, glid => $GLID2, color=>"#409D27"});
    foreach my $i (0..$#sp )  {
        my @accs = GeneLynx::ResourceLinks->new
                    (-glid        =>$sp[$i]->{glid},
                     -resourceobj => $cdna_res,
                     -db          => $env->genelynx_db($sp[$i]->{sp}))->_raw_id_list;

        push @rows, ($q->td({-colspan=>2},uc($sp[$i]->{sp})
                            )

                     ,
                    $q->td({-colspan=>2}, $env->genelynx_db($sp[$i]->{sp})->get_description($sp[$i]->{glid})
                    ." [GeneLynx "
                                              .$glinfo{$sp[$i]->{sp}}->{short}. "#"
                                              .$q->a({ -href => $glinfo{$sp[$i]->{sp}}->{url}->($sp[$i]->{glid})},
                                                    $sp[$i]->{glid}
                                                    )."]\n"));

        my $checked = "CHECKED";


        my @maps;
        my @accrows;
        foreach my $acc (@accs)  {
            if(my $map = ($env->mapping_db($sp[$i]->{sp})->get_mappings_for_acc($acc))[0])
            {
               push @maps, $map;
            }
        }
        my $image_file = "IMG/".$sp[$i]->{sp}."_".$sp[$i]->{glid}.".png";
        # unlesss (#-e ABS_TMP_DIR."/".$image_file)  {
        my $trp = CONSNP::BioGraphics::TranscriptPanel->new_from_mappings
                        ({border_color => $sp[$i]->{color},
                          fill_color   => $sp[$i]->{color}
                         }, @maps);
        if (my $panel = $trp->panel(1))  {
            open F, ">".$env->param("ABS_TMP_DIR")."/".$image_file
                or die "Could not open ".$env->param("ABS_TMP_DIR")."/$image_file";
            print F $trp->_panel1->gd->png;
            close F;
            my $chr;
            foreach my $map (@{$trp->blocks->[0]}) {
                push @accrows, ($q->td({-align=>"left"},["<input type=radio name=acc".($i+1)." value=".$map->qName." $checked>".$map->qName
                                    #.$map->tName.":".$map->tStart."-".$map->tEnd."[".$map->strand."]"
                                    ]));
                $checked= "";
                $chr = $map->tName;

            }
            $chr =~ s/chr//;
            push @rows,
                ($q->td({width=>150},
                        $q->table({-width=>150}, $q->Tr(\@accrows)))
                .$q->td({-width=>450, -align =>"center"},
                        $q->a({ -href => sprintf("http://genome.ucsc.edu".
                                                  "/cgi-bin/hgTracks".
                                                  "?position=chr%s:%d-%d&Submit=Submit&db=%s",
                                                  $chr,
                                                  $trp->panel(1)->start,
                                                  $trp->panel(1)->end,
                                                  $env->assembly($sp[$i]->{sp})
                                                 ),
                               -target => "_blank"
                              },$q->img({-src=>$env->param("REL_TMP_DIR")."/".$image_file,
                                -alt=> "Transcript map for "}))));
            push @rows, ($q->td("").
                         $q->td({-align=>"center"},
                                $q->font({-size=>-2},
                                         "Click on image to view region in UCSC Genome Browser")));
        }
        else {
            push @rows,($q->td({-colspan => 2}, "Could not find cDNA mappings for @maps"));
        }
    }

    push @rows, $q->td({-colspan=>3},"
                       Extract region".
                       $q->textfield(-name => "seqstart",
                                     -size=> 8,
                                     -value => "-3000").
                       "to".
                       $q->textfield(-name => "seqend",
                                     -size => 8,
                                     -value => "500").
                       "relative to".
                       $q->radio_group(-name => 'relative_to',
                                       -values => ['5prime', '3prime'],
                                       -default => '5prime',
                                       -labels => {'5prime' => '5\' end (TSS)',
                                                   '3prime' => '3\' end'
                                                  }
                                       ),
                    );


    $template->param
        ( 'TITLE' => "RAVEN reference sequence selection",
          'CONTENT' =>  $q->start_form({-action=>"a", -method=>"GET"}).
                        $q->hidden(-name => "rm",
                                    -value =>"view",
                                    -override=>1).
                        $q->hidden(-name => "type",
                                    -value =>($self->state ? $self->state->{current_view} : "graphical"),
                                    -override=>1).
                        $q->hidden(-name => "event",
                                    -value =>"genomic_seqs_selected",
                                    -override=>1).
                        $q->hidden( -name=>"sp1",
                                    -value=>$SP1,
                                    -override=>1).
                        $q->hidden(-name=>"sp2",
                                -value=>$SP2,
                                -override=>1).
                        $q->table({-bgcolor => $env->param('COLOR2')},
                                    $q->Tr(\@rows)).
                        $q->submit.
                        $q->end_form
          );
    return $template->output();


}

1;
