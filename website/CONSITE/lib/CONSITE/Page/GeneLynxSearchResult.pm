package CONSITE::Page::GeneLynxSearchResult;

use ConSNPWeb::Page::PageI;
use GeneLynx::Search::MultiSpecies::Quick;
use strict;

use vars qw'@ISA';
@ISA = qw'ConSNPWeb::Page::PageI';


sub output  {
    my ($self) = @_;
    my $q = $self->query;
    my $env = $self->environment;
    my $template = $self->template;
    if (!$self->can("row_color")) {
        $self->{'row_color'} = "#ffffdf";
    }
    my $csearch = GeneLynx::Search::MultiSpecies::Quick->new
                         (-raw_query_string => $q->param("querystring")."",
                          -raw_query_boolean => $q->param("querybool")."",
                          -species_dbs => $env->genelynx_db);

    my @rows;

    my $state_clause = "";
    if ($self->can("state") and $self->state)  {
        $state_clause = "&state_id".$self->state->id;
    }
    foreach my $hit ($csearch->hit_list)  {
        my %mh = $hit->all_matches_hash;
        next unless ($mh{"human"} and $mh{"mouse"});
        my $url = $env->param("REL_CGI_BIN_DIR").
                    "/a?rm=selectrelseq&sp1=human&glid1=".
                    $mh{"human"}->glid.
                    "&sp2=mouse&glid2=".
                    $mh{"mouse"}->glid.
                    "&event=".$state_clause;
        my $img1 = $env->param("REL_IMG_DIR")."/select_pair1.gif";
        my $img2 = $env->param("REL_IMG_DIR")."/select_pair2.gif";
        my $image_name = "im".int(rand(1000));
        push @rows, ($q->td(["Human HS#". $mh{"human"}->glid . $q->br.
                             $env->genelynx_db('human')->get_description($mh{"human"}->glid),
                             $q->a({-style=>"gl",
                                    -href=> $url,
                                    -onMouseOver => "MM_swapImage('$image_name','','$img2',1)",
                                    -onMouseOut => "MM_swapImgRestore()",
                                   }, $q->img({-src=>$img1, -border=>0, -name=>$image_name})),


                             "Mouse MM#" . $mh{"mouse"}->glid. $q->br.
                                 $env->genelynx_db('mouse')->get_description($mh{"mouse"}->glid)."\n"
                            ]));

        push @rows, ($q->td({-colspan=>3}, $q->hr));

    }
    unless (@rows) {
        push @rows, ($q->td({-colspan=>3}, $q->p("The search returned no hits.")));
    }

    $template->param
        ('HIDDEN' => $q->hidden(-name=>"rm",
                                -value=>"t_selectseq",
                                -override=>1)
                    .$q->hidden(-name=>"event",
                                -value=>"gene_pair_selected",
                                -override=>1)
                    ,
          'TITLE' => "RAVEN Search Results",
          'CONTENT' =>  $q->start_form({-action=>"a", -method=>"GET"}).
                        $q->hidden({-name => "rm", -value =>"factorselect"}).
                        $q->table({-bgcolor => $env->param('COLOR2')},
                                    $q->Tr({-bgcolor=>$self->row_color},\@rows)).
                        $q->end_form
          );
    return $template->output();



}

1;
