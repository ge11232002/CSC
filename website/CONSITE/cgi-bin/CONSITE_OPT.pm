# server name

use constant WEBSERVER        => 'http://asp.ii.uib.no:8090';

# consite perl library

use constant CONSITE_PERLLIB       => '/opt/www/CONSITE/lib';

# RAVEN perl library

use constant CONSNP_PERLLIB       => '/opt/www/CONSNP_OLD_TEMP/lib';
use constant GL_PERLLIB       => '/opt/www/GENELYNX2/lib';
use constant AT_PERLLIB       => qw'/opt/www/AT/lib /opt/www/AT/Malin';


# Absolute paths

use constant ABS_HTML_DIR     => '/opt/www/CONSITE/html'    ;
use constant ABS_CGI_BIN_DIR  => '/opt/www/CONSITE/cgi-bin' ;
use constant ABS_IMAGES_DIR   => '/opt/www/CONSITE/html/IMG';
use constant ABS_TMP_DIR      => '/opt/www/CONSITE/html/TEMP'    ;
use constant ABS_TEMPLATE_DIR => '/opt/www/CONSITE/TEMPLATES'    ;
use constant LOG_DIR          => '/opt/www/CONSITE/html/TEMP';


# HTML page (RELATIVE) paths

use constant REL_CGI_BIN_DIR  => '/cgi-bin/CONSITE';
use constant REL_IMG_DIR      => '/CONSITE/IMG';
use constant REL_HTML_DIR     => '/CONSITE';
use constant REL_TMP_DIR      => '/CONSITE/TEMP';

# which alignment program to use

use constant ALIGNMENT_PROGRAM => "ORCA";
use constant DEFAULT_ALIGNMENT_METHOD => "orca";

# DPB settings
use constant DPB_CGI => "http://pegasus.cgb.ki.se/cgi-bin/DPB.cgi";
use constant MAX_SEQ_LENGTH => 10000;
use constant DPB_TIMEOUT => 600;

# orca settings

use constant ORCA_BINARY => "/opt/www/orca/orca.pl";
use constant ORCA_LIB    => '/opt/www/orca/Orca';

# Layout

use constant ADMIN_EMAIL   => 'Boris.Lenhard@bccs.uib.no';

# TMP directory cleaning
use constant CLEAN_TEMPFILES_OLDER_THAN => 3; # in days

# development options

use constant DEBUGGING_ON => 1;

# from CONSNP

# human GeneLynx db access

use constant HUMAN_GENELYNX_DB_HOST => "nautilus.cgb.ki.se";
use constant HUMAN_GENELYNX_DB_NAME => "GLHUMAN_1_99";
use constant HUMAN_GENELYNX_DB_USER => "borisl";
use constant HUMAN_GENELYNX_DB_PASS => "chongabrdja";

# mouse GeneLynx db access

use constant MOUSE_GENELYNX_DB_HOST => "nautilus.cgb.ki.se";
use constant MOUSE_GENELYNX_DB_NAME => "GLMOUSE_1_99";
use constant MOUSE_GENELYNX_DB_USER => "borisl";
use constant MOUSE_GENELYNX_DB_PASS => "chongabrdja";

# human assembly db access

use constant HUMAN_ASSEMBLY_DB_HOST => "localhost";
use constant HUMAN_ASSEMBLY_DB_NAME => "HS_MAY04";
use constant HUMAN_ASSEMBLY_DB_USER => "borisl";
use constant HUMAN_ASSEMBLY_DB_PASS => "chongabrdja";

# human mapping db access

use constant HUMAN_MAPPING_DB_HOST => "localhost";
use constant HUMAN_MAPPING_DB_NAME => "AT_HS_MAY04";
use constant HUMAN_MAPPING_DB_USER => "borisl";
use constant HUMAN_MAPPING_DB_PASS => "chongabrdja";


# mouse assembly db access

use constant MOUSE_ASSEMBLY_DB_HOST => "localhost";
use constant MOUSE_ASSEMBLY_DB_NAME => "MM_MAR05";
use constant MOUSE_ASSEMBLY_DB_USER => "borisl";
use constant MOUSE_ASSEMBLY_DB_PASS => "chongabrdja";

# mouse mapping db access

use constant MOUSE_MAPPING_DB_HOST => "localhost";
use constant MOUSE_MAPPING_DB_NAME => "AT_MM_MAR05";
use constant MOUSE_MAPPING_DB_USER => "borisl";
use constant MOUSE_MAPPING_DB_PASS => "chongabrdja";

# JASPAR access

use constant JASPAR_DB_HOST => "localhost";
use constant JASPAR_DB_NAME => "JASPAR2_1";
use constant JASPAR_DB_USER => "borisl";
use constant JASPAR_DB_PASS => "chongabrdja";

# UCSC assemblies
use constant MOUSE_ASSEMBLY => 'mm6';
use constant HUMAN_ASSEMBLY => 'hg17';
use constant RAT_ASSEMBLY   => 'rn3';

# layout parameters

%LAYOUT_PARAMS = (  'COLOR1'            => '#2f2f2f',
                    'COLOR2'            => '#ffffff',
                    'COLOR3'            => '#dfdfdf',
                    'REL_TMP_DIR'       =>  REL_TMP_DIR,
                    'ABS_TMP_DIR'       =>  ABS_TMP_DIR,
                    'REL_CGI_BIN_DIR'   =>  REL_CGI_BIN_DIR,
                    'REL_IMG_DIR'       =>  REL_IMG_DIR,
                    'REL_HTML_DIR'      =>  REL_HTML_DIR,
                    'IMAGE_DIR'         =>  REL_IMG_DIR,
                    'ABS_TEMPLATE_DIR'  =>  ABS_TEMPLATE_DIR,
                    'HUMAN_ASSEMBLY'    =>  HUMAN_ASSEMBLY,
                    'MOUSE_ASSEMBLY'    =>  MOUSE_ASSEMBLY,
                    'RAT_ASSEMBLY'      =>  RAT_ASSEMBLY,
                    

                  );


1;
