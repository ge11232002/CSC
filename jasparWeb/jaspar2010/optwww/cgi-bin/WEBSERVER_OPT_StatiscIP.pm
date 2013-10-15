# settings for WEBINTERFACE using JASPAR DB

#HTML AND CGI PATHS
use constant BASE_DIR => "/opt/www/jaspar_2010/";
use constant BASE_URL => "http://jaspar.genereg.net/";
#use constant BASE_URL => "http://people.binf.ku.dk/albin/jaspar_2010/";

#HTML base 
#use constant HTML_BASE_URL => BASE_URL."html/";
use constant HTML_BASE_URL => BASE_URL;
#use constant HTML_BASE_DIR => BASE_DIR."html/";
use constant HTML_BASE_DIR => BASE_DIR."html/";

# CGI web 
#use constant CGI_BASE_URL => BASE_URL."cgi-bin/";
#use constant CGI_BASE_DIR => BASE_DIR."cgi-bin/";
use constant CGI_BASE_URL => BASE_URL."cgi-bin/jaspar2010/";
use constant CGI_BASE_DIR => BASE_DIR."cgi-bin/";

# rest of things are relative to these - just for simplifiction in the main script
use constant ACTION => CGI_BASE_URL."jaspar_db.pl";
use constant TEMPLATES => HTML_BASE_DIR."TEMPLATES/";
use constant TEMPLATES_URL => HTML_BASE_URL."TEMPLATES/";
use constant TEMP_DIR => HTML_BASE_DIR."/TEMP/"; 
use constant TEMP_URL => HTML_BASE_URL."/TEMP/"; 
use constant ERROR => TEMP_DIR. "error.txt";
use constant ABS_HTML_PATH => BASE_DIR."html";
use constant DOWNLOAD_PATH => HTML_BASE_DIR."DOWNLOAD/";

use constant STAMP_DIR => BASE_DIR."../jaspar/bin/STAMP/";
use constant PWMRANDOM_DIR => BASE_DIR."bin/PWMrandomization/";

use constant SUBMIT_PATH => ACTION."?rm=submission";
use constant SUBMISSION_DIRECTORY => TEMP_DIR."SUBMISSIONS/";


1;
