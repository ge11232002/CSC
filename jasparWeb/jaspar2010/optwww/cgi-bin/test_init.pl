use lib "/opt/www/jaspar_2010/";

use InitDB;
my ($db_info)=InitDB->init;
  require Data::Dumper;
    print Data::Dumper->Dump([$db_info]);

#print $db_info->collections;

print  $db_info->description("CORE");
