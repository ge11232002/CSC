#Ancora technical documentation
First version by Pär Engström, April 2008. 
This version by Ge Tan, May 2013. 
This document was produced by the 
[GitHub Flavored Markdown](https://help.github.com/articles/github-flavored-markdown "Markdown").

This document describes the implementation the Ancora web resource. 
For a general description of Ancora and its user interface, 
see the [Ancora publication](http://genomebiology.com/2008/9/2/R34 "Ancora publication").

**Table of Contents**  *generated with [DocToc](http://doctoc.herokuapp.com/)*

- [Ancora technical documentation](#ancora-technical-documentation)
	- [Installation](#installation)
		- [Software dependecies](#software-dependecies)
		- [GBrowse2](#gbrowse2)
		- [GBrowse2 Advanced Installation](#gbrowse2-advanced-installation)
			- [Running GBrowse under FastCGI](#running-gbrowse-under-fastcgi)
			- [User Account Database](#user-account-database)
			- [Displaying Next Generation Sequencing Data](#displaying-next-generation-sequencing-data)
			- [Configuring the Uploaded Track Database](#configuring-the-uploaded-track-database)
		- [UCSC Genome Browser source and utilities](#ucsc-genome-browser-source-and-utilities)
		- [ProServer DAS server](#proserver-das-server)
		- [Software components](#software-components)
	- [Configuration](#configuration)
		- [Ancora Web Resource](#ancora-web-resource)
		- [Apache Configuration](#apache-configuration)
		- [Data Files](#data-files)
		- [The CNE database](#the-cne-database)
		- [MySQL account for Perl scripts](#mysql-account-for-perl-scripts)
	- [Generating CNEs](#generating-cnes)
		- [Obtaining alignment files](#obtaining-alignment-files)
		- [Creating filter files](#creating-filter-files)
		- [Scanning for CNEs](#scanning-for-cnes)
		- [Removing unannotated repeats with BLAT](#removing-unannotated-repeats-with-blat)
		- [Loading CNEs into the cne database](#loading-cnes-into-the-cne-database)
	- [Setting up a genome browser](#setting-up-a-genome-browser)
		- [Creating an annotation database for GBrowse](#creating-an-annotation-database-for-gbrowse)
		- [Obtaining genome annotations](#obtaining-genome-annotations)
		- [Creating and loading GFF files](#creating-and-loading-gff-files)
		- [Creating a configuration file for GBrowse](#creating-a-configuration-file-for-gbrowse)

## Installation
 This documentation focuses on the Ancora installation on Olifant at MRC CSC. 
 Olifant is running the CentOS release 5.8 (Final). 
 However, ancora is supposed to run on other Linux/Unix distribution without too much difficulty.

### Software dependecies
 The following softwares are essential for the whole implementation.
 
 *  A httpd server (Currently [Apache](http://httpd.apache.org/docs/2.2/ "Apache") 2.2.3 is used on olifant)
 *  [MySQL](http://dev.mysql.com/downloads/mysql/5.0.html "MySQL") client if you have a dedicated MySQL server. if not, MySQL server is alse required. (MySQL 5.0.95 is used on olifant)
 *  [BioPerl](http://www.bioperl.org/wiki/Installing_BioPerl_on_Unix "BioPerl"). Installation is quite straightforward following the instructions. (BioPerl 1.6.1 is used on olifant.)

### GBrowse2
This implementation of Ancora utilizes the latest version of GBrowse 2.54. 
GBrowse 2.X is a complete rewrite of GBrowse 1.X version. 
There are several advantages compared to the GBrowse 1.X. 
See the [details](http://gmod.org/wiki/GBrowse_2.0_HOWTO#Installation_from_Source_Code "GBrowse").

Since the GBrowse 2.X is perl-based and 
all the modules are hosted on [CPAN](http://search.cpan.org/~lds/GBrowse-2.54/), 
the easiest way to install GBrowse 2 is using the standard Perl module 
building procedure.

For a smooth installation, 
please install some prerequisites before launching the GBrowse2 installation. 
Details are mentioned on the [page](http://gmod.org/wiki/GBrowse_2.0_Prerequisites). 
It is highly recommened to use the common package managers, Debian(DEB) 
or RedHat Package Manager(RPM).

After installing the prerequisites, run this command in terminal

```sh
sudo perl -MCPAN -e 'install Bio::Graphics::Browser2'
```

If not all the necessary prerequisites are installed, it will be notified during the test step. 
During the installation, there are some prompted configuration questions about the location of files.
In principle, the default paths are fine and you can customize it to fit your needs.
On olifant, we configrued as follows:

* Directory for GBrowse's config and support files? [/etc/gbrowse2] /opt/www/gbrowse2/conf
* Directory for GBrowse's static images & HTML files? [/var/www/html/gbrowse2] /opt/www/gbrowse2/html
* Directory for GBrowse's temporary data [/var/tmp/gbrowse2] /opt/www/gbrowse2/tmp
* Directory for GBrowse's sessions, uploaded tracks and other persistent data [/var/lib/gbrowse2] /opt/www/gbrowse2/lib
* Directory for GBrowse's example databases [/var/lib/gbrowse2/databases] /opt/www/gbrowse2/databases
* Directory for GBrowse's CGI script executables? [/var/www/cgi-bin/gb2] /opt/www/gbrowse2/cgi-bin
* Internet port to run demo web site on (for demo)? [8000]
* Apache loadable module directory (for demo)? [/etc/httpd/modules]
* User account under which Apache daemon runs? [apache]
* Automatically update Apache config files to run GBrowse? [y]
* Automatically update system config files to run gbrowse-slave? [y]

After the installation, the GBrowse2 demos should be accessable now. 
Type ```YourServersIP/gbrowse2``` in the internet browser and 
try the example database *yeast*.

However, because of some bugs in ```Bio::Graphics 2.33```, 
the tracks generated by plugins might complain some errors like 

```sh
Can't locate object method "height" via package "gdTinyFont"
(perhaps you forgot to load "gdTinyFont"?) at /usr/local/share/perl/5.8.8/Bio/Graphics/Panel.pm line 641. at
/usr/local/lib/perl/5.8.8/Bio/Graphics/Browser2/Render.pm line 3678.
```

As long as Lincoln D. Stein does not fix the bug, just downgrade to the ```Bio::Graphics 2.32```. Here is the trick how to do it.

```sh
wget http://backpan.perl.org/authors/id/L/LD/LDS/Bio-Graphics-2.32.tar.gz
tar xzf Bio-Graphics-2.32.tar.gz
perl Build.PL
./Build
./Build test
sudo ./Build install --uninst 1
```

The "--uninst 1" will make sure that the files for the
newer (2.33) version are removed in case they are installed 
in another directory.
**Note**: For some unknown reasons, the tracks are not displayed on Mac OS 10.8.3 with Firefox 20.0, Safari 6.0.4 and Opera 12.15. However, Chrome 26.0.1410.65 works fine. So far, it works on Windows 7 with all major browsers.

### GBrowse2 Advanced Installation
This section is optional and incomplete. The advance installation should cover the topics such as *FastCGI*, *User Account Database*, *Displaying Next Generation Sequencing Data*, *Configuring the Uploaded Track Database*.

#### Running GBrowse under FastCGI
The idea of FastCGI is to make the script as a long-running process at the first time of the script is requested. 
By eliminating the startup time for later use of script, 
the responsiveness of a FastCGI is significantly improved.

To use this facility, you will need a version of Apache equipped with FastCGI support, 
the module [mod_fcgid](http://httpd.apache.org/mod_fcgid/). 
For a Debian (DEB) based system,

```sh
sudo apt-get install libapache2-mod-fastcgi libfcgi-perl
```

or for Cent OS,

```sh
cd /etc/yum.repos.d/
wget http://centos.karan.org/kbsingh-CentOS-Extras.repo
sudo vim /etc/yum.repos.d/kbsingh-CentOS-Extras.repo 
# set gpgcheck to 0 and enabled to 1 in the [kbs-CentOS-Testing] section.
yum install mod_fcgid
# Down the latest atomic-release rpm from http://www6.atomicorp.com/channels/atomic/centos/5/x86_64/RPMS/
# Install atomic-release rpm:
rpm -Uvh atomic-release*rpm
yum install fcgi-perl
```

**Note**: Because the olifant has a old version of *Apache* and *mod_fcgid*. 
You may need to adapt the GBrowse2 conf file for Apache 
(/etc/httpd/conf.d/gbrowse2.conf on olifant) 
to the older version of *mod_fcgid* with the names mapping 
described on the [page](http://httpd.apache.org/mod_fcgid/mod/mod_fcgid.html).

#### User Account Database
To be implemented

#### Displaying Next Generation Sequencing Data
To be implemented

#### Configuring the Uploaded Track Database
To be implemented

### UCSC Genome Browser source and utilities
The UCSC Genome Browser source code is needed to compile the C program 
for detecting CNEs, which makes use of several library functions 
from the UCSC source. 
The UCSC Genome Browser source also contains a collection of useful utilities, 
some of which are needed at steps below. 
It is recommended to compile all the UCSC utilities by following the instructions in the README file in the UCSC source package. 
Short story on olifant,

```sh
export MACHTYPE=x86_64
export MYSQLINC=/usr/include/mysql
export MYSQLLIBS="/usr/lib64/mysql/libmysqlclient.a -lz"
export PATH=$HOME/bin/$MACHTYPE:$PATH
mkdir -p $HOME/bin/$MACHTYPE.
wegt http://hgdownload.cse.ucsc.edu/admin/jksrc.zip
unzip jksrc.zip
cd kent/src
make libs
cd utils/
make
```

### ProServer DAS server
To be implemented

### Software components
Ancora requires two software packages developed in Boris Lenhard’s group: 
cne and AT. The cne package contains Perl libraries, 
C source code and Perl and bash scripts for detecting and 
working with CNEs, as well as Perl libraries and CGI scripts 
directly underlying the Ancora web resource. 
The AT package is a collection of Perl modules related to genome 
and transcriptome sequence analysis. 
They were largely implemented by Pär Engström during his PhD studies 
but also include contributions from several other past 
and present members of the Lenhard group. 
Ancora only uses a few well-tested modules from the AT package, 
which also includes many modules in early stages of development.

The cne and AT packages are maintained in the CSC repository. 
To retrieve the most recent versions, do:

```sh
git clone git@github.com:ge11232002/CSC.git
```

One copy of the packages should be placed under ```/opt/www/cne``` and ```/opt/www/AT```, respectively, 
where the Ancora genome browser will access them from.
This could be done by creating soft links for the ```cne``` and ```AT``` under ```/opt/www/```.
To run scripts that use the modules, 
you will have to add the installation paths to your PERL5LIB environment variable. 
If you use bash as your shell, 
you can do this by adding the following line to the file ```~/.bashrc```:

```sh
export PERL5LIB=/opt/www/AT/lib:/opt/www/cne/perl_lib:$PERL5LIB
```

(If you want scripts to use a working copy of the modules in a different location, modify the above line accordingly.)
The cne/tools directory in the cne package contains source code for several C programs. 
Ancora currently only needs one of these programs: ceScan, 
the program we use to detect CNE. 
Instructions for how to compile the C programs are in cne/tools/README.txt.

## Configuration
This section will describe all the necessary configurations step by step
for setting up the Ancora web resource.

### Ancora Web Resource
Under the repository, there are two other folders *ancora* and *gbrowse2*,
besides the *cne* and *AT* packages.
The *ancora* directory contains html and cgi files 
related to the Ancora web resource,
while the *gbrowse2* directory contains the *conf* files for GBrowse2.
If you are setting up the Ancora web resource, 
you will need to copy some files to the ancora and gbrowse2 directory on your server. 
On olifant, they are ```/opt/www/ancora``` and ```/opt/www/gbrowse2``` 
(As we installed GBrowse2 to this directory).
The following commands create the links there:

```sh
cd ancora
## the html files
ln -s html /opt/www/ancora/
## the cgi files
ln -s cgi-bin /opt/www/ancora/
## the das server settings
ln -s Bio-Das-ProServer /opt/www/ancora/
## put the cgi files to GBrowse cgi-bin folder
cd gbrowse2
ln -s cgi-bin/* /opt/www/gbrowse2/cgi-bin/
## put the gbrowse conf files to the GBrowse's conf folder 
ln -s conf/*.conf /opt/www/gbrowse2/conf/
## put the gbrowse plugins to the GBrowse's plugins folder
ln -s conf/plugins/* /opt/www/gbrowse2/conf/plugins/
```

### Apache Configuration
The Apache httpd server is configured to look for html files 
under ```/var/www/html``` and CGI scripts in ```/var/www/cgi-bin```.
However, if we use the technique called "VirtualHost",
the html files can be left under ```/opt/www/ancora```.
The mapping of the domain http://ancora.olifant.cscdom.csc.mrc.ac.uk/ 
to this directory is configured in Apache configuration file ```/etc/httpd/conf```.

```
<VirtualHost *:80>
	ServerName ancora.olifant.cscdom.csc.mrc.ac.uk
	ScriptAlias /cgi-bin /opt/www/gbrowse2/cgi-bin
	DocumentRoot /opt/www/ancora/html
    <Directory /opt/www/gbrowse2/cgi-bin>
    	Options FollowSymLinks Indexes Includes ExecCGI
    	AllowOverride Indexes
    </Directory>
</VirtualHost>
```

The last three lines is for the usage of symlinks of CGIs. 
By default, apache does not allow running the symlinks in the cgi-bin folder.

**Note**: When the CGI script is called by ```/cgi-bin/gbrowse```,
the CGI will run in the traditional mode.
If you want to take the advantage of FastCGI, 
the CGI will be executed in ```/fgb2/gbrowse```.
The settings of FastCGI for GBrowse2 is in ```/etc/httpd/conf.d/gbrowse2.conf```.
By default, Ancora will call the CGI in FastCGI mode.

### Data Files
Ancora expects genome assembly sequences in binary 2bit format to 
be present under ```/export/data/goldenpath```. 
There should be one subdirectory for each assembly, 
named using a UCSC assembly identifier (e.g. hg18, mm9) 
or some other suitable identifier for assemblies not served by UCSC. 
Each such subdirectory should contain a file named assembly.2bit 
containing the soft-masked assembly sequence. 
Two-bit files can be made by downloading the assembly 
in soft-masked fasta format from http://genome.ucsc.edu/ 
and converting it to 2bit format using the program faToTwoBit, 
which is part of the UCSC utilities. E.g. 
for the latest human assembly:

```sh
cd /export/data/goldenpath
mkdir hg18
cd hg18
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/bigZips/chromFa.zip
unip chromFa.zip
faToTwoBit *.fa assembly.2bit
rm *.fa chromFa.zip
```
We use additional files (annotation and aligment files) from UCSC 
in the CNE detection pipeline and 
in constructing the other annotation tracks shown in the genome browser. 
How to retrieve these files is described in the corresponding sections below.

### The CNE database
Ancora reads information about genome assemblies and CNEs 
from a database called cne served by the MySQL server on the same host. 
Programmatic access to this database is provided 
by the Perl module CNE::DB in the cne package.

If this database does not exist, create it:

```sh
mysql -u username -p -e 'create database cne'
mysql -u username -p -e 'grant select on cne.* to nobody@localhost'
```
(Replacing username with your MySQL username, or root if necessary.) 
The second command gives user nobody read access to the database. 
This is necessary because Ancora accesses the database as user nobody.

Presently there is only one table that must exist in the database. 
This table is called assembly and holds information 
about the genome assemblies for which comparisons are available in Ancora. 
The cne package contains a file with SQL commands 
that create the table and fills it with sample data:

```sh
mysql -u username -p cne < cne/scripts/cne_pipeline/create_assembly_table.sql
```

If the assembly you are setting up a browser for is not already present in this table, 
you need to execute an SQL INSERT statement
to add a row describing the assembly. 
Most of the fields in the table are self-explanatory. 
The following may not be:
*	ensembl_ver – Currently only used by the DAS service. Specifies which Ensembl version that DAS tracks should be added to. Leave blank for the most recent Ensembl release, or add a version in the form of an archive name (e.g. “apr2007”) for an older release.
*	default_ensembl_loc – Currently only used by the DAS service. Specifies which location the user should be taken to in Ensembl when tracks are added. Specify as an Ensembl location string (e.g. “7:8541098-8656549”).
*	ucsc_db – Deprecated, so can be left blank. Earlier releases of Ancora required an UCSC annotation database to be present for some assemblies.

To add one new row into the table or update one row,

```sql
INSERT INTO assembly 
(assembly_id, assembly_name, organism_common, organism_latin, ensembl_ver, default_ensembl_loc) 
VALUES 
("hg19", "NCBI Build 37", "human", "Homo sapiens", "feb2009", "11:31766034-31797085");
UPDATE assembly
SET ensembl_ver="nov2009", default_ensembl_loc="GL172646.1:2369612-2552500"
WHERE assembly_id="xenTro3";
/* Finally output the table for backup*/
mysqldump -u username -p cne assembly > create_assembly_table.sql 
```

### MySQL account for Perl scripts
Generating CNEs and setting up the annotation database for Ancora 
involves running several Perl scripts that need to access local MySQL databases. 
Some of these scripts require that a MySQL username and password 
are specified in a file called MyPerlVars.pm 
located in a directory mentioned in your PERL5LIB environment variable. 
I have such a file in a directory called .perl under my home directory, 
and have added the following line to my .bashrc to 
make Perl look for libraries in this directory:

```sh
export PERL5LIB=~/.perl:$PERL5LIB
```

The file MyPerlVars.pm should contain the following lines:

```perl
package MyPerlVars;
use warnings;
use strict;

use Exporter;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw($sqlUser $sqlPass);

our $sqlUser = "username"; # replace username with your username
our $sqlPass = "password"; # replace password with your password

1;
```

For security reasons, 
it is a good idea to change the permissions of MyPerlVars.pm file 
so that nobody but you can read it.

## Generating CNEs
To scan for CNEs, 
you need to have the cne package, the AT package and 
genome assembly sequences installed as outlined above. 
In addition, you need to obtain alignment files for the genomes 
that are to be compared, and create filter files that specify 
which regions (typically exons and repeats) 
that are to be ignored in the scanning. 
The alignment and filter files are not required to run the web resource, 
so they can be removed after CNE generation.

### Obtaining alignment files
The program that scans for CNEs (ceScan) takes alignments in axt format as input. 
We typically obtain alignment files in axt format from UCSC and 
keep them under ```/export/downloads/ucsc/axtNet```.
This directory contains one subdirectory for each genome assembly. 
The directories are named using UCSC assembly identifiers, 
as for the genome assembly directories (see above). 
Each directory contains axt files for pairwise net alignments 
with the corresponding assembly as the reference; 
e.g. directory ```/export/downloads/ucsc/axtNet/hg18``` 
only contains net alignments where hg18 is reference. 
To detect CNEs between two genomes, 
you need two sets of net alignments - 
one where each genome is the reference 
(see the [Ancora](http://genomebiology.com/2008/9/2/R34) 
paper for an explanation of the rationale behind this).
The files we download from UCSC are typically compressed with gzip. 
There is no need to decompress them, 
because ceScan can read gzip-compressed files.

As an example, to detect CNEs between the current human and mouse genome assemblies, 
the required alignment files can be obtained as follows:
```sh
rsync -avzP \
  rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/vsMm10/axtNet/* \
  /export/downloads/ucsc/axtNet/hg19/
rsync -avzP \
  rsync://hgdownload.cse.ucsc.edu/goldenPath/mm10/vsHg19/axtNet/* \
  /export/downloads/ucsc/axtNet/mm10/
```
I prepared a script ```cne/scripts/cne_pipeline/downloadAlignments.r``` 
to download the pairwise alignments.

### Creating filter files
Filer files list regions to be ignored in the scan. 
One filter file should be created for each assembly.
The current file for human assembly hg19 is ```/export/data/CNEs/hg19/filters/filter_regions.hg19.bed```. 
and files for other assemblies are in corresponding locations. 
Filter files should adhere to the three-column bed format. 
The filter files do not have to be sorted or non-redundant 
(i.e. they may contain overlapping features).

We produce filter files from exon and repeat annotation.
Please refer to the section “Obtaining genome annotations” 
below for information about how to download exon and repeat annotation 
from UCSC and load it into a local MySQL database or from other sources.
Here's genes and repeats information we used.

| Assembly   |   name    |        Exon  |        Repeat  |
| -----------|-----------|------------- | -------------- |
| hg19       | Human     |RefSeq Genes, Ensembl Genes, UCSC Known Genes | RepeatMasker |
| mm10       | Mouse     |RefSeq Genes, Ensembl Genes, UCSC Known Genes | RepeatMasker |
| xenTro3    | Frog      |RefSeq Genes, Ensembl Genes                   | RepeatMasker |
| tetNig2    | Tetraodon |Ensembl Genes                                 |              |
| canFam3    | Dog       |RefSeq Genes, Ensembl Genes                   | RepeatMasker |
| galGal4    | Chicken   |RefSeq Genes, Ensembl Genes                   | RepeatMasker |
| danRer7    | Zebrafish |RefSeq Genes, Ensembl Genes                   | RepeatMasker |
| fr3        | Fugu      |RefSeq Genes                                  | RepeatMasker |
| anoCar2    | Lizard    |Ensembl Genes                                 | RepeatMasker |
| equCab2    | Horse     |RefSeq Genes, Ensembl Genes                   | RepeatMasker |
| oryLat2    | Medaka    |RefSeq Genes, Ensembl Genes                   | RepeatMasker |
| monDom5    | Opossum   |RefSeq Genes, Ensembl Genes                   | RepeatMasker |
| gasAcu1    |Stickleback|RefSeq Genes, Ensembl Genes                   | RepeatMasker |
| rn5        | Rat       |RefSeq Genes, Ensembl Genes                   | RepeatMasker |
| dm3        | D. melanogaster| RefSeq Genes, Ensembl Genes             | RepeatMasker |
| droAna2    | D. ananassae   |                                         | RepeatMasker |
| dp3        |D. pseudoobscura|                                         | RepeatMasker |
| ce4        | C. elegans     | RefSeq Genes                            | RepeatMasker |
| cb3        | C. briggsae    |                                         | RepeatMasker |
| caeRem2    | C. remanei     |                                         | RepeatMasker |
| caePb1     | C. brenneri    |                                         | RepeatMasker |
| C_albicans_SC5314_A21 |     | CDG exon                               | CGD repeat_region, long_terminal_repeat  |
| C_dubliniensis_CD36   |     | CDG exon                               | CGD repeat_region, long_terminal_repeat  |
| C_glabrata_CBS138     |     | CDG exon                               | CGD repeat_region, long_terminal_repeat  |
| C_parapsilosis_CDC317 |     | CDG exon                               | CGD repeat_region, long_terminal_repeat  |

To prepare the exon file from ensembl gff, use the script ```cne/scripts/cne_pipeline/ens4filters.r```.

Once the data is in a local database and other exons files are ready,
the script ```cne/scripts/cne_pipeline/create_filter.pl``` in the cne package 
can be used to create a filter file from it. 
Run the script without arguments for usage information. 
The username and password that the script should use to connect to the database must be declared in MyPerlVars.pm (see above).

```sh
perl /opt/www/cne/scripts/cne_pipeline/create_filter.pl \
-g refGene,knownGene -b hg19 -r -a /export/data/goldenpath/hg19/assembly.2bit \
-f ensembleExons.txt \
> /export/data/CNEs/hg19/filters/filter_regions.hg19.bed &
```

### Scanning for CNEs
The cne package contains a Perl script cne/scripts/cne_pipeline/detect_cnes.pl 
that we run to detect CNEs. 
The simplest way to run the script is as follows:

```sh
perl /opt/www/cne/scripts/cne_pipeline/detect_cnes.pl all
```

If executed in this way, 
the script will scan for CNEs for all pairwise comparisons and thresholds 
listed in the configuration file ```comparisons.txt``` 
located in the same directory as the script itself. 
Running the script with a pair of assembly ids as arguments 
causes it to detect CNEs only for that particular pairwise comparison 
(which must be listed in the configuration file):

```sh
perl /opt/www/cne/scripts/cne_pipeline/detect_cnes.pl hg18 mm9
```

If thresholds are also specified on the command line, 
the script will not attempt to read a configuration file, 
but simply generate CNEs for the indicated comparison and thresholds:

```sh
perl /opt/www/cne/scripts/cne_pipeline/detect_cnes.pl hg18 mm9 45,50 48,50 49,50
```

Each threshold is specified as two values separated by a comma: 
identities,columns - e.g. 45,50 for 45 (95%) 
identities over 50 alignment columns. 
Any number of thresholds can be given. 
The configuration file is formatted in the same manner; 
see the default file comparisons.txt for an example. 
Some additional usage information can be obtained 
by running the script without any arguments.

Briefly, the script does the following:
For each pair of genome assemblies a and b that are to be compared,
* Run ceScan to detect CNEs based on net alignments with a as reference.
* Run ceScan to detect CNEs based on net alignments with b as reference.
* Run the script merge_twoway_cne_sets.pl to combine the output 
	from the two ceScan runs. 
    In this process, any CNEs that overlap on both genomes are merged.

Output files are placed in the current working directory. 
The files are named according the parameters of the comparison. 
For example, a search for CNEs between assemblies hg19 and mm10 
at a threshold of 49 identities over 50 columns generates a file 
named cne2w_hg19_mm10_49_50. 
Assembly ids are ordered alphanumerically in the filenames. 
We currently keep these files under /export/data/CNEs/pre-blatFilter/.

The output files contain one line for each CNE. 
There are nine columns:
*	Sequence name in first assembly
*	Start coordinate in first assembly
*	End coordinate in first assembly
*	Sequence name in second assembly
*	Start coordinate in second assembly
*	End coordinate in second assembly
*	Alignment orientation (+ or -)
*	Percent identity (percentage of alignment columns with identities)
*	CIGAR string, according to the DAS specification for alignments (see http://www.dasregistry.org/extension_alignment.jsp).

The coordinate ranges are half-open, zero-based as in the bed format, 
so that three-column bed files can be generated by cutting out columns 1-3 or 4-6.

### Removing unannotated repeats with BLAT
This step has been implemented as a separate script because it takes long to run. 
The goal of this step is to remove CNEs that correspond to repeated sequence 
not annotated as such in the RepeatMasker track from UCSC. 
The strategy is simple: use blat to align all CNEs back to 
the genome and remove those having more than N hits 
with percent identity of at least *I*. 
Based on some testing by David Fredman we use I = 90% for all genomes. 
We use N = 4 for all genomes except the fish genomes, 
for which we use N = 8. 
The largest family of duplicated CNEs we have found is 
linked to the IRX genes and comprises 4 homologous CNEs 
in the human genome, 
motivating the chosen values of N for land vertebrates. 
The two-fold higher value for fish is motivated 
by the extra gene duplication that occurred in the teleost lineage. 
We are currently using N = 4 for nematodes and flies as well, 
although we have not explored CNE duplication in those genomes. 
Note that N should not be set to anything lower than 2 
because many assemblies contain some artificially duplicated sequence 
(this means that N=4 may in fact be too strict for assemblies 
that may contain artificially duplicated sequence from the IRX clusters, 
e.g. an IRX cluster contig not assigned to a chromosome; 
we have not explored this further).

The blat filter is implemented in the script ```cne/scripts/cne_pipeline/blat_filter.pl```.
It can simply be run with the cne files as arguments, e.g.:

```sh
perl /opt/www/cne/scripts/cne_pipeline/blat_filter.pl cne2w_hg18_mm9_45_50 cne2w_hg18_mm9_48_50 cne2w_hg18_mm9_49_50
```

Any number of CNE files can be given as arguments. 
For each input file, a corresponding file with filtered CNEs will be created in the working directory. 
This file will have the same name, 
but with the prefix cne2wBf instead of cne2w. 
The format is also the same, 
except that four columns with diagnostic information are added at the end (see the script for details). 
Future versions of the script might not output these extra columns. 
We currently keep the output files under ```/export/data/CNEs/blatFiltered```.

The script reads the N-values for different assemblies 
from a configuration file. 
By default, the script looks for a file named ```blat_hit_cutoffs.txt``` 
located in the same directory as the script itself. 
Run the script without any arguments for more information 
about how it can be configured.

### Loading CNEs into the cne database
The script ```cne/scripts/cne_pipeline/load_cnes.pl``` loads CNEs into a MySQL database. 
For example, the following command loads two files 
with hg19-mm10 CNEs into the database cne:

```sh
perl /opt/www/cne/scripts/cne_pipeline/load_cnes.pl \
	-d cne -r cne2wBf_hg19_mm10_29_30 cne2wBf_hg19_mm10_30_30
```

## Setting up a genome browser
Setting up a genome browser for a new assembly involves 
creating and loading CNEs for the assembly, 
as described in the previous section, 
as well as a number of additional steps detailed in this section.

### Creating an annotation database for GBrowse
Create a MySQL database to store the annotation (genes etc) 
that is to be shown alongside the CNEs. 
For the main Ancora installation on olifant, 
we name these databases gbrowse_gff_assembly_id, 
e.g. gbrowse_gff_hg19 for the latest human assembly. 
Make sure the user nobody@localhost has SELECT permission on the database
as described above. 

```sh
mysql -u root -p -e 'create database gbrowse_gff_hg19'
mysql -u root -p -e 'grant select on gbrowse_gff_hg19.* to nobody@localhost'
```

Note: please make sure the work in database name can not start with capital, i.e., gbrowse_gff_Hg19 will not work.

### Obtaining genome annotations
This step usually involves more manual work than the other steps, 
because annotations are available in diverse databases and 
in different and changing formats. 
To keep things simple, 
we have tried to retrieve as much annotation as possible 
from the UCSC Genome Browser database.
Begin by creating a local database to hold UCSC annotations for the assembly, 
if one does not already exist. 
We have named these databases UCSC_id, 
replacing id with the UCSC assembly id, e.g. UCSC_hg18, UCSC_mm9.

The annotations available for download from UCSC can be browsed 
at http://hgdownload.cse.ucsc.edu/downloads.html. 
Annotations files can be downloaded in batch with rsync.

```sh
rsync -avzP \
  rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/*_gap.* \
  rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/*_rmsk.* \
  /export/data/CNEs/hg19/annotation/
```

The UCSC annotations come as MySQL table dumps, 
with one table definition file (ending in .sql) 
and one data file (ending in .txt.gz) for each table. 
The following series of commands can be used to 
load a set of downloaded files into a local MySQL database 
(named UCSC_hg19 in this example).

```sh
## first create the database and grant the select to nobody
mysql -u root -p -e 'create database UCSC_hg19' 
mysql -u root -p -e 'grant select on UCSC_hg19.* to nobody@localhost' 
## unzip and load the dump
gunzip *.txt.gz
cat *.sql | mysql -u root -p UCSC_hg19
mysqlimport -u root -Lp UCSC_hg19 *.txt
```

**Annotation from UCSC used in Ancora browsers**

| Annotation type |old UCSC tables | new UCSC table |
| --------------- |--------------- | -------------- |
| Gaps            | \*_gap.\*      | gap.\*         |
| RepeatMasker    | \*_rmsk.\*     | rmsk.\*        |
| RefSeq genes    | refGene.\*, refLink.\* | refGene.\*, refLink.\*|
| UCSC Known Genes| knownGene.\*, knownIsoforms.\*, kgXref.\* | knownGene.\*, knownIsoforms.\*, kgXref.\* |
| Aligned UCSC Known Gene peptides from other species |  blastKG\*, kgXref.\* | kgXref.\*|
| ZFIN Genes      | vegaGene.\* | vegaGene.\* | 
| CpG islands     | cpgIslandExt.\*| cpgIslandExt.\*|
| ORegAnno Regulatory Elements | oreganno.\*, oregannoAttr.\*|  oreganno.\*, oregannoAttr.\*  |
| Synteny blocks  | netDanRer7.\*  | netDanRer7.\* |

The Ancora browsers also use annotation from sources other than UCSC. 
For each of these sources, 
we download flat files as detailed below.

NOTE: Synteny blocks are not included currently. Not sure the settings of parameters.

**Synteny blocks**

Download the net files (For example, netDanRer7.sql and netDanRer7.txt.gz)
from UCSC annotation database and load them into mySQL UCSC database (UCSC_hg19).

```sh
## for hg19 with danRer7
cd /export/data/CNEs/synteny
perl /opt/www/cne/scripts/synteny/join_nets.pl danRer7 hg19 UCSC_danRer7 100000 300000
```

Please adapt the database username , password, host and port in ```join_nets.pl``` for your environment.

**Ensembl genes**

Tab-separated text file obtained using BioMart. 
Choose the most recent Ensembl version for your assembly of interest.
Under Dataset, choose genes for the organism of interest (e.g. "Homo sapiens genes").
Do not apply any filters. 
Choose to get the following attributes (categorized under "Structures"):
* Ensembl Gene ID
* Chromosome Name
* Ensembl Transcript ID
* Gene Start (bp)
* Gene End (bp)
* Transcript Start (bp)
* Transcript End (bp)
* Strand
* Associated Gene Name
* Exon Chr Start (bp)
* Exon Chr End (bp)
* Genomic coding start
* Genomic coding end

To convert the downloaded file into gff, 

```sh
perl /opt/www/cne/scripts/gbrowse_db/ens2gff.pl \
    /export/data/goldenpath/hg19/assembly.2bit \
    /export/data/CNEs/hg19/annotation/ensembl_genes.txt \
    >>/export/data/CNEs/hg19/gff/hg19.gff
```

**miRBase microRNAs**

GFF3 file obtained from http://microrna.sanger.ac.uk/sequences/ftp.shtml.
To convert the downloaded gff3 to gff,

```sh
  perl /opt/www/cne/scripts/gbrowse_db/mirbase2gff.pl \
    /export/data/goldenpath/hg19/assembly.2bit \
    /export/data/CNEs/hg19/annotation/hsa.gff3 \
    >>/export/data/CNEs/hg19/gff/hg19.gff
```
Note: Because of the old gff version of c. elegans 4 (WS170),
the line 61 should be uncommented and line 63 should be commented 
in ```mirbase2gff.pl```.
This is not good. Should be fixed later.

**MGI genes**

This is just for Mouse.
Tab-delimited MGI coordinate file 
(currently named MGI_Gene_Model_Coord.rpt) 
obtained from ftp://ftp.informatics.jax.org/pub/reports/index.html 
```sh
  perl /opt/www/cne/scripts/gbrowse_db/mgi2gff.pl \
    /export/data/goldenpath/mm10/assembly.2bit \
    /export/data/CNEs/mm10/annotation/MGI_Gene_Model_Coord.rpt \
    >>/export/data/CNEs/gff/mm10.gff
```

**FlyBase genes**

This is just for fruitfly.
GFF file from ftp://ftp.flybase.net/. 
Currently used file is named ```dmel-all-r5.51.gff.gz```.
The conversion can be done by the script ```cne/scripts/gbrowse_db/extract_genes_from_gff3.r```.

**RedFly**

This is just for fruitfly.

RedFly CRMs: GFF file from http://redfly.ccr.buffalo.edu/ (obtained by clicking “Search”, then “Download all CRM”) 
RedFly TFBSs: GFF file from http://redfly.ccr.buffalo.edu/ (obtained by clicking “Search”, then “Download all TFBS”)

```sh
perl /opt/www/cne/scripts/gbrowse_db/redfly2gff.pl \
	/export/data/goldenpath/dm3/assembly.2bit \
    /export/data/CNEs/dm3/annotation/redfly_CRM.gff CRM \
    >>/export/data/CNEs/gff/dm3.gff
perl /opt/www/cne/scripts/gbrowse_db/redfly2gff.pl \
	/export/data/goldenpath/dm3/assembly.2bit \
    /export/data/CNEs/dm3/annotation/redfly_TFBS.gff TFBS \
    >>/export/data/CNEs/gff/dm3.gff
```

**WormBase genes**

GFF file from ftp://ftp.wormbase.org/pub/wormbase/releases/WS220/species/c_elegans/.
File mapping gene and transcript IDs to symbols are downloaded from the same place.
The conversion can be done by the script ```cne/scripts/gbrowse_db/extract_genes_from_gff3.r```.

**CGD annotations**

The annotation for Candida is downloaded from [CGD](http://www.candidagenome.org/download/).
They provide ```tab``` format and ```gff3``` format. 
For display in GBrowse, we use the ```gff3``` format 
with lasted version [C_albicans_SC5314_version_A19-s01-m04-r03_features.gff](http://www.candidagenome.org/download/gff/C_albicans_SC5314/C_albicans_SC5314_version_A19-s01-m04-r03_features.gff).

Name conventions for candida:

| full name | short name |
| --------------- |--------------- | 
| C_albicans_SC5314A21 | SC5314A21 |
| C_dubliniensis_CD36 | CD36 |
| C_glabrata_CBS138 | CBS138 |
| C_parapsilosis_CDC317 |  CDC317 | 
| C_tropicalis_MYA-3404 | MYA3404 |
| C_guilliermondii_ATCC_6260 | ATCC6260 |
| C_albicans_WO-1 | WO1 |

The original reference genome from CGD is not repeat masked, 
i.e., all the bases in the upper case.
To mask the repeats into lower case, ```cne/scripts/gbrowse_db/candida/repeatmaskReference.r```
can be used to do that.

The evolutionary distance between these species can be found from 
[Broad Institute](http://www.broadinstitute.org/annotation/genome/candida_group/MultiHome.html).

**Aspergillus annotation**

The reference genome is downloaded from [Aspergillus Genome Database](http://www.aspgd.org/).

| full name | short name |
| --------- | ---------- |
| A_nidulans_FGSC_A4 | FGSCA4    |
| A_oryzae_RIB40  | RIB40     |
| A_fumigatus_Af293 | AF293    |
| A_niger_CBS_513_88    | CBS51388 |
| A_flavus_NRRL_3357 | NRRL3357 |

The evolutionary distance between these species can be found from 
[Broad Institute](http://www.broadinstitute.org/annotation/genome/aspergillus_group/MultiHome.html).

**S. accharomyces**
The reference genome from [SGD](http://yeastgenome.org/).

| full name | short name |
| --------- | ---------- |
| S_bayanus/623-6C | Sb6236C | 

**Pucciniales**
The reference genome is downloaded from Fungi Ensembl.

| full name | short name |
| --------- | ---------- |
| Puccinia graminis | ASM14992v1 |
| Puccinia triticina| ASM15152v1 |
| Melampsora larici-populina | GCA0002040551 |

**Magnaporthales**

The reference genome is downloaded from Fungi Ensembl.

| full name | short name |
| --------- | ---------- |
| Magnaporthe oryzae | MG8 |
| Gaeumannomyces graminis | GaegraminisV2 |
| Magnaporthe poae   | MagpoaeATCC64411V1 |

**Tsetse Fly**

The reference genome (Glossina morsitans) is from [Sanger Institute](ftp://ftp.sanger.ac.uk/pub/pathogens/Glossina/morsitans/Assemblies/tsetseGenome-v1.tar.gz).
The Glossina morsitans assembly was released by the Sanger Center in December 2010 and is in the process of being released. It consists of 13,807 supercontigs.
The annotation is from [VectorBase](https://www.vectorbase.org/downloads) with 
Gene set: GmorY1.3 25 Feb 2014.
Assembly ID: GmorY1.
CNEs are generated between Drosophila (dm3) vs GmorY1.


### Creating and loading GFF files
The annotations in the UCSC MySQL tables and 
the flat files must be converted to a GFF format suitable for 
import into the GBrowse annotation database. 
A number of scripts that carry out those conversions are available 
in the directory ```cne/scripts/gbrowse_db``` in the cne package. 
The script ```ucsc2gff.pl``` extracts data 
from UCSC annotation tables stored in a MySQL database and 
outputs GFF format suitable for import into the GBrowse database. 
The scripts that process flat files to generate GBrowse GFF files are listed in Table above.

```sh
perl /opt/www/cne/scripts/gbrowse_db/ucsc2gff.pl \
	-a /export/data/goldenpath/hg18/assembly.2bit \
    -d UCSC_hg18   assembly rmsk gap refGene \
    >hg18.gff
```

The script ```ucsc2gff.pl``` can also read genome assembly information from 2bit files 
and output a GFF feature for each sequence in the assembly 
(the chromosomes, or supercontigs if the genome is not assembled into chromosomes). 
These GFF lines are also needed by GBrowse.

In the same directory as the Perl conversion scripts, 
there is a bash script ```make_gff.sh``` that runs the conversion scripts 
to generate GFF files for all assemblies currently in Ancora. 
The script creates GFF files named by assembly id 
(e.g. hg19.gff, mm10.gff). 
It is recommended carry out the conversions by running this bash script 
instead of running the Perl scripts directly, 
as it also provides a record of which annotation is being used for each assembly. 
The bash script does not require any arguments, 
but it expects the UCSC databases with relevant annotation tables to be available on the local host, 
as well as relevant flat files to be present in the correct subdirectories. 
Please see the script for where to place the flat files so that it can find them.

When GFF files have been created 
they can be loaded into GBrowse annotation databases 
using the script ```load_gff.sh```. 
The script does not take any arguments, 
but prompts for a MySQL username and password. 
NOTE: use this script with caution as 
it erases all data stored in the GBrowse annotation databases 
before it loads the GFF files.

In details, for gff2 files,

```sh
bp_bulk_load_gff.pl -u $USER -p $PASS -c -d gbrowse_gff_hg19 --maxfeature 1000000000 /export/data/CNEs/gff/hg19.gff
```

And for gff3 files,
```sh
bp_seqfeature_load.pl -u $USER -p $PASS -d gbrowse_gff_hg19 -f -c genome1.gff3 genome1.fasta
```

### Creating a configuration file for GBrowse
GBrowse requires a configuration file for each assembly. 
We keep the live versions of these files in the repository under ```ancora/gbrowse2/conf```.
The files are named by assembly id and have the suffix .conf 
(e.g. hg19.conf for the current human assembly). 
The configuration files are fairly complex, 
so the easiest way to set up a file for a new assembly is 
to copy one written for a different assembly and modify it. 
The file hg19.conf begins with a few comment lines listing what should be changed 
when adapting a file for a new assembly. 
It is typically best to start from a file for a different assembly 
for the same species, if available, 
and otherwise from a file for a similar species. 
The reason for this is that the types of annotation available differs 
between species and some parameters, 
in particular the semantic zoom levels, 
depend on genome compactness. 
The general format of GBrowse configuration files is well described 
in the [GBrowse 2.0 HOWTO](http://gmod.org/wiki/GBrowse_2.0_HOWTO), 
which should be read before editing a configuration file.

The configuration of plugins that display CNE locations and densities 
deserves a more detailed description. 
These plugins are configured in three places in the configuration file:

(1) The line starting with “plugins =” specifies 
what plugins are to be loaded. 
For example, if CNEs from three comparisons are to be shown, 
this should be specified as:

```
plugins = CNEPlotInstance1 CNEPlotInstance2 CNEPlotInstance3
```

Each of the plugins listed on this line corresponds to a file 
in the subdirectory plugins under the configuration directory. 
The files CNEPlotInstance\*.pm in this directory are 
identical subclasses of the same base class ```Bio::Graphics::Browser2::Plugin::CNEPlot```.
This is a trick to be able to load several instances of the same plugin. 
There are currently thirteen CNE plugin files (up to CNEPlotInstance13.pm) 
in the plugins directory, 
so up to thirteen comparisons can be shown in any Ancora browser. 
If you wish to show more comparisons, 
simply add additional plugin files by copying one of the existing files and changing number on the first line in the new file.

(2) Parameters for each plugin are specified in the section following the comment

```
# Plugin config
```

One parameter must be specified for each plugin:
*	asm_id – Id for the other assembly in the comparison (e.g. hg19)

The following general parameters can also be specified:
*	asm_name – Name to be displayed for the other assembly in the comparison. 
	Usually the common name for the species (e.g. Human). 
    Defaults to asm_id.
*	window_size – Window size, in kb, to use for density computations. 
	Defaults to 50 kb.
*	max_score – Height of density plot. Defaults to 20.
*	nr_cne_sets – Number of CNE sets from this comparison 
	that the user can define and simultaneously view. Defaults to 3.
*	nr_graphs – Separate densities based on chromosome in other genome into this number of graphs. 
	Defaults to 1 (meaning don’t separate).
* 	rank_un_low – Give low priority to poorly assembled sequence 
	and alternative haplotypes when separating densities. Defaults to 0 (no). 
    Set this option to a true value to activate it.

In addition, the following parameters can also be used to define CNE sets, 
replacing i with a set number (1,2,3,…):
*	min_cne_idi - Similarity threshold for the CNE set. 
	This parameter specifies which CNE set to retrieve from the database, 
    so it must correspond to a CNE set present in the database. 
    Specified as: percent-identity/columns. 
    E.g.: “0.98/50” for 98% identity over 50 alignment columns. 
    Default: undefined.
*	min_cne_leni - Required shortest length for CNEs, in bp. Default: 50 bp.
*	show_cnesi - Show CNE locations. Default: 0 (no). Set to a true value to activate.
*	show_graphi - Show density plot. Default: 0 (no). Set to a true value to activate.

All of the above parameters except asm_id and asm_name can be modified by the user via the web interface. GBrowse uses cookies to save user settings.

(3) At the end of the configuration file, 
there is a track definition for each CNE plugin, e.g.:

```
[plugin:Mouse HCNEs]
citation = HCNEs conserved in mouse.
See <a href="/methods.html">Methods</a> for a description of how these were detected.
```

The order in which the track definitions are listed 
in the configuration file specifies the default order 
of the CNE tracks in the browser. 
This order can be changed by the user, 
who can also restore the default order by clicking the red “[Reset]” link in the genome browser.

### Making download files
Ancora has a download section where files with CNE locations (bed format) 
and CNE densities (wig format) can be downloaded. 
The files are located under ```/opt/www/ancora/html/downloads```. 
Download files can be created using the Perl script ```cne/scripts/cne_pipeline/make_download_files.pl``` 
in the cne package. 
Run the script without arguments for usage information. 
There is a wrapper bash script ```make_download_files.sh``` 
that runs the Perl script to create files for the assemblies currently 
shown in Ancora. 
It is therefore recommended to use this bash script 
for updating the Ancora download files. 
The bash script does not require any arguments.

The output files are placed in subdirectories 
created under the directory from which the scripts are run. 
The directory structure created is identical to that 
in the Ancora download section, 
so that the directories created by the scripts simply 
can be moved to the Ancora download directory ```/opt/www/ancora/html/downloads``` 
after checking that the files seem OK.

## The HTML pages and CGI scripts
The HTML pages are located under ```cne/ancora/html``` in the cne package. 
These files should be identical to the live HTML pages 
that we keep under ```/opt/www/ancora/html```. 
The style of the pages is defined in the stylesheet ```cne/ancora/html/main.css```. 
A stylesheet for GBrowse (```cne/ancora/html/gbrowse/gbrowse-cne.css```) 
is also part of the cne package, 
to give GBrowse a look and feel similar to the Ancora pages.

There are three CGI scripts for Ancora, 
located under ```/opt/www/ancora/cgi-bin``` (live version) 
and under ```cne/ancora/cgi-bin``` in the cne package:
*	cne_color_legend - Shows available color legends for chromosome colors. Ancora runs this script in a popup window opened from the CNE track configuration page.
*	cne_das_loader - Handles addition of DAS tracks to Ensembl.
*	cne_details - Shows details of CNEs or synteny block. Invoked if user clicks on one of these features in the genome browser.


**TO DO**
switch to gff3!!!!!



