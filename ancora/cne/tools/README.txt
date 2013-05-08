Several of the C programs here use libraries from the the UCSC Genome Browser
source code. To compile these programs, you need to set the environment variable
JKSRC to point to the location of the UCSC source. On our servers shire and
angband this is /usr/local/jksrc/kent/src, so you would do:

export JKSRC=/usr/local/jksrc/kent/src

Then go into the directory for the program you wish to compile and type

make

Executables will be put in your ~/bin directory.

If you are installing on a machine where the UCSC source is not available,
you first need to obtain it and compile the libraries
The source code can be downloaded from
http://genome.ucsc.edu/FAQ/FAQlicense
It comes with a README file explaining how to compile it. 
Instructions are also available at
http://genome.ucsc.edu/admin/jk-install.html
On our server shire this currently boils down to:

export MACHTYPE=x86_64
export MYSQLINC=/usr/include/mysql
export MYSQLLIBS= "-L/usr/lib64/mysql -lmysqlclient -lz -lcrypt -lnsl -lm -L/usr/lib64 -lssl -lcrypto"
mkdir ~/bin/$MACHTYPE
cd kent/src
(make clean)
make libs

Note that the the you should set the MACHTYPE, MYSQLINC and MYSQLLIBS
variables according to the configuration of your system. See the README file
in the UCSC source package for more details.

