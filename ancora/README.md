#Ancora technical documentation
First version by Pär Engström, April 2008. This version by Ge Tan, May 2013. This document was produced by the [GitHub Flavored Markdown](https://help.github.com/articles/github-flavored-markdown "Markdown").

This document describes the implementation the Ancora web resource. For a general description of Ancora and its user interface, see the [Ancora publication](http://genomebiology.com/2008/9/2/R34 "Ancora publication").

**Table of contents**

*    [Installation](#installation)
	*	[Software dependencies](#dependencies)
	*	[GBrowse2](#GBrowse2)
	*	[UCSC Genome Browser source and utilities](#UCSC)
	*	[ProServer DAS server](#DAS)
	
<h2 id="installation">Installation</h2>
 This documentation focuses on the Ancora installation on Olifant at csc. Olifant is running the CentOS release 5.8 (Final). However, ancora is supposed to run on other Linux/Unix distribution without too much difficulty.

<h3 id="dependencies">Software dependecies</h3>
The following softwares are essential for the whole implementation.
 
 *	A httpd server (Currently [Apache](http://httpd.apache.org/docs/2.2/ "Apache") 2.2.3 is used on olifant)
 *	[MySQL](http://dev.mysql.com/downloads/mysql/5.0.html "MySQL") client if you have a dedicated MySQL server. if not, MySQL server is alse required. (MySQL 5.0.95 is used on olifant)
 *	[BioPerl](http://www.bioperl.org/wiki/Installing_BioPerl_on_Unix "BioPerl"). Installation is quite straightforward following the instructions. (BioPerl 1.6.1 is used on olifant.)

<h3 id="GBrowse2">GBrowse2</h3>
This implementation of Ancora utilizes the latest version of GBrowse 2.54. GBrowse 2.X is a complete rewrite of GBrowse 1.X version. There are several advantages compared to the GBrowse 1.X. See the [details](http://gmod.org/wiki/GBrowse_2.0_HOWTO#Installation_from_Source_Code "GBrowse").

Since the GBrowse 2.X is perl-based and all the modules are hosted on [CPAN](http://search.cpan.org/~lds/GBrowse-2.54/), the easiest way to install GBrowse 2 is using the standard Perl module building procedure.








