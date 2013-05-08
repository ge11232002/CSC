--
-- Table structure for table 'HSP'
--

CREATE TABLE HSP (
  hsp_id int(10) unsigned NOT NULL auto_increment,
  mapping_id mediumint(8) unsigned NOT NULL default '0',
  qStart int(10) unsigned NOT NULL default '0',
  qEnd int(10) unsigned NOT NULL default '0',
  tStart int(10) unsigned NOT NULL default '0',
  tEnd int(10) unsigned NOT NULL default '0',
  blockSize int(10) unsigned NOT NULL default '0',
  qSeq mediumtext,
  tSeq mediumtext,
  PRIMARY KEY  (hsp_id)
) TYPE=MyISAM;

--
-- Table structure for table 'MAPPING'
--

CREATE TABLE MAPPING (
  mapping_id mediumint(8) unsigned NOT NULL auto_increment,
  bin smallint(5) unsigned NOT NULL default '0',
  target_db varchar(64) default NULL,
  matches int(10) unsigned NOT NULL default '0',
  misMatches int(10) unsigned NOT NULL default '0',
  repMatches int(10) unsigned NOT NULL default '0',
  nCount int(10) unsigned NOT NULL default '0',
  qNumInsert int(10) unsigned NOT NULL default '0',
  qBaseInsert int(10) unsigned NOT NULL default '0',
  tNumInsert int(10) unsigned NOT NULL default '0',
  tBaseInsert int(10) unsigned NOT NULL default '0',
  strand char(2) NOT NULL default '',
  qName varchar(40) NOT NULL default '',
  qType enum('EST','mRNA', 'unknown') NOT NULL default 'unknown',
  qSize int(10) unsigned NOT NULL default '0',
  qStart int(10) unsigned NOT NULL default '0',
  qEnd int(10) unsigned NOT NULL default '0',
  tName varchar(40) NOT NULL default '',
  tSize int(10) unsigned NOT NULL default '0',
  tStart int(10) unsigned NOT NULL default '0',
  tEnd int(10) unsigned NOT NULL default '0',
  PRIMARY KEY  (mapping_id)
) TYPE=MyISAM;

--
-- Table structure for table 'CONFIG'
--

CREATE TABLE CONFIG (
  id varchar(80) NOT NULL,
  value text NOT NULL
) TYPE=MyISAM;

INSERT CONFIG VALUES ('use_binning','1');
INSERT CONFIG VALUES ('use_mrnainfo','1');

---
--- Table strcuture for table 'LIB_ORIANN_ACCURACY'
---

CREATE TABLE LIB_ORIANN_ACCURACY (
  libid int unsigned not null,
  nr_ests_tot mediumint unsigned not null,
  nr_ests_splori mediumint unsigned not null,
  nr_ests_annori mediumint unsigned not null,
  correct_ppt smallint unsigned not null,
  PRIMARY KEY(libid),
  KEY(nr_ests_annori),
  KEY(correct_ppt)
);
