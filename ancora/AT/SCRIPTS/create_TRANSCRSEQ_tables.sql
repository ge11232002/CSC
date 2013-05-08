

CREATE TABLE SEQ (
  acc char(12) NOT NULL default '',
  version tinyint(2) unsigned NOT NULL default 0,
  division char(3) NOT NULL default '',
  replaced tinyint(1) unsigned NOT NULL default 0,
  locus varchar(20) NOT NULL default '',
  type enum('EST','mRNA', 'unknown') NOT NULL default 'unknown',
  organism enum('Homo sapiens', 'Mus musculus', 'Rattus norvegicus', 'unknown')
    NOT NULL default 'unknown',
  data_source_id int unsigned,
  description text NOT NULL default '',
  seq mediumtext NOT NULL default ''
) TYPE=MyISAM max_rows=20000000 avg_row_length=700;


CREATE TABLE DATA_SOURCE (
  data_source_id int unsigned NOT NULL auto_increment primary key,
  name varchar(80) NOT NULL default '',
  version varchar(40) NOT NULL default 0,
  component varchar(80) NOT NULL default '',
  load_start_time datetime,
  load_end_time datetime,
  unique key (name, version, component)
) TYPE=MyISAM;


CREATE TABLE CONFIG (
  id varchar(80) NOT NULL,
  value text NOT NULL
) TYPE=MyISAM;


CREATE TABLE RIKEN_ID (
  clone_id char(10) NOT NULL default '',
  rearray_id char(10) NOT NULL default '',
  seq_id mediumint unsigned,
  acc char(12) NOT NULL default ''
) TYPE=MyISAM;
