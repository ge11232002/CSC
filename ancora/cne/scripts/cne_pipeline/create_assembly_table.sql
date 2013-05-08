-- MySQL dump 10.10
--
-- Host: localhost    Database: cne
-- ------------------------------------------------------
-- Server version	5.0.22

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `assembly`
--

DROP TABLE IF EXISTS `assembly`;
CREATE TABLE `assembly` (
  `assembly_id` varchar(255) NOT NULL default '',
  `assembly_name` varchar(255) NOT NULL default '',
  `organism_common` varchar(255) NOT NULL default '',
  `organism_latin` varchar(255) NOT NULL default '',
  `ensembl_ver` varchar(255) NOT NULL default '',
  `default_ensembl_loc` varchar(255) NOT NULL default '',
  `ucsc_db` varchar(255) NOT NULL default '',
  PRIMARY KEY  (`assembly_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

--
-- Dumping data for table `assembly`
--


/*!40000 ALTER TABLE `assembly` DISABLE KEYS */;
LOCK TABLES `assembly` WRITE;
INSERT INTO `assembly` VALUES ('hg18','NCBI Build 36','human','Homo sapiens','','11:31766034-31797085','UCSC_HS_MAR06'),('danRer4','Zv6','zebrafish','Danio rerio','apr2007','3:51040000-51139999','UCSC_DR_MAR06'),('mm8','NCBI Build 36','mouse','Mus musculus','aug2007','6:98800000-99299999','UCSC_MM_MAR06'),('monDom4','monDom4','opossum','Monodelphis domestica','','',''),('galGal3','v2.1','chicken','Gallus gallus','','',''),('xenTro2','v4.1','frog','Xenopus tropicalis','','',''),('dm2','BDGP Release 4','','Drosophila melanogaster','','','UCSC_DM_APR04'),('dm3','BDGP Release 5','','Drosophila melanogaster','','','UCSC_dm3'),('tetNig1','V7','tetraodon','Tetraodon nigroviridis','','',''),('fr1','v3.0','fugu','Takifugu rubripes','','','UCSC_fr1'),('fr2','v4.0','fugu','Takifugu rubripes','','',''),('gasAcu1','v1.0','stickleback','Gasterosteus aculeatus','','',''),('oryLat1','v1.0','medaka','Oryzias latipes','','',''),('mm9','NCBI Build 37','mouse','Mus musculus','','',''),('danRer5','Zv7','zebrafish','Danio rerio','','7:8541098-8656549','UCSC_danRer5'),('droAna3','Feb 2006','','Drosophila ananassae','','',''),('dp4','Release 2','','Drosophila pseudoobscura','','',''),('droVir3','Feb 2006','','Drosophila virilis','','',''),('droMoj3','Feb 2006','','Drosophila mojavensis','','','');
UNLOCK TABLES;
/*!40000 ALTER TABLE `assembly` ENABLE KEYS */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

