-- MySQL dump 10.11
--
-- Host: localhost    Database: cne
-- ------------------------------------------------------
-- Server version	5.0.95

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
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
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
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `assembly`
--

LOCK TABLES `assembly` WRITE;
/*!40000 ALTER TABLE `assembly` DISABLE KEYS */;
INSERT INTO `assembly` VALUES ('hg18','NCBI Build 36','human','Homo sapiens','nov2005','11:31766034-31797085','UCSC_HS_MAR06'),('danRer4','Zv6','zebrafish','Danio rerio','apr2007','3:51040000-51139999','UCSC_DR_MAR06'),('galGal4','v4.0','chicken','Gallus gallus','nov2011','Z:9897818-9987257',''),('monDom4','monDom4','opossum','Monodelphis domestica','','',''),('galGal3','v2.1','chicken','Gallus gallus','may2006','Z:9907346-9929706',''),('xenTro2','v4.1','frog','Xenopus tropicalis','aug2005','GL172646.1:2369612-2552500',''),('dm2','BDGP Release 4','fruitfly','Drosophila melanogaster','','2L:21650001-21700000','UCSC_DM_APR04'),('dm3','BDGP Release 5','fruitfly','Drosophila melanogaster','apr2006','2L:21650001-21700000','UCSC_dm3'),('tetNig1','V7','tetraodon','Tetraodon nigroviridis','','',''),('fr1','v3.0','fugu','Takifugu rubripes','','','UCSC_fr1'),('fr2','v4.0','fugu','Takifugu rubripes','jun2005','scaffold_611:1-67046',''),('gasAcu1','v1.0','stickleback','Gasterosteus aculeatus','feb2006','groupI:18287526-18588610',''),('oryLat1','v1.0','medaka','Oryzias latipes','oct2005','14:31116527-31129337',''),('mm9','NCBI Build 37','mouse','Mus musculus','apr2007','6:98800000-99299999',''),('danRer5','Zv7','zebrafish','Danio rerio','','7:8541098-8656549','UCSC_danRer5'),('droAna3','Feb 2006','Drosophila ananassae','Drosophila ananassae','','',''),('dp4','Release 2','','Drosophila pseudoobscura','','',''),('droVir3','Feb 2006','','Drosophila virilis','','',''),('droMoj3','Feb 2006','','Drosophila mojavensis','','',''),('hg19','NCBI Build 37','human','Homo sapiens','feb2009','11:31766034-31797085',''),('monDom5','monDom5','opossum','Monodelphis domestica','oct2006',' 6:36217906-36423052',''),('mm10','NCBI Build 38','mouse','Mus musculus','jan2012',' 6:98800000-99299999',''),('xenTro3','v4.2','frog','Xenopus tropicalis','nov2009','GL172646.1:2369612-2552500',''),('tetNig2','v8.0','tetraodon','Tetraodon nigroviridis','mar2007','1:783604-1032531',''),('fr3','v5.0','fugu','Takifugu rubripes','','',''),('oryLat2','v2.0','medaka','Oryzias latipes','oct2005','14:31116527-31129337',''),('danRer7','Zv9','zebrafish','Danio rerio','apr2010','10:10138322-10349251',''),('danRer6','Zv8','zebrafish','Danio rerio','dec2008','10:10138322-10349251',''),('equCab2','EquCab2','horse','Equus caballus','sep2007','X:14898464-14976430',''),('canFam3','CanFam3.1','dog','Canis lupus familiaris','sep2011','10:13152725-13327708',''),('rn5','Rnor_5.0','rat','Rattus norvegicus','mar2012','5:162132007-162167596',''),('rn4','v3.4','rat','Rattus norvegicus','dec2004','5:162132007-162167596',''),('canFam2','canFam2','dog','Canis lupus familiaris','may2006','X:42152150-42169554',''),('ce10','WS220','Caenorhabditis elegans','Caenorhabditis elegans','','',''),('ornAna1','5.0.1','platypus','Ornithorhynchus anatinus','dec2005','Ultra274:1048851-1231454',''),('cb3','cb3','C. briggsae','C. briggsae','','',''),('ce6','WS190','Caenorhabditis elegans','Caenorhabditis elegans','','',''),('ce4','WS170','Caenorhabditis elegans','Caenorhabditis elegans','','',''),('caeRem2','v1.0','C. remanei','C. remanei','','',''),('caePb1','v4.0','C. brenneri','C. brenneri','','',''),('anoCar2','v2.0','lizard','Anolis carolinensis','may2010','GL343273.1:846138-853994',''),('droAna2','Aug 2005','Drosophila ananassae','Drosophila ananassae','','',''),('spur2.5.16','spur2.5.16','Purple sea urchin','Strongylocentrotus purpuratus','','',''),('lytVar0.4','lytVar04','Green sea urchin','Lytechinus variegatus','','',''),('petMar2','WUGSC 7.0','lamprey','Petromyzon marinus','jan2011','GL476598:172728-528162','');
/*!40000 ALTER TABLE `assembly` ENABLE KEYS */;
UNLOCK TABLES;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2013-06-13 15:58:33
