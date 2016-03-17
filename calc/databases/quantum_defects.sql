BEGIN TRANSACTION;


DROP TABLE IF EXISTS `model_potential`;
CREATE TABLE `model_potential` ( `element` text,`L` int,`ac` real,`Z` int,`a1` real,`a2` real,`a3` real,`a4` real,`rc` real);

INSERT INTO `model_potential` VALUES('H',0,'0.0',1,'0.0','0.0','0.0','0.0','1.0');

-- Phys. Rev. A 49, 982 (1994)
INSERT INTO `model_potential` VALUES('Li',0,'0.1923',3,'2.47718079','1.84150932','-0.02169712','-0.11988362','0.61340824');
INSERT INTO `model_potential` VALUES('Li',1,'0.1923',3,'3.45414648','2.55151080','-0.21646561','-0.06990078','0.61566441');
INSERT INTO `model_potential` VALUES('Li',2,'0.1923',3,'2.51909839','2.43712450','0.32505524','0.10602430','2.34126273');
INSERT INTO `model_potential` VALUES('Li',3,'0.1923',3,'2.51909839','2.43712450','0.32505524','0.10602430','2.34126273');

-- Phys. Rev. A 49, 982 (1994)
INSERT INTO `model_potential` VALUES('Na',0,'0.9448',11,'4.82223117','2.45449865','-1.12255048','-1.42631393','0.45489422');
INSERT INTO `model_potential` VALUES('Na',1,'0.9448',11,'5.08382502','2.18226881','-1.19534623','-1.03142861','0.45798739');
INSERT INTO `model_potential` VALUES('Na',2,'0.9448',11,'3.53324124','2.48697936','-0.75688448','-1.27852357','0.71875312');
INSERT INTO `model_potential` VALUES('Na',3,'0.9448',11,'1.11056646','1.05458759','1.73203428','-0.09265696','28.6735059');

-- Phys. Rev. A 49, 982 (1994)
INSERT INTO `model_potential` VALUES('K',0,'5.3310',19,'3.56079437','1.83909642','-1.74701102','-1.03237313','0.83167545');
INSERT INTO `model_potential` VALUES('K',1,'5.3310',19,'3.65670429','1.67520788','-2.07416615','-0.89030421','0.85235381');
INSERT INTO `model_potential` VALUES('K',2,'5.3310',19,'4.12713694','1.79837462','-1.69935174','-0.98913582','0.83216907');
INSERT INTO `model_potential` VALUES('K',3,'5.3310',19,'1.42310446','1.27861156','4.77441476','-0.94829262','6.50294371');

-- Phys. Rev. A 49, 982 (1994)
INSERT INTO `model_potential` VALUES('Rb',0,'9.076',37,'3.69628474','1.64915255','-9.86069196','0.19579987','1.66242117');
INSERT INTO `model_potential` VALUES('Rb',1,'9.076',37,'4.44088978','1.92828831','-16.79597770','-0.81633314','1.50195124');
INSERT INTO `model_potential` VALUES('Rb',2,'9.076',37,'3.78717363','1.57027864','-11.6558897','0.52942835','4.86851938');
INSERT INTO `model_potential` VALUES('Rb',3,'9.076',37,'2.39848933','1.76810544','-12.0710678','0.77256589','4.79831327');

-- Phys. Rev. A 49, 982 (1994)
INSERT INTO `model_potential` VALUES('Cs',0,'15.6440',55,'3.49546309','1.47533800','-9.72143084','0.02629242','1.92046930');
INSERT INTO `model_potential` VALUES('Cs',1,'15.6440',55,'4.69366096','1.71398344','-24.65624280','-0.09543125','2.13383095');
INSERT INTO `model_potential` VALUES('Cs',2,'15.6440',55,'4.32466196','1.61365288','-6.70128850','-0.74095193','0.93007296');
INSERT INTO `model_potential` VALUES('Cs',3,'15.6440',55,'3.01048361','1.40000001','-3.20036138','0.00034538','1.99969677');


DROP TABLE IF EXISTS `rydberg_ritz`;
CREATE TABLE `rydberg_ritz` ( `element` text,`L` int,`J` real,`d0` real,`d2` real,`d4` real,`d6` real,`d8` real,`Ry` real);

INSERT INTO `rydberg_ritz` VALUES('H',0,'0.5','0.0','0.0','0.0','0.0','0.0','109737.31568525');

-- T. F. Gallagher, ``Rydberg Atoms'', Cambridge University Press (2005), ISBN: 978-0-52-102166-1
INSERT INTO `rydberg_ritz` VALUES('Li',0,'0.5','0.399468','0.030233','-0.0028','0.0115','0.0','109728.64');
INSERT INTO `rydberg_ritz` VALUES('Li',1,'0.5','0.47263','-0.02613','0.0221','-0.0683','0.0','109728.64');
INSERT INTO `rydberg_ritz` VALUES('Li',1,'1.5','0.47263','-0.02613','0.0221','-0.0683','0.0','109728.64');
INSERT INTO `rydberg_ritz` VALUES('Li',2,'1.5','0.002129','-0.01491','0.1759','-0.8507','0.0','109728.64');
INSERT INTO `rydberg_ritz` VALUES('Li',2,'2.5','0.002129','-0.01491','0.1759','-0.8507','0.0','109728.64');
INSERT INTO `rydberg_ritz` VALUES('Li',3,'2.5','0.000305','-0.00126','0.0','0.0','0.0','109728.64');
INSERT INTO `rydberg_ritz` VALUES('Li',3,'3.5','0.000305','-0.00126','0.0','0.0','0.0','109728.64');

-- T. F. Gallagher, ``Rydberg Atoms'', Cambridge University Press (2005), ISBN: 978-0-52-102166-1
INSERT INTO `rydberg_ritz` VALUES('Na',0,'0.5','1.347969','0.06137','0.0','0.0','0.0','109734.69');
INSERT INTO `rydberg_ritz` VALUES('Na',1,'0.5','0.855424','0.1222','0.0','0.0','0.0','109734.69');
INSERT INTO `rydberg_ritz` VALUES('Na',1,'1.5','0.854608','0.122','0.0','0.0','0.0','109734.69');
INSERT INTO `rydberg_ritz` VALUES('Na',2,'1.5','0.015543','-0.08535','0.7958','-4.0513','0.0','109734.69');
INSERT INTO `rydberg_ritz` VALUES('Na',2,'2.5','0.015543','-0.08535','0.7958','-4.0513','0.0','109734.69');
INSERT INTO `rydberg_ritz` VALUES('Na',3,'2.5','0.001663','-0.0098','0.0','0.0','0.0','109734.69');

-- T. F. Gallagher, ``Rydberg Atoms'', Cambridge University Press (2005), ISBN: 978-0-52-102166-1
INSERT INTO `rydberg_ritz` VALUES('K',0,'0.5','2.180197','0.136','0.0759','0.117','-0.206','109735.774');
INSERT INTO `rydberg_ritz` VALUES('K',1,'0.5','1.713892','0.2332','0.16137','0.5345','-0.234','109735.774');
INSERT INTO `rydberg_ritz` VALUES('K',1,'1.5','1.710848','0.2354','0.11551','1.105','-2.0356','109735.774');
INSERT INTO `rydberg_ritz` VALUES('K',2,'1.5','0.27697','-1.0249','-0.709174','11.839','-26.689','109735.774');
INSERT INTO `rydberg_ritz` VALUES('K',2,'2.5','0.277158','-1.0256','-0.59201','10.0053','-19.0244','109735.774');
INSERT INTO `rydberg_ritz` VALUES('K',3,'2.5','0.010098','-0.100224','1.56334','-12.6851','0.0','109735.774');
INSERT INTO `rydberg_ritz` VALUES('K',3,'3.5','0.010098','-0.100224','1.56334','-12.6851','0.0','109735.774');

-- [1] Phys. Rev. A 83, 052515 (2011) - Rb87
-- [2] Phys. Rev. A 67, 052502 (2003) - Rb85
-- [3] Phys. Rev. A 74, 054502 (2006) - Rb85
INSERT INTO `rydberg_ritz` VALUES('Rb',0,'0.5','3.1311807','0.1787','0.0','0.0','0.0','109736.62301604665'); -- [1]
INSERT INTO `rydberg_ritz` VALUES('Rb',1,'0.5','2.6548849','0.29','0.0','0.0','0.0','109736.605'); -- [2]
INSERT INTO `rydberg_ritz` VALUES('Rb',1,'1.5','2.6416737','0.295','0.0','0.0','0.0','109736.605'); -- [2]
INSERT INTO `rydberg_ritz` VALUES('Rb',2,'1.5','1.3480948','-0.6054','0.0','0.0','0.0','109736.62301604665'); -- [1]
INSERT INTO `rydberg_ritz` VALUES('Rb',2,'2.5','1.3464622','-0.594','0.0','0.0','0.0','109736.62301604665'); -- [1]
INSERT INTO `rydberg_ritz` VALUES('Rb',3,'2.5','0.0165192','-0.085','0.0','0.0','0.0','109736.605'); -- [3]
INSERT INTO `rydberg_ritz` VALUES('Rb',3,'3.5','0.0165437','-0.086','0.0','0.0','0.0','109736.605'); -- [3]

/*-- T. F. Gallagher, ``Rydberg Atoms'', Cambridge University Press (2005), ISBN: 978-0-52-102166-1
INSERT INTO `rydberg_ritz` VALUES('Rb',0,'0.5','3.13109','0.204','-1.8','0.0','0.0','109736.605');
INSERT INTO `rydberg_ritz` VALUES('Rb',1,'0.5','2.65456','0.388','-7.904','116.437','-405.907','109736.605');
INSERT INTO `rydberg_ritz` VALUES('Rb',1,'1.5','2.64145','0.33','-0.97495','14.6001','-44.7265','109736.605');
INSERT INTO `rydberg_ritz` VALUES('Rb',2,'1.5','1.347157','-0.59553','-1.50517','-2.4206','19.736','109736.605');
INSERT INTO `rydberg_ritz` VALUES('Rb',2,'2.5','1.347157','-0.59553','-1.50517','-2.4206','19.736','109736.605');
INSERT INTO `rydberg_ritz` VALUES('Rb',3,'2.5','0.016312','-0.064007','-0.36005','3.239','0.0','109736.605');
INSERT INTO `rydberg_ritz` VALUES('Rb',3,'3.5','0.016312','-0.064007','-0.36005','3.239','0.0','109736.605');*/

-- [1] Phys. Rev. A 93, 013424 (2016)
-- [2] Phys. Rev. A 26, 2733 (1982)
INSERT INTO `rydberg_ritz` VALUES('Cs',0,'0.5','4.0493532','0.2391','0.06','11','-209','109736.8627339'); -- [1]
INSERT INTO `rydberg_ritz` VALUES('Cs',1,'0.5','3.5915871','0.36273','0.0','0.0','0.0','109736.8627339'); -- [1]
INSERT INTO `rydberg_ritz` VALUES('Cs',1,'1.5','3.5590676','0.37469','0.0','0.0','0.0','109736.8627339'); -- [1]
INSERT INTO `rydberg_ritz` VALUES('Cs',2,'1.5','2.475365','0.5554','0.0','0.0','0.0','109736.85999132106'); -- [2]
INSERT INTO `rydberg_ritz` VALUES('Cs',2,'2.5','2.4663144','0.01381','−0.392','−1.9','0.0','109736.8627339'); -- [1]
INSERT INTO `rydberg_ritz` VALUES('Cs',3,'2.5','0.033392','-0.191','0.0','0.0','0.0','109736.85999132106'); -- [2]
INSERT INTO `rydberg_ritz` VALUES('Cs',3,'3.5','0.033537','-0.191','0.0','0.0','0.0','109736.85999132106'); -- [2]

/*-- Phys. Rev. A 26, 2733 (1982)
INSERT INTO `rydberg_ritz` VALUES('Cs',0,'0.5','4.049325','0.2462','0.0','0.0','0.0','109736.85999132106');
INSERT INTO `rydberg_ritz` VALUES('Cs',1,'0.5','3.591556','0.3714','0.0','0.0','0.0','109736.85999132106');
INSERT INTO `rydberg_ritz` VALUES('Cs',1,'1.5','3.559058','0.374','0.0','0.0','0.0','109736.85999132106');
INSERT INTO `rydberg_ritz` VALUES('Cs',2,'1.5','2.475365','0.5554','0.0','0.0','0.0','109736.85999132106');
INSERT INTO `rydberg_ritz` VALUES('Cs',2,'2.5','2.466210','0.067','0.0','0.0','0.0','109736.85999132106');
INSERT INTO `rydberg_ritz` VALUES('Cs',3,'2.5','0.033392','-0.191','0.0','0.0','0.0','109736.85999132106');
INSERT INTO `rydberg_ritz` VALUES('Cs',3,'3.5','0.033537','-0.191','0.0','0.0','0.0','109736.85999132106');*/

/*-- T. F. Gallagher, ``Rydberg Atoms'', Cambridge University Press (2005), ISBN: 978-0-52-102166-1
INSERT INTO `rydberg_ritz` VALUES('Cs',0,'0.5','4.049325','0.246','0.0','0.0','0.0','109736.86');
INSERT INTO `rydberg_ritz` VALUES('Cs',1,'0.5','3.591556','0.3714','0.0','0.0','0.0','109736.86');
INSERT INTO `rydberg_ritz` VALUES('Cs',1,'1.5','3.559058','0.374','0.0','0.0','0.0','109736.86');
INSERT INTO `rydberg_ritz` VALUES('Cs',2,'1.5','2.475365','0.5554','0.0','0.0','0.0','109736.86');
INSERT INTO `rydberg_ritz` VALUES('Cs',2,'2.5','2.46621','0.0167','0.0','0.0','0.0','109736.86');
INSERT INTO `rydberg_ritz` VALUES('Cs',3,'2.5','0.033392','-0.191','0.0','0.0','0.0','109736.86');
INSERT INTO `rydberg_ritz` VALUES('Cs',3,'3.5','0.033537','-0.191','0.0','0.0','0.0','109736.86');*/


COMMIT;

