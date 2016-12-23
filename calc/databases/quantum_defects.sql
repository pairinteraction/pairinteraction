-- Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
--
-- Licensed under the Apache License, Version 2.0 (the "License");
-- you may not use this file except in compliance with the License.
-- You may obtain a copy of the License at
--
--     http://www.apache.org/licenses/LICENSE-2.0
--
-- Unless required by applicable law or agreed to in writing, software
-- distributed under the License is distributed on an "AS IS" BASIS,
-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
-- See the License for the specific language governing permissions and
-- limitations under the License.

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

-- [1] Phys. Rev. A 34, 2889 (1986) (Li 7)
-- [2] T. F. Gallagher, ``Rydberg Atoms'', Cambridge University Press (2005), ISBN: 978-0-52-102166-1 
-- [3] Johansson I 1958 Ark. Fysik 15 169
INSERT INTO `rydberg_ritz` VALUES('Li',0,'0.5','0.3995101','0.029','0.0','0.0','0.0','109728.64'); -- [1]
INSERT INTO `rydberg_ritz` VALUES('Li',1,'0.5','0.0471780','-0.024','0.0','0.0','0.0','109728.64'); -- [1]
INSERT INTO `rydberg_ritz` VALUES('Li',1,'1.5','0.0471665','-0.024','0.0','0.0','0.0','109728.64'); -- [1]
INSERT INTO `rydberg_ritz` VALUES('Li',2,'1.5','0.002129','-0.01491','0.1759','-0.8507','0.0','109728.64'); -- [2,3]
INSERT INTO `rydberg_ritz` VALUES('Li',2,'2.5','0.002129','-0.01491','0.1759','-0.8507','0.0','109728.64'); -- [2,3]
INSERT INTO `rydberg_ritz` VALUES('Li',3,'2.5','0.000305','-0.00126','0.0','0.0','0.0','109728.64'); -- [2,3]
INSERT INTO `rydberg_ritz` VALUES('Li',3,'3.5','0.000305','-0.00126','0.0','0.0','0.0','109728.64'); -- [2,3]
INSERT INTO `rydberg_ritz` VALUES('Li',4,'3.5','0.0','0.0','0.0','0.0','0.0','109728.64');

-- [1] Phys. Rev. A 45, 4720 (1992)
-- [2] Quantum Electron. 25 914 (1995)
-- [3] J. Phys. B: At. Mol. Opt. Phys. 30 2345 (1997)
INSERT INTO `rydberg_ritz` VALUES('Na',0,'0.5','1.34796938','0.0609892','0.0196743','-0.001045','0.0','109734.69'); -- [1,2]
INSERT INTO `rydberg_ritz` VALUES('Na',1,'0.5','0.85544502','0.112067','0.0479','0.0457','0.0','109734.69'); -- [2] 
INSERT INTO `rydberg_ritz` VALUES('Na',1,'1.5','0.85462615','0.112344','0.0497','0.0406','0.0','109734.69'); -- [2] 
INSERT INTO `rydberg_ritz` VALUES('Na',2,'1.5','0.014909286','-0.042506','0.00840','0.0','0.0','109734.69'); -- [2,3]
INSERT INTO `rydberg_ritz` VALUES('Na',2,'2.5','0.01492422','-0.042585','0.00840','0.0','0.0','109734.69'); -- [2,3]
INSERT INTO `rydberg_ritz` VALUES('Na',3,'2.5','0.001632977','-0.0069906','0.00423','0.0','0.0','109734.69'); -- [3]
INSERT INTO `rydberg_ritz` VALUES('Na',3,'3.5','0.001630875','-0.0069824','0.00352','0.0','0.0','109734.69'); -- [3]
INSERT INTO `rydberg_ritz` VALUES('Na',4,'3.5','0.00043825','-0.00283','0.0','0.0','0.0','109734.69'); -- [3]
INSERT INTO `rydberg_ritz` VALUES('Na',4,'4.5','0.00043740','-0.00297','0.0','0.0','0.0','109734.69'); -- [3]
INSERT INTO `rydberg_ritz` VALUES('Na',5,'4.5','0.00016114','-0.00185','0.0','0.0','0.0','109734.69'); -- [3]
INSERT INTO `rydberg_ritz` VALUES('Na',5,'5.5','0.00015796','-0.00148','0.0','0.0','0.0','109734.69'); -- [3]
INSERT INTO `rydberg_ritz` VALUES('Na',6,'5.5','0.0','0.0','0.0','0.0','0.0','109734.69');

-- [1] Phys. Scr. 27, 300 (1983) 
-- [2] Opt. Commun. 39, 370 (1981)
-- [3] Ark. Fys., 10 p.583 (1956) 
INSERT INTO `rydberg_ritz` VALUES('K',0,'0.5','2.180197','0.136','0.0759','0.117','-0.206','109735.774'); -- [1,2]
INSERT INTO `rydberg_ritz` VALUES('K',1,'0.5','1.713892','0.2332','0.16137','0.5345','-0.234','109735.774'); -- [1]
INSERT INTO `rydberg_ritz` VALUES('K',1,'1.5','1.710848','0.2354','0.11551','1.105','-2.0356','109735.774'); -- [1]
INSERT INTO `rydberg_ritz` VALUES('K',2,'1.5','0.27697','-1.0249','-0.709174','11.839','-26.689','109735.774'); -- [1,2]
INSERT INTO `rydberg_ritz` VALUES('K',2,'2.5','0.277158','-1.0256','-0.59201','10.0053','-19.0244','109735.774'); -- [1,2]
INSERT INTO `rydberg_ritz` VALUES('K',3,'2.5','0.010098','-0.100224','1.56334','-12.6851','0.0','109735.774'); -- [1,3]
INSERT INTO `rydberg_ritz` VALUES('K',3,'3.5','0.010098','-0.100224','1.56334','-12.6851','0.0','109735.774'); -- [1,3]
INSERT INTO `rydberg_ritz` VALUES('K',4,'3.5','0.0','0.0','0.0','0.0','0.0','109735.774');

-- [1] Phys. Rev. A 83, 052515 (2011) - Rb87
-- [2] Phys. Rev. A 67, 052502 (2003) - Rb85
-- [3] Phys. Rev. A 74, 054502 (2006) - Rb85
-- [4] Phys. Rev. A 74, 062712 (2006) - Rb85
INSERT INTO `rydberg_ritz` VALUES('Rb',0,'0.5','3.1311807','0.1787','0.0','0.0','0.0','109736.62301604665'); -- [1]
INSERT INTO `rydberg_ritz` VALUES('Rb',1,'0.5','2.6548849','0.29','0.0','0.0','0.0','109736.605'); -- [2]
INSERT INTO `rydberg_ritz` VALUES('Rb',1,'1.5','2.6416737','0.295','0.0','0.0','0.0','109736.605'); -- [2]
INSERT INTO `rydberg_ritz` VALUES('Rb',2,'1.5','1.3480948','-0.6054','0.0','0.0','0.0','109736.62301604665'); -- [1]
INSERT INTO `rydberg_ritz` VALUES('Rb',2,'2.5','1.3464622','-0.594','0.0','0.0','0.0','109736.62301604665'); -- [1]
INSERT INTO `rydberg_ritz` VALUES('Rb',3,'2.5','0.0165192','-0.085','0.0','0.0','0.0','109736.605'); -- [3]
INSERT INTO `rydberg_ritz` VALUES('Rb',3,'3.5','0.0165437','-0.086','0.0','0.0','0.0','109736.605'); -- [3]
INSERT INTO `rydberg_ritz` VALUES('Rb',4,'3.5','0.004','0.0','0.0','0.0','0.0','109736.605'); -- [4]
INSERT INTO `rydberg_ritz` VALUES('Rb',4,'4.5','0.004','0.0','0.0','0.0','0.0','109736.605'); -- [4]
INSERT INTO `rydberg_ritz` VALUES('Rb',5,'4.5','0.0','0.0','0.0','0.0','0.0','109736.605');

-- [1] Phys. Rev. A 93, 013424 (2016)
-- [2] Phys. Rev. A 26, 2733 (1982)
INSERT INTO `rydberg_ritz` VALUES('Cs',0,'0.5','4.0493532','0.2391','0.06','11','-209','109736.8627339'); -- [1]
INSERT INTO `rydberg_ritz` VALUES('Cs',1,'0.5','3.5915871','0.36273','0.0','0.0','0.0','109736.8627339'); -- [1]
INSERT INTO `rydberg_ritz` VALUES('Cs',1,'1.5','3.5590676','0.37469','0.0','0.0','0.0','109736.8627339'); -- [1]
INSERT INTO `rydberg_ritz` VALUES('Cs',2,'1.5','2.475365','0.5554','0.0','0.0','0.0','109736.85999132106'); -- [2]
INSERT INTO `rydberg_ritz` VALUES('Cs',2,'2.5','2.4663144','0.01381','−0.392','−1.9','0.0','109736.8627339'); -- [1]
INSERT INTO `rydberg_ritz` VALUES('Cs',3,'2.5','0.033392','-0.191','0.0','0.0','0.0','109736.85999132106'); -- [2]
INSERT INTO `rydberg_ritz` VALUES('Cs',3,'3.5','0.033537','-0.191','0.0','0.0','0.0','109736.85999132106'); -- [2]
INSERT INTO `rydberg_ritz` VALUES('Cs',4,'3.5','0.0','0.0','0.0','0.0','0.0','109736.85999132106');

COMMIT;

