// * CONSTANTS
let margin = { top: 40, right: 200, bottom: 180, left: 130 };
let w = 1200 - margin.left - margin.right,
  h = 900 - margin.top - margin.bottom;

// * DATA
// let dataset = interactions;
let dataset = {
  '3sn6' : [
  ['R','F',139,'34.51x51','A','V',217, ["hydrophobic", "van-der-waals"]],
  ['R','Q',229,'5.68x68','A','R',385, ["hydrophobic", "polar-backbone-sidechain", "polar-sidechain-sidechain", "h-bond"]],
  ['R','R',63,'12.49x49','B','D',312, ["polar-backbone-sidechain"]],
  ['R','E',225,'5.64x64','A','Q',384, ["polar-sidechain-sidechain", "van-der-waals"]],
  ['R','E',62,'12.48x48','B','D',312, ["polar-sidechain-sidechain"]],
  ['R','R',239,'-','A','D',343, ["polar-sidechain-backbone"]],
  ['R','F',139,'34.51x51','A','I',383, ["hydrophobic"]],
  ['R','R',239,'-','A','L',346, ["hydrophobic"]],
  ['R','S',143,'34.55x55','A','A',39, ["hydrophobic", "van-der-waals"]],
  ['R','F',139,'34.51x51','A','C',379, ["hydrophobic"]],
  ['R','I',233,'5.72x72','A','L',394, ["hydrophobic"]],
  ['R','L',275,'6.37x37','A','L',393, ["hydrophobic", "van-der-waals"]],
  ['R','Q',229,'5.68x68','A','L',388, ["hydrophobic"]],
  ['R','I',135,'3.54x54','A','L',388, ["hydrophobic", "van-der-waals"]],
  ['R','T',274,'6.36x36','A','L',393, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['R','I',135,'3.54x54','A','Q',384, ["polar-backbone-sidechain", "van-der-waals"]],
  ['R','E',225,'5.64x64','A','R',380, ["polar-sidechain-sidechain"]],
  ['R','K',232,'5.71x71','A','D',381, ["hydrophobic", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
  ['R','I',233,'5.72x72','A','R',385, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['R','Q',229,'5.68x68','A','D',381, ["polar-sidechain-backbone", "polar-sidechain-sidechain"]],
  ['R','R',239,'-','A','T',350, ["hydrophobic", "polar-sidechain-sidechain"]],
  ['R','I',135,'3.54x54','A','Y',391, ["hydrophobic"]],
  ['R','K',235,'5.74x74','A','D',323, ["polar-backbone-sidechain"]],
  ['R','A',271,'6.33x33','A','L',393, ["hydrophobic", "van-der-waals"]],
  ['R','A',226,'5.65x65','A','L',388, ["hydrophobic", "van-der-waals"]],
  ['R','I',135,'3.54x54','A','H',387, ["hydrophobic", "polar-backbone-sidechain"]],
  ['R','F',139,'34.51x51','A','R',380, ["hydrophobic", "van-der-waals"]],
  ['R','T',274,'6.36x36','A','Y',391, ["polar-sidechain-backbone"]],
  ['R','Q',142,'34.54x54','A','I',383, ["hydrophobic"]],
  ['R','R',228,'5.67x67','A','D',381, ["polar-sidechain-sidechain"]],
  ['R','P',138,'34.50x50','A','R',380, ["hydrophobic"]],
  ['R','T',136,'3.55x55','A','R',380, ["polar-backbone-sidechain", "van-der-waals"]],
  ['R','F',139,'34.51x51','A','F',376, ["edge-to-face", "face-to-edge", "hydrophobic", "van-der-waals"]],
  ['R','L',230,'5.69x69','A','L',394, ["hydrophobic", "van-der-waals"]],
  ['R','P',138,'34.50x50','A','I',383, ["hydrophobic", "van-der-waals"]],
  ['R','R',228,'5.67x67','A','Q',384, ["polar-sidechain-sidechain"]],
  ['R','I',233,'5.72x72','A','Y',358, ["hydrophobic", "van-der-waals"]],
  ['R','I',135,'3.54x54','A','L',393, ["hydrophobic", "van-der-waals"]],
  ['R','Q',142,'34.54x54','A','H',387, ["polar-sidechain-sidechain"]],
  ['R','P',138,'34.50x50','A','H',387, ["hydrophobic"]],
  ['R','Y',141,'34.53x53','A','H',387, ["edge-to-face", "face-to-edge", "pi-cation", "hydrophobic"]],
  ['R','D',130,'3.49x49','A','Y',391, ["polar-sidechain-sidechain"]],
  ['R','Q',229,'5.68x68','A','Q',384, ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','R',131,'3.50x50','A','Y',391, ["cation-pi", "hydrophobic", "van-der-waals"]],
  ['R','A',134,'3.53x53','A','H',387, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['R','T',274,'6.36x36','A','E',392, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['R','R',239,'-','A','R',347, ["hydrophobic", "polar-sidechain-backbone"]],
  ['R','K',140,'34.52x52','A','R',380, ["polar-sidechain-sidechain"]],
  ['R','V',222,'5.61x61','A','L',393, ["hydrophobic", "van-der-waals"]],
  ['R','F',139,'34.51x51','A','H',41, ["edge-to-face", "face-to-edge", "pi-cation", "hydrophobic"]],
  ['R','P',138,'34.50x50','A','Q',384, ["hydrophobic", "van-der-waals"]],
  ],
  '4x1h' : [
  ['A','A',246,'6.29x29','C','F',350, ["hydrophobic", "van-der-waals"]],
  ['A','V',139,'3.54x54','C','V',340, ["hydrophobic", "water-mediated"]],
  ['A','V',138,'3.53x53','C','D',343, ["hydrophobic", "polar-backbone-sidechain", "water-mediated"]],
  ['A','T',242,'6.25x25','C','F',350, ["hydrophobic", "van-der-waals"]],
  ['A','M',309,'7.56x56','C','G',348, ["water-mediated", "water-mediated"]],
  ['A','I',305,'7.52x52','C','G',348, ["water-mediated"]],
  ['A','Q',312,'8.49x49','C','S',346, ["water-mediated"]],
  ['A','K',141,'3.56x56','C','V',340, ["water-mediated"]],
  ['A','N',73,'2.40x40','C','G',348, ["water-mediated"]],
  ['A','T',229,'5.64x64','C','V',340, ["hydrophobic"]],
  ['A','E',249,'6.32x32','C','F',350, ["hydrophobic"]],
  ['A','V',250,'6.33x33','C','L',349, ["hydrophobic"]],
  ['A','A',233,'5.68x68','C','V',340, ["hydrophobic"]],
  ['A','N',73,'2.40x40','C','C',347, ["water-mediated", "water-mediated"]],
  ['A','V',139,'3.54x54','C','L',344, ["hydrophobic"]],
  ['A','L',72,'2.39x39','C','C',347, ["water-mediated", "water-mediated"]],
  ['A','T',70,'2.37x37','C','S',346, ["water-mediated"]],
  ['A','A',246,'6.29x29','C','L',344, ["hydrophobic"]],
  ['A','N',73,'2.40x40','C','S',346, ["water-mediated", "water-mediated"]],
  ['A','L',72,'2.39x39','C','S',346, ["hydrophobic", "van-der-waals", "water-mediated", "water-mediated"]],
  ['A','K',245,'6.28x28','C','F',350, ["hydrophobic"]],
  ['A','V',230,'5.65x65','C','L',344, ["hydrophobic"]],
  ['A','N',310,'8.47x47','C','S',346, ["water-mediated"]],
  ['A','A',246,'6.29x29','C','L',341, ["hydrophobic", "van-der-waals"]],
  ['A','N',310,'8.47x47','C','C',347, ["polar-sidechain-backbone", "water-mediated"]],
  ['A','K',311,'8.48x48','C','F',350, ["polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
  ['A','T',242,'6.25x25','C','L',341, ["hydrophobic"]],
  ['A','M',253,'6.36x36','C','G',348, ["water-mediated"]],
  ['A','R',135,'3.50x50','C','C',347, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals", "water-mediated", "water-mediated"]],
  ['A','R',135,'3.50x50','C','S',346, ["water-mediated"]],
  ['A','R',135,'3.50x50','C','L',349, ["hydrophobic", "van-der-waals"]],
  ['A','K',141,'3.56x56','C','D',343, ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals", "water-mediated"]],
  ['A','N',310,'8.47x47','C','G',348, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals", "water-mediated", "water-mediated", "water-mediated", "water-mediated"]],
  ['A','L',226,'5.61x61','C','L',344, ["hydrophobic"]],
  ['A','V',139,'3.54x54','C','D',343, ["hydrophobic", "water-mediated"]],
  ['A','V',250,'6.33x33','C','L',344, ["hydrophobic"]],
  ['A','L',226,'5.61x61','C','L',349, ["hydrophobic"]],
  ['A','V',138,'3.53x53','C','V',340, ["water-mediated"]],
  ['A','T',243,'6.26x26','C','L',341, ["hydrophobic"]],
  ['A','P',71,'2.38x38','C','S',346, ["water-mediated"]],
  ['A','Y',306,'7.53x53','C','G',348, ["water-mediated", "water-mediated"]],
  ['A','Q',312,'8.49x49','C','C',347, ["water-mediated"]],
  ['A','Q',312,'8.49x49','C','G',348, ["water-mediated"]],
  ['A','R',135,'3.50x50','C','G',348, ["water-mediated", "water-mediated"]],
  ],
  '5g53' : [
  ['A','K',227,'6.29x29','C','L',394, ["hydrophobic", "polar-sidechain-sidechain"]],
  ['A','L',235,'6.37x37','C','L',393, ["hydrophobic", "van-der-waals"]],
  ['A','Q',207,'5.68x68','C','L',388, ["hydrophobic"]],
  ['A','M',211,'5.72x72','C','Y',358, ["hydrophobic", "polar-sidechain-sidechain"]],
  ['A','Q',207,'5.68x68','C','Q',384, ["polar-sidechain-sidechain", "h-bond"]],
  ['A','A',203,'5.64x64','C','Q',384, ["hydrophobic", "polar-backbone-sidechain"]],
  ['A','L',110,'34.51x51','C','V',217, ["hydrophobic", "van-der-waals"]],
  ['A','I',106,'3.54x54','C','Q',384, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['A','I',200,'5.61x61','C','L',388, ["hydrophobic"]],
  ['A','P',109,'34.50x50','C','R',380, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['A','R',293,'8.48x48','C','E',392, ["hydrophobic", "polar-backbone-sidechain"]],
  ['A','Q',207,'5.68x68','C','Y',360, ["polar-sidechain-sidechain"]],
  ['A','I',200,'5.61x61','C','L',393, ["hydrophobic"]],
  ['A','L',208,'5.69x69','C','L',394, ["hydrophobic", "van-der-waals"]],
  ['A','A',204,'5.65x65','C','L',388, ["hydrophobic"]],
  ['A','L',110,'34.51x51','C','R',380, ["hydrophobic"]],
  ['A','P',109,'34.50x50','C','Q',384, ["hydrophobic"]],
  ['A','Q',207,'5.68x68','C','R',385, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['A','Q',207,'5.68x68','C','D',381, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['A','R',111,'34.52x52','C','R',380, ["polar-sidechain-sidechain"]],
  ['A','I',106,'3.54x54','C','L',388, ["hydrophobic", "van-der-waals"]],
  ['A','Y',112,'34.53x53','C','H',387, ["edge-to-face", "face-to-edge", "pi-cation", "hydrophobic"]],
  ['A','L',110,'34.51x51','C','I',383, ["hydrophobic"]],
  ['A','I',106,'3.54x54','C','Y',391, ["hydrophobic"]],
  ['A','R',107,'3.55x55','C','R',380, ["polar-backbone-sidechain", "van-der-waals"]],
  ['A','Q',210,'5.71x71','C','D',381, ["polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
  ['A','L',110,'34.51x51','C','C',379, ["hydrophobic"]],
  ['A','L',110,'34.51x51','C','H',41, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['A','A',203,'5.64x64','C','L',388, ["hydrophobic"]],
  ['A','L',110,'34.51x51','C','F',219, ["hydrophobic"]],
  ['A','R',291,'7.56x56','C','E',392, ["polar-backbone-sidechain"]],
  ['A','R',102,'3.50x50','C','Y',391, ["cation-pi", "hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['A','I',108,'3.56x56','C','R',380, ["polar-backbone-sidechain"]],
  ['A','A',231,'6.33x33','C','L',393, ["hydrophobic"]],
  ['A','L',110,'34.51x51','C','F',376, ["hydrophobic", "van-der-waals"]],
  ['A','P',109,'34.50x50','C','I',383, ["hydrophobic", "van-der-waals"]],
  ['A','R',296,'8.51x51','C','E',392, ["polar-sidechain-sidechain", "van-der-waals"]],
  ['A','R',111,'34.52x52','C','V',217, ["hydrophobic", "van-der-waals"]],
  ['A','R',291,'7.56x56','C','Y',391, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['A','I',106,'3.54x54','C','H',387, ["hydrophobic"]],
  ['A','A',105,'3.53x53','C','H',387, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['A','R',111,'34.52x52','C','D',215, ["polar-sidechain-backbone", "van-der-waals"]],
  ['B','R',111,'34.52x52','D','D',215, ["polar-sidechain-backbone", "van-der-waals"]],
  ['B','I',106,'3.54x54','D','L',388, ["hydrophobic", "van-der-waals"]],
  ['B','I',106,'3.54x54','D','Y',391, ["hydrophobic"]],
  ['B','L',110,'34.51x51','D','I',383, ["hydrophobic"]],
  ['B','L',110,'34.51x51','D','H',41, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['B','L',110,'34.51x51','D','F',376, ["hydrophobic", "van-der-waals"]],
  ['B','P',109,'34.50x50','D','Q',384, ["hydrophobic"]],
  ['B','G',114,'34.55x55','D','A',39, ["hydrophobic", "van-der-waals"]],
  ['B','R',296,'8.51x51','D','E',392, ["polar-sidechain-sidechain", "van-der-waals"]],
  ['B','Y',112,'34.53x53','D','H',387, ["edge-to-face", "face-to-edge", "pi-cation", "hydrophobic"]],
  ['B','L',110,'34.51x51','D','R',380, ["hydrophobic"]],
  ['B','A',231,'6.33x33','D','L',393, ["hydrophobic", "van-der-waals"]],
  ['B','I',106,'3.54x54','D','Q',384, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['B','Q',207,'5.68x68','D','L',388, ["hydrophobic"]],
  ['B','A',204,'5.65x65','D','L',393, ["hydrophobic", "van-der-waals"]],
  ['B','A',203,'5.64x64','D','Q',384, ["hydrophobic", "polar-backbone-sidechain"]],
  ['B','R',102,'3.50x50','D','Y',391, ["hydrophobic", "polar-backbone-sidechain"]],
  ['B','A',105,'3.53x53','D','H',387, ["hydrophobic", "polar-backbone-sidechain"]],
  ['B','I',200,'5.61x61','D','L',393, ["hydrophobic", "van-der-waals"]],
  ['B','A',203,'5.64x64','D','L',388, ["hydrophobic"]],
  ['B','Q',207,'5.68x68','D','R',385, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['B','L',110,'34.51x51','D','V',217, ["hydrophobic", "van-der-waals"]],
  ['B','L',235,'6.37x37','D','L',393, ["hydrophobic", "van-der-waals"]],
  ['B','A',204,'5.65x65','D','L',388, ["hydrophobic"]],
  ['B','P',109,'34.50x50','D','I',383, ["hydrophobic", "van-der-waals"]],
  ['B','A',105,'3.53x53','D','Y',391, ["van-der-waals"]],
  ['B','I',106,'3.54x54','D','H',387, ["hydrophobic"]],
  ['B','R',111,'34.52x52','D','V',217, ["hydrophobic", "van-der-waals"]],
  ['B','I',108,'3.56x56','D','R',380, ["polar-backbone-sidechain"]],
  ['B','R',291,'7.56x56','D','E',392, ["polar-backbone-sidechain"]],
  ['B','L',110,'34.51x51','D','F',219, ["hydrophobic"]],
  ['B','R',111,'34.52x52','D','R',380, ["polar-sidechain-sidechain"]],
  ['B','P',109,'34.50x50','D','R',380, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['B','N',113,'34.54x54','D','A',39, ["hydrophobic", "van-der-waals"]],
  ['B','L',110,'34.51x51','D','C',379, ["hydrophobic"]],
  ['B','R',293,'8.48x48','D','E',392, ["polar-backbone-sidechain", "van-der-waals"]],
  ['B','R',107,'3.55x55','D','R',380, ["polar-backbone-sidechain", "van-der-waals"]],
  ['B','I',200,'5.61x61','D','L',388, ["hydrophobic"]],
  ['B','Q',207,'5.68x68','D','Q',384, ["polar-sidechain-sidechain"]],
  ['B','Q',207,'5.68x68','D','D',381, ["polar-sidechain-backbone"]],
  ],
  '5uz7' : [
  ['R','Y',243,'3.53x53','A','Y',391, ["van-der-waals"]],
  ['R','N',396,'8.48x48','A','E',392, ["polar-backbone-sidechain"]],
  ['R','T',330,'-','A','Y',358, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['R','R',180,'2.46x46','A','Q',390, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
  ['R','L',348,'6.45x45','A','L',393, ["hydrophobic"]],
  ['R','Q',415,'8.67x67','B','V',307, ["polar-sidechain-backbone", "van-der-waals"]],
  ['R','F',253,'-','A','V',217, ["hydrophobic"]],
  ['R','T',254,'-','A','H',387, ["polar-sidechain-sidechain"]],
  ['R','N',396,'8.48x48','A','R',356, ["polar-sidechain-sidechain"]],
  ['R','F',253,'-','A','H',41, ["polar-backbone-sidechain", "van-der-waals"]],
  ['R','V',249,'3.59x59','A','Q',384, ["polar-backbone-sidechain"]],
  ['R','R',180,'2.46x46','A','Y',391, ["cation-pi", "hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['R','I',248,'3.58x58','A','Q',384, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['R','H',184,'2.50x50','A','Y',391, ["hydrophobic"]],
  ['R','Q',415,'8.67x67','B','Q',44, ["hydrophobic", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
  ['R','V',252,'-','A','R',380, ["hydrophobic", "van-der-waals"]],
  ['R','K',326,'5.64x64','A','R',385, ["polar-backbone-sidechain"]],
  ['R','Q',408,'8.60x60','B','H',311, ["polar-sidechain-backbone", "van-der-waals"]],
  ['R','I',411,'8.63x63','B','V',307, ["hydrophobic"]],
  ['R','I',248,'3.58x58','A','H',387, ["hydrophobic"]],
  ['R','V',252,'-','A','I',383, ["hydrophobic", "van-der-waals"]],
  ['R','H',331,'-','A','Y',358, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['R','K',326,'5.64x64','A','R',380, ["polar-sidechain-sidechain"]],
  ['R','C',394,'7.60x60','A','E',392, ["polar-backbone-sidechain"]],
  ['R','K',326,'5.64x64','A','Q',384, ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','L',244,'3.54x54','A','Y',391, ["hydrophobic"]],
  ['R','L',247,'3.57x57','A','H',387, ["hydrophobic", "van-der-waals"]],
  ['R','R',404,'8.56x56','B','D',312, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','V',252,'-','A','Q',384, ["hydrophobic", "van-der-waals"]],
  ['R','T',345,'6.42x42','A','L',393, ["hydrophobic", "van-der-waals"]],
  ['R','E',329,'-','A','R',385, ["polar-backbone-sidechain"]],
  ['R','I',248,'3.58x58','A','L',388, ["hydrophobic"]],
  ['R','L',323,'5.61x61','A','L',388, ["hydrophobic", "van-der-waals"]],
  ['R','M',327,'5.65x65','A','L',394, ["hydrophobic", "van-der-waals"]],
  ['R','Q',408,'8.60x60','B','A',309, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['R','L',247,'3.57x57','A','Y',391, ["hydrophobic"]],
  ['R','L',348,'6.45x45','A','E',392, ["hydrophobic"]],
  ['R','Q',408,'8.60x60','B','G',310, ["polar-sidechain-backbone"]],
  ],
  '5vai' : [
  ['R','Y',402,'7.57x57','A','E',392, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['R','K',415,'8.56x56','B','D',291, ["polar-sidechain-backbone"]],
  ['R','N',338,'-','A','Y',360, ["hydrophobic", "van-der-waals"]],
  ['R','L',255,'3.58x58','A','R',380, ["polar-backbone-sidechain"]],
  ['R','R',419,'8.60x60','B','N',293, ["polar-sidechain-sidechain"]],
  ['R','R',419,'8.60x60','B','H',311, ["polar-sidechain-backbone"]],
  ['R','L',251,'3.54x54','A','Y',391, ["hydrophobic", "van-der-waals"]],
  ['R','H',180,'2.50x50','A','Y',391, ["cation-pi", "hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','R',170,'12.48x48','B','R',52, ["polar-sidechain-sidechain"]],
  ['R','L',255,'3.58x58','A','Q',384, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['R','H',171,'12.49x49','B','D',312, ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','K',334,'5.64x64','A','R',385, ["polar-backbone-sidechain", "van-der-waals"]],
  ['R','R',419,'8.60x60','B','G',310, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['R','L',339,'-','A','L',394, ["hydrophobic"]],
  ['R','R',264,'4.40x41','A','Q',35, ["polar-backbone-sidechain"]],
  ['R','N',406,'8.47x47','A','E',392, ["hydrophobic", "polar-backbone-sidechain", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','A',256,'3.59x59','A','Q',384, ["polar-backbone-sidechain"]],
  ['R','S',352,'6.41x41','A','E',392, ["polar-sidechain-backbone"]],
  ['R','R',419,'8.60x60','B','A',309, ["hydrophobic", "van-der-waals"]],
  ['R','L',359,'6.48x48','A','Y',391, ["hydrophobic", "van-der-waals"]],
  ['R','E',408,'8.49x49','A','Q',390, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
  ['R','K',342,'-','A','T',350, ["polar-sidechain-backbone", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
  ['R','L',254,'3.57x57','A','H',387, ["hydrophobic", "van-der-waals"]],
  ['R','N',407,'8.48x48','A','E',392, ["polar-backbone-sidechain"]],
  ['R','E',262,'4.38x39','A','Q',35, ["polar-backbone-sidechain", "van-der-waals"]],
  ['R','L',356,'6.45x45','A','L',393, ["hydrophobic"]],
  ['R','R',419,'8.60x60','B','F',292, ["polar-sidechain-backbone"]],
  ['R','V',405,'7.60x60','A','E',392, ["polar-backbone-sidechain", "van-der-waals"]],
  ['R','S',352,'6.41x41','A','L',394, ["polar-sidechain-backbone"]],
  ['R','N',338,'-','A','C',359, ["polar-sidechain-backbone"]],
  ['R','R',176,'2.46x46','A','Q',390, ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','R',176,'2.46x46','A','Y',391, ["hydrophobic"]],
  ['R','K',415,'8.56x56','B','F',292, ["cation-pi", "hydrophobic", "van-der-waals"]],
  ['R','L',356,'6.45x45','A','Y',391, ["hydrophobic"]],
  ['R','L',339,'-','A','Y',358, ["hydrophobic", "van-der-waals"]],
  ['R','K',334,'5.64x64','A','D',381, ["hydrophobic", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
  ['R','R',419,'8.60x60','B','D',312, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','E',412,'8.53x53','B','D',312, ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','N',406,'8.47x47','A','Q',390, ["polar-sidechain-backbone"]],
  ['R','T',353,'6.42x42','A','L',393, ["hydrophobic"]],
  ['R','A',256,'3.59x59','A','R',380, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['R','S',261,'4.37x38','A','Q',35, ["hydrophobic", "polar-backbone-sidechain", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','Q',263,'4.39x40','A','Q',35, ["polar-backbone-sidechain", "van-der-waals"]],
  ['R','V',331,'5.61x61','A','L',388, ["hydrophobic", "van-der-waals"]],
  ['R','Q',263,'4.39x40','A','Q',31, ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','L',339,'-','A','R',385, ["polar-backbone-sidechain"]],
  ['R','S',352,'6.41x41','A','L',393, ["polar-sidechain-backbone", "van-der-waals"]],
  ['R','L',401,'7.56x56','A','E',392, ["polar-backbone-sidechain", "van-der-waals"]],
  ],
  '6b3j' : [
  ['R','N',407,'8.48x48','A','E',392, ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','K',334,'5.64x64','A','L',388, ["hydrophobic"]],
  ['R','K',415,'8.56x56','B','D',312, ["polar-sidechain-sidechain", "h-bond"]],
  ['R','R',176,'2.46x46','A','Y',391, ["cation-pi", "hydrophobic"]],
  ['R','L',255,'3.58x58','A','Q',384, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['R','V',331,'5.61x61','A','L',393, ["hydrophobic"]],
  ['R','E',423,'8.64x64','B','Q',44, ["hydrophobic"]],
  ['R','S',261,'-','A','Q',35, ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','L',251,'3.54x54','A','Y',391, ["hydrophobic", "van-der-waals"]],
  ['R','V',331,'5.61x61','A','L',394, ["hydrophobic"]],
  ['R','L',255,'3.58x58','A','L',388, ["hydrophobic"]],
  ['R','S',258,'-','A','I',383, ["hydrophobic"]],
  ['R','R',176,'2.46x46','A','Q',390, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain"]],
  ['R','V',259,'-','A','V',217, ["hydrophobic"]],
  ['R','K',334,'5.64x64','A','R',385, ["hydrophobic", "polar-backbone-sidechain"]],
  ['R','R',348,'6.37x37','A','L',394, ["hydrophobic"]],
  ['R','Y',250,'3.53x53','A','Y',391, ["van-der-waals"]],
  ['R','S',352,'6.41x41','A','L',394, ["polar-sidechain-backbone"]],
  ['R','K',334,'5.64x64','A','Q',384, ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','V',331,'5.61x61','A','L',388, ["hydrophobic"]],
  ['R','H',171,'12.49x49','B','D',312, ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','R',419,'8.60x60','B','A',309, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['R','L',255,'3.58x58','A','H',387, ["hydrophobic"]],
  ['R','E',262,'4.38x39','A','R',38, ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','V',327,'5.57x57','A','L',393, ["hydrophobic"]],
  ['R','L',254,'3.57x57','A','H',387, ["hydrophobic", "van-der-waals"]],
  ['R','H',180,'2.50x50','A','Y',391, ["hydrophobic"]],
  ['R','K',334,'5.64x64','A','D',381, ["hydrophobic", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
  ['R','L',356,'6.45x45','A','L',393, ["hydrophobic"]],
  ['R','S',352,'6.41x41','A','L',393, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['R','E',262,'4.38x39','A','K',34, ["polar-sidechain-sidechain"]],
  ],
  '6cmo' : [
  ['R','M',309,'7.56x56','A','G',352, ["van-der-waals"]],
  ['R','V',139,'3.54x54','A','L',348, ["hydrophobic"]],
  ['R','K',141,'3.56x56','A','D',193, ["polar-sidechain-sidechain"]],
  ['R','S',240,'-','A','E',318, ["hydrophobic", "polar-backbone-sidechain", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
  ['R','K',66,'12.48x48','B','D',312, ["polar-sidechain-sidechain"]],
  ['R','V',139,'3.54x54','A','N',347, ["hydrophobic"]],
  ['R','K',311,'8.48x48','A','K',349, ["polar-sidechain-backbone"]],
  ['R','T',242,'6.25x25','A','D',315, ["polar-sidechain-backbone"]],
  ['R','E',239,'-','A','E',318, ["polar-backbone-sidechain"]],
  ['R','T',242,'6.25x25','A','F',354, ["hydrophobic", "van-der-waals"]],
  ['R','Q',237,'5.72x72','A','D',341, ["hydrophobic"]],
  ['R','R',135,'3.50x50','A','L',353, ["hydrophobic"]],
  ['R','R',147,'34.55x55','A','R',32, ["hydrophobic", "polar-sidechain-backbone"]],
  ['R','A',246,'6.29x29','A','F',354, ["hydrophobic", "van-der-waals"]],
  ['R','E',249,'6.32x32','A','L',353, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['R','A',246,'6.29x29','A','L',348, ["hydrophobic"]],
  ['R','K',311,'8.48x48','A','L',353, ["polar-sidechain-backbone"]],
  ['R','S',240,'-','A','K',345, ["polar-sidechain-sidechain", "h-bond"]],
  ['R','Q',237,'5.72x72','A','Y',320, ["polar-sidechain-sidechain"]],
  ['R','E',239,'-','A','Y',320, ["hydrophobic", "van-der-waals"]],
  ['R','E',249,'6.32x32','A','F',354, ["hydrophobic"]],
  ['R','K',311,'8.48x48','A','F',354, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','A',241,'6.24x24','A','E',318, ["polar-backbone-sidechain"]],
  ['R','K',245,'6.28x28','A','F',354, ["hydrophobic"]],
  ['R','N',310,'8.47x47','A','G',352, ["hydrophobic", "polar-sidechain-backbone"]],
  ['R','M',253,'6.36x36','A','L',353, ["hydrophobic", "van-der-waals"]],
  ['R','K',311,'8.48x48','A','G',352, ["polar-sidechain-backbone"]],
  ['R','R',147,'34.55x55','A','A',31, ["polar-sidechain-backbone"]],
  ['R','V',250,'6.33x33','A','L',353, ["hydrophobic"]],
  ['R','T',243,'6.26x26','A','D',341, ["polar-sidechain-sidechain"]],
  ],
  '6d9h' : [
  ['R','Q',210,'5.68x68','A','I',345, ["hydrophobic", "van-der-waals"]],
  ['R','K',224,'6.25x25','A','E',319, ["polar-sidechain-sidechain"]],
  ['R','R',291,'7.56x56','A','G',353, ["hydrophobic"]],
  ['R','I',207,'5.65x65','A','L',349, ["hydrophobic"]],
  ['R','Q',38,'12.48x48','B','F',335, ["hydrophobic", "van-der-waals"]],
  ['R','K',231,'6.32x32','A','F',355, ["polar-sidechain-backbone"]],
  ['R','L',113,'34.51x51','A','L',195, ["hydrophobic"]],
  ['R','R',108,'3.53x53','A','N',348, ["hydrophobic", "polar-backbone-sidechain", "polar-sidechain-sidechain"]],
  ['R','Q',210,'5.68x68','A','D',342, ["polar-sidechain-sidechain", "van-der-waals"]],
  ['R','P',112,'34.50x50','A','N',348, ["polar-backbone-sidechain"]],
  ['R','L',211,'5.69x69','A','K',346, ["hydrophobic"]],
  ['R','K',213,'5.71x71','A','D',342, ["polar-sidechain-sidechain"]],
  ['R','K',228,'6.29x29','A','F',355, ["cation-pi", "hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','L',236,'6.37x37','A','L',354, ["hydrophobic", "van-der-waals"]],
  ['R','V',203,'5.61x61','A','L',354, ["hydrophobic"]],
  ['R','K',224,'6.25x25','A','D',316, ["hydrophobic", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
  ['R','I',207,'5.65x65','A','L',354, ["hydrophobic"]],
  ['R','V',109,'3.54x54','A','L',349, ["hydrophobic", "van-der-waals"]],
  ['R','R',108,'3.53x53','A','D',351, ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','L',113,'34.51x51','A','T',341, ["hydrophobic"]],
  ['R','K',294,'8.49x49','A','D',351, ["hydrophobic", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
  ['R','P',112,'34.50x50','A','I',344, ["hydrophobic", "van-der-waals"]],
  ['R','P',112,'34.50x50','A','I',345, ["hydrophobic", "van-der-waals"]],
  ['R','L',113,'34.51x51','A','I',344, ["hydrophobic"]],
  ['R','R',105,'3.50x50','A','L',354, ["hydrophobic"]],
  ['R','I',232,'6.33x33','A','F',355, ["hydrophobic"]],
  ['R','I',292,'8.47x47','A','G',353, ["hydrophobic"]],
  ['R','R',108,'3.53x53','A','C',352, ["polar-sidechain-sidechain", "van-der-waals"]],
  ['R','K',228,'6.29x29','A','D',316, ["polar-sidechain-sidechain"]],
  ['R','N',37,'1.60x60','B','D',312, ["polar-sidechain-sidechain"]],
  ['R','F',45,'2.40x40','A','D',351, ["van-der-waals"]],
  ['R','L',113,'34.51x51','A','F',337, ["hydrophobic"]],
  ['R','D',42,'2.37x37','A','D',351, ["polar-sidechain-sidechain"]],
  ['R','R',105,'3.50x50','A','C',352, ["hydrophobic", "polar-sidechain-backbone"]],
  ['R','K',301,'8.56x56','B','D',312, ["polar-sidechain-backbone"]],
  ['R','I',292,'8.47x47','A','C',352, ["hydrophobic", "van-der-waals"]],
  ['R','I',232,'6.33x33','A','L',354, ["hydrophobic", "van-der-waals"]],
  ['R','K',301,'8.56x56','B','F',292, ["cation-pi", "hydrophobic", "van-der-waals"]],
  ],
  '6dde' : [
  ['R','R',263,'-','A','I',319, ["polar-sidechain-backbone"]],
  ['R','R',263,'-','A','Y',320, ["cation-pi", "hydrophobic", "van-der-waals"]],
  ['R','V',173,'34.51x51','A','F',336, ["hydrophobic", "van-der-waals"]],
  ['R','D',177,'34.55x55','A','R',32, ["polar-backbone-sidechain", "polar-sidechain-sidechain"]],
  ['R','L',259,'5.65x65','A','I',344, ["hydrophobic"]],
  ['R','V',169,'3.54x54','A','L',348, ["hydrophobic", "van-der-waals"]],
  ['R','L',176,'34.54x54','A','I',343, ["hydrophobic", "van-der-waals"]],
  ['R','L',176,'34.54x54','A','L',194, ["hydrophobic"]],
  ['R','V',173,'34.51x51','A','L',194, ["hydrophobic"]],
  ['R','M',264,'-','A','D',341, ["hydrophobic", "polar-backbone-sidechain"]],
  ['R','P',172,'34.50x50','A','I',343, ["hydrophobic", "van-der-waals"]],
  ['R','L',176,'34.54x54','A','R',32, ["hydrophobic", "polar-backbone-sidechain"]],
  ['R','R',277,'6.32x32','A','L',353, ["polar-sidechain-backbone"]],
  ['R','T',103,'2.39x39','A','C',351, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','V',173,'34.51x51','A','D',193, ["hydrophobic", "van-der-waals"]],
  ['R','L',259,'5.65x65','A','L',348, ["hydrophobic"]],
  ['R','D',340,'8.47x47','A','G',352, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['R','M',281,'6.36x36','A','L',353, ["hydrophobic", "van-der-waals"]],
  ['R','I',278,'6.33x33','A','F',354, ["hydrophobic", "van-der-waals"]],
  ['R','E',270,'6.25x25','A','D',315, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','M',255,'5.61x61','A','L',353, ["hydrophobic"]],
  ['R','T',103,'2.39x39','A','D',350, ["polar-sidechain-backbone"]],
  ['R','V',262,'5.68x68','A','D',341, ["hydrophobic", "van-der-waals"]],
  ['R','K',271,'6.26x26','A','K',314, ["polar-sidechain-backbone", "van-der-waals"]],
  ['R','R',165,'3.50x50','A','L',353, ["hydrophobic"]],
  ['R','R',165,'3.50x50','A','C',351, ["polar-sidechain-sidechain", "van-der-waals"]],
  ['R','R',179,'34.57x57','A','N',347, ["polar-sidechain-sidechain", "van-der-waals"]],
  ['R','V',262,'5.68x68','A','I',344, ["hydrophobic"]],
  ['R','P',172,'34.50x50','A','T',340, ["hydrophobic"]],
  ['R','A',168,'3.53x53','A','N',347, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['R','I',278,'6.33x33','A','L',348, ["hydrophobic"]],
  ['R','M',264,'-','A','T',316, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['R','K',271,'6.26x26','A','D',315, ["hydrophobic", "van-der-waals"]],
  ['R','I',278,'6.33x33','A','L',353, ["hydrophobic"]],
  ['R','P',172,'34.50x50','A','I',344, ["hydrophobic", "van-der-waals"]],
  ['R','M',264,'-','A','K',345, ["hydrophobic"]],
  ['R','S',268,'6.23x23','A','D',315, ["polar-backbone-sidechain"]],
  ['R','R',258,'5.64x64','A','I',344, ["hydrophobic", "van-der-waals"]],
  ['R','R',263,'-','A','D',341, ["polar-backbone-sidechain"]],
  ['R','R',182,'4.40x40','A','R',24, ["hydrophobic", "polar-sidechain-sidechain"]],
  ],
  '6ddf' : [
  ['R','E',341,'8.48x48','A','F',354, ["hydrophobic"]],
  ['R','R',263,'-','A','I',319, ["polar-sidechain-backbone", "van-der-waals"]],
  ['R','V',173,'34.51x51','A','F',336, ["hydrophobic", "van-der-waals"]],
  ['R','L',176,'34.54x54','A','I',343, ["hydrophobic"]],
  ['R','L',176,'34.54x54','A','L',194, ["hydrophobic"]],
  ['R','M',264,'-','A','D',341, ["polar-backbone-sidechain"]],
  ['R','R',277,'6.32x32','A','F',354, ["polar-sidechain-backbone", "van-der-waals"]],
  ['R','P',172,'34.50x50','A','I',343, ["hydrophobic", "van-der-waals"]],
  ['R','R',277,'6.32x32','A','L',353, ["polar-sidechain-backbone", "van-der-waals"]],
  ['R','V',173,'34.51x51','A','L',194, ["hydrophobic"]],
  ['R','V',173,'34.51x51','A','D',193, ["hydrophobic"]],
  ['R','L',259,'5.65x65','A','I',344, ["hydrophobic"]],
  ['R','D',340,'8.47x47','A','G',352, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['R','K',174,'34.52x52','A','D',193, ["polar-sidechain-sidechain"]],
  ['R','E',341,'8.48x48','A','L',353, ["polar-sidechain-backbone"]],
  ['R','R',165,'3.50x50','A','C',351, ["hydrophobic", "polar-sidechain-backbone"]],
  ['R','V',169,'3.54x54','A','L',348, ["hydrophobic", "van-der-waals"]],
  ['R','P',172,'34.50x50','A','N',347, ["polar-backbone-sidechain"]],
  ['R','P',172,'34.50x50','A','T',340, ["hydrophobic", "van-der-waals"]],
  ['R','A',168,'3.53x53','A','N',347, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['R','I',278,'6.33x33','A','L',348, ["hydrophobic"]],
  ['R','K',271,'6.26x26','A','D',315, ["hydrophobic", "van-der-waals"]],
  ['R','M',255,'5.61x61','A','L',353, ["hydrophobic"]],
  ['R','K',271,'6.26x26','A','K',317, ["polar-sidechain-backbone"]],
  ['R','R',263,'-','A','D',341, ["hydrophobic", "polar-backbone-sidechain"]],
  ['R','I',278,'6.33x33','A','L',353, ["hydrophobic"]],
  ['R','V',262,'5.68x68','A','I',344, ["hydrophobic"]],
  ['R','R',263,'-','A','Y',320, ["cation-pi", "hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','D',177,'34.55x55','A','R',32, ["polar-sidechain-sidechain"]],
  ['R','E',341,'8.48x48','A','G',352, ["polar-sidechain-backbone", "van-der-waals"]],
  ['R','P',172,'34.50x50','A','I',344, ["hydrophobic", "van-der-waals"]],
  ['R','A',168,'3.53x53','A','C',351, ["van-der-waals"]],
  ['R','D',164,'3.49x49','A','C',351, ["polar-sidechain-sidechain"]],
  ['R','L',176,'34.54x54','A','R',32, ["hydrophobic"]],
  ['R','T',103,'2.39x39','A','C',351, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain"]],
  ['R','L',259,'5.65x65','A','L',348, ["hydrophobic", "van-der-waals"]],
  ['R','R',263,'-','A','E',318, ["van-der-waals"]],
  ['R','M',281,'6.36x36','A','L',353, ["hydrophobic", "van-der-waals"]],
  ['R','I',278,'6.33x33','A','F',354, ["hydrophobic", "van-der-waals"]],
  ['R','E',270,'6.25x25','A','D',315, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','T',103,'2.39x39','A','D',350, ["polar-sidechain-backbone", "van-der-waals"]],
  ['R','V',262,'5.68x68','A','D',341, ["hydrophobic", "van-der-waals"]],
  ['R','K',271,'6.26x26','A','K',314, ["polar-sidechain-backbone", "van-der-waals"]],
  ['R','R',165,'3.50x50','A','L',353, ["hydrophobic"]],
  ['R','R',179,'34.57x57','A','N',347, ["polar-sidechain-sidechain"]],
  ['R','D',340,'8.47x47','A','L',353, ["polar-sidechain-backbone"]],
  ['R','R',179,'34.57x57','A','C',351, ["polar-sidechain-sidechain", "van-der-waals"]],
  ['R','M',264,'-','A','T',316, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['R','M',264,'-','A','K',345, ["hydrophobic"]],
  ['R','S',268,'6.23x23','A','D',315, ["polar-backbone-sidechain", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','R',258,'5.64x64','A','I',344, ["hydrophobic", "van-der-waals"]],
  ],
  '6g79' : [
  ['S','L',316,'6.37x37','A','L',353, ["hydrophobic"]],
  ['S','A',235,'5.65x65','A','L',348, ["hydrophobic"]],
  ['S','R',238,'5.68x68','A','D',341, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
  ['S','R',308,'6.29x29','A','N',316, ["polar-sidechain-backbone"]],
  ['S','T',315,'6.36x36','A','L',353, ["hydrophobic", "polar-sidechain-backbone"]],
  ['S','V',155,'34.51x51','A','T',340, ["hydrophobic"]],
  ['S','R',147,'3.50x50','A','C',351, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
  ['S','A',150,'3.53x53','A','N',347, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['S','S',372,'7.56x56','A','C',351, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['S','I',151,'3.54x54','A','N',347, ["hydrophobic"]],
  ['S','K',311,'6.32x32','A','G',352, ["polar-sidechain-backbone"]],
  ['S','A',154,'34.50x50','A','I',344, ["hydrophobic", "van-der-waals"]],
  ['S','I',151,'3.54x54','A','I',344, ["hydrophobic"]],
  ['S','A',312,'6.33x33','A','L',353, ["hydrophobic", "van-der-waals"]],
  ['S','R',238,'5.68x68','A','L',348, ["hydrophobic"]],
  ['S','T',315,'6.36x36','A','G',352, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['S','I',239,'5.69x69','A','Y',354, ["hydrophobic"]],
  ['S','I',231,'5.61x61','A','L',353, ["hydrophobic", "van-der-waals"]],
  ['S','R',238,'5.68x68','A','I',344, ["hydrophobic", "van-der-waals"]],
  ['S','S',372,'7.56x56','A','G',352, ["hydrophobic"]],
  ['S','V',155,'34.51x51','A','F',336, ["hydrophobic"]],
  ['S','R',238,'5.68x68','A','Y',354, ["hydrophobic"]],
  ['S','A',154,'34.50x50','A','I',343, ["hydrophobic"]],
  ['S','R',308,'6.29x29','A','Y',354, ["cation-pi", "hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
  ['S','K',311,'6.32x32','A','Y',354, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
  ['S','A',235,'5.65x65','A','L',353, ["hydrophobic"]],
  ['S','R',161,'34.57x57','A','N',347, ["polar-sidechain-sidechain"]],
  ['S','N',373,'8.47x47','A','G',350, ["polar-sidechain-backbone"]],
  ['S','I',151,'3.54x54','A','L',348, ["hydrophobic", "van-der-waals"]],
  ],
  '6gdg' : [
  ['A','Q',207,'5.68x68','D','L',384, ["hydrophobic"]],
  ['A','L',110,'34.51x51','D','I',373, ["hydrophobic"]],
  ['A','Y',112,'34.53x53','D','H',377, ["edge-to-face", "face-to-edge", "pi-cation", "hydrophobic"]],
  ['A','Q',207,'5.68x68','D','Q',374, ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
  ['A','R',102,'3.50x50','D','Y',381, ["cation-pi", "hydrophobic", "van-der-waals"]],
  ['A','N',36,'12.49x49','B','D',312, ["hydrophobic", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
  ['A','R',111,'34.52x52','D','K',216, ["hydrophobic", "polar-sidechain-backbone"]],
  ['A','M',211,'5.72x72','D','Y',348, ["hydrophobic"]],
  ['A','L',110,'34.51x51','D','H',41, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['A','P',109,'34.50x50','D','R',370, ["hydrophobic"]],
  ['A','I',106,'3.54x54','D','H',377, ["hydrophobic", "van-der-waals"]],
  ['A','S',35,'12.48x48','B','D',333, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
  ['A','L',110,'34.51x51','D','V',217, ["hydrophobic"]],
  ['A','S',234,'6.36x36','D','L',383, ["hydrophobic", "van-der-waals"]],
  ['A','E',294,'8.49x49','D','Q',380, ["polar-backbone-sidechain"]],
  ['A','L',110,'34.51x51','D','C',369, ["hydrophobic"]],
  ['A','I',106,'3.54x54','D','Q',374, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['A','S',35,'12.48x48','B','D',312, ["hydrophobic"]],
  ['A','S',35,'12.48x48','B','F',335, ["hydrophobic", "van-der-waals"]],
  ['A','I',106,'3.54x54','D','Y',381, ["hydrophobic"]],
  ['A','I',200,'5.61x61','D','L',383, ["hydrophobic", "van-der-waals"]],
  ['A','N',34,'1.60x60','B','D',312, ["polar-sidechain-sidechain"]],
  ['A','L',208,'5.69x69','D','L',384, ["hydrophobic"]],
  ['A','Q',210,'5.71x71','D','D',371, ["polar-sidechain-sidechain"]],
  ['A','S',234,'6.36x36','D','E',382, ["polar-sidechain-backbone"]],
  ['A','L',110,'34.51x51','D','R',370, ["hydrophobic"]],
  ['A','R',111,'34.52x52','D','D',215, ["polar-sidechain-backbone"]],
  ['A','R',291,'7.56x56','D','Y',381, ["polar-sidechain-backbone"]],
  ['A','R',293,'8.48x48','D','Q',380, ["polar-backbone-sidechain"]],
  ['A','L',110,'34.51x51','D','F',366, ["hydrophobic", "van-der-waals"]],
  ['A','Q',207,'5.68x68','D','D',371, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
  ['A','R',111,'34.52x52','D','V',217, ["hydrophobic", "van-der-waals"]],
  ['A','Q',207,'5.68x68','D','L',378, ["hydrophobic"]],
  ['A','L',235,'6.37x37','D','L',383, ["hydrophobic", "van-der-waals"]],
  ['A','A',204,'5.65x65','D','L',378, ["hydrophobic"]],
  ['A','I',106,'3.54x54','D','L',378, ["hydrophobic", "van-der-waals"]],
  ['A','A',231,'6.33x33','D','L',383, ["hydrophobic", "van-der-waals"]],
  ['A','H',230,'6.32x32','D','E',382, ["polar-sidechain-backbone"]],
  ['A','N',113,'34.54x54','D','A',39, ["hydrophobic"]],
  ['A','P',109,'34.50x50','D','Q',374, ["hydrophobic", "van-der-waals"]],
  ['A','R',293,'8.48x48','D','E',382, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['A','A',105,'3.53x53','D','H',377, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['A','G',114,'34.55x55','D','A',39, ["hydrophobic"]],
  ['A','A',203,'5.64x64','D','L',378, ["hydrophobic"]],
  ['A','Q',207,'5.68x68','D','R',375, ["polar-sidechain-backbone"]],
  ['A','R',107,'3.55x55','D','R',370, ["polar-backbone-sidechain", "polar-sidechain-sidechain", "van-der-waals"]],
  ['A','R',291,'7.56x56','D','E',382, ["polar-backbone-sidechain", "polar-sidechain-backbone", "van-der-waals"]],
  ['A','N',113,'34.54x54','D','R',38, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['A','P',109,'34.50x50','D','I',373, ["hydrophobic", "van-der-waals"]],
  ['A','Q',38,'12.51x51','B','R',52, ["polar-sidechain-sidechain", "van-der-waals"]],
  ],
  '6e3y' : [
  ['R','V',243,'3.60x60','A','R',380, ["polar-backbone-sidechain"]],
  ['R','L',240,'3.57x57','A','Y',391, ["hydrophobic"]],
  ['R','V',242,'3.59x59','A','R',380, ["polar-backbone-sidechain", "van-der-waals"]],
  ['R','I',241,'3.58x58','A','H',387, ["hydrophobic"]],
  ['R','G',389,'8.48x48','A','E',392, ["polar-backbone-sidechain", "van-der-waals"]],
  ['R','E',248,'-','A','A',39, ["hydrophobic"]],
  ['R','K',319,'5.64x64','A','Q',384, ["hydrophobic", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
  ['R','T',323,'-','A','Y',358, ["hydrophobic"]],
  ['R','L',240,'3.57x57','A','H',387, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['R','F',246,'-','A','C',379, ["hydrophobic"]],
  ['R','L',316,'5.61x61','A','L',388, ["hydrophobic", "van-der-waals"]],
  ['R','K',333,'6.37x37','A','L',394, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','I',241,'3.58x58','A','L',388, ["hydrophobic"]],
  ['R','R',336,'6.40x40','A','L',393, ["polar-sidechain-backbone"]],
  ['R','I',241,'3.58x58','A','Y',391, ["hydrophobic"]],
  ['R','R',336,'6.40x40','A','L',394, ["polar-sidechain-sidechain"]],
  ['R','K',319,'5.64x64','A','D',381, ["hydrophobic", "van-der-waals"]],
  ['R','F',246,'-','A','R',380, ["pi-cation", "hydrophobic", "van-der-waals"]],
  ['R','N',388,'8.47x47','A','E',392, ["polar-sidechain-sidechain", "van-der-waals"]],
  ['R','K',319,'5.64x64','A','L',388, ["hydrophobic"]],
  ['R','S',168,'12.49x49','B','D',312, ["polar-sidechain-sidechain"]],
  ['R','R',173,'2.46x46','A','Y',391, ["cation-pi", "hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['R','R',173,'2.46x46','A','Q',390, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','E',248,'-','A','R',38, ["hydrophobic", "van-der-waals"]],
  ['R','F',246,'-','A','I',383, ["hydrophobic", "van-der-waals"]],
  ['R','R',397,'8.56x56','B','H',311, ["polar-sidechain-backbone"]],
  ['R','R',397,'8.56x56','B','D',312, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
  ['R','H',177,'2.50x50','A','Y',391, ["hydrophobic"]],
  ['R','S',168,'12.49x49','B','F',335, ["hydrophobic"]],
  ['R','L',320,'5.65x65','A','L',394, ["hydrophobic", "van-der-waals"]],
  ['R','V',242,'3.59x59','A','Q',384, ["polar-backbone-sidechain"]],
  ['R','V',245,'-','A','I',383, ["hydrophobic", "van-der-waals"]],
  ['R','F',387,'7.60x60','A','E',392, ["polar-backbone-sidechain"]],
  ['R','I',241,'3.58x58','A','L',393, ["hydrophobic"]],
  ['R','R',336,'6.40x40','A','E',392, ["polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','L',237,'3.54x54','A','Y',391, ["hydrophobic", "van-der-waals"]],
  ['R','F',246,'-','A','F',219, ["hydrophobic"]],
  ['R','I',241,'3.58x58','A','Q',384, ["polar-backbone-sidechain", "van-der-waals"]],
  ['R','K',319,'5.64x64','A','R',385, ["polar-backbone-sidechain", "van-der-waals"]],
  ['R','F',246,'-','A','H',41, ["pi-cation", "hydrophobic"]],
  ['R','L',341,'6.45x45','A','L',393, ["hydrophobic", "van-der-waals"]],
  ['R','I',312,'5.57x57','A','L',393, ["hydrophobic"]],
  ['R','V',245,'-','A','H',387, ["polar-backbone-sidechain"]],
  ['R','L',316,'5.61x61','A','L',393, ["hydrophobic"]],
  ['R','F',246,'-','A','V',217, ["hydrophobic", "van-der-waals"]],
  ['R','K',167,'12.48x48','B','R',52, ["hydrophobic", "polar-backbone-sidechain"]],
  ['R','V',245,'-','A','R',380, ["hydrophobic", "van-der-waals"]],
  ['R','F',246,'-','A','F',376, ["hydrophobic", "van-der-waals"]],
  ['R','R',397,'8.56x56','B','F',292, ["cation-pi"]],
  ['R','V',245,'-','A','Q',384, ["hydrophobic", "van-der-waals"]],
  ['R','L',316,'5.61x61','A','L',394, ["hydrophobic", "van-der-waals"]],
  ]
  }

let pdb_ids = Object.keys(dataset);
let data_t = Object.keys(dataset).map(key => dataset[key]);

for (let i = 0; i < pdb_ids.length; i++) {
  const pdb = pdb_ids[i];
  for (let j = 0; j < data_t[i].length; j++) {
    data_t[i][j].push(pdb);
}};

// https://stackoverflow.com/questions/10865025/merge-flatten-an-array-of-arrays-in-javascript/25804569#comment50580016_10865042
const flattenOnce = (array) => {return [].concat(...array)};
data_t = flattenOnce(data_t);

let keys = [
  "rec_chain",
  "rec_aa",
  "rec_pos",
  "rec_gn",
  "sig_chain",
  "sig_aa",
  // "sig_pos",
  "sig_gn",
  "int_ty",
  "pdb_id"
];

data_t = data_t.map(function(e) {
  let obj = {};
  keys.forEach(function(key, i) {
    obj[key] = e[i];
  });
  return obj;
});

// * DEFINE ADDITIONAL DATASETS
// These are used as the x and y axis tick mark labels
let data_t_rec = _.uniqBy(data_t, (t) => [t.rec_gn, t.pdb_id].join());
let data_t_sig = _.uniqBy(data_t, (t) => [t.sig_gn, t.pdb_id].join());

let int_ty = []
for (let i = 0; i < data_t.length; i++) {
  const int_arr = data_t[i].int_ty;
  int_ty.push(int_arr)
}
int_ty = flattenOnce(int_ty).filter((v, i, a) => a.indexOf(v) === i);
let rm_index = int_ty.indexOf("undefined");
if (rm_index > -1) {
  int_ty.splice(rm_index, 1);
}

// * SETTING UP SVG FOR OUTPUT
let svg = d3
  .select("body")
  .select("div#content")
  .append("div")
  .classed("svg-container", true) //container class to make it responsive
  .append("svg")
  .attr("preserveAspectRatio", "xMinYMin meet")
  .attr(
    "viewBox",
    "0 0 " +
      (w + margin.left + margin.right) +
      " " +
      (h + 200 + margin.top + margin.bottom)
  )
  .classed("svg-content", true) //class to make it responsive
  .append("g")
  .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

// * SETTING THE X/Y SCALE
let xScale = d3
  .scaleBand()
  .domain(
    d3
      .map(data_t, (d: any) => d.rec_gn)
      .keys()
      .sort(d3.ascending)
  )
  .range([0, w])
  // .round(true)
  .padding(1);
let yScale = d3
  .scaleBand()
  .domain(
    d3
      .map(data_t, (d: any) => d.sig_gn)
      .keys()
      .sort(d3.descending)
  )
  .range([h, 0])
  // .round(true)
  .padding(1);

// * SETTING THE PDB/SIG-PROT SCALE
// TODO: DEFINE SCALE FOR PDB ID AND SIGPROT ID
let pdbScale = d3
  .scaleBand()
  .domain(
    d3
      .map(data_t, (d: any) => d.pdb_id)
      .keys()
      .sort(d3.descending)
  )
  .range([120,0])
  .padding(1);

let sigScale = d3
  .scaleBand()
  .domain(
    d3
      .map(data_t, (d: any) => d.pdb_id)
      .keys()
      .sort(d3.descending)
  )
  .range([120,0])
  .padding(1);


// * SETTING THE COLOR SCALE
let colScale = d3
  .scaleOrdinal()
  .domain(int_ty)
  .range(d3.schemeSet3);

// * DEFINING AXIS FOR X/Y AND GRID
let xAxis = d3
  .axisBottom(xScale)
  .tickSize(0)
  .tickPadding(8);
let yAxis = d3
  .axisRight(yScale)
  .tickSize(0)
  .tickPadding(8);
let xAxisGrid = d3
  .axisTop(xScale)
  .tickSize(h - yScale.step())
  .tickFormat(d => "");
let yAxisGrid = d3
  .axisRight(yScale)
  .tickSize(w - xScale.step())
  .tickFormat(d => "");

// * ADD TOOLTIP FUNCTIONALITY
let tip = d3
  .tip()
  .attr("class", "d3-tip")
  .html(function(d) {
    return d.rec_gn + "<br>" + d.sig_gn + "<br>" + d.int_ty;
  });
svg.call(tip);

// * RENDER DATA
let shift_left: number = 7 / 8;
let shift_top: number = 1 / 8;
let scale_size: number = shift_left - shift_top;
let offset: number = 1;

// array for data in infobox
let info_data = [];

svg
  .append("g")
  .attr("id", "interact")
  .selectAll("rects")
  .data(data_t)
  .enter()
  .append("rect")
  .attr("x", function(d: any) {
    return xScale(d.rec_gn) - shift_left * xScale.step() + offset;
  })
  .attr("y", function(d: any) {
    return yScale(d.sig_gn) + shift_top * yScale.step() + offset;
  })
  .attr("rx", function() {
    if (data_t.length < 15) {
      return 5;
    } else {
      return 3;
    }
  })
  .attr("ry", function() {
    if (data_t.length < 15) {
      return 5;
    } else {
      return 3;
    }
  })
  .attr("width", xScale.step() * scale_size)
  .attr("height", yScale.step() * scale_size)
  .attr("fill", function(d: any) {
    if (d.int_ty === undefined) {
      return "none";
    } else {
      return colScale(d.int_ty[0]);
    }
  })
  .on("mouseover", function(d) {
    tip.show(d);
  })
  .on("mouseout", function(d) {
    tip.hide();
  })
  .on("click", function(d) {
    let index;
    // let rect_x = d3.event.target.getAttribute('x')
    // let rect_y = d3.event.target.getAttribute('y')
    // console.log(rect_x, rect_y)

    // https://stackoverflow.com/a/20251369/8160230
    // select the rect under cursor
    let curr = d3.select(this);

    // Determine if current rect was clicked before
    let active = d.active ? false : true;

    // Update whether or not the elements are active
    d.active = active;

    // set style in regards to active
    if (d.active) {
      curr.style("stroke", "black").style("stroke-width", 2);
      info_data.push(d);
    } else {
      curr.style("stroke", "none").style("stroke-width", 2);
      index = info_data.indexOf(d);
      info_data.splice(index, 1);
    }
    infoBoxUpdate();
  });

// * ADD INFOBOX ELEMENT
svg
  .append("g")
  .attr("id", "infobox")
  .attr("transform", "translate(-15," + (int_ty.length + 2) * 20 + ")");

function infoBoxUpdate() {
  // create selection and bind data
  let info_box = d3
    .select("g#infobox")
    .selectAll("text")
    .data(info_data);

  // update existing nodes
  info_box
    .attr("y", function(d, i) {
      return i * 15;
    })
    .attr("text-anchor", "end")
    .attr("class", "legend");

  // create nodes for new data
  info_box
    .enter()
    .append("text")
    .attr("y", function(d, i) {
      return i * 15;
    })
    .attr("text-anchor", "end")
    .attr("class", "legend")
    .text(function(d) {
      return d.rec_gn + " : " + d.sig_gn;
    });

  // discard removed nodes
  info_box.exit().remove();

  // print the data again in case it changed
  info_box.text(function(d) {
    return d.rec_gn + " : " + d.sig_gn;
  });
}

// * ADDING COLOR LEGEND
svg
  .append("g")
  .attr("class", "legendOrdinal")
  .attr("transform", "translate(-30," + yScale.step() + ")");

let legendOrdinal = d3
  .legendColor()
  .cells(int_ty.length)
  .scale(colScale)
  // .cellFilter(function (d) { return d.label !== "undefined" })
  .orient("vertical")
  .labelOffset(-20);

svg
  .select(".legendOrdinal")
  .call(legendOrdinal)
  .selectAll("rect")
  .attr("rx", 3)
  .attr("ry", 3);
svg
  .select(".legendOrdinal")
  .selectAll("text")
  .attr("class", "legend")
  .attr("text-anchor", "end");


// * APPENDING COL TICK ANNOTATION FOR RECEPTOR GNs
svg
  .append("g")
  .attr("id", "recPDB")
  .attr("transform", "translate(" + 0 + "," + h + ")")
  .selectAll("text")
  .data(Object.keys(dataset))
  .enter()
  .append("text")
  .attr("class", "x label")
  .attr("x", 0)
  .attr("y", function(d: any) {
    return pdbScale(d);
  })
  .attr("text-anchor", "end")
  .attr("dy", 75)
  .text(function(d: any) {
    return d;
  });


// * APPENDING ROW TICK ANNOTATION FOR SIGPROT GNs
svg
  .append("g")
  .attr("id", "sigPDB")
  .attr("transform", "translate(" + w + "," + 0 + ")rotate(-45)")
  .selectAll("text")
  .data(Object.keys(dataset))
  .enter()
  .append("text")
  .attr("class", "x label")
  .attr("x", function(d: any) {
    return sigScale(d);
  })
  .attr("y", function(d: any) {
    return sigScale(d);
  })
  .attr("text-anchor", "begin")
  .attr("dx", 45)
  .attr("dy", 45)
  .text(function(d: any) {
    return d;
  });


// * APPENDING AMINOACID SEQUENCE [RECEPTOR]
svg
  .append("g")
  .attr("id", "recAA")
  .attr("transform", "translate(" + -xScale.step() / 2 + "," + h + ")")
  .selectAll("text")
  .data(data_t_rec)
  .enter()
  .append("text")
  .attr("class", "res_label")
  .attr("x", function(d: any) {
    return xScale(d.rec_gn);
  })
  .attr("y", function(d: any) {
    return pdbScale(d.pdb_id);
  })
  .attr("text-anchor", "middle")
  .attr("dy", 75)
  .text(function(d: any) {
    return d.rec_aa;
  });

// * APPENDING AMINOACID SEQUENCE [SIGPROT]
svg
  .append("g")
  .attr("id", "sigAA")
  .attr(
    "transform",
    "translate(" + (w + (1 / 3) * margin.right) + "," + yScale.step() / 2 + ")"
  )
  .selectAll("text")
  .data(data_t_sig)
  .enter()
  .append("text")
  .attr("class", "res_label")
  .attr("x", function(d: any) {
    return sigScale(d.pdb_id);
  })
  .attr("y", function(d: any) {
    return yScale(d.sig_gn);
  })
  .attr("text-anchor", "middle")
  .attr("dy", 5)
  .text(function(d: any) {
    return d.sig_aa;
  });

// * AMINOACID SEQUENCE BOX
let seq_rect_h = 20;
// d3.select("g#recAA")
//   .append("rect")
//   .style("stroke", "black")
//   .style("fill", "none")
//   .attr("x", yScale.step() / 2)
//   .attr("y", 60)
//   .attr("width", w - xScale.step())
//   .attr("height", seq_rect_h);

// d3.select("g#sigAA")
//   .append("rect")
//   .style("stroke", "black")
//   .style("fill", "none")
//   .attr("x", -seq_rect_h / 2)
//   .attr("y", yScale.step() / 2)
//   .attr("width", seq_rect_h)
//   .attr("height", h - yScale.step());

// * DRAWING AXES
svg
  .append("g")
  .attr("class", "x axis")
  .attr("transform", "translate(" + -xScale.step() / 2 + "," + h + ")")
  .call(xAxis)
  .selectAll("text")
  .attr("text-anchor", "end")
  .attr("font-size", "12px")
  .attr("dx", "-5px")
  .attr("dy", "-5px")
  .attr("transform", "rotate(-90)");

svg
  .append("g")
  .attr("class", "y axis")
  .attr(
    "transform",
    "translate(" + (w - xScale.step()) + "," + yScale.step() / 2 + ")"
  )
  .call(yAxis)
  .selectAll("text")
  .attr("font-size", "12px");

// * DRAWING GRIDLINES
svg
  .append("g")
  .attr("class", "x grid")
  .attr("transform", "translate(" + 0 + "," + h + ")")
  .call(xAxisGrid);

svg
  .append("g")
  .attr("class", "y grid")
  .attr("transform", "translate(" + 0 + "," + yScale.step() + ")")
  .call(yAxisGrid);

// * ADDITIONAL FIGURE LINES
// top x line
svg
  .append("line")
  .style("stroke", "black")
  .attr("x1", 0)
  .attr("y1", yScale.step())
  .attr("x2", w - xScale.step())
  .attr("y2", yScale.step());

// left y line
svg
  .append("line")
  .style("stroke", "black")
  .attr("x1", 0)
  .attr("y1", yScale.step())
  .attr("x2", 0)
  .attr("y2", h);

// * ADD AXIS LABELS
svg
  .append("text")
  .attr("class", "x label")
  .attr("text-anchor", "end")
  .attr("x", 0)
  .attr("y", h + 15)
  .text("GPCR");

svg
  .append("text")
  .attr("class", "y label")
  .attr("text-anchor", "begin")
  .attr("x", w - 0.8 * xScale.step())
  .attr("y", 0.8 * yScale.step())
  .text("G-Protein");


 $(document).ready( function () {
  $('#table_data').DataTable({
    data: data_t,
    columns: [
      // { data: 'int_ty' },
      { data: 'pdb_id' },
      { data: 'rec_aa' },
      { data: 'rec_gn' },
      { data: 'rec_pos' },
      { data: 'sig_aa' },
      { data: 'sig_gn' }
    ]
  });
});
