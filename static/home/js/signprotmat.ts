// * CONSTANTS
let margin = { top: 0, right: 150, bottom: 180, left: 130 };
let w = 650 - margin.left - margin.right,
  h = 550 - margin.top - margin.bottom;

// * DATA
// let dataset = interactions;
let dataset = {
  '3sn6' : [
  ['R','139','34.51x51','A','376', ["edge-to-face", "face-to-edge", "hydrophobic", "van-der-waals"]],
  ['R','233','5.72x72','A','394', ["hydrophobic"]],
  ['R','233','5.72x72','A','385', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['R','135','3.54x54','A','388', ["hydrophobic", "van-der-waals"]],
  ['R','135','3.54x54','A','393', ["hydrophobic", "van-der-waals"]],
  ['R','62','12.48x48','B','312', ["polar-sidechain-sidechain"]],
  ['R','229','5.68x68','A','388', ["hydrophobic"]],
  ['R','63','12.49x49','B','312', ["polar-backbone-sidechain"]],
  ['R','142','34.54x54','A','387', ["polar-sidechain-sidechain"]],
  ['R','274','6.36x36','A','391', ["polar-sidechain-backbone"]],
  ['R','239','-','A','350', ["hydrophobic", "polar-sidechain-sidechain"]],
  ['R','232','5.71x71','A','381', ["hydrophobic", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
  ['R','135','3.54x54','A','384', ["polar-backbone-sidechain", "van-der-waals"]],
  ['R','139','34.51x51','A','383', ["hydrophobic"]],
  ['R','142','34.54x54','A','383', ["hydrophobic"]],
  ['R','235','5.74x74','A','323', ["polar-backbone-sidechain"]],
  ['R','271','6.33x33','A','393', ["hydrophobic", "van-der-waals"]],
  ['R','225','5.64x64','A','380', ["polar-sidechain-sidechain"]],
  ['R','138','34.50x50','A','380', ["hydrophobic"]],
  ['R','233','5.72x72','A','358', ["hydrophobic", "van-der-waals"]],
  ['R','139','34.51x51','A','217', ["hydrophobic", "van-der-waals"]],
  ['R','138','34.50x50','A','384', ["hydrophobic", "van-der-waals"]],
  ['R','228','5.67x67','A','381', ["polar-sidechain-sidechain"]],
  ['R','239','-','A','347', ["hydrophobic", "polar-sidechain-backbone"]],
  ['R','139','34.51x51','A','379', ["hydrophobic"]],
  ['R','135','3.54x54','A','387', ["hydrophobic", "polar-backbone-sidechain"]],
  ['R','141','34.53x53','A','387', ["edge-to-face", "face-to-edge", "pi-cation", "hydrophobic"]],
  ['R','140','34.52x52','A','380', ["polar-sidechain-sidechain"]],
  ['R','239','-','A','346', ["hydrophobic"]],
  ['R','136','3.55x55','A','380', ["polar-backbone-sidechain", "van-der-waals"]],
  ['R','274','6.36x36','A','392', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['R','225','5.64x64','A','384', ["polar-sidechain-sidechain", "van-der-waals"]],
  ['R','230','5.69x69','A','394', ["hydrophobic", "van-der-waals"]],
  ['R','229','5.68x68','A','384', ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','239','-','A','343', ["polar-sidechain-backbone"]],
  ['R','138','34.50x50','A','387', ["hydrophobic"]],
  ['R','135','3.54x54','A','391', ["hydrophobic"]],
  ['R','139','34.51x51','A','41', ["edge-to-face", "face-to-edge", "pi-cation", "hydrophobic"]],
  ['R','130','3.49x49','A','391', ["polar-sidechain-sidechain"]],
  ['R','138','34.50x50','A','383', ["hydrophobic", "van-der-waals"]],
  ['R','275','6.37x37','A','393', ["hydrophobic", "van-der-waals"]],
  ['R','228','5.67x67','A','384', ["polar-sidechain-sidechain"]],
  ['R','229','5.68x68','A','385', ["hydrophobic", "polar-backbone-sidechain", "polar-sidechain-sidechain", "h-bond"]],
  ['R','134','3.53x53','A','387', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['R','226','5.65x65','A','388', ["hydrophobic", "van-der-waals"]],
  ['R','274','6.36x36','A','393', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['R','139','34.51x51','A','380', ["hydrophobic", "van-der-waals"]],
  ['R','143','34.55x55','A','39', ["hydrophobic", "van-der-waals"]],
  ['R','222','5.61x61','A','393', ["hydrophobic", "van-der-waals"]],
  ['R','229','5.68x68','A','381', ["polar-sidechain-backbone", "polar-sidechain-sidechain"]],
  ['R','131','3.50x50','A','391', ["cation-pi", "hydrophobic", "van-der-waals"]],
  ],
  '4x1h' : [
  ['A','135','3.50x50','C','347', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals", "water-mediated", "water-mediated"]],
  ['A','72','2.39x39','C','347', ["water-mediated", "water-mediated"]],
  ['A','310','8.47x47','C','346', ["water-mediated"]],
  ['A','246','6.29x29','C','341', ["hydrophobic", "van-der-waals"]],
  ['A','135','3.50x50','C','346', ["water-mediated"]],
  ['A','309','7.56x56','C','348', ["water-mediated", "water-mediated"]],
  ['A','305','7.52x52','C','348', ["water-mediated"]],
  ['A','138','3.53x53','C','343', ["hydrophobic", "polar-backbone-sidechain", "water-mediated"]],
  ['A','139','3.54x54','C','340', ["hydrophobic", "water-mediated"]],
  ['A','73','2.40x40','C','346', ["water-mediated", "water-mediated"]],
  ['A','306','7.53x53','C','348', ["water-mediated", "water-mediated"]],
  ['A','312','8.49x49','C','346', ["water-mediated"]],
  ['A','139','3.54x54','C','344', ["hydrophobic"]],
  ['A','250','6.33x33','C','344', ["hydrophobic"]],
  ['A','250','6.33x33','C','349', ["hydrophobic"]],
  ['A','135','3.50x50','C','348', ["water-mediated", "water-mediated"]],
  ['A','243','6.26x26','C','341', ["hydrophobic"]],
  ['A','141','3.56x56','C','343', ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals", "water-mediated"]],
  ['A','312','8.49x49','C','348', ["water-mediated"]],
  ['A','141','3.56x56','C','340', ["water-mediated"]],
  ['A','229','5.64x64','C','340', ["hydrophobic"]],
  ['A','73','2.40x40','C','347', ["water-mediated", "water-mediated"]],
  ['A','71','2.38x38','C','346', ["water-mediated"]],
  ['A','242','6.25x25','C','350', ["hydrophobic", "van-der-waals"]],
  ['A','246','6.29x29','C','344', ["hydrophobic"]],
  ['A','311','8.48x48','C','350', ["polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
  ['A','226','5.61x61','C','344', ["hydrophobic"]],
  ['A','70','2.37x37','C','346', ["water-mediated"]],
  ['A','249','6.32x32','C','350', ["hydrophobic"]],
  ['A','139','3.54x54','C','343', ["hydrophobic", "water-mediated"]],
  ['A','226','5.61x61','C','349', ["hydrophobic"]],
  ['A','230','5.65x65','C','344', ["hydrophobic"]],
  ['A','310','8.47x47','C','348', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals", "water-mediated", "water-mediated", "water-mediated", "water-mediated"]],
  ['A','233','5.68x68','C','340', ["hydrophobic"]],
  ['A','246','6.29x29','C','350', ["hydrophobic", "van-der-waals"]],
  ['A','138','3.53x53','C','340', ["water-mediated"]],
  ['A','310','8.47x47','C','347', ["polar-sidechain-backbone", "water-mediated"]],
  ['A','135','3.50x50','C','349', ["hydrophobic", "van-der-waals"]],
  ['A','73','2.40x40','C','348', ["water-mediated"]],
  ['A','72','2.39x39','C','346', ["hydrophobic", "van-der-waals", "water-mediated", "water-mediated"]],
  ['A','245','6.28x28','C','350', ["hydrophobic"]],
  ['A','253','6.36x36','C','348', ["water-mediated"]],
  ['A','312','8.49x49','C','347', ["water-mediated"]],
  ['A','242','6.25x25','C','341', ["hydrophobic"]],
  ],
  '5g53' : [
  ['A','106','3.54x54','C','388', ["hydrophobic", "van-der-waals"]],
  ['A','111','34.52x52','C','215', ["polar-sidechain-backbone", "van-der-waals"]],
  ['A','207','5.68x68','C','384', ["polar-sidechain-sidechain", "h-bond"]],
  ['A','200','5.61x61','C','388', ["hydrophobic"]],
  ['A','110','34.51x51','C','41', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['A','109','34.50x50','C','384', ["hydrophobic"]],
  ['A','107','3.55x55','C','380', ["polar-backbone-sidechain", "van-der-waals"]],
  ['A','293','8.48x48','C','392', ["hydrophobic", "polar-backbone-sidechain"]],
  ['A','200','5.61x61','C','393', ["hydrophobic"]],
  ['A','204','5.65x65','C','388', ["hydrophobic"]],
  ['A','210','5.71x71','C','381', ["polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
  ['A','291','7.56x56','C','392', ["polar-backbone-sidechain"]],
  ['A','110','34.51x51','C','379', ["hydrophobic"]],
  ['A','110','34.51x51','C','376', ["hydrophobic", "van-der-waals"]],
  ['A','106','3.54x54','C','391', ["hydrophobic"]],
  ['A','207','5.68x68','C','385', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['A','207','5.68x68','C','381', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['A','105','3.53x53','C','387', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['A','110','34.51x51','C','217', ["hydrophobic", "van-der-waals"]],
  ['A','203','5.64x64','C','384', ["hydrophobic", "polar-backbone-sidechain"]],
  ['A','110','34.51x51','C','380', ["hydrophobic"]],
  ['A','106','3.54x54','C','384', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['A','111','34.52x52','C','217', ["hydrophobic", "van-der-waals"]],
  ['A','291','7.56x56','C','391', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['A','106','3.54x54','C','387', ["hydrophobic"]],
  ['A','111','34.52x52','C','380', ["polar-sidechain-sidechain"]],
  ['A','110','34.51x51','C','219', ["hydrophobic"]],
  ['A','211','5.72x72','C','358', ["hydrophobic", "polar-sidechain-sidechain"]],
  ['A','110','34.51x51','C','383', ["hydrophobic"]],
  ['A','108','3.56x56','C','380', ["polar-backbone-sidechain"]],
  ['A','207','5.68x68','C','388', ["hydrophobic"]],
  ['A','231','6.33x33','C','393', ["hydrophobic"]],
  ['A','109','34.50x50','C','380', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['A','208','5.69x69','C','394', ["hydrophobic", "van-der-waals"]],
  ['A','207','5.68x68','C','360', ["polar-sidechain-sidechain"]],
  ['A','109','34.50x50','C','383', ["hydrophobic", "van-der-waals"]],
  ['A','112','34.53x53','C','387', ["edge-to-face", "face-to-edge", "pi-cation", "hydrophobic"]],
  ['A','235','6.37x37','C','393', ["hydrophobic", "van-der-waals"]],
  ['A','203','5.64x64','C','388', ["hydrophobic"]],
  ['A','102','3.50x50','C','391', ["cation-pi", "hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['A','296','8.51x51','C','392', ["polar-sidechain-sidechain", "van-der-waals"]],
  ['A','227','6.29x29','C','394', ["hydrophobic", "polar-sidechain-sidechain"]],
  ],
  '5g53' : [
  ['B','200','5.61x61','D','388', ["hydrophobic"]],
  ['B','203','5.64x64','D','388', ["hydrophobic"]],
  ['B','235','6.37x37','D','393', ["hydrophobic", "van-der-waals"]],
  ['B','114','34.55x55','D','39', ["hydrophobic", "van-der-waals"]],
  ['B','207','5.68x68','D','384', ["polar-sidechain-sidechain"]],
  ['B','109','34.50x50','D','384', ["hydrophobic"]],
  ['B','293','8.48x48','D','392', ["polar-backbone-sidechain", "van-der-waals"]],
  ['B','110','34.51x51','D','217', ["hydrophobic", "van-der-waals"]],
  ['B','203','5.64x64','D','384', ["hydrophobic", "polar-backbone-sidechain"]],
  ['B','207','5.68x68','D','381', ["polar-sidechain-backbone"]],
  ['B','106','3.54x54','D','384', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['B','113','34.54x54','D','39', ["hydrophobic", "van-der-waals"]],
  ['B','296','8.51x51','D','392', ["polar-sidechain-sidechain", "van-der-waals"]],
  ['B','207','5.68x68','D','385', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['B','102','3.50x50','D','391', ["hydrophobic", "polar-backbone-sidechain"]],
  ['B','111','34.52x52','D','380', ["polar-sidechain-sidechain"]],
  ['B','110','34.51x51','D','41', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['B','231','6.33x33','D','393', ["hydrophobic", "van-der-waals"]],
  ['B','106','3.54x54','D','388', ["hydrophobic", "van-der-waals"]],
  ['B','110','34.51x51','D','379', ["hydrophobic"]],
  ['B','204','5.65x65','D','393', ["hydrophobic", "van-der-waals"]],
  ['B','204','5.65x65','D','388', ["hydrophobic"]],
  ['B','105','3.53x53','D','387', ["hydrophobic", "polar-backbone-sidechain"]],
  ['B','110','34.51x51','D','219', ["hydrophobic"]],
  ['B','291','7.56x56','D','392', ["polar-backbone-sidechain"]],
  ['B','110','34.51x51','D','376', ["hydrophobic", "van-der-waals"]],
  ['B','111','34.52x52','D','217', ["hydrophobic", "van-der-waals"]],
  ['B','106','3.54x54','D','387', ["hydrophobic"]],
  ['B','207','5.68x68','D','388', ["hydrophobic"]],
  ['B','110','34.51x51','D','380', ["hydrophobic"]],
  ['B','111','34.52x52','D','215', ["polar-sidechain-backbone", "van-der-waals"]],
  ['B','109','34.50x50','D','383', ["hydrophobic", "van-der-waals"]],
  ['B','200','5.61x61','D','393', ["hydrophobic", "van-der-waals"]],
  ['B','109','34.50x50','D','380', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['B','112','34.53x53','D','387', ["edge-to-face", "face-to-edge", "pi-cation", "hydrophobic"]],
  ['B','105','3.53x53','D','391', ["van-der-waals"]],
  ['B','106','3.54x54','D','391', ["hydrophobic"]],
  ['B','107','3.55x55','D','380', ["polar-backbone-sidechain", "van-der-waals"]],
  ['B','110','34.51x51','D','383', ["hydrophobic"]],
  ['B','108','3.56x56','D','380', ["polar-backbone-sidechain"]],
  ],
  '5uz7' : [
  ['R','326','5.64x64','A','384', ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','415','8.67x67','B','307', ["polar-sidechain-backbone", "van-der-waals"]],
  ['R','408','8.60x60','B','311', ["polar-sidechain-backbone", "van-der-waals"]],
  ['R','331','-','A','358', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['R','253','-','A','41', ["polar-backbone-sidechain", "van-der-waals"]],
  ['R','348','6.45x45','A','393', ["hydrophobic"]],
  ['R','329','-','A','385', ["polar-backbone-sidechain"]],
  ['R','408','8.60x60','B','309', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['R','243','3.53x53','A','391', ["van-der-waals"]],
  ['R','247','3.57x57','A','387', ["hydrophobic", "van-der-waals"]],
  ['R','249','3.59x59','A','384', ["polar-backbone-sidechain"]],
  ['R','411','8.63x63','B','307', ["hydrophobic"]],
  ['R','394','7.60x60','A','392', ["polar-backbone-sidechain"]],
  ['R','348','6.45x45','A','392', ["hydrophobic"]],
  ['R','323','5.61x61','A','388', ["hydrophobic", "van-der-waals"]],
  ['R','404','8.56x56','B','312', ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','254','-','A','387', ["polar-sidechain-sidechain"]],
  ['R','180','2.46x46','A','391', ["cation-pi", "hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['R','252','-','A','384', ["hydrophobic", "van-der-waals"]],
  ['R','415','8.67x67','B','44', ["hydrophobic", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
  ['R','396','8.48x48','A','356', ["polar-sidechain-sidechain"]],
  ['R','326','5.64x64','A','385', ["polar-backbone-sidechain"]],
  ['R','252','-','A','380', ["hydrophobic", "van-der-waals"]],
  ['R','252','-','A','383', ["hydrophobic", "van-der-waals"]],
  ['R','408','8.60x60','B','310', ["polar-sidechain-backbone"]],
  ['R','247','3.57x57','A','391', ["hydrophobic"]],
  ['R','253','-','A','217', ["hydrophobic"]],
  ['R','327','5.65x65','A','394', ["hydrophobic", "van-der-waals"]],
  ['R','345','6.42x42','A','393', ["hydrophobic", "van-der-waals"]],
  ['R','180','2.46x46','A','390', ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
  ['R','248','3.58x58','A','387', ["hydrophobic"]],
  ['R','184','2.50x50','A','391', ["hydrophobic"]],
  ['R','244','3.54x54','A','391', ["hydrophobic"]],
  ['R','330','-','A','358', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['R','248','3.58x58','A','388', ["hydrophobic"]],
  ['R','396','8.48x48','A','392', ["polar-backbone-sidechain"]],
  ['R','326','5.64x64','A','380', ["polar-sidechain-sidechain"]],
  ['R','248','3.58x58','A','384', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ],
  '5vai' : [
  ['R','415','8.56x56','B','292', ["cation-pi", "hydrophobic", "van-der-waals"]],
  ['R','262','4.38x39','A','35', ["polar-backbone-sidechain", "van-der-waals"]],
  ['R','356','6.45x45','A','393', ["hydrophobic"]],
  ['R','264','4.40x41','A','35', ["polar-backbone-sidechain"]],
  ['R','176','2.46x46','A','390', ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','261','4.37x38','A','35', ["hydrophobic", "polar-backbone-sidechain", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','419','8.60x60','B','310', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['R','263','4.39x40','A','31', ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','338','-','A','359', ["polar-sidechain-backbone"]],
  ['R','405','7.60x60','A','392', ["polar-backbone-sidechain", "van-der-waals"]],
  ['R','256','3.59x59','A','384', ["polar-backbone-sidechain"]],
  ['R','339','-','A','385', ["polar-backbone-sidechain"]],
  ['R','353','6.42x42','A','393', ["hydrophobic"]],
  ['R','419','8.60x60','B','293', ["polar-sidechain-sidechain"]],
  ['R','339','-','A','358', ["hydrophobic", "van-der-waals"]],
  ['R','408','8.49x49','A','390', ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
  ['R','419','8.60x60','B','312', ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','255','3.58x58','A','384', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['R','171','12.49x49','B','312', ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','352','6.41x41','A','393', ["polar-sidechain-backbone", "van-der-waals"]],
  ['R','406','8.47x47','A','390', ["polar-sidechain-backbone"]],
  ['R','256','3.59x59','A','380', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['R','339','-','A','394', ["hydrophobic"]],
  ['R','419','8.60x60','B','309', ["hydrophobic", "van-der-waals"]],
  ['R','331','5.61x61','A','388', ["hydrophobic", "van-der-waals"]],
  ['R','407','8.48x48','A','392', ["polar-backbone-sidechain"]],
  ['R','412','8.53x53','B','312', ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','254','3.57x57','A','387', ["hydrophobic", "van-der-waals"]],
  ['R','352','6.41x41','A','392', ["polar-sidechain-backbone"]],
  ['R','263','4.39x40','A','35', ["polar-backbone-sidechain", "van-der-waals"]],
  ['R','342','-','A','350', ["polar-sidechain-backbone", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
  ['R','176','2.46x46','A','391', ["hydrophobic"]],
  ['R','334','5.64x64','A','385', ["polar-backbone-sidechain", "van-der-waals"]],
  ['R','359','6.48x48','A','391', ["hydrophobic", "van-der-waals"]],
  ['R','255','3.58x58','A','380', ["polar-backbone-sidechain"]],
  ['R','334','5.64x64','A','381', ["hydrophobic", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
  ['R','352','6.41x41','A','394', ["polar-sidechain-backbone"]],
  ['R','180','2.50x50','A','391', ["cation-pi", "hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','356','6.45x45','A','391', ["hydrophobic"]],
  ['R','338','-','A','360', ["hydrophobic", "van-der-waals"]],
  ['R','415','8.56x56','B','291', ["polar-sidechain-backbone"]],
  ['R','419','8.60x60','B','311', ["polar-sidechain-backbone"]],
  ['R','406','8.47x47','A','392', ["hydrophobic", "polar-backbone-sidechain", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','402','7.57x57','A','392', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['R','401','7.56x56','A','392', ["polar-backbone-sidechain", "van-der-waals"]],
  ['R','251','3.54x54','A','391', ["hydrophobic", "van-der-waals"]],
  ['R','419','8.60x60','B','292', ["polar-sidechain-backbone"]],
  ['R','170','12.48x48','B','52', ["polar-sidechain-sidechain"]],
  ],
  '6b3j' : [
  ['R','356','6.45x45','A','393', ["hydrophobic"]],
  ['R','262','4.38x39','A','38', ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','250','3.53x53','A','391', ["van-der-waals"]],
  ['R','334','5.64x64','A','384', ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','327','5.57x57','A','393', ["hydrophobic"]],
  ['R','258','-','A','383', ["hydrophobic"]],
  ['R','261','-','A','35', ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','176','2.46x46','A','391', ["cation-pi", "hydrophobic"]],
  ['R','352','6.41x41','A','393', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['R','334','5.64x64','A','385', ["hydrophobic", "polar-backbone-sidechain"]],
  ['R','262','4.38x39','A','34', ["polar-sidechain-sidechain"]],
  ['R','423','8.64x64','B','44', ["hydrophobic"]],
  ['R','259','-','A','217', ["hydrophobic"]],
  ['R','334','5.64x64','A','381', ["hydrophobic", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
  ['R','352','6.41x41','A','394', ["polar-sidechain-backbone"]],
  ['R','348','6.37x37','A','394', ["hydrophobic"]],
  ['R','254','3.57x57','A','387', ["hydrophobic", "van-der-waals"]],
  ['R','180','2.50x50','A','391', ["hydrophobic"]],
  ['R','331','5.61x61','A','393', ["hydrophobic"]],
  ['R','334','5.64x64','A','388', ["hydrophobic"]],
  ['R','415','8.56x56','B','312', ["polar-sidechain-sidechain", "h-bond"]],
  ['R','255','3.58x58','A','384', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['R','171','12.49x49','B','312', ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','251','3.54x54','A','391', ["hydrophobic", "van-der-waals"]],
  ['R','331','5.61x61','A','388', ["hydrophobic"]],
  ['R','176','2.46x46','A','390', ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain"]],
  ['R','331','5.61x61','A','394', ["hydrophobic"]],
  ['R','255','3.58x58','A','387', ["hydrophobic"]],
  ['R','419','8.60x60','B','309', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['R','255','3.58x58','A','388', ["hydrophobic"]],
  ['R','407','8.48x48','A','392', ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
  ],
  '6cmo' : [
  ['R','246','6.29x29','A','348', ["hydrophobic"]],
  ['R','141','3.56x56','A','193', ["polar-sidechain-sidechain"]],
  ['R','237','5.72x72','A','320', ["polar-sidechain-sidechain"]],
  ['R','241','6.24x24','A','318', ["polar-backbone-sidechain"]],
  ['R','245','6.28x28','A','354', ["hydrophobic"]],
  ['R','242','6.25x25','A','354', ["hydrophobic", "van-der-waals"]],
  ['R','249','6.32x32','A','354', ["hydrophobic"]],
  ['R','246','6.29x29','A','354', ["hydrophobic", "van-der-waals"]],
  ['R','240','-','A','345', ["polar-sidechain-sidechain", "h-bond"]],
  ['R','239','-','A','318', ["polar-backbone-sidechain"]],
  ['R','135','3.50x50','A','353', ["hydrophobic"]],
  ['R','66','12.48x48','B','312', ["polar-sidechain-sidechain"]],
  ['R','240','-','A','318', ["hydrophobic", "polar-backbone-sidechain", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
  ['R','139','3.54x54','A','347', ["hydrophobic"]],
  ['R','237','5.72x72','A','341', ["hydrophobic"]],
  ['R','243','6.26x26','A','341', ["polar-sidechain-sidechain"]],
  ['R','311','8.48x48','A','353', ["polar-sidechain-backbone"]],
  ['R','311','8.48x48','A','352', ["polar-sidechain-backbone"]],
  ['R','242','6.25x25','A','315', ["polar-sidechain-backbone"]],
  ['R','253','6.36x36','A','353', ["hydrophobic", "van-der-waals"]],
  ['R','249','6.32x32','A','353', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['R','147','34.55x55','A','31', ["polar-sidechain-backbone"]],
  ['R','147','34.55x55','A','32', ["hydrophobic", "polar-sidechain-backbone"]],
  ['R','309','7.56x56','A','352', ["van-der-waals"]],
  ['R','139','3.54x54','A','348', ["hydrophobic"]],
  ['R','311','8.48x48','A','349', ["polar-sidechain-backbone"]],
  ['R','250','6.33x33','A','353', ["hydrophobic"]],
  ['R','310','8.47x47','A','352', ["hydrophobic", "polar-sidechain-backbone"]],
  ['R','239','-','A','320', ["hydrophobic", "van-der-waals"]],
  ['R','311','8.48x48','A','354', ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
  ],
  '6d9h' : [
  ['R','45','2.40x40','A','351', ["van-der-waals"]],
  ['R','38','12.48x48','B','335', ["hydrophobic", "van-der-waals"]],
  ['R','207','5.65x65','A','354', ["hydrophobic"]],
  ['R','232','6.33x33','A','354', ["hydrophobic", "van-der-waals"]],
  ['R','113','34.51x51','A','337', ["hydrophobic"]],
  ['R','228','6.29x29','A','316', ["polar-sidechain-sidechain"]],
  ['R','294','8.49x49','A','351', ["hydrophobic", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
  ['R','203','5.61x61','A','354', ["hydrophobic"]],
  ['R','37','1.60x60','B','312', ["polar-sidechain-sidechain"]],
  ['R','113','34.51x51','A','344', ["hydrophobic"]],
  ['R','228','6.29x29','A','355', ["cation-pi", "hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','210','5.68x68','A','345', ["hydrophobic", "van-der-waals"]],
  ['R','210','5.68x68','A','342', ["polar-sidechain-sidechain", "van-der-waals"]],
  ['R','207','5.65x65','A','349', ["hydrophobic"]],
  ['R','292','8.47x47','A','353', ["hydrophobic"]],
  ['R','113','34.51x51','A','341', ["hydrophobic"]],
  ['R','224','6.25x25','A','319', ["polar-sidechain-sidechain"]],
  ['R','236','6.37x37','A','354', ["hydrophobic", "van-der-waals"]],
  ['R','224','6.25x25','A','316', ["hydrophobic", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
  ['R','112','34.50x50','A','345', ["hydrophobic", "van-der-waals"]],
  ['R','105','3.50x50','A','354', ["hydrophobic"]],
  ['R','109','3.54x54','A','349', ["hydrophobic", "van-der-waals"]],
  ['R','113','34.51x51','A','195', ["hydrophobic"]],
  ['R','292','8.47x47','A','352', ["hydrophobic", "van-der-waals"]],
  ['R','213','5.71x71','A','342', ["polar-sidechain-sidechain"]],
  ['R','231','6.32x32','A','355', ["polar-sidechain-backbone"]],
  ['R','112','34.50x50','A','344', ["hydrophobic", "van-der-waals"]],
  ['R','291','7.56x56','A','353', ["hydrophobic"]],
  ['R','301','8.56x56','B','312', ["polar-sidechain-backbone"]],
  ['R','232','6.33x33','A','355', ["hydrophobic"]],
  ['R','42','2.37x37','A','351', ["polar-sidechain-sidechain"]],
  ['R','105','3.50x50','A','352', ["hydrophobic", "polar-sidechain-backbone"]],
  ['R','112','34.50x50','A','348', ["polar-backbone-sidechain"]],
  ['R','108','3.53x53','A','352', ["polar-sidechain-sidechain", "van-der-waals"]],
  ['R','108','3.53x53','A','351', ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','108','3.53x53','A','348', ["hydrophobic", "polar-backbone-sidechain", "polar-sidechain-sidechain"]],
  ['R','211','5.69x69','A','346', ["hydrophobic"]],
  ['R','301','8.56x56','B','292', ["cation-pi", "hydrophobic", "van-der-waals"]],
  ],
  '6dde' : [
  ['R','263','-','A','341', ["polar-backbone-sidechain"]],
  ['R','176','34.54x54','A','343', ["hydrophobic", "van-der-waals"]],
  ['R','173','34.51x51','A','194', ["hydrophobic"]],
  ['R','264','-','A','341', ["hydrophobic", "polar-backbone-sidechain"]],
  ['R','258','5.64x64','A','344', ["hydrophobic", "van-der-waals"]],
  ['R','268','6.23x23','A','315', ["polar-backbone-sidechain"]],
  ['R','264','-','A','316', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['R','281','6.36x36','A','353', ["hydrophobic", "van-der-waals"]],
  ['R','271','6.26x26','A','314', ["polar-sidechain-backbone", "van-der-waals"]],
  ['R','177','34.55x55','A','32', ["polar-backbone-sidechain", "polar-sidechain-sidechain"]],
  ['R','263','-','A','319', ["polar-sidechain-backbone"]],
  ['R','165','3.50x50','A','353', ["hydrophobic"]],
  ['R','278','6.33x33','A','354', ["hydrophobic", "van-der-waals"]],
  ['R','165','3.50x50','A','351', ["polar-sidechain-sidechain", "van-der-waals"]],
  ['R','277','6.32x32','A','353', ["polar-sidechain-backbone"]],
  ['R','259','5.65x65','A','348', ["hydrophobic"]],
  ['R','172','34.50x50','A','343', ["hydrophobic", "van-der-waals"]],
  ['R','255','5.61x61','A','353', ["hydrophobic"]],
  ['R','103','2.39x39','A','351', ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','176','34.54x54','A','194', ["hydrophobic"]],
  ['R','169','3.54x54','A','348', ["hydrophobic", "van-der-waals"]],
  ['R','103','2.39x39','A','350', ["polar-sidechain-backbone"]],
  ['R','264','-','A','345', ["hydrophobic"]],
  ['R','262','5.68x68','A','344', ["hydrophobic"]],
  ['R','172','34.50x50','A','344', ["hydrophobic", "van-der-waals"]],
  ['R','278','6.33x33','A','348', ["hydrophobic"]],
  ['R','182','4.40x40','A','24', ["hydrophobic", "polar-sidechain-sidechain"]],
  ['R','176','34.54x54','A','32', ["hydrophobic", "polar-backbone-sidechain"]],
  ['R','278','6.33x33','A','353', ["hydrophobic"]],
  ['R','271','6.26x26','A','315', ["hydrophobic", "van-der-waals"]],
  ['R','259','5.65x65','A','344', ["hydrophobic"]],
  ['R','179','34.57x57','A','347', ["polar-sidechain-sidechain", "van-der-waals"]],
  ['R','340','8.47x47','A','352', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['R','263','-','A','320', ["cation-pi", "hydrophobic", "van-der-waals"]],
  ['R','168','3.53x53','A','347', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['R','173','34.51x51','A','336', ["hydrophobic", "van-der-waals"]],
  ['R','270','6.25x25','A','315', ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','172','34.50x50','A','340', ["hydrophobic"]],
  ['R','173','34.51x51','A','193', ["hydrophobic", "van-der-waals"]],
  ['R','262','5.68x68','A','341', ["hydrophobic", "van-der-waals"]],
  ],
  '6ddf' : [
  ['R','263','-','A','341', ["hydrophobic", "polar-backbone-sidechain"]],
  ['R','176','34.54x54','A','343', ["hydrophobic"]],
  ['R','271','6.26x26','A','317', ["polar-sidechain-backbone"]],
  ['R','263','-','A','319', ["polar-sidechain-backbone", "van-der-waals"]],
  ['R','259','5.65x65','A','348', ["hydrophobic", "van-der-waals"]],
  ['R','271','6.26x26','A','314', ["polar-sidechain-backbone", "van-der-waals"]],
  ['R','177','34.55x55','A','32', ["polar-sidechain-sidechain"]],
  ['R','341','8.48x48','A','354', ["hydrophobic"]],
  ['R','165','3.50x50','A','353', ["hydrophobic"]],
  ['R','165','3.50x50','A','351', ["hydrophobic", "polar-sidechain-backbone"]],
  ['R','271','6.26x26','A','315', ["hydrophobic", "van-der-waals"]],
  ['R','103','2.39x39','A','351', ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain"]],
  ['R','263','-','A','318', ["van-der-waals"]],
  ['R','179','34.57x57','A','347', ["polar-sidechain-sidechain"]],
  ['R','103','2.39x39','A','350', ["polar-sidechain-backbone", "van-der-waals"]],
  ['R','164','3.49x49','A','351', ["polar-sidechain-sidechain"]],
  ['R','278','6.33x33','A','348', ["hydrophobic"]],
  ['R','340','8.47x47','A','352', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['R','173','34.51x51','A','194', ["hydrophobic"]],
  ['R','278','6.33x33','A','353', ["hydrophobic"]],
  ['R','172','34.50x50','A','347', ["polar-backbone-sidechain"]],
  ['R','174','34.52x52','A','193', ["polar-sidechain-sidechain"]],
  ['R','179','34.57x57','A','351', ["polar-sidechain-sidechain", "van-der-waals"]],
  ['R','263','-','A','320', ["cation-pi", "hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','172','34.50x50','A','343', ["hydrophobic", "van-der-waals"]],
  ['R','172','34.50x50','A','340', ["hydrophobic", "van-der-waals"]],
  ['R','169','3.54x54','A','348', ["hydrophobic", "van-der-waals"]],
  ['R','264','-','A','341', ["polar-backbone-sidechain"]],
  ['R','258','5.64x64','A','344', ["hydrophobic", "van-der-waals"]],
  ['R','176','34.54x54','A','32', ["hydrophobic"]],
  ['R','264','-','A','316', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['R','281','6.36x36','A','353', ["hydrophobic", "van-der-waals"]],
  ['R','277','6.32x32','A','353', ["polar-sidechain-backbone", "van-der-waals"]],
  ['R','173','34.51x51','A','336', ["hydrophobic", "van-der-waals"]],
  ['R','264','-','A','345', ["hydrophobic"]],
  ['R','277','6.32x32','A','354', ["polar-sidechain-backbone", "van-der-waals"]],
  ['R','255','5.61x61','A','353', ["hydrophobic"]],
  ['R','341','8.48x48','A','353', ["polar-sidechain-backbone"]],
  ['R','168','3.53x53','A','351', ["van-der-waals"]],
  ['R','176','34.54x54','A','194', ["hydrophobic"]],
  ['R','262','5.68x68','A','344', ["hydrophobic"]],
  ['R','172','34.50x50','A','344', ["hydrophobic", "van-der-waals"]],
  ['R','340','8.47x47','A','353', ["polar-sidechain-backbone"]],
  ['R','278','6.33x33','A','354', ["hydrophobic", "van-der-waals"]],
  ['R','341','8.48x48','A','352', ["polar-sidechain-backbone", "van-der-waals"]],
  ['R','173','34.51x51','A','193', ["hydrophobic"]],
  ['R','270','6.25x25','A','315', ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
  ['R','168','3.53x53','A','347', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['R','259','5.65x65','A','344', ["hydrophobic"]],
  ['R','262','5.68x68','A','341', ["hydrophobic", "van-der-waals"]],
  ['R','268','6.23x23','A','315', ["polar-backbone-sidechain", "polar-sidechain-sidechain", "van-der-waals"]],
  ],
  '6g79' : [
  ['S','308','6.29x29','A','306', ["polar-sidechain-backbone"]],
  ['S','373','8.47x47','A','340', ["polar-sidechain-backbone"]],
  ['S','151','3.54x54','A','338', ["hydrophobic", "van-der-waals"]],
  ['S','311','6.32x32','A','342', ["polar-sidechain-backbone"]],
  ['S','312','6.33x33','A','343', ["hydrophobic", "van-der-waals"]],
  ['S','154','34.50x50','A','333', ["hydrophobic"]],
  ['S','372','7.56x56','A','342', ["hydrophobic"]],
  ['S','147','3.50x50','A','341', ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
  ['S','239','5.69x69','A','344', ["hydrophobic"]],
  ['S','150','3.53x53','A','337', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['S','155','34.51x51','A','330', ["hydrophobic"]],
  ['S','372','7.56x56','A','341', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['S','315','6.36x36','A','342', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['S','161','34.57x57','A','337', ["polar-sidechain-sidechain"]],
  ['S','231','5.61x61','A','343', ["hydrophobic", "van-der-waals"]],
  ['S','151','3.54x54','A','334', ["hydrophobic"]],
  ['S','315','6.36x36','A','343', ["hydrophobic", "polar-sidechain-backbone"]],
  ['S','316','6.37x37','A','343', ["hydrophobic"]],
  ['S','238','5.68x68','A','338', ["hydrophobic"]],
  ['S','238','5.68x68','A','334', ["hydrophobic", "van-der-waals"]],
  ['S','238','5.68x68','A','344', ["hydrophobic"]],
  ['S','235','5.65x65','A','343', ["hydrophobic"]],
  ['S','311','6.32x32','A','344', ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
  ['S','154','34.50x50','A','334', ["hydrophobic", "van-der-waals"]],
  ['S','235','5.65x65','A','338', ["hydrophobic"]],
  ['S','238','5.68x68','A','331', ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
  ['S','155','34.51x51','A','326', ["hydrophobic"]],
  ['S','308','6.29x29','A','344', ["cation-pi", "hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
  ['S','151','3.54x54','A','337', ["hydrophobic"]],
  ],
  '6gdg' : [
  ['A','106','3.54x54','D','378', ["hydrophobic", "van-der-waals"]],
  ['A','112','34.53x53','D','377', ["edge-to-face", "face-to-edge", "pi-cation", "hydrophobic"]],
  ['A','207','5.68x68','D','374', ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
  ['A','106','3.54x54','D','374', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['A','110','34.51x51','D','373', ["hydrophobic"]],
  ['A','203','5.64x64','D','378', ["hydrophobic"]],
  ['A','109','34.50x50','D','374', ["hydrophobic", "van-der-waals"]],
  ['A','102','3.50x50','D','381', ["cation-pi", "hydrophobic", "van-der-waals"]],
  ['A','207','5.68x68','D','384', ["hydrophobic"]],
  ['A','113','34.54x54','D','39', ["hydrophobic"]],
  ['A','105','3.53x53','D','377', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['A','204','5.65x65','D','378', ["hydrophobic"]],
  ['A','200','5.61x61','D','383', ["hydrophobic", "van-der-waals"]],
  ['A','110','34.51x51','D','369', ["hydrophobic"]],
  ['A','106','3.54x54','D','381', ["hydrophobic"]],
  ['A','109','34.50x50','D','370', ["hydrophobic"]],
  ['A','211','5.72x72','D','348', ["hydrophobic"]],
  ['A','293','8.48x48','D','380', ["polar-backbone-sidechain"]],
  ['A','207','5.68x68','D','375', ["polar-sidechain-backbone"]],
  ['A','291','7.56x56','D','381', ["polar-sidechain-backbone"]],
  ['A','110','34.51x51','D','370', ["hydrophobic"]],
  ['A','113','34.54x54','D','38', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
  ['A','106','3.54x54','D','377', ["hydrophobic", "van-der-waals"]],
  ['A','294','8.49x49','D','380', ["polar-backbone-sidechain"]],
  ['A','110','34.51x51','D','217', ["hydrophobic"]],
  ['A','293','8.48x48','D','382', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['A','207','5.68x68','D','371', ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
  ['A','111','34.52x52','D','217', ["hydrophobic", "van-der-waals"]],
  ['A','207','5.68x68','D','378', ["hydrophobic"]],
  ['A','34','1.60x60','B','312', ["polar-sidechain-sidechain"]],
  ['A','114','34.55x55','D','39', ["hydrophobic"]],
  ['A','230','6.32x32','D','382', ["polar-sidechain-backbone"]],
  ['A','107','3.55x55','D','370', ["polar-backbone-sidechain", "polar-sidechain-sidechain", "van-der-waals"]],
  ['A','234','6.36x36','D','383', ["hydrophobic", "van-der-waals"]],
  ['A','109','34.50x50','D','373', ["hydrophobic", "van-der-waals"]],
  ['A','110','34.51x51','D','366', ["hydrophobic", "van-der-waals"]],
  ['A','36','12.49x49','B','312', ["hydrophobic", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
  ['A','208','5.69x69','D','384', ["hydrophobic"]],
  ['A','35','12.48x48','B','333', ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
  ['A','35','12.48x48','B','335', ["hydrophobic", "van-der-waals"]],
  ['A','35','12.48x48','B','312', ["hydrophobic"]],
  ['A','291','7.56x56','D','382', ["polar-backbone-sidechain", "polar-sidechain-backbone", "van-der-waals"]],
  ['A','231','6.33x33','D','383', ["hydrophobic", "van-der-waals"]],
  ['A','111','34.52x52','D','215', ["polar-sidechain-backbone"]],
  ['A','111','34.52x52','D','216', ["hydrophobic", "polar-sidechain-backbone"]],
  ['A','110','34.51x51','D','41', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
  ['A','235','6.37x37','D','383', ["hydrophobic", "van-der-waals"]],
  ['A','210','5.71x71','D','371', ["polar-sidechain-sidechain"]],
  ['A','234','6.36x36','D','382', ["polar-sidechain-backbone"]],
  ['A','38','12.51x51','B','52', ["polar-sidechain-sidechain", "van-der-waals"]],
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
  "rec_aa",
  "rec_gn",
  "rec_pos",
  "sig_aa",
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
// the next two datasets could also be generated by _.uniqBy lodash
// data that has unique receptor__generic_name entries
let data_t_rec = data_t.filter(
  (thing, index, self) =>
    index === self.findIndex(t => t.rec_gn === thing.rec_gn)
);
// data that has unique sigprot__generic_name entries
let data_t_sig = data_t.filter(
  (thing, index, self) =>
    index === self.findIndex(t => t.sig_gn === thing.sig_gn)
);
// unique interaction types withou undefined
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
      (h + margin.top + margin.bottom)
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


// * SETTING THE COLOR SCALE
let colScale = d3
  .scaleOrdinal()
  .domain(int_ty)
  .range(d3.schemeDark2);

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

// * APPENDING AMINOACID SEQUENCE [RECEPTOR]
svg
  .append("g")
  .attr("id", "recAA")
  .attr("transform", "translate(" + -xScale.step() / 2 + "," + h + ")")
  .selectAll("text")
  .data(data_t_rec)
  .enter()
  .append("text")
  .attr("x", function(d: any) {
    return xScale(d.rec_gn);
  })
  .attr("y", function(d: any, i) {
    return d * 20;
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
d3.select("g#recAA")
  .append("rect")
  .style("stroke", "black")
  .style("fill", "none")
  .attr("x", yScale.step() / 2)
  .attr("y", 60)
  .attr("width", w - xScale.step())
  .attr("height", seq_rect_h);

d3.select("g#sigAA")
  .append("rect")
  .style("stroke", "black")
  .style("fill", "none")
  .attr("x", -seq_rect_h / 2)
  .attr("y", yScale.step() / 2)
  .attr("width", seq_rect_h)
  .attr("height", h - yScale.step());

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

// TODO: ADD THE PDB ID AS LABEL



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
