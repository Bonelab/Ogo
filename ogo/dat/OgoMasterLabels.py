# /------------------------------------------------------------------------------+
# | 22-AUG-2022                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+

from collections import OrderedDict

labels_hdr = \
    '################################################\n' + \
    '# ITK-SnAP Label Description File\n' + \
    '# File format: \n' + \
    '# IDX   -R-  -G-  -B-  -A--  VIS MSH  LABEL\n' + \
    '# Fields: \n' + \
    '#    IDX:   Zero-based index \n' + \
    '#    -R-:   Red color component (0..255)\n' + \
    '#    -G-:   Green color component (0..255)\n' + \
    '#    -B-:   Blue color component (0..255)\n' + \
    '#    -A-:   Label transparency (0.00 .. 1.00)\n' + \
    '#    VIS:   Label visibility (0 or 1)\n' + \
    '#    IDX:   Label mesh visibility (0 or 1)\n' + \
    '#  LABEL:   Label description \n' + \
    '#\n' + \
    '# Label Ranges\n' + \
    '#               Bones:   1 to 80\n' + \
    '#       Miscellaneous:  81 to 90\n' + \
    '#           Int Calib:  91 to 110\n' + \
    '# Mindways Model 3 CT: 111 to 120\n' + \
    '#            B-MAS200: 121 to 130\n' + \
    '# Mindways Model 3 QA: 131 to 140\n' + \
    '#       QRM-BDC 3-rod: 141 to 150\n' + \
    '#       QRM-BDC 6-rod: 151 to 160\n' + \
    '# Img Ana QCT-3D Plus: 161 to 170\n' + \
    '#\n' + \
    '# Hint:\n' + \
    '# If you need to add new labels, add it into the \n' + \
    '# correct range and assign it a unique color. To\n' + \
    '# see possible colors, ogoPrintLabels --colors.\n' + \
    '# \n' + \
    '# August 3, 2022\n' + \
    '################################################\n' + \
    '# IND    R    G    B    A  VIS  IDX  LABEL\n'  + \
    '################################################'

labels_dict = {
    0: {'RGB': [    0,    0,    0], 'A': 0.0, 'VIS': 0, 'IDX': 0, 'LABEL': 'Clear Label'},
    1: {'RGB': [  255,    0,    0], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Femur Right'}, # femur_right
    2: {'RGB': [    0,  255,    0], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Femur Left'}, # femur_left
    3: {'RGB': [    0,    0,  255], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Pelvis Right'}, # hip_right
    4: {'RGB': [  255,  255,    0], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Pelvis Left'}, # hip_left
    5: {'RGB': [    0,  255,  255], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Sacrum'}, # sacrum
    6: {'RGB': [  255,    0,  255], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'L5'}, # vertebrae_L5
    7: {'RGB': [  255,  239,  213], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'L4'}, # vertebrae_L4
    8: {'RGB': [    0,    0,  205], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'L3'}, # vertebrae_L3
    9: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'L2'}, # vertebrae_L2
   10: {'RGB': [  210,  180,  140], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'L1'}, # vertebrae_L1
   11: {'RGB': [    0,  139,  139], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'T12'}, # vertebrae_T12
   12: {'RGB': [  224,  255,  255], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'T11'}, # vertebrae_T11
   13: {'RGB': [   72,  209,  204], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'T10'}, # vertebrae_T10
   14: {'RGB': [  176,  224,  230], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'T9'}, # vertebrae_T9
   15: {'RGB': [  100,  149,  237], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'T8'}, # vertebrae_T8
   16: {'RGB': [  173,  216,  230], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Radius Right'},
   17: {'RGB': [    0,    0,  139], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Radius Left'},
   18: {'RGB': [  139,   69,   19], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Ulna Right'},
   19: {'RGB': [   72,   61,  139], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Ulna Left'},
   20: {'RGB': [  100,  149,  237], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Humerus Right'}, # humerus_right
   21: {'RGB': [  153,   50,  204], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Humerus Left'}, # humerus_left
   22: {'RGB': [  147,  112,  219], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Tibia Right'},
   23: {'RGB': [  153,    0,  204], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Tibia Left'},
   24: {'RGB': [  216,  191,  216], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Fibula Right'},
   25: {'RGB': [  255,    0,  255], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Fibula Left'},
   26: {'RGB': [  219,  112,  147], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Patella Right'},
   27: {'RGB': [  255,  182,  193], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Patella Left'},
   28: {'RGB': [  245,  245,  220], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Calcaneus Right'},
   29: {'RGB': [  154,  205,   50], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Calcaneus Left'},
   30: {'RGB': [  224,  255,  255], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Navicular Right'},
   31: {'RGB': [    0,  206,  209], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Navicular Left'},
   32: {'RGB': [  138,   43,  226], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Talus Right'},
   33: {'RGB': [  176,  224,  230], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Talus Left'},
   34: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'L6'},
   35: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'T7'}, # vertebrae_T7
   36: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'T6'}, # vertebrae_T6
   37: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'T5'}, # vertebrae_T5
   38: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'T4'}, # vertebrae_T4
   39: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'T3'}, # vertebrae_T3
   40: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'T2'}, # vertebrae_T2
   41: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'T1'}, # vertebrae_T1
   42: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'C7'}, # vertebrae_C7
   43: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'C6'}, # vertebrae_C6
   44: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'C5'}, # vertebrae_C5
   45: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'C4'}, # vertebrae_C4
   46: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'C3'}, # vertebrae_C3
   47: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'C2'}, # vertebrae_C2
   48: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'C1'}, # vertebrae_C1
   49: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rib Left 1'}, # rib_left_1
   50: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rib Left 2'}, # rib_left_2
   51: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rib Left 3'}, # rib_left_3
   52: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rib Left 4'}, # rib_left_4
   53: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rib Left 5'}, # rib_left_5
   54: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rib Left 6'}, # rib_left_6
   55: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rib Left 7'}, # rib_left_7
   56: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rib Left 8'}, # rib_left_8
   57: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rib Left 9'}, # rib_left_9
   58: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rib Left 10'}, # rib_left_10
   59: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rib Left 11'}, # rib_left_11
   60: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rib Left 12'}, # rib_left_12
   61: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rib Right 1'}, # rib_right_1
   62: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rib Right 2'}, # rib_right_2
   63: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rib Right 3'}, # rib_right_3
   64: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rib Right 4'}, # rib_right_4
   65: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rib Right 5'}, # rib_right_5
   66: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rib Right 6'}, # rib_right_6
   67: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rib Right 7'}, # rib_right_7
   68: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rib Right 8'}, # rib_right_8
   69: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rib Right 9'}, # rib_right_9
   70: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rib Right 10'}, # rib_right_10
   71: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rib Right 11'}, # rib_right_11
   72: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rib Right 12'}, # rib_right_12
   73: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Scapula Left'}, # scapula_left
   74: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Scapula Right'}, # scapula_right
   75: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Clavicula Left'}, # clavicula_left
   76: {'RGB': [  152,  251,  152], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Clavicula Right'}, # clavicula_right
   81: {'RGB': [  192,  192,  192], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Not bone'},
   91: {'RGB': [  128,    0,    0], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Adipose'},
   92: {'RGB': [  128,  128,    0], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Air'},
   93: {'RGB': [    0,  128,    0], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Blood'},
   94: {'RGB': [  128,    0,  128], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Cortical Bone'},
   95: {'RGB': [    0,  128,  128], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Skeletal Muscle'},
  111: {'RGB': [  255,    0,    0], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rod A - Mindways Model 3 CT - Low K2HPO4'},
  112: {'RGB': [    0,  255,    0], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rod B - Mindways Model 3 CT'},
  113: {'RGB': [    0,    0,  255], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rod C - Mindways Model 3 CT'},
  114: {'RGB': [  255,    0,  255], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rod D - Mindways Model 3 CT'},
  115: {'RGB': [    0,  255,  255], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rod E - Mindways Model 3 CT - High K2HPO4'},
  121: {'RGB': [  255,    0,    0], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rod A - B-MAS200 -   0 mg/cc CHA'},
  122: {'RGB': [    0,  255,    0], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rod B - B-MAS200 -  50 mg/cc CHA'},
  123: {'RGB': [    0,    0,  255], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rod C - B-MAS200 - 100 mg/cc CHA'},
  124: {'RGB': [  255,  255,    0], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rod D - B-MAS200 - 150 mg/cc CHA'},
  125: {'RGB': [    0,  255,  255], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rod E - B-MAS200 - 200 mg/cc CHA'},
  131: {'RGB': [  219,  112,  147], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rod A - Mindways Model 3 QA - Low K2HPO4'},
  132: {'RGB': [  255,  182,  193], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rod B - Mindways Model 3 QA'},
  133: {'RGB': [  245,  245,  220], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rod C - Mindways Model 3 QA'},
  134: {'RGB': [  154,  205,   50], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rod D - Mindways Model 3 QA - High K2HPO4'},
  141: {'RGB': [  100,  149,  237], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rod A - QRM-BDC 3-rod -   0 mg/cc CHA'},
  142: {'RGB': [  173,  216,  230], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rod B - QRM-BDC 3-rod - 400 mg/cc CHA'},
  143: {'RGB': [    0,    0,  139], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rod C - QRM-BDC 3-rod - 800 mg/cc CHA'},
  151: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rod A - QRM-BDC 6-rod -   0 mg/cc CHA'},
  152: {'RGB': [  210,  180,  140], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rod B - QRM-BDC 6-rod - 100 mg/cc CHA'},
  153: {'RGB': [    0,  139,  139], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rod C - QRM-BDC 6-rod - 200 mg/cc CHA'},
  154: {'RGB': [  224,  255,  255], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rod D - QRM-BDC 6-rod - 400 mg/cc CHA'},
  155: {'RGB': [   72,  209,  204], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rod E - QRM-BDC 6-rod - 600 mg/cc CHA'},
  156: {'RGB': [  176,  224,  230], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rod F - QRM-BDC 6-rod - 800 mg/cc CHA'},
  161: {'RGB': [    0,    0,  189], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rod A - Image Analysis QCT-3D Plus -   0 mg/cc CHA'},
  162: {'RGB': [  139,   69,   19], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rod B - Image Analysis QCT-3D Plus -  75 mg/cc CHA'},
  163: {'RGB': [   72,   61,  139], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Rod C - Image Analysis QCT-3D Plus - 150 mg/cc CHA'},
  200: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Spleen'}, # spleen # organs
  201: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Kidney Right'}, # kidney_right
  202: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Kidney Left'}, # kidney_left
  203: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Gallbladder'}, # gallbladder
  204: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Liver'}, # liver
  205: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Stomach'}, # stomach
  206: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Pancreas'}, # pancreas
  207: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Adrenal Gland Right'}, # adrenal_gland_right
  208: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Adrenal Gland Left'}, # adrenal_gland_left
  209: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Esophagus'}, # esophagus
  210: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Trachea'}, # trachea
  211: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Brain'}, # brain
  212: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Small Bowel'}, # small_bowel
  213: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Duodenum'}, # duodenum
  214: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Colon'}, # colon
  215: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Face'}, # face
  216: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Urinary Bladder'}, # urinary_bladder
  217: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Lung Upper Lobe Left'}, # lung_upper_lobe_left
  218: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Lung Lower Lobe Left'}, # lung_lower_lobe_left
  219: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Lung Upper Lobe Right'}, # lung_upper_lobe_right
  220: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Lung Middle Lobe Right'}, # lung_middle_lobe_right
  221: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Lung Lower Lobe Right'}, # lung_lower_lobe_right
  230: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Heart Myocardium'}, # heart_myocardium # heart
  231: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Heart Atrium Left'}, # heart_atrium_left
  232: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Heart Ventricle Left'}, # heart_ventricle_left
  233: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Heart Atrium Right'}, # heart_atrium_right
  234: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Heart Ventricle Right'}, # heart_ventricle_right
  240: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Aorta'}, # aorta # vasculature
  241: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Inferior Vena Cava'}, # inferior_vena_cava
  242: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Portal Vein and Splenic Vein'}, # portal_vein_and_splenic_vein
  243: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Pulmonary Artery'}, # pulmonary_artery
  244: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Iliac Artery Left'}, # iliac_artery_left
  245: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Iliac Artery Right'}, # iliac_artery_right
  246: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Iliac Vena Left'}, # iliac_vena_left
  247: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Iliac Vena Right'}, # iliac_vena_right
  250: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Gluteus Maximus Left # muscle'}, # gluteus_maximus_left # muscle
  251: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Gluteus Maximus Right'}, # gluteus_maximus_right
  252: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Gluteus Medius Left'}, # gluteus_medius_left
  253: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Gluteus Medius Right'}, # gluteus_medius_right
  254: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Gluteus Minimus Left'}, # gluteus_minimus_left
  255: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Gluteus Minimus Right'}, # gluteus_minimus_right
  256: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Autochthon Left'}, # autochthon_left
  257: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Autochthon Right'}, # autochthon_right
  258: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Iliopsoas Left'}, # iliopsoas_left
  259: {'RGB': [  205,  133,   63], 'A': 1.0, 'VIS': 1, 'IDX': 1, 'LABEL': 'Iliopsoas right'}  # iliopsoas_right
}

#print(labels_hdr)
#print("\n".join("{:5d} {:4d} {:4d} {:4d}  {:.1f} {:4d} {:4d}  \"{}\"".format(k,v['RGB'][0],v['RGB'][1],v['RGB'][2],v['A'],v['VIS'],v['IDX'],v['LABEL']) \
#                                                   for k, v in labels_dict.items()))


colors128_hdr = \
    '################################################\n' + \
    '# IND    R     G     B   DESC  - 128 colors\n' + \
    '################################################'

colors128_dict = {
       1: {'RGB':[128,   0,   0], 'DESC': 'maroon'},
       2: {'RGB':[139,   0,   0], 'DESC': 'dark red'},
       3: {'RGB':[165,  42,  42], 'DESC': 'brown'},
       4: {'RGB':[178,  34,  34], 'DESC': 'firebrick'},
       5: {'RGB':[220,  20,  60], 'DESC': 'crimson'},
       6: {'RGB':[255,   0,   0], 'DESC': 'red'},
       7: {'RGB':[255,  99,  71], 'DESC': 'tomato'},
       8: {'RGB':[255, 127,  80], 'DESC': 'coral'},
       9: {'RGB':[205,  92,  92], 'DESC': 'indian red'},
      10: {'RGB':[240, 128, 128], 'DESC': 'light coral'},
      11: {'RGB':[233, 150, 122], 'DESC': 'dark salmon'},
      12: {'RGB':[250, 128, 114], 'DESC': 'salmon'},
      13: {'RGB':[255, 160, 122], 'DESC': 'light salmon'},
      14: {'RGB':[255,  69,   0], 'DESC': 'orange red'},
      15: {'RGB':[255, 140,   0], 'DESC': 'dark orange'},
      16: {'RGB':[255, 165,   0], 'DESC': 'orange'},
      17: {'RGB':[255, 215,   0], 'DESC': 'gold'},
      18: {'RGB':[184, 134,  11], 'DESC': 'dark golden rod'},
      19: {'RGB':[218, 165,  32], 'DESC': 'golden rod'},
      20: {'RGB':[238, 232, 170], 'DESC': 'pale golden rod'},
      21: {'RGB':[189, 183, 107], 'DESC': 'dark khaki'},
      22: {'RGB':[240, 230, 140], 'DESC': 'khaki'},
      23: {'RGB':[128, 128,   0], 'DESC': 'olive'},
      24: {'RGB':[255, 255,   0], 'DESC': 'yellow'},
      25: {'RGB':[154, 205,  50], 'DESC': 'yellow green'},
      26: {'RGB':[ 85, 107,  47], 'DESC': 'dark olive green'},
      27: {'RGB':[107, 142,  35], 'DESC': 'olive drab'},
      28: {'RGB':[127, 255,   0], 'DESC': 'chartreuse'},
      29: {'RGB':[173, 255,  47], 'DESC': 'green yellow'},
      30: {'RGB':[  0, 100,   0], 'DESC': 'dark green'},
      31: {'RGB':[  0, 128,   0], 'DESC': 'green'},
      32: {'RGB':[  0, 255,   0], 'DESC': 'lime'},
      33: {'RGB':[ 50, 205,  50], 'DESC': 'lime green'},
      34: {'RGB':[144, 238, 144], 'DESC': 'light green'},
      35: {'RGB':[152, 251, 152], 'DESC': 'pale green'},
      36: {'RGB':[143, 188, 143], 'DESC': 'dark sea green'},
      37: {'RGB':[  0, 250, 154], 'DESC': 'medium spring green'},
      38: {'RGB':[  0, 255, 127], 'DESC': 'spring green'},
      39: {'RGB':[ 46, 139,  87], 'DESC': 'sea green'},
      40: {'RGB':[102, 205, 170], 'DESC': 'medium aqua marine'},
      41: {'RGB':[ 60, 179, 113], 'DESC': 'medium sea green'},
      42: {'RGB':[ 32, 178, 170], 'DESC': 'light sea green'},
      43: {'RGB':[ 47,  79,  79], 'DESC': 'dark slate gray'},
      44: {'RGB':[  0, 128, 128], 'DESC': 'teal'},
      45: {'RGB':[  0, 139, 139], 'DESC': 'dark cyan'},
      46: {'RGB':[  0, 255, 255], 'DESC': 'aqua'},
      47: {'RGB':[  0, 255, 255], 'DESC': 'cyan'},
      48: {'RGB':[224, 255, 255], 'DESC': 'light cyan'},
      49: {'RGB':[  0, 206, 209], 'DESC': 'dark turquoise'},
      50: {'RGB':[ 64, 224, 208], 'DESC': 'turquoise'},
      51: {'RGB':[ 72, 209, 204], 'DESC': 'medium turquoise'},
      52: {'RGB':[175, 238, 238], 'DESC': 'pale turquoise'},
      53: {'RGB':[127, 255, 212], 'DESC': 'aqua marine'},
      54: {'RGB':[176, 224, 230], 'DESC': 'powder blue'},
      55: {'RGB':[ 95, 158, 160], 'DESC': 'cadet blue'},
      56: {'RGB':[ 70, 130, 180], 'DESC': 'steel blue'},
      57: {'RGB':[100, 149, 237], 'DESC': 'corn flower blue'},
      58: {'RGB':[  0, 191, 255], 'DESC': 'deep sky blue'},
      59: {'RGB':[ 30, 144, 255], 'DESC': 'dodger blue'},
      60: {'RGB':[173, 216, 230], 'DESC': 'light blue'},
      61: {'RGB':[135, 206, 235], 'DESC': 'sky blue'},
      62: {'RGB':[135, 206, 250], 'DESC': 'light sky blue'},
      63: {'RGB':[  0,   0, 139], 'DESC': 'dark blue'},
      64: {'RGB':[  0,   0, 205], 'DESC': 'medium blue'},
      65: {'RGB':[  0,   0, 255], 'DESC': 'blue'},
      66: {'RGB':[ 65, 105, 225], 'DESC': 'royal blue'},
      67: {'RGB':[138,  43, 226], 'DESC': 'blue violet'},
      68: {'RGB':[ 75,   0, 130], 'DESC': 'indigo'},
      69: {'RGB':[ 72,  61, 139], 'DESC': 'dark slate blue'},
      70: {'RGB':[106,  90, 205], 'DESC': 'slate blue'},
      71: {'RGB':[123, 104, 238], 'DESC': 'medium slate blue'},
      72: {'RGB':[147, 112, 219], 'DESC': 'medium purple'},
      73: {'RGB':[139,   0, 139], 'DESC': 'dark magenta'},
      74: {'RGB':[148,   0, 211], 'DESC': 'dark violet'},
      75: {'RGB':[153,  50, 204], 'DESC': 'dark orchid'},
      76: {'RGB':[186,  85, 211], 'DESC': 'medium orchid'},
      77: {'RGB':[128,   0, 128], 'DESC': 'purple'},
      78: {'RGB':[216, 191, 216], 'DESC': 'thistle'},
      79: {'RGB':[221, 160, 221], 'DESC': 'plum'},
      80: {'RGB':[238, 130, 238], 'DESC': 'violet'},
      81: {'RGB':[255,   0, 255], 'DESC': 'magenta'},
      82: {'RGB':[218, 112, 214], 'DESC': 'orchid'},
      83: {'RGB':[199,  21, 133], 'DESC': 'medium violet red'},
      84: {'RGB':[219, 112, 147], 'DESC': 'pale violet red'},
      85: {'RGB':[255,  20, 147], 'DESC': 'deep pink'},
      86: {'RGB':[255, 105, 180], 'DESC': 'hot pink'},
      87: {'RGB':[255, 182, 193], 'DESC': 'light pink'},
      88: {'RGB':[255, 192, 203], 'DESC': 'pink'},
      89: {'RGB':[250, 235, 215], 'DESC': 'antique white'},
      90: {'RGB':[245, 245, 220], 'DESC': 'beige'},
      91: {'RGB':[255, 228, 196], 'DESC': 'bisque'},
      92: {'RGB':[255, 235, 205], 'DESC': 'blanched almond'},
      93: {'RGB':[245, 222, 179], 'DESC': 'wheat'},
      94: {'RGB':[255, 248, 220], 'DESC': 'corn silk'},
      95: {'RGB':[255, 250, 205], 'DESC': 'lemon chiffon'},
      96: {'RGB':[250, 250, 210], 'DESC': 'light golden rod yellow'},
      97: {'RGB':[255, 255, 224], 'DESC': 'light yellow'},
      98: {'RGB':[139,  69,  19], 'DESC': 'saddle brown'},
      99: {'RGB':[160,  82,  45], 'DESC': 'sienna'},
     100: {'RGB':[210, 105,  30], 'DESC': 'chocolate'},
     101: {'RGB':[205, 133,  63], 'DESC': 'peru'},
     102: {'RGB':[244, 164,  96], 'DESC': 'sandy brown'},
     103: {'RGB':[222, 184, 135], 'DESC': 'burly wood'},
     104: {'RGB':[210, 180, 140], 'DESC': 'tan'},
     105: {'RGB':[188, 143, 143], 'DESC': 'rosy brown'},
     106: {'RGB':[255, 228, 225], 'DESC': 'misty rose'},
     107: {'RGB':[255, 240, 245], 'DESC': 'lavender blush'},
     108: {'RGB':[250, 240, 230], 'DESC': 'linen'},
     109: {'RGB':[253, 245, 230], 'DESC': 'old lace'},
     110: {'RGB':[255, 239, 213], 'DESC': 'papaya whip'},
     111: {'RGB':[255, 245, 238], 'DESC': 'sea shell'},
     112: {'RGB':[245, 255, 250], 'DESC': 'mint cream'},
     113: {'RGB':[112, 128, 144], 'DESC': 'slate gray'},
     114: {'RGB':[119, 136, 153], 'DESC': 'light slate gray'},
     115: {'RGB':[176, 196, 222], 'DESC': 'light steel blue'},
     116: {'RGB':[230, 230, 250], 'DESC': 'lavender'},
     117: {'RGB':[255, 250, 240], 'DESC': 'floral white'},
     118: {'RGB':[240, 248, 255], 'DESC': 'alice blue'},
     119: {'RGB':[248, 248, 255], 'DESC': 'ghost white'},
     120: {'RGB':[  0,   0,   0], 'DESC': 'black'},
     121: {'RGB':[105, 105, 105], 'DESC': 'dim gray'},
     122: {'RGB':[128, 128, 128], 'DESC': 'gray'},
     123: {'RGB':[169, 169, 169], 'DESC': 'dark gray'},
     124: {'RGB':[192, 192, 192], 'DESC': 'silver'},
     125: {'RGB':[211, 211, 211], 'DESC': 'light gray'},
     126: {'RGB':[220, 220, 220], 'DESC': 'gainsboro'},
     127: {'RGB':[245, 245, 245], 'DESC': 'white smoke'},
     128: {'RGB':[255, 255, 255], 'DESC': 'white'}
}

#print(colors128_hdr)
#print("\n".join("{:5d} {:4d}  {:4d}  {:4d}   {:s}".format(k,v['RGB'][0],v['RGB'][1],v['RGB'][2],v['DESC']) \
#                                                   for k, v in colors128_dict.items()))

colors16_hdr = \
    '################################################\n' + \
    '# IND    R     G     B   DESC  - 16 colors\n' + \
    '################################################'

colors16_dict = {
     1: {'RGB':[  0,   0,   0], 'DESC': 'Black'},
     2: {'RGB':[255, 255, 255], 'DESC': 'White'},
     3: {'RGB':[255,   0,   0], 'DESC': 'Red'},
     4: {'RGB':[  0, 255,   0], 'DESC': 'Lime'},
     5: {'RGB':[  0,   0, 255], 'DESC': 'Blue'},
     6: {'RGB':[255, 255,   0], 'DESC': 'Yellow'},
     7: {'RGB':[  0, 255, 255], 'DESC': 'Cyan'},
     8: {'RGB':[255,   0, 255], 'DESC': 'Magenta'},
     9: {'RGB':[192, 192, 192], 'DESC': 'Silver'},
    10: {'RGB':[128, 128, 128], 'DESC': 'Gray'},
    11: {'RGB':[128,   0,   0], 'DESC': 'Maroon'},
    12: {'RGB':[128, 128,   0], 'DESC': 'Olive'},
    13: {'RGB':[  0, 128,   0], 'DESC': 'Green'},
    14: {'RGB':[128,   0, 128], 'DESC': 'Purple'},
    15: {'RGB':[  0, 128, 128], 'DESC': 'Teal'},
    16: {'RGB':[  0,   0, 128], 'DESC': 'Navy'}
    }

#print(colors16_hdr)
#print("\n".join("{:5d} {:4d}  {:4d}  {:4d}   {:s}".format(k,v['RGB'][0],v['RGB'][1],v['RGB'][2],v['DESC']) \
#                                                   for k, v in colors16_dict.items()))
