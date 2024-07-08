from astropy import constants as const
from astropy import units as u
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from math import pi
import seaborn as sns
import scipy.stats as stats
import MDM_targets as md
import pandas as pd

#Declare constants 
EARTH_MASS = const.M_earth # Value = 5.97216787e+24 kg
SUN_MASS = const.M_sun # Value = 1.98840987e+30 kg


#Declare arrays
q_array = np.random.uniform(0,1, 200)
inclination_angle_array = md.sin_distn()
theta_array = np.random.uniform(0, 2 * pi, 100) * u.radian

rv_times = md.rv_times()



'''
    In this part of the code, I wrote the measured (known) radial velocity curves
    for each given KIC target. I used the rc_times table and wrote the measured
    velocities in the form of an array, which we will use after in our 
    Chi-square statistics test
'''
KIC_1570924_time_of_observation = np.array([7904.8034, 7906.8315, 7907.81672, 7908.79553, 7908.93257, 7909.78197, 
                                            7909.92973, 7910.88811, 7911.81382, 7912.74221, 7913.77538, 7913.93496,
                                            7914.74206, 7914.87532, 7915.74782, 7915.88332, 7916.75207, 7916.89091,
                                            7917.80579, 7917.92246])
KIC_1570924_measured_RVs = np.array([-6.01430605, -2.79444638, -14.74919465, -7.28997981, -13.76927981, 4.03473816,
                                     -6.98376184, -9.17476184, -11.5024829, -5.30175823, 5.976315, -10.399385, -8.75004155,
                                     -9.37544155, -10.45104242,  -8.79274242, -10.12132926, -9.03902926, -17.65446889, -8.37066889])
'''
This file contains a tuple consisting of these instructions as the first \nelement, and a dictionary mapping
 the KIC name of the targets to a 2-tuple\nconaining the HJD times of the observation and the measured RVs.
'''
KIC_3219623_time_of_observation = np.array([7906.74891, 7906.88693, 7907.74473, 7907.8944, 7908.8798, 7909.80402, 7909.89322, 7910.90191,
                                            7911.75498, 7912.89384, 7913.8424, 7913.90054, 7914.85325, 7915.80738, 7915.8654, 7916.73832,
                                            7916.86271, 7917.75048, 7917.85351])
KIC_3219623_measured_RVs = np.array([-15.87344638, -17.86364638, -23.58099465, -29.28889465, -20.81007981, -22.86716184, -22.59036184,
                                     -25.59416184, -33.5015829, -29.04095823, -14.142685, -8.786385, -13.62394155, -28.07554242, -25.85554242,
                                     -25.05692926, -21.77162926, -26.15586889, -20.69166889])
                                            
KIC_3248885_time_of_observation = np.array([7904.93806, 7906.78638, 7906.89583, 7906.92054, 7907.84185, 7908.81424, 7908.94525,
                                            7909.74815, 7909.89836, 7910.86986, 7911.804, 7912.78572, 7913.88223, 7914.78595,
                                            7914.91465, 7915.7823, 7916.81653, 7917.823])
KIC_3248885_measured_RVs = np.array([-12.22140605, -9.73524638, -14.33824638, -20.28854638, -53.74639465, -52.42627981, -42.58767981,
                                     -34.74106184, -14.28746184, -12.77616184, -50.5761829, -62.79715823, -27.935585, -15.91764155,
                                     -9.57944155, -23.60114242, -42.07872926, -46.38186889])

KIC_3539632_time_of_observation = np.array([7904.85205, 7906.74374, 7907.81084, 7908.86453, 7909.81217, 7910.85126,
                                            7911.79575, 7912.90926, 7913.7851, 7913.94311, 7914.7203, 7914.86145, 
                                            7915.72984, 7916.7192, 7916.87223, 7917.71446, 7917.86266])
KIC_3539632_measured_RVs = np.array([-20.31600605, -11.72344638, -26.76719465, -21.42027981, -15.08956184, -20.94156184,
                                     -21.8344829, -12.10005823, -23.335585, -24.812085, -21.88834155, -17.28414155,
                                     -21.51404242, -19.32972926, -15.64572926, -19.98696889, -23.03696889])

KIC_3540728_time_of_observation = np.array([7904.90137, 7906.87733, 7907.79253, 7907.9322, 7908.79114, 7908.92829, 7909.83896,
                                            7910.81033, 7911.82078, 7912.79304, 7913.79768, 7914.76933, 7914.90127, 7915.73822,
                                            7915.87363, 7916.78197, 7916.92125, 7917.81911, 7917.94817])
KIC_3540728_measured_RVs = np.array([-12.42480605, -23.32224638, -19.87099465, -25.45679465, -13.14327981, -28.71907981, -21.84416184,
                                     -25.01726184, -25.5392829, -12.46075823, -18.877185, -15.31004155, -19.28754155, -23.60094242,
                                     -24.00164242, -15.02022926, -25.13312926, -20.30196889, -21.42556889])

KIC_4036736_time_of_observation = np.array([7904.92132, 7906.80539, 7907.86827, 7908.75913, 7908.89712, 7909.76613, 7909.91541,
                                            7910.85774, 7911.90764, 7912.92477, 7913.89507, 7914.849, 7915.80306, 7916.84835,
                                            7917.74252, 7917.87489])
KIC_4036736_measured_RVs = np.array([-9.29480605, -16.47574638, -28.15809465, -17.80147981, -19.79557981, -14.75566184, -29.25616184,
                                     -23.15046184, -42.9573829, -33.20495823, -12.023685, -21.18474155, -23.45484242, -21.91502926,
                                     -23.71666889, -24.16346889])

KIC_4249702_time_of_observation = np.array([7904.81764, 7906.8647, 7906.91604, 7907.79696, 7907.9369, 7908.78247,
                                            7908.91986, 7909.86145, 7910.82122, 7911.7844, 7912.86998, 7913.83345, 
                                            7914.83925, 7915.76014, 7915.89556, 7916.83044, 7917.81514, 7917.9399])
KIC_4249702_measured_RVs = np.array([-17.88240605, -22.78654638, -4.12864638, -24.86959465, -21.41739465, -23.07587981,
                                     -29.44687981, -12.23746184, -28.22876184, -23.5693829, -18.70285823, -22.445585,
                                     -17.28464155, -27.32264242, -28.12754242, -22.06542926, -18.08606889, -25.37006889])

KIC_4454890_time_of_observation = np.array([7904.884, 7906.7669, 7906.90207, 7907.7778, 7907.91858, 7908.81855, 7909.78596,
                                            7909.93382, 7910.87547, 7911.9227, 7912.75946, 7913.81974, 7914.72951, 7914.87097,
                                            7915.83919, 7915.93134, 7916.7901, 7916.92958, 7917.7636, 7917.88425])
KIC_4454890_measured_RVs = np.array([-16.41120605, -14.40754638, -26.97814638, -25.23209465, -30.02249465, -21.19097981,
                                     -13.98346184, -23.49896184, -17.40156184, -29.3014829, -14.77405823, -13.735685,
                                     -20.76644155, -25.42434155, -33.79684242, -23.65304242, -25.42952926, -26.69092926,
                                     -17.09226889, -20.12816889])

KIC_4480434_time_of_observation = np.array([7904.92696, 7906.88189, 7906.9112, 7907.8582, 7908.84131, 7909.84723, 7910.82645,
                                            7911.85814, 7912.84291, 7913.82901,7914.79956, 7914.92737, 7915.77767, 7916.764, 
                                            7916.85732, 7916.90276, 7917.72403])
KIC_4480434_measured_RVs = np.array([-29.18680605, -24.91314638, -19.59454638, -32.55709465, -33.93467981, -30.07406184, -19.20516184,
                                     -49.6316829, -27.76265823, -25.504585, -16.93914155, -25.86194155, -35.66924242, -38.95552926,
                                     -32.58922926, -27.72612926, -18.80686889])

KIC_5213142_time_of_observation = np.array([7904.84317, 7906.81054, 7907.80564, 7908.77768, 7908.91511, 7909.8802, 7910.9096,
                                            7911.84796, 7912.79726, 7913.87337, 7914.79521, 7914.92282, 7915.76919, 7915.9048,
                                            7916.7557, 7916.89437, 7917.73487, 7917.87018])
KIC_5213142_measured_RVs = np.array([6.86469395, 2.56215362, -12.99579465, -0.86957981, -3.90687981, -1.30566184, -2.09396184,
                                     -15.3734829, -9.62675823, 10.111515, 1.58475845, -1.34674155, -12.60304242, -11.62614242,
                                     -4.42222926, -1.24402926, -9.08116889, 15.47733111])

KIC_5553362_time_of_observation = np.array([7904.83654, 7906.79506, 7906.93986, 7907.76058, 7907.90325, 7908.74989, 
                                            7908.88779, 7909.85255, 7910.89165, 7911.89752, 7912.77561, 7913.88618,
                                            7914.82161, 7914.94996, 7915.8346, 7915.9267, 7916.72857, 7916.88027, 7917.84085])
KIC_5553362_measured_RVs = np.array([-2.42270605, -11.23634638, -13.62114638, -24.37919465, -26.05539465, -30.60757981,
                                     -14.52377981, -14.61346184, -12.52606184, -19.8595829, -23.07085823, -16.523685,
                                     -17.02544155, -16.24854155, -23.81354242, -16.36494242, -33.33472926, -24.18012926, -25.79256889])

KIC_5609753_time_of_observation = np.array([7904.85842, 7904.9067, 7906.75807, 7906.78164, 7906.85506, 7906.93088, 7907.73989,
                                            7907.76736, 7907.83741, 7907.8896, 7907.92776, 7908.76744, 7908.80942, 7908.86854,
                                            7908.90708, 7909.75639, 7909.79009, 7909.8431, 7909.88477, 7909.9254, 7910.93649,
                                            7911.76624, 7911.83823, 7911.8799, 7911.93628, 7912.7375 , 7912.78051, 7912.87694,
                                            7912.94217, 7913.76644, 7913.81128, 7913.8517, 7913.91706, 7913.93879, 7914.72516,
                                            7914.75133, 7914.80822, 7914.86553, 7914.89702, 7914.94616, 7915.73397, 7915.7737,
                                            7915.83064, 7915.89133, 7915.93973, 7916.73389, 7916.76001, 7916.80341, 7916.85216,
                                            7916.89884, 7916.9395, 7917.7192, 7917.74625, 7917.78927, 7917.84927, 7917.89653, 7917.93572])
KIC_5609753_measured_RVs = np.array([88.35779395, 60.06299395, 57.05585362, 55.96955362, 53.60535362, 60.02045362, 70.73020535,
                                     55.13850535, 59.60150535, 56.70350535, 47.04270535, 28.40652019, 26.03952019, 35.73442019,
                                     19.20312019, 50.73633816, 54.15343816, 54.01643816, 50.43793816, 66.67833816, 66.46263816,
                                     17.5891171, 25.8704171, 16.5195171, 34.0407171, 54.37394177, 52.42094177, 62.73644177, 56.52204177,
                                     74.019215, 78.593615, 82.005615, 71.724215, 68.153415, 35.52315845, 35.86155845, 31.75875845,
                                     19.88545845, 27.75855845, 22.99795845, 37.55375758, 33.28545758, 44.53685758, 48.22515758,
                                     46.53185758, 66.01777074, 71.91737074, 68.11697074, 71.07547074, 62.90917074, 57.55557074,
                                     44.67213111, 40.63783111, 33.31863111, 34.71333111, 26.12433111, 29.00703111])

KIC_6425783_time_of_observation = np.array([7904.8887, 7906.92516, 7907.78766, 7907.94604, 7908.76322, 7908.90148, 7909.87632,
                                            7910.87966, 7911.88985, 7912.86191, 7913.89064, 7914.84447, 7915.84862, 7915.94809,
                                            7916.80828, 7917.81084, 7917.93171])
KIC_6425783_measured_RVs = np.array([-28.05890605, -34.18694638, -30.46419465, -35.81789465, -4.65527981, -29.79827981, -23.18536184,
                                     -25.09826184, -18.4113829, -27.86195823, -26.229685, -29.46924155, -27.65064242, -31.35614242,
                                     -24.12362926, -25.27676889, -28.49196889])

KIC_6780052_time_of_observation = np.array([7904.91181, 7906.77136, 7906.90633, 7907.85259, 7908.85996, 7909.82989, 7910.86128, 7911.775,
                                            7912.85001, 7913.80269, 7914.77304, 7914.90503, 7915.82195, 7915.91728, 7916.77751, 7916.91671,
                                            7917.77947, 7917.90026])
KIC_6780052_measured_RVs = np.array([-28.59170605, -21.09614638, -15.96734638, -26.54909465, -22.01137981, -18.97836184, -22.79616184,
                                     -32.0563829, -6.87685823, -13.454485, -13.83964155, -19.21594155, -27.63884242, -25.89664242, 
                                     -22.04352926, -27.88682926, -16.92156889, -24.64956889])

KIC_6844101_time_of_observation = np.array([7904.79316, 7904.79635, 7906.87265, 7907.78218, 7907.92313, 7908.82822,
                                            7909.86975, 7910.89758, 7911.78752, 7912.88333, 7913.81547, 7914.75558,
                                            7914.88745, 7915.75182, 7915.88724, 7916.77373, 7916.91291, 7917.79332, 7917.90986])
KIC_6844101_measured_RVs = np.array([-16.89540605, -15.55850605, -27.24284638, -26.05909465, -26.23149465, -19.53657981,
                                     -18.79466184, -20.09076184, -23.8857829 , -24.71235823, -20.338785  , -17.11314155,
                                     -26.14964155, -24.22194242, -21.44814242, -12.78592926, -31.32262926, -24.20626889, -33.07126889])

KIC_7421325_time_of_observation = np.array([7904.82371, 7906.84493, 7907.86347, 7908.84585, 7909.75236, 7909.90257, 7910.86565,
                                            7911.88449, 7912.88842, 7913.86826, 7914.82613, 7915.74338, 7915.87881, 7916.83854, 7917.84523])
KIC_7421325_measured_RVs = np.array([0.79679395, -19.90874638, -51.88679465, 7.47052019, 11.97933816, 2.18773816, 6.15463816, -8.4374829,
                                     -0.22235823, -0.909985, -6.98304155, -26.63044242, -38.78844242, -22.28242926, -2.53136889])

KIC_7919763_time_of_observation = np.array([7904.93187, 7906.94403, 7907.82039, 7908.80016, 7908.93607, 7909.75982, 7909.91002, 7910.91908,
                                            7911.91468, 7912.91562, 7913.86342, 7914.80322, 7914.93102, 7915.82563, 7915.92098, 7916.83369,
                                            7917.79676, 7917.91327])
KIC_7919763_measured_RVs = np.array([-20.68860605, -4.25804638, -31.42499465, -17.81817981, -18.07787981, -22.76526184, -14.96096184,
                                     -8.67786184, -23.7222829, -12.84095823, -10.611385, -20.13224155, -14.09324155, -33.20974242, -16.34704242,
                                     -22.84832926, -21.40366889, -20.79466889])

KIC_8442720_time_of_observation = np.array([7904.87483, 7906.85886, 7907.83015, 7908.85109, 7909.82535, 7910.92741, 7911.84334,
                                            7912.83551, 7913.78857, 7913.94654, 7914.74657, 7914.88263, 7915.84268, 7915.93496,
                                            7916.74716, 7916.886, 7917.76704, 7917.88757])
KIC_8442720_measured_RVs = np.array([-54.01940605, -24.71004638, -65.73309465, -24.01667981, -31.58826184, -59.02446184, -27.5864829,
                                     -39.41395823, -45.927185, -45.820085, -38.48294155, -19.36314155, -46.39864242, -46.19314242, -57.10062926,
                                     -48.62502926, -22.43456889, -20.40256889])

KIC_8651471_time_of_observation = np.array([7904.82794, 7906.81498, 7907.77245, 7907.9131, 7908.85507, 7909.79331, 7909.93705,
                                            7911.83017, 7912.81732, 7913.8242, 7914.82958, 7915.78573, 7916.81158, 7917.75745, 7917.87799])
KIC_8651471_measured_RVs = np.array([-4.88300605, -0.91854638, -22.13019465, 0.30780535, -5.06087981, -7.54076184, 3.96863816,
                                     -18.0607829, -6.76695823, 0.644615, -3.47934155, -11.63844242, -13.07772926, -9.90196889, -0.65956889])

KIC_9151271_time_of_observation = np.array([7907.74843, 7907.80106, 7907.8473, 7907.90902, 7907.94086, 7908.75409, 7908.78625, 7908.83181, 7908.8921,
                                            7908.92352, 7909.76946, 7909.82011, 7909.85677, 7909.90621, 7909.94549, 7910.93184, 7911.76135, 7911.79817,
                                            7911.85233, 7911.93117, 7912.75442, 7912.81158, 7912.93663, 7913.75499, 7913.79282, 7913.83746, 7913.87736,
                                            7913.93012, 7914.73309, 7914.78103, 7914.81705, 7914.87866, 7914.93652, 7915.72496, 7915.75542, 7915.81074,
                                            7915.8564, 7915.90875, 7915.94352, 7916.72349, 7916.78551, 7916.82559, 7916.87574, 7916.92488, 7916.94313,
                                            7917.72979, 7917.77523, 7917.82672, 7917.86605, 7917.92581, 7917.94357])
KIC_9151271_measured_RVs = np.array([-33.29739465, -38.80649465, -37.59259465, -40.18919465, -39.38629465, -45.13497981, -47.23507981, -37.86297981,
                                     -47.04837981, -37.97757981, -34.18986184, -41.82726184, -45.26106184, -39.48636184, -36.46236184, -46.04006184,
                                     -34.7138829, -31.0943829, -38.3196829, -40.0322829, -35.31925823, -41.72975823, -36.21785823, -37.460285,
                                     -30.326285, -32.464485, -34.987285, -32.807085, -35.21214155, -29.64634155, -38.48054155, -41.50434155,
                                     -37.44094155, -28.51274242, -36.39874242, -34.56724242, -30.40434242, -29.57184242, -35.58114242, -27.80672926,
                                     -27.16182926, -25.37782926, -32.42572926, -33.29202926, -35.66862926, -31.79356889, -25.80156889, -28.44836889,
                                     -26.38696889, -34.28156889, -40.04866889])

KIC_9653110_time_of_observation = np.array([7904.81111, 7906.79016, 7908.88311, 7909.79751, 7909.94122, 7910.88306,
                                            7911.82488, 7912.83031, 7913.8064, 7914.77675, 7914.9087 , 7915.71978, 
                                            7915.86862, 7916.76901, 7916.90804, 7917.78323, 7917.90386])
KIC_9653110_measured_RVs = np.array([-39.85270605, -18.99964638, -34.99267981, -26.15876184, -11.46366184, -39.87286184,
                                     -46.8367829, -20.28775823, -41.893785, -40.41324155, -36.41084155, -30.41554242,
                                     -37.92994242, -32.65932926, -35.79462926, -39.51856889, -48.50476889])

KIC_9655045_time_of_observation = np.array([7912.90315, 7913.74957, 7913.91211, 7914.7590])
KIC_9655045_measured_RVs = np.array([-6.64595823, -27.973585, 0.080715, -9.54344155])

KIC_9710336_time_of_observation = np.array([7904.94352, 7906.79936, 7906.94845, 7907.73473, 7907.8838, 7908.77307, 7908.91071,
                                            7909.86481, 7910.846, 7911.8734, 7912.72712, 7913.84566, 7914.81167, 7914.94109,
                                            7915.85211, 7916.79471, 7916.93474, 7917.80098, 7917.91747])
KIC_9710336_measured_RVs = np.array([-26.92980605, -1.15464638, -9.64574638, -18.52349465, -28.93639465, -29.75657981, -33.31617981,
                                     -33.27546184, -22.18866184, -13.3846829, -25.07025823, -20.528885, -15.96414155, -17.53464155, 
                                     3.05135758, -17.15142926, -29.47732926, -32.36416889, -26.07046889])

KIC_9964938_time_of_observation = np.array([7906.75231, 7906.89041, 7907.75248, 7907.89772, 7908.82316, 7909.83359, 7910.92305,
                                            7911.92629, 7912.85512, 7913.74172, 7913.90775, 7914.83371, 7915.7958, 7916.7985,
                                            7916.94687, 7917.83057])
KIC_9964938_measured_RVs = np.array([ 5.57495362, -3.00884638, -12.58459465, 3.91660535, -6.59847981, -1.98196184, 3.63523816,
                                     -0.0317829, -3.13645823, -3.761185, -5.924885, -5.30984155, -8.53854242, 0.85977074,
                                     -6.15582926, -2.40106889])

KIC_10153521_time_of_observation = np.array([7904.89447, 7906.82581, 7907.87221, 7908.8724, 7909.81551, 7909.88805, 7910.83138,
                                             7911.7455, 7912.74676, 7913.77998, 7913.90372, 7914.73685, 7914.85624, 7915.7911,
                                             7915.86017, 7916.74143, 7916.86719, 7917.75341, 7917.85644])
KIC_10153521_measured_RVs = np.array([3.63409395, 0.52775362, 4.58720535, 5.92142019, 6.58253816, 7.77143816, -5.59406184, -9.0161829,
                                      7.77824177, 2.434815, 16.147215, -3.91394155, 9.74515845, -6.03184242, 5.06405758, -5.28622926,
                                      2.31767074, 0.91543111, 1.24343111])

KIC_10802309_time_of_observation = np.array([7912.89795, 7913.8551])
KIC_10802309_measured_RVs = np.array([-38.30475823, -31.852785])

KIC_11073910_time_of_observation = np.array([7912.93131, 7913.77032, 7913.92475])
KIC_11073910_measured_RVs = np.array([9.29234177, 10.867515, -4.347985])

KIC_11819949_time_of_observation = np.array([7904.86279, 7906.83648, 7907.82518, 7908.80419, 7908.94008, 7909.77312, 7909.92,
                                             7910.90473, 7911.76974, 7912.82445, 7913.75953, 7913.92019, 7914.7626, 7914.89177,
                                             7915.81434, 7915.91231, 7916.82087, 7917.7707, 7917.8912])
KIC_11819949_measured_RVs = np.array([-0.12440605, -1.83214638, 1.97710535, 1.25442019, 2.49032019, -4.18236184, 13.76743816,
                                      6.62913816, -19.1659829, 24.92514177, 3.872715, -4.712785, 14.05745845, 10.70645845, 
                                      -5.00974242, -1.34854242, 5.50417074, 18.23933111, 9.28503111])

KIC_12736892_time_of_observation = np.array([7904.86891, 7906.84933, 7907.87733, 7908.83552, 7909.8069, 7910.91476, 7911.86527, 
                                             7912.80223, 7913.85889, 7914.79027, 7914.91789, 7915.76303, 7915.89845, 7916.84284, 7917.83597]) 
KIC_12736892_measured_RVs = np.array([19.03479395, -17.33624638, 23.09790535, -9.94997981, 7.75093816, -8.53056184, -30.1987829,
                                      22.34844177, -6.493385, 5.50325845, 10.94235845, 8.83165758, 3.28695758, -19.84362926, 19.48293111])

time_array = KIC_6425783_time_of_observation *u.day

#This is a function to calculate the semi amplitude of RV - note need to researhch what this function actually is
def calculate_amplitude(mass1, q_array , inclination_angle_array, period): #m/s
    """ 
    Here is a function to calculate an array of amplitudes for a binary system, 
    the size of the array depends on the size of the input arrays
    The inputs consist of:
        mass of star1 (bigger object) in kg.
        'q' which is the mass2:mass1 ratio 
        angle of inclination
        and period of orbit in days (the code will change it to seconds)
        the amplitude output has units of m/sec
    """
    mass1 = (mass1)# mass1 is in kg
    mass2 = q_array * (mass1) # mass 2 is in kg
    Gravitational_constant = const.G # Value = 6.67 x 10^-11 m^3/(kg s^2)
    orbital_period = period * u.day
    orbital_period = orbital_period.to(u.second)#change the period from days to seconds
    total_mass = (mass1 + mass2) # in kg
    amplitude = (
                    ( ((mass2**3) / (total_mass**2)) 
                    * (np.sin(inclination_angle_array)**3)
                  *((2 * pi * Gravitational_constant) / orbital_period))**(1/3)
                )
    return amplitude #units of m/s

'''Here is a function to calculate the radial velocity of a star
    q_array is the mass2:mass1 ratio
    radial_velocity = Acos(wt+o),  where
    A=amplitude in m/s
    w=2*pi/period in days
    t=time_array
    o=theta_array
    period should be entered in days and the code will transform it to seconds,
    to get units of m/s
'''
def calculate_radial_velocity(mass1, q_array, time_array, period , theta_array, inclination_angle_array):#m/s
    #Specify inputs as quantity objects
    time_array = time_array.reshape(len(time_array), 1, 1, 1) # time array 
    theta_array = theta_array.reshape(1, len(theta_array), 1, 1) #Phase angle array
    q_array = q_array.reshape(1, 1, len(q_array),1)
    inclination_angle_array = inclination_angle_array.reshape(1, 1, 1, len(inclination_angle_array))
    amplitude = calculate_amplitude(mass1, q_array, inclination_angle_array, period) #m/s
    period = period * u.day
    omega=(2*pi) / (period.to(u.second)) #period is in seconds
    radial_velocity = ((amplitude*np.cos(((omega*time_array).to(u.radian, equivalencies = u.dimensionless_angles())) + theta_array))) #units of m/s
    return radial_velocity # in m/s

""" Here is a function to calculate an array of velocities with noise
    period should be entered in days to get units of m/s
    standard deviation needs to have units of m/sec
    'q'is the mass2:mass1 ratio
"""
def calculate_velocity_with_noise(mass1, q_array, time_array, period, std, theta_array,
                                  inclination_angle_array):#km/s
    radial_velocity_array = calculate_radial_velocity(mass1, q_array, time_array, period,theta_array,inclination_angle_array) #m/s
    noise = np.random.normal(0, std, size = radial_velocity_array.shape)  * (u.km/u.second)
    velocity_with_noise = radial_velocity_array + noise
    return (velocity_with_noise.to(u.km/u.s)) #km/sec


'''
    This function generates an array of radial velocity values using known quantities about the stars
    and some of the modeled quantities. This array is pretty big and we will use it in a chi-square
    test to find the best fit radial velocity curve for a given KIC target

9653110
9655045
9710336
9964938
10153521
11819949
12736892

'''
def generate_rv_array(KICnumber):
    std = 10 #our standard deviation will be 10 km/s
    time_array = KIC_6425783_time_of_observation*u.day # this gives us days of observation over a certain number of periods
    mass = md.target_mass(KICnumber) * SUN_MASS
    target_period = md.target_period(KICnumber) #* u.day
    inclination_angle_array = md.sin_distn()
    rv_array = calculate_velocity_with_noise(mass, q_array, time_array, target_period, std, theta_array, inclination_angle_array) #m/sec
    return rv_array # output in km/s


""" In this function we will be running a Chi Square statistical test to see if
    if any of our target stars fall within the range of a binary star distribution.
    In other words, this function will generate a chi-square value for each radial velocity curve
    in the array, then we will use the function after to find the minimum Chi-squre value
    which corresponds to the best-fit curve
    

9653110
9655045
9710336
9964938
10153521
11819949
12736892

""" #need to work out the units
def Chi_square_distribution(KICnumber):
    radial_velocities = generate_rv_array(KICnumber)
    measured_rv = KIC_6425783_measured_RVs * (u.km/u.s) #measured rv
    measured_rv = measured_rv.reshape(len(measured_rv), 1, 1, 1)
    chi_square = ((radial_velocities - measured_rv)**2) / ((10* u.km/u.s)**2)
    chi_values_array = np.sum(chi_square, axis = 0)
    return chi_values_array

      
'''
    This function is to find the minimum Chi-square value which corresponds tp the 
    best-fit radial velocity curve for a given KIC target with the given parameters
'''          
def minimum(KICnumber):
    s = Chi_square_distribution(KICnumber)
    min_value = s.min() #returns the minimum value in the multidimensional array s
    max_value = s.max()
    print("max value is ", max_value)
    print("min value is ", min_value)
    min_index = np.argwhere(s == min_value)
    # print("The index of the minimum value is ", min_index)
    theta_value = theta_array[min_index[0][0]]#0,0 specifies the value of theta and all rows
    q_value = q_array[min_index[0][1]]#0,1 specifies the value of q and all rows
    i_value = inclination_angle_array[min_index[0][2]]#0,1 specifies the value of q and all rows
    print("theta value is ", theta_value)
    print("q value is ", q_value)
    print("inclination angle value is ", i_value)
    best_fit_parameters = [theta_value,q_value,i_value,min_value]
    return best_fit_parameters#This returns an array for best-fit parameters of our model
      
def hist():
    df = theta_array.value
    plt.hist(df, bins = 500, color = 'slategrey', normed = 1)
    plt.title("Phase Angle Probability Distribution")
    

    

    
#np.argwhere(s) returns the shape of the matrix
# fig.savefig('my_figure.png') saves figures might use it to save plots to specific locations
# from IPython.display import Image
# Image('my_figure.png') displays a saved image or plot as defined above

            
            
# Have a graph that will output a single line for any single combination of inputs minus time time will always be constant
# The Astropy project is a library > go to documentation > from astropy import units as u.
# binary_semimajor_axis = 5 * u.AU
# if you have an array you have to go to speed.to()value and it will take the array out of the function
# When you plot you need to pull out the value as we did in the comment above.
# Adjust it so we can use the astropy framework to have the units
# have a file that calls the relevant functions as we need them.    