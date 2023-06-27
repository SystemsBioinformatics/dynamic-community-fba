import pytest
import cbmpy
import numpy as np
import endPointFBA.KineticModel as km
import endPointFBA.dynamic_fba as df
import endPointFBA.helpers.build_community_matrix as bm
from unittest import TestCase


@pytest.fixture
def model():
    model = bm.combine_models(
        [cbmpy.loadModel("data/bigg_models/e_coli_core.xml")]
    )
    model.createObjectiveFunction("R_BIOMASS_Ecoli_core_w_GAM")

    return model


@pytest.fixture
def kinetic_model(model):
    kinetic_model: km.KineticModel = km.KineticModel(
        model, {"R_GLCpts": [10, 5]}
    )

    return kinetic_model


# def test_create_external_species_reaction_dict(model):
#     dd = df.create_species_reactions_dict(model)
#     actual_dict = {
#         "M_akg_e": ["R_AKGt2r"],
#         "M_h_e": [
#             "R_AKGt2r",
#             "R_PIt2r",
#             "R_ACt2r",
#             "R_ATPS4r",
#             "R_PYRt2",
#             "R_SUCCt2_2",
#             "R_CYTBD",
#             "R_D_LACt2",
#             "R_SUCCt3",
#             "R_ETOHt2r",
#             "R_THD2",
#             "R_FORt2",
#             "R_FUMt2_2",
#             "R_GLUt2r",
#             "R_MALt2_2",
#             "R_NADH16",
#         ],
#         "M_pi_e": ["R_PIt2r"],
#         "M_acald_e": ["R_ACALDt"],
#         "M_ac_e": ["R_ACt2r"],
#         "M_pyr_e": ["R_PYRt2"],
#         "M_co2_e": ["R_CO2t"],
#         "M_succ_e": ["R_SUCCt2_2", "R_SUCCt3"],
#         "M_lac__D_e": ["R_D_LACt2"],
#         "M_etoh_e": ["R_ETOHt2r"],
#         "M_for_e": ["R_FORt2", "R_FORt"],
#         "M_fru_e": ["R_FRUpts2"],
#         "M_fum_e": ["R_FUMt2_2"],
#         "M_glc__D_e": ["R_GLCpts"],
#         "M_gln__L_e": ["R_GLNabc"],
#         "M_glu__L_e": ["R_GLUt2r"],
#         "M_h2o_e": ["R_H2Ot"],
#         "M_mal__L_e": ["R_MALt2_2"],
#         "M_nh4_e": ["R_NH4t"],
#         "M_o2_e": ["R_O2t"],
#     }
#     TestCase().assertDictEqual(dd, actual_dict)


# def test_if_core_model_is_correct(kinetic_model):
#     actual_dict = {
#         "M_glc__D_e": [
#             10,
#             9.8989898989899,
#             9.789653002983036,
#             9.671377964334022,
#             9.543523913323515,
#             9.405422557532395,
#             9.256381355663502,
#             9.095688023281728,
#             8.922616677699489,
#             8.73643598638444,
#             8.536419747903434,
#             8.321860405953542,
#             8.092086073712155,
#             7.846481723889006,
#             7.584515272741775,
#             7.3057693423959975,
#             7.009979506248541,
#             6.69707977725918,
#             6.367255942852327,
#             6.021007014852384,
#             5.659214450226936,
#             5.283217773161022,
#             4.894893616563057,
#             4.496732799918911,
#             4.0919066872096765,
#             3.6843096680988117,
#             3.2785594670731903,
#             2.8799321012988752,
#             2.4942058405118073,
#             2.1273921376924996,
#             1.7853459551050337,
#             1.4732768986050175,
#             1.1952241771573626,
#             0.9535997935489413,
#             0.7489214798464696,
#             0.5798262053202592,
#             0.4433788799757659,
#             0.3355825087662514,
#             0.25195294215427333,
#             0.18803033086169624,
#         ],
#         "biomass": [
#             0.1,
#             0.10861170230885982,
#             0.11793092919638619,
#             0.12800913943723263,
#             0.1389001648144169,
#             0.15066000231874085,
#             0.16334650920447258,
#             0.17701897792916652,
#             0.19173756359332036,
#             0.20756253143421777,
#             0.22455328621339843,
#             0.24276713902683425,
#             0.26225776032581066,
#             0.28307326110436803,
#             0.3052538379044956,
#             0.3288289125663649,
#             0.35381369622455133,
#             0.3802051116301716,
#             0.40797702261286745,
#             0.4370747505629889,
#             0.46740891410291147,
#             0.4988487218616511,
#             0.5312149952282383,
#             0.5642734162708527,
#             0.5977288020066721,
#             0.631221604392436,
#             0.6643282986626843,
#             0.6965677595630279,
#             0.7274159376114352,
#             0.7563308011452538,
#             0.7827881720959953,
#             0.8063264221614483,
#             0.8265941955842866,
#             0.8433915685959,
#             0.8566935495052769,
#             0.866626161507264,
#             0.8734572194144682,
#             0.8775791680886169,
#             0.8794270817451405,
#             0.8794308385200111,
#         ],
#         "M_atp_c": [0],
#         "M_f6p_c": [0],
#         "M_adp_c": [0],
#         "M_fdp_c": [0],
#         "M_h_c": [0],
#         "M_coa_c": [0],
#         "M_pyr_c": [0],
#         "M_accoa_c": [0],
#         "M_for_c": [0],
#         "M_g6p_c": [0],
#         "M_3pg_c": [0],
#         "M_13dpg_c": [0],
#         "M_6pgl_c": [0],
#         "M_h2o_c": [0],
#         "M_6pgc_c": [0],
#         "M_acald_c": [0],
#         "M_nad_c": [0],
#         "M_nadh_c": [0],
#         "M_akg_e": [
#             0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#         ],
#         "M_h_e": [
#             0,
#             -1.3754587658213804,
#             -2.8648960422063166,
#             -4.476801105618535,
#             -6.220100994992932,
#             -8.104139084350436,
#             -10.13864040130067,
#             -12.333660439179972,
#             -14.699513575320793,
#             -17.24667647669416,
#             -19.98566104473923,
#             -22.92685052960816,
#             -26.080291450458567,
#             -29.455432936316136,
#             -33.06080413232475,
#             -36.90361953836667,
#             -40.989301790168525,
#             -45.32091182041966,
#             -49.8984781082393,
#             -54.71822067485238,
#             -59.77167281787442,
#             -65.04471596308349,
#             -70.51656261036368,
#             -76.15875166928464,
#             -81.9342618378867,
#             -87.79690288157781,
#             -93.69120838718996,
#             -99.55311496906431,
#             -105.31174595122341,
#             -110.8925777943354,
#             -116.22209687147046,
#             -121.23370663174082,
#             -125.87414164572891,
#             -130.10913795374196,
#             -133.92688800012897,
#             -137.34943656948522,
#             -140.41871707184958,
#             -143.1746719354223,
#             -145.6634631488012,
#             -147.93173989968471,
#         ],
#         "M_akg_c": [0],
#         "M_2pg_c": [0],
#         "M_pi_e": [
#             0,
#             -0.03167986928360292,
#             -0.06596250923474664,
#             -0.10303722124774833,
#             -0.1431020363027964,
#             -0.1863629505299529,
#             -0.23303280341049465,
#             -0.28332971410802704,
#             -0.3374749751907509,
#             -0.39569028438706,
#             -0.4581941739932321,
#             -0.5251974743380198,
#             -0.5968976229105658,
#             -0.6734716056246461,
#             -0.7550672934992759,
#             -0.841792920657896,
#             -0.9337044443012671,
#             -1.0307905441539236,
#             -1.132955073085967,
#             -1.239996884896078,
#             -1.3515871723103938,
#             -1.4672447931124701,
#             -1.586310602946134,
#             -1.7079226164356018,
#             -1.830994943941961,
#             -1.9542049160784725,
#             -2.0759945122904355,
#             -2.1945938171045323,
#             -2.308075009691209,
#             -2.414444118173069,
#             -2.5117728486895623,
#             -2.5983630092053454,
#             -2.67292206729594,
#             -2.734714563393763,
#             -2.7836485605650885,
#             -2.8201876603367984,
#             -2.8453170730600306,
#             -2.8604804856476216,
#             -2.867278405615875,
#             -2.8672922256635913,
#         ],
#         "M_pi_c": [0],
#         "M_etoh_c": [0],
#         "M_acald_e": [
#             0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#         ],
#         "M_ac_c": [0],
#         "M_actp_c": [0],
#         "M_co2_c": [0],
#         "M_pep_c": [0],
#         "M_oaa_c": [0],
#         "M_cit_c": [0],
#         "M_acon_C_c": [0],
#         "M_icit_c": [0],
#         "M_ac_e": [
#             0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#         ],
#         "M_amp_c": [0],
#         "M_succoa_c": [0],
#         "M_e4p_c": [0],
#         "M_g3p_c": [0],
#         "M_gln__L_c": [0],
#         "M_glu__L_c": [0],
#         "M_nadph_c": [0],
#         "M_r5p_c": [0],
#         "M_nadp_c": [0],
#         "M_pyr_e": [
#             0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#         ],
#         "M_co2_e": [
#             0,
#             0.23958530845592246,
#             0.4990223246849814,
#             0.7797892806747166,
#             1.083440556298991,
#             1.4116029261304526,
#             1.765969493568055,
#             2.148290745044986,
#             2.5603620463075156,
#             3.0040067752445023,
#             3.481054141125114,
#             3.993310579422296,
#             4.542523438182031,
#             5.130335493739025,
#             5.758228664604586,
#             6.427455156906074,
#             7.138954212824856,
#             7.893252708467151,
#             8.690348157084156,
#             9.529573363302418,
#             10.409443254531933,
#             11.327486577850376,
#             12.280068571186362,
#             13.26221583437226,
#             14.267461842947007,
#             15.287741005684822,
#             16.313370283820994,
#             17.333168100122315,
#             18.334766023905726,
#             19.305161765708156,
#             20.231532211738674,
#             21.10226455007838,
#             21.908075646868713,
#             22.643001341323526,
#             23.30499877494937,
#             23.89688315205607,
#             24.424868019353,
#             24.89623465980735,
#             25.319373169869028,
#             25.702748966191486,
#         ],
#         "M_ru5p__D_c": [0],
#         "M_xu5p__D_c": [0],
#         "M_succ_e": [
#             0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#         ],
#         "M_succ_c": [0],
#         "M_o2_c": [0],
#         "M_q8h2_c": [0],
#         "M_q8_c": [0],
#         "M_lac__D_e": [
#             0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#         ],
#         "M_lac__D_c": [0],
#         "M_etoh_e": [
#             0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#         ],
#         "M_fum_c": [0],
#         "M_s7p_c": [0],
#         "M_dhap_c": [0],
#         "M_for_e": [
#             0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#         ],
#         "M_fru_e": [
#             0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#         ],
#         "M_fum_e": [
#             0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#         ],
#         "M_gln__L_e": [
#             0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#         ],
#         "M_glu__L_e": [
#             0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#         ],
#         "M_h2o_e": [
#             0,
#             0.3023163927545811,
#             0.6296383853231375,
#             0.9838190559912946,
#             1.36680491687313,
#             1.7806306470210886,
#             2.2274108052171155,
#             2.709327787872208,
#             3.228615154546701,
#             3.78753527922392,
#             4.3883500992179965,
#             5.033283526949372,
#             5.7244738674993725,
#             6.463914356927691,
#             7.253379721436101,
#             8.09433648760451,
#             8.987834701602987,
#             9.934378823625982,
#             10.933775980605338,
#             11.984960676303466,
#             13.085796748423194,
#             14.232860207379401,
#             15.421211082426954,
#             16.644169107855674,
#             17.893117528284428,
#             19.157371660721104,
#             20.424163342599478,
#             21.678806287883262,
#             22.90511467984249,
#             24.086137853570673,
#             25.20523437255477,
#             26.247428739671264,
#             27.20087840518292,
#             28.05816288360353,
#             28.81705726696564,
#             29.481294762939616,
#             30.059039788455785,
#             30.560432351832105,
#             30.997031804133364,
#             31.380434966306687,
#         ],
#         "M_mal__L_e": [
#             0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#             0.0,
#         ],
#         "M_nh4_e": [
#             0,
#             -0.04695789034975076,
#             -0.09777377072205454,
#             -0.15272823552334203,
#             -0.21211481870005233,
#             -0.27623886064363,
#             -0.34541584539014797,
#             -0.4199690828521591,
#             -0.5002265867616573,
#             -0.5865169714045027,
#             -0.6791641590644191,
#             -0.7784806556855218,
#             -0.8847591155045804,
#             -0.9982618781498982,
#             -1.1192081273256338,
#             -1.2477582944418746,
#             -1.3839953227732336,
#             -1.5279024326969999,
#             -1.6793371089034437,
#             -1.838001199869866,
#             -2.003407326820356,
#             -2.1748423105672114,
#             -2.3513291259805382,
#             -2.531590084241706,
#             -2.7140156115819822,
#             -2.8966451644310753,
#             -3.077169346947885,
#             -3.2529646793452787,
#             -3.421173624607634,
#             -3.57884059248484,
#             -3.723107344805043,
#             -3.851456714761946,
#             -3.9619728296819985,
#             -4.053565545239724,
#             -4.1260985867423745,
#             -4.1802591334668096,
#             -4.217507526023213,
#             -4.239983687753611,
#             -4.250059991339903,
#             -4.2500804762819175,
#         ],
#         "M_o2_e": [
#             0,
#             -0.22962931941664966,
#             -0.47829237744103925,
#             -0.7474079145713322,
#             -1.0384680757570437,
#             -1.3530348974497561,
#             -1.692734594276764,
#             -2.0592491046610766,
#             -2.454304249037278,
#             -2.879653732653403,
#             -3.3370580869338045,
#             -3.8282574899933732,
#             -4.3549372414693615,
#             -4.918684496576265,
#             -5.520934702603198,
#             -6.162906051088099,
#             -6.845520198619651,
#             -7.569307578911508,
#             -8.334295921241418,
#             -9.139881244176546,
#             -9.984681808937557,
#             -10.86637757050612,
#             -11.781540915202994,
#             -12.725469337821524,
#             -13.692037574947092,
#             -14.673595708846724,
#             -15.660950337737065,
#             -16.6434761132915,
#             -17.609410458433146,
#             -18.54637772650413,
#             -19.442160805978492,
#             -20.28568057341753,
#             -21.06806009735372,
#             -21.783566348869805,
#             -22.43018536236632,
#             -23.010586646737522,
#             -23.530674127987936,
#             -23.997275383580103,
#             -24.418277520663473,
#             -24.8016489737785,
#         ],
#         "M_mal__L_c": [0],
#         "M_nh4_c": [0],
#         "M_glx_c": [0],
#     }

#     ts = np.linspace(0, 15, 100)
#     y = df.dynamic_fba(
#         kinetic_model,
#         "R_BIOMASS_Ecoli_core_w_GAM",
#         ts,
#         {"M_glc__D_e": [10], "biomass": [0.1]},
#     )

#     TestCase().assertDictEqual(y, actual_dict)


# def test_if_community_model_is_correct():
#     pass
