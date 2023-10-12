import numpy as np

cur_dict = {
            "well_diameter": 0.25,
            "fluid_velocity": 2.0,
            "well_depth": 1500.0,
            "turb_outlet_temp": 60.0,
            "capacity_factor": 90.0,
            "fluid_density": 1000.0,
            "overall_pump_efficiency": 70.0
        }

# Defaults
geothermal_class = 1
plant_type = "binary"
if geothermal_class in [1,2]:
    plant_type = "flash"
t_amb = 50. # fahrenheit
p_amb = 14.7 # psia
turbine_eff = 0.8 # D32
t_wetbulb = 60. # fahrenheit
cooling_water_temp_rise  = 25.
pinchPt_cond = 7.5 # fahrenheit
pinchPt_twr = 5.0 # fahrenheit
tCond = t_wetbulb + cooling_water_temp_rise + pinchPt_cond + pinchPt_twr
# NCG Removal
mw_ncg = 44.0 # Mol Wt D148
mw_H2O = 18 # Mol Wt D149
ncg_level = 2000. # ppm D22
ncg_flow = 1000.0*ncg_level/1000000. # lb/hour D156
p_ncg_stage1 = (0.0000825*2000.+0.15)*0.49 # D163


steam_properties = [[0.0, -2.55175E-12, 2.41218E-08, -9.19096E-06, 0.001969537, -0.197885257, 8.089410675],             # Pressure psia
                    [1.01226E-14, -1.88058E-11, 1.49248E-08, -5.97605E-06, 0.001346286, 0.8382772, -24.1139345],        # Enthalpy btu/lb
                    [7.39915E-18, -1.29452E-14, 8.84301E-12, -1.84191E-09, -1.20262E-06, 0.002032431, -0.060089552],    # Entropy btu/lb-R
                    [1.40682E-18, -2.69957E-15, 2.17758E-12, -9.15282E-10, 2.2418E-07, -2.3968E-05, 0.017070952],       # Specific Volume ft^3/lb
                    [1.38098E-16, -2.18187E-13, 1.42033E-10, -4.65228E-08, 8.65998E-06, -0.000806161, 1.02617],         # Cp btu/lb-F
                    [-9.0287E-10, 3.4638E-07, -5.4475E-05, 0.00456759, -0.226287, 7.7497, 134.575],                     # T ft^3/lb
                    [-5.62597E-18, 1.09949E-14, -9.73487E-12, 5.19477E-09, -2.64887E-06, 0.000829554, 0.30074]]         # k btu/hr-ft-F

well_temp = {1:572, 2:392, 3:320., 4:0, 5:0, 999:0} 
tSource = well_temp[geothermal_class] # fahrenheit
t_c = (tSource-32)/1.8

pSource_sat = steam_properties[0][0]*tSource**6 + steam_properties[0][1]*tSource**5 + steam_properties[0][2]*tSource**4\
             + steam_properties[0][3]*tSource**3 + steam_properties[0][4]*tSource**2 + steam_properties[0][5]*tSource + steam_properties[0][6]# psia
hSource = steam_properties[1][0]*tSource**6 + steam_properties[1][1]*tSource**5 + steam_properties[1][2]*tSource**4\
             + steam_properties[1][3]*tSource**3 + steam_properties[1][4]*tSource**2 + steam_properties[1][5]*tSource + steam_properties[1][6]# btu/lb
sSource = steam_properties[2][0]*tSource**6 + steam_properties[2][1]*tSource**5 + steam_properties[2][2]*tSource**4\
             + steam_properties[2][3]*tSource**3 + steam_properties[2][4]*tSource**2 + steam_properties[2][5]*tSource + steam_properties[2][6] # btu/lb-R
vSource = steam_properties[3][0]*tSource**6 + steam_properties[3][1]*tSource**5 + steam_properties[3][2]*tSource**4\
             + steam_properties[3][3]*tSource**3 + steam_properties[3][4]*tSource**2 + steam_properties[3][5]*tSource + steam_properties[3][6] # ft^3/lb

sio2 = -0.0000001334837*t_c**4 + 0.0000706584462*t_c**3 - 0.0036294799613*t_c**2 + 0.3672417729236*t_c + 4.205944351495
tAmphSiO2 = 0.0000000000249634*sio2**4 - 0.00000000425191*sio2**3 - 0.000119669*sio2**2 + 0.307616*sio2 - 0.294394 
tSi = tAmphSiO2*1.8 + 32 + 1
tSink = max((tSi-32)/1.8, (t_amb-32)/1.8+15+5+5+2)*1.8 + 32
# in the GETEM code they compare this to the calculated binary outlet temp at the plant 
# brine effectiveness but we can't estimate the 2nd law eff in the same way so 
# we go with this number
hSink = steam_properties[1][0]*tSink**6 + steam_properties[1][1]*tSink**5 + steam_properties[1][2]*tSink**4\
             + steam_properties[1][3]*tSink**3 + steam_properties[1][4]*tSink**2 + steam_properties[1][5]*tSink + steam_properties[1][6]# btu/lb
sSink = steam_properties[2][0]*tSink**6 + steam_properties[2][1]*tSink**5 + steam_properties[2][2]*tSink**4\
             + steam_properties[2][3]*tSink**3 + steam_properties[2][4]*tSink**2 + steam_properties[2][5]*tSink + steam_properties[2][6] # btu/lb-R
vSink = steam_properties[3][0]*tSink**6 + steam_properties[3][1]*tSink**5 + steam_properties[3][2]*tSink**4\
             + steam_properties[3][3]*tSink**3 + steam_properties[3][4]*tSink**2 + steam_properties[3][5]*tSink + steam_properties[3][6] # ft^3/lb

availableEnergy = ((hSource-hSink)-(tSink+460)*(sSource-sSink))/3.413 # w-hr/lb  | *3.413 -> btu/lb
            # maximum power that could be produced with IDEAL power cycle

''' Flash Plant Performance '''
# Fraction of inlet GF Injection 
if plant_type == "binary":
    frac_inlet = 1
    second_law_eff = 0.0000000550046*t_c**3 - 0.0000407813*\
                      t_c**2 + 0.0101746*t_c - 0.413853
else:
    second_law_eff = -0.0000000000866018*t_c**4 + 0.000000179053*t_c**3 - \
        0.000105902*t_c**2 + 0.0248737*t_c - 1.66615
    bin1 = [125, 325, 675]
    bin2 = [2, 20, 200, 1000]
    tHpEst = tSource-0.5*(tSource - 212) # T52
    pHP = -0.0000000000000280188*tHpEst**6 + 0.0000000000434771*tHpEst**5 - \
          0.00000000675074*tHpEst**4 + 0.00000158988*tHpEst**3 - \
          0.0000909892*tHpEst**2 + 0.00600322*tHpEst - 0.0607916 # T55
    pLP = -0.0000000000000280188*212.**6 + 0.0000000000434771*212.**5 - \
          0.00000000675074*212.**4 + 0.00000158988*212.**3 - \
          0.0000909892*212.**2 + 0.00600322*212. - 0.0607916 # T54
    t57 = pHP + 1
    qq = np.searchsorted(bin2, t57, side='right')
    flash_temp1_vals = [-7.1864066799*t57**6 + 63.304761377*t57**5 - 222.30982965*t57**4 + 400.57269432*t57**3 - 403.56297354*t57**2 + 255.85632577*t57 + 14.788238833,
                        -0.000012178121702*t57**6 + 0.00094358038872*t57**5 - 0.029734376328*t57**4 + 0.49468791547*t57**3 - 4.8016701723*t57**2 + 31.491049082*t57 + 78.871966537,
                        -1.2878309875E-11*t57**6 + 0.00000001053096688*t57**5 - 0.0000034988475881*t57**4 + 0.00061292292067*t57**3 - 0.062604066919*t57**2 + 4.3688747745*t57 + 161.40853789,
                        -4.3371351471E-16*t57**6 + 1.8867165569E-12*t57**5 - 0.0000000034275245432*t57**4 + 0.0000034048164769*t57**3 - 0.0020724712921*t57**2 + 0.93056131917*t57 + 256.29706201,
                        -3.2288675486E-19*t57**6 + 4.589696886E-15*t57**5 - 2.7823504188E-11*t57**4 + 0.000000094407417758*t57**3 - 0.00020256473758*t57**2 + 0.33345911089*t57 + 342.90613285] # bin2
    flash_temp1 = flash_temp1_vals[qq]
    flash_pressure1 = -0.0000000000000280188*flash_temp1**6 + \
                        0.0000000000434771*flash_temp1**5 - 0.00000000675074*flash_temp1**4 + \
                        0.00000158988*flash_temp1**3 - 0.0000909892*flash_temp1**2 + \
                        0.00600322*flash_temp1 - 0.0607916 + 1
    flash_pressure2 = -0.0000000000000280188*212**6 + \
                        0.0000000000434771*212**5 - 0.00000000675074*212**4 + \
                        0.00000158988*212**3 - 0.0000909892*212**2 + \
                        0.00600322*212 - 0.0607916 + 1
    
    # Define flash conditions
    hf1_vals = [ -0.000000004480902*flash_temp1**4 + 0.0000020320904*flash_temp1**3 - 0.00034115062*flash_temp1**2 + 1.0234315*flash_temp1 - 32.479184,
                 0.00000000025563678*flash_temp1**4 + 0.000000073480055*flash_temp1**3 - 0.000027703224*flash_temp1**2 + 0.9998551*flash_temp1 - 31.760088,
                 5.8634263518E-11*flash_temp1**5 - 0.00000013378773724*flash_temp1**4 + 0.00012227602697*flash_temp1**3 - 0.055373746094*flash_temp1**2 + 13.426933583*flash_temp1 - 1137.0718729,
                 0.000037021613128*flash_temp1**5 - 0.12714518982*flash_temp1**4 + 174.6587566*flash_temp1**3 - 119960.00955*flash_temp1**2 + 41194401.715*flash_temp1 - 5658291651.7]
    hg1_vals = [ -7.2150559138E-10*flash_temp1**4 - 0.00000015844186585*flash_temp1**3 - 0.000030268712038*flash_temp1**2 + 0.44148580795*flash_temp1 + 1061.0996074,
                 -0.0000000005035389718*flash_temp1**4 - 0.00000051596852593*flash_temp1**3 + 0.000099006018886*flash_temp1**2 + 0.42367961566*flash_temp1 + 1061.9537518,
                 -4.9118123157E-13*flash_temp1**6 + 0.0000000013698021251*flash_temp1**5 - 0.0000015842735401*flash_temp1**4 + 0.00096963380389*flash_temp1**3 - 0.33157805684*flash_temp1**2 + 60.38391862*flash_temp1 - 3413.791688,
                 -0.000048138033984*flash_temp1**5 + 0.16531315908*flash_temp1**4 - 227.07686319*flash_temp1**3 + 155953.29919*flash_temp1**2 - 53551582.984*flash_temp1 + 7355226428.1]
    qq = np.searchsorted(bin1, flash_temp1, side='right') 

    hf1 = hf1_vals[qq]
    hg1 = hg1_vals[qq]
    
    qq = np.searchsorted(bin2, flash_pressure2, side='right')
    flash_temp2_vals = [-7.1864066799*flash_pressure2**6 + 63.304761377*flash_pressure2**5 - 222.30982965*flash_pressure2**4 + 400.57269432*flash_pressure2**3 - 403.56297354*flash_pressure2**2 + 255.85632577*flash_pressure2 + 14.788238833,
                        -0.000012178121702*flash_pressure2**6 + 0.00094358038872*flash_pressure2**5 - 0.029734376328*flash_pressure2**4 + 0.49468791547*flash_pressure2**3 - 4.8016701723*flash_pressure2**2 + 31.491049082*flash_pressure2 + 78.871966537,
                        -1.2878309875E-11*flash_pressure2**6 + 0.00000001053096688*flash_pressure2**5 - 0.0000034988475881*flash_pressure2**4 + 0.00061292292067*flash_pressure2**3 - 0.062604066919*flash_pressure2**2 + 4.3688747745*flash_pressure2 + 161.40853789,
                        -4.3371351471E-16*flash_pressure2**6 + 1.8867165569E-12*flash_pressure2**5 - 0.0000000034275245432*flash_pressure2**4 + 0.0000034048164769*flash_pressure2**3 - 0.0020724712921*flash_pressure2**2 + 0.93056131917*flash_pressure2 + 256.29706201,
                        -3.2288675486E-19*flash_pressure2**6 + 4.589696886E-15*flash_pressure2**5 - 2.7823504188E-11*flash_pressure2**4 + 0.000000094407417758*flash_pressure2**3 - 0.00020256473758*flash_pressure2**2 + 0.33345911089*flash_pressure2 + 342.90613285]
    flash_temp2 = flash_temp2_vals[qq]

    
    qq = np.searchsorted(bin1, flash_temp1, side='right') 
    hf2_vals = [-0.000000004480902*flash_temp2**4 + 0.0000020320904*flash_temp2**3 - 0.00034115062*flash_temp2**2 + 1.0234315*flash_temp2 - 32.479184,
                0.00000000025563678*flash_temp2**4 + 0.000000073480055*flash_temp2**3 - 0.000027703224*flash_temp2**2 + 0.9998551*flash_temp2 - 31.760088,
                5.8634263518E-11*flash_temp2**5 - 0.00000013378773724*flash_temp2**4 + 0.00012227602697*flash_temp2**3 - 0.055373746094*flash_temp2**2 + 13.426933583*flash_temp2 - 1137.0718729,
                0.000037021613128*flash_temp2**5 - 0.12714518982*flash_temp2**4 + 174.6587566*flash_temp2**3 - 119960.00955*flash_temp2**2 + 41194401.715*flash_temp2 - 5658291651.7]
    hg2_vals = [-7.2150559138E-10*flash_temp2**4 - 0.00000015844186585*flash_temp2**3 - 0.000030268712038*flash_temp2**2 + 0.44148580795*flash_temp2 + 1061.0996074,
                -0.0000000005035389718*flash_temp2**4 - 0.00000051596852593*flash_temp2**3 + 0.000099006018886*flash_temp2**2 + 0.42367961566*flash_temp2 + 1061.9537518,
                -4.9118123157E-13*flash_temp2**6 + 0.0000000013698021251*flash_temp2**5 - 0.0000015842735401*flash_temp2**4 + 0.00096963380389*flash_temp2**3 - 0.33157805684*flash_temp2**2 + 60.38391862*flash_temp2 - 3413.791688,
                -0.000048138033984*flash_temp2**5 + 0.16531315908*flash_temp2**4 - 227.07686319*flash_temp2**3 + 155953.29919*flash_temp2**2 - 53551582.984*flash_temp2 + 7355226428.1]

    hf2 = hf2_vals[qq]
    hg2 = hg2_vals[qq]
    x1 = (hSource-hf1)/(hg1-hf1)
    x2 = (hf1-hf2)/(hg2-hf2)
    m_steam1 = x1*1000.0
    m_steam2 = x2*1000.0*(1-x1)

    # Define look up tables for Pressure, enthalpy, entropy, and temperature for a variety of things
    p_CondSat_table = [lambda temp: 9.97153E-15*temp**6 + 0.0000000000168375*temp**5 + 0.00000000113502*temp**4 + 0.000000341917*temp**3 + 0.0000184319*temp**2 + 0.00111108*temp + 0.0212481,
                       lambda temp: -0.0000000000000280188*temp**6 + 0.0000000000434771*temp**5 - 0.00000000675074*temp**4 + 0.00000158988*temp**3 - 0.0000909892*temp**2 + 0.00600322*temp - 0.0607916,
                       lambda temp: 0.00000000000024303*temp**6 - 0.000000000662939*temp**5 + 0.000000768058*temp**4 - 0.000453695*temp**3 + 0.150475*temp**2 - 26.49*temp + 1934.47,
                       lambda temp: 0.000000000000722829*temp**6 - 0.00000000201101*temp**5 + 0.00000225787*temp**4 - 0.00125136*temp**3 + 0.347131*temp**2 - 38.4978*temp] # Bin1  D153->U67
    hGF_in_table = [lambda temp: -0.000000004480902*temp**4 + 0.0000020320904*temp**3 - 0.00034115062*temp**2 + 1.0234315*temp - 32.479184,
                    lambda temp: 0.00000000025563678*temp**4 + 0.000000073480055*temp**3 - 0.000027703224*temp**2 + 0.9998551*temp - 31.760088,
                    lambda temp: 5.8634263518E-11*temp**5 - 0.00000013378773724*temp**4 + 0.00012227602697*temp**3 - 0.055373746094*temp**2 + 13.426933583*temp - 1137.0718729,
                    lambda temp: 0.000037021613128*temp**5 - 0.12714518982*temp**4 + 174.6587566*temp**3 - 119960.00955*temp**2 + 41194401.715*temp - 5658291651.7] # Bin1 D70->V66
    flash1_hf_table = [lambda temp: -0.000000004480902*temp**4 + 0.0000020320904*temp**3 - 0.00034115062*temp**2 + 1.0234315*temp - 32.479184,
                       lambda temp: 0.00000000025563678*temp**4 + 0.000000073480055*temp**3 - 0.000027703224*temp**2 + 0.9998551*temp - 31.760088,
                       lambda temp: 5.8634263518E-11*temp**5 - 0.00000013378773724*temp**4 + 0.00012227602697*temp**3 - 0.055373746094*temp**2 + 13.426933583*temp - 1137.0718729,
                       lambda temp: 0.000037021613128*temp**5 - 0.12714518982*temp**4 + 174.6587566*temp**3 - 119960.00955*temp**2 + 41194401.715*temp - 5658291651.7] # Bin1  D82->V73
    flash1_hg_table = [lambda temp: -7.2150559138E-10*temp**4 - 0.00000015844186585*temp**3 - 0.000030268712038*temp**2 + 0.44148580795*temp + 1061.0996074,
                       lambda temp: -0.0000000005035389718*temp**4 - 0.00000051596852593*temp**3 + 0.000099006018886*temp**2 + 0.42367961566*temp + 1061.9537518,
                       lambda temp: -4.9118123157E-13*temp**6 + 0.0000000013698021251*temp**5 - 0.0000015842735401*temp**4 + 0.00096963380389*temp**3 - 0.33157805684*temp**2 + 60.38391862*temp - 3413.791688,
                       lambda temp: -0.000048138033984*temp**5 + 0.16531315908*temp**4 - 227.07686319*temp**3 + 155953.29919*temp**2 - 53551582.984*temp + 7355226428.1] # Bin1  D83->V74
    flash2_hf_table = [lambda temp: -0.000000004480902*temp**4 + 0.0000020320904*temp**3 - 0.00034115062*temp**2 + 1.0234315*temp - 32.479184,
                       lambda temp: 0.00000000025563678*temp**4 + 0.000000073480055*temp**3 - 0.000027703224*temp**2 + 0.9998551*temp - 31.760088,
                       lambda temp: 5.8634263518E-11*temp**5 - 0.00000013378773724*temp**4 + 0.00012227602697*temp**3 - 0.055373746094*temp**2 + 13.426933583*temp - 1137.0718729,
                       lambda temp: 0.000037021613128*temp**5 - 0.12714518982*temp**4 + 174.6587566*temp**3 - 119960.00955*temp**2 + 41194401.715*temp - 5658291651.7] # Bin1 D90->V77
    flash2_hg_table = [lambda temp: -7.2150559138E-10*temp**4 - 0.00000015844186585*temp**3 - 0.000030268712038*temp**2 + 0.44148580795*temp + 1061.0996074,
                       lambda temp: -0.0000000005035389718*temp**4 - 0.00000051596852593*temp**3 + 0.000099006018886*temp**2 + 0.42367961566*temp + 1061.9537518,
                       lambda temp: -4.9118123157E-13*temp**6 + 0.0000000013698021251*temp**5 - 0.0000015842735401*temp**4 + 0.00096963380389*temp**3 - 0.33157805684*temp**2 + 60.38391862*temp - 3413.791688,
                       lambda temp: -0.000048138033984*temp**5 + 0.16531315908*temp**4 - 227.07686319*temp**3 + 155953.29919*temp**2 - 53551582.984*temp + 7355226428.1] # Bin1 D91 -> V78
    turb_perf_hTin_table = [lambda temp: -7.2150559138E-10*temp**4 - 0.00000015844186585*temp**3 - 0.000030268712038*temp**2 + 0.44148580795*temp + 1061.0996074,
                            lambda temp: -0.0000000005035389718*temp**4 - 0.00000051596852593*temp**3 + 0.000099006018886*temp**2 + 0.42367961566*temp + 1061.9537518,
                            lambda temp: -4.9118123157E-13*temp**6 + 0.0000000013698021251*temp**5 - 0.0000015842735401*temp**4 + 0.00096963380389*temp**3 - 0.33157805684*temp**2 + 60.38391862*temp - 3413.791688,
                            lambda temp: -0.000048138033984*temp**5 + 0.16531315908*temp**4 - 227.07686319*temp**3 + 155953.29919*temp**2 - 53551582.984*temp + 7355226428.1] # Bin1  L67->V76
    turb_perf_sCondLiq_table = [lambda temp: 0.000000002964634*temp**3 - 0.000002499642*temp**2 + 0.002195168*temp- 0.06778459,
                                lambda temp: -0.0000000000009593628*temp**4 + 0.000000002223552*temp**3 - 0.00000216312*temp**2 + 0.00215356*temp - 0.06615222,
                                lambda temp: 5.297508E-16*temp**6 - 0.000000000001495185*temp**5 + 0.000000001739088*temp**4 - 0.000001064881*temp**3 + 0.0003612208*temp**2 - 0.06296176*temp + 4.729245,
                                lambda temp: 0.0000000001100136*temp**6 - 0.0000004450574*temp**5 + 0.0007501037*temp**4 - 0.674174*temp**3 + 340.7939*temp**2 - 91866.44*temp + 10317090] # Bin1 L70->W69
    turb_perf_sCondVap_table = [lambda temp: 0.0000000000238219*temp**4 - 0.00000002213415*temp**3 + 0.00001113945*temp**2 - 0.004205795*temp + 2.312154,
                                lambda temp: 0.000000000008392618*temp**4 - 0.00000001361597*temp**3 + 0.000009334758*temp**2 - 0.004032959*temp + 2.305898,
                                lambda temp: -8.283605E-16*temp**6 + 0.000000000002333203*temp**5 - 0.000000002708314*temp**4 + 0.00000165401*temp**3 - 0.0005589582*temp**2 + 0.09784205*temp - 5.19791,
                                lambda temp: -0.000000000141418*temp**6 + 0.0000005720121*temp**5 - 0.0009639254*temp**4 + 0.8662219*temp**3 - 437.8104*temp**2 + 118002.2*temp - 13250460] # Bin1 L69->W71
    turb_perf_sTin_table = [lambda temp: 0.0000000000238219*temp**4 - 0.00000002213415*temp**3 + 0.00001113945*temp**2 - 0.004205795*temp + 2.312154,
                            lambda temp: 0.000000000008392618*temp**4 - 0.00000001361597*temp**3 + 0.000009334758*temp**2 - 0.004032959*temp + 2.305898,
                            lambda temp: -8.283605E-16*temp**6 + 0.000000000002333203*temp**5 - 0.000000002708314*temp**4 + 0.00000165401*temp**3 - 0.0005589582*temp**2 + 0.09784205*temp - 5.19791,
                            lambda temp: -0.000000000141418*temp**6 + 0.0000005720121*temp**5 - 0.0009639254*temp**4 + 0.8662219*temp**3 - 437.8104*temp**2 + 118002.2*temp - 13250460] # Bin1 L68->W75
    flash_hCond_g_table = [ lambda temp: -7.2150559138E-10*temp**4 - 0.00000015844186585*temp**3 - 0.000030268712038*temp**2 + 0.44148580795*temp + 1061.0996074,
                            lambda temp: -0.0000000005035389718*temp**4 - 0.00000051596852593*temp**3 + 0.000099006018886*temp**2 + 0.42367961566*temp + 1061.9537518,
                            lambda temp: -4.9118123157E-13*temp**6 + 0.0000000013698021251*temp**5 - 0.0000015842735401*temp**4 + 0.00096963380389*temp**3 - 0.33157805684*temp**2 + 60.38391862*temp - 3413.791688,
                            lambda temp: -0.000048138033984*temp**5 + 0.16531315908*temp**4 - 227.07686319*temp**3 + 155953.29919*temp**2 - 53551582.984*temp + 7355226428.1] # Bin1 D78->V70
    flash_hCond_f_table = [ lambda temp: -0.000000004480902*temp**4 + 0.0000020320904*temp**3 - 0.00034115062*temp**2 + 1.0234315*temp - 32.479184,
                            lambda temp: 0.00000000025563678*temp**4 + 0.000000073480055*temp**3 - 0.000027703224*temp**2 + 0.9998551*temp - 31.760088,
                            lambda temp: 5.8634263518E-11*temp**5 - 0.00000013378773724*temp**4 + 0.00012227602697*temp**3 - 0.055373746094*temp**2 + 13.426933583*temp - 1137.0718729,
                            lambda temp: 0.0000000000007261395*temp**6 - 0.00000000216551*temp**5 + 0.000002698676*temp**4 - 0.001765175*temp**3 + 0.6471745*temp**2 - 125.9218*temp + 10153.58] # Bin1 D77->V68
    lpTurb_sTin_table = [lambda temp: 0.0000000000238219*temp**4 - 0.00000002213415*temp**3 + 0.00001113945*temp**2 - 0.004205795*temp + 2.312154,
                         lambda temp: 0.000000000008392618*temp**4 - 0.00000001361597*temp**3 + 0.000009334758*temp**2 - 0.004032959*temp + 2.305898,
                         lambda temp: -8.283605E-16*temp**6 + 0.000000000002333203*temp**5 - 0.000000002708314*temp**4 + 0.00000165401*temp**3 - 0.0005589582*temp**2 + 0.09784205*temp - 5.19791,
                         lambda temp: -0.000000000141418*temp**6 + 0.0000005720121*temp**5 - 0.0009639254*temp**4 + 0.8662219*temp**3 - 437.8104*temp**2 + 118002.2*temp - 13250460] # Bin1 L85->W79
    lpTurb_hTin_table = [lambda temp: -7.2150559138E-10*temp**4 - 0.00000015844186585*temp**3 - 0.000030268712038*temp**2 + 0.44148580795*temp + 1061.0996074,
                         lambda temp: -0.0000000005035389718*temp**4 - 0.00000051596852593*temp**3 + 0.000099006018886*temp**2 + 0.42367961566*temp + 1061.9537518,
                         lambda temp: -4.9118123157E-13*temp**6 + 0.0000000013698021251*temp**5 - 0.0000015842735401*temp**4 + 0.00096963380389*temp**3 - 0.33157805684*temp**2 + 60.38391862*temp - 3413.791688,
                         lambda temp: -0.000048138033984*temp**5 + 0.16531315908*temp**4 - 227.07686319*temp**3 + 155953.29919*temp**2 - 53551582.984*temp + 7355226428.1] # Bin1 L84->V80    
    u90_table = [lambda temp: -7.1864066799*temp**6 + 63.304761377*temp**5 - 222.30982965*temp**4 + 400.57269432*temp**3 - 403.56297354*temp**2 + 255.85632577*temp + 14.788238833,
                 lambda temp: -0.000012178121702*temp**6 + 0.00094358038872*temp**5 - 0.029734376328*temp**4 + 0.49468791547*temp**3 - 4.8016701723*temp**2 + 31.491049082*temp + 78.871966537,
                 lambda temp: -1.2878309875E-11*temp**6 + 0.00000001053096688*temp**5 - 0.0000034988475881*temp**4 + 0.00061292292067*temp**3 - 0.062604066919*temp**2 + 4.3688747745*temp + 161.40853789,
                 lambda temp: -4.3371351471E-16*temp**6 + 1.8867165569E-12*temp**5 - 0.0000000034275245432*temp**4 + 0.0000034048164769*temp**3 - 0.0020724712921*temp**2 + 0.93056131917*temp + 256.29706201,
                 lambda temp: -3.2288675486E-19*temp**6 + 4.589696886E-15*temp**5 - 2.7823504188E-11*temp**4 + 0.000000094407417758*temp**3 - 0.00020256473758*temp**2 + 0.33345911089*temp + 342.90613285] # Bin2  U90
    u93_table = [lambda temp: -7.1864066799*temp**6 + 63.304761377*temp**5 - 222.30982965*temp**4 + 400.57269432*temp**3 - 403.56297354*temp**2 + 255.85632577*temp + 14.788238833,
                 lambda temp: -0.000012178121702*temp**6 + 0.00094358038872*temp**5 - 0.029734376328*temp**4 + 0.49468791547*temp**3 - 4.8016701723*temp**2 + 31.491049082*temp + 78.871966537,
                 lambda temp: -1.2878309875E-11*temp**6 + 0.00000001053096688*temp**5 - 0.0000034988475881*temp**4 + 0.00061292292067*temp**3 - 0.062604066919*temp**2 + 4.3688747745*temp + 161.40853789,
                 lambda temp: -4.3371351471E-16*temp**6 + 1.8867165569E-12*temp**5 - 0.0000000034275245432*temp**4 + 0.0000034048164769*temp**3 - 0.0020724712921*temp**2 + 0.93056131917*temp + 256.29706201,
                 lambda temp: -3.2288675486E-19*temp**6 + 4.589696886E-15*temp**5 - 2.7823504188E-11*temp**4 + 0.000000094407417758*temp**3 - 0.00020256473758*temp**2 + 0.33345911089*temp + 342.90613285] # Bin2 U93
        
    qq = np.searchsorted(bin1, tCond, side='right')
    h_condF = flash_hCond_f_table[qq](tCond) # D77
    h_condG = flash_hCond_g_table[qq](tCond) # D78
    qq = np.searchsorted(bin1, tSource-1, side='right') 
    hGF_in = hGF_in_table[qq](tSource-1) # btu/lb D70    
    qq = np.searchsorted(bin1, tCond, side='right')        
    ncg_Ptotal = ((0.0000825*ncg_level+0.15)*0.49) + p_CondSat_table[qq](tCond) # D154/D169
    ncg_Pratio = np.exp(np.log(p_amb/ncg_Ptotal)/3) # per stage D155

    # NCG Removal Stage 1 160-175
    ncg_stage1_Pinter1 = ncg_Pratio*ncg_Ptotal # D160
    ncg_stage1_PrJet = ncg_stage1_Pinter1/ncg_Ptotal   # D161
    # D162 = tCond + 460
    ncg_stage1_h2o_vent_flow = (p_CondSat_table[qq](tCond)/p_ncg_stage1)*(mw_H2O/mw_ncg)*ncg_flow # lb/hour D164
    ncg_stage1_totalVentFlow = ncg_flow + ncg_stage1_h2o_vent_flow # D165
    ncg_stage1_amntNCG = ncg_flow/mw_ncg # D166
    ncg_stage1_amntH2O = ncg_stage1_h2o_vent_flow/mw_H2O # D167
    ncg_stage1_MWvent = ncg_stage1_totalVentFlow/(ncg_stage1_amntNCG+ncg_stage1_amntH2O) # D168    
    ncg_stage1_PsucPstm = ncg_Ptotal/(pHP+1.) # D171
    ncg_stage1_area_ratio = ((3.5879*(ncg_stage1_PrJet)**-2.1168)+0.1)*ncg_stage1_PsucPstm**(-1.155*ncg_stage1_PrJet**-0.0453) # D172
    entrainmentRatioDesign = (1.0035*ncg_stage1_area_ratio+8.9374)*ncg_stage1_PsucPstm**(2.9594*ncg_stage1_area_ratio**-0.8458+0.99) # D173
    entrainmentRatio = entrainmentRatioDesign*(((460+flash_temp1)*ncg_stage1_MWvent)/((tCond + 460)*mw_H2O))**0.5 # D174
    ncg_stage1_supplySteam = ncg_stage1_totalVentFlow/entrainmentRatio # D175

    stage1_pNCG = (ncg_Pratio*ncg_Ptotal)-p_CondSat_table[qq](tCond) # psia D187
    stage1_moleRatio = stage1_pNCG/p_CondSat_table[qq](tCond) # ncg/h2o D188
    stage1_flowNCG = ncg_flow/mw_ncg # mole/h D189
    steam_leaving_stage1 = (stage1_flowNCG/stage1_moleRatio)*mw_H2O # lb/h D191
    stage1_steamCond = ncg_stage1_supplySteam+ncg_stage1_h2o_vent_flow-steam_leaving_stage1 # lb/h D192
    stage1_dhCondSteam = h_condG-h_condF # D193
    stage1_Q_rejected = stage1_steamCond*stage1_dhCondSteam # btu/h D194

    # NCG Removal Stage 2 209-224
    ncg_stage2_Pinter2 = ncg_Pratio*ncg_stage1_Pinter1 # D209  P,inter-2
    ncg_stage2_PrJet = ncg_stage2_Pinter2/ncg_stage1_Pinter1 # D210  Pr,jet
    # D211  T = tCond+460
    ncg_stage2_h2o_vent_flow = ((ncg_flow/mw_ncg)/(((ncg_Pratio*ncg_Ptotal)-
                                p_CondSat_table[qq](tCond))/p_CondSat_table[qq](tCond)))*mw_H2O # D213  vent flow,h2o
    ncg_stage2_totalVentFlow = ncg_flow + ncg_stage2_h2o_vent_flow # D214  total vent flow
    ncg_stage2_amntNCG = ncg_flow/mw_ncg # D215  amount ncg
    ncg_stage2_amntH2O = (ncg_flow/mw_ncg)/(((ncg_Pratio*ncg_Ptotal)-
                            p_CondSat_table[qq](tCond))/p_CondSat_table[qq](tCond)) # D216  amount h2o
    ncg_stage2_MWvent = ncg_stage2_totalVentFlow/(ncg_stage2_amntNCG+ncg_stage2_amntH2O) # D217  MW,vent
    ncg_stage2_PsucPstm = ncg_stage1_Pinter1/(pHP + 1) # D220  Psuc/Pstm
    ncg_stage2_area_ratio = ((3.5879*(ncg_stage2_PrJet)**-2.1168)+0.1)*ncg_stage2_PsucPstm**(-1.155*ncg_stage2_PrJet**-0.0453) # D221  Area Ratio
    entrainmentRatioDesign2 = (1.0035*ncg_stage2_area_ratio+8.9374)*ncg_stage2_PsucPstm**(2.9594*ncg_stage2_area_ratio**-0.8458+0.99) # D222  Entrainment Ratio,design
    entrainmentRatio2 = entrainmentRatioDesign2*(((460+flash_temp1)*ncg_stage2_MWvent)/((tCond+460)*mw_H2O))**0.5 # D223  Entrainment Ratio
    ncg_stage2_supplySteam = ncg_stage2_totalVentFlow/entrainmentRatio2 # D224  supply stm flow_2nd stage

    stage2_pNCG = ncg_stage2_Pinter2-p_CondSat_table[qq](tCond) # psia D236
    stage2_moleRatio = stage2_pNCG/p_CondSat_table[qq](tCond) # ncg/h2o D237
    stage2_flowNCG = ncg_flow/mw_ncg # mole/h D238
    steam_leaving_stage2 = stage2_flowNCG/stage2_moleRatio # lb/h D239
    stage2_steamCond = ncg_stage2_supplySteam+ncg_stage2_h2o_vent_flow-\
        (steam_leaving_stage2*mw_H2O) # lb/h D241
    stage2_Q_rejected = stage2_steamCond*stage1_dhCondSteam # btu/h D243

    stage2_pNCG = p_amb-p_CondSat_table[qq](tCond) # D284
    stage2_moleRatio = stage2_pNCG/p_CondSat_table[qq](tCond) # D285
    steam_leaving_stage2 = stage2_flowNCG/stage2_moleRatio # mole/hour D287 

    totSteamFlow = ncg_stage1_supplySteam + ncg_stage2_supplySteam # D306
    qRejected_NCG_system = stage1_Q_rejected+stage2_Q_rejected # btu/hour D307
    
    # Turbine Perfomance 
    qq = np.searchsorted(bin1, flash_temp1, side='right')
    x_flash1 = (hGF_in-flash1_hf_table[qq](flash_temp1))/(flash1_hg_table[qq](flash_temp1)-flash1_hf_table[qq](flash_temp1))# D84
    mass_steam_flash1 = x_flash1*1000.0 # D86

    qq2 = np.searchsorted(bin1, flash_temp2, side='right')
    x_flash2 = (flash1_hf_table[qq](flash_temp1)-flash2_hf_table[qq2](flash_temp2))/\
                (flash2_hg_table[qq2](flash_temp2)-flash2_hf_table[qq2](flash_temp2)) # D92
    mass_steam_flash2 = x_flash2*1000.0*(1-x_flash1) # D94/L92

    hpTurbine_mSteam = mass_steam_flash1-totSteamFlow # lb/hour L77
    lpTurbine_mSteam = mass_steam_flash2 # lb/hour L92

    qq = np.searchsorted(bin2, pHP, side='right')
    u90_temp = u90_table[qq](pHP)
    qq = np.searchsorted(bin1, u90_temp, side='right')
    hpTurbine_sTin = turb_perf_sTin_table[qq](u90_temp) # L68
    qq = np.searchsorted(bin1, tCond, side='right')
    hpTurbine_sCondVap = turb_perf_sCondVap_table[qq](tCond) # L69
    hpTurbine_sCondLiq = turb_perf_sCondLiq_table[qq](tCond) # L70
    hpTurbine_hExIsent = h_condF + (h_condG-h_condF)*\
        (hpTurbine_sTin-hpTurbine_sCondLiq)/(hpTurbine_sCondVap-hpTurbine_sCondLiq) # L71
    qq = np.searchsorted(bin1, u90_temp, side='right')
    hpTurbine_dhIsent = turb_perf_hTin_table[qq](u90_temp)-hpTurbine_hExIsent # L72
    hpTurbine_A = 0.5*turbine_eff*hpTurbine_dhIsent # L73
    hpTurbine_hexPrime =(turb_perf_hTin_table[qq](u90_temp)-hpTurbine_A*
                         (1-h_condF/(h_condG-h_condF)))/(1+hpTurbine_A/(h_condG-h_condF)) # L74

    qq = np.searchsorted(bin2, pLP, side='right')
    u93_temp = u93_table[qq](pLP)
    qq = np.searchsorted(bin1, u93_temp, side='right')
    lpTurbine_sTin = lpTurb_sTin_table[qq](u93_temp) # L85
    lpTurbine_hExIsent = h_condF+(h_condG-h_condF)*(lpTurbine_sTin-
                        hpTurbine_sCondLiq)/(hpTurbine_sCondVap-hpTurbine_sCondLiq) # L86
    lpTurbine_dhIsent = lpTurb_hTin_table[qq](u93_temp)-lpTurbine_hExIsent # L87
    lpTurbine_A = 0.5*turbine_eff*lpTurbine_dhIsent # L88
    lpTurbine_hexPrime = (lpTurb_hTin_table[qq](u93_temp)-lpTurbine_A*
                          (1-h_condF/(h_condG-h_condF)))/(1+lpTurbine_A/(h_condG-h_condF)) # L89

    # Heat Rejected from Tower
    mSteam_rejected = hpTurbine_mSteam+lpTurbine_mSteam # lb/hour D97
    hCondIn_rejected = (hpTurbine_mSteam*hpTurbine_hexPrime+lpTurbine_mSteam*
                        lpTurbine_hexPrime)/(hpTurbine_mSteam+lpTurbine_mSteam) # D98
    q_Cond = mSteam_rejected*(hCondIn_rejected-h_condF) # D100
    qRejectedTower = q_Cond+qRejected_NCG_system # btu/hr  D102
    cwFlow = qRejectedTower/cooling_water_temp_rise # lb/hour D115
    drift = 0.001*cwFlow # lb/hour D131

    # Evaporative Water loss   
    a = -0.0001769*np.log(cooling_water_temp_rise) + 0.0011083
    b = 0.0657628*np.log(cooling_water_temp_rise) - 0.4091309
    c = -6.7041142*np.log(cooling_water_temp_rise) + 44.3438937
    d = -0.0325112*cooling_water_temp_rise**2 + 6.831236*cooling_water_temp_rise - 64.6250943 
    evaporativeWaterLoss = (a*t_wetbulb**3+b*t_wetbulb**2+c*t_wetbulb+d)*qRejectedTower/1000000 # D130    
    loss_in_NCG_removal = steam_leaving_stage2*mw_H2O # lb/hour D132
    fraction_gf_flow_not_injected = (evaporativeWaterLoss+drift+loss_in_NCG_removal)/1000.0 # D136
    frac_inlet = 1-fraction_gf_flow_not_injected     

''' Total Pumping Power '''

rho_amb = 1/(steam_properties[3][0]*t_amb**6 + steam_properties[3][1]*t_amb**5 + steam_properties[3][2]*t_amb**4\
            + steam_properties[3][3]*t_amb**3 + steam_properties[3][4]*t_amb**2 + steam_properties[3][5]*t_amb + steam_properties[3][6]) # lb/cft
g = 32.174 # ft/s^2  G17
production_depth = cur_dict['well_depth']*3.28084 # G13/G15
dT = (((tSource-t_amb)-32)/1.8)/(cur_dict['well_depth'])# F/ft  
t_correction = 30.796*t_c**-0.552

p_hydrostaticProd = ((p_amb/14.50377)+(1/0.000000000464)*(np.exp(((rho_amb*16.01846)*
                    9.807*0.000000000464*(cur_dict['well_depth']-0.5*(0.0009/t_correction)*
                    dT*cur_dict['well_depth']**2))/100000)-1))*14.50377# G14/G16/G195
flow_per_Production_Well = 873029.52 # lb/hour G19/G347
productivityIndex_ProdWell = 2500. # lb/hour-psi G20/G26
drawdown = flow_per_Production_Well/productivityIndex_ProdWell # psi G21

# flow_ration_Prod2Inj_Well = 0.75 # G23
flow_per_InjWell = flow_per_Production_Well/0.75 # lb/hour G25/G65
buildup = flow_per_InjWell/productivityIndex_ProdWell # psi G27

# Production well head pressure/pump suction
prodWellHead_Temp = tSource-(1.23/production_depth)*production_depth # F G30/G198
prodWellHead_Psat = steam_properties[0][0]*prodWellHead_Temp**6 + steam_properties[0][1]*prodWellHead_Temp**5 +\
                    steam_properties[0][2]*prodWellHead_Temp**4 + steam_properties[0][3]*prodWellHead_Temp**3 + \
                    steam_properties[0][4]*prodWellHead_Temp**2 + steam_properties[0][5]*prodWellHead_Temp + steam_properties[0][6]# psia G32/G200
prodWellHead_excess_pres = 50.0 # psia G33
prodWellHead_Pres = prodWellHead_Psat+prodWellHead_excess_pres # psia G34

# Injection well head pressure/pump suction
injWell_dP_Binary = 40. # psid G37
injWell_Pflash = 15.71 # psia G38
if plant_type == "binary":
    injWell_P = prodWellHead_Pres-injWell_dP_Binary
    tR_bin = 1+(-0.002954*t_c - 0.121503)*(availableEnergy/(((hSource-hSink)-(tSink+460)*(sSource-sSink))/3.413))
    injWell_T = max(((t_c+273)*tR_bin)-273, ((tAmphSiO2*1.8+32 +1)-32)/1.8)*1.8+32 # F G40
else:
    injWell_P = injWell_Pflash # psia G39
    injWell_T = max(tAmphSiO2, (0.897*(t_c+273)-57)-273)*1.8+32 # F G40  OUT!G701  ->  Geofluid H96
injWell_rho = 1/(steam_properties[3][0]*injWell_T**6 + steam_properties[3][1]*injWell_T**5 + steam_properties[3][2]*injWell_T**4\
             + steam_properties[3][3]*injWell_T**3 + steam_properties[3][4]*injWell_T**2 + steam_properties[3][5]*injWell_T + steam_properties[3][6]) # lb/cft G42

# Friction injection casing
diameter = 12.5/12 # ft G68/G206
fric_casing_A = np.pi*diameter**2/4 # sqft G69
length_interval = 3937.0 # ft G70 !The lookup table for this is unreachable which means we can't change well depth
casing_surf_roughness = 0.00015 # ft G71
injWell_dT = 1.13/length_interval # F/ft G45
fric_casing_T = injWell_T+0.5*injWell_dT*0.75*length_interval # F G72
fric_casing_Psat = steam_properties[0][0]*fric_casing_T**6 + steam_properties[0][1]*fric_casing_T**5 +\
                   steam_properties[0][2]*fric_casing_T**4 + steam_properties[0][3]*fric_casing_T**3 + \
                   steam_properties[0][4]*fric_casing_T**2 + steam_properties[0][5]*fric_casing_T + steam_properties[0][6]# psia G73
fric_casing_visc = 407.22*fric_casing_T**-1.194/3600 # lb/ft-s G74
fric_casing_rho = 1/(steam_properties[3][0]*fric_casing_T**6 + steam_properties[3][1]*fric_casing_T**5 + steam_properties[3][2]*fric_casing_T**4\
             + steam_properties[3][3]*fric_casing_T**3 + steam_properties[3][4]*fric_casing_T**2 + steam_properties[3][5]*fric_casing_T + steam_properties[3][6]) # lb/cft G75
rhoXhead = length_interval*fric_casing_rho/144 # psi G76
pstar_psat = 0.5*(rhoXhead+injWell_P)/fric_casing_Psat # G77
fric_casing_flow = flow_per_InjWell/(3600*fric_casing_rho) # cft/s G80
fric_casing_vel = fric_casing_flow/fric_casing_A # fps G81
fric_casing_rho*= (1+(7.15037E-19*fric_casing_T**5.91303)*(pstar_psat-1))
fric_casing_visc*= (1+(4.02401E-18*fric_casing_T**5.736882)*(pstar_psat-1))
fric_casing_Rew = diameter*fric_casing_vel*fric_casing_rho/fric_casing_visc # G82
serghide_a = -2*np.log10((casing_surf_roughness/diameter)/3.7+12/fric_casing_Rew) # G84
serghide_v = -2*np.log10((casing_surf_roughness/diameter)/3.7+2.51*serghide_a/fric_casing_Rew) # G85
serghide_c = -2*np.log10((casing_surf_roughness/diameter)/3.7+2.51*serghide_v/fric_casing_Rew) # G86
friction_factor = (serghide_a-(serghide_v-serghide_a)**2/(serghide_c-2*serghide_v+serghide_a))**-2 # G87/G88
friction_head_loss = (((friction_factor*1/diameter)*(fric_casing_vel**2)/(2*g))*length_interval)*fric_casing_rho/144 # bar G91
pBottomUpperInterval = injWell_P+fric_casing_rho*length_interval/144-friction_head_loss # psi G92

# injection pump
intermediate_int2_T = injWell_T+injWell_dT*0.75*(length_interval+0.5*0) # F G99
intermediate_int2_psat = steam_properties[0][0]*intermediate_int2_T**6 + steam_properties[0][1]*intermediate_int2_T**5 +\
                         steam_properties[0][2]*intermediate_int2_T**4 + steam_properties[0][3]*intermediate_int2_T**3 + \
                         steam_properties[0][4]*intermediate_int2_T**2 + steam_properties[0][5]*intermediate_int2_T + steam_properties[0][6] # G100
intermediate_int2_rho = 1/(steam_properties[3][0]*intermediate_int2_T**6 + steam_properties[3][1]*intermediate_int2_T**5 +\
                           steam_properties[3][2]*intermediate_int2_T**4 + steam_properties[3][3]*intermediate_int2_T**3 + \
                           steam_properties[3][4]*intermediate_int2_T**2 + steam_properties[3][5]*intermediate_int2_T + \
                           steam_properties[3][6]) # lb/cft G102
intermediate_int2_rho_corr = 1+(7.15037E-19*intermediate_int2_T**5.91303)*\
                             (((pBottomUpperInterval+0.5*0)/intermediate_int2_psat)-1)# G105
bottom_intermediate_intervalP2 = pBottomUpperInterval+0*intermediate_int2_rho*\
                                 intermediate_int2_rho_corr/144 # G119

# Injection Interval 
bottom_intermediate_intervalP3 = bottom_intermediate_intervalP2 # G146
injInterval_D = 12.25/12 # ft G149
injInterval_T = injWell_T+injWell_dT*0.75*(length_interval+0+0+0.5*984.2) # F G153
injInterval_psat = steam_properties[0][0]*injInterval_T**6 + steam_properties[0][1]*injInterval_T**5 +\
                   steam_properties[0][2]*injInterval_T**4 + steam_properties[0][3]*injInterval_T**3 + \
                   steam_properties[0][4]*injInterval_T**2 + steam_properties[0][5]*injInterval_T + steam_properties[0][6] # psia G154
injInterval_visc = 407.22*injInterval_T**-1.194/3600 # lb/ft-s G155
injInterval_rho = 1/(steam_properties[3][0]*injInterval_T**6 + steam_properties[3][1]*injInterval_T**5 +\
                     steam_properties[3][2]*injInterval_T**4 + steam_properties[3][3]*injInterval_T**3 + \
                     steam_properties[3][4]*injInterval_T**2 + steam_properties[3][5]*injInterval_T + \
                     steam_properties[3][6]) # lb/cft G156
injInterval_rho_head = 984.2*injInterval_rho/144 # G157
injInterval_P_Psat = (bottom_intermediate_intervalP3+0.5*injInterval_rho_head)/injInterval_psat # psia G158
injInterval_rho_corr = 1+(7.15037E-19*injInterval_T**5.91303)*(
                        ((bottom_intermediate_intervalP3+0.5*(984.2*injInterval_rho/144))/
                         injInterval_psat)-1) # G159
injInterval_mu_corr = 1+(4.02401E-18*injInterval_T**5.736882)*(injInterval_P_Psat-1) # G160
injInterval_flow = flow_per_InjWell/(3600*injInterval_rho*injInterval_rho_corr) # fps G161
injInterval_velo = injInterval_flow/(np.pi*injInterval_D**2/4) # fps G162
injInterval_Rew = injInterval_D*injInterval_velo*injInterval_rho*injInterval_rho_corr/(injInterval_visc*injInterval_mu_corr) # G163
serghide_a = -2*np.log10((0.02/injInterval_D)/3.7+12/injInterval_Rew) # G165
serghide_v = -2*np.log10((0.02/injInterval_D)/3.7+2.51*serghide_a/injInterval_Rew) # G166
serghide_c = -2*np.log10((0.02/injInterval_D)/3.7+2.51*serghide_v/injInterval_Rew) # G167
injInterval_friction_factor = (serghide_a-(serghide_v-serghide_a)**2/(serghide_c-2*serghide_v+serghide_a))**-2 # G169
injInterval_friction_head_loss = (((injInterval_friction_factor*1/injInterval_D)*\
                                   (injInterval_velo**2)/(2*g))*984.2)*injInterval_rho*injInterval_rho_corr/144 # psi G172
bottom_hole_pressure_injection_well = bottom_intermediate_intervalP3+injInterval_rho*injInterval_rho_corr*984.2/144-injInterval_friction_head_loss # psi G174
injectPump_excessP = bottom_hole_pressure_injection_well-p_hydrostaticProd # psi G180
injection_pump_head_used = -injectPump_excessP+buildup+1 # G183


# Production Pumping Bottom of Hole
prodWell_BottomP = p_hydrostaticProd-drawdown # psia G196
prodWell_BottomDP = prodWell_BottomP-(prodWellHead_Psat+prodWellHead_excess_pres) # psia G202
t_avg = 320-0.5*injWell_dT*984.2 # F G211
below_pump_rho = 1/(steam_properties[3][0]*t_avg**6 + steam_properties[3][1]*t_avg**5 + steam_properties[3][2]*t_avg**4\
             + steam_properties[3][3]*t_avg**3 + steam_properties[3][4]*t_avg**2 + steam_properties[3][5]*t_avg + steam_properties[3][6]) # lb/cft G213
below_pump_Psat = steam_properties[0][0]*t_avg**6 + steam_properties[0][1]*t_avg**5 +\
                  steam_properties[0][2]*t_avg**4 + steam_properties[0][3]*t_avg**3 + \
                  steam_properties[0][4]*t_avg**2 + steam_properties[0][5]*t_avg + steam_properties[0][6] # psia G212
below_pump_P_Psat = (prodWell_BottomP-0.5*below_pump_rho*984.2/144)/below_pump_Psat # G214
below_pump_rho *= 1+(7.15037E-19*t_avg**5.91303)*(below_pump_P_Psat-1) # lb/cft G217

blow_pump_diameter = 12.25/12 # ft G206
below_pump_flow = (flow_per_Production_Well/below_pump_rho)/3600 # cft/s G218
below_pump_vel = below_pump_flow/(np.pi*blow_pump_diameter**2/4) # fps G219
below_pump_visc = (1+(4.02401E-18*t_avg**5.736882)*(below_pump_P_Psat-1))*407.22*t_avg**-1.194/3600 # fps G220
below_pump_Rew = blow_pump_diameter*below_pump_vel*below_pump_rho/below_pump_visc # G221
serghide_a = -2*np.log10((0.02/blow_pump_diameter)/3.7+12/below_pump_Rew) # G223
serghide_v = -2*np.log10((0.02/blow_pump_diameter)/3.7+2.51*serghide_a/below_pump_Rew) # G224
serghide_c = -2*np.log10((0.02/blow_pump_diameter)/3.7+2.51*serghide_v/below_pump_Rew) # G225
below_friction_factor = (serghide_a-(serghide_v-serghide_a)**2/(serghide_c-2*serghide_v+serghide_a))**-2 # G227
below_friction_head_loss = (((below_friction_factor*1/blow_pump_diameter)*(below_pump_vel**2)/(2*g))*984.2)*below_pump_rho/144 # bar G230
pTopProdZone = prodWell_BottomP-below_friction_head_loss-below_pump_rho*984.2/144 # psi G231

# Upper Interval Below Pump # 288
uibp_diameter = 12.5/12 # ft G290
uibp_T = tSource-984.2*injWell_dT # F G292
distance2surface = production_depth-984.2 # ft G293
uibp_tavg = uibp_T-0.5*injWell_dT*distance2surface # G295
uibp_Psat = steam_properties[0][0]*uibp_tavg**6 + steam_properties[0][1]*uibp_tavg**5 +\
            steam_properties[0][2]*uibp_tavg**4 + steam_properties[0][3]*uibp_tavg**3 + \
            steam_properties[0][4]*uibp_tavg**2 + steam_properties[0][5]*uibp_tavg + steam_properties[0][6] # psia G296
uibp_rho = 1/(steam_properties[3][0]*uibp_tavg**6 + steam_properties[3][1]*uibp_tavg**5 + steam_properties[3][2]*uibp_tavg**4\
           + steam_properties[3][3]*uibp_tavg**3 + steam_properties[3][4]*uibp_tavg**2 + steam_properties[3][5]*uibp_tavg + steam_properties[3][6]) # lb/cft G297 
uibp_P_Psat = 0.5*(uibp_Psat+pTopProdZone)/uibp_Psat # G298
uibp_rho *= 1+(7.15037E-19*uibp_tavg**5.91303)*(uibp_P_Psat-1) # G301
uibp_flowRate = (flow_per_Production_Well/uibp_rho)/3600 # cfs G302
uibp_velo = uibp_flowRate/(np.pi*uibp_diameter**2/4) # lb/ft-s G303
uibp_visc = (1+(4.02401E-18*uibp_tavg**5.736882)*(uibp_P_Psat-1))* 407.22*uibp_tavg**-1.194/3600 # lb/ft-s G304 
uibp_Rew = uibp_diameter*uibp_velo*uibp_rho/uibp_visc # G305
serghide_a = -2*np.log10((0.00015/uibp_diameter)/3.7+12/uibp_Rew) # G307
serghide_v = -2*np.log10((0.00015/uibp_diameter)/3.7+2.51*serghide_a/uibp_Rew) # G308
serghide_c = -2*np.log10((0.00015/uibp_diameter)/3.7+2.51*serghide_v/uibp_Rew) # G309
uibp_friction_factor = (serghide_a-(serghide_v-serghide_a)**2/(serghide_c-2*serghide_v+serghide_a))**-2 # G311
uibp_friction_head_loss = (uibp_friction_factor*1/uibp_diameter)*(uibp_velo**2)/(2*g) # ft/ft G312

# Friction production casing 346
prodCas_Diameter = 8.681/12 # ft G349
prodCas_velo = uibp_flowRate/(np.pi*prodCas_Diameter**2/4) # fps G351
prodCas_Rec = prodCas_velo*prodCas_Diameter*uibp_rho/uibp_visc # G354
serghide_a = -2*np.log10((0.00015/prodCas_Diameter)/3.7+12/prodCas_Rec) # G357
serghide_v = -2*np.log10((0.00015/prodCas_Diameter)/3.7+2.51*serghide_a/prodCas_Rec) # G358
serghide_c = -2*np.log10((0.00015/prodCas_Diameter)/3.7+2.51*serghide_v/prodCas_Rec) # G359
prodCas_friction_factor = (serghide_a-(serghide_v-serghide_a)**2/(serghide_c-2*serghide_v+serghide_a))**-2 # G361
prodCas_friction_head_loss = (prodCas_friction_factor*949.1351129/prodCas_Diameter)*(prodCas_velo**2)/(2*g)# ft G363

if plant_type == "binary":
    pumpLift = 949.1351129+prodCas_friction_head_loss # ft G368
    production_flow_rate = 794119.5420575 # lb/hr  https://publications.mygeoenergynow.org/grc/1033909.pdf
else:
    pumpLift = 0
    production_flow_rate = 439651.42196051 # lb/hr  https://publications.mygeoenergynow.org/grc/1033909.pdf

production_pump_eff = 0.675 # G372
production_pump_work = ((((pumpLift*1)/(60*33000))/production_pump_eff)*0.7457)*1000 # w-h/lb G386

total_pump_power = production_pump_work + (((((((injection_pump_head_used*144/injWell_rho)*1.)
                                               /(60*33000))/production_pump_eff)*0.7457)*1000)*frac_inlet)  # w-h/lb

''' Total Plant Power '''
brine_effectiveness = availableEnergy*second_law_eff
total_power = (brine_effectiveness*production_flow_rate/1000)-(total_pump_power*production_flow_rate/1000) # kW

