import json
import pathlib
import os
import pandas as pd
import struct
import datetime
import numpy as np
import multiprocessing

PROJECT_PATH = pathlib.Path(__file__).parent
DATA_PATH = PROJECT_PATH / "data"

class solar_pv:
    """
    A class representing a solar photovoltaic system.

    Parameters:
    - lat (float): Latitude of the location (in degrees).
    - lon (float): Longitude of the location (in degrees).
    - rooftop_area (float): Area of the rooftop available for the solar PV system (in square feet).
    - system_loss (float, optional): System loss as a percentage (default is 14.08%).
    - dc_ac (float, optional): DC to AC ratio (default is 1.2).
    - invert_eff (float, optional): Inverter efficiency as a percentage (default is 96%).
    - per_area (float, optional): Percentage of the rooftop area covered by solar panels (default is 25%).
    - tilt (float, optional): Tilt angle of the solar panels (in degrees, default is 20°).
    - azimuth (float, optional): Azimuth angle of the solar panels (in degrees, default is 180°).
    - mod_type (int, optional): Type of solar module (0 for Standard, 1 for Premium, 2 for Thin Film, default is 0).

    Attributes:
    - lat (float): Latitude of the location.
    - lon (float): Longitude of the location.
    - system_loss (float): System loss as a decimal.
    - dc_ac (float): DC to AC ratio.
    - invert_eff (float): Inverter efficiency as a decimal.
    - tilt (float): Tilt angle of the solar panels.
    - azimuth (float): Azimuth angle of the solar panels.
    - dc_nameplate (float): DC system nameplate capacity (in kW).
    - solar_data (DataFrame): Solar data for the location.
    - metadata (dict): Metadata about the location and solar data.
    - tz (int): Time zone of the location.

    Methods:
    - read_bin_file(id): Read solar data from a binary file for the specified location ID.

    Note: The class initializes by reading solar data from a binary file based on the given latitude and longitude.

    Example usage:
    >>> pv_system = solar_pv(lat=37.7749, lon=-122.4194, rooftop_area=100, tilt=30, azimuth=180)
    >>> print(pv_system.solar_data.head())
    >>> print(pv_system.metadata)
    """
    def __init__(self, lat:float, lon:float, rooftop_area:float, system_loss=14.08, dc_ac=1.2, 
                 invert_eff=96, per_area=25, tilt=20, azimuth=180, mod_type=0, ground=False) -> None:
        # Latitude and Longitude
        self.lat = lat
        self.lon = lon
        lat = str(lat)
        lon = str(lon)
        # Look up what the station ID for the Lat/Lon
        # key = f'{lat[:lat.find(".")+3]}&{lon[:lon.find(".")+3]}' GSA Data only
        key = f'{lat[:lat.find(".")]}&{lon[:lon.find(".")]}'
        with open(os.path.join(DATA_PATH, "Lat_Lon_keyMap.json")) as json_data:
            map_dict = json.load(json_data)
        if key in map_dict.keys():
            station_id = map_dict[key][0]
        else:
            station_id = self.find_closest_station(map_dict)
        # Other system attributes
        self.system_loss = system_loss/100
        self.dc_ac = dc_ac
        self.invert_eff = invert_eff/100
        self.tilt = tilt 
        self.azimuth = azimuth
        if mod_type == 0:
            self.nom_eff = 0.16 # Standard
            self.gamma = -0.0047
        elif mod_type == 1:
            self.nom_eff = 0.18 # Premium
            self.gamma = -0.0035
        else:
            self.nom_eff = 0.11 # Thin Film
            self.gamma = -0.0020
        # DC System Size
        # Size (kW) = Array Area (m²) × 1 kW/m² × Module Efficiency (%)
        self.dc_nameplate = ((rooftop_area*0.092903)*(per_area/100)) * 1 * self.nom_eff
        # installed nominal operating temperature 
        if ground:
            self.inoct = 45
        else:
            self.inoct = 50
        # Read the solar data from the station binary
        self.solar_data = None
        self.metadata = None
        self.solar_power = None
        self.tz = None # Time zone
        self.monthly_ac = np.zeros(12)
        self.monthly_rad = np.zeros(12)
        self.total = 0
        self.read_bin_file(station_id)

    def read_bin_file(self, id:int):
        """
        Read solar data from a binary file for the specified location ID.

        Parameters:
        - id (int): Location ID used to retrieve the corresponding binary file.

        This method reads solar data from a binary file and populates the 'solar_data' and 'metadata' attributes
        of the solar_pv instance.
        """
        metadata = {}

        with open(os.path.join(DATA_PATH, 'solar_data', f"{id}.bin"), "rb") as f:
            bin_data = f.read()

        # Write the metadata
        data_pair = [x.decode("utf-8") for x in struct.unpack("6s5s", bin_data[:11])]
        metadata[data_pair[0]] = data_pair[1] # Source
        data_pair = [x.decode("utf-8") for x in struct.unpack("11s10s", bin_data[11:32])]
        metadata[data_pair[0]] = data_pair[1].rstrip("\x00") # Location ID
        data_pair = [x.decode("utf-8") for x in struct.unpack("4s30s", bin_data[32:66])]
        metadata[data_pair[0]] = data_pair[1].rstrip("\x00") # City
        data_pair = [x.decode("utf-8") for x in struct.unpack("5s15s", bin_data[66:86])]
        metadata[data_pair[0]] = data_pair[1].rstrip("\x00") # State
        data_pair = [x.decode("utf-8") for x in struct.unpack("7s20s", bin_data[86:113])]
        metadata[data_pair[0]] = data_pair[1].rstrip("\x00") # Country
        data_pair = struct.unpack("8sd", bin_data[113:129])
        metadata[data_pair[0].decode("utf-8")] = data_pair[1] # Latitude
        data_pair = struct.unpack("9sd", bin_data[129:153])
        metadata[data_pair[0].decode("utf-8")] = data_pair[1] # Longitude
        data_pair = struct.unpack("9si", bin_data[153:169])
        metadata[data_pair[0].decode("utf-8")] = data_pair[1] # Time Zone
        data_pair = struct.unpack("9sd", bin_data[169:193])
        metadata[data_pair[0].decode("utf-8")] = data_pair[1] # Elevation
        data_pair = struct.unpack("15si", bin_data[193:213])
        metadata[data_pair[0].decode("utf-8")] = data_pair[1] # Local Time Zone

        # Write data to dictionary for conversion to pd.DataFrame
        tz = datetime.timezone(datetime.timedelta(hours=metadata['Time Zone']))
        self.tz = metadata['Time Zone']
        headers = [x.decode("utf-8") for x in struct.unpack("4s5s3s4s6s3s3s3s11s", bin_data[213:255])]
        data_dict = {k: [] for k in headers}
        data_dict['datetime'] = []

        index = 255
        for _ in range(8760):
            temp = struct.unpack("9d", bin_data[index:index+72])
            for qq, key in enumerate(headers):
                if qq <= 4:
                    data_dict[key].append(int(temp[qq]))
                else:
                    data_dict[key].append(temp[qq])
            data_dict['datetime'].append(datetime.datetime(int(temp[0]), int(temp[1]), 
                                        int(temp[2]), int(temp[3]), int(temp[4]), tzinfo=tz))
            index += 72

        data = pd.DataFrame.from_dict(data_dict)
        data = data.set_index('datetime')
        
        self.solar_data = data
        self.metadata = metadata

    def find_closest_station(self, map_dict):
        min_dist = 1.0e20
        station_id = ''
        for key, value in map_dict.items():
            split = key.split('&')
            lat = float(split[0])
            lon = float(split[1])
            dist = np.arccos( np.sin(np.radians(self.lat))*np.sin(np.radians(lat)) + \
                             np.cos(np.radians(self.lat))*np.cos(np.radians(lat))*\
                                np.cos(np.radians(lon)-np.radians(self.lon)) ) * 6371000
            if dist < min_dist:
                station_id = value[0]
                min_dist = dist
        return station_id

    def calculate_solar_position(self):
        '''
        Calculates the position of the sun at every hour for an entire year
        These formulas are taken from NOAA's solar calculator 
        https://gml.noaa.gov/grad/solcalc/
        The calculations in the NOAA Sunrise/Sunset and Solar Position Calculators are based on equations from 
        Astronomical Algorithms, by Jean Meeus. The sunrise and sunset results are theoretically accurate to 
        within a minute for locations between +/- 72° latitude, and within 10 minutes outside of those latitudes. 
        However, due to variations in atmospheric composition, temperature, pressure and conditions, observed 
        values may vary from calculations.
        '''   
        self.solar_data['Julian Day'] = pd.DatetimeIndex(self.solar_data.index).to_julian_date()
        self.solar_data['Julian Century'] = (self.solar_data['Julian Day']-2451545)/36525
        self.solar_data['Geom Mean Long Sun (deg)'] = (280.46646+self.solar_data['Julian Century']*
                                                       (36000.76983 + self.solar_data['Julian Century']*0.0003032))%360.
        self.solar_data['Geom Mean Anom Sun (deg)'] = 357.52911+self.solar_data['Julian Century']*\
                                                        (35999.05029 - 0.0001537*self.solar_data['Julian Century'])
        self.solar_data['Eccent Earth Orbit'] = 0.016708634-self.solar_data['Julian Century']*\
                                                (0.000042037+0.0000001267*self.solar_data['Julian Century'])
        self.solar_data['Sun Eq of Ctr'] = np.sin(np.radians(self.solar_data['Geom Mean Anom Sun (deg)']))*\
                                            (1.914602-self.solar_data['Julian Century']*(0.004817+0.000014*self.solar_data['Julian Century']))+\
                                            np.sin(np.radians(2*self.solar_data['Geom Mean Anom Sun (deg)']))*(0.019993-0.000101*self.solar_data['Julian Century'])+\
                                            np.sin(np.radians(3*self.solar_data['Geom Mean Anom Sun (deg)']))*0.000289
        self.solar_data['Sun True Long (deg)'] = self.solar_data['Geom Mean Long Sun (deg)'] + self.solar_data['Sun Eq of Ctr']
        self.solar_data['Sun True Anom (deg)'] = self.solar_data['Geom Mean Anom Sun (deg)'] + self.solar_data['Sun Eq of Ctr']
        self.solar_data['Sun Rad Vector (AUs)'] = (1.000001018*(1-self.solar_data['Eccent Earth Orbit']**2))/\
                                            (1+self.solar_data['Eccent Earth Orbit']*np.cos(np.radians(self.solar_data['Sun True Anom (deg)'])))
        self.solar_data['Sun App Long (deg)'] = self.solar_data['Sun True Long (deg)']-0.00569-0.00478*\
                                            np.sin(np.radians(125.04-1934.136*self.solar_data['Julian Century']))
        self.solar_data['Mean Obliq Ecliptic (deg)'] = 23+(26+((21.448-self.solar_data['Julian Century']*(46.815+
                                            self.solar_data['Julian Century']*(0.00059-self.solar_data['Julian Century']*0.001813))))/60)/60
        self.solar_data['Obliq Corr (deg)'] = self.solar_data['Mean Obliq Ecliptic (deg)']+0.00256*\
                                            np.cos(np.radians(125.04-1934.136*self.solar_data['Julian Day'])) 
        self.solar_data['Sun Rt Ascen (deg)'] = np.rad2deg(np.arctan2(np.cos(np.radians(self.solar_data['Sun App Long (deg)'])),
                                            np.cos(np.radians(self.solar_data['Obliq Corr (deg)']))*np.sin(np.radians(self.solar_data['Sun App Long (deg)']))))
        self.solar_data['Sun Declin (deg)'] = np.rad2deg(np.arcsin(np.sin(np.radians(self.solar_data['Obliq Corr (deg)']))*\
                                                                  np.sin(np.radians(self.solar_data['Sun App Long (deg)']))))
        self.solar_data['var y'] = np.tan(np.radians(self.solar_data['Obliq Corr (deg)']/2))*\
                                            np.tan(np.radians(self.solar_data['Obliq Corr (deg)']/2))
        self.solar_data['Eq of Time (minutes)']	= 4*np.rad2deg(self.solar_data['var y']*np.sin(2*np.radians(self.solar_data['Geom Mean Long Sun (deg)']))\
                                            -2*self.solar_data['Eccent Earth Orbit']*np.sin(np.radians(self.solar_data['Geom Mean Anom Sun (deg)']))+4*\
                                            self.solar_data['Eccent Earth Orbit']*self.solar_data['var y']*np.sin(np.radians(self.solar_data['Geom Mean Anom Sun (deg)']))\
                                            *np.cos(2*np.radians(self.solar_data['Geom Mean Long Sun (deg)']))-0.5*self.solar_data['var y']**2*\
                                            np.sin(4*np.radians(self.solar_data['Geom Mean Long Sun (deg)']))-1.25*self.solar_data['Eccent Earth Orbit']**2*\
                                            np.sin(2*np.radians(self.solar_data['Geom Mean Anom Sun (deg)'])))
        self.solar_data['HA Sunrise (deg)']	= np.rad2deg(np.arccos(np.cos(np.radians(90.833))/(np.cos(np.radians(self.lat))*
                                            np.cos(np.radians(self.solar_data['Sun Declin (deg)'])))-np.tan(np.radians(self.lat))*
                                            np.tan(np.radians(self.solar_data['Sun Declin (deg)']))))
        self.solar_data['Solar Noon (LST)']	= (720-4*self.lon-self.solar_data['Eq of Time (minutes)']+self.metadata['Local Time Zone']*60)/1440
        self.solar_data['Sunrise Time (LST)'] = (self.solar_data['Solar Noon (LST)']*1440-self.solar_data['HA Sunrise (deg)']*4)/1440
        self.solar_data['Sunset Time (LST)'] = (self.solar_data['Solar Noon (LST)']*1440+self.solar_data['HA Sunrise (deg)']*4)/1440
        self.solar_data['Sunlight Duration (minutes)'] = 8*self.solar_data['HA Sunrise (deg)']
        self.solar_data['True Solar Time (min)'] = ((self.solar_data.index.hour+0.5)/24*1440+self.solar_data['Eq of Time (minutes)']+4*self.lon-60*self.metadata['Local Time Zone'])%1440        
        if np.sum(self.solar_data['True Solar Time (min)'] < 0) > 0:
            dmy = np.zeros(self.solar_data.shape[0]) 
            dmy[self.solar_data['True Solar Time (min)'] < 0] = self.solar_data['True Solar Time (min)']/4 + 180
            dmy[self.solar_data['True Solar Time (min)'] > 0] = self.solar_data['True Solar Time (min)']/4 - 180
            self.solar_data['Hour Angle (deg)'] = dmy.tolist()	
        else:
            self.solar_data['Hour Angle (deg)'] = self.solar_data['True Solar Time (min)']/4 - 180        
        self.solar_data['Solar Zenith Angle (deg)'] = np.rad2deg(np.arccos(np.sin(np.radians(self.lat))*np.sin(np.radians(self.solar_data['Sun Declin (deg)']))\
                                                        +np.cos(np.radians(self.lat))*np.cos(np.radians(self.solar_data['Sun Declin (deg)']))*\
                                                        np.cos(np.radians(self.solar_data['Hour Angle (deg)']))))

    def calculate_solar_power(self):
        '''

        Citations
        [1] Kalogirou, S. A. (2014). Solar energy engineering : processes and systems. Academic Press.
        [2] “Perez Sky Diffuse Model” GitHub, NREL, 10 June 1998, https://github.com/NREL/ssc/blob/6cddbd76efe4ab20f1670718f82005a5e306cd2e/shared/lib_irradproc.cpp#L1575. Accessed 25 Sept. 2023.
        [3] “PVWatts Power Out” GitHub, NREL, 10 June 1998, https://github.com/NREL/ssc/blob/6cddbd76efe4ab20f1670718f82005a5e306cd2e/ssc/cmod_pvwattsv5.cpp#L327. Accessed 25 Sept. 2023.        
        '''
        # Incident Angle
        # [1] equation 2.18 pg 60
        # cos(theta) = sin(L)sin(delta)cos(beta) - cos(L)sin(delta)sin(beta)cos(Zs)
        #               + cos(L)cos(delta)cos(h)cos(beta)
        #               + sin(L)cos(delta)cos(h)sin(beta)cos(Zs)
        #               + cos(delta)sin(h)sin(beta)sin(Zs)
        # Where L = Latitude, Delta = solar declination, Beta = PV array tilt angle from horizontal
        # h = hour angle, Zs = surface azimuth angle
        inc_angle = np.zeros(self.solar_data.shape[0])
        L = np.radians(self.lat)
        delta = np.radians(self.solar_data['Sun Declin (deg)'].values)   
        beta = np.radians(self.tilt) 
        h = np.radians(self.solar_data['Hour Angle (deg)'])
        Zs = np.radians(180.-self.azimuth)    
        inc_angle = np.rad2deg(np.arccos(np.sin(L)*np.sin(delta)*np.cos(beta) -
                                         np.cos(L)*np.sin(delta)*np.sin(beta)*np.cos(Zs) +
                                         np.cos(L)*np.cos(delta)*np.cos(h)*np.cos(beta) +
                                         np.sin(L)*np.cos(delta)*np.cos(h)*np.sin(beta)*np.cos(Zs) +
                                         np.cos(delta)*np.sin(h)*np.sin(beta)*np.sin(Zs)))
        self.solar_power = pd.DataFrame({'Incident Angle': inc_angle})
        
        # Perez Sky Diffuse calculation [2]       
        radiation = np.zeros((self.solar_data.shape[0],3)) 
        for qq in range(self.solar_data.shape[0]):
            rad = self.perez(qq)
            radiation[qq,0] = rad[0]
            radiation[qq,1] = rad[1]
            radiation[qq,2] = rad[2]   
        self.solar_power['Beam'] = radiation[:,0].tolist()
        self.solar_power['Total Sky Diffuse'] = radiation[:,1].tolist()
        self.solar_power['Ground Diffuse'] = radiation[:,2].tolist()
        self.solar_power['Absorbed Solar Radiation (W/m^2)'] = self.solar_power['Beam'] + \
            self.solar_power['Total Sky Diffuse'] + self.solar_power['Ground Diffuse'] 
        # Cell Temperature 
        # [1] equation 9.35 pg 504
        self.solar_power['Cell Temp'] = (self.inoct-20)*(self.solar_power['Absorbed Solar Radiation (W/m^2)']/800)*\
            (1-((self.nom_eff/100)/0.9))+self.solar_data['Temperature']
        # DC Power [3]
        self.solar_power['DC Power'] =(self.dc_nameplate*(1+ self.gamma*(self.solar_power['Cell Temp']-25))*\
                                       (self.solar_power['Absorbed Solar Radiation (W/m^2)']/1000))*(1-self.system_loss)
        plr = np.zeros(self.solar_data.shape[0])
        plr = self.solar_power['DC Power'].values/((self.dc_nameplate/self.dc_ac)/(self.invert_eff))
        eta = np.zeros(self.solar_data.shape[0])
        ac_pow = np.zeros(self.solar_data.shape[0])
        filter_ = plr > 0
        if np.sum(filter_) > 0:
            eta[filter_] = (-0.0162*plr[filter_] + -0.0059/plr[filter_] + 0.9858)*self.invert_eff/0.9637
            ac_pow[filter_] = self.solar_power[filter_]['DC Power']*eta[filter_]*0.9
            #                                                                   |-> from here over is wrong and added to account
            #                                                                       for extra losses that PVWatts has that I don't
            ac_pow[np.where(ac_pow > self.dc_nameplate/self.dc_ac)] = self.dc_nameplate/self.dc_ac
        ac_pow[np.where(ac_pow < 0)] = 0.0
        self.solar_power['AC Power'] = ac_pow.tolist()               

    def perez(self, i)->list:
        # Perez Sky Diffuse calculation [2]
        F11R = [-0.0083117, 0.1299457, 0.3296958, 0.5682053, 0.8730280, 1.1326077, 1.0601591, 0.6777470]
        F12R = [0.5877285, 0.6825954, 0.4868735, 0.1874525, -0.3920403, -1.2367284, -1.5999137, -0.3272588]
        F13R = [-0.0620636, -0.1513752, -0.2210958, -0.2951290, -0.3616149, -0.4118494, -0.3589221, -0.2504286]
        F21R = [-0.0596012, -0.0189325, 0.0554140, 0.1088631, 0.2255647, 0.2877813, 0.2642124, 0.1561313]
        F22R = [0.0721249, 0.0659650, -0.0639588, -0.1519229, -0.4620442, -0.8230357, -1.1272340, -1.3765031]
        F23R = [-0.0220216, -0.0288748, -0.0260542, -0.0139754, 0.0012448, 0.0558651, 0.1310694, 0.2506212]
        EPSBINS = [1.065, 1.23, 1.5, 1.95, 2.8, 4.5, 6.2]
        out = [0.0, 0.0, 0.0]
        solar_data = self.solar_data.iloc[i]
        if solar_data['DNI'] < 0.0:
            solar_data['DNI'] = 0.0
        
        if solar_data['Solar Zenith Angle (deg)'] < 0.0 or np.radians(solar_data['Solar Zenith Angle (deg)']) > 1.5271631:
            # Zenith not between 0 and 87.5 deg 
            if solar_data['DHI'] < 0.0:
                solar_data['DHI'] = 0.0
            if np.cos(np.radians(self.solar_power.iloc[i]['Incident Angle'])) > 0.0 and np.radians(solar_data['Solar Zenith Angle (deg)']) < 1.5707963:
                # Zen between 87.5 and 90 and incident < 90 deg 
                out[0] = solar_data['DNI'] * np.cos(np.radians(self.solar_power.iloc[i]['Incident Angle']))
                out[2] = solar_data['DHI'] * (1.0 + np.cos(np.radians(self.tilt))) / 2.0
                return out
            else:
                # Isotropic diffuse only
                out[2] = solar_data['DHI'] * (1.0 + np.cos(np.radians(self.tilt))) / 2.0
                return out
        else:
            # Zen between 0 and 87.5 deg
            cz = np.cos(np.radians(solar_data['Solar Zenith Angle (deg)']))
            zh = cz
            if cz < 0.0871557: zh = 0.0871557
            if solar_data['DHI'] <= 0.0:
                # Diffuse is zero or less
                if np.cos(np.radians(self.solar_power.iloc[i]['Incident Angle'])) > 0.0:
                    # Incident < 90 deg
                    out[0] = solar_data['DNI'] * np.cos(np.radians(self.solar_power.iloc[i]['Incident Angle']))
                    return out
                else:
                    return out
            else:
                # Diffuse is greater than zero
                air_mass = 1.0 / (cz + 0.15 * pow(93.9 - solar_data['Solar Zenith Angle (deg)'], -1.253))
                delta = solar_data['DHI'] * air_mass / 1367.0
                t = pow(solar_data['Solar Zenith Angle (deg)'], 3.0)
                eps = (solar_data['DNI'] + solar_data['DHI']) / solar_data['DHI']
                eps = (eps + t * 0.000005534) / (1.0 + t * 0.000005534)
                qq = 0
                while qq < 7 and eps > EPSBINS[qq]:
                    qq += 1
                x = F11R[qq] + F12R[qq] * delta + F13R[qq] * np.radians(solar_data['Solar Zenith Angle (deg)'])
                f1 = max(0, x)
                f2 = F21R[qq] + F22R[qq] * delta + F23R[qq] * np.radians(solar_data['Solar Zenith Angle (deg)'])
                if np.cos(np.radians(self.solar_power.iloc[i]['Incident Angle'])) < 0.0:
                    zc = 0.0
                else:
                    zc = np.cos(np.radians(self.solar_power.iloc[i]['Incident Angle']))

                a = solar_data['DHI'] * (1 - f1) * (1.0 + np.cos(np.radians(self.tilt))) / 2.0 # isotropic diffuse
                b = solar_data['DHI'] * f1 * zc / zh # circumsolar diffuse
                c = solar_data['DHI'] * f2 * np.sin(np.radians(self.tilt)) # horizon brightness term

                out[0] = solar_data['DNI'] * zc # Beam
                out[1] = a + b + c # total sky diffuse
                out[2] = 0.2 * (solar_data['DNI'] * cz + solar_data['DHI']) * \
                    (1.0 - np.cos(np.radians(self.tilt))) / 2.0  # ground diffuse    
                return out

    def post_process(self):
        n_days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        for i in range(12):
            self.monthly_ac[i] = self.solar_power.loc[self.solar_power.index.month==i+1]['AC Power'].sum()
            self.monthly_rad[i] = (self.solar_power.loc[self.solar_power.index.month==i+1]['Absorbed Solar Radiation (W/m^2)'].sum()/n_days[i])/1000
        self.total = self.solar_power['AC Power'].sum()

    def analyze(self):
        '''
        '''
        self.calculate_solar_position()
        self.calculate_solar_power()
        self.post_process()


if __name__ == "__main__":    
    # Read required Lat/Lon
    solar_calc = solar_pv(38.886777, -77.02997, 374686.66666)
    solar_calc.analyze()
    
   