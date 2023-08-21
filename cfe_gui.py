from PyQt6.QtWidgets import *
from PyQt6.QtGui import *
from PyQt6.uic import loadUi

import sys
import pathlib
import json
import os
import pandas as pd
import numpy as np
import geopandas as gpd
import rasterio
from pyproj import Transformer
import time

import rc_cfe

PROJECT_PATH = pathlib.Path(__file__).parent
DATA_PATH = PROJECT_PATH / "data"
PROJECT_UI = PROJECT_PATH / "cfe_model.ui"
ASSUMP_DLG = PROJECT_PATH / "model_assumptions.ui"
DEFAULT_WHITE = u"background-color: rgb(255, 255, 255);"
DEFAULT_ERROR = u"background-color: rgb(255, 103, 103);"

class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super(MainWindow, self).__init__()

        # Load the user interface file
        loadUi(PROJECT_UI, self)

        # add status bar
        self.statusbar = self.statusBar()
        self.statusbar.showMessage("Ready", 50000)

        # Set model default assumptions from data folder
        assert os.path.exists(os.path.join(DATA_PATH, "model_assumptions.json")), \
            f"The following file could not be located:\n\t {os.path.join(DATA_PATH, 'model_assumptions.json')}\n\nThe program can't start without it."
        with open(os.path.join(DATA_PATH, 'model_assumptions.json')) as json_data:
            self.model_assump = json.load(json_data)

        # Load default FRPP Dataset
        self.frpp_path.setText(str(os.path.join(DATA_PATH, 'frpp_public_dataset_fy21_final_02242023.csv')))

        # Establish function connections
        ''' Landing Page '''
        self.frpp_SetNewPath_pushButton.pressed.connect(self.set_newFRPP_path)
        self.load_frpp_pushButton.pressed.connect(self.process_FRPP)
        ''' Build CFE Model Page '''
        self.model_assumpt_pushButton.pressed.connect(self.launch_model_assump_dlg)
        self.build_CFE_pushButton.pressed.connect(self.build_cfe_model)

        # Establish class variables
        self.frpp_df = None


    def set_newFRPP_path(self):
        """
        Sets a new path to a user supplied FRPP worksheet
        This can be either a csv or excel but it must contain
        the established headers. This will be checked when the
        user hits the `Load and Process Data` button.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        #options = QFileDialog.options()
        file, _ = QFileDialog.getOpenFileName(self,"Select an FRPP data set to Load", str(DATA_PATH), "CSV Files (*.csv);;Excel Files (*.xlsx)")
        if file:
            self.frpp_path.setText(str(file))
            if str(self.frpp_path.text()).split(".")[-1] == "xlsx":
                out_msg = "You can utilize excel files but it is recommended that\nyou save the sheet as a csv.\n\n"
                out_msg += "This will speed up the processing time by roughly\n150% (from roughly 3 minutes to 1 seocond)"
                msg = QMessageBox(self)
                msg.setIcon(QMessageBox.Icon.Information)
                msg.setText(out_msg)
                msg.setWindowTitle("Consider Converting to CSV")
                msg.setStandardButtons(QMessageBox.StandardButton.Ok)
                msg.exec()
        else:
            self.frpp_path.setText(str(os.path.join(DATA_PATH, 'frpp_public_dataset_fy21_final_02242023.csv')))

    def process_FRPP(self):
        """
        This function takes the path from `self.frpp_path` and 
        attempts to load the data from the file at that location.
        If the data is not an excel or csv file, or if the data
        does not match the anticipated headers, then an error will
        be raised.

        Parameters
        ----------
        str
            path to FRPP data

        Returns
        -------
        pd.dataframe
            self.frpp_df
        """
        # assure that the page is valid before proceeding
        valid, errors = self.validate_page1()
        if not valid:
            msg = QMessageBox(self)
            msg.setIcon(QMessageBox.Icon.Warning)
            msg.setText("\n\n".join(errors))
            msg.setWindowTitle("Warning")
            msg.setStandardButtons(QMessageBox.StandardButton.Ok)
            msg.exec()
            return
        
        # Proceed to reading the data
        self.statusbar.showMessage("Reading FRPP Data", 50000)
        self.landing_page.update()
        self.landing_page.repaint()
        # ---
        # we only want the following columns  
        column_names = ['Reporting Agency', 'Real Property Unique Identifier', 'State Name', 'US/Foreign', 'Zip Code',
                        'Latitude', 'Longitude', 'Real Property Type', 'Real Property Use', 'Asset Status','Acres', 
                        'Square Feet (Buildings)', 'Square Feet Unit of Measure', 'Unit Of Measure',
                        'Asset Height Range', 'Asset Height']
        
        column_dtype = {'Reporting Agency': str, 'Real Property Unique Identifier': str, 'State Name': str, 'US/Foreign': str,
                        'Zip Code':pd.Int64Dtype(),'Latitude': np.double, 'Longitude': np.double, 'Real Property Type': str, 
                        'Real Property Use': str, 'Asset Status': str, 'Acres': np.double, 'Square Feet (Buildings)': np.double, 
                        'Square Feet Unit of Measure': str, 'Unit Of Measure': str,
                        'Asset Height Range': str, 'Asset Height': np.double}
        # check if the path is a csv or excel
        if str(self.frpp_path.text()).split(".")[-1] == "csv":
            self.frpp_df = pd.read_csv(str(self.frpp_path.text()), header=0, 
                                       usecols=column_names, dtype=column_dtype, encoding = 'Latin')
        else:
            self.frpp_df = pd.read_excel(str(self.frpp_path.text()),usecols=column_names,
                                          dtype=column_dtype)
            
        # Apply agreed filters 
        # - only US and US Territories; No "Disposed" assets
        # - all buildings
        # - land: vacant, r&d, office building locations
        # - structures: parking structures
        self.frpp_df = self.frpp_df.loc[(self.frpp_df['Reporting Agency'] == str(self.agency_comboBox.currentText()))
            & (self.frpp_df['US/Foreign'].isin(['UNITED STATES', 'US TERRITORIES']))
            & (~self.frpp_df['Asset Status'].isin(['Disposed']))
            & ((self.frpp_df['Real Property Type'] == 'Building')
            | ((self.frpp_df['Real Property Type'] == 'Land') & (self.frpp_df['Real Property Use'].isin(['Vacant', 'Office Building Locations', 'Research and Development'])))
            | ((self.frpp_df['Real Property Type'] == 'Structure') & (self.frpp_df['Real Property Use'] == 'Parking Structures')))]
        
        # Clean Data
        self.frpp_df['Zip Code'] = self.frpp_df['Zip Code'].apply(lambda x: int(str(int(x))[0:5]) if len(str(int(x))) > 5 else int(x))
        self.pg1_progressBar.setValue(int(1/6*100))
        self.statusbar.showMessage("Estimating Rooftop Area", 50000)
        self.landing_page.update()
        self.landing_page.repaint()

        # Estimate Rooftop Area for Buildings
        '''        
        Adapted from PVWattsv8_Implementation notebook
        This code estimates rooftop SIZE only, and does not apply any assumptions 
        of how much of that rooftop area is available for solar.

        Steps:
        - Subset for BUILDINGS only
        - Based on available data, creates two engineered features:
            - est_num_stories
            - est_rooftop_area_sqft
        - Note: FRPP column "Square Feet Unit of Measure" distinguishes whether Sq Ft is reported in gross or rentable sq ft.
        '''
        height_to_stories = {'Height > 0 feet and <= 30 feet above ground level': 1,
                             'Height > 30 feet and <= 100 feet above ground level': 6, 
                             'Height > 100 feet and < 200 feet above ground level': 13, 
                             'Height >= 200 feet above ground level': 22}
        
        self.frpp_df['est_num_stories'] = list(np.zeros((self.frpp_df.shape[0],1)))
        self.frpp_df['est_rooftop_area_sqft'] = list(np.zeros((self.frpp_df.shape[0],1)))
        for index, row in self.frpp_df.iterrows():
            if row['Real Property Type'] == 'Building':
                if not np.isnan(row['Asset Height']):
                    self.frpp_df.at[index, 'est_num_stories'] = round(row['Asset Height']/12, 0)
                elif isinstance(row['Asset Height Range'], str):
                    self.frpp_df.at[index, 'est_num_stories'] = height_to_stories[row['Asset Height Range']]
                else:
                    self.frpp_df.at[index, 'est_num_stories'] = 2
                self.frpp_df.at[index, 'est_rooftop_area_sqft'] = row['Square Feet (Buildings)'] / self.frpp_df.at[index, 'est_num_stories']
            else:
                self.frpp_df.at[index, 'est_num_stories'] = np.nan
                self.frpp_df.at[index, 'est_rooftop_area_sqft'] = np.nan
        
        self.pg1_progressBar.setValue(int(2/6*100))
        self.statusbar.showMessage("Adding Geothermal Resource Data", 50000)
        self.landing_page.update()
        self.landing_page.repaint()

        '''
        Add Geothermal Resource Data
        **Geothermal Resource Data**
        - NREL Map/Legend for lower 48 contiguous: https://www.nrel.gov/gis/geothermal.html
        - NREL Geothermal Shapefiles: https://www.nrel.gov/gis/assets/data/egs.zip
        Define subset of RP that we will consider for geothermal (vacant land)
        '''
        assets_geo = self.frpp_df.loc[(self.frpp_df['Real Property Type'] == 'Land') & (self.frpp_df['Real Property Use'] == 'Vacant')].copy()        
        geo_path = os.path.join(DATA_PATH,'egs','lower-48-geothermal-high-resolution')
        # load 
        # - .shp: file that contains the geometry for all features
        # - .shx: file that indexes the geometry
        # - .dbf: file that stores feature attributes in tabular format
        # - .prj: contains info on projection format including cooridnate system and projection info
        # - (.sbn adn .sbx spatial index of features)
        # - (.shp.xml: geospatial metadata in XML format)

        # read in shape file (geopandas uses fiona to read in all layers in the folder)
        # Required before spatial join: Ensure CRS systems align between the geometry fields in the two datasets
        # -FRPP Data: lat/lon specifies WSG84, which corresponds to EPSG4326 using EPSG.
        # - Note: Re-projecting FRPP data to the geothermal data CRS (EPSG:5070) does not work. 
        # Try re-projecting geothermal to EPSG4326, which is the FRPP CRS
        geothermal = gpd.read_file(geo_path).to_crs('EPSG:4326')
        geo_df = gpd.GeoDataFrame(assets_geo, geometry = gpd.points_from_xy(x=assets_geo.Longitude, y=assets_geo.Latitude), crs='EPSG:4326')
        # Do a spatial join between geo shape files and frpp latlon
        merged_df = geo_df.sjoin(geothermal, how = 'left')
        merged_df.rename(columns = {'CLASS': 'Geothermal_CLASS'}, inplace = True)
        merged_df.drop(['geometry', 'index_right'], axis = 1, inplace = True)
        self.frpp_df = self.frpp_df.merge(merged_df[['Real Property Unique Identifier', 'Geothermal_CLASS']], how = 'left', on = 'Real Property Unique Identifier')

        self.pg1_progressBar.setValue(int(3/6*100))
        self.statusbar.showMessage("Adding Wind Resource Data", 50000)
        self.landing_page.update()
        self.landing_page.repaint()
        '''
        Add Wind Resource Data
        These static maps illustrate multiyear average wind speeds at various heights derived from NREL's WIND Toolkit.
        https://www.nrel.gov/gis/wind-resource-maps.html
        Downloaded as raster files - elements ("pixels") mapped to regions on earth's surface, pixel per spatial bounding box

        **Citation**
        Draxl, C., B.M. Hodge, A. Clifton, and J. McCaa. 2015. Overview and Meteorological Validation of the Wind Integration National Dataset Toolkit (Technical Report, NREL/TP-5000-61740). Golden, CO: National Renewable Energy Laboratory.
        
        Draxl, C., B.M. Hodge, A. Clifton, and J. McCaa. 2015. "The Wind Integration National Dataset (WIND) Toolkit." Applied Energy 151: 355366.
        
        Lieberman-Cribbin, W., C. Draxl, and A. Clifton. 2014. Guide to Using the WIND Toolkit Validation Code (Technical Report, NREL/TP-5000-62595). Golden, CO: National Renewable Energy Laboratory.
        
        King, J., A. Clifton, and B.M. Hodge. 2014. Validation of Power Output for the WIND Toolkit (Technical Report, NREL/TP-5D00-61714). Golden, CO: National Renewable Energy Laboratory.
        '''
        # Gather wind data for Land - Vacant, R&D - RP categories only
        assets_wind = self.frpp_df.loc[(self.frpp_df['Real Property Type'] == 'Land') & (self.frpp_df['Real Property Use'].isin(['Vacant', 'Research and Development']))].copy()
        # Read in wind file (100m as per previous discussions)
        # Note metadata documents state spatial resolution is 2km x 2km
        wind = rasterio.open(os.path.join(DATA_PATH,'us-wind-data','us-wind-data','wtk_conus_100m_mean_masked.tif'))
        wind_values = []
        for x, y in assets_wind[['Latitude', 'Longitude']].iterrows():
        #     print(y[0], y[1])
            lat = y[0]
            lon = y[1]
            transformer = Transformer.from_crs('EPSG:4326', wind.crs, always_xy = True)
            xx, yy = transformer.transform(lon, lat)
            
            value = list(wind.sample([(xx, yy)]))[0]
            wind_values.append(value)
        assets_wind['Annual Average Wind Speed (m/s) at 100m above surface level'] = wind_values[0]
        self.frpp_df = self.frpp_df.merge(assets_wind[['Real Property Unique Identifier', 'Annual Average Wind Speed (m/s) at 100m above surface level']], how = 'left', on = 'Real Property Unique Identifier')

        self.pg1_progressBar.setValue(int(4/6*100))
        self.statusbar.showMessage("Adding Solar Resource Data", 50000)
        self.landing_page.update()
        self.landing_page.repaint()
        '''
         Add Solar Resource Data
         **Solar data used in existing tools**:
         - **REopt**: Solar production data is taken from the PVWatts dataset, which includes many international  locations. 
                      The REopt web tool will use the closest available location that is found to have resource data, so the 
                      user should independently confirm that PVWatts includes data for a  location that is acceptably close to 
                      their site location. The available resource data locations can be found using NREL's PV Watts. Users who 
                      have access to hourly custom solar production data for their site can upload it in the Advanced Inputs 
                      section, and it will be used instead of PVWatts data." p. 23, https://reopt.nrel.gov/tool/reopt-user-manual.pdf
        
         - **PVWatts**: Solar resource data is solar irradiance and meteorological data that describe the conditions at the system's 
                        location. PVWatts® uses a set of weather data prepared from the NREL National Solar Radiation Database (NSRDB) 
                        where it is available, and a collection of data from other sources for the rest of the world. PVWatts® uses 
                        hourly typical meteorological year (TMY) data, which is one year's worth of data that represents the solar resource 
                        over a multi-year period. For more details about TMY data, see this NSRDB article on TMY data. When you type a 
                        street address or latitude and longitude for the system's location, PVWatts®shows the latitude and longitude of the 
                        NSRDB grid cell for the location and the distance in miles between the location and the grid cell center. For 
                        locations outside of the NSRDB area, it shows the coordinates or name of the nearest available solar resource data 
                        site. PVWatts® also displays a map with a rectangle representing the NSRDB grid cell and a red pin indicating the system 
                        location. For locations outside of the NSRDB area, the map shows pins for the nearest available data sites.
        
         - Note. It is not possible to run PVWatts® using your own solar resource data file or a data from a source other than those discussed here. 
           If you want to run PVWatts® simulations with your own solar resource data file, you can use the version of PVWatts® in NREL's System Advisor Model.
        
         - The NSRDB solar resource data for PVWatts is a special set of files from the NSRDB. These files were collected from the following NSRDB datasets:
             - PSM V3 TMY (tmy-2020)
             - Himarawi PSM V3 TMY (tmy-2020)
         - Solar resource data sources for locations not covered by the NSRDB include:
             - Solar and Wind Energy Resource Assessment Programme (SWERA)
             - The ASHRAE International Weather for Energy Calculations Version 1.1 (IWEC)
             - Canadian Weather for Energy Calculations (CWEC)
        '''
        solar = pd.read_csv(os.path.join(DATA_PATH,'20230524_FRPP_Solar_Part6.csv'))
        # Remove 0 zipcode - that happened bc an asset had a latlon and not a zip. Not relevant for GSA. 
        # 60 zipcodes returned 'Unprocessable Entity' error. We will drop those as well for now
        solar = solar.loc[(solar['Solar Value']!="['Unprocessable Entity']") & (solar['Zip Code'] != 0)]

        # Take the median value per zipcode
        clean_solar = pd.DataFrame(solar.groupby('Zip Code')['Solar Value'].median())
        clean_solar.reset_index(inplace = True)
        clean_solar.columns = ['Zip Code', 'Solar Value']
        self.frpp_df = self.frpp_df.merge(clean_solar, how = 'left', on = 'Zip Code')

        self.frpp_df.rename(columns = {'Solar Value': 'Annual Solar Radiation kWh/m2/day'}, inplace= True)
        self.frpp_df = self.frpp_df.drop(['Square Feet Unit of Measure', 'Unit Of Measure', 'Asset Height Range', 'Asset Height'], axis=1)

        self.pg1_progressBar.setValue(int(5/6*100))
        self.statusbar.showMessage("FRPP Data Read Successfully: Building Table", 50000)
        self.landing_page.update()
        self.landing_page.repaint()
        # Set the values for the next page
        self.agency_table_label.setText(f"FRPP Data for {str(self.agency_comboBox.currentText())}")
        self.nAssets_lineEdit.setText(str(self.frpp_df.shape[0]))
        self.total_land_lineEdit.setText(f'{self.frpp_df["Acres"].sum():.2f}')
        self.total_build_lineEdit.setText(f'{self.frpp_df["est_rooftop_area_sqft"].sum():.2f}')
        self.vacant_lineEdit.setText(f"{self.frpp_df.loc[(self.frpp_df['Real Property Type'] == 'Land') & (self.frpp_df['Real Property Use'] == 'Vacant')]['Acres'].sum():.2f}")
        self.RnD_lineEdit.setText(f"{self.frpp_df.loc[(self.frpp_df['Real Property Type'] == 'Land') & (self.frpp_df['Real Property Use'] == 'Research and Development')]['Acres'].sum():.2f}")
        self.office_lineEdit.setText(f"{self.frpp_df.loc[(self.frpp_df['Real Property Type'] == 'Land') & (self.frpp_df['Real Property Use'] == 'Office Building Locations')]['Acres'].sum():.2f}")
        # Build table from frpp data to show user
        self.agency_tableWidget.setRowCount(self.frpp_df.shape[0])
        self.agency_tableWidget.setColumnCount(self.frpp_df.shape[1])
        self.agency_tableWidget.setHorizontalHeaderLabels(self.frpp_df.columns.values)
        for kk in range(self.frpp_df.shape[0]):
            for qq, header in enumerate(self.frpp_df.columns.values):
                self.agency_tableWidget.setItem(kk,qq, QTableWidgetItem(str(self.frpp_df.iloc[kk][header])))

        self.pg1_progressBar.setValue(int(6/6*100))
        self.landing_page.update()
        self.landing_page.repaint()
        time.sleep(1)
        self.stackedWidget.setCurrentIndex(1)

    def launch_model_assump_dlg(self):
        dlg = model_assump_dlg(self.model_assump)
        if dlg.exec():
            if dlg.overwrite_checkBox.isChecked():
                print("write new json")


    def build_cfe_model(self):
        a = 1

    # Functions not tied to any button
    def validate_page1(self)->bool | list:
        """
        Goes through all of page one and makes sure that all
        entries are valid.

        Parameters
        ----------
        None

        Returns
        -------
        bool
            valid
        list
            errors
        """

        errors = ["Please address the following issues to continue:"]
        valid = True

        # Reset all page options to valid background
        self.frpp_path.setStyleSheet(DEFAULT_WHITE)
        self.agency_comboBox.setStyleSheet(DEFAULT_WHITE)

        # Check the path
        # The only way the user could mess with this is to set the path then delete it before hitting proceed
        if not os.path.exists(str(self.frpp_path.text())):
            self.frpp_path.setStyleSheet(DEFAULT_ERROR)
            valid = False
            errors.append(f"The path specified {self.frpp_path.text()}\nDoes not exist. Please select a valid FRPP dataset to continue.")

        # Check that an agency has been selected
        if self.agency_comboBox.currentIndex() == 0:
            self.agency_comboBox.setStyleSheet(DEFAULT_ERROR)
            valid = False
            errors.append("You must select an agency to continue")

        return valid, errors

        

class model_assump_dlg(QDialog):
    def __init__(self, assump_json):
        super(model_assump_dlg, self).__init__()
        loadUi(ASSUMP_DLG, self)
        # Load the values from the assumption json
        ####   Wind  ####
        self.wind_frac_lineEdit.setText(str(assump_json['wind']['system']['fraction_of_average_wind_speed']))        
        self.wind_turb_diam_lineEdit.setText(str(assump_json['wind']['system']['turbine_diameter']))
        self.wind_CapFac_lineEdit.setText(str(assump_json['wind']['system']['capacity_factor']))
        self.wind_Nturb_lineEdit.setText(str(assump_json['wind']['system']['N_turbines']))
        self.wind_spacing_lineEdit.setText(str(assump_json['wind']['system']['turbine_spacing'])) 
        self.wind_pow_coef_lineEdit.setText(str(assump_json['wind']['system']['power_coefficient']))
        self.wind_eff_lineEdit.setText(str(assump_json['wind']['system']['overall_efficiency']))	
        self.wg_disc_rate_lineEdit.setText(str(assump_json['wind']['gov_rates']['discount_rate']))
        self.wg_infl_rate_lineEdit.setText(str(assump_json['wind']['gov_rates']['inflation_rate']))
        self.wg_ece_LB_lineEdit.setText(str(assump_json['wind']['gov_rates']['electric_cost_escalation_lb'])) 
        self.wg_ece_UB_lineEdit.setText(str(assump_json['wind']['gov_rates']['electric_cost_escalation_ub']))
        self.wg_tax_lineEdit.setText(str(assump_json['wind']['gov_rates']['taxation']))
        self.wg_proj_life_lineEdit.setText(str(assump_json['wind']['gov_rates']['project_life']))
        self.wt_disc_rate_lineEdit.setText(str(assump_json['wind']['third_party_rates']['discount_rate']))
        self.wt_infl_rate_lineEdit.setText(str(assump_json['wind']['third_party_rates']['inflation_rate'])) 
        self.wt_ece_LB_lineEdit.setText(str(assump_json['wind']['third_party_rates']['electric_cost_escalation_lb'])) 
        self.wt_ece_UB_lineEdit.setText(str(assump_json['wind']['third_party_rates']['electric_cost_escalation_ub']))
        self.wt_tax_lineEdit.setText(str(assump_json['wind']['third_party_rates']['taxation'])) 
        self.wt_proj_life_lineEdit.setText(str(assump_json['wind']['third_party_rates']['project_life']))
        ####   Solar  ####
        self.per_roof_lineEdit.setText(str(assump_json['solar']['system']['precent_roof_avail']))
        self.rad_2_elec_lineEdit.setText(str(assump_json['solar']['system']['per_rad_to_elec']))
        self.fract_avg_SolRad_lineEdit.setText(str(assump_json['solar']['system']['avg_radiation_frac'])) 
        self.perc_land_used_lineEdit.setText(str(assump_json['solar']['system']['perc_land_used'])) 
        self.sun_hours_lineEdit.setText(str(assump_json['solar']['system']['avg_sun_hours'])) 
        self.dc_ac_lineEdit.setText(str(assump_json['solar']['system']['dc_to_ac_ratio'])) 
        self.cap_fac_lineEdit.setText(str(assump_json['solar']['system']['capacity_factor']))
        self.sg_disc_rate_lineEdit.setText(str(assump_json['solar']['gov_rates']['discount_rate'])) 
        self.sg_infl_rate_lineEdit.setText(str(assump_json['solar']['gov_rates']['inflation_rate'])) 
        self.sg_ece_LB_lineEdit.setText(str(assump_json['solar']['gov_rates']['electric_cost_escalation_lb'])) 
        self.sg_ece_UB_lineEdit.setText(str(assump_json['solar']['gov_rates']['electric_cost_escalation_ub'])) 
        self.sg_tax_lineEdit.setText(str(assump_json['solar']['gov_rates']['taxation'])) 
        self.sg_proj_life_lineEdit.setText(str(assump_json['solar']['gov_rates']['project_life'])) 
        self.st_disc_rate_lineEdit.setText(str(assump_json['solar']['third_party_rates']['discount_rate'])) 
        self.st_infl_rate_lineEdit.setText(str(assump_json['solar']['third_party_rates']['inflation_rate']))
        self.st_ece_LB_lineEdit.setText(str(assump_json['solar']['third_party_rates']['electric_cost_escalation_lb'])) 
        self.st_ece_UB_lineEdit.setText(str(assump_json['solar']['third_party_rates']['electric_cost_escalation_ub'])) 
        self.st_tax_lineEdit.setText(str(assump_json['solar']['third_party_rates']['taxation'])) 
        self.st_proj_life_lineEdit.setText(str(assump_json['solar']['third_party_rates']['project_life'])) 
        ####   Concentrating Solar  ####
        self.opt_eff_lineEdit.setTextstr(assump_json['conc_solar']['system']['optical_eff']) 
        self.elec_eff_lineEdit.setText(assump_json['conc_solar']['system']['elec_eff']) 
        self.therm_stor_lineEdit.setText(assump_json['conc_solar']['system']['thermal_storage'])
        self.avg_sunhrs_lineEdit.setText(assump_json['conc_solar']['system']['avg_sun_hours']) 
        self.avg_sol_rad_lineEdit.setText(assump_json['conc_solar']['system']['avg_solar_rad'])
        self.per_agency_land_lineEdit.setText(assump_json['conc_solar']['system']['perc_land_used']) 		
        self.cs_cap_fac_lineEdit.setText(assump_json['conc_solar']['system']['capacity_factor']) 
        self.csg_disc_rate_lineEdit.setText(str(assump_json['conc_solar']['gov_rates']['discount_rate'])) 
        self.csg_infl_rate_lineEdit.setText(str(assump_json['conc_solar']['gov_rates']['inflation_rate'])) 
        self.csg_ece_LB_lineEdit.setText(str(assump_json['conc_solar']['gov_rates']['electric_cost_escalation_lb'])) 
        self.csg_ece_UB_lineEdit.setText(str(assump_json['conc_solar']['gov_rates']['electric_cost_escalation_ub'])) 
        self.csg_tax_lineEdit.setText(str(assump_json['conc_solar']['gov_rates']['taxation'])) 
        self.csg_proj_life_lineEdit.setText(str(assump_json['conc_solar']['gov_rates']['project_life']))
        self.cst_disc_rate_lineEdit.setText(str(assump_json['conc_solar']['third_party_rates']['discount_rate']))
        self.cst_infl_rate_lineEdit.setText(str(assump_json['conc_solar']['third_party_rates']['inflation_rate'])) 
        self.cst_ece_LB_lineEdit.setText(str(assump_json['conc_solar']['third_party_rates']['electric_cost_escalation_lb'])) 
        self.cst_ece_UB_lineEdit.setText(str(assump_json['conc_solar']['third_party_rates']['electric_cost_escalation_ub'])) 
        self.cst_tax_lineEdit.setText(str(assump_json['conc_solar']['third_party_rates']['taxation'])) 
        self.cst_proj_life_lineEdit.setText(str(assump_json['conc_solar']['third_party_rates']['project_life']))
        ####   Hydrogen  ####	
        self.h_cap_fac_lineEdit.setText(str(assump_json['hydrogen']['system']['capacity_factor'])) 
        self.hrs_oper_lineEdit.setText(str(assump_json['hydrogen']['system']['daily_operation'])) 
        for ind, entry in enumerate([self.fuel_emp_comboBox.itemText(i) for i in range(self.fuel_emp_comboBox.count())]):
            if entry == str(assump_json['hydrogen']['system']['fuel_used']):
                self.fuel_emp_comboBox.setCurrentIndex(ind)
        self.cur_den_lineEdit.setText(str(assump_json['hydrogen']['system']['curr_density'])) 
        self.cell_volt_lineEdit.setText(str(assump_json['hydrogen']['system']['cell_active_area'])) 
        self.cell_act_area_lineEdit.setText(str(assump_json['hydrogen']['system']['cell_voltage'])) 
        self.Nstacks_lineEdit.setText(str(assump_json['hydrogen']['system']['N_stacks'])) 
        # If we ever add the other option we will have to change this
        self.prot_exch_comboBox.setCurrentIndex(0)
        self.heat_rec_comboBox.setCurrentIndex(0)
        self.heat_val_lineEdit.setText(str(assump_json['hydrogen']['system']['heating_value'])) 
        self.stack_replace_lineEdit.setText(str(assump_json['hydrogen']['system']['stack_replacement'])) 
        self.annual_degrade_lineEdit.setText(str(assump_json['hydrogen']['system']['annual_degredation'])) 
        self.hg_disc_rate_lineEdit.setText(str(assump_json['hydrogen']['gov_rates']['discount_rate'])) 
        self.hg_infl_rate_lineEdit.setText(str(assump_json['hydrogen']['gov_rates']['inflation_rate'])) 
        self.hg_ece_LB_lineEdit.setText(str(assump_json['hydrogen']['gov_rates']['electric_cost_escalation_lb'])) 
        self.hg_ece_UB_lineEdit.setText(str(assump_json['hydrogen']['gov_rates']['electric_cost_escalation_ub'])) 
        self.hg_tax_lineEdit.setText(str(assump_json['hydrogen']['gov_rates']['taxation'])) 
        self.hg_proj_life_lineEdit.setText(str(assump_json['hydrogen']['gov_rates']['project_life'])) 
        self.ht_disc_rate_lineEdit.setText(str(assump_json['hydrogen']['third_party_rates']['discount_rate'])) 
        self.ht_infl_rate_lineEdit.setText(str(assump_json['hydrogen']['third_party_rates']['inflation_rate'])) 
        self.ht_ece_LB_lineEdit.setText(str(assump_json['hydrogen']['third_party_rates']['electric_cost_escalation_lb'])) 
        self.ht_ece_UB_lineEdit.setText(str(assump_json['hydrogen']['third_party_rates']['electric_cost_escalation_ub'])) 
        self.ht_tax_lineEdit.setText(str(assump_json['hydrogen']['third_party_rates']['taxation'])) 
        self.ht_proj_life_lineEdit.setText(str(assump_json['hydrogen']['third_party_rates']['project_life'])) 	
        ####   Geothermal  ####	
        self.well_diam_lineEdit.setText(str(assump_json['geo_therm']['system']['well_diameter'])) 
        self.fluid_vel_lineEdit.setText(str(assump_json['geo_therm']['system']['fluid_velocity'])) 
        self.well_depth_lineEdit.setText(str(assump_json['geo_therm']['system']['well_depth'])) 
        self.turb_outlet_lineEdit.setText(str(assump_json['geo_therm']['system']['turb_outlet_temp'])) 
        self.geo_cap_fac_lineEdit.setText(str(assump_json['geo_therm']['system']['capacity_factor'])) 
        self.fluid_den_lineEdit.setText(str(assump_json['geo_therm']['system']['fluid_density'])) 
        self.pump_eff_lineEdit.setText(str(assump_json['geo_therm']['system']['overall_pump_efficiency'])) 	
        self.gg_disc_rate_lineEdit.setText(str(assump_json['geo_therm']['gov_rates']['discount_rate'])) 
        self.gg_infl_rate_lineEdit.setText(str(assump_json['geo_therm']['gov_rates']['inflation_rate']))
        self.gg_ece_LB_lineEdit.setText(str(assump_json['geo_therm']['gov_rates']['electric_cost_escalation_lb'])) 
        self.gg_ece_UB_lineEdit.setText(str(assump_json['geo_therm']['gov_rates']['electric_cost_escalation_ub'])) 
        self.gg_tax_lineEdit.setText(str(assump_json['geo_therm']['gov_rates']['taxation'])) 
        self.gg_proj_life_lineEdit.setText(str(assump_json['geo_therm']['gov_rates']['project_life']))
        self.gt_disc_rate_lineEdit.setText(str(assump_json['geo_therm']['third_party_rates']['discount_rate']))
        self.gt_infl_rate_lineEdit.setText(str(assump_json['geo_therm']['third_party_rates']['inflation_rate'])) 
        self.gt_ece_LB_lineEdit.setText(str(assump_json['geo_therm']['third_party_rates']['electric_cost_escalation_lb'])) 
        self.gt_ece_UB_lineEdit.setText(str(assump_json['geo_therm']['third_party_rates']['electric_cost_escalation_ub'])) 
        self.gt_tax_lineEdit.setText(str(assump_json['geo_therm']['third_party_rates']['taxation'])) 
        self.gt_proj_life_lineEdit.setText(str(assump_json['geo_therm']['third_party_rates']['project_life']))     
        # Set functions
        self.ok_button = self.buttonBox.button(QDialogButtonBox.StandardButton.Ok)
        self.ok_button.clicked.connect(self.validate)

    def validate(self)->bool | list:
        """
        Goes through all of page one and makes sure that all
        entries are valid.

        Parameters
        ----------
        None

        Returns
        -------
        bool
            valid
        list
            errors
        """

        errors = ["Please address the following issues to continue:"]
        valid = True

        # Reset all page options to valid background
        self.wind_frac_lineEdit.setStyleSheet(DEFAULT_WHITE)        
        self.wind_turb_diam_lineEdit.setStyleSheet(DEFAULT_WHITE)
        self.wind_CapFac_lineEdit.setStyleSheet(DEFAULT_WHITE)
        self.wind_Nturb_lineEdit.setStyleSheet(DEFAULT_WHITE)
        self.wind_spacing_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.wind_pow_coef_lineEdit.setStyleSheet(DEFAULT_WHITE)
        self.wind_eff_lineEdit.setStyleSheet(DEFAULT_WHITE)		
        self.wg_disc_rate_lineEdit.setStyleSheet(DEFAULT_WHITE)
        self.wg_infl_rate_lineEdit.setStyleSheet(DEFAULT_WHITE)
        self.wg_ece_LB_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.wg_ece_UB_lineEdit.setStyleSheet(DEFAULT_WHITE)
        self.wg_tax_lineEdit.setStyleSheet(DEFAULT_WHITE)
        self.wg_proj_life_lineEdit.setStyleSheet(DEFAULT_WHITE)
        self.wt_disc_rate_lineEdit.setStyleSheet(DEFAULT_WHITE)
        self.wt_infl_rate_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.wt_ece_LB_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.wt_ece_UB_lineEdit.setStyleSheet(DEFAULT_WHITE)
        self.wt_tax_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.wt_proj_life_lineEdit.setStyleSheet(DEFAULT_WHITE)		
        self.per_roof_lineEdit.setStyleSheet(DEFAULT_WHITE)
        self.rad_2_elec_lineEdit.setStyleSheet(DEFAULT_WHITE)
        self.fract_avg_SolRad_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.perc_land_used_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.sun_hours_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.dc_ac_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.cap_fac_lineEdit.setStyleSheet(DEFAULT_WHITE)		
        self.sg_disc_rate_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.sg_infl_rate_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.sg_ece_LB_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.sg_ece_UB_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.sg_tax_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.sg_proj_life_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.st_disc_rate_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.st_infl_rate_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.st_ece_LB_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.st_ece_UB_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.st_tax_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.st_proj_life_lineEdit.setStyleSheet(DEFAULT_WHITE) 		
        self.opt_eff_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.elec_eff_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.therm_stor_lineEdit.setStyleSheet(DEFAULT_WHITE)
        self.avg_sunhrs_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.avg_sol_rad_lineEdit.setStyleSheet(DEFAULT_WHITE)
        self.per_agency_land_lineEdit.setStyleSheet(DEFAULT_WHITE) 		
        self.cs_cap_fac_lineEdit.setStyleSheet(DEFAULT_WHITE) 		
        self.csg_disc_rate_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.csg_infl_rate_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.csg_ece_LB_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.csg_ece_UB_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.csg_tax_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.csg_proj_life_lineEdit.setStyleSheet(DEFAULT_WHITE)
        self.cst_disc_rate_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.cst_infl_rate_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.cst_ece_LB_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.cst_ece_UB_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.cst_tax_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.cst_proj_life_lineEdit.setStyleSheet(DEFAULT_WHITE)		
        self.h_cap_fac_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.hrs_oper_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.fuel_emp_comboBox.setStyleSheet(DEFAULT_WHITE) 
        self.cur_den_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.cell_volt_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.cell_act_area_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.Nstacks_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.prot_exch_comboBox.setStyleSheet(DEFAULT_WHITE)  
        self.heat_rec_comboBox.setStyleSheet(DEFAULT_WHITE)  
        self.heat_val_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.stack_replace_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.annual_degrade_lineEdit.setStyleSheet(DEFAULT_WHITE) 		
        self.hg_disc_rate_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.hg_infl_rate_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.hg_ece_LB_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.hg_ece_UB_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.hg_tax_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.hg_proj_life_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.ht_disc_rate_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.ht_infl_rate_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.ht_ece_LB_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.ht_ece_UB_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.ht_tax_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.ht_proj_life_lineEdit.setStyleSheet(DEFAULT_WHITE) 		
        self.well_diam_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.fluid_vel_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.well_depth_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.turb_outlet_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.geo_cap_fac_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.fluid_den_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.pump_eff_lineEdit.setStyleSheet(DEFAULT_WHITE) 		
        self.gg_disc_rate_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.gg_infl_rate_lineEdit.setStyleSheet(DEFAULT_WHITE)
        self.gg_ece_LB_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.gg_ece_UB_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.gg_tax_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.gg_proj_life_lineEdit.setStyleSheet(DEFAULT_WHITE)
        self.gt_disc_rate_lineEdit.setStyleSheet(DEFAULT_WHITE)
        self.gt_infl_rate_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.gt_ece_LB_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.gt_ece_UB_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.gt_tax_lineEdit.setStyleSheet(DEFAULT_WHITE) 
        self.gt_proj_life_lineEdit.setStyleSheet(DEFAULT_WHITE)

        # Check wind
        diam = np.nan
        try:
            val = float(self.wind_frac_lineEdit.text()) 
            if val < 0:
                self.wind_frac_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The Fraction of Averate Wind Speed must be greater than zero")
            elif val > 3:
                self.wind_frac_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The Fraction of Averate Wind Speed must be less than 3")   
        except:
            self.wind_frac_lineEdit.setStyleSheet(DEFAULT_ERROR)
            valid = False
            errors.append(f"The Fraction of Averate Wind Speed must be numeric")

        try:
            val = float(self.self.wind_turb_diam_lineEdit.text()) 
            diam = val
            if val < 9:
                self.wind_turb_diam_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The Turbine Rotor Diameter must be greater than or equal to 9 meters")
            elif val > 252:
                self.wind_turb_diam_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The Turbine Rotor Diameter must be less than or equal to 252 meters")   
        except:
            self.wind_turb_diam_lineEdit.setStyleSheet(DEFAULT_ERROR)
            valid = False
            errors.append(f"The Turbine Rotor Diameter must be numeric")

        try:
            val = float(self.self.wind_CapFac_lineEdit.text()) 
            if val <= 0:
                self.wind_CapFac_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The Wind Capacity Factor must be greater than 0%")
            elif val > 100:
                self.wind_CapFac_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The Wind Capacity Factor must be less than or equal to 100%")   
        except:
            self.wind_CapFac_lineEdit.setStyleSheet(DEFAULT_ERROR)
            valid = False
            errors.append(f"The Wind Capacity Factor must be numeric")

        # I think the spacing is wrong. Nothing supports a 2*diameter. 5 was the lowest https://www.renewablesfirst.co.uk/renewable-energy-technologies/windpower/community-windpower/location-size-no-of-wind-turbines/#:~:text=The%20number%20of%20wind%20turbines,turbine%20it%20is%20410%20metres.
        # I saw and SAM has it as 8 https://github.com/NREL/SAM/blob/85b6974a7f9f02e931b9013b7d250a7eccd87b1c/deploy/runtime/ui/Wind%20Farm%20Specifications.txt#L1631C30-L1631C30
        try:
            val = int(self.self.wind_Nturb_lineEdit.text())              
        except:
            self.wind_Nturb_lineEdit.setStyleSheet(DEFAULT_ERROR)
            valid = False
            errors.append(f"The number of turbines must be numeric")
        
        try:
            val = int(self.self.wind_spacing_lineEdit.text())              
        except:
            self.wind_spacing_lineEdit.setStyleSheet(DEFAULT_ERROR)
            valid = False
            errors.append(f"The distance between turbines must be numeric")

        if not np.nan(diam):
            try:
                nturb = int(self.self.wind_Nturb_lineEdit.text())
                req_spa = int(self.self.wind_spacing_lineEdit.text())  
                spacing = (req_spa/nturb)/diam
                if spacing < 5:
                    self.wind_Nturb_lineEdit.setStyleSheet(DEFAULT_ERROR)
                    self.wind_spacing_lineEdit.setStyleSheet(DEFAULT_ERROR)
                    valid = False
                    errors.append(f"The distance between turbines must be greater than or equal to 5 turbine diameters. Currently: {spacing}")
            except:
                pass
        
        try:
            val = float(self.self.wind_pow_coef_lineEdit.text()) 
            if val <= 0:
                self.wind_pow_coef_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The Overall Electrical Efficiency must be greater than 0%")
            elif val > 100:
                self.wind_pow_coef_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The Overall Electrical Efficiency must be less than or equal to 100%")   
        except:
            self.wind_pow_coef_lineEdit.setStyleSheet(DEFAULT_ERROR)
            valid = False
            errors.append(f"The Overall Electrical Efficiency must be numeric")
        
        try:
            val = float(self.self.wind_eff_lineEdit.text()) 
            if val <= 0:
                self.wind_eff_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The wind power coefficient must be greater than 0")
            elif val > 1:
                self.wind_eff_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The wind power coefficient must be less than or equal to 1")   
        except:
            self.wind_eff_lineEdit.setStyleSheet(DEFAULT_ERROR)
            valid = False
            errors.append(f"The wind power coefficient must be numeric")
        # Check Solar
        try:
            val = float(self.self.per_roof_lineEdit.text()) 
            if val <= 0:
                self.per_roof_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The percent of rooftop available for solar must be greater than 0%")
            elif val > 100:
                self.per_roof_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The percent of rooftop available for solar must be less than or equal to 100%")   
        except:
            self.per_roof_lineEdit.setStyleSheet(DEFAULT_ERROR)
            valid = False
            errors.append(f"The percent of rooftop available for solar must be numeric")

        try:
            val = float(self.self.rad_2_elec_lineEdit.text()) 
            if val <= 0:
                self.rad_2_elec_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The percent of radiation converted to electricity for solar must be greater than 0%")
            elif val > 100:
                self.rad_2_elec_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The percent of radiation converted to electricity for solar must be less than or equal to 100%")   
        except:
            self.rad_2_elec_lineEdit.setStyleSheet(DEFAULT_ERROR)
            valid = False
            errors.append(f"The percent of radiation converted to electricity for solar must be numeric")

        try:
            val = float(self.self.fract_avg_SolRad_lineEdit.text()) 
            if val <= 0:
                self.fract_avg_SolRad_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The fraction of average solar radiation must be greater than 0")
            elif val > 1:
                self.fract_avg_SolRad_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The fraction of average solar radiation must be less than or equal to 1")   
        except:
            self.fract_avg_SolRad_lineEdit.setStyleSheet(DEFAULT_ERROR)
            valid = False
            errors.append(f"The fraction of average solar radiation must be numeric")
        
        try:
            val = float(self.self.perc_land_used_lineEdit.text()) 
            if val <= 0:
                self.perc_land_used_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The percent of land used for solar must be greater than 0%")
            elif val > 100:
                self.perc_land_used_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The percent of land used for solar must be less than or equal to 100%")   
        except:
            self.perc_land_used_lineEdit.setStyleSheet(DEFAULT_ERROR)
            valid = False
            errors.append(f"The percent of land used for solar must be numeric")

        try:
            val = float(self.self.sun_hours_lineEdit.text()) 
            if val < 1:
                self.sun_hours_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The average number of sun hours per day must be greater than or equal to 1")
            elif val > 24:
                self.sun_hours_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The average number of sun hours per day must be less than or equal to 24")   
        except:
            self.sun_hours_lineEdit.setStyleSheet(DEFAULT_ERROR)
            valid = False
            errors.append(f"The average number of sun hours per day must be numeric")

        try:
            val = float(self.self.dc_ac_lineEdit.text()) 
            if val < 1:
                self.dc_ac_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The D.C. to A.C. ratio must be greater than or equal to 1")
            elif val > 1.35:
                self.dc_ac_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The D.C. to A.C. ratio must be less than or equal to 1.35")   
        except:
            self.dc_ac_lineEdit.setStyleSheet(DEFAULT_ERROR)
            valid = False
            errors.append(f"The D.C. to A.C. ratio must be numeric")

        try:
            val = float(self.self.cap_fac_lineEdit.text()) 
            if val <= 0:
                self.cap_fac_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The solar capacity factor must be greater than 0%")
            elif val > 100:
                self.cap_fac_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The solar capacity factor must be less than or equal to 100%")   
        except:
            self.cap_fac_lineEdit.setStyleSheet(DEFAULT_ERROR)
            valid = False
            errors.append(f"The solar capacity factor must be numeric")  
        # Check Concentrating Solar
        try:
            val = float(self.self.opt_eff_lineEdit.text()) 
            if val <= 0:
                self.opt_eff_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The optical efficiency for concentrating solar must be greater than 0%")
            elif val > 100:
                self.opt_eff_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The optical efficiency for concentrating solar must be less than or equal to 100%")   
        except:
            self.opt_eff_lineEdit.setStyleSheet(DEFAULT_ERROR)
            valid = False
            errors.append(f"The optical efficiency for concentrating solar must be numeric")  

        try:
            val = float(self.self.elec_eff_lineEdit.text()) 
            if val <= 0:
                self.elec_eff_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The electrical efficiency for concentrating solar must be greater than 0%")
            elif val > 100:
                self.elec_eff_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The electrical efficiency for concentrating solar must be less than or equal to 100%")   
        except:
            self.elec_eff_lineEdit.setStyleSheet(DEFAULT_ERROR)
            valid = False
            errors.append(f"The electrical efficiency for concentrating solar must be numeric") 
        
        try:
            val = float(self.self.therm_stor_lineEdit.text()) 
            if val < 0:
                self.therm_stor_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The thermal storage for concentrating solar must be greater than or equal to 0")
            elif val > 30:
                self.therm_stor_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The thermal storage for concentrating solar must be less than or equal to 30")   
        except:
            self.therm_stor_lineEdit.setStyleSheet(DEFAULT_ERROR)
            valid = False
            errors.append(f"The thermal storage for concentrating solar must be numeric")
        
        try:
            val = float(self.self.avg_sunhrs_lineEdit.text()) 
            if val < 0:
                self.avg_sunhrs_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The average hours of sun per day for concentrating solar must be greater than or equal to 0")
            elif val > 30:
                self.avg_sunhrs_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The average hours of sun per day for concentrating solar must be less than or equal to 24")   
        except:
            self.avg_sunhrs_lineEdit.setStyleSheet(DEFAULT_ERROR)
            valid = False
            errors.append(f"The average hours of sun per day for concentrating solar must be numeric")
        
        try:
            val = float(self.self.avg_sol_rad_lineEdit.text()) 
            if val < 0.1:
                self.avg_sol_rad_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The average solar radiation must be greater than or equal to 0.1")
            elif val > 1.115:
                self.avg_sol_rad_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The average solar radiation must be less than or equal to 1.115")   
        except:
            self.avg_sol_rad_lineEdit.setStyleSheet(DEFAULT_ERROR)
            valid = False
            errors.append(f"The average solar radiation must be numeric")
        
        try:
            val = float(self.self.per_agency_land_lineEdit.text()) 
            if val <= 0:
                self.per_agency_land_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The precent of land used for concentrating solar must be greater than 0%")
            elif val > 100:
                self.per_agency_land_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The precent of land used for concentrating solar must be less than or equal to 100%")   
        except:
            self.per_agency_land_lineEdit.setStyleSheet(DEFAULT_ERROR)
            valid = False
            errors.append(f"The precent of land used for concentrating solar must be numeric") 
        		
        try:
            val = float(self.self.cs_cap_fac_lineEdit.text()) 
            if val <= 0:
                self.cs_cap_fac_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The precent of land used for concentrating solar must be greater than 0%")
            elif val > 100:
                self.cs_cap_fac_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The precent of land used for concentrating solar must be less than or equal to 100%")   
        except:
            self.cs_cap_fac_lineEdit.setStyleSheet(DEFAULT_ERROR)
            valid = False
            errors.append(f"The precent of land used for concentrating solar must be numeric") 
        # Check Hydrogen Cells
        try:
            val = float(self.self.h_cap_fac_lineEdit.text()) 
            if val <= 0:
                self.h_cap_fac_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The capacity factor for hydrogen fuel cells must be greater than 0%")
            elif val > 100:
                self.h_cap_fac_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The capacity factor for hydrogen fuel cells must be less than or equal to 100%")   
        except:
            self.h_cap_fac_lineEdit.setStyleSheet(DEFAULT_ERROR)
            valid = False
            errors.append(f"The capacity factor for hydrogen fuel cells must be numeric") 
        
        try:
            val = float(self.self.hrs_oper_lineEdit.text()) 
            if val < 0:
                self.hrs_oper_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The operational hours per day for hydrogen fuel cells must be greater than or equal to 0")
            elif val > 30:
                self.hrs_oper_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The operational hours per day for hydrogen fuel cells must be less than or equal to 24")   
        except:
            self.hrs_oper_lineEdit.setStyleSheet(DEFAULT_ERROR)
            valid = False
            errors.append(f"The operational hours per day for hydrogen fuel cells must be numeric")
        
        try:
            val = float(self.self.cur_den_lineEdit.text()) 
            if val < 0:
                self.cur_den_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The hydrogen fuel cell current density must be greater than or equal to 0")
            elif val > 15000:
                self.cur_den_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The hydrogen fuel cell current density must be less than or equal to 15,000")   
        except:
            self.cur_den_lineEdit.setStyleSheet(DEFAULT_ERROR)
            valid = False
            errors.append(f"The hydrogen fuel cell current density must be numeric")
        
        try:
            val = float(self.self.cell_volt_lineEdit.text()) 
            if val < 0:
                self.cell_volt_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The hydrogen fuel cell active area must be greater than or equal to 0")
            elif val > 5:
                self.cell_volt_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The hydrogen fuel cell active area must be less than or equal to 5")   
        except:
            self.cell_volt_lineEdit.setStyleSheet(DEFAULT_ERROR)
            valid = False
            errors.append(f"The hydrogen fuel cell active area must be numeric")

        try:
            val = float(self.self.cell_act_area_lineEdit.text()) 
            if val < 0.6:
                self.cell_act_area_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The hydrogen fuel cell voltage must be greater than or equal to 0.6")
            elif val > 0.8:
                self.cell_act_area_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The hydrogen fuel cell voltage must be less than or equal to 0.8")   
        except:
            self.cell_act_area_lineEdit.setStyleSheet(DEFAULT_ERROR)
            valid = False
            errors.append(f"The hydrogen fuel cell voltage must be numeric")

        try:
            val = int(self.self.Nstacks_lineEdit.text()) 
            if val < 1:
                self.Nstacks_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The number of stacks for hydrogen fuel cell must be greater than or equal to 1")
            elif val > 100:
                self.Nstacks_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The number of stacks for hydrogen fuel cell must be less than or equal to 100")   
        except:
            self.Nstacks_lineEdit.setStyleSheet(DEFAULT_ERROR)
            valid = False
            errors.append(f"The number of stacks for hydrogen fuel cell must be numeric")
        
        try:
            val = float(self.self.stack_replace_lineEdit.text()) 
            if val < 1:
                self.stack_replace_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The time to replace hydrogen fuel cell must be greater than or equal to 1")
            elif val > 10:
                self.stack_replace_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The time to replace hydrogen fuel cell must be less than or equal to 10")   
        except:
            self.stack_replace_lineEdit.setStyleSheet(DEFAULT_ERROR)
            valid = False
            errors.append(f"The time to replace for hydrogen fuel cell must be numeric")       
        
        try:
            val = float(self.self.annual_degrade_lineEdit.text()) 
            if val < 3:
                self.annual_degrade_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The hydrogen fuel cell degradation must be greater than or equal to 1")
            elif val > 20:
                self.annual_degrade_lineEdit.setStyleSheet(DEFAULT_ERROR)
                valid = False
                errors.append(f"The hydrogen fuel cell degradation must be less than or equal to 10")   
        except:
            self.annual_degrade_lineEdit.setStyleSheet(DEFAULT_ERROR)
            valid = False
            errors.append(f"The hydrogen fuel cell degradation must be numeric") 
        # Check Geothermal
        self.well_diam_lineEdit.text()
        self.fluid_vel_lineEdit.text()
        self.well_depth_lineEdit.text()
        self.turb_outlet_lineEdit.text()
        self.geo_cap_fac_lineEdit.text()
        self.fluid_den_lineEdit.text()
        self.pump_eff_lineEdit.text()
        # Check Goverment Rates
        for entry in [self.wg_disc_rate_lineEdit, self.sg_disc_rate_lineEdit,
                      self.csg_disc_rate_lineEdit, self.hg_disc_rate_lineEdit,
                      self.gg_disc_rate_lineEdit]:
        
        for entry in [self.wg_infl_rate_lineEdit, self.sg_infl_rate_lineEdit,
                      self.csg_infl_rate_lineEdit, self.hg_infl_rate_lineEdit,
                      self.gg_infl_rate_lineEdit]:
        
        for entry in [self.wg_ece_LB_lineEdit, self.sg_ece_LB_lineEdit,
                      self.csg_ece_LB_lineEdit, self.hg_ece_LB_lineEdit,
                      self.gg_ece_LB_lineEdit]:
        
        for entry in [self.wg_ece_UB_lineEdit, self.sg_ece_UB_lineEdit,
                      self.csg_ece_UB_lineEdit, self.hg_ece_UB_lineEdit,
                      self.gg_ece_UB_lineEdit]:
        
        for entry in [self.wg_tax_lineEdit, self.sg_tax_lineEdit,
                      self.csg_tax_lineEdit, self.hg_tax_lineEdit,
                      self.gg_tax_lineEdit]:
        
        for entry in [self.wg_proj_life_lineEdit, self.sg_proj_life_lineEdit,
                      self.csg_proj_life_lineEdit, self.hg_proj_life_lineEdit,
                      self.gg_proj_life_lineEdit]:
        

        # Check Third-Party Rates
        for entry in [self.wt_disc_rate_lineEdit, self.st_disc_rate_lineEdit,
                      self.cst_disc_rate_lineEdit, self.ht_disc_rate_lineEdit,
                      self.gt_disc_rate_lineEdit]:
        
        for entry in [self.wt_infl_rate_lineEdit, self.st_infl_rate_lineEdit,
                      self.cst_infl_rate_lineEdit, self.ht_infl_rate_lineEdit,
                      self.gt_infl_rate_lineEdit]:
        
        for entry in [self.wt_ece_LB_lineEdit, self.st_ece_LB_lineEdit,
                      self.cst_ece_LB_lineEdit, self.ht_ece_LB_lineEdit,
                      self.gt_ece_LB_lineEdit]:
        
        for entry in [self.wt_ece_UB_lineEdit, self.st_ece_UB_lineEdit,
                      self.cst_ece_UB_lineEdit, self.ht_ece_UB_lineEdit,
                      self.gt_ece_UB_lineEdit]:
        
        for entry in [self.wt_tax_lineEdit, self.st_tax_lineEdit
                      self.cst_tax_lineEdit, self.ht_tax_lineEdit,
                      self.gt_tax_lineEdit]:
        
        for entry in [self.wt_proj_life_lineEdit, self.st_proj_life_lineEdit,
                      self.cst_proj_life_lineEdit, self.ht_proj_life_lineEdit,
                      self.gt_proj_life_lineEdit]:
        
        
        # Check the path
        # The only way the user could mess with this is to set the path then delete it before hitting proceed
        if not os.path.exists(str(self.frpp_path.text())):
            self.frpp_path.setStyleSheet(DEFAULT_ERROR)
            valid = False
            errors.append(f"The path specified {self.frpp_path.text()}\nDoes not exist. Please select a valid FRPP dataset to continue.")

        # Check that an agency has been selected
        if self.agency_comboBox.currentIndex() == 0:
            self.agency_comboBox.setStyleSheet(DEFAULT_ERROR)
            valid = False
            errors.append("You must select an agency to continue")

        return valid, errors


if __name__ == "__main__":
    app = QApplication(sys.argv)
    ui = MainWindow()
    ui.show()
    sys.exit(app.exec())