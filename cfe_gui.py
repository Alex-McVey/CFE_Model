from PyQt6.QtWidgets import *
from PyQt6.QtGui import *
from PyQt6.QtCore import Qt, QCoreApplication
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
from datetime import datetime
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

from model_assump_dlg import model_assump_dlg
import rc_cfe

PROJECT_PATH = pathlib.Path(__file__).parent
DATA_PATH = PROJECT_PATH / "data"
PROJECT_UI = PROJECT_PATH / "cfe_model.ui"
ASSUMP_DLG = PROJECT_PATH / "model_assumptions.ui"
DEFAULT_WHITE = u"background-color: rgb(255, 255, 255);"
DEFAULT_ERROR = u"background-color: rgb(255, 103, 103);"


class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)


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

        # Add matplot objects to layouts
        self.cfe_bar_canvas = MplCanvas(self, width=5, height=4, dpi=100)
        self.toolbar_pg3 = NavigationToolbar(self.cfe_bar_canvas, self)
        self.cfe_barchrt_layout.addWidget(self.toolbar_pg3)
        self.cfe_barchrt_layout.addWidget(self.cfe_bar_canvas)
        
        # Establish function connections
        ''' Landing Page '''
        self.frpp_SetNewPath_pushButton.pressed.connect(self.set_newFRPP_path)
        self.load_frpp_pushButton.pressed.connect(self.process_FRPP)
        ''' Build CFE Model Page '''
        self.model_assumpt_pushButton.pressed.connect(self.launch_model_assump_dlg)
        self.build_CFE_pushButton.pressed.connect(self.build_cfe_model)

        # Establish class variables
        self.assmpt_dlf_errors = []
        self.frpp_df = None
        self.cfe_use_df = pd.DataFrame()

        # Set defaults
        self.stackedWidget.setCurrentIndex(0)


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
        assets_wind['Annual Average Wind Speed (m/s)'] = wind_values
        self.frpp_df = self.frpp_df.merge(assets_wind[['Real Property Unique Identifier', 'Annual Average Wind Speed (m/s)']], how = 'left', on = 'Real Property Unique Identifier')

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
        dlg = model_assump_dlg(self.model_assump, color_red=self.assmpt_dlf_errors)
        if dlg.exec():
            self.assmpt_dlf_errors = dlg.color_red.copy()
            self.model_assump = dlg.assump_json.copy()
            if dlg.valid:
                if dlg.overwrite_checkBox.isChecked():
                    json_string = json.dumps(self.model_assump, indent=4)
                    with open(os.path.join(DATA_PATH, "model_assumptions.json"), 'w') as outfile:
                        outfile.write(json_string)
            else:
                self.assmpt_dlf_errors = dlg.color_red.copy()
                msg = QMessageBox(self)
                msg.setIcon(QMessageBox.Icon.Warning)
                msg.setText("\n\n".join(dlg.errors))
                msg.setWindowTitle("Warning")
                msg.setStandardButtons(QMessageBox.StandardButton.Ok)
                msg.exec()
                self.launch_model_assump_dlg()


    def build_cfe_model(self):
        # Check to see if the agency is in the energy data csv
        if not os.path.exists(os.path.join(DATA_PATH, "Agency_Energy_Data.csv")):
            msg = QMessageBox(self)
            msg.setIcon(QMessageBox.Icon.Warning)
            msg.setText(f"The agency energy data csv is missing {os.path.join(DATA_PATH, 'Agency_Energy_Data.csv')}\nThis data is needed to continue")
            msg.setWindowTitle("Warning")
            msg.setStandardButtons(QMessageBox.StandardButton.Ok)
            msg.exec()
            return
        else:
            agency_energy_data = pd.read_csv(os.path.join(DATA_PATH, "Agency_Energy_Data.csv"), header=0)
            if str(self.agency_comboBox.currentText()) not in agency_energy_data['Agency'].values:
                msg = QMessageBox(self)
                msg.setIcon(QMessageBox.Icon.Warning)
                msg.setText(f"The data for the requested agency: {str(self.agency_comboBox.currentText())}\nIs not present in the data set\n\nThe calculation can't continue")
                msg.setWindowTitle("Warning")
                msg.setStandardButtons(QMessageBox.StandardButton.Ok)
                msg.exec()
                return
            else:
                agency_energy_data = agency_energy_data.loc[(agency_energy_data['Agency'] == str(self.agency_comboBox.currentText()))]

        # These variables are just for the report at the end
        n_solar_roof, n_solar_roof_built = 0, 0
        n_solar_grnd, n_solar_grnd_built = 0, 0
        n_wind, n_wind_built = 0, 0
        n_geot, n_geot_built = 0, 0
        n_fuelcell, n_fuelcell_built = 0, 0
        n_conc_sol, n_conc_sol_built = 0, 0

        ''' Solar Power '''
        cur_dict = self.model_assump['solar']['system']
        # Calculate the A.C. Roof-top solar power  
        dmy = np.zeros(self.frpp_df.shape[0])
        dmy[:] = np.nan      
        filter_ = ~self.frpp_df['est_rooftop_area_sqft'].isna()
        n_solar_roof = filter_.sum()
        dmy[filter_] = cur_dict['dc_to_ac_ratio'] *(cur_dict['per_rad_to_elec']/100.) * \
            (0.0929031299*self.frpp_df[filter_]['Annual Solar Radiation kWh/m2/day']) * \
            self.frpp_df[filter_]['est_rooftop_area_sqft'] * (cur_dict['precent_roof_avail']/100.) / \
            cur_dict['avg_sun_hours']
            # 0.0929031299 is m2 to ft2
        self.frpp_df['Rooftop Solar Power'] = dmy.tolist()            
        # Calculate the A.C. Ground mounted Solar
        dmy = np.zeros(self.frpp_df.shape[0])
        dmy[:] = np.nan  
        filter_ = ~self.frpp_df['Acres'].isna()
        n_solar_grnd = filter_.sum()
        dmy[filter_] = cur_dict['dc_to_ac_ratio'] * (cur_dict['per_rad_to_elec']/100.) * (cur_dict['perc_land_used']/100.) * \
            (0.0929031299*self.frpp_df[filter_]['Annual Solar Radiation kWh/m2/day']) * self.frpp_df[filter_]['Acres'] * (43560.) / \
            cur_dict['avg_sun_hours']
        self.frpp_df['Ground Solar Power'] = dmy.tolist()
        # Calculate annual energy rate for solar
        self.frpp_df['Annual Rooftop Solar Power (kWh/yr)'] = self.frpp_df['Rooftop Solar Power'] * 24 * 365 * (cur_dict['capacity_factor']/100.)
        self.frpp_df['Annual Ground Solar Power (kWh/yr)'] = self.frpp_df['Ground Solar Power'] * 24 * 365 * (cur_dict['capacity_factor']/100.)

        ''' Wind Power '''
        cur_dict = self.model_assump['wind']['system']
        # Number of turbines on land
        #                           (*N Defined by user between 5 and 8ish)
        # --->  Primary         X <----  N*Turb Diam -----> X
        #                       -
        # --->  Wind            |
        #                      2*Turb Diam                 
        # --->  Direction       |  
        #                       -
        # --->                  X                           X <- Turbine
        dmy = np.zeros(self.frpp_df.shape[0])
        dmy[:] = np.nan         
        filter_ = ~self.frpp_df['Annual Average Wind Speed (m/s)'].isna()
        dmy[filter_] = np.floor(self.frpp_df[filter_]['Acres'] / (2*cur_dict['turbine_diameter']*\
                                                    cur_dict['turbine_spacing'] * 0.000247105)) # convert to acres
        dmy[np.where(dmy == 0)] = 1 # If it gets rounded down to 0 I think it should still be able to support 1
        self.frpp_df['N Wind Turbines'] = dmy.tolist()
        # Power (W) = (0.5 * Cp * rho * A * v^3) * N_Turbines
        #   - Cp coefficient of performance
        #   - rho density of the air in kg/m3
        #   - A cross-sectional area of the wind in m2
        #   - v velocity of the wind in m/s
        wind_mod = cur_dict['fraction_of_average_wind_speed']
        dmy = np.zeros(self.frpp_df.shape[0])
        dmy[:] = np.nan          
        filter_ = np.logical_and(~self.frpp_df['Annual Average Wind Speed (m/s)'].isna(), \
                  np.logical_or(self.frpp_df['Annual Average Wind Speed (m/s)']*wind_mod > 25 , \
                                self.frpp_df['Annual Average Wind Speed (m/s)']*wind_mod <= 3))  # Cut in/out speeds
        dmy[filter_] = 0.0 
        filter_ = np.logical_and(~self.frpp_df['Annual Average Wind Speed (m/s)'].isna(), \
                  np.logical_and(self.frpp_df['Annual Average Wind Speed (m/s)']*wind_mod <= 25, \
                                 self.frpp_df['Annual Average Wind Speed (m/s)']*wind_mod > 12))  # Power output is constant between these values (ie it cant spin any faster)
        dmy[filter_] = 0.5 * cur_dict['power_coefficient'] * 1.293 * (np.pi * cur_dict['turbine_diameter']**2 / 4.) * \
                                              12**3 * self.frpp_df[filter_]['N Wind Turbines'] * 0.001 
        n_wind = filter_.sum()
        filter_ = np.logical_and(~self.frpp_df['Annual Average Wind Speed (m/s)'].isna(), \
                  np.logical_and(self.frpp_df['Annual Average Wind Speed (m/s)']*wind_mod <= 12, \
                                 self.frpp_df['Annual Average Wind Speed (m/s)']*wind_mod > 3)) # Between these speeds the turbine speed is proportional to wind speed
        dmy[filter_] = 0.5 * cur_dict['power_coefficient'] * 1.293 * (np.pi * cur_dict['turbine_diameter']**2 / 4.) * \
                       (self.frpp_df[filter_]['Annual Average Wind Speed (m/s)']*wind_mod)**3 *\
                       self.frpp_df[filter_]['N Wind Turbines'] * 0.001  
        n_wind += filter_.sum()      
        self.frpp_df['Wind Power (kW)'] = dmy.tolist()
        self.frpp_df['Annual Wind Power (kWh/yr)'] = self.frpp_df['Wind Power (kW)'] * 24 * 365 * (cur_dict['capacity_factor']/100.)
        
        ''' Concentrating Solar '''
        cur_dict = self.model_assump['conc_solar']['system']
        # Explain
        dmy = np.zeros(self.frpp_df.shape[0])
        dmy[:] = np.nan 
        filter_ = np.logical_and(self.frpp_df['Real Property Type'] == 'Land', self.frpp_df['Real Property Use'] == 'Vacant')
        n_conc_sol = filter_.sum()
        # calculate
        self.frpp_df['Concentrating Solar Power (kW)'] = dmy.tolist()
        self.frpp_df['Annual Concentrating Solar Power (kWh)'] = self.frpp_df['Concentrating Solar Power (kW)'] * 24 * 365 * (cur_dict['capacity_factor']/100.)

        ''' GeoThermal '''
        cur_dict = self.model_assump['geo_therm']['system']
        dmy = np.zeros(self.frpp_df.shape[0])
        dmy[:] = np.nan  
        filter_ = np.logical_and(~self.frpp_df['Geothermal_CLASS'].isna(), self.frpp_df['Geothermal_CLASS'] <=3)
        n_geot = filter_.sum()
        enthalpy_diff, inlet_temp, outlet_temp = self.geothermal_enthalpy(self.frpp_df[filter_]['Geothermal_CLASS']) 
        dmy[filter_] = (((np.pi*cur_dict['well_diameter']**2/4)*cur_dict['fluid_velocity']*cur_dict['fluid_density'])*\
                         enthalpy_diff* (1-(outlet_temp+273)/(inlet_temp+273))) - \
                         ((((np.pi*cur_dict['well_diameter']**2/4)*cur_dict['fluid_velocity']*cur_dict['fluid_density']\
                         *9.81 *cur_dict['well_depth']/(cur_dict['overall_pump_efficiency']/100))*2)/1000)
        self.frpp_df['Geothermal Power (kW)'] = dmy.tolist()
        self.frpp_df['Annual Geothermal Power (kWh/yr)'] = self.frpp_df['Geothermal Power (kW)'] * 24 * 365 * (cur_dict['capacity_factor']/100.)

        ''' Hydrogen Fuel Cells '''
        cur_dict = self.model_assump['hydrogen']['system']
        self.frpp_df['Fuel Cell (kW)'] = (cur_dict['curr_density'] * cur_dict['cell_active_area'] * \
                                          cur_dict['cell_voltage']) * cur_dict['N_stacks'] * 0.001
        self.frpp_df['Annual Fuel Cell (kW/yr)'] =  self.frpp_df['Fuel Cell (kW)'] * cur_dict['daily_operation'] * \
            365. * (cur_dict['capacity_factor']/100.)
        n_fuelcell = self.frpp_df.shape[0]

        ''' Finally Sum the various renewables to find the total power per site '''
        # The following stipulations are used when calculating renewable potential (property use listed)
        # Rooftop Solar
        #   - buildings 
        # Ground mounted solar
        #   - Vacant 
        #   - R&D
        #   - Office land
        # Wind
        #   - Vacant 
        #   - R&D
        # Concentrating Solar
        #   - Vacant
        # Geothermal
        #   - Vacant
        # Fuel Cell
        #   - No restrictions
        #
        # This means that for vacant land we will potentiallys have data for Wind, Ground Solar, Concentrating Solar, and Geothermal 
        # in this instance we will choose the one with the greatest power potential (though in the future that should
        # probably be weighed against its cost)
        power = np.array([self.frpp_df['Ground Solar Power'].to_numpy(),
                         self.frpp_df['Concentrating Solar Power (kW)'].to_numpy(), 
                         self.frpp_df['Wind Power (kW)'].to_numpy(), 
                         self.frpp_df['Geothermal Power (kW)'].to_numpy() ]).transpose()
        energy = np.array([self.frpp_df['Annual Ground Solar Power (kWh/yr)'].to_numpy(), 
                          self.frpp_df['Annual Concentrating Solar Power (kWh)'].to_numpy(), 
                          self.frpp_df['Annual Wind Power (kWh/yr)'].to_numpy(), 
                          self.frpp_df['Annual Geothermal Power (kWh/yr)'].to_numpy()]).transpose()
        # Vacant Land
        filter_ = np.logical_and(self.frpp_df['Real Property Type'] == 'Land', self.frpp_df['Real Property Use'] == 'Vacant')
        indicies = np.where(filter_)[0]
        for qq in indicies:
            if not np.isnan(power[qq,:]).all():
                mx_val = np.nanargmax(power[qq,:])
                chng_ind = [x for x in [0, 1, 2, 3] if x != mx_val]
                for kk in chng_ind:
                    power[qq,kk] = 0.0
                    energy[qq,kk] = 0.0
        # R&D Land
        filter_ = np.logical_and(self.frpp_df['Real Property Type'] == 'Land', self.frpp_df['Real Property Use'] == 'Research and Development')
        indicies = np.where(filter_)[0]
        for qq in indicies:
            if not np.isnan(power[qq,:]).all():
                mx_val = np.nanargmax(power[qq,:])
                chng_ind = [x for x in [0, 1, 2, 3] if x != mx_val]
                for kk in chng_ind:
                    power[qq,kk] = 0.0
                    energy[qq,kk] = 0.0
        # Office Land
        filter_ = np.logical_and(self.frpp_df['Real Property Type'] == 'Land', self.frpp_df['Real Property Use'] == 'Office Building Locations')
        indicies = np.where(filter_)[0]
        for qq in indicies:
            if not np.isnan(power[qq,:]).all():
                mx_val = np.nanargmax(power[qq,:])
                chng_ind = [x for x in [0, 1, 2, 3] if x != mx_val]
                for kk in chng_ind:
                    power[qq,kk] = 0.0
                    energy[qq,kk] = 0.0


        power = power.transpose()
        energy = energy.transpose()
        self.frpp_df['Ground Solar Power'] = power[0,:]
        self.frpp_df['Concentrating Solar Power (kW)'] = power[1,:]
        self.frpp_df['Wind Power (kW)'] = power[2,:]
        self.frpp_df['Geothermal Power (kW)'] = power[3,:]
        self.frpp_df['Annual Ground Solar Power (kWh/yr)'] = energy[0,:]
        self.frpp_df['Annual Concentrating Solar Power (kWh)'] = energy[1,:]
        self.frpp_df['Annual Wind Power (kWh/yr)'] = energy[2,:]
        self.frpp_df['Annual Geothermal Power (kWh/yr)'] = energy[3,:]
        n_solar_roof_built = np.sum(~self.frpp_df['Rooftop Solar Power'].isna().values)
        n_solar_grnd_built = np.sum(np.logical_and(~self.frpp_df['Ground Solar Power'].isna().values,
                                                   self.frpp_df['Ground Solar Power'] != 0))
        n_wind_built = np.sum(np.logical_and(~self.frpp_df['Wind Power (kW)'].isna().values,
                                             self.frpp_df['Wind Power (kW)'] != 0))
        n_geot_built = np.sum(np.logical_and(~self.frpp_df['Geothermal Power (kW)'].isna().values,
                                             self.frpp_df['Geothermal Power (kW)'] != 0))
        n_fuelcell_built = np.sum(~self.frpp_df['Fuel Cell (kW)'].isna().values)
        n_conc_sol_built = np.sum(np.logical_and(~self.frpp_df['Concentrating Solar Power (kW)'].isna().values,
                                                 self.frpp_df['Concentrating Solar Power (kW)'] != 0))

        # Finalize
        self.frpp_df['Total Power (kW)'] = self.frpp_df[['Rooftop Solar Power','Ground Solar Power',
                                                         'Wind Power (kW)','Fuel Cell (kW)','Geothermal Power (kW)',
                                                         'Concentrating Solar Power (kW)']].sum(axis=1)
           
        self.frpp_df['Total Energy (kWh)'] = self.frpp_df[['Annual Rooftop Solar Power (kWh/yr)','Annual Ground Solar Power (kWh/yr)',
                                                          'Annual Wind Power (kWh/yr)','Annual Fuel Cell (kW/yr)',
                                                          'Annual Geothermal Power (kWh/yr)','Annual Concentrating Solar Power (kWh)']].sum(axis=1)
        

        ''' Build CFE Energy Usage and Transition Table '''
        cur_year = datetime.now().year
        year_ind = list(range(cur_year,2036))
        
        pp = -(1000*agency_energy_data['Electricity (MWh)'].values[0] * (float(self.energy_proj_lineEdit.text())+100)) /\
            (float(self.oper_days_lineEdit.text())*(float(self.CB_oper_lineEdit.text())*\
            (float(self.cur_cfe_lineEdit.text())-100) - float(self.cur_cfe_lineEdit.text()) * \
            float(self.ren_oper_lineEdit.text())))

        self.cfe_use_df['year'] = year_ind
        # Power requirement with annual percent energy growth
        dmy = [pp]
        for year in year_ind[1:]:
            dmy.append(dmy[0]*(1. + float(self.AEG_lineEdit.text())/100.)**(year-2023))
        self.cfe_use_df['Power in year N (kW)'] = dmy
         
        # Carbon based vs Renewable Power percent per year
        cfe = [dmy[0]*(float(self.cur_cfe_lineEdit.text())/100.)]
        cb = [dmy[0]-cfe[0]]
        ann_dec_in_carbon_based = cb[0]/(2030-cur_year)
        for qq in range(len(year_ind[1:])):
            cb.append(max(cb[qq]-ann_dec_in_carbon_based,0.))
            cfe.append(dmy[qq+1]-cb[-1])
        self.cfe_use_df['Carbon-Based Power in year N (kW)'] = cb
        self.cfe_use_df['Renewable Power in year N (kW)'] = cfe

        # Carbon based vs Renewable Energy percent per year
        self.cfe_use_df['Carbon-Based Energy in year N (kWh)'] = self.cfe_use_df['Carbon-Based Power in year N (kW)'] *\
              float(self.CB_oper_lineEdit.text()) * float(self.oper_days_lineEdit.text())
        self.cfe_use_df['Renewable Energy in year N (kWh)'] = self.cfe_use_df['Renewable Power in year N (kW)'] *\
              float(self.ren_oper_lineEdit.text()) * float(self.oper_days_lineEdit.text())
        
        # Total Energy in year N
        self.cfe_use_df['Total Energy in year N (MWh)'] = (self.cfe_use_df['Renewable Energy in year N (kWh)'] + \
            self.cfe_use_df['Carbon-Based Energy in year N (kWh)']) * 0.001
        dmy = []
        for year in year_ind:
            dmy.append(agency_energy_data['Electricity (MWh)'].values[0] * \
                (1. + float(self.energy_proj_lineEdit.text())/100.) * \
                (1 + float(self.AEG_lineEdit.text())/100.)**(year-2023))
        self.cfe_use_df['Projected Total Energy in year N (MWh)'] = dmy

        ''' Populate the next Page '''        
        self.cfe_summary_textEdit.setHtml(QCoreApplication.translate("MainWindow", u"<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
            "<html><head><meta name=\"qrichtext\" content=\"1\" /><meta charset=\"utf-8\" /><style type=\"text/css\">\n"
            "p, li { white-space: pre-wrap; }\n"
            "hr { height: 1px; border-width: 0; }\n"
            "li.unchecked::marker { content: \"\\2610\"; }\n"
            "li.checked::marker { content: \"\\2612\"; }\n"
            "</style></head><body style=\" font-family:'Curier'; font-size:9pt; font-weight:400; font-style:normal;\">\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:24pt; font-weight:700; text-decoration: underline; color:#000086;\">CFE Potential for {str(self.agency_comboBox.currentText())}</span></p>\n"
            "<p style=\"-qt-paragraph-type:empty;  -qt-block-indent:0; text-indent:0px; font-size:12pt; color:#000086;\"><br /></p>\n"
            "<p style=\"  -qt-block-indent"
                                    f":0; text-indent:0px;\"><span style=\" font-size:14pt;\">Total Energy Required: 		{agency_energy_data['Electricity (MWh)'].values[0]:,.0f} MWh</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt;\">Total Estimated Production Capacity:	{self.frpp_df['Total Energy (kWh)'].sum()*0.001:,.0f} MWh</span></p>\n"
            "<p style=\"-qt-paragraph-type:empty;  -qt-block-indent:0; text-indent:0px; font-size:14pt;\"><br /></p>\n"
            "<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:22pt; color:#990000;\">Wind</span></p>\n"
            "<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:10pt; color:#990000;\">===================================================</span></p>\n"
            "<p style=\"  -qt-block-i"
                                    f"ndent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Number of Sites avaialable:{n_wind:>12}</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Number of sites built:{n_wind_built:>22}</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Total Potential:{self.frpp_df['Wind Power (kW)'].sum()*0.001:>33,.0f} MW</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Total Annual Generation:{self.frpp_df['Annual Wind Power (kWh/yr)'].sum()*0.001:>17,.0f} MWh</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Percent of Agency Demand:{(self.frpp_df['Annual Wind Power (kWh/yr)'].sum()*0.001)/agency_energy_data['Electricity (MWh)'].values[0]*100:>10.1f}%</span></p>\n"
            "<p style=\"-qt-paragraph-type:empty;  -qt-block-indent:0; text-indent:0px;  font-size:14pt;\"><br /></p>\n"
            "<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:22pt; color:#990000;\">Rooftop Solar PV</span></p>\n"
            "<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:10pt; color:#990000;\">===================================================</span></p>\n"
            "<p style=\"  -qt-block-i"
                                    f"ndent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Number of Sites avaialable:{n_solar_roof:>15}</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Number of sites built:{n_solar_roof_built:>26}</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Total Potential:{self.frpp_df['Rooftop Solar Power'].sum()*0.001:>37,.0f} MW</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Total Annual Generation:{self.frpp_df['Annual Rooftop Solar Power (kWh/yr)'].sum()*0.001:>19,.0f} MWh</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Percent of Agency Demand:{(self.frpp_df['Annual Rooftop Solar Power (kWh/yr)'].sum()*0.001)/agency_energy_data['Electricity (MWh)'].values[0]*100:>17.1f}%</span></p>\n"
            "<p style=\"-qt-paragraph-type:empty;  -qt-block-indent:0; text-indent:0px; font-size:14pt;\"><br /></p>\n"
            "<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:22pt; color:#990000;\">Ground Mounted Solar PV</span></p>\n"
            "<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:10pt; color:#990000;\">===================================================</span></p>\n"
            "<p style=\"  -qt-block-i"
                                    f"ndent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Number of Sites avaialable:{n_solar_grnd:>16}</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Number of sites built{n_solar_grnd_built:>27}</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Total Potential:{self.frpp_df['Ground Solar Power'].sum()*0.001:>36,.0f} MW</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Total Annual Generation:{self.frpp_df['Annual Ground Solar Power (kWh/yr)'].sum()*0.001:>17,.0f} MWh</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Percent of Agency Demand:{(self.frpp_df['Annual Ground Solar Power (kWh/yr)'].sum()*0.001)/agency_energy_data['Electricity (MWh)'].values[0]*100:>17.1f}%</span></p>\n" 
            "<p style=\"-qt-paragraph-type:empty;  -qt-block-indent:0; text-indent:0px; font-size:14pt;\"><br /></p>\n"
            "<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:22pt; color:#990000;\">Hydrogen Fuel Cell</span></p>\n"
            "<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:10pt; color:#990000;\">===================================================</span></p>\n"
            "<p style=\"  -qt-block-i"
                                    f"ndent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Number of Sites avaialable:{n_fuelcell:>17}</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Number of sites built:{n_fuelcell_built:>27}</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Total Potential:{self.frpp_df['Fuel Cell (kW)'].sum()*0.001:>38,.0f} MW</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Total Annual Generation:{self.frpp_df['Annual Fuel Cell (kW/yr)'].sum()*0.001:>19,.0f} MWh</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Percent of Agency Demand:{(self.frpp_df['Annual Fuel Cell (kW/yr)'].sum()*0.001)/agency_energy_data['Electricity (MWh)'].values[0]*100:>17.1f}%</span></p>\n"
            "<p style=\"-qt-paragraph-type:empty;  -qt-block-indent:0; text-indent:0px; font-size:14pt;\"><br /></p>\n"
            "<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:22pt; color:#990000;\">Geothermal Power</span></p>\n"
            "<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:10pt; color:#990000;\">===================================================</span></p>\n"
            "<p style=\"  -qt-block-i"
                                    f"ndent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Number of Sites avaialable:{n_geot:>14}</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Number of sites built:{n_geot_built:>24}</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Total Potential:{self.frpp_df['Geothermal Power (kW)'].sum()*0.001:>35,.0f} MW</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Total Annual Generation:{self.frpp_df['Annual Geothermal Power (kWh/yr)'].sum()*0.001:>19,.0f} MWh</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Percent of Agency Demand:{(self.frpp_df['Annual Geothermal Power (kWh/yr)'].sum()*0.001)/agency_energy_data['Electricity (MWh)'].values[0]*100:>13.1f}%</span></p>\n"
            "<p style=\"-qt-paragraph-type:em"
                                    "pty;  -qt-block-indent:0; text-indent:0px; font-size:14pt; color:#000000;\"><br /></p>\n"
            "<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:22pt; color:#990000;\">Concentrating Solar</span></p>\n"
            "<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:10pt; color:#990000;\">===================================================</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Number of Sites avaialable:{n_conc_sol:>14}</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Number of sites built:{n_conc_sol_built:>25}</span></p>\n"
            "<p "
                                    f"style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Total Potential:{self.frpp_df['Concentrating Solar Power (kW)'].sum()*0.001:>35,.0f} MW</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Total Annual Generation:{self.frpp_df['Annual Concentrating Solar Power (kWh)'].sum()*0.001:>19,.0f} MWh</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Percent of Agency Demand:{(self.frpp_df['Annual Concentrating Solar Power (kWh)'].sum()*0.001)/agency_energy_data['Electricity (MWh)'].values[0]*100:>13.1f}%</span></p></body></html>", None))
        
        bar_data = {"Carbon-Based Energy": self.cfe_use_df['Carbon-Based Energy in year N (kWh)'].values,
                    "Renewable Energy": self.cfe_use_df['Renewable Energy in year N (kWh)'].values}
        x = np.arange(len(year_ind))  # the label locations
        width = 0.25  # the width of the bars
        multiplier = 0

        for attribute, measurement in bar_data.items():
            offset = width * multiplier
            rects = self.cfe_bar_canvas.axes.bar(x + offset, measurement, width, label=attribute)
            multiplier += 1

        # Add some text for labels, title and custom x-axis tick labels, etc.
        self.cfe_bar_canvas.axes.set_ylabel('Energy (MWh)')
        self.cfe_bar_canvas.axes.set_title('Energy Source by Fiscal Year')
        self.cfe_bar_canvas.axes.set_xticks(x + width, year_ind)
        self.cfe_bar_canvas.axes.legend(loc='upper right', ncols=2)

        self.stackedWidget.setCurrentIndex(2)


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
    
    def geothermal_enthalpy(self, geo_class:np.array)->np.array:
        '''
        Takes in an array of geothermal classes and
        interpolates the inlet and outlet enthaly and 
        returns delta(enthalpy).

        Parameters
        ----------
        geo_class: np.array
            - geothermal class

        Returns
        -------
        delta_h: np.array
            - the difference between inlet and outlet enthalpy
        '''
        well_temp = {1:300, 2:200, 3:160, 4:0, 5:0, 999:0}
        enth_dict = {50:209, 60:251, 70:293, 80:334.9, 90:376.9, 100:419, 150:2746, 160:2757.6, 200:2804, 300:2749}

        delta_h = np.zeros_like(geo_class)
        out_enth = np.zeros_like(geo_class)
        in_enth = np.zeros_like(geo_class)
        in_temp = np.zeros_like(geo_class)
        
        out_temp = self.model_assump['geo_therm']['system']['turb_outlet_temp']
        if out_temp in enth_dict.keys():
            out_enth[:] = enth_dict[out_temp]
        else:
            # linear interplolate the enthalpy (outlet is bound between 50 and 100)
            out_enth[:] = np.interp(out_temp, list(enth_dict.keys()), list(enth_dict.values()))

        for qq,val in enumerate(geo_class):
            in_enth[qq] = enth_dict[well_temp[val]]
            in_temp[qq] = well_temp[val]

        delta_h = in_enth - out_enth
        out_temp =  np.ones_like(geo_class)
        out_temp[:] = self.model_assump['geo_therm']['system']['turb_outlet_temp']
        return delta_h, in_temp, out_temp
        



if __name__ == "__main__":
    app = QApplication(sys.argv)
    ui = MainWindow()
    ui.show()
    sys.exit(app.exec())