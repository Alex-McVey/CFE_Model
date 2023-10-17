from PyQt6.QtWidgets import *
from PyQt6.QtGui import *
from PyQt6.QtCore import *#Qt, QCoreApplication, pyqtSignal
from PyQt6.uic import loadUi

import sys
import pathlib
import json
import os
from solar import solar_pv
from geothermal import getem_routine
import pandas as pd
import numpy as np
import numpy_financial as npf
import geopandas as gpd
import rasterio
from pyproj import Transformer
import time
from datetime import datetime
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import logging
# import timeit
# import debugpy

from model_assump_dlg import model_assump_dlg
import rc_cfe


PROJECT_PATH = pathlib.Path(__file__).parent
DATA_PATH = PROJECT_PATH / "data"
PROJECT_UI = PROJECT_PATH / "cfe_model.ui"
ASSUMP_DLG = PROJECT_PATH / "model_assumptions.ui"
DEFAULT_WHITE = u"background-color: rgb(255, 255, 255);"
DEFAULT_ERROR = u"background-color: rgb(255, 103, 103);"
logging.basicConfig(format="%(message)s", level=logging.INFO)

class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        super(MplCanvas, self).__init__(self.fig)    


class buildCFEWorker(QObject):
    progress_signal = pyqtSignal(int, str)
    result_signal = pyqtSignal(object)  

    def __init__(self, frpp_df, agency_energy_data, agency_price_data, model_assump):
        super().__init__()
        self.frpp_df = frpp_df
        self.agency_energy_data = agency_energy_data
        self.agency_price_data = agency_price_data
        self.model_assump = model_assump
        self.count = 0
        self.solar_count = 0
        self.dmy = None
        self.dmy2 = None

    def build_cfe_model(self):
        # debugpy.debug_this_thread()        
        ''' Solar Power '''        
        cur_dict = self.model_assump['solar']['system']        
        # Calculate the A.C. Roof-top solar power  
        dmy = np.zeros(self.frpp_df.shape[0])
        dmy[:] = np.nan     
        dmy2 = np.zeros(self.frpp_df.shape[0])
        dmy2[:] = np.nan          
        filter_ = ~self.frpp_df['est_rooftop_area_sqft'].isna()
        solar_count = filter_.sum()
        count = 0                
        for qq in np.where(filter_)[0]: 
            if cur_dict['use_lat_tilt']:
                solar_calc = solar_pv(self.frpp_df.iloc[qq]['Latitude'], self.frpp_df.iloc[qq]['Longitude'], 
                                    self.frpp_df.iloc[qq]['est_rooftop_area_sqft'], system_loss=cur_dict['system_losses'],
                                    dc_ac=cur_dict['dc_to_ac_ratio'], invert_eff=cur_dict['invert_eff'],
                                    per_area=cur_dict['precent_roof_avail'], tilt=self.frpp_df.iloc[qq]['Latitude'],
                                    azimuth=cur_dict['azimuth'], mod_type=cur_dict['module_type'])       
            else:
                solar_calc = solar_pv(self.frpp_df.iloc[qq]['Latitude'], self.frpp_df.iloc[qq]['Longitude'], 
                                    self.frpp_df.iloc[qq]['est_rooftop_area_sqft'], system_loss=cur_dict['system_losses'],
                                    dc_ac=cur_dict['dc_to_ac_ratio'], invert_eff=cur_dict['invert_eff'],
                                    per_area=cur_dict['precent_roof_avail'], tilt=cur_dict['tilt'],
                                    azimuth=cur_dict['azimuth'], mod_type=cur_dict['module_type'])
            solar_calc.analyze()
            dmy[qq] = solar_calc.total
            dmy2[qq] = solar_calc.dc_nameplate
            count += 1 
            self.progress_signal.emit(int(count/solar_count*100), "Rooftop Solar")  
        self.frpp_df['Annual Rooftop Solar Power (kWh/yr)'] = dmy.tolist()
        self.frpp_df['Rooftop Solar Power'] = dmy2.tolist()
                         
        # Calculate the A.C. Ground mounted Solar
        dmy = np.zeros(self.frpp_df.shape[0])
        dmy[:] = np.nan  
        dmy2 = np.zeros(self.frpp_df.shape[0])
        dmy2[:] = np.nan          
        filter_ = ~self.frpp_df['Acres'].isna()
        solar_count = filter_.sum()
        count = 0
        for qq in np.where(filter_)[0]:            
            if cur_dict['use_lat_tilt']:
                solar_calc = solar_pv(self.frpp_df.iloc[qq]['Latitude'], self.frpp_df.iloc[qq]['Longitude'], 
                                    self.frpp_df.iloc[qq]['Acres']*43560, system_loss=cur_dict['system_losses'],
                                    dc_ac=cur_dict['dc_to_ac_ratio'], invert_eff=cur_dict['invert_eff'],
                                    per_area=cur_dict['perc_land_used'], tilt=self.frpp_df.iloc[qq]['Latitude'],
                                    azimuth=cur_dict['azimuth'], mod_type=cur_dict['module_type'], ground=True)       
            else:
                solar_calc = solar_pv(self.frpp_df.iloc[qq]['Latitude'], self.frpp_df.iloc[qq]['Longitude'], 
                                    self.frpp_df.iloc[qq]['Acres']*43560, system_loss=cur_dict['system_losses'],
                                    dc_ac=cur_dict['dc_to_ac_ratio'], invert_eff=cur_dict['invert_eff'],
                                    per_area=cur_dict['perc_land_used'], tilt=cur_dict['tilt'],
                                    azimuth=cur_dict['azimuth'], mod_type=cur_dict['module_type'], ground=True)
            solar_calc.analyze()
            dmy[qq] = solar_calc.total
            dmy2[qq] = solar_calc.dc_nameplate
            count += 1 
            self.progress_signal.emit(int(count/solar_count*100), "Ground Mounted Solar") 

        self.frpp_df['Ground Solar Power'] = dmy2.tolist()        
        self.frpp_df['Annual Ground Solar Power (kWh/yr)'] = dmy.tolist()

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
        
        # https://windexchange.energy.gov/maps-data/332
        cap_factor_by_state = { "ALABAMA" : cur_dict['capacity_factor'], "ALASKA" : 28.6, "ARIZONA" : 24.3,
                                "ARKANSAS" : cur_dict['capacity_factor'], "CALIFORNIA" : 26.8, "COLORADO" : 35.6,
                                "CONNECTICUT" : cur_dict['capacity_factor'], "DELAWARE" : cur_dict['capacity_factor'], "DISTRICT OF COLUMBIA": cur_dict['capacity_factor'],
                                "FLORIDA" : cur_dict['capacity_factor'], "GEORGIA" : cur_dict['capacity_factor'], "HAWAII" : 30.5,
                                "IDAHO" : 30.1, "ILLINOIS" : 34.5, "INDIANA" : 29.7,
                                "IOWA" : 35.3, "KANSAS" : 42.5, "KENTUCKY" : cur_dict['capacity_factor'],
                                "LOUISIANA" : cur_dict['capacity_factor'], "MAINE" : 29.4, "MARYLAND" : 33.7,
                                "MASSACHUSETTS" : 28.7, "MICHIGAN" : 33.1, "MINNESOTA" : 35.9,
                                "MISSISSIPPI" : cur_dict['capacity_factor'], "MISSOURI" : 31.7, "MONTANA" : 35.5,
                                "NEBRASKA" : 43.8, "NEVADA" : 27.2, "NEW HAMPSHIRE" : 25.4,
                                "NEW JERSEY" : cur_dict['capacity_factor'], "NEW MEXICO" : 37.6, "NEW YORK" : 25.7,
                                "NORTH CAROLINA" : cur_dict['capacity_factor'], "NORTH DAKOTA" : 42.9, "OHIO" : 33.5,
                                "OKLAHOMA" : 41.1, "OREGON" : 22.1, "PENNSYLVANIA" : 29.9,
                                "RHODE ISLAND" : cur_dict['capacity_factor'], "SOUTH CAROLINA" : cur_dict['capacity_factor'],
                                "SOUTH DAKOTA" : 39.6, "TENNESSEE" : 18.3, "TEXAS" : 36.0,
                                "UTAH" : 25.2, "VERMONT" : 28.5, "VIRGINIA" : cur_dict['capacity_factor'],
                                "WASHINGTON" : 25.8, "WEST VIRGINIA" : 28.0, "WISCONSIN" : 27.4, "WYOMING" : 33.2}

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
        capFac = np.zeros(self.frpp_df.shape[0])
        capFac[:] = [cap_factor_by_state[state]/100. if type(state) == str else cur_dict['capacity_factor']/100. for state in self.frpp_df['State Name'].values]
        self.frpp_df['Annual Wind Power (kWh/yr)'] = self.frpp_df['Wind Power (kW)'] * 24 * 365 * capFac
        
        ''' Concentrating Solar '''
        cur_dict = self.model_assump['conc_solar']['system']
        # Explain
        dmy = np.zeros(self.frpp_df.shape[0])
        dmy[:] = np.nan 
        filter_ = np.logical_and(self.frpp_df['Real Property Type'] == 'Land', self.frpp_df['Real Property Use'] == 'Vacant')
        dmy[filter_] = (1.0 *(1000/cur_dict['avg_sun_hours'])*((cur_dict['perc_land_used']/100)*self.frpp_df[filter_]['Acres'] * 43560.))*\
            (cur_dict['optical_eff']/100)*(cur_dict['elec_eff']/100)/1000.  # The 1.0 is 1 kWh/ft2/day of radiation. No idea why
        # calculate
        self.frpp_df['Concentrating Solar Power (kW)'] = dmy.tolist()
        self.frpp_df['Annual Concentrating Solar Power (kWh)'] = self.frpp_df['Concentrating Solar Power (kW)'] * 24 * 365 * (cur_dict['capacity_factor']/100.)

        ''' GeoThermal '''        
        cur_dict = self.model_assump['geo_therm']['system']
        # There are three possible values. calculate all three then just assign
        geoth_power = {1: getem_routine(1,cur_dict),
                       2: getem_routine(2,cur_dict),
                       3: getem_routine(3,cur_dict)}
        dmy = np.zeros(self.frpp_df.shape[0])
        dmy[:] = np.nan  
        filter_ = np.logical_and(~self.frpp_df['Geothermal_CLASS'].isna(), self.frpp_df['Geothermal_CLASS'] <=3)  
        for qq in np.where(filter_)[0]:
            dmy[qq] = geoth_power[self.frpp_df.iloc[qq]['Geothermal_CLASS']]              
        self.frpp_df['Geothermal Power (kW)'] = dmy.tolist()
        self.frpp_df['Annual Geothermal Power (kWh/yr)'] = self.frpp_df['Geothermal Power (kW)'] * 24 * 365 * (cur_dict['capacity_factor']/100.)

        ''' Hydrogen Fuel Cells '''
        cur_dict = self.model_assump['hydrogen']['system']
        self.frpp_df['Fuel Cell (kW)'] = (cur_dict['curr_density'] * cur_dict['cell_active_area'] * \
                                          cur_dict['cell_voltage']) * cur_dict['N_stacks'] * 0.001
        self.frpp_df['Annual Fuel Cell (kW/yr)'] =  self.frpp_df['Fuel Cell (kW)'] * cur_dict['daily_operation'] * \
            365. * (cur_dict['capacity_factor']/100.)
                
        self.result_signal.emit(self.frpp_df)    


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

        self.econ_canvas = MplCanvas(self, width=8, height=4, dpi=100)
        self.econ_toolbar = NavigationToolbar(self.econ_canvas, self)
        self.econ_NPV_layout.addWidget(self.econ_toolbar)
        self.econ_NPV_layout.addWidget(self.econ_canvas)
        
        # Establish function connections
        ''' Menu - File '''
        self.actionNew.triggered.connect(self.reset)
        self.actionLoad.triggered.connect(self.load)
        self.actionExport_Results.triggered.connect(self.export)
        ''' Landing Page '''
        self.frpp_SetNewPath_pushButton.pressed.connect(self.set_newFRPP_path)
        self.load_frpp_pushButton.pressed.connect(self.process_FRPP)
        ''' Build CFE Model Page '''
        self.model_assumpt_pushButton.pressed.connect(self.launch_model_assump_dlg)
        self.build_CFE_pushButton.pressed.connect(self.build_cfe_model)
        ''' Econ Tab '''
        self.econ_cfe_comboBox.currentIndexChanged.connect(self.econ_val_changed)
        self.gov_radioButton.toggled.connect(self.owner_switch)
        self.third_party_radioButton.toggled.connect(self.owner_switch)
        self.econ_energy_storage_lineEdit.editingFinished.connect(self.econ_val_changed)
        self.econ_bat_rep_lineEdit.editingFinished.connect(self.econ_val_changed)
        self.econ_agency_cons_lineEdit.editingFinished.connect(self.econ_val_changed)
        self.econ_agency_elec_cost_lineEdit.editingFinished.connect(self.econ_val_changed)
        self.econ_interest_rate_lineEdit.editingFinished.connect(self.econ_val_changed)
        self.econ_proj_life_lineEdit.editingFinished.connect(self.econ_val_changed)
        self.econ_marcs_comboBox.currentIndexChanged.connect(self.econ_val_changed)
        self.econ_ann_deg_lineEdit.editingFinished.connect(self.econ_val_changed)
        self.econ_rec_lineEdit.editingFinished.connect(self.econ_val_changed)
        self.econ_cont_per_lineEdit.editingFinished.connect(self.econ_val_changed)
        self.econ_tax_rate_lineEdit.editingFinished.connect(self.econ_val_changed)
        self.econ_insur_rate_lineEdit.editingFinished.connect(self.econ_val_changed)
        self.econ_inflat_rate_lineEdit.editingFinished.connect(self.econ_val_changed)
        self.econ_fed_tax_cred_lineEdit.editingFinished.connect(self.econ_val_changed)
        self.econ_tax_cred_per_lineEdit.editingFinished.connect(self.econ_val_changed)
        self.econ_price_esc_lineEdit.editingFinished.connect(self.econ_val_changed)
        self.econ_cont_per_Slider.valueChanged.connect(self.contract_period_slider)
        self.econ_tax_rate_Slider.valueChanged.connect(self.tax_rate_slider)
        self.econ_insur_rate_Slider.valueChanged.connect(self.insurance_rate_slider)
        self.econ_inflat_rate_Slider.valueChanged.connect(self.inflation_rate_slider)
        self.econ_fed_tax_cred_Slider.valueChanged.connect(self.fed_tax_credit_slider)
        self.econ_tax_cred_per_Slider.valueChanged.connect(self.tax_cred_period_slider)
        self.econ_price_esc_Slider.valueChanged.connect(self.price_escalat_slider)
        self.econ_cont_per_Slider.sliderReleased.connect(self.econ_val_changed)
        self.econ_tax_rate_Slider.sliderReleased.connect(self.econ_val_changed)
        self.econ_insur_rate_Slider.sliderReleased.connect(self.econ_val_changed)
        self.econ_inflat_rate_Slider.sliderReleased.connect(self.econ_val_changed)
        self.econ_fed_tax_cred_Slider.sliderReleased.connect(self.econ_val_changed)
        self.econ_tax_cred_per_Slider.sliderReleased.connect(self.econ_val_changed)
        self.econ_price_esc_Slider.sliderReleased.connect(self.econ_val_changed)        

        # Establish class variables
        self.assmpt_dlf_errors = []
        self.frpp_df = None
        self.econ_df = pd.DataFrame()
        self.cfe_use_df = pd.DataFrame()
        self.agency_energy_data = None
        self.agency_price_data = None
        self.output_folder = DATA_PATH
        #self.threadpool = QThreadPool()
        self.thread = QThread()
        self.solar_power = None
        self.solar_energy = None
        self.solar_count = 0
        self.solar_max = 0

        # Set defaults
        self.actionExport_Results.setEnabled(False)
        self.stackedWidget.setCurrentIndex(0)
        self.tabWidget.setCurrentIndex(0)
        self.build_cfe_progressBar.setVisible(False)
        self.time_left_label_label.setVisible(False)

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
        
        if not os.path.exists(os.path.join(DATA_PATH,"BA_eGrid_Zip.csv")):
            msg = QMessageBox(self)
            msg.setIcon(QMessageBox.Icon.Warning)
            msg.setText("The BA_eGrid_Zip.csv seems to be missing from the data folder.\nThe calculation can continue but if you would like to stop and find the file, please hit cancel.")
            msg.setWindowTitle("Warning")
            msg.setStandardButtons(QMessageBox.StandardButton.Ok | QMessageBox.StandardButton.Cancel)
            ret = msg.exec()
            if ret == QMessageBox.StandardButton.Cancel:
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
        #pd.Int64Dtype()
        column_dtype = {'Reporting Agency': str, 'Real Property Unique Identifier': str, 'State Name': str, 'US/Foreign': str,
                        'Zip Code':str,'Latitude': np.double, 'Longitude': np.double, 'Real Property Type': str, 
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
        self.frpp_df['Zip Code'] = self.frpp_df['Zip Code'].apply(lambda x: x[0:5] if len(x) > 5 else int(x))
        self.frpp_df = self.frpp_df.astype({'Zip Code': int})
        self.pg1_progressBar.setValue(int(1/7*100))
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
        
        self.pg1_progressBar.setValue(int(2/7*100))
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

        self.pg1_progressBar.setValue(int(3/7*100))
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
            
            value = list(wind.sample([(xx, yy)]))[0][0]
            wind_values.append(value)
        assets_wind['Annual Average Wind Speed (m/s)'] = wind_values
        self.frpp_df = self.frpp_df.merge(assets_wind[['Real Property Unique Identifier', 'Annual Average Wind Speed (m/s)']], how = 'left', on = 'Real Property Unique Identifier')

        self.pg1_progressBar.setValue(int(4/7*100))
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

        # add balancing authority and egrid data
        if os.path.exists(os.path.join(DATA_PATH,"BA_eGrid_Zip.csv")):
            self.pg1_progressBar.setValue(int(5/7*100))
            self.statusbar.showMessage("Adding Egrid and Balancing Authority", 50000)
            self.landing_page.update()
            self.landing_page.repaint()
            column_dtype = {'Zip Code': pd.Int64Dtype(), 'eGRID Subregion': str,
                            'Balancing Authority ID': str, 'Balancing Authority': str}
            egrid = pd.read_csv(os.path.join(DATA_PATH,'BA_eGrid_Zip.csv'), dtype=column_dtype)
            # egrid.drop(['Agency Code', 'Country', 'State'], axis = 1, inplace = True)
            # egrid = egrid.loc[(egrid['Agency'] == self.agency_code())]
            # egrid.drop(['Agency'], axis = 1, inplace = True)
            self.frpp_df = self.frpp_df.merge(egrid, how = 'left', on = 'Zip Code')

        self.pg1_progressBar.setValue(int(6/7*100))
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

        self.pg1_progressBar.setValue(int(7/7*100))
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
        # Check to see if the agency is in the energy data csv and energy cost data csv
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
                self.agency_energy_data = agency_energy_data.loc[(agency_energy_data['Agency'] == str(self.agency_comboBox.currentText()))]

        if not os.path.exists(os.path.join(DATA_PATH, "Agency_Energy_Price_Data.csv")):
            msg = QMessageBox(self)
            msg.setIcon(QMessageBox.Icon.Warning)
            msg.setText(f"The agency energy price data csv is missing {os.path.join(DATA_PATH, 'Agency_Energy_Price_Data.csv')}\nThis data is needed to continue")
            msg.setWindowTitle("Warning")
            msg.setStandardButtons(QMessageBox.StandardButton.Ok)
            msg.exec()
            return
        else:
            agency_price_data = pd.read_csv(os.path.join(DATA_PATH, "Agency_Energy_Price_Data.csv"), header=0)
            if str(self.agency_comboBox.currentText()) not in agency_price_data['Agency'].values:
                msg = QMessageBox(self)
                msg.setIcon(QMessageBox.Icon.Warning)
                msg.setText(f"The data for the requested agency: {str(self.agency_comboBox.currentText())}\nIs not present in the price data set\n\nThe calculation can't continue")
                msg.setWindowTitle("Warning")
                msg.setStandardButtons(QMessageBox.StandardButton.Ok)
                msg.exec()
                return
            else:
                self.agency_price_data = agency_price_data.loc[(agency_price_data['Agency'] == str(self.agency_comboBox.currentText()))]

        self.toggle_cfe_running(True)
        
        build_model = buildCFEWorker(self.frpp_df.copy(), self.agency_energy_data.copy(), 
                                   self.agency_price_data.copy(), self.model_assump)
        build_model.moveToThread(self.thread)
        self.thread.started.connect(build_model.build_cfe_model)
        build_model.result_signal.connect(self.thread.quit)
        build_model.result_signal.connect(build_model.deleteLater)
        self.thread.finished.connect(self.thread.deleteLater)
        build_model.result_signal.connect(self.handle_cfeTask_results)
        build_model.progress_signal.connect(self.update_cfeTask_progress)
        self.thread.finished.connect(self.finalize_frpp_df) 
        self.thread.start()
        QApplication.processEvents()
        print("Thread starting")

    def finalize_frpp_df(self):
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
        self.frpp_df['Ground Solar Power Built'] = power[0,:]
        self.frpp_df['Concentrating Solar Power Built(kW)'] = power[1,:]
        self.frpp_df['Wind Power Built (kW)'] = power[2,:]
        self.frpp_df['Geothermal Power Built (kW)'] = power[3,:]
        self.frpp_df['Ground Solar Energy Built (kWh/yr)'] = energy[0,:]
        self.frpp_df['Concentrating Solar Energy Built (kWh)'] = energy[1,:]
        self.frpp_df['Wind Energy Built (kWh/yr)'] = energy[2,:]
        self.frpp_df['Geothermal Energy Built (kWh/yr)'] = energy[3,:]        
        # Finalize
        self.frpp_df['Total Power (kW)'] = self.frpp_df[['Rooftop Solar Power','Ground Solar Power Built',
                                                         'Wind Power Built (kW)','Fuel Cell (kW)','Geothermal Power Built (kW)',
                                                         'Concentrating Solar Power Built(kW)']].sum(axis=1)
           
        self.frpp_df['Total Energy (kWh)'] = self.frpp_df[['Annual Rooftop Solar Power (kWh/yr)','Ground Solar Energy Built (kWh/yr)',
                                                          'Wind Energy Built (kWh/yr)','Annual Fuel Cell (kW/yr)',
                                                          'Geothermal Energy Built (kWh/yr)','Concentrating Solar Energy Built (kWh)']].sum(axis=1)
        

        ''' Build CFE Energy Usage and Transition Table '''
        cur_year = datetime.now().year
        year_ind = list(range(cur_year,2036))
        
        pp = -(1000*self.agency_energy_data['Electricity (MWh)'].values[0] * (float(self.energy_proj_lineEdit.text())+100)) /\
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
            dmy.append(self.agency_energy_data['Electricity (MWh)'].values[0] * \
                (1. + float(self.energy_proj_lineEdit.text())/100.) * \
                (1 + float(self.AEG_lineEdit.text())/100.)**(year-2023))
        self.cfe_use_df['Projected Total Energy in year N (MWh)'] = dmy

        ''' Populate the next Page '''  
        self.populate_report()
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

        self.setup_econ_page()
        self.build_econ_model()
        self.toggle_cfe_running(False)
        self.actionExport_Results.setEnabled(True)
        self.stackedWidget.setCurrentIndex(2)

    def owner_switch(self):
        # Block the signals
        self.toggle_econ_signals(True)
        
        if self.gov_radioButton.isChecked():
            # Government Options
            cur_dict = self.model_assump[self.assump_from_econ_cfe()]['gov_rates']
            self.econ_rate_frame.setStyleSheet(u"font: 12pt \"Segoe UI\";\n"
                                    "background-color: rgb(200, 200, 200);")  
            self.econ_interest_rate_lineEdit.setText("2.0")            
            self.econ_cont_per_lineEdit.setText("2")
            self.econ_cont_per_Slider.setValue(2)            
            self.econ_proj_life_lineEdit.setText(str(cur_dict['project_life']))             
            self.econ_inflat_rate_lineEdit.setText(str(cur_dict['inflation_rate']))   
            self.econ_inflat_rate_Slider.setValue(max(1,int(cur_dict['inflation_rate']*10))) # increments of 0.1            
            self.econ_fed_tax_cred_lineEdit.setText("0.0")            
            self.econ_fed_tax_cred_Slider.setValue(0) 
            self.econ_fed_tax_cred_lineEdit.setEnabled(False)
            self.econ_fed_tax_cred_Slider.setEnabled(False)
            self.econ_tax_cred_per_lineEdit.setText("0")            
            self.econ_tax_cred_per_Slider.setValue(1) 
            self.econ_tax_cred_per_lineEdit.setEnabled(False)
            self.econ_tax_cred_per_Slider.setEnabled(False) 
            self.econ_tax_rate_lineEdit.setText(str(cur_dict['taxation']))  
            self.econ_tax_rate_Slider.setValue(int(cur_dict['taxation'])) 
            esc_val = 0.5 * (cur_dict['electric_cost_escalation_lb'] + cur_dict['electric_cost_escalation_ub']) 
            self.econ_price_esc_lineEdit.setText(str(esc_val))
            self.econ_price_esc_Slider.setValue(24)                   
            
        else:
            # Third Party Options
            cur_dict = self.model_assump[self.assump_from_econ_cfe()]['third_party_rates']
            self.econ_rate_frame.setStyleSheet(u"font: 12pt \"Segoe UI\";\n"
                                    "background-color: rgb(255, 146, 124);")            
            self.econ_interest_rate_lineEdit.setText("12.0")
            self.econ_cont_per_lineEdit.setText("25")
            self.econ_cont_per_Slider.setValue(25)
            self.econ_proj_life_lineEdit.setText(str(cur_dict['project_life']))
            self.econ_inflat_rate_lineEdit.setText(str(cur_dict['inflation_rate']))   
            self.econ_inflat_rate_Slider.setValue(max(1,int(cur_dict['inflation_rate']*10))) # increments of 0.1
            self.econ_fed_tax_cred_lineEdit.setText("30.0")            
            self.econ_fed_tax_cred_Slider.setValue(30) 
            self.econ_fed_tax_cred_lineEdit.setEnabled(True)
            self.econ_fed_tax_cred_Slider.setEnabled(True)            
            self.econ_tax_cred_per_lineEdit.setText("1")            
            self.econ_tax_cred_per_Slider.setValue(1) 
            self.econ_tax_cred_per_lineEdit.setEnabled(True)
            self.econ_tax_cred_per_Slider.setEnabled(True)
            self.econ_tax_rate_lineEdit.setText(str(cur_dict['taxation']))  
            self.econ_tax_rate_Slider.setValue(int(cur_dict['taxation']))
            esc_val = 0.5 * (cur_dict['electric_cost_escalation_lb'] + cur_dict['electric_cost_escalation_ub']) 
            self.econ_price_esc_lineEdit.setText(str(esc_val))
            self.econ_price_esc_Slider.setValue(24)       

        # Unblock the signals
        self.toggle_econ_signals(False)      
        self.econ_val_changed()

    def econ_val_changed(self):
        # Block Signals
        self.toggle_econ_signals(True)        

        # Begin by setting the derived quantaties   
        tot_eng = self.frpp_df[self.df_header_from_econ_cfe()].sum()*0.001  
        self.econ_power_produced_lineEdit.setText(f"{tot_eng:.2f}")
        
        self.econ_calc_engstor_lineEdit.setText(f"{(float(self.econ_energy_storage_lineEdit.text())/100.)*1000*tot_eng:.2f}")  

        cons = float(self.econ_agency_cons_lineEdit.text())
        tot_price = float(self.econ_agency_elec_cost_lineEdit.text())
        self.econ_unit_cost_lineEdit.setText(f"{tot_price/(cons*1000):.6f}") 

        if "Solar" in self.econ_cfe_comboBox.currentText() or "Hydrogen" in self.econ_cfe_comboBox.currentText():
            self.econ_ann_deg_lineEdit.setEnabled(True)
            try:
                float(self.econ_ann_deg_lineEdit.text())
            except ValueError:
                self.econ_ann_deg_lineEdit.setText("0.5")
        else:
            self.econ_ann_deg_lineEdit.setEnabled(False)

        valid, errors = self.validate_econ()
        if not valid:
            msg = QMessageBox(self)
            msg.setIcon(QMessageBox.Icon.Warning)
            msg.setText("\n\n".join(errors))
            msg.setWindowTitle("Warning")
            msg.setStandardButtons(QMessageBox.StandardButton.Ok)
            msg.exec()

        self.set_sliders()
        self.toggle_econ_signals(False)
        self.build_econ_model()

    def contract_period_slider(self):
        self.toggle_econ_signals(True)
        current_value = self.econ_cont_per_Slider.value()
        self.econ_cont_per_lineEdit.setText(f"{current_value}")
        self.toggle_econ_signals(False)

    def tax_rate_slider(self):
        self.toggle_econ_signals(True)
        current_value = self.econ_tax_rate_Slider.value()
        self.econ_tax_rate_lineEdit.setText(f"{current_value}")
        self.toggle_econ_signals(False)

    def insurance_rate_slider(self):
        self.toggle_econ_signals(True)
        current_value = float(self.econ_insur_rate_Slider.value())/10.
        self.econ_insur_rate_lineEdit.setText(f"{current_value}")
        self.toggle_econ_signals(False)

    def inflation_rate_slider(self):
        self.toggle_econ_signals(True)
        current_value = float(self.econ_inflat_rate_Slider.value())/10.
        self.econ_inflat_rate_lineEdit.setText(f"{current_value}")
        self.toggle_econ_signals(False)

    def fed_tax_credit_slider(self):
        self.toggle_econ_signals(True)
        current_value = self.econ_fed_tax_cred_Slider.value()
        self.econ_fed_tax_cred_lineEdit.setText(f"{current_value}")
        self.toggle_econ_signals(False)

    def tax_cred_period_slider(self):
        self.toggle_econ_signals(True)
        current_value = self.econ_tax_cred_per_Slider.value()
        self.econ_tax_cred_per_lineEdit.setText(f"{current_value}")
        self.toggle_econ_signals(False)

    def price_escalat_slider(self):
        self.toggle_econ_signals(True)
        if self.gov_radioButton.isChecked():
            cur_dict = self.model_assump[self.assump_from_econ_cfe()]['gov_rates']
        else:
            cur_dict = self.model_assump[self.assump_from_econ_cfe()]['third_party_rates']

        esc_val = cur_dict['electric_cost_escalation_lb'] + self.econ_price_esc_Slider.value()*\
            ((cur_dict['electric_cost_escalation_ub']-cur_dict['electric_cost_escalation_lb'])/50)
        self.econ_price_esc_lineEdit.setText(f"{esc_val}")  
        self.toggle_econ_signals(False) 

    def reset(self):
        '''
        Starts the analysis over from the beginning
        '''
        # Reset class variables
        self.assmpt_dlf_errors = []
        self.frpp_df = None
        self.econ_df = pd.DataFrame()
        self.cfe_use_df = pd.DataFrame()
        self.agency_energy_data = None
        self.agency_price_data = None
        self.cfe_bar_canvas.axes.cla() 
        self.econ_canvas.axes.cla() 

        # Restore default choices
        self.agency_comboBox.setCurrentIndex(0)
        self.pg1_progressBar.setValue(0)
        self.build_cfe_progressBar.setVisible(False)
        self.time_left_label_label.setVisible(False)

        # Move to the first page
        self.actionExport_Results.setEnabled(False)
        self.stackedWidget.setCurrentIndex(0)
        self.tabWidget.setCurrentIndex(0)
        self.statusbar.showMessage("Model Reset", 30000)

    def load(self):
        '''
        User will be asked to choose    
        '''        
        # Ask for a folder to save
        file, _ = QFileDialog.getOpenFileName(self,"Select a CFE Model output of the form Agency_Energy_Data.csv", str(self.output_folder), "CSV Files (*.csv);")
        if file:            
            error = self.agency_set(file.split("/")[-1].split("_")[0])
            if error == -1:
                msg = QMessageBox(self)
                msg.setIcon(QMessageBox.Icon.Warning)
                msg.setText("The provided csv did not have an agency code that matched our data.\n\nMake sure you only load data that has been saved from the file/export function of this program.")
                msg.setWindowTitle("Warning")
                msg.setStandardButtons(QMessageBox.StandardButton.Ok)
                msg.exec()
                self.statusbar.showMessage("Import Canceled", 30000)  
                return
            dummy = pd.read_csv(file, header=0)
            required_headers = ['Reporting Agency', 'Real Property Unique Identifier', 'US/Foreign', 'State Name', 'Zip Code', 
                                'Latitude', 'Longitude', 'Real Property Type', 'Real Property Use', 'Asset Status', 'Acres', 
                                'Square Feet (Buildings)', 'est_num_stories', 'est_rooftop_area_sqft', 'Geothermal_CLASS', 
                                'Annual Average Wind Speed (m/s)', 'Annual Solar Radiation kWh/m2/day', 'Rooftop Solar Power', 
                                'Ground Solar Power', 'Annual Rooftop Solar Power (kWh/yr)', 'Annual Ground Solar Power (kWh/yr)', 
                                'N Wind Turbines', 'Wind Power (kW)', 'Annual Wind Power (kWh/yr)', 'Concentrating Solar Power (kW)', 
                                'Annual Concentrating Solar Power (kWh)', 'Geothermal Power (kW)', 'Annual Geothermal Power (kWh/yr)', 
                                'Fuel Cell (kW)', 'Annual Fuel Cell (kW/yr)', 'Ground Solar Power Built', 'Concentrating Solar Power Built(kW)', 
                                'Wind Power Built (kW)', 'Geothermal Power Built (kW)', 'Ground Solar Energy Built (kWh/yr)', 
                                'Concentrating Solar Energy Built (kWh)', 'Wind Energy Built (kWh/yr)', 'Geothermal Energy Built (kWh/yr)', 
                                'Total Power (kW)', 'Total Energy (kWh)']
            if list(dummy.columns.values[1:]).sort() != required_headers.sort():
                msg = QMessageBox(self)
                msg.setIcon(QMessageBox.Icon.Warning)
                msg.setText("The provided csv did not match the required headers.\n\nMake sure you only load data that has been saved from the file/export function of this program.")
                msg.setWindowTitle("Warning")
                msg.setStandardButtons(QMessageBox.StandardButton.Ok)
                msg.exec()
                self.statusbar.showMessage("Import Canceled", 30000)   
                if self.frpp_df is None:
                    self.stackedWidget.setCurrentIndex(0)
                else:
                    self.stackedWidget.setCurrentIndex(1)
            else:
                self.frpp_df = dummy.copy()

                ''' Build CFE Energy Usage and Transition Table '''
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
                        self.agency_energy_data = agency_energy_data.loc[(agency_energy_data['Agency'] == str(self.agency_comboBox.currentText()))]

                if not os.path.exists(os.path.join(DATA_PATH, "Agency_Energy_Price_Data.csv")):
                    msg = QMessageBox(self)
                    msg.setIcon(QMessageBox.Icon.Warning)
                    msg.setText(f"The agency energy price data csv is missing {os.path.join(DATA_PATH, 'Agency_Energy_Price_Data.csv')}\nThis data is needed to continue")
                    msg.setWindowTitle("Warning")
                    msg.setStandardButtons(QMessageBox.StandardButton.Ok)
                    msg.exec()
                    return
                else:
                    agency_price_data = pd.read_csv(os.path.join(DATA_PATH, "Agency_Energy_Price_Data.csv"), header=0)
                    if str(self.agency_comboBox.currentText()) not in agency_price_data['Agency'].values:
                        msg = QMessageBox(self)
                        msg.setIcon(QMessageBox.Icon.Warning)
                        msg.setText(f"The data for the requested agency: {str(self.agency_comboBox.currentText())}\nIs not present in the price data set\n\nThe calculation can't continue")
                        msg.setWindowTitle("Warning")
                        msg.setStandardButtons(QMessageBox.StandardButton.Ok)
                        msg.exec()
                        return
                    else:
                        self.agency_price_data = agency_price_data.loc[(agency_price_data['Agency'] == str(self.agency_comboBox.currentText()))]
                
                cur_year = datetime.now().year
                year_ind = list(range(cur_year,2036))
                
                pp = -(1000*self.agency_energy_data['Electricity (MWh)'].values[0] * (float(self.energy_proj_lineEdit.text())+100)) /\
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
                    dmy.append(self.agency_energy_data['Electricity (MWh)'].values[0] * \
                        (1. + float(self.energy_proj_lineEdit.text())/100.) * \
                        (1 + float(self.AEG_lineEdit.text())/100.)**(year-2023))
                self.cfe_use_df['Projected Total Energy in year N (MWh)'] = dmy

                ''' Populate the next Page '''        
                self.populate_report()
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

                self.setup_econ_page()
                self.build_econ_model()

                self.actionExport_Results.setEnabled(True)
                self.stackedWidget.setCurrentIndex(2)
        else:
            self.statusbar.showMessage("Import Canceled", 30000)

    def export(self):
        '''
        User will be asked to provide an export folder then 
        the following tables will be exported as csv.
        - self.frpp_df
        - self.econ_df
        - self.cfe_use_df    
        '''

        # Ask for a folder to save
        folder = QFileDialog.getExistingDirectory(self,"Select a directory to export CFE Results", str(self.output_folder))
        if folder:
            self.output_folder = folder
            now = datetime.now()
            new_folder = os.path.join(folder, f"CFE_Export_{now.year}{now.month}{now.day}{now.hour}{now.minute}")
            if not os.path.exists(new_folder):
                os.mkdir(new_folder)

            self.frpp_df.to_csv(os.path.join(new_folder, f"{self.agency_code()}_Energy_Data.csv"))
            self.econ_df.to_csv(os.path.join(new_folder, f"{self.agency_code()}_Econ_Data.csv"))
            self.cfe_use_df.to_csv(os.path.join(new_folder, f"{self.agency_code()}_Energy_Transition.csv"))
        else:
            self.statusbar.showMessage("Export Canceled", 30000)

    # Functions not tied to any button
    def build_econ_model(self):
        self.econ_df = pd.DataFrame()
        eng_key = self.df_header_from_econ_cfe(req="Energy")
        # Average install cost in ($/kW)
        # Fixed O&M cost in ($/kW/year)
        cost = {"Rooftop Solar": {"install_cost": 1592., "fixed_om": 18.},
            "Ground-Mounted Solar": {"install_cost": 1327., "fixed_om": 18.},
            "Wind": {"install_cost": 1718., "fixed_om": 42.},
            "Geothermal": {"install_cost": 3076., "fixed_om": 143.22},
            "Hydrogen Fuel Cells": {"install_cost": 6639., "fixed_om": 30.65},
            "Concentrating Solar": {"install_cost": 2383., "fixed_om": 67.32}}
        cur_year = datetime.now().year
        year_ind = list(range(cur_year + 1,cur_year + 1 + int(float(self.econ_proj_life_lineEdit.text()))))
        self.econ_df['Year'] = year_ind

        strg_cost = 0
        try:
            eng_strg_pec = float(self.econ_energy_storage_lineEdit.text())
            strg_cost = float(self.econ_calc_engstor_lineEdit.text())*1383. + \
                ((0.*self.frpp_df[eng_key].sum()/1000.)*(eng_strg_pec/100))
            # The ^ 0 is listed as the fuel cost ($/MWh)
        except: 
            pass

        investment = cost[self.econ_cfe_comboBox.currentText()]['install_cost']* \
            float(self.econ_power_produced_lineEdit.text())*1000 + strg_cost        
        
        # Annual Loan Payments
        dmy = np.zeros(len(year_ind))
        ind = int(float(self.econ_cont_per_lineEdit.text()))
        dmy[:ind] = -1*(investment * (float(self.econ_interest_rate_lineEdit.text())/100.)) / \
            (1-(1+(float(self.econ_interest_rate_lineEdit.text())/100.))**\
             (-float(self.econ_cont_per_lineEdit.text())))
        self.econ_df['Anual Loan Payments'] = dmy.tolist() 
        
        # MACRS depreciation
        # https://xplaind.com/370120/macrs#google_vignette
        dmy = np.zeros(len(year_ind))
        if self.third_party_radioButton.isChecked():
            macrs_methods = {0: [3, 2.0], 1: [5, 2.0], 2: [7, 2.0],
                            3: [10,2.0], 4: [15,1.5], 5: [20,1.5]}
            macrs = macrs_methods[self.econ_marcs_comboBox.currentIndex()]
            if self.econ_cfe_comboBox.currentText() == "Wind":
                dereciable = investment * 0.6
            else:
                dereciable = ((100-float(self.econ_fed_tax_cred_lineEdit.text())/2.)/100.) * investment * 0.6
            dmy[0] = dereciable * (1/macrs[0]) * macrs[1] * 0.5
            for qq in range(1,min(macrs[0]+1, len(year_ind))):
                macr = (dereciable - np.sum(dmy[:qq])) * (1/macrs[0]) * macrs[1]
                remaining_life = (macrs[0]-qq)+0.5                
                strt_line = (dereciable - np.sum(dmy[:qq]))*(1/remaining_life)
                if remaining_life == 0.5:
                    strt_line = dmy[qq-1]*remaining_life
                dmy[qq] = max(macr,strt_line)
        
        self.econ_df['MACRS Depreciation'] = dmy.tolist() 

        # Energy by year
        dmy = np.zeros(len(year_ind))
        dmy[0] = self.frpp_df[eng_key].sum()
        if self.econ_cfe_comboBox.currentText() in ["Rooftop Solar", "Ground-Mounted Solar","Hydrogen Fuel Cells"]:
            for qq in range(1,len(dmy)):
                dmy[qq] = dmy[qq-1] *(1-float(self.econ_ann_deg_lineEdit.text())/100.)
        else:
            dmy[:] = self.frpp_df[eng_key].sum()

        self.econ_df['Energy (kWh)'] = dmy.tolist() 

        # Discounted Energy
        self.econ_df['Discounted Energy'] = self.econ_df['Energy (kWh)'] / \
            (1+float(self.econ_interest_rate_lineEdit.text())/100.)**(self.econ_df['Year']-cur_year)

        # Price 
        dmy = np.zeros(len(year_ind))
        dmy[0] = float(self.econ_unit_cost_lineEdit.text())
        for qq in range(1,len(dmy)):
            dmy[qq] = dmy[0]*(1+float(self.econ_price_esc_lineEdit.text())/100.)**(year_ind[qq]-cur_year)
        self.econ_df['Price ($/kWh)'] = dmy.tolist() 

        # Revenue
        dmy = np.zeros(len(year_ind))
        if self.third_party_radioButton.isChecked():
            for qq in range(len(dmy)):
                if qq < int(float(self.econ_cont_per_lineEdit.text())):
                    dmy[qq] = self.econ_df.iloc[qq]['Energy (kWh)'] * self.econ_df.iloc[qq]['Price ($/kWh)']
                else:
                    dmy[qq] = self.econ_df.iloc[qq]['Energy (kWh)'] * float(self.econ_unit_cost_lineEdit.text())
        self.econ_df['Revenue ($)'] = dmy.tolist()

        # Cost without CFE Contract
        self.econ_df['Cost without CFE ($)'] = (1 + float(self.econ_price_esc_lineEdit.text())/100)**(self.econ_df['Year']-cur_year)*\
            float(self.econ_unit_cost_lineEdit.text())*self.econ_df['Energy (kWh)'] 
        
        # Government Annual Savings ($)
        self.econ_df['Government Annual Savings ($)'] = self.econ_df['Cost without CFE ($)'] - self.econ_df['Revenue ($)']

        # Battery Replacement 
        n_years = int(float(self.econ_bat_rep_lineEdit.text()))
        dmy = np.zeros(len(year_ind))
        for qq in range(1,len(dmy)):
            if qq%n_years == 0:
                dmy[qq] = strg_cost
        self.econ_df['Battery Replacement Cost ($)'] = dmy.tolist()

        # Insurance 
        self.econ_df['Insurance ($)'] = investment*(float(self.econ_insur_rate_lineEdit.text())/100.)*\
            (1 + float(self.econ_inflat_rate_lineEdit.text())/100)**(self.econ_df['Year']-(cur_year+1))
        
        # O&M Cost
        om_cost = cost[self.econ_cfe_comboBox.currentText()]['fixed_om']*\
            float(self.econ_power_produced_lineEdit.text()) * 1000.
        dmy = np.zeros(len(year_ind))
        dmy[0] = om_cost
        for qq in range(1,len(dmy)):
            dmy[qq] = om_cost*(1+float(self.econ_inflat_rate_lineEdit.text())/100.)**(qq)
        self.econ_df['OM Cost ($)'] = dmy.tolist()

        # Government Expenses 
        dmy = np.zeros(len(year_ind))
        self.econ_df['Government Expenses ($)'] = dmy.tolist()
        if self.gov_radioButton.isChecked():
            self.econ_df['Government Expenses ($)'] = 0 - self.econ_df['Battery Replacement Cost ($)'] -\
                 self.econ_df['Insurance ($)'] - self.econ_df['OM Cost ($)']

        # After-tax Profit + Depreciation ($)
        if self.third_party_radioButton.isChecked():
            self.econ_df['After-tax Profit + Depreciation ($)'] = (self.econ_df['Revenue ($)'] - self.econ_df['OM Cost ($)'] - \
                           self.econ_df['MACRS Depreciation'] - self.econ_df['Insurance ($)'] - self.econ_df['Battery Replacement Cost ($)']) *\
                           (1-float(self.econ_tax_rate_lineEdit.text())*0.01) + self.econ_df['MACRS Depreciation']
        else:
            self.econ_df['After-tax Profit + Depreciation ($)'] = -1*(self.econ_df['Battery Replacement Cost ($)'] +\
                                                                       self.econ_df['Insurance ($)'] + self.econ_df['OM Cost ($)'])

        # Net Taxes ($)
        if self.third_party_radioButton.isChecked():
            dmy = np.zeros(len(year_ind))
            self.econ_df['Net Taxes ($)'] = dmy.tolist()
        else:
            self.econ_df['Net Taxes ($)'] = self.econ_df['Revenue ($)'] - self.econ_df['After-tax Profit + Depreciation ($)']

        # Cash Flow
        tax_cred = np.zeros(len(year_ind))
        prod_ben = np.zeros(len(year_ind))
        if self.third_party_radioButton.isChecked():
            for qq in range(int(float(self.econ_tax_cred_per_lineEdit.text()))):
                tax_cred[qq] = investment * float(self.econ_fed_tax_cred_lineEdit.text())/100.
            for qq in range(10):
                prod_ben[qq] = 0.015 * (1+float(self.econ_price_esc_lineEdit.text())/100.)**(qq+1) * self.econ_df.iloc[qq]['Energy (kWh)']
        rec = np.zeros(len(year_ind))
        for qq in range(len(rec)):
            rec[qq] = float(self.econ_rec_lineEdit.text()) * self.econ_df.iloc[qq]['Energy (kWh)'] / 1000.
        #### PTC tax credit (only counts for the first 10 years)
        ptc = np.zeros(len(year_ind))
        ptc[0] = investment * 0.4

        dmy = np.zeros(len(year_ind))
        for qq in range(len(dmy)):
            if self.gov_radioButton.isChecked():
                dmy[qq] = self.econ_df.iloc[qq]['After-tax Profit + Depreciation ($)'] + self.econ_df.iloc[qq]['MACRS Depreciation']\
                      + self.econ_df.iloc[qq]['Cost without CFE ($)'] 
            else:
                dmy[qq] = self.econ_df.iloc[qq]['After-tax Profit + Depreciation ($)']
            if self.econ_cfe_comboBox.currentText() in ["Rooftop Solar", "Ground-Mounted Solar", "Concentrating Solar", "Geothermal"]:
                dmy[qq] += tax_cred[qq] + rec[qq]
            elif self.econ_cfe_comboBox.currentText() in ["Wind"]:
                dmy[qq] += ptc[qq] + prod_ben[qq] + rec[qq]
        self.econ_df['Cash Flow ($)'] = dmy.tolist()

        # Discounted Annualized
        self.econ_df['Discounted Annualized ($)'] = (1/(1+float(self.econ_interest_rate_lineEdit.text())/ \
                                                      100)**(self.econ_df['Year']-cur_year))*self.econ_df['Cash Flow ($)']
        
        # Cumulative Cash Flow
        dmy = np.zeros(len(year_ind))
        dmy[0] = -1*investment + self.econ_df.iloc[0]['Discounted Annualized ($)']
        for qq in range(1,len(dmy)):
            dmy[qq] = dmy[qq-1] + self.econ_df.iloc[qq]['Discounted Annualized ($)']
        self.econ_df['Cumulative Cash Flow ($)'] = dmy.tolist()

        # Build the bar chart    
        self.econ_canvas.axes.cla()    
        x = np.arange(len(year_ind)+1)  # the label locations
        self.econ_canvas.axes.bar(x, np.concatenate(([-investment],
                                  self.econ_df['Cumulative Cash Flow ($)'].values[:])), 0.8, label="NPV") 
        self.econ_canvas.axes.grid(True, axis="y", color = "grey", linewidth = "1.4")          

        # Add some text for labels, title and custom x-axis tick labels, etc.
        self.econ_canvas.axes.set_ylabel('Net Present Value ($)')
        self.econ_canvas.axes.set_title('Cummulative Cash Flow by Fiscal Year')
        self.econ_canvas.axes.set_xticks(x, [cur_year] + year_ind)
        xlabels = ['%i' % i for i in [cur_year]+ year_ind]
        self.econ_canvas.axes.set_xticklabels(xlabels, rotation=45)
        self.econ_canvas.draw()

        # Fill in the output textedits
        ltZero = self.econ_df['Anual Loan Payments'] < 0
        eqZero = self.econ_df['Anual Loan Payments'] == 0

        self.econ_OM_lineEdit.setText(f"{om_cost:,.2f}")
        self.econ_net_cap_cost_lineEdit.setText(f"{investment:,.2f}")
        self.econ_irr_lineEdit.setText(f"{npf.irr(np.concatenate(([-investment],self.econ_df['Cash Flow ($)'].values[:])))*100:.2f}")
        self.econ_dev_npv_lineEdit.setText(f"{self.econ_df[ltZero]['Discounted Annualized ($)'].sum() - investment:,.2f}")
        self.econ_total_energy_lineEdit.setText(f"{self.econ_df['Energy (kWh)'].sum():,.2f}")
        dev_lcc = abs(self.econ_df['Anual Loan Payments'].sum()) + self.econ_df[ltZero]['Insurance ($)'].sum() +\
              self.econ_df[ltZero]['OM Cost ($)'].sum()
        self.econ_developer_LCC_lineEdit.setText(f"{dev_lcc:,.2f}")
        if self.third_party_radioButton.isChecked():
            total = self.econ_df[ltZero]['Revenue ($)'].sum() + self.econ_df[eqZero]['Battery Replacement Cost ($)'].sum() + \
                self.econ_df[eqZero]['Insurance ($)'].sum() + self.econ_df[eqZero]['OM Cost ($)'].sum()            
        else: 
            total = investment + self.econ_df['Battery Replacement Cost ($)'].sum() + \
                self.econ_df['Insurance ($)'].sum() + self.econ_df['OM Cost ($)'].sum()
        self.econ_owner_LCC_lineEdit.setText(f"{total:,.2f}")
        self.econ_annual_payments_lineEdit.setText(f"{total/float(self.econ_cont_per_lineEdit.text()):,.2f}")
        self.econ_owner_savings_lineEdit.setText(f"{self.econ_df['Cost without CFE ($)'].sum()-total:,.2f}")
        if self.third_party_radioButton.isChecked():
            lvl_cst = (npf.npv(float(self.econ_interest_rate_lineEdit.text())/100, \
                               np.concatenate(([0.],self.econ_df['Revenue ($)'].values[:]))))*\
                                1000/self.econ_df['Discounted Energy'].sum()
        else:
            lvl_cst = (-npf.npv(float(self.econ_interest_rate_lineEdit.text())/100, \
                np.concatenate(([-investment],self.econ_df['Government Expenses ($)'].values[:]))))*\
                1000/self.econ_df['Discounted Energy'].sum()
        self.lvl_cost_lineEdit.setText(f"{lvl_cst:.2f}")
        self.cost_benefit_lineEdit.setText(f"{self.econ_df['Cost without CFE ($)'].sum()/total:.2f}")
        smpl_pybk = total / (self.econ_df['Cost without CFE ($)'].sum()/ \
                             float(self.econ_proj_life_lineEdit.text()))
        self.simp_payback_lineEdit.setText(f"{smpl_pybk:.2f}")

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

    def validate_econ(self)->bool | list:
        """
        Goes through all of the econ tab and makes sure that all
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

        errors = ["There was an issue with the requested input:"]
        valid = True

        if self.gov_radioButton.isChecked():
            cur_dict = self.model_assump[self.assump_from_econ_cfe()]['gov_rates']
            int_rate = 2
        else:
            cur_dict = self.model_assump[self.assump_from_econ_cfe()]['third_party_rates']
            int_rate = 12

        # Check the Energy Storage Percentage 
        try:
            val = float(self.econ_energy_storage_lineEdit.text())   
            if val < 0:
                valid = False
                errors.append(f"The requested percent energy storage {val} is outside the acceptable bounds of 0 - 100%.\nThe value has been changed to 0")
                self.econ_energy_storage_lineEdit.setText("0")
            elif val > 100:
                valid = False
                errors.append(f"The requested percent energy storage {val} is outside the acceptable bounds of 0 - 100%.\nThe value has been changed to 100")
                self.econ_energy_storage_lineEdit.setText("100")
        except ValueError:            
            valid = False
            errors.append(f"{self.econ_energy_storage_lineEdit.text()} was not a valid choice for energy storage.\nThe value has been changed to 0")
            self.econ_energy_storage_lineEdit.setText("0")

        # Check the Battery Replacement
        try:
            val = float(self.econ_bat_rep_lineEdit.text())   
            if val < 1:
                valid = False
                errors.append(f"The requested battery replacement time: {val}, is outside the acceptable bounds.\nThe value has been changed to 1")
                self.econ_bat_rep_lineEdit.setText("1")
        except ValueError:            
            valid = False
            errors.append(f"{self.econ_bat_rep_lineEdit.text()} was not a valid choice for battery replacement.\nThe value has been changed to 5")
            self.econ_bat_rep_lineEdit.setText("5")

        # Check the Agency Consumption
        try:
            val = float(self.econ_agency_cons_lineEdit.text())   
            if val < 1:
                valid = False
                errors.append(f"The requested annual agency energy consumption: {val}, must be greater than 0.\nThe value has been changed to the default 2022 value")
                self.econ_agency_cons_lineEdit.setText(f"{self.agency_energy_data['Electricity (MWh)'].values[0]}")
        except ValueError:            
            valid = False
            errors.append(f"{self.econ_agency_cons_lineEdit.text()} was not a valid choice for agency consumption.\nThe value has been changed to the default 2022 value")
            self.econ_agency_cons_lineEdit.setText(f"{self.agency_energy_data['Electricity (MWh)'].values[0]}")

        # Check the Agency Electricity Costs
        try:
            val = float(self.econ_agency_elec_cost_lineEdit.text())   
            if val < 1:
                valid = False
                errors.append(f"The requested annual agency energy consumption: {val}, must be greater than 0.\nThe value has been changed to the default 2022 value")
                self.econ_agency_elec_cost_lineEdit.setText(f"{self.agency_price_data['Electricity'].values[0]}")
        except ValueError:            
            valid = False
            errors.append(f"{self.econ_agency_elec_cost_lineEdit.text()} was not a valid choice for agency consumption.\nThe value has been changed to the default 2022 value")
            self.econ_agency_elec_cost_lineEdit.setText(f"{self.agency_price_data['Electricity'].values[0]}")

        # Check the Interest Rate
        try:
            val = float(self.econ_interest_rate_lineEdit.text())               
        except ValueError:            
            valid = False
            errors.append(f"{self.econ_interest_rate_lineEdit.text()} was not a valid choice for interest rate.\nThe value has been changed to {int_rate}")
            self.econ_interest_rate_lineEdit.setText(f"{int_rate}")

        # Check the project life
        try:
            val = float(self.econ_proj_life_lineEdit.text())   
            if val < 1:
                valid = False
                errors.append(f"The requested project life: {val}, is outside the acceptable bounds.\nThe value has been changed to {str(cur_dict['project_life'])}")
                self.econ_proj_life_lineEdit.setText({str(cur_dict['project_life'])})
        except ValueError:            
            valid = False
            errors.append(f"{self.econ_proj_life_lineEdit.text()} was not a valid choice for project life.\nThe value has been changed to {str(cur_dict['project_life'])}")
            self.econ_proj_life_lineEdit.setText(str(cur_dict['project_life'])) 

        # Check annual degradation
        if "Solar" in self.econ_cfe_comboBox.currentText() or "Hydrogen" in self.econ_cfe_comboBox.currentText():
            try:
                val = float(self.econ_ann_deg_lineEdit.text())   
                if val < 0:
                    valid = False
                    errors.append(f"The requested annual degradation: {val}, is outside the acceptable bounds.\nThe value has been changed to 0.5")
                    self.econ_ann_deg_lineEdit.setText("0.5")
            except ValueError:            
                valid = False
                errors.append(f"{self.econ_ann_deg_lineEdit.text()} was not a valid choice for annual degradation.\nThe value has been changed to 0.5")
                self.econ_ann_deg_lineEdit.setText("0.5") 

        # Check renewable energy certification
        try:
            val = float(self.econ_rec_lineEdit.text())   
            if val < 0:
                valid = False
                errors.append(f"The requested renewable energy certification: {val}, is outside the acceptable bounds.\nThe value has been changed to 0")
                self.econ_rec_lineEdit.setText("0")
        except ValueError:            
            valid = False
            errors.append(f"{self.econ_rec_lineEdit.text()} was not a valid choice for renewable energy certification.\nThe value has been changed to 19")
            self.econ_rec_lineEdit.setText(f"{19}") 

        # Check contract period
        current_value = self.econ_cont_per_Slider.value()
        try:
            val = float(self.econ_cont_per_lineEdit.text())   
            if val < 1:
                valid = False
                errors.append(f"The requested contract period: {val}, is outside the acceptable bounds of 1 - 30.\nThe value has been changed to {current_value}")
                self.econ_cont_per_lineEdit.setText(f"{current_value}")
            elif val > 30:
                valid = False
                errors.append(f"The requested contract period: {val}, is outside the acceptable bounds of 1 - 30.\nThe value has been changed to {current_value}")
                self.econ_cont_per_lineEdit.setText(f"{current_value}")
        except ValueError:            
            valid = False
            errors.append(f"{self.econ_cont_per_lineEdit.text()} was not a valid choice for contract period.\nThe value has been changed to {current_value}")
            self.econ_cont_per_lineEdit.setText(f"{current_value}") 

        # Check tax rate
        current_value = self.econ_tax_rate_Slider.value()
        try:
            val = float(self.econ_tax_rate_lineEdit.text())   
            if val < 0:
                valid = False
                errors.append(f"The requested tax rate: {val}, is outside the acceptable bounds of 0 - 50%.\nThe value has been changed to {current_value}")
                self.econ_tax_rate_lineEdit.setText(f"{current_value}")
            elif val > 50:
                valid = False
                errors.append(f"The requested tax rate: {val}, is outside the acceptable bounds of 0 - 50%.\nThe value has been changed to {current_value}")
                self.econ_tax_rate_lineEdit.setText(f"{current_value}")
        except ValueError:            
            valid = False
            errors.append(f"{self.econ_tax_rate_lineEdit.text()} was not a valid choice for tax rate.\nThe value has been changed to {current_value}")
            self.econ_tax_rate_lineEdit.setText(f"{current_value}") 

        # Check insurance rate
        current_value = self.econ_insur_rate_Slider.value()/10.
        try:
            val = float(self.econ_insur_rate_lineEdit.text())   
            if val < 0.1:
                valid = False
                errors.append(f"The requested insurance rate: {val}, is outside the acceptable bounds of 0.1 - 5.0%.\nThe value has been changed to {current_value}")
                self.econ_insur_rate_lineEdit.setText(f"{current_value}")
            elif val > 5:
                valid = False
                errors.append(f"The requested insurance rate: {val}, is outside the acceptable bounds of 0.1 - 5.0%.\nThe value has been changed to {current_value}")
                self.econ_insur_rate_lineEdit.setText(f"{current_value}")
        except ValueError:            
            valid = False
            errors.append(f"{self.econ_insur_rate_lineEdit.text()} was not a valid choice for insurance rate.\nThe value has been changed to {current_value}")
            self.econ_insur_rate_lineEdit.setText(f"{current_value}") 

        if self.third_party_radioButton.isChecked():
            # Check inflation rate
            current_value = self.econ_inflat_rate_Slider.value()/10
            try:
                val = float(self.econ_inflat_rate_lineEdit.text())   
                if val < 0.1:
                    valid = False
                    errors.append(f"The requested inflation rate: {val}, is outside the acceptable bounds of 0.1 - 10.0%.\nThe value has been changed to {current_value}")
                    self.econ_inflat_rate_lineEdit.setText(f"{current_value}")
                elif val > 10.0:
                    valid = False
                    errors.append(f"The requested inflation rate: {val}, is outside the acceptable bounds of 0.1 - 10.0%.\nThe value has been changed to {current_value}")
                    self.econ_inflat_rate_lineEdit.setText(f"{current_value}")
            except ValueError:            
                valid = False
                errors.append(f"{self.econ_inflat_rate_lineEdit.text()} was not a valid choice for inflation rate.\nThe value has been changed to {current_value}")
                self.econ_inflat_rate_lineEdit.setText(f"{current_value}") 
            
            # Check federal tax credit        
            current_value = self.econ_fed_tax_cred_Slider.value()
            try:
                val = float(self.econ_fed_tax_cred_lineEdit.text())   
                if val < 0:
                    valid = False
                    errors.append(f"The requested federal tax credit: {val}, is outside the acceptable bounds of 0 - 50%.\nThe value has been changed to {current_value}")
                    self.econ_fed_tax_cred_lineEdit.setText(f"{current_value}")
                elif val > 50.0:
                    valid = False
                    errors.append(f"The requested federal tax credit: {val}, is outside the acceptable bounds of 0 - 50%.\nThe value has been changed to {current_value}")
                    self.econ_fed_tax_cred_lineEdit.setText(f"{current_value}")
            except ValueError:            
                valid = False
                errors.append(f"{self.econ_fed_tax_cred_lineEdit.text()} was not a valid choice for federal tax credit.\nThe value has been changed to {current_value}")
                self.econ_fed_tax_cred_lineEdit.setText(f"{current_value}") 

            # Check federal tax credit period        
            current_value = self.econ_tax_cred_per_Slider.value()
            try:
                val = float(self.econ_tax_cred_per_lineEdit.text())   
                if val < 1:
                    valid = False
                    errors.append(f"The requested federal tax credit period: {val}, is outside the acceptable bounds of 1 - 10.\nThe value has been changed to {current_value}")
                    self.econ_tax_cred_per_lineEdit.setText(f"{current_value}")
                elif val > 10:
                    valid = False
                    errors.append(f"The requested federal tax credit period: {val}, is outside the acceptable bounds of 1 - 10.\nThe value has been changed to {current_value}")
                    self.econ_tax_cred_per_lineEdit.setText(f"{current_value}")
            except ValueError:            
                valid = False
                errors.append(f"{self.econ_tax_cred_per_lineEdit.text()} was not a valid choice for federal tax creditperiod.\nThe value has been changed to {current_value}")
                self.econ_tax_cred_per_lineEdit.setText(f"{current_value}") 
            
        # Check price escalation
        current_value = self.econ_price_esc_Slider.value()
        lb = cur_dict['electric_cost_escalation_lb']
        ub = cur_dict['electric_cost_escalation_ub']
        esc_val = lb + current_value*((ub-lb)/50)
        try:
            val = float(self.econ_price_esc_lineEdit.text())   
            if val < lb:
                valid = False
                errors.append(f"The requested price escalation: {val}, is outside the acceptable bounds of {lb} - {ub}%/year.\nThe value has been changed to {esc_val}")
                self.econ_price_esc_lineEdit.setText(f"{esc_val}")
            elif val > ub:
                valid = False
                errors.append(f"The requested price escalation: {val}, is outside the acceptable bounds of {lb} - {ub}%/year.\nThe value has been changed to {esc_val}")
                self.econ_price_esc_lineEdit.setText(f"{esc_val}")
        except ValueError:            
            valid = False
            errors.append(f"{self.econ_price_esc_lineEdit.text()} was not a valid choice for price escalation.\nThe value has been changed to {esc_val}")
            self.econ_price_esc_lineEdit.setText(f"{esc_val}") 

        return valid, errors

    def update_cfeTask_progress(self, val, run_label):           
        self.build_cfe_progressBar.setValue(val)        
        self.time_left_label_label.setText(f"{run_label} is being calculated")

    def handle_cfeTask_results(self, frpp_df):        
        self.frpp_df = frpp_df.copy()

    def toggle_cfe_running(self, running:bool):
        if running:
            self.energy_proj_lineEdit.setReadOnly(True)
            self.cur_cfe_lineEdit.setReadOnly(True)
            self.AEG_lineEdit.setReadOnly(True)
            self.CB_oper_lineEdit.setReadOnly(True)
            self.ren_oper_lineEdit.setReadOnly(True)
            self.oper_days_lineEdit.setReadOnly(True)
            self.model_assumpt_pushButton.setEnabled(False)
            self.build_CFE_pushButton.setEnabled(False)
            self.build_cfe_progressBar.setVisible(True)
            self.time_left_label_label.setVisible(True)
        else:
            self.energy_proj_lineEdit.setReadOnly(False)
            self.cur_cfe_lineEdit.setReadOnly(False)
            self.AEG_lineEdit.setReadOnly(False)
            self.CB_oper_lineEdit.setReadOnly(False)
            self.ren_oper_lineEdit.setReadOnly(False)
            self.oper_days_lineEdit.setReadOnly(False)
            self.model_assumpt_pushButton.setEnabled(True)
            self.build_CFE_pushButton.setEnabled(True)
            self.build_cfe_progressBar.setVisible(False)
            self.time_left_label_label.setVisible(False)
 
    def populate_report(self):
        n_solar_roof = np.sum(~self.frpp_df['est_rooftop_area_sqft'].isna())
        n_solar_grnd = np.sum(~self.frpp_df['Acres'].isna())
        n_wind = np.sum(~self.frpp_df["Wind Power (kW)"].isna())
        n_geot = np.sum(np.logical_and(~self.frpp_df['Geothermal_CLASS'].isna(), self.frpp_df['Geothermal_CLASS'] <=3))
        n_fuelcell = self.frpp_df.shape[0]
        n_conc_sol = np.sum(np.logical_and(self.frpp_df['Real Property Type'] == 'Land', self.frpp_df['Real Property Use'] == 'Vacant'))

        n_solar_roof_built = np.sum(~self.frpp_df['Rooftop Solar Power'].isna().values)
        n_solar_grnd_built = np.sum(np.logical_and(~self.frpp_df['Ground Solar Power Built'].isna().values,
                                                   self.frpp_df['Ground Solar Power Built'] != 0))
        n_wind_built = np.sum(np.logical_and(~self.frpp_df['Wind Power Built (kW)'].isna().values,
                                             self.frpp_df['Wind Power Built (kW)'] != 0))
        n_geot_built = np.sum(np.logical_and(~self.frpp_df['Geothermal Power Built (kW)'].isna().values,
                                             self.frpp_df['Geothermal Power Built (kW)'] != 0))
        n_fuelcell_built = np.sum(~self.frpp_df['Fuel Cell (kW)'].isna().values)
        n_conc_sol_built = np.sum(np.logical_and(~self.frpp_df['Concentrating Solar Power Built(kW)'].isna().values,
                                                 self.frpp_df['Concentrating Solar Power Built(kW)'] != 0))

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
                                    f":0; text-indent:0px;\"><span style=\" font-size:14pt;\">Total Energy Required: 		{self.agency_energy_data['Electricity (MWh)'].values[0]:,.0f} MWh</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt;\">Total Estimated Production Capacity:	{self.frpp_df['Total Energy (kWh)'].sum()*0.001:,.0f} MWh</span></p>\n"
            "<p style=\"-qt-paragraph-type:empty;  -qt-block-indent:0; text-indent:0px; font-size:14pt;\"><br /></p>\n"
            "<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:22pt; color:#990000;\">Wind</span></p>\n"
            "<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:10pt; color:#990000;\">===================================================</span></p>\n"
            "<p style=\"  -qt-block-i"
                                    f"ndent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Number of Sites available: {n_wind}</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Number of sites built: {n_wind_built}</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Total Potential: {self.frpp_df['Wind Power (kW)'].sum()*0.001:,.0f} MW</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Total Annual Generation: {self.frpp_df['Annual Wind Power (kWh/yr)'].sum()*0.001:,.0f} MWh</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Percent of Agency Demand: {(self.frpp_df['Annual Wind Power (kWh/yr)'].sum()*0.001)/self.agency_energy_data['Electricity (MWh)'].values[0]*100:.1f}%</span></p>\n"
            "<p style=\"-qt-paragraph-type:empty;  -qt-block-indent:0; text-indent:0px;  font-size:14pt;\"><br /></p>\n"
            "<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:22pt; color:#990000;\">Rooftop Solar PV</span></p>\n"
            "<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:10pt; color:#990000;\">===================================================</span></p>\n"
            "<p style=\"  -qt-block-i"
                                    f"ndent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Number of Sites available: {n_solar_roof}</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Number of sites built: {n_solar_roof_built}</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Total Potential: {self.frpp_df['Rooftop Solar Power'].sum()*0.001:,.0f} MW</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Total Annual Generation: {self.frpp_df['Annual Rooftop Solar Power (kWh/yr)'].sum()*0.001:,.0f} MWh</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Percent of Agency Demand: {(self.frpp_df['Annual Rooftop Solar Power (kWh/yr)'].sum()*0.001)/self.agency_energy_data['Electricity (MWh)'].values[0]*100:.1f}%</span></p>\n"
            "<p style=\"-qt-paragraph-type:empty;  -qt-block-indent:0; text-indent:0px; font-size:14pt;\"><br /></p>\n"
            "<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:22pt; color:#990000;\">Ground Mounted Solar PV</span></p>\n"
            "<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:10pt; color:#990000;\">===================================================</span></p>\n"
            "<p style=\"  -qt-block-i"
                                    f"ndent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Number of Sites available: {n_solar_grnd}</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Number of sites built: {n_solar_grnd_built}</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Total Potential: {self.frpp_df['Ground Solar Power'].sum()*0.001:,.0f} MW</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Total Annual Generation: {self.frpp_df['Annual Ground Solar Power (kWh/yr)'].sum()*0.001:,.0f} MWh</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Percent of Agency Demand: {(self.frpp_df['Annual Ground Solar Power (kWh/yr)'].sum()*0.001)/self.agency_energy_data['Electricity (MWh)'].values[0]*100:.1f}%</span></p>\n" 
            "<p style=\"-qt-paragraph-type:empty;  -qt-block-indent:0; text-indent:0px; font-size:14pt;\"><br /></p>\n"
            "<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:22pt; color:#990000;\">Hydrogen Fuel Cell</span></p>\n"
            "<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:10pt; color:#990000;\">===================================================</span></p>\n"
            "<p style=\"  -qt-block-i"
                                    f"ndent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Number of Sites available: {n_fuelcell}</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Number of sites built: {n_fuelcell_built}</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Total Potential: {self.frpp_df['Fuel Cell (kW)'].sum()*0.001:,.0f} MW</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Total Annual Generation: {self.frpp_df['Annual Fuel Cell (kW/yr)'].sum()*0.001:,.0f} MWh</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Percent of Agency Demand: {(self.frpp_df['Annual Fuel Cell (kW/yr)'].sum()*0.001)/self.agency_energy_data['Electricity (MWh)'].values[0]*100:.1f}%</span></p>\n"
            "<p style=\"-qt-paragraph-type:empty;  -qt-block-indent:0; text-indent:0px; font-size:14pt;\"><br /></p>\n"
            "<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:22pt; color:#990000;\">Geothermal Power</span></p>\n"
            "<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:10pt; color:#990000;\">===================================================</span></p>\n"
            "<p style=\"  -qt-block-i"
                                    f"ndent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Number of Sites available: {n_geot}</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Number of sites built: {n_geot_built}</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Total Potential: {self.frpp_df['Geothermal Power (kW)'].sum()*0.001:,.0f} MW</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Total Annual Generation: {self.frpp_df['Annual Geothermal Power (kWh/yr)'].sum()*0.001:,.0f} MWh</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Percent of Agency Demand: {(self.frpp_df['Annual Geothermal Power (kWh/yr)'].sum()*0.001)/self.agency_energy_data['Electricity (MWh)'].values[0]*100:.1f}%</span></p>\n"
            "<p style=\"-qt-paragraph-type:em"
                                    "pty;  -qt-block-indent:0; text-indent:0px; font-size:14pt; color:#000000;\"><br /></p>\n"
            "<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:22pt; color:#990000;\">Concentrating Solar</span></p>\n"
            "<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:10pt; color:#990000;\">===================================================</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Number of Sites available: {n_conc_sol}</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Number of sites built: {n_conc_sol_built}</span></p>\n"
            "<p "
                                    f"style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Total Potential: {self.frpp_df['Concentrating Solar Power (kW)'].sum()*0.001:,.0f} MW</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Total Annual Generation: {self.frpp_df['Annual Concentrating Solar Power (kWh)'].sum()*0.001:,.0f} MWh</span></p>\n"
            f"<p style=\"  -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; color:#000000;\">Percent of Agency Demand: {(self.frpp_df['Annual Concentrating Solar Power (kWh)'].sum()*0.001)/self.agency_energy_data['Electricity (MWh)'].values[0]*100:.1f}%</span></p></body></html>", None))     

    def setup_econ_page(self)->None:  
        # Block Signals
        self.toggle_econ_signals(True)
             
        tot_eng = self.agency_energy_data['Electricity (MWh)'].values[0]                             
        tot_price = self.agency_price_data['Electricity'].values[0]
        # Establish CFE dropdown menu   
        for key, value in {'Wind Power (kW)': 'Wind', 'Rooftop Solar Power': "Rooftop Solar", 
                    'Ground Solar Power': "Ground-Mounted Solar",
                    'Fuel Cell (kW)': "Hydrogen Fuel Cells", 'Geothermal Power (kW)': "Geothermal", 
                    'Concentrating Solar Power (kW)': "Concentrating Solar"}.items():
            if self.frpp_df[key].sum() > 0:
                self.econ_cfe_comboBox.addItem(value)        
        self.econ_cfe_comboBox.setCurrentIndex(0)
                
        cur_dict = self.model_assump[self.assump_from_econ_cfe()]['gov_rates']

        # set ownership to government
        self.gov_radioButton.setChecked(True)

        # set input values
        # Top input box
        self.econ_power_produced_lineEdit.setText(f"{self.frpp_df[self.df_header_from_econ_cfe()].sum()*0.001:.2f}")
        self.econ_energy_storage_lineEdit.setText("0")
        self.econ_bat_rep_lineEdit.setText("5")
        self.econ_calc_engstor_lineEdit.setText("0.0")
        self.econ_agency_cons_lineEdit.setText(str(tot_eng))
        self.econ_agency_elec_cost_lineEdit.setText(str(tot_price))
        self.econ_unit_cost_lineEdit.setText(f"{tot_price/(tot_eng*1000):.6f}")

        # Bottom input box
        self.econ_interest_rate_lineEdit.setText("2.0")
        self.econ_proj_life_lineEdit.setText(str(cur_dict['project_life']))
        if "Solar" in self.econ_cfe_comboBox.currentText() or "Hydrogen" in self.econ_cfe_comboBox.currentText():
            self.econ_ann_deg_lineEdit.setText("0.5")
        else:
            self.econ_ann_deg_lineEdit.setEnabled(False)
        self.econ_rec_lineEdit.setText("19")        
        self.econ_cont_per_lineEdit.setText("2")
        self.econ_cont_per_Slider.setValue(2)
        self.econ_tax_rate_lineEdit.setText(str(cur_dict['taxation']))
        self.econ_tax_rate_Slider.setValue(int(cur_dict['taxation']))
        self.econ_insur_rate_lineEdit.setText("0.5")
        self.econ_insur_rate_Slider.setValue(5) # increments of 0.1
        self.econ_inflat_rate_lineEdit.setText(str(cur_dict['inflation_rate']))
        self.econ_inflat_rate_Slider.setValue(max(1,int(cur_dict['inflation_rate']*10))) # increments of 0.1
        self.econ_fed_tax_cred_lineEdit.setText("0.0")
        self.econ_fed_tax_cred_lineEdit.setEnabled(False)
        self.econ_fed_tax_cred_Slider.setValue(0) 
        self.econ_fed_tax_cred_Slider.setEnabled(False)
        self.econ_tax_cred_per_lineEdit.setText("0.0")
        self.econ_tax_cred_per_lineEdit.setEnabled(False)
        self.econ_tax_cred_per_Slider.setValue(1) 
        self.econ_tax_cred_per_Slider.setEnabled(False)
        esc_val = 0.5 * (cur_dict['electric_cost_escalation_lb'] + cur_dict['electric_cost_escalation_ub'])
        self.econ_price_esc_lineEdit.setText(str(esc_val))
        self.econ_price_esc_Slider.setValue(24) 

        self.toggle_econ_signals(False)
   
    def df_header_from_econ_cfe(self, req="Power"):        
        cfe_text = self.econ_cfe_comboBox.currentText()
        if req == "Power":
            if cfe_text == "Rooftop Solar":
                return 'Rooftop Solar Power'
            elif cfe_text == "Ground-Mounted Solar": 
                return 'Ground Solar Power'
            elif cfe_text == "Wind": 
                return 'Wind Power (kW)'
            elif cfe_text == "Geothermal":
                return 'Geothermal Power (kW)'
            elif cfe_text == "Hydrogen Fuel Cells":
                return 'Fuel Cell (kW)'
            elif cfe_text == "Concentrating Solar":
                return 'Concentrating Solar Power (kW)'
        elif req == "Energy":
            if cfe_text == "Rooftop Solar":
                return 'Annual Rooftop Solar Power (kWh/yr)'
            elif cfe_text == "Ground-Mounted Solar": 
                return 'Annual Ground Solar Power (kWh/yr)'
            elif cfe_text == "Wind": 
                return 'Annual Wind Power (kWh/yr)'
            elif cfe_text == "Geothermal":
                return 'Annual Geothermal Power (kWh/yr)'
            elif cfe_text == "Hydrogen Fuel Cells":
                return 'Annual Fuel Cell (kW/yr)'
            elif cfe_text == "Concentrating Solar":
                return 'Annual Concentrating Solar Power (kWh)'
        
    def assump_from_econ_cfe(self):        
        cfe_text = self.econ_cfe_comboBox.currentText()
        if cfe_text == "Rooftop Solar":
            return 'solar'
        elif cfe_text == "Ground-Mounted Solar": 
            return 'solar'
        elif cfe_text == "Wind": 
            return 'wind'
        elif cfe_text == "Geothermal":
            return 'geo_therm'
        elif cfe_text == "Hydrogen Fuel Cells":
            return 'hydrogen'
        elif cfe_text == "Concentrating Solar":
            return 'conc_solar'

    def toggle_econ_signals(self, toggle:bool):
        self.econ_cfe_comboBox.blockSignals(toggle)
        self.gov_radioButton.blockSignals(toggle)
        self.third_party_radioButton.blockSignals(toggle)
        self.econ_energy_storage_lineEdit.blockSignals(toggle)
        self.econ_bat_rep_lineEdit.blockSignals(toggle)
        self.econ_agency_cons_lineEdit.blockSignals(toggle)
        self.econ_agency_elec_cost_lineEdit.blockSignals(toggle)
        self.econ_interest_rate_lineEdit.blockSignals(toggle)
        self.econ_proj_life_lineEdit.blockSignals(toggle)
        self.econ_marcs_comboBox.blockSignals(toggle)
        self.econ_ann_deg_lineEdit.blockSignals(toggle)
        self.econ_rec_lineEdit.blockSignals(toggle)
        self.econ_cont_per_lineEdit.blockSignals(toggle)
        self.econ_tax_rate_lineEdit.blockSignals(toggle)
        self.econ_insur_rate_lineEdit.blockSignals(toggle)
        self.econ_inflat_rate_lineEdit.blockSignals(toggle)
        self.econ_fed_tax_cred_lineEdit.blockSignals(toggle)
        self.econ_tax_cred_per_lineEdit.blockSignals(toggle)
        self.econ_price_esc_lineEdit.blockSignals(toggle)
        self.econ_cont_per_Slider.blockSignals(toggle)
        self.econ_tax_rate_Slider.blockSignals(toggle)
        self.econ_insur_rate_Slider.blockSignals(toggle)
        self.econ_inflat_rate_Slider.blockSignals(toggle)
        self.econ_fed_tax_cred_Slider.blockSignals(toggle)
        self.econ_tax_cred_per_Slider.blockSignals(toggle)
        self.econ_price_esc_Slider.blockSignals(toggle)

    def set_sliders(self):
        self.toggle_econ_signals(True)

        cont_per = int(float(self.econ_cont_per_lineEdit.text()))
        self.econ_cont_per_Slider.setValue(cont_per)

        tax_rate = int(float(self.econ_tax_rate_lineEdit.text()))
        self.econ_tax_rate_Slider.setValue(tax_rate)

        insur_rate = int(10*float(self.econ_insur_rate_lineEdit.text()))
        self.econ_insur_rate_Slider.setValue(insur_rate) # increments of 0.1

        inflat_rate = int(10*float(self.econ_inflat_rate_lineEdit.text()))
        self.econ_inflat_rate_Slider.setValue(max(1,inflat_rate)) # increments of 0.1

        if self.econ_fed_tax_cred_lineEdit.isEnabled():
            fed_tax_cred = int(float(self.econ_fed_tax_cred_lineEdit.text()))        
            self.econ_fed_tax_cred_Slider.setValue(fed_tax_cred) 

        if self.econ_tax_cred_per_lineEdit.isEnabled():
            tax_cred_per = int(float(self.econ_tax_cred_per_lineEdit.text()))
            self.econ_tax_cred_per_Slider.setValue(max(1,tax_cred_per)) 
              
        if self.gov_radioButton.isChecked():
            cur_dict = self.model_assump[self.assump_from_econ_cfe()]['gov_rates']
        else:
            cur_dict = self.model_assump[self.assump_from_econ_cfe()]['third_party_rates']

        esc_val = int((float(self.econ_price_esc_lineEdit.text()) - cur_dict['electric_cost_escalation_lb'])/
                      ((cur_dict['electric_cost_escalation_ub']-cur_dict['electric_cost_escalation_lb'])/50))         
        self.econ_price_esc_Slider.setValue(esc_val) 

        self.toggle_econ_signals(False)

    def agency_code(self):
        ind = self.agency_comboBox.currentIndex()
        code = {0: 'NA', 1: "USDA", 2: "DOC", 3: "USACE",
                4: "CSOSA", 5: "DOE", 6: "EPA", 7: "EOP",
                8: "FRTIB", 9: "GSA", 10: "HHS", 11: "DHS",
                12: "IND", 13: "DOI", 14: "DOJ", 15: "DOL",
                16: "MSPB", 17: "NASA", 18: "NCUA", 19: "PSA",
                20: "FRCA", 21: "SMITH", 22: "DOS", 23: "USAID",
                24: "TVA", 25: "DOT", 26: "DFC", 27: "USHMM",
                28: "USAGM", 29: "VA"}
        return code[ind]
    
    def agency_set(self, agency):        
        code = {"NA":0, "USDA":1, "DOC": 2, "USACE": 3,
                "CSOSA":4, "DOE":5, "EPA": 6, "EOP": 7,
                "FRTIB":8, "GSA":9, "HHS":10, "DHS":11,
                 "IND":12, "DOI":13, "DOJ":14, "DOL":15,
                 "MSPB":16, "NASA":17, "NCUA":18,"PSA":19,
                 "FRCA":20, "SMITH":21, "DOS":22, "USAID":23,
                 "TVA":24,  "DOT":25, "DFC":26, "USHMM":27,
                 "USAGM":28, "VA":29 }
        if agency in code.keys():
            self.agency_comboBox.setCurrentIndex(code[agency])
            return 0
        else:
            return -1
        


if __name__ == "__main__":
    app = QApplication(sys.argv)
    ui = MainWindow()
    ui.show()
    sys.exit(app.exec())