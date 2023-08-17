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
        
        print("test")


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

        




if __name__ == "__main__":
    app = QApplication(sys.argv)
    ui = MainWindow()
    ui.show()
    sys.exit(app.exec())