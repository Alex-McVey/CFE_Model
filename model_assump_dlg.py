from PyQt6.QtWidgets import *
from PyQt6.QtGui import *
from PyQt6.uic import loadUi
import pathlib
import numpy as np

PROJECT_PATH = pathlib.Path(__file__).parent
ASSUMP_DLG = PROJECT_PATH / "model_assumptions.ui"
DEFAULT_WHITE = u'font: 10pt "Segoe UI";background-color: rgb(255, 255, 255);'
DEFAULT_ERROR = u'font: 10pt "Segoe UI";background-color: rgb(255, 103, 103);'

class model_assump_dlg(QDialog):
    def __init__(self, assump_json, color_red=[]):
        super(model_assump_dlg, self).__init__()
        loadUi(ASSUMP_DLG, self)
        self.assump_json = assump_json
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
        self.opt_eff_lineEdit.setText(str(assump_json['conc_solar']['system']['optical_eff']))
        self.elec_eff_lineEdit.setText(str(assump_json['conc_solar']['system']['elec_eff'])) 
        self.therm_stor_lineEdit.setText(str(assump_json['conc_solar']['system']['thermal_storage']))
        self.avg_sunhrs_lineEdit.setText(str(assump_json['conc_solar']['system']['avg_sun_hours'])) 
        self.avg_sol_rad_lineEdit.setText(str(assump_json['conc_solar']['system']['avg_solar_rad']))
        self.per_agency_land_lineEdit.setText(str(assump_json['conc_solar']['system']['perc_land_used'])) 		
        self.cs_cap_fac_lineEdit.setText(str(assump_json['conc_solar']['system']['capacity_factor'])) 
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
        self.cell_volt_lineEdit.setText(str(assump_json['hydrogen']['system']['cell_voltage'])) 
        self.cell_act_area_lineEdit.setText(str(assump_json['hydrogen']['system']['cell_active_area'])) 
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
        # valid and message class values
        self.valid = True
        self.errors = []
        # other class var
        self.color_red = color_red.copy()
        self.line_edits = [self.wind_frac_lineEdit, self.wind_turb_diam_lineEdit, self.wind_CapFac_lineEdit, self.wind_Nturb_lineEdit, 
                           self.wind_spacing_lineEdit, self.wind_pow_coef_lineEdit, self.wind_eff_lineEdit, self.wg_disc_rate_lineEdit, 
                           self.wg_infl_rate_lineEdit, self.wg_ece_LB_lineEdit, self.wg_ece_UB_lineEdit, self.wg_tax_lineEdit, 
                           self.wg_proj_life_lineEdit, self.wt_disc_rate_lineEdit, self.wt_infl_rate_lineEdit, self.wt_ece_LB_lineEdit, 
                           self.wt_ece_UB_lineEdit, self.wt_tax_lineEdit, self.wt_proj_life_lineEdit, self.per_roof_lineEdit, 
                           self.rad_2_elec_lineEdit, self.fract_avg_SolRad_lineEdit, self.perc_land_used_lineEdit, self.sun_hours_lineEdit, 
                           self.dc_ac_lineEdit, self.cap_fac_lineEdit, self.sg_disc_rate_lineEdit, self.sg_infl_rate_lineEdit, 
                           self.sg_ece_LB_lineEdit, self.sg_ece_UB_lineEdit, self.sg_tax_lineEdit, self.sg_proj_life_lineEdit, 
                           self.st_disc_rate_lineEdit, self.st_infl_rate_lineEdit, self.st_ece_LB_lineEdit, self.st_ece_UB_lineEdit, 
                           self.st_tax_lineEdit, self.st_proj_life_lineEdit, self.opt_eff_lineEdit, self.elec_eff_lineEdit, self.therm_stor_lineEdit, 
                           self.avg_sunhrs_lineEdit, self.avg_sol_rad_lineEdit, self.per_agency_land_lineEdit, self.cs_cap_fac_lineEdit, 
                           self.csg_disc_rate_lineEdit, self.csg_infl_rate_lineEdit, self.csg_ece_LB_lineEdit, self.csg_ece_UB_lineEdit, 
                           self.csg_tax_lineEdit, self.csg_proj_life_lineEdit, self.cst_disc_rate_lineEdit, self.cst_infl_rate_lineEdit, 
                           self.cst_ece_LB_lineEdit, self.cst_ece_UB_lineEdit, self.cst_tax_lineEdit, self.cst_proj_life_lineEdit, self.h_cap_fac_lineEdit, 
                           self.hrs_oper_lineEdit, self.fuel_emp_comboBox, self.cur_den_lineEdit, self.cell_volt_lineEdit, self.cell_act_area_lineEdit, 
                           self.Nstacks_lineEdit, self.prot_exch_comboBox,  self.heat_rec_comboBox,  self.heat_val_lineEdit, self.stack_replace_lineEdit, 
                           self.annual_degrade_lineEdit, self.hg_disc_rate_lineEdit, self.hg_infl_rate_lineEdit, self.hg_ece_LB_lineEdit, 
                           self.hg_ece_UB_lineEdit, self.hg_tax_lineEdit, self.hg_proj_life_lineEdit, self.ht_disc_rate_lineEdit, self.ht_infl_rate_lineEdit, 
                           self.ht_ece_LB_lineEdit, self.ht_ece_UB_lineEdit, self.ht_tax_lineEdit, self.ht_proj_life_lineEdit, self.well_diam_lineEdit, 
                           self.fluid_vel_lineEdit, self.well_depth_lineEdit, self.turb_outlet_lineEdit, self.geo_cap_fac_lineEdit, self.fluid_den_lineEdit, 
                           self.pump_eff_lineEdit, self.gg_disc_rate_lineEdit, self.gg_infl_rate_lineEdit, self.gg_ece_LB_lineEdit, self.gg_ece_UB_lineEdit, 
                           self.gg_tax_lineEdit, self.gg_proj_life_lineEdit, self.gt_disc_rate_lineEdit, self.gt_infl_rate_lineEdit, self.gt_ece_LB_lineEdit, 
                           self.gt_ece_UB_lineEdit, self.gt_tax_lineEdit, self.gt_proj_life_lineEdit]
        # Reset all page options to valid background
        for qq, obj in enumerate(self.line_edits):
            if qq in self.color_red:
                obj.setStyleSheet(DEFAULT_ERROR)  
            else:
                obj.setStyleSheet(DEFAULT_WHITE)
        # Set functions
        self.ok_button = self.buttonBox.button(QDialogButtonBox.StandardButton.Ok)
        self.ok_button.clicked.connect(self.validate)
    
    def validate(self):
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
        color_red = []
       
        # Check wind
        diam = np.nan
        try:
            val = float(self.wind_frac_lineEdit.text()) 
            self.assump_json['wind']['system']['fraction_of_average_wind_speed'] = val
            if val < 0:
                color_red.append(self.wind_frac_lineEdit)
                valid = False
                errors.append(f"The Fraction of Averate Wind Speed must be greater than zero")
            elif val > 3:
                color_red.append(self.wind_frac_lineEdit)
                valid = False
                errors.append(f"The Fraction of Averate Wind Speed must be less than 3")   
        except:
            color_red.append(self.wind_frac_lineEdit)
            valid = False
            errors.append(f"The Fraction of Averate Wind Speed must be numeric")

        try:
            val = float(self.wind_turb_diam_lineEdit.text()) 
            self.assump_json['wind']['system']['turbine_diameter'] = val
            diam = val
            if val < 9:
                color_red.append(self.wind_turb_diam_lineEdit)
                valid = False
                errors.append(f"The Turbine Rotor Diameter must be greater than or equal to 9 meters")
            elif val > 252:
                color_red.append(self.wind_turb_diam_lineEdit)
                valid = False
                errors.append(f"The Turbine Rotor Diameter must be less than or equal to 252 meters")   
        except:
            color_red.append(self.wind_turb_diam_lineEdit)
            valid = False
            errors.append(f"The Turbine Rotor Diameter must be numeric")

        try:
            val = float(self.wind_CapFac_lineEdit.text()) 
            self.assump_json['wind']['system']['capacity_factor'] = val
            if val <= 0:
                color_red.append(self.wind_CapFac_lineEdit)
                valid = False
                errors.append(f"The Wind Capacity Factor must be greater than 0%")
            elif val > 100:
                color_red.append(self.wind_CapFac_lineEdit)
                valid = False
                errors.append(f"The Wind Capacity Factor must be less than or equal to 100%")   
        except:
            color_red.append(self.wind_CapFac_lineEdit)
            valid = False
            errors.append(f"The Wind Capacity Factor must be numeric")

        # I think the spacing is wrong. Nothing supports a 2*diameter. 5 was the lowest https://www.renewablesfirst.co.uk/renewable-energy-technologies/windpower/community-windpower/location-size-no-of-wind-turbines/#:~:text=The%20number%20of%20wind%20turbines,turbine%20it%20is%20410%20metres.
        # I saw and SAM has it as 8 https://github.com/NREL/SAM/blob/85b6974a7f9f02e931b9013b7d250a7eccd87b1c/deploy/runtime/ui/Wind%20Farm%20Specifications.txt#L1631C30-L1631C30
        try:
            val = float(self.wind_Nturb_lineEdit.text())  
            self.assump_json['wind']['system']['N_turbines']            
        except:
            color_red.append(self.wind_Nturb_lineEdit)
            valid = False
            errors.append(f"The number of turbines must be numeric")
        
        try:
            val = float(self.wind_spacing_lineEdit.text())  
            self.assump_json['wind']['system']['turbine_spacing'] = val             
        except:
            color_red.append(self.wind_spacing_lineEdit)
            valid = False
            errors.append(f"The distance between turbines must be numeric")

        if not np.isnan(diam):
            try:
                nturb = float(self.wind_Nturb_lineEdit.text())
                req_spa = float(self.wind_spacing_lineEdit.text())  
                spacing = (req_spa/nturb)/diam
                if spacing < 5:
                    color_red.append(self.wind_Nturb_lineEdit)
                    color_red.append(self.wind_spacing_lineEdit)
                    valid = False
                    errors.append(f"The distance between turbines must be greater than or equal to 5 turbine diameters. Currently: {spacing}")
            except:
                pass
        
        try:
            val = float(self.wind_eff_lineEdit.text()) 
            self.assump_json['wind']['system']['overall_efficiency'] = val
            if val <= 0:
                color_red.append(self.wind_eff_lineEdit)
                valid = False
                errors.append(f"The Overall Electrical Efficiency must be greater than 0%")
            elif val > 100:
                color_red.append(self.wind_eff_lineEdit)
                valid = False
                errors.append(f"The Overall Electrical Efficiency must be less than or equal to 100%")   
        except:
            color_red.append(self.wind_eff_lineEdit)
            valid = False
            errors.append(f"The Overall Electrical Efficiency must be numeric")
        
        try:
            val = float(self.wind_pow_coef_lineEdit.text()) 
            self.assump_json['wind']['system']['power_coefficient'] = val
            if val <= 0:
                color_red.append(self.wind_pow_coef_lineEdit)
                valid = False
                errors.append(f"The wind power coefficient must be greater than 0")
            elif val > 1:
                color_red.append(self.wind_pow_coef_lineEdit)
                valid = False
                errors.append(f"The wind power coefficient must be less than or equal to 1")   
        except:
            color_red.append(self.wind_pow_coef_lineEdit)
            valid = False
            errors.append(f"The wind power coefficient must be numeric")
        # Check Solar
        try:
            val = float(self.per_roof_lineEdit.text()) 
            self.assump_json['solar']['system']['precent_roof_avail'] = val
            if val <= 0:
                color_red.append(self.per_roof_lineEdit)
                valid = False
                errors.append(f"The percent of rooftop available for solar must be greater than 0%")
            elif val > 100:
                color_red.append(self.per_roof_lineEdit)
                valid = False
                errors.append(f"The percent of rooftop available for solar must be less than or equal to 100%")   
        except:
            color_red.append(self.per_roof_lineEdit)
            valid = False
            errors.append(f"The percent of rooftop available for solar must be numeric")

        try:
            val = float(self.rad_2_elec_lineEdit.text()) 
            self.assump_json['solar']['system']['per_rad_to_elec'] = val
            if val <= 0:
                color_red.append(self.rad_2_elec_lineEdit)
                valid = False
                errors.append(f"The percent of radiation converted to electricity for solar must be greater than 0%")
            elif val > 100:
                color_red.append(self.rad_2_elec_lineEdit)
                valid = False
                errors.append(f"The percent of radiation converted to electricity for solar must be less than or equal to 100%")   
        except:
            color_red.append(self.rad_2_elec_lineEdit)
            valid = False
            errors.append(f"The percent of radiation converted to electricity for solar must be numeric")

        try:
            val = float(self.fract_avg_SolRad_lineEdit.text()) 
            self.assump_json['solar']['system']['avg_radiation_frac'] = val 
            if val <= 0:
                color_red.append(self.fract_avg_SolRad_lineEdit)
                valid = False
                errors.append(f"The fraction of average solar radiation must be greater than 0")
            elif val > 1:
                color_red.append(self.fract_avg_SolRad_lineEdit)
                valid = False
                errors.append(f"The fraction of average solar radiation must be less than or equal to 1")   
        except:
            color_red.append(self.fract_avg_SolRad_lineEdit)
            valid = False
            errors.append(f"The fraction of average solar radiation must be numeric")
        
        try:
            val = float(self.perc_land_used_lineEdit.text()) 
            self.assump_json['solar']['system']['perc_land_used'] = val 
            if val <= 0:
                color_red.append(self.perc_land_used_lineEdit)
                valid = False
                errors.append(f"The percent of land used for solar must be greater than 0%")
            elif val > 100:
                color_red.append(self.perc_land_used_lineEdit)
                valid = False
                errors.append(f"The percent of land used for solar must be less than or equal to 100%")   
        except:
            color_red.append(self.perc_land_used_lineEdit)
            valid = False
            errors.append(f"The percent of land used for solar must be numeric")

        try:
            val = float(self.sun_hours_lineEdit.text()) 
            self.assump_json['solar']['system']['avg_sun_hours'] = val
            if val < 1:
                color_red.append(self.sun_hours_lineEdit)
                valid = False
                errors.append(f"The average number of sun hours per day must be greater than or equal to 1")
            elif val > 24:
                color_red.append(self.sun_hours_lineEdit)
                valid = False
                errors.append(f"The average number of sun hours per day must be less than or equal to 24")   
        except:
            color_red.append(self.sun_hours_lineEdit)
            valid = False
            errors.append(f"The average number of sun hours per day must be numeric")

        try:
            val = float(self.dc_ac_lineEdit.text()) 
            self.assump_json['solar']['system']['dc_to_ac_ratio'] = val 
            if val < 1:
                color_red.append(self.dc_ac_lineEdit)
                valid = False
                errors.append(f"The D.C. to A.C. ratio must be greater than or equal to 1")
            elif val > 1.35:
                color_red.append(self.dc_ac_lineEdit)
                valid = False
                errors.append(f"The D.C. to A.C. ratio must be less than or equal to 1.35")   
        except:
            color_red.append(self.dc_ac_lineEdit)
            valid = False
            errors.append(f"The D.C. to A.C. ratio must be numeric")

        try:
            val = float(self.cap_fac_lineEdit.text()) 
            self.assump_json['solar']['system']['capacity_factor'] = val 
            if val <= 0:
                color_red.append(self.cap_fac_lineEdit)
                valid = False
                errors.append(f"The solar capacity factor must be greater than 0%")
            elif val > 100:
                color_red.append(self.cap_fac_lineEdit)
                valid = False
                errors.append(f"The solar capacity factor must be less than or equal to 100%")   
        except:
            color_red.append(self.cap_fac_lineEdit)
            valid = False
            errors.append(f"The solar capacity factor must be numeric")  
        # Check Concentrating Solar
        try:
            val = float(self.opt_eff_lineEdit.text()) 
            self.assump_json['conc_solar']['system']['optical_eff'] = val
            if val <= 0:
                color_red.append(self.opt_eff_lineEdit)
                valid = False
                errors.append(f"The optical efficiency for concentrating solar must be greater than 0%")
            elif val > 100:
                color_red.append(self.opt_eff_lineEdit)
                valid = False
                errors.append(f"The optical efficiency for concentrating solar must be less than or equal to 100%")   
        except:
            color_red.append(self.opt_eff_lineEdit)
            valid = False
            errors.append(f"The optical efficiency for concentrating solar must be numeric")  

        try:
            val = float(self.elec_eff_lineEdit.text())
            self.assump_json['conc_solar']['system']['elec_eff'] = val 
            if val <= 0:
                color_red.append(self.elec_eff_lineEdit)
                valid = False
                errors.append(f"The electrical efficiency for concentrating solar must be greater than 0%")
            elif val > 100:
                color_red.append(self.elec_eff_lineEdit)
                valid = False
                errors.append(f"The electrical efficiency for concentrating solar must be less than or equal to 100%")   
        except:
            color_red.append(self.elec_eff_lineEdit)
            valid = False
            errors.append(f"The electrical efficiency for concentrating solar must be numeric") 
        
        try:
            val = float(self.therm_stor_lineEdit.text()) 
            self.assump_json['conc_solar']['system']['thermal_storage'] = val 
            if val < 0:
                color_red.append(self.therm_stor_lineEdit)
                valid = False
                errors.append(f"The thermal storage for concentrating solar must be greater than or equal to 0")
            elif val > 30:
                color_red.append(self.therm_stor_lineEdit)
                valid = False
                errors.append(f"The thermal storage for concentrating solar must be less than or equal to 30")   
        except:
            color_red.append(self.therm_stor_lineEdit)
            valid = False
            errors.append(f"The thermal storage for concentrating solar must be numeric")
        
        try:
            val = float(self.avg_sunhrs_lineEdit.text()) 
            self.assump_json['conc_solar']['system']['avg_sun_hours'] = val
            if val < 0:
                color_red.append(self.avg_sunhrs_lineEdit)
                valid = False
                errors.append(f"The average hours of sun per day for concentrating solar must be greater than or equal to 0")
            elif val > 30:
                color_red.append(self.avg_sunhrs_lineEdit)
                valid = False
                errors.append(f"The average hours of sun per day for concentrating solar must be less than or equal to 24")   
        except:
            color_red.append(self.avg_sunhrs_lineEdit)
            valid = False
            errors.append(f"The average hours of sun per day for concentrating solar must be numeric")
        
        try:
            val = float(self.avg_sol_rad_lineEdit.text()) 
            self.assump_json['conc_solar']['system']['avg_solar_rad'] = val 
            if val < 0.1:
                color_red.append(self.avg_sol_rad_lineEdit)
                valid = False
                errors.append(f"The average solar radiation must be greater than or equal to 0.1")
            elif val > 1.115:
                color_red.append(self.avg_sol_rad_lineEdit)
                valid = False
                errors.append(f"The average solar radiation must be less than or equal to 1.115")   
        except:
            color_red.append(self.avg_sol_rad_lineEdit)
            valid = False
            errors.append(f"The average solar radiation must be numeric")
        
        try:
            val = float(self.per_agency_land_lineEdit.text()) 
            self.assump_json['conc_solar']['system']['perc_land_used'] = val 
            if val <= 0:
                color_red.append(self.per_agency_land_lineEdit)
                valid = False
                errors.append(f"The precent of land used for concentrating solar must be greater than 0%")
            elif val > 100:
                color_red.append(self.per_agency_land_lineEdit)
                valid = False
                errors.append(f"The precent of land used for concentrating solar must be less than or equal to 100%")   
        except:
            color_red.append(self.per_agency_land_lineEdit)
            valid = False
            errors.append(f"The precent of land used for concentrating solar must be numeric") 
        		
        try:
            val = float(self.cs_cap_fac_lineEdit.text()) 
            self.assump_json['conc_solar']['system']['capacity_factor'] = val 
            if val <= 0:
                color_red.append(self.cs_cap_fac_lineEdit)
                valid = False
                errors.append(f"The precent of land used for concentrating solar must be greater than 0%")
            elif val > 100:
                color_red.append(self.cs_cap_fac_lineEdit)
                valid = False
                errors.append(f"The precent of land used for concentrating solar must be less than or equal to 100%")   
        except:
            color_red.append(self.cs_cap_fac_lineEdit)
            valid = False
            errors.append(f"The precent of land used for concentrating solar must be numeric") 
        # Check Hydrogen Cells
        try:
            val = float(self.h_cap_fac_lineEdit.text())
            self.assump_json['hydrogen']['system']['capacity_factor'] = val 
            if val <= 0:
                color_red.append(self.h_cap_fac_lineEdit)
                valid = False
                errors.append(f"The capacity factor for hydrogen fuel cells must be greater than 0%")
            elif val > 100:
                color_red.append(self.h_cap_fac_lineEdit)
                valid = False
                errors.append(f"The capacity factor for hydrogen fuel cells must be less than or equal to 100%")   
        except:
            color_red.append(self.h_cap_fac_lineEdit)
            valid = False
            errors.append(f"The capacity factor for hydrogen fuel cells must be numeric") 
        
        try:
            val = float(self.hrs_oper_lineEdit.text()) 
            self.assump_json['hydrogen']['system']['daily_operation'] = val
            if val < 0:
                color_red.append(self.hrs_oper_lineEdit)
                valid = False
                errors.append(f"The operational hours per day for hydrogen fuel cells must be greater than or equal to 0")
            elif val > 30:
                color_red.append(self.hrs_oper_lineEdit)
                valid = False
                errors.append(f"The operational hours per day for hydrogen fuel cells must be less than or equal to 24")   
        except:
            color_red.append(self.hrs_oper_lineEdit)
            valid = False
            errors.append(f"The operational hours per day for hydrogen fuel cells must be numeric")
        
        try:
            val = float(self.cur_den_lineEdit.text()) 
            self.assump_json['hydrogen']['system']['curr_density'] = val
            if val < 0:
                color_red.append(self.cur_den_lineEdit)
                valid = False
                errors.append(f"The hydrogen fuel cell current density must be greater than or equal to 0")
            elif val > 15000:
                color_red.append(self.cur_den_lineEdit)
                valid = False
                errors.append(f"The hydrogen fuel cell current density must be less than or equal to 15,000")   
        except:
            color_red.append(self.cur_den_lineEdit)
            valid = False
            errors.append(f"The hydrogen fuel cell current density must be numeric")
        
        try:
            val = float(self.cell_act_area_lineEdit.text()) 
            self.assump_json['hydrogen']['system']['cell_active_area'] = val            
            if val < 0:
                color_red.append(self.cell_act_area_lineEdit)
                valid = False
                errors.append(f"The hydrogen fuel cell active area must be greater than or equal to 0")
            elif val > 5:
                color_red.append(self.cell_act_area_lineEdit)
                valid = False
                errors.append(f"The hydrogen fuel cell active area must be less than or equal to 5")   
        except:
            color_red.append(self.cell_act_area_lineEdit)
            valid = False
            errors.append(f"The hydrogen fuel cell active area must be numeric")

        try:
            val = float(self.cell_volt_lineEdit.text()) 
            self.assump_json['hydrogen']['system']['cell_voltage'] = val
            if val < 0.6:
                color_red.append(self.cell_volt_lineEdit)
                valid = False
                errors.append(f"The hydrogen fuel cell voltage must be greater than or equal to 0.6")
            elif val > 0.8:
                color_red.append(self.cell_volt_lineEdit)
                valid = False
                errors.append(f"The hydrogen fuel cell voltage must be less than or equal to 0.8")   
        except:
            color_red.append(self.cell_volt_lineEdit)
            valid = False
            errors.append(f"The hydrogen fuel cell voltage must be numeric")

        try:
            val = int(self.Nstacks_lineEdit.text()) 
            self.assump_json['hydrogen']['system']['N_stacks'] = val
            if val < 1:
                color_red.append(self.Nstacks_lineEdit)
                valid = False
                errors.append(f"The number of stacks for hydrogen fuel cell must be greater than or equal to 1")
            elif val > 100:
                color_red.append(self.Nstacks_lineEdit)
                valid = False
                errors.append(f"The number of stacks for hydrogen fuel cell must be less than or equal to 100")   
        except:
            color_red.append(self.Nstacks_lineEdit)
            valid = False
            errors.append(f"The number of stacks for hydrogen fuel cell must be numeric")
        
        try:
            val = float(self.stack_replace_lineEdit.text()) 
            self.assump_json['hydrogen']['system']['stack_replacement'] = val 
            if val < 1:
                color_red.append(self.stack_replace_lineEdit)
                valid = False
                errors.append(f"The time to replace hydrogen fuel cell must be greater than or equal to 1")
            elif val > 10:
                color_red.append(self.stack_replace_lineEdit)
                valid = False
                errors.append(f"The time to replace hydrogen fuel cell must be less than or equal to 10")   
        except:
            color_red.append(self.stack_replace_lineEdit)
            valid = False
            errors.append(f"The time to replace for hydrogen fuel cell must be numeric")       
        
        try:
            val = float(self.annual_degrade_lineEdit.text()) 
            self.assump_json['hydrogen']['system']['annual_degredation'] = val
            if val < 3:
                color_red.append(self.annual_degrade_lineEdit)
                valid = False
                errors.append(f"The hydrogen fuel cell degradation must be greater than or equal to 1")
            elif val > 20:
                color_red.append(self.annual_degrade_lineEdit)
                valid = False
                errors.append(f"The hydrogen fuel cell degradation must be less than or equal to 10")   
        except:
            color_red.append(self.annual_degrade_lineEdit)
            valid = False
            errors.append(f"The hydrogen fuel cell degradation must be numeric") 
        # Check Geothermal
        try:
            # https://workingincaes.inl.gov/SiteAssets/CAES%20Files/FORGE/inl_ext-16-38751%20GETEM%20User%20Manual%20Final.pdf
            val = float(self.well_diam_lineEdit.text()) 
            self.assump_json['geo_therm']['system']['well_diameter'] = val
            if val < 0.1778:
                color_red.append(self.well_diam_lineEdit)
                valid = False
                errors.append(f"The geothermal well diameter must be greater than or equal to 0.1778 m (7 inches)")
            elif val > 0.3397:
                color_red.append(self.well_diam_lineEdit)
                valid = False
                errors.append(f"The geothermal well diameter must be less than or equal to 0.3397 m (13 3/8 inches)")   
        except:
            color_red.append(self.well_diam_lineEdit)
            valid = False
            errors.append(f"The geothermal well diameter must be numeric") 
        
        try:
            val = float(self.fluid_vel_lineEdit.text()) 
            self.assump_json['geo_therm']['system']['fluid_velocity'] = val
            if val < 0.1:
                color_red.append(self.fluid_vel_lineEdit)
                valid = False
                errors.append(f"The geothermal well fluid velocity must be greater than or equal to 0.1 m/s")
            elif val > 50:
                color_red.append(self.fluid_vel_lineEdit)
                valid = False
                errors.append(f"The geothermal well fluid velocity must be less than or equal to 50 m/s")   
        except:
            color_red.append(self.fluid_vel_lineEdit)
            valid = False
            errors.append(f"The geothermal well fluid velocity must be numeric") 

        try:
            # https://utahforge.com/2021/03/29/did-you-know-geothermal-wells-can-be-highly-deviated-too/
            val = float(self.well_depth_lineEdit.text()) 
            self.assump_json['geo_therm']['system']['well_depth'] = val
            if val < 1000:
                color_red.append(self.well_depth_lineEdit)
                valid = False
                errors.append(f"The geothermal well depth must be greater than or equal to 1000 m")
            elif val > 6000:
                color_red.append(self.well_depth_lineEdit)
                valid = False
                errors.append(f"The geothermal well depth must be less than or equal to 6000 m")   
        except:
            color_red.append(self.well_depth_lineEdit)
            valid = False
            errors.append(f"The geothermal well depth must be numeric") 
        
        try:
            val = float(self.turb_outlet_lineEdit.text()) 
            self.assump_json['geo_therm']['system']['turb_outlet_temp'] = val
            if val < 50:
                color_red.append(self.turb_outlet_lineEdit)
                valid = False
                errors.append(f"The geothermal turbine outlet temperature must be greater than or equal to 50 C")
            elif val > 100:
                color_red.append(self.turb_outlet_lineEdit)
                valid = False
                errors.append(f"The geothermal turbine outlet temperature must be less than or equal to 100 C")   
        except:
            color_red.append(self.turb_outlet_lineEdit)
            valid = False
            errors.append(f"The geothermal turbine outlet temperature must be numeric") 
        
        try:
            val = float(self.geo_cap_fac_lineEdit.text()) 
            self.assump_json['geo_therm']['system']['capacity_factor'] = val 
            if val <= 0:
                color_red.append(self.geo_cap_fac_lineEdit)
                valid = False
                errors.append(f"The capacity factor for geothermal wells must be greater than 0%")
            elif val > 100:
                color_red.append(self.geo_cap_fac_lineEdit)
                valid = False
                errors.append(f"The capacity factor for geothermal wells must be less than or equal to 100%")   
        except:
            color_red.append(self.geo_cap_fac_lineEdit)
            valid = False
            errors.append(f"The capacity factor for geothermal wells must be numeric") 
        
        try:
            val = float(self.fluid_den_lineEdit.text()) 
            self.assump_json['geo_therm']['system']['fluid_density'] = val 
            if val <= 0:
                color_red.append(self.fluid_den_lineEdit)
                valid = False
                errors.append(f"The geothermal turbine fluid density must be greater than 0")             
        except:
            color_red.append(self.fluid_den_lineEdit)
            valid = False
            errors.append(f"The geothermal turbine fluid density must be numeric") 
        
        try:
            val = float(self.pump_eff_lineEdit.text()) 
            self.assump_json['geo_therm']['system']['overall_pump_efficiency'] = val
            if val <= 0:
                color_red.append(self.pump_eff_lineEdit)
                valid = False
                errors.append(f"The pumping efficieny for geothermal wells must be greater than 0%")
            elif val > 100:
                color_red.append(self.pump_eff_lineEdit)
                valid = False
                errors.append(f"The pumping efficieny for geothermal wells must be less than or equal to 100%")   
        except:
            color_red.append(self.pump_eff_lineEdit)
            valid = False
            errors.append(f"The pumping efficieny for geothermal wells must be numeric") 
        # Check Goverment Rates
        for entry in [(self.wg_disc_rate_lineEdit, "wind"), (self.sg_disc_rate_lineEdit, "solar"),
                      (self.csg_disc_rate_lineEdit, "conc_solar"), (self.hg_disc_rate_lineEdit, "hydrogen"),
                      (self.gg_disc_rate_lineEdit, "geo_therm")]:
            try:
                val = float(entry[0].text())    
                self.assump_json[entry[1]]['gov_rates']['discount_rate'] = val             
            except:
                color_red.append(entry)
                valid = False
                errors.append(f"Government owned discounted rate must be numeric")
        
        for entry in [(self.wg_infl_rate_lineEdit, "wind"), (self.sg_infl_rate_lineEdit, "solar"),
                      (self.csg_infl_rate_lineEdit, "conc_solar"), (self.hg_infl_rate_lineEdit, "hydrogen"),
                      (self.gg_infl_rate_lineEdit, "geo_therm")]:
            try:
                val = float(entry[0].text())    
                self.assump_json[entry[1]]['gov_rates']['inflation_rate'] = val                   
            except:
                color_red.append(entry)
                valid = False
                errors.append(f"Government owned inflation rate must be numeric")
        
        for entry in [(self.wg_ece_LB_lineEdit, "wind"), (self.sg_ece_LB_lineEdit, "solar"),
                      (self.csg_ece_LB_lineEdit, "conc_solar"), (self.hg_ece_LB_lineEdit, "hydrogen"),
                      (self.gg_ece_LB_lineEdit, "geo_therm")]:
            try:
                val = float(entry[0].text())    
                self.assump_json[entry[1]]['gov_rates']['electric_cost_escalation_lb'] = val                  
            except:
                color_red.append(entry)
                valid = False
                errors.append(f"Government owned electric cost escalation lower bound must be numeric")
        
        for entry in [(self.wg_ece_UB_lineEdit, "wind"), (self.sg_ece_UB_lineEdit, "solar"),
                      (self.csg_ece_UB_lineEdit, "conc_solar"), (self.hg_ece_UB_lineEdit, "hydrogen"),
                      (self.gg_ece_UB_lineEdit, "geo_therm")]:
            try:
                val = float(entry[0].text())    
                self.assump_json[entry[1]]['gov_rates']['electric_cost_escalation_ub'] = val                  
            except:
                color_red.append(entry)
                valid = False
                errors.append(f"Government owned electric cost escalation upper bound must be numeric")

        for entry in [[self.wg_ece_LB_lineEdit,self.wg_ece_UB_lineEdit], [self.sg_ece_LB_lineEdit,self.sg_ece_UB_lineEdit],
                      [self.csg_ece_LB_lineEdit,self.csg_ece_UB_lineEdit], [self.hg_ece_LB_lineEdit,self.hg_ece_UB_lineEdit],
                      [self.gg_ece_LB_lineEdit,self.gg_ece_UB_lineEdit]]:       
            try:
                val_lb = float(entry[0].text())
                val_ub = float(entry[1].text())
                if val_ub <= val_lb:
                    color_red.append(entry[0])
                    color_red.append(entry[1])
                    valid = False
                    errors.append(f"The goverment owned upper bound on electric cost escalation must be larger than the lower bound") 
            except:
                pass
        
        for entry in [(self.wg_tax_lineEdit, "wind"), (self.sg_tax_lineEdit, "solar"),
                      (self.csg_tax_lineEdit, "conc_solar"), (self.hg_tax_lineEdit, "hydrogen"),
                      (self.gg_tax_lineEdit, "geo_therm")]:
            try:
                val = float(entry[0].text())    
                self.assump_json[entry[1]]['gov_rates']['taxation'] = val                   
            except:
                color_red.append(entry)
                valid = False
                errors.append(f"overnment owned taxation rate must be numeric")
        
        for entry in [(self.wg_proj_life_lineEdit, "wind"), (self.sg_proj_life_lineEdit, "solar"),
                      (self.csg_proj_life_lineEdit, "conc_solar"), (self.hg_proj_life_lineEdit, "hydrogen"),
                      (self.gg_proj_life_lineEdit, "geo_therm")]:
            try:
                val = float(entry[0].text())    
                self.assump_json[entry[1]]['gov_rates']['project_life'] = val 
                if val <= 0:
                    color_red.append(entry)
                    valid = False
                    errors.append(f"Government owned project life must be greater than 0") 
            except:
                color_red.append(entry)
                valid = False
                errors.append(f"Government owned project life must be numeric")
        

        # Check Third-Party Rates
        for entry in [(self.wt_disc_rate_lineEdit, "wind"), (self.st_disc_rate_lineEdit, "solar"),
                      (self.cst_disc_rate_lineEdit, "conc_solar"), (self.ht_disc_rate_lineEdit, "hydrogen"),
                      (self.gt_disc_rate_lineEdit, "geo_therm")]:
            try:
                val = float(entry[0].text())    
                self.assump_json[entry[1]]['third_party_rates']['discount_rate'] = val                 
            except:
                color_red.append(entry)
                valid = False
                errors.append(f"Third party owned discounted rate must be numeric")
        
        for entry in [(self.wt_infl_rate_lineEdit, "wind"), (self.st_infl_rate_lineEdit, "solar"),
                      (self.cst_infl_rate_lineEdit, "conc_solar"), (self.ht_infl_rate_lineEdit, "hydrogen"),
                      (self.gt_infl_rate_lineEdit, "geo_therm")]:
            try:
                val = float(entry[0].text())    
                self.assump_json[entry[1]]['third_party_rates']['inflation_rate'] = val                
            except:
                color_red.append(entry)
                valid = False
                errors.append(f"Third party owned inflation rate must be numeric")
        
        for entry in [(self.wt_ece_LB_lineEdit, "wind"), (self.st_ece_LB_lineEdit, "solar"),
                      (self.cst_ece_LB_lineEdit, "conc_solar"), (self.ht_ece_LB_lineEdit, "hydrogen"),
                      (self.gt_ece_LB_lineEdit, "geo_therm")]:
            try:
                val = float(entry[0].text())    
                self.assump_json[entry[1]]['third_party_rates']['electric_cost_escalation_lb'] = val                
            except:
                color_red.append(entry)
                valid = False
                errors.append(f"Third party owned electric cost escalation lower bound must be numeric")
        
        for entry in [(self.wt_ece_UB_lineEdit, "wind"), (self.st_ece_UB_lineEdit, "solar"),
                      (self.cst_ece_UB_lineEdit, "conc_solar"), (self.ht_ece_UB_lineEdit, "hydrogen"),
                      (self.gt_ece_UB_lineEdit, "geo_therm")]:
            try:
                val = float(entry[0].text())    
                self.assump_json[entry[1]]['third_party_rates']['electric_cost_escalation_ub'] = val                  
            except:
                color_red.append(entry)
                valid = False
                errors.append(f"Third party owned electric cost escalation upper bound must be numeric")
        
        for entry in [[self.wt_ece_LB_lineEdit,self.wt_ece_UB_lineEdit], [self.st_ece_LB_lineEdit,self.st_ece_UB_lineEdit],
                      [self.cst_ece_LB_lineEdit,self.cst_ece_UB_lineEdit], [self.ht_ece_LB_lineEdit,self.ht_ece_UB_lineEdit],
                      [self.gt_ece_LB_lineEdit,self.gt_ece_UB_lineEdit]]:       
            try:
                val_lb = float(entry[0].text())
                val_ub = float(entry[1].text())
                if val_ub <= val_lb:
                    color_red.append(entry[0])
                    color_red.append(entry[1])
                    valid = False
                    errors.append(f"The third party owned upper bound on electric cost escalation must be larger than the lower bound") 
            except:
                pass
        
        for entry in [(self.wt_tax_lineEdit, "wind"), (self.st_tax_lineEdit, "solar"),
                      (self.cst_tax_lineEdit, "conc_solar"), (self.ht_tax_lineEdit, "hydrogen"),
                      (self.gt_tax_lineEdit, "geo_therm")]:
            try:
                val = float(entry[0].text())    
                self.assump_json[entry[1]]['third_party_rates']['taxation'] = val                 
            except:
                color_red.append(entry)
                valid = False
                errors.append(f"Third party owned taxation rate must be numeric")
        
        for entry in [(self.wt_proj_life_lineEdit, "wind"), (self.st_proj_life_lineEdit, "solar"),
                      (self.cst_proj_life_lineEdit, "conc_solar"), (self.ht_proj_life_lineEdit, "hydrogen"),
                      (self.gt_proj_life_lineEdit, "geo_therm")]:
            try:
                val = float(entry[0].text())    
                self.assump_json[entry[1]]['third_party_rates']['project_life'] = val  
                if val <= 0:
                    color_red.append(entry)
                    valid = False
                    errors.append(f"Third party owned project life must be greater than 0") 
            except:
                color_red.append(entry)
                valid = False
                errors.append(f"Third party owned project life must be numeric")
        
        errors.append("")
        errors.append("*If any values were found to be non-numeric, they were reset to their default values")
        self.valid = valid
        self.errors = errors
        self.color_red = []
        for qq, obj in enumerate(self.line_edits):
            if obj in color_red:
                self.color_red.append(qq)