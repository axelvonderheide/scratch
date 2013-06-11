#--------------------------------------------------------
# Definition of the general input parameters used in the 
# calculation of stresses and required number of 
# reinforcement layers
#--------------------------------------------------------

#--------------------------------------------------------
# Geometry
#--------------------------------------------------------
# thickness at border [m]
D_b = 0.040

# thickness at top [m]
D_t = 0.030

#--------------------------------------------------------
# GdG
#--------------------------------------------------------
# tensile strength of the concrete[MPa]
f_ctk = 1.6

# flexural tensile trength of the concrete [MPa]
f_m = 10.5

#--------------------------------------------------------
# GdT
#--------------------------------------------------------

# INDEX l: longitudinal direction of the textile (MAG-02-02-06a)
# tensile strength of the tensile specimen [N/mm2]
# containing a gamma-factor of 1.5 and d long term reduction factor of 0.7
f_td_l = 251

# cross sectional area of the reinforcement [mm2/m]
a_t_l = 71.65

# INDEX q: orthogonal direction of the textile (MAG-02-02-06a)
# tensile strength of the tensile specimen [kN/m]
f_td_q = 238

# cross sectional area of the reinforcement [mm2/m]
a_t_q = 53.31

# tensile strength of the textile reinforcement [kN/m]
F_Rtex_l = a_t_l * f_td_l / 1000.  

# tensile strength of the textile reinforcement [kN/m]
F_Rtex_q = a_t_q * f_td_q / 1000.

# print defined values
print '---'
print 'Defined input params in file INPUTPARAMS'
print '---'
print ' f_ctk    [N/mm2] = ', f_ctk
print ' f_m      [N/mm2] = ', f_m
print '---'
print ' F_Rtex_l [kN/m]  = ', F_Rtex_l
print ' F_Rtex_q [kN/m]  = ', F_Rtex_q
print '---\n'

#        print ' f_ctk    [N/mm2] = ', self.selected_data_file.f_ctk
#        print ' f_m      [N/mm2] = ', self.selected_data_file.f_m
#        print ' F_Rtex_l [kN/m]  = ', self.selected_data_file.F_Rtex_l
#        print ' F_Rtex_q [kN/m]  = ', self.selected_data_file.F_Rtex_q

