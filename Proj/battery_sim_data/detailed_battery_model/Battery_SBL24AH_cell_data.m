function data = Battery_SBL24AH_cell_data()
% A123:AMP20M1HD datasheet

AH      = 24;
SOC_vec = 0:10:100;
OCV_vec = fliplr([4.1798 4.0718 3.9978 3.926 3.857 3.7546 3.6832 3.6418 3.6063 3.5299 3.4397]);
SOC_R0_vec = 10:10:90;
R0_vec  = [1.70 1.52 1.47 1.48 1.56 1.46 1.50 1.53 1.48] * 1e-3;
T_vec   = 25 + 273.15;

data.SOC_vec = SOC_vec;
data.OCV_vec = OCV_vec;
data.T_vec   = T_vec;
data.R0_vec  = R0_vec;
data.AH      = AH;
data.SOC_R0_vec = SOC_R0_vec;
end