function parameters_values = loadParameters(file_name)
% Read all lines & collect in cell array
fid = fopen(file_name);
txt = textscan(fid,'%s','delimiter','\n'); 
fclose(fid);
% Convert string to numerical value
global Lambda_M Lambda_S  Lambda_m1 Lambda_m2 ...
       beta_M beta_H b ...
       mu_M mu_H ...
       alpha_c alpha_h ...   
       p q ...
       Ms0 M10 M20 ... 
       Is0 I10 I20 ...
       Sm10 Sm20 Ym10 ...
       Ym20 Rec0 z0 w0 ...
       h T

par = str2double(txt{1}); 
Lambda_M = par(1); 
Lambda_S = par(2);
Lambda_m1 = par(3); 
Lambda_m2 = par(4);
beta_M = par(5);
beta_H = par(6); 
b = par(7);
mu_M = par(8);
mu_H = par(9);
alpha_c = par(10); 
alpha_h = par(11);
p = par(12);
q = par(13);
Ms0 = par(14);
M10 = par(15);
M20 = par(16);
Is0 = par(17);
I10 = par(18);
I20 = par(19);
Sm10 = par(20);
Sm20 = par(21);
Ym10 = par(22);
Ym20 = par(23);
Rec0 = par(24);
z0 = par(25);
w0 = par(26);
h = par(27); 
T = par(28);
%
parameters_values = true;