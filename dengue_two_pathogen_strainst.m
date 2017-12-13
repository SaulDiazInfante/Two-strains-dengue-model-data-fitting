%{ 
Scritp to solve our two strain model.
Here we use ode45 to integrate the proposed model
by Diaz-Infante and Olmos.
OUTPUT from ode45 is of the form
[t, X]
where
X = [x1(t) x2(t)]
%}

clear
clc
global Ms0 M10 M20 Is0 I10 I20 Sm10 Sm20 Ym10 Ym20 Rec0 z0 w0 h T 

%% Set up model parameters. 
% Uncoment to load an especific parameters set from a data file.
 
% Load parameters file data
%file_name = 'parameters_11-Dec-201718:19:08.dat';
% par = load_parameters(file_name);
par = newParameters();
saveParameters(par);
%
% Time interval
tspan = 0 : h :T;
%
%% Set options for the ODE solver ode45
% 
X0 = [Ms0; M10; M20; Is0; I10; I20; Sm10; Sm20; Ym10; Ym20; Rec0; z0; w0];
options = odeset('RelTol', 0.0000001, 'AbsTol', 0.000000001);

% Now solve the first order system
[t, X] = ode45('dengue_twostrains', tspan, X0,  options);
%% Basic Reproductive Number $\mathcal{R_0}$
R0 = reproductiveBasicNumber();
text = ['R0: max{', num2str(R0(1)),', ' num2str(R0(2)),'}'];
disp(text);
%% Ploting Results
popultationsPlot(t, X);
