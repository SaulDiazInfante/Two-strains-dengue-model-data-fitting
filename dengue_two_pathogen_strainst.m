%{ 
Scritp to solve our two strain model.
Here we use ode45 to integrate the proposed model
by Diaz-Infante and Olmos.
OUTPUT from ode45 is of the form
[t, X]
where
X = [x1(t) x2(t)]
%}
%
%
% Global parameters.
%
%
clear
clc
global Lambda_M beta_M beta_H b mu_M Lambda_S mu_H alpha_c alpha_h ...
        Lambda_m1 Lambda_m2 p q

% Mosquitoes constants
Lambda_M = 20000;
beta_M = 0.011;
b = 3.8;
mu_M = 0.04;

%Human constants
Lambda_S = 1.4;
beta_H = 0.01;
mu_H = 0.00003653;  % 75 years.
alpha_c = 0.1;      % 10 days for recovering
alpha_h = 0.5;      % 2 days for recovering   
%
%
Lambda_m1 = .015;
Lambda_m2 = 0.065;
%
% Data fit parameters
p = 0.05;
q = 0.9;

% Numerical constants
h=0.1; T = 250;
%
% Time interval
tspan = 0 : h :T;
%
% Set options for the ODE solver ode45
% 
% Initial condition mosquitoes
Ms0 = 120000; M10 = 20; M20 = 30;

% Initial condition humans
Is0 = 35600; I10 = 1; I20=20; 
Sm10 = 4400; Sm20 = 0;
Ym10 = 0; Ym20 = 0; 
Rec0 = 0;

z0 = p * (I10 + I20 + (1-q) * Ym10 + Ym20);
w0 = p * q * Ym10;
    
X0 = [Ms0; M10; M20; Is0; I10; I20; Sm10; Sm20; Ym10; Ym20; Rec0; z0; w0];
par = [Lambda_M; Lambda_S; ...
       Lambda_m1; Lambda_m2; ...
       beta_M; beta_H; b; ...
       mu_M; mu_H; ...
       alpha_c ; alpha_h; ...   
       p; q; ...
       Ms0; M10; M20; ... 
       Is0; I10; I20; ...
       Sm10; Sm20; Ym10; ...
       Ym20; Rec0; z0; w0; ...
       h; T];
saveParameters(par);

% Now solve the first order system
options = odeset('RelTol', 0.0000001, 'AbsTol', 0.000000001);
[t, X] = ode45('dengue_twostrains', tspan, X0,  options);
% Solutions in time
z = X(:,12); w = X(:,13);

R0 = reproductiveBasicNumber();
text = ['R0: max{', num2str(R0(1)),', ' num2str(R0(2)),'}'];
disp(text);
popultationsPlot(t, X);
%%%%
figure(2)
plot(t, z,'b',t, w,'r')
legend('CD','HDF')
xlabel('t')

figure(3)
hold off

% Load data sets of DF and HDF for 2010 in Hermosillo, Sonora

datFD=load('dengue_c_her2010.dat');
diaFD=datFD(:,1)-40386;
casosFD=datFD(:,2);

plot(diaFD, casosFD,'b*');
grid on
hold on

datFHD=load('dengue_h_her2010.dat');
diaFHD=datFHD(:,1)-40384;
casosFHD=datFHD(:,2);

plot(diaFHD,casosFHD,'ro')
hold on
plot(t, z,'b',t, w, 'r')
legend('CD','HDF')
xlabel('t')
axis([0 160 0 70])

title('Confirmed DF and DHF cases, Hermosillo 2010')
xlabel('date')
ylabel('Incidence')

legend('FD','FHD')
