function R0 = reproductiveBasicNumber() 
global Lambda_M beta_M beta_H b mu_M Lambda_S mu_H alpha_c alpha_h ...
        Lambda_m1 Lambda_m2 

N_H = (Lambda_S + Lambda_m1 + Lambda_m2) / mu_H;
N_M = Lambda_M / mu_M;
Pi_r =(beta_H * beta_M * b^2 * N_M) / (mu_M * mu_H * N_H ^ 2);
R01sq = Pi_r * (Lambda_m1 / (alpha_h + mu_H) ...
                + Lambda_S / (alpha_c + mu_H));
R02sq = Pi_r * (Lambda_m2 / (alpha_h + mu_H) ...
                + Lambda_S / (alpha_c  +mu_H));

R01 = sqrt(R01sq);
R02 = sqrt(R02sq);
R0 = [R01, R02];