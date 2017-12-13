function Fp = dengue_twostrains(~, F)
%{
    This function computes the rhds of our two strains model.
    zp reprsents the reported cases and wq the new hemorragic dengue
    cases.
%}

global Lambda_M beta_M beta_H b mu_M Lambda_S mu_H alpha_c alpha_h ...
        Lambda_m1 Lambda_m2 Sigma p q

Ms = F(1); M1 = F(2); M2 = F(3); Is = F(4);
I1 = F(5); I2 = F(6); Sm1 = F(7); Sm2 = F(8);
Ym1 = F(9); Ym2 = F(10); R = F(11);

N_H = Is + I1 + I2 + Sm1 + Sm2 + Ym1 + Ym2 + R;
cM= (beta_M * b / N_H);
% Infection forces
A1 = cM * I1; A2 = cM * I2; A3 = cM * Ym1; A4 = cM * Ym2;
%
%
% rhs funcuntions
%
Msp = Lambda_M - (A1 + A2 + A3 + A4) * Ms - mu_M * Ms;
M1p = (A1 + A4) * Ms - mu_M * M1;
M2p = (A2 + A3) * Ms - mu_M * M2;

cH = (beta_H * b / N_H);
B1 = cH * M1; B2 = cH * M2;

Isp = Lambda_S - (B1 + B2) * Is - mu_H * Is;
I1p = B1 * Is - (alpha_c + mu_H) * I1;
I2p = B2 * Is - (alpha_c + mu_H) * I2;
Sm1p = Lambda_m1 - Sigma * B2 * Sm1 - mu_H * Sm1;
Sm2p = Lambda_m2 - B1 * Sm2 - mu_H * Sm2;
Ym1p= Sigma * B2 * Sm1 - (alpha_h + mu_H) * Ym1;
Ym2p = B1 * Sm2 -(alpha_h + mu_H) * Ym2;
Rp = alpha_c * (I1+I2) + alpha_h * (Ym1 + Ym2) - mu_H * R;
%
% Data fit equations
%
zp = p * (I1p + I2p + (1-q) * Ym1p + Ym2p); 
wp = p * q * Ym1p;


Fp=[Msp; M1p; M2p; Isp; I1p; I2p; Sm1p; Sm2p; Ym1p; Ym2p; Rp; zp; wp];
