% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 8: Nonlinear Robust Optimization
%% q has normal distribution with uncertain mean and sigma
% mean in [-0.1, 0.1]; sigma in [0.1 0.3]
% moments of q are polynomial function of mean and sigma
%Polynomial moments of uncertain normal distribution;
% m: mean in A1=[-0.1, 0.1],  % ss: sigma in A2=[0.1, 0.3]
%% Polynomial moment of order a : int q^a pr(q)dq : parametric integral
% List of Polynomial moments up to order 16
% Provided by Tillmann Weisser
%Jean B. Lasserre, Tillmann Weisser,”Distributionally robust polynomial chance-constraints under mixture ambiguity sets”, Mathematical Programming, pp 1–45, 2019
%%
clc
sdpvar mm3 ss3
py=[
   1
    mm3,
    mm3^2+ss3^2,
    mm3^3+3*mm3*ss3^2,
    mm3^4+6*mm3^2*ss3^2+3*ss3^4,
    mm3^5+10*mm3^3*ss3^2+15*mm3*ss3^4,
    mm3^6+15*mm3^4*ss3^2+45*mm3^2*ss3^4+15*ss3^6,
    mm3^7+21*mm3^5*ss3^2+105*mm3^3*ss3^4+105*mm3*ss3^6,
    mm3^8+28*mm3^6*ss3^2+210*mm3^4*ss3^4+420*mm3^2*ss3^6+105*ss3^8,
    mm3*(mm3^8 + 36*mm3^6*ss3^2 + 378*mm3^4*ss3^4 + 1260*mm3^2*ss3^6 + 945*ss3^8),
    mm3^10 + 45*mm3^8*ss3^2 + 630*mm3^6*ss3^4 + 3150*mm3^4*ss3^6 + 4725*mm3^2*ss3^8 + 945*ss3^10,
    mm3*(mm3^10 + 55*mm3^8*ss3^2 + 990*mm3^6*ss3^4 + 6930*mm3^4*ss3^6 + 17325*mm3^2*ss3^8 + 10395*ss3^10),
    mm3^12 + 66*mm3^10*ss3^2 + 1485*mm3^8*ss3^4 + 13860*mm3^6*ss3^6 + 51975*mm3^4*ss3^8 + 62370*mm3^2*ss3^10 + 10395*ss3^12,
    mm3*(mm3^12 + 78*mm3^10*ss3^2 + 2145*mm3^8*ss3^4 + 25740*mm3^6*ss3^6 + 135135*mm3^4*ss3^8 + 270270*mm3^2*ss3^10 + 135135*ss3^12),
    mm3^14 + 91*mm3^12*ss3^2 + 3003*mm3^10*ss3^4 + 45045*mm3^8*ss3^6 + 315315*mm3^6*ss3^8 + 945945*mm3^4*ss3^10 + 945945*mm3^2*ss3^12 + 135135*ss3^14,
    mm3*(mm3^14 + 105*mm3^12*ss3^2 + 4095*mm3^10*ss3^4 + 75075*mm3^8*ss3^6 + 675675*mm3^6*ss3^8 + 2837835*mm3^4*ss3^10 + 4729725*mm3^2*ss3^12 + 2027025*ss3^14),
    mm3^16 + 120*mm3^14*ss3^2 + 5460*mm3^12*ss3^4 + 120120*mm3^10*ss3^6 + 1351350*mm3^8*ss3^8 + 7567560*mm3^6*ss3^10 + 18918900*mm3^4*ss3^12 + 16216200*mm3^2*ss3^14 + 2027025*ss3^16];