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
sdpvar mm1 ss1
py=[
   1
    mm1,
    mm1^2+ss1^2,
    mm1^3+3*mm1*ss1^2,
    mm1^4+6*mm1^2*ss1^2+3*ss1^4,
    mm1^5+10*mm1^3*ss1^2+15*mm1*ss1^4,
    mm1^6+15*mm1^4*ss1^2+45*mm1^2*ss1^4+15*ss1^6,
    mm1^7+21*mm1^5*ss1^2+105*mm1^3*ss1^4+105*mm1*ss1^6,
    mm1^8+28*mm1^6*ss1^2+210*mm1^4*ss1^4+420*mm1^2*ss1^6+105*ss1^8,
    mm1*(mm1^8 + 36*mm1^6*ss1^2 + 378*mm1^4*ss1^4 + 1260*mm1^2*ss1^6 + 945*ss1^8),
    mm1^10 + 45*mm1^8*ss1^2 + 630*mm1^6*ss1^4 + 3150*mm1^4*ss1^6 + 4725*mm1^2*ss1^8 + 945*ss1^10,
    mm1*(mm1^10 + 55*mm1^8*ss1^2 + 990*mm1^6*ss1^4 + 6930*mm1^4*ss1^6 + 17325*mm1^2*ss1^8 + 10395*ss1^10),
    mm1^12 + 66*mm1^10*ss1^2 + 1485*mm1^8*ss1^4 + 13860*mm1^6*ss1^6 + 51975*mm1^4*ss1^8 + 62370*mm1^2*ss1^10 + 10395*ss1^12,
    mm1*(mm1^12 + 78*mm1^10*ss1^2 + 2145*mm1^8*ss1^4 + 25740*mm1^6*ss1^6 + 135135*mm1^4*ss1^8 + 270270*mm1^2*ss1^10 + 135135*ss1^12),
    mm1^14 + 91*mm1^12*ss1^2 + 3003*mm1^10*ss1^4 + 45045*mm1^8*ss1^6 + 315315*mm1^6*ss1^8 + 945945*mm1^4*ss1^10 + 945945*mm1^2*ss1^12 + 135135*ss1^14,
    mm1*(mm1^14 + 105*mm1^12*ss1^2 + 4095*mm1^10*ss1^4 + 75075*mm1^8*ss1^6 + 675675*mm1^6*ss1^8 + 2837835*mm1^4*ss1^10 + 4729725*mm1^2*ss1^12 + 2027025*ss1^14),
    mm1^16 + 120*mm1^14*ss1^2 + 5460*mm1^12*ss1^4 + 120120*mm1^10*ss1^6 + 1351350*mm1^8*ss1^8 + 7567560*mm1^6*ss1^10 + 18918900*mm1^4*ss1^12 + 16216200*mm1^2*ss1^14 + 2027025*ss1^16];