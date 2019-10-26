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
sdpvar mm ss
py=[
   1
    mm,
    mm^2+ss^2,
    mm^3+3*mm*ss^2,
    mm^4+6*mm^2*ss^2+3*ss^4,
    mm^5+10*mm^3*ss^2+15*mm*ss^4,
    mm^6+15*mm^4*ss^2+45*mm^2*ss^4+15*ss^6,
    mm^7+21*mm^5*ss^2+105*mm^3*ss^4+105*mm*ss^6,
    mm^8+28*mm^6*ss^2+210*mm^4*ss^4+420*mm^2*ss^6+105*ss^8,
    mm*(mm^8 + 36*mm^6*ss^2 + 378*mm^4*ss^4 + 1260*mm^2*ss^6 + 945*ss^8),
    mm^10 + 45*mm^8*ss^2 + 630*mm^6*ss^4 + 3150*mm^4*ss^6 + 4725*mm^2*ss^8 + 945*ss^10,
    mm*(mm^10 + 55*mm^8*ss^2 + 990*mm^6*ss^4 + 6930*mm^4*ss^6 + 17325*mm^2*ss^8 + 10395*ss^10),
    mm^12 + 66*mm^10*ss^2 + 1485*mm^8*ss^4 + 13860*mm^6*ss^6 + 51975*mm^4*ss^8 + 62370*mm^2*ss^10 + 10395*ss^12,
    mm*(mm^12 + 78*mm^10*ss^2 + 2145*mm^8*ss^4 + 25740*mm^6*ss^6 + 135135*mm^4*ss^8 + 270270*mm^2*ss^10 + 135135*ss^12),
    mm^14 + 91*mm^12*ss^2 + 3003*mm^10*ss^4 + 45045*mm^8*ss^6 + 315315*mm^6*ss^8 + 945945*mm^4*ss^10 + 945945*mm^2*ss^12 + 135135*ss^14,
    mm*(mm^14 + 105*mm^12*ss^2 + 4095*mm^10*ss^4 + 75075*mm^8*ss^6 + 675675*mm^6*ss^8 + 2837835*mm^4*ss^10 + 4729725*mm^2*ss^12 + 2027025*ss^14),
    mm^16 + 120*mm^14*ss^2 + 5460*mm^12*ss^4 + 120120*mm^10*ss^6 + 1351350*mm^8*ss^8 + 7567560*mm^6*ss^10 + 18918900*mm^4*ss^12 + 16216200*mm^2*ss^14 + 2027025*ss^16];