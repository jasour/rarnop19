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
sdpvar mm2 ss2
py=[
   1
    mm2,
    mm2^2+ss2^2,
    mm2^3+3*mm2*ss2^2,
    mm2^4+6*mm2^2*ss2^2+3*ss2^4,
    mm2^5+10*mm2^3*ss2^2+15*mm2*ss2^4,
    mm2^6+15*mm2^4*ss2^2+45*mm2^2*ss2^4+15*ss2^6,
    mm2^7+21*mm2^5*ss2^2+105*mm2^3*ss2^4+105*mm2*ss2^6,
    mm2^8+28*mm2^6*ss2^2+210*mm2^4*ss2^4+420*mm2^2*ss2^6+105*ss2^8,
    mm2*(mm2^8 + 36*mm2^6*ss2^2 + 378*mm2^4*ss2^4 + 1260*mm2^2*ss2^6 + 945*ss2^8),
    mm2^10 + 45*mm2^8*ss2^2 + 630*mm2^6*ss2^4 + 3150*mm2^4*ss2^6 + 4725*mm2^2*ss2^8 + 945*ss2^10,
    mm2*(mm2^10 + 55*mm2^8*ss2^2 + 990*mm2^6*ss2^4 + 6930*mm2^4*ss2^6 + 17325*mm2^2*ss2^8 + 10395*ss2^10),
    mm2^12 + 66*mm2^10*ss2^2 + 1485*mm2^8*ss2^4 + 13860*mm2^6*ss2^6 + 51975*mm2^4*ss2^8 + 62370*mm2^2*ss2^10 + 10395*ss2^12,
    mm2*(mm2^12 + 78*mm2^10*ss2^2 + 2145*mm2^8*ss2^4 + 25740*mm2^6*ss2^6 + 135135*mm2^4*ss2^8 + 270270*mm2^2*ss2^10 + 135135*ss2^12),
    mm2^14 + 91*mm2^12*ss2^2 + 3003*mm2^10*ss2^4 + 45045*mm2^8*ss2^6 + 315315*mm2^6*ss2^8 + 945945*mm2^4*ss2^10 + 945945*mm2^2*ss2^12 + 135135*ss2^14,
    mm2*(mm2^14 + 105*mm2^12*ss2^2 + 4095*mm2^10*ss2^4 + 75075*mm2^8*ss2^6 + 675675*mm2^6*ss2^8 + 2837835*mm2^4*ss2^10 + 4729725*mm2^2*ss2^12 + 2027025*ss2^14),
    mm2^16 + 120*mm2^14*ss2^2 + 5460*mm2^12*ss2^4 + 120120*mm2^10*ss2^6 + 1351350*mm2^8*ss2^8 + 7567560*mm2^6*ss2^10 + 18918900*mm2^4*ss2^12 + 16216200*mm2^2*ss2^14 + 2027025*ss2^16];