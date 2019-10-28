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
sdpvar mm4 ss4
py=[
   1
    mm4,
    mm4^2+ss4^2,
    mm4^3+3*mm4*ss4^2,
    mm4^4+6*mm4^2*ss4^2+3*ss4^4,
    mm4^5+10*mm4^3*ss4^2+15*mm4*ss4^4,
    mm4^6+15*mm4^4*ss4^2+45*mm4^2*ss4^4+15*ss4^6,
    mm4^7+21*mm4^5*ss4^2+105*mm4^3*ss4^4+105*mm4*ss4^6,
    mm4^8+28*mm4^6*ss4^2+210*mm4^4*ss4^4+420*mm4^2*ss4^6+105*ss4^8,
    mm4*(mm4^8 + 36*mm4^6*ss4^2 + 378*mm4^4*ss4^4 + 1260*mm4^2*ss4^6 + 945*ss4^8),
    mm4^10 + 45*mm4^8*ss4^2 + 630*mm4^6*ss4^4 + 3150*mm4^4*ss4^6 + 4725*mm4^2*ss4^8 + 945*ss4^10,
    mm4*(mm4^10 + 55*mm4^8*ss4^2 + 990*mm4^6*ss4^4 + 6930*mm4^4*ss4^6 + 17325*mm4^2*ss4^8 + 10395*ss4^10),
    mm4^12 + 66*mm4^10*ss4^2 + 1485*mm4^8*ss4^4 + 13860*mm4^6*ss4^6 + 51975*mm4^4*ss4^8 + 62370*mm4^2*ss4^10 + 10395*ss4^12,
    mm4*(mm4^12 + 78*mm4^10*ss4^2 + 2145*mm4^8*ss4^4 + 25740*mm4^6*ss4^6 + 135135*mm4^4*ss4^8 + 270270*mm4^2*ss4^10 + 135135*ss4^12),
    mm4^14 + 91*mm4^12*ss4^2 + 3003*mm4^10*ss4^4 + 45045*mm4^8*ss4^6 + 315315*mm4^6*ss4^8 + 945945*mm4^4*ss4^10 + 945945*mm4^2*ss4^12 + 135135*ss4^14,
    mm4*(mm4^14 + 105*mm4^12*ss4^2 + 4095*mm4^10*ss4^4 + 75075*mm4^8*ss4^6 + 675675*mm4^6*ss4^8 + 2837835*mm4^4*ss4^10 + 4729725*mm4^2*ss4^12 + 2027025*ss4^14),
    mm4^16 + 120*mm4^14*ss4^2 + 5460*mm4^12*ss4^4 + 120120*mm4^10*ss4^6 + 1351350*mm4^8*ss4^8 + 7567560*mm4^6*ss4^10 + 18918900*mm4^4*ss4^12 + 16216200*mm4^2*ss4^14 + 2027025*ss4^16];