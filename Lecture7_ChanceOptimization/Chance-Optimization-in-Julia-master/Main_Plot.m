clc;
Co = dlmread('Sol.txt');

d=10;n=2;
vpow=[];for k = 2*d:-1:0; vpow = [vpow;genpow(n,k)]; end; 

syms x1 q1 
P=sum(Co.*x1.^vpow(:,1).*q1.^vpow(:,2)); % polynomial indicator fuction

% plot polynomial indicator fuction
[x1,q1] = meshgrid(-0.96:0.01:0.96,-0.96:0.01:0.96); P=eval(P);
close all;surfc(x1,q1,P,'FaceColor','red','EdgeColor','none','FaceAlpha',0.7);
camlight; lighting gouraud; hold on;grid on;set(gca,'fontsize',25)
title('dual SDP solution');

