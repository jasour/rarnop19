% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 9: Algorithms for Large Scale Semidefinite Programs (SDPs)
%% SEMIDEFINITE PROGRAMMING FOR CHANCE CONSTRAINED OPTIMIZATION OVER SEMIALGEBRAIC SETS
% SIAM J. OPTIM. Vol. 25, No. 3, pp. 1411–1440
% A. M. JASOUR, N. S. AYBAT, AND C. M. LAGOA
% Section 5.2.  First-order augmented Lagrangian algorithm.
% Simple Example in Section 3.4
% The optimal values of the SDP in (3.7)
%% Simple Exampl 3.4 in the paper

clc;clear;close all;

addpath ./PROPACK_SVT
rng(13,'twister')
global Mind_y Mind_x Mind_xq AInd_y AInd_x Ny Nx Dc1 Cc1 d1 Yx Yq nx Lmo Lco AInd_L


%% Problem Setting
nvar=2; % Number of Decision and Uncertain varibales
nx=1; % Number of Decision Variables
d=10; % Relaxation Order
Ny= round(factorial(nvar+2*d)/(factorial(2*d)*factorial(nvar))); % Number of y variables : Moments of Measure Mu  
Nx= round(factorial(nx+2*d)/(factorial(2*d)*factorial(nx))); % Number of x variables : Moments of Measure Mu_x  

y=rand(Ny,1); y(1)=1; % Moments of measure Mu
x=rand(Nx,1); x(1)=1; % Moments of measure Mu_x

%% Semialgebraic Set

pvar x1 z1
pc1=1*(0.5*z1*(z1^2+(x1-0.5)^2)-(z1^4+z1^2*(x1-0.5)^2+(x1-0.5)^4)) ;


%%
Cpc1=pc1.coef; Cc1=0+Cpc1(:,:); % Cc1 : coefficent, 
Dpc1=pc1.deg;  Dc1=0+Dpc1(:,:); % Bc1 : degrees

d1=floor(abs(2*d- max(sum(Dc1,2)) )/2);


%% Pre-Calculation
[AInd_y,Mind_y,AInd_x,Mind_x]=MomentIndex(nvar,nx,Ny,Nx,d); % Mind: Moment Matrix Indx, Mindx: Moment Matrix Indx
[Yx,Yq,Mind_xq]=MomentIndex_xq(nvar,nx,Ny,Nx,2*d); % Upper Bound
[Lco,Lmo,AInd_L] = LMoment2(Dc1,Cc1,y,d1);
%% Constants

LM=sum(Lco(:,:,:).*y(Lmo(:,:,:)),3);
disp('Lipschitze Constant')
L_My = lansvd('MO','MOt',length(Mind_y(:)),Ny,1,'L')^2;
L_LM = lansvd('LO','LOt',length(LM(:)),Ny,1,'L')^2;
L_Mx = lansvd('MxO','MxOt',length(Mind_x(:)),Nx,1,'L')^2;
L_MAx = lansvd('UPO','UPOt',length(Mind_y(:)),Nx,1,'L')^2;

%% Cost Function
cr =zeros(size(x)); cr(diag(Mind_x))=1; % cr'*x = Trace(Mx) : Rank
cp=zeros(size(y));cp(1)=-1;% cp'y : Probability

%% Algorithm Parameters *******************************************
tol=1e-3;
scale=0.8; 
Ly =  scale*(2*L_My + L_LM);
Lx =  scale*(L_Mx + L_MAx);

eta_y=inf; eta_x=inf;
nrm_GPx = inf; nrm_GPy = inf;
mu=2*(Ly+Lx);
lambdar=1e-5; lambdap=1;
coef_mu=1/0.95; delta=0.9;
coef_eta_y = 0.5; coef_eta_y0 = 0.5;
coef_eta_x=0.5; coef_eta_x0 = 0.5;

%% Dual Variables
Y_My=zeros(size(Mind_y)); 
Y_L=zeros(size(Lmo(:,:,1)));  
Y_Up=zeros(size(Mind_y));
Y_Mx=zeros(size(Mind_x));

%% ALLC
k = 0;total_iter=0;opt_flag = 0;
yp = y; xp = x;
clc; disp('start');
  
while ~opt_flag
k = k+1;

%% APG
y1 = y; y2 = y; y1p = y;
x1 = x; x2 = x; x1p = x;
t=1; flag = 0;iter=0;
    
while ~flag
   iter=iter+1;
 clc
disp('-----------------------------------------------')   
disp(['Moment Number :',num2str(Ny)])
disp(['Moment x Number :',num2str(Nx)])
disp(['Moment Matrix : ',num2str(size(Mind_y,1)),'*',num2str(size(Mind_y,1))])
disp(['Localization Matrix : ',num2str(size(Lmo,1)),'*',num2str(size(Lmo,2))])
disp(['Upper Bound Moment Matrix : ',num2str(size(Mind_y,1)),'*',num2str(size(Mind_y,1))])
disp(['Moment Matrix x : ',num2str(size(Mind_x,1)),'*',num2str(size(Mind_x,1))])
disp(['Ly=',num2str(Ly),', Lx=', num2str(Lx)])
disp('-----------------------------------------------')   

str1=sprintf('Outer Loop : %d ---------- Inner Loop : %d',k,iter);disp(str1);disp(['Total iter: ',num2str(total_iter)])
str2=sprintf('Probability : %d -- %d',y(1) , y2(1));disp(str2);
disp('x_outer   x_inner');disp([x(2:1+nx), x2(2:1+nx)])
   
        
    %% Operators

    M=y2(Mind_y);
    [V,D] = eig(M - Y_My);diagD=diag(D);ind=find(diagD > 0);diagD(ind)=0;

    LM=sum(Lco(:,:,:).*y2(Lmo(:,:,:)),3);
    [Vc,Dc] = eig(LM - Y_L);diagDc=diag(Dc);indc=find(diagDc > 0);
    if isempty(indc)~=1;diagDc(indc)=0;end
    
    Mx=x2(Mind_x);
    [Vx,Dx] = eig(Mx - Y_Mx);diagDx=diag(Dx);indx=find(diagDx > 0);diagDx(indx)=0;
    
    MAx = MomentCross(x2,Yx,Yq);
    Mup=MAx-M;
    [Vup,Dup] = eig(Mup - Y_Up);diagDup=diag(Dup);indup=find(diagDup > 0);diagDup(indup)=0;


    %% Adjoint Operators
    AdjM=AdjMoment(V*diag(diagD)*V',AInd_y); 
    AdjL=AdjLMoment(Vc*diag(diagDc)*Vc',AInd_L); 
    AdjMup=AdjMoment(Vup*diag(diagDup)*Vup',AInd_y);
    AdjMx= AdjMoment(Vx*diag(diagDx)*Vx',AInd_x);
    AdjMAx = AdjMoment_xq(Vup*diag(diagDup)*Vup',Yq,Mind_xq,AInd_y,Nx);


    %% Gradients

    GPy = lambdap*cp/mu + AdjM + AdjL - AdjMup;
    GPx = lambdar*cr/mu + AdjMx + AdjMAx;
    
    y1p = y1; y1= y2 - GPy/Ly;
    nrm_GPy = norm(GPy);
    
    x1p = x1; x1= x2 - GPx/Lx; x1(1)=1;
    nrm_GPx = norm(x2-x1)*Lx;

    if iter==1;
        eta_bar_y = nrm_GPy;
        eta_bar_x = nrm_GPx;
        if k==1;
            eta_y = eta_bar_y*coef_eta_y0;
            eta_x = eta_bar_x*coef_eta_x0;
        else
            eta_y = eta_bar_y*coef_eta_y;
            eta_x = eta_bar_x*coef_eta_x;
        end
    else
        if nrm_GPy <= eta_y && nrm_GPx <= eta_x 
            flag=1;
        end
    end

    tp=t; t=(1+sqrt(1+4*t^2))/2; 
    y2=y1 + (tp-1)*(y1-y1p)/t;
    x2=x1 + (tp-1)*(x1-x1p)/t;

end % End iter

       
   total_iter = total_iter+iter;
   y=y1; x=x1;
   
   
   if norm(y-yp)<tol && norm(x-xp)<tol 
       opt_flag=1; 
   end

  
%%  Y and mu Update

       M=y2(Mind_y);
       [V,D] = eig(M - Y_My);
       diagD=diag(D);
       ind=find(diagD > 0);
       diagD(ind)=0;
       Y_My=-V*diag(diagD/coef_mu)*V'; 
       
       LM=sum(Lco(:,:,:).*y2(Lmo(:,:,:)),3);
       [Vc,Dc] = eig(LM - Y_L);
       diagDc=diag(Dc);
       indc=find(diagDc > 0);
       if isempty(indc) == 0
            diagDc(indc)=0;
       end
       Y_L=-Vc*diag(diagDc/coef_mu)*Vc';
       
       Mx=x2(Mind_x);
       [Vx,Dx] = eig(Mx - Y_Mx);
       diagDx=diag(Dx);
       indx=find(diagDx > 0);
       diagDx(indx)=0;
       Y_Mx=-Vx*diag(diagDx/coef_mu)*Vx';
     
       MAx = MomentCross(x2,Yx,Yq);
       Mup=MAx-M;       
       [Vup,Dup] = eig(Mup - Y_Up);
       diagDup=diag(Dup);
       indup=find(diagDup > 0);
       diagDup(indup)=0;
       Y_Up=-Vup*diag(diagDup/coef_mu)*Vup';
        
  mu=coef_mu*mu; 
  eta_y = (k/(k+1))^(2*(1+delta))*eta_y/(coef_mu)^2 ;
  eta_x = (k/(k+1))^(2*(1+delta))*eta_x/(coef_mu)^2 ;
  yp = y;
  xp = x;

end
%%
disp('----------------');
disp('End')


disp('----------------');
str1=sprintf('Probability : %d ',y(1));disp(str1)
str2=sprintf('Decision Parameter: %d', x(2:1+nx)); disp(str2)
disp('----------------');

