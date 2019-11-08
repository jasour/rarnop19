%% SEMIDEFINITE PROGRAMMING FOR CHANCE CONSTRAINED OPTIMIZATION OVER SEMIALGEBRAIC SETS
% SIAM J. OPTIM. Vol. 25, No. 3, pp. 1411–1440
% A. M. JASOUR, N. S. AYBAT, AND C. M. LAGOA
%% Section 3.5. Orthogonal basis.
% The optimal values of the SDP in (3.7) formed by Chebyshev polynomials of the first kind

clc;clear all
addpath ./chebfun
global x1 x2 Mind_y

d=4;
tic
%% Genreate Powers
disp('Generate Powers')
vpowd=[]; for k = 0:d; vpowd = [vpowd;genpow(2,k)]; end
vpow2d=[];for k = 0:2*d; vpow2d = [vpow2d;genpow(2,k)]; end
vpowdx=[]; for k = 0:d; vpowdx = [vpowdx;genpow(1,k)]; end
vpow2dx=[];for k = 0:2*d; vpow2dx = [vpow2dx;genpow(1,k)]; end

%% Vector of Orthogonal Moments mov, mox 
mov=[];moxv=[];
for i=size(vpow2d,1):-1:1; 
eval(['syms mo',num2str(vpow2d(i,1)),'mo',num2str(vpow2d(i,2)),' real']);eval(['mov=[mo',num2str(vpow2d(i,1)),'mo',num2str(vpow2d(i,2)),';mov];']);
end
for i=size(vpow2dx,1):-1:1; 
    eval(['syms mox',num2str(vpow2dx(i,1)),' real']); eval(['moxv=[mox',num2str(vpow2dx(i,1)),';moxv];']);
end
moxv=[moxv];

%% Monomial Basis in Terms of Chebyshev Basis : x^ix^j = Tij
disp('Generate Monomial Basis in Terms of Chebyshev Basis')
pvar x1 x2
T=[]; for i=1:size(vpow2d);clc;
    disp({'T',i,size(vpow2d,1)});T=[T;SubCheb(x1^vpow2d(i,1)*x2^vpow2d(i,2),mov) ];
end
Tx=[];
for i=1:size(vpow2dx)
    clc;
    disp({'Tx',i,size(vpow2dx,1)});
    Tx=[Tx;SubChebX(x1^vpow2dx(i,1),moxv) ]; 
end
Tx(1)=mox0;
%% C: Univariant Chebyshev Basis  T_i = C* x^i
disp('Generate Univariant Chebyshev Basis')
N=2*d+1;C=zeros(N,N);
for i=0:N-1; cc=flip(ChebyshevPoly(i))'; C(i+1,1:size(cc,2))=cc; end

%% P: 2D Chebyshev Basis T_ij = C* x^ix^j
disp('Generate 2D Chebyshev Basis')

for k=size(vpow2d,1):-1:1
clc;disp({'P',k})
k1=vpow2d(k,1)+1;k2=vpow2d(k,2)+1;
cs1=sparse(C(k1,:)); cs2=sparse(C(k2,:));CS=cs1'*cs2;[I,J,V]=find(CS); 
if size(I,1)==1; I=I';J=J';V=V';end
ind=[I-1,J-1,V];
for j=1:size(ind,1);P(k,glex2num(ind(j,1:2)))=V(j); end
end
%%
clc
%% Moment Matrix x1-x2
disp('Generate Moment Matrix Mxq')
[Mind_y,Mind_x]=MomentIndex(2,1,d);
Mm1=T(Mind_y);
Mxq=P(1:size(Mm1,1),1:size(Mm1,1))*Mm1*P(1:size(Mm1,1),1:size(Mm1,1))';

AA1=[];Po1=1000*sum(mov);
for i=1:size(Mxq,1)
    for j=1:i;     clc;disp({'Linear Operator x-q',i,j,size(Mxq,1)});
     Po=Mxq(i,j)+Po1;
     [Co1,Do2]=coeffs(Po);
     AA1{i,j}=double(Co1)-1000;
     AA1{j,i}= AA1{i,j};
    end
end
AA1=cell2mat(AA1(:));
%% Localization Matrix

disp('Generate Localization Matrix LMxq')
pc1= 0.5*x2*(x2^2+(x1-0.5)^2)-(x2^4+x2^2*(x1-0.5)^2+(x1-0.5)^4);%0.5*((0.5)^2 - (x1-0.5)^2 - (z1-0)^2); 
Cpc1=pc1.coef; Cc1=0+Cpc1(:,:); % Cc1 :coefficent, 
Dpc1=pc1.deg;  Dc1=0+Dpc1(:,:); % Bc1: pwers
d1=floor(abs(2*d- max(sum(Dc1,2)) )/2);
LMm=LMoment(Dc1,Cc1,T,d1);
LMxq=P(1:size(LMm,1),1:size(LMm,1))*LMm*P(1:size(LMm,1),1:size(LMm,1))';

AA2=[];Po1=1000*sum(mov);
for i=1:size(LMxq,1);
    for j=1:i;clc;disp({'Linear Operator Localization',i,j,size(LMxq,1)});
     Po=LMxq(i,j)+Po1;
     [Co1,Do2]=coeffs(Po);
     %AA2=[AA2;double(Co1)-1000];     
     AA2{i,j}=double(Co1)-1000;
     AA2{j,i}= AA2{i,j};
    end
end
AA2=cell2mat(AA2(:));

%% Moment Matrix x1
disp('Generate Moment Matrix Mx')
Mmx1=Tx(Mind_x);
Mx=C(1:d+1,1:d+1)*Mmx1*C(1:d+1,1:d+1)';

AA3=[];Po1=1000*sum(moxv);
for i=1:size(Mx,1)
    for j=1:i;clc;disp({'Linear Operator x',i,j,size(Mx,1)});
     Po=Mx(i,j)+Po1;
     [Co1,Do2]=coeffs(Po);  
         AA3{i,j}=double(Co1)-1000;
     AA3{j,i}= AA3{i,j};
    end
end
AA3=cell2mat(AA3(:));
%% Upper Moment Matrix
disp('Generate Upper Moment Matrix MUxq')
[Yx,Yq]=MomentIndex_xq(2,1,2*d);
MAx = MomentCross(Tx',Yx,Yq);
Mmup=MAx; %
MUxq= P(1:size(Mmup,1),1:size(Mmup,1))*Mmup*P(1:size(Mmup,1),1:size(Mmup,1))';

AA4=[];Po1=1000*sum(moxv);
for i=1:size(MUxq,1)
    for j=1:i;clc;disp({'Linear Operator Upper',i,j,size(MUxq,1)});
     Po=MUxq(i,j)+Po1;
     [Co1,Do2]=coeffs(Po); 
     AA4{i,j}=double(Co1)-1000;
     AA4{j,i}= AA4{i,j};
    end
end
AA4=cell2mat(AA4(:));

%% Sedumi Matrixes

A1=[AA1,zeros(size(AA1,1),size(moxv,1)-1)];
A2=[AA2,zeros(size(AA2,1),size(moxv,1)-1)];
A3=[zeros(size(AA3,1),size(mov,1)),AA3(:,2:end)];
A4=[zeros(size(AA4,1),size(mov,1)),AA4(:,2:end)]-A1;

A=[-A1;-A2;-A3;-A4];

b=[1;zeros(size(A,2)-1,1)];

c=[zeros(size(A1,1),1);zeros(size(A2,1),1);AA3(:,1);AA4(:,1)];


%%
K.f=[];
K.s=[size(Mxq,1),size(LMxq,1),size(Mx,1),size(MUxq,1)];
pars = [];
pars.fid = 0;

clc;disp('start')
[x,y,info] = sedumi(A,b,c,K,pars);

Pro=y(1)
Dec= y(size(AA1,2)+1)

toc
