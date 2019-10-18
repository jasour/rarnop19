

% number of uncertainty samples
N=1000000;
% Obtained Control input
u=Decision_x;

% uncertainties
x1=random('Uniform',-0.1,0.1,1,N);
x2=random('Uniform',-0.1,0.1,1,N);
x3=random('Uniform',-0.1,0.1,1,N);
w=random('Beta',2,5,1,N);

% x(k+1)
f1=0.2*w.*x2;
f2=x1.*x3;
f3=1.2*x1-0.5*x2+x3+2*u;

% success set
K=1^2-(f1-0).^2/0.03^2-(f2-0).^2/0.02^2-(f3-0.9).^2/0.4^2;

% estimated probability
Pro=size(find(K>=0),2)/N

