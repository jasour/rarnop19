
Prob_List=[];
for u=-1:0.01:1

% number of uncertainty samples
N=1000000;

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
a = 0.1;b = 0.1;c=0.4;
K=1^2-(f1-0).^2/a^2-(f2-0).^2/b^2-(f3-0.9).^2/c^2;


Pro=size(find(K>=0),2)/N;
Prob_List=[Prob_List;u,Pro];

end
plot([-1:0.01:1],Prob_List(:,2),'--','LineWidth',3)

 ylim([0 1])