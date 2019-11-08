function PM=SubCheb(P,mov)
global x1 x2
p1=P - x1^1000^1000*x2^1000^1000;%
Cpc=p1.coef; Cc=0+Cpc(:,:); Cc(1,:)=[];% Cc1 :coefficent, 
Dpc=p1.deg;  Dc=0+Dpc(:,:); Dc(1,:)=[];% Bc1: pwers

P2H=[];PmH=[];
pm=0;PM=0;
z1=chebfun('z1');z2=chebfun('z2');

for i=1:size(Dc,1)
    
%P2H=[P2H;Cc(i)*x1^Dc(i,1)*x2^Dc(i,2)];

pp1=Cc(i,1)*z1.^Dc(i,1); pc1=rot90(chebpoly(pp1),-1); 
spc1=sparse(pc1);
pp2=z2.^Dc(i,2); pc2=rot90(chebpoly(pp2),-1); 
spc2=sparse(pc2);

sT=spc1'*spc2; [I,J,V]=find(sT); 
if size(I,1)==1; I=I';J=J';V=V';end
ind=[I-1,J-1,V];

pm=0;
for j=1:size(ind,1)
    co=ind(j,3);% if co<10^(-10); co=0;end
pm=pm+ co*mov(glex2num(ind(j,1:2)));
end
%PmH=[PmH;pm];
PM=PM+pm;% Output
end
