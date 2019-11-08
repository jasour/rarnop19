function PM=SubChebX(P,moxv)

global x1
p1=P - x1^1000^1000;
Dpc=p1.deg;  Dc=0+Dpc(:,:); Dc(1,:)=[];% Bc1: pwers
Cpc=p1.coef; Cc=0+Cpc(:,:); Cc(1,:)=[];% Cc1 :coefficent, 

if isempty(Cc)==1
    
    PM=0;
else
    
if Dc==0
   PM = Cc; 
else
    
pm=0;PM=0;
z1=chebfun('z1');

for i=1:size(Dc,1)

pp1=Cc(i,1)*z1.^Dc(i,1); pc1=rot90(chebpoly(pp1),-1); 
spc1=sparse(pc1);
[I,J,V]=find(spc1); ind=[[J-1]',V'];

pm=0;
for j=1:size(ind,1)
pm=pm+ind(j,2)*moxv(glex2num(ind(j,1)));
end
PM=PM + pm;% Output

end
end
end