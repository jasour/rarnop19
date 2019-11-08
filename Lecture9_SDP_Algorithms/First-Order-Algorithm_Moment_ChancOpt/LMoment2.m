function [Lco,Lmo,AInd_L]=LMoment2(Dc,Cc,x,d1)
global Ny
n2=size(Dc,2);
vpow=[];
for k = 0:d1
    vpow = [vpow;genpow(n2,k)];
end
size(vpow,1)
for i=1:size(vpow,1)   
    for j=1:i
        clc;disp('Localization Matrix');disp([i,j,size(vpow,1)])
        a=vpow(i,:)+vpow(j,:);
%         b=0;
        for k=1:size(Dc)
%            b=b+Cc(k)*x(glex2num(a+Dc(k,:)));
            Lco(i,j,k)= Cc(k);
            Lmo(i,j,k)= glex2num(a+Dc(k,:));
       
           Lco(j,i,k)=Lco(i,j,k);         
          Lmo(j,i,k)=Lmo(i,j,k);
        end
          
          %LM(i,j)=b;
%         LM(j,i)=LM(i,j);
    end
end

AInd_L=zeros(Ny,1);
for i=1:Ny
    clc;disp('Adjoint Localization Matrix');disp([i,Ny])
[ind]=find(Lmo==i);
AInd_L(i,1:size(ind,1))=ind';
end