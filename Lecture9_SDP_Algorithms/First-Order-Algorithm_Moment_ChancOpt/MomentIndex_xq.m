function [Yx,Yq,indxq]=MomentIndex_xq(nvar,nx,Ny,Nx,d)
clc;disp('Upper Bound Moment Matrix');
vpow=[];
for k = 0:d
    vpow = [vpow;genpow(nvar,k)];
end

%% moment q
yq1=[1];for i=1:d ;yq1(i+1,1)=(1/2)*(((1)^(i+1) - (-1)^(i+1))/(i+1));end % Moments of Mu_q
yq=[];
for i=1:nvar-nx
yq=[yq,yq1];
end

%%
Yx=[];
for i=1:size(vpow,1)
    g=glex2num([vpow(i,1:nx)]);
    Yx=[Yx;g];
  
        
    Yq(i,:)=1;
    for j=1:nvar-nx
        Yq(i,:)=Yq(i,:)*yq(vpow(i,nx+j)+1,j);
    end
end
Yx=Yx+Ny*ones(size(Yx));
Yx=Yx-Ny;
%ADDED BY SERHAT
indxq=cell(Nx,1);
for i=1:Nx
    indxq{i}=find(Yx==i);
end