function [AInd_y,Mind_y,AInd_x,Mind_x]=MomentIndex(nvar,nx,Ny,Nx,d)
vpow=[];
for k = 0:d
    vpow = [vpow;genpow(nvar,k)];
end

for i=1:size(vpow,1)   
    for j=1:i
        clc
        disp('Moment Matrix')
        disp([i,j,size(vpow,1)])
        
       Mind_y(i,j)=(glex2num(vpow(i,:)+vpow(j,:)));
       Mind_y(j,i)=Mind_y(i,j);
    end
end

AInd_y=zeros(Ny,size(vpow,1));
for i=1:Ny
    clc;disp('Adjoint Moment Matrix');disp([i,Ny])
    ind=find(Mind_y==i);
AInd_y(i,1:size(ind,1))=ind';
end
%%
vpow=[];
for k = 0:d
    vpow = [vpow;genpow(nx,k)];
end

Mind_x=[];
for i=1:size(vpow,1)   
    for j=1:i
        clc;disp('Moment Matrix x');disp([i,j,size(vpow,1)])
       Mind_x(i,j)=(glex2num(vpow(i,:)+vpow(j,:)));   
       Mind_x(j,i)=Mind_x(i,j);
    end
end

AInd_x=zeros(Nx,size(vpow,1));
for i=1:Nx
    clc;disp('Adjoint Moment Matrix x');disp([i,Nx])
    ind=find(Mind_x==i);
AInd_x(i,1:size(ind,1))=[ind'];
end

Mind_x=Mind_x+Ny*ones(size(Mind_x));
Mind_x=Mind_x-Ny;