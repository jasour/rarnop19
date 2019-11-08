function [Mind_x]=MomentIndex(nx,d)

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
