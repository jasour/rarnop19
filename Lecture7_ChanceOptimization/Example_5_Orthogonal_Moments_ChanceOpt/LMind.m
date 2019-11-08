
function [LMm]=LMind(Dc,Cc,d,dc,m)

n2=size(Dc,2);
vpow=[];
for k = 0:2*d
    vpow = [vpow;genpow(n2,k)];
end

x=[];
for i=1:m
    in=[];
    for j=1:size(Dc,2)
        in=[in,num2str(vpow(i,j))];
    end
    eval(['syms',' ','m',in]);
    x=[x,eval(['m',in])];
end

Ccc=[];
for i=1:size(Cc,1)
    in=[];
    for j=1:size(Dc,2)
        in=[in,num2str(Dc(i,j))];
    end
    eval(['syms',' ','a',in]);
    Ccc=[Ccc;eval(['a',in])];
end

mBc=[];
for i=1:size(Dc,1)
    mBc=[mBc;glex2num(Dc(i,:))-1];
end
[Dc,Cc,mBc,Ccc];

LM=LMoment(Dc,Ccc,x,dc);
LMm=subs(LM,Ccc,Cc)
