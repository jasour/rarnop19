
function [LMInd,LMm]=LMind(Dc,Cc,d,dc,m)

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


Yl=[];
syms y
for i=1:size(LM,1)*size(LM,2)
yl=y^i;
    Yl=[Yl,yl'];
end
Yl=reshape(Yl,size(LM));
Yl=Yl';
T=trace(Yl*LM);

A=[];
for i=1: m
    in=[];
    for j=1:size(Dc,2)
        in=[in,num2str(vpow(i,j))];
    end
    mm=eval(['m',in]);
    [C1,C2]=coeffs(T,mm) ;
    if size(C2,2)==1
        A=[A;0];
    else
        A=[A;(C1(1))];
    end
end
A2=subs(A,Ccc,Cc);

LMInd=zeros(size(A2,1),size(LM,1)*size(LM,2));
for i=1:size(A2,1)
    S=[];
    S=rot90(sym2poly(A2(i)),2);
    S=S(1,2:end);
   LMInd(i,1:size(S,2))= S  ;
end 