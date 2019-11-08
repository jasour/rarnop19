function LM=LMoment(Dc,Cc,x,d1)

n2=size(Dc,2);
vpow=[];
for k = 0:d1
    vpow = [vpow;genpow(n2,k)];
end

for i=1:size(vpow,1)   
    for j=1:i
        clc;disp({'LM',j,i,size(vpow,1)})
        a=vpow(i,:)+vpow(j,:);
        b=0;
        for k=1:size(Dc)
            b=b+Cc(k)*x(glex2num(a+Dc(k,:)));
        end
        LM(i,j)=b;
        clear b a
        LM(j,i)=LM(i,j);
    end
end