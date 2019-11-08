function  AtMtZ=AdjMoment_xq(Z,Yq,Mind_xq,AInd_y,Nx)

AtMtZ=zeros(Nx,1);
MtZ=AdjMoment(Z,AInd_y);
MtZ=MtZ.*Yq;
for i=1:Nx
    AtMtZ(i)=sum(MtZ(Mind_xq{i}));
end