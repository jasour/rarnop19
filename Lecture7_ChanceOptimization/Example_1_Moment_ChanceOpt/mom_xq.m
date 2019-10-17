function yxq = mom_xq(nx,nq,d,yq,yx)
 
%vpow: powers of (nx+nq) dimensional monomials of order k
vpow=[];for k = 0:2*d; vpow = [vpow;genpow(nx+nq,k)]; end
Yx=[];
for i=1:size(vpow,1)
Yx=[Yx;glex2num([vpow(i,1:nx)])];
Yq(i,:)=1; for j=1:nq; Yq(i,:)=Yq(i,:)*yq(vpow(i,nx+j)+1,j); end
end
yxq=yx(Yx).*Yq; 

end