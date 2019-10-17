function xx = extract_mom(yx,nx,d)
 
mset clear; warning('off','YALMIP:strict')
mset('yalmip',true);mset(sdpsettings('solver','mosek')); 
mpol('x',1,nx)
p=sum(x.^2);
y=mom(mmon([x],2*d));
ys=yx;
P = msdp(min(p),y==ys);
[status,obj] = msol(P);
if status==1; clc; xx=double([x]); else
 disp('Increase the Relaxation Order to extract the solution');
end

end