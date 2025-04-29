#############################################################################
# Chance Constrained Set = {all x: probability( g(x,w)>=0 ) >= 1-Delta }
#  where w is a random vector, x is design vector, g(x,w)>=0 is a nonlinear safety constraint, and Delta is acceptable risk level.
#############################################################################
# Lecture 7: Nonlinear Chance Constrained and Chance Optimization, 
# rarnop.mit.edu 
# Ashkan Jasour, Research Scientist, MIT 2020
# jasour.mit.edu 
###########################################################################

using MomentOpt
using DynamicPolynomials
using MosekTools
using DelimitedFiles # to save the resuults in text file


# relaxation order: SDP uses 2*d number of the moments of uncertainties.
d = 10

# polnomial variables
@polyvar x1 w1

# nonlinear safety constraint
K = @set(0.5*w1*(w1^2+(x1-0.5)^2)-(w1^4+w1^2*(x1-0.5)^2+(x1-0.5)^4)>=0)

# Lebesgue Measure over [-1,1]^2 to calculate the integral
B = @set(1-x1^2>=0 && 1-w1^2>=0)


#SDP solver
gmp = GMPModel(with_optimizer(Mosek.Optimizer))
# gmp = GMPModel(with_optimizer(ProxSDP.Optimizer,tol_primal=0.001, max_iter=10^4 ,log_verbose=true))

# measure mu, slack measure mu_s, and Lebesgue measure mu_Leb
@variable gmp mu Meas([x1,w1], support = K)
@variable gmp mu_s Meas([x1,w1], support = B)
mu_Leb = lebesgue_measure_box(variable_box(x1 => [-1, 1], w1 => [-1, 1]); normalize = true)

# obj: maximize the volume of mu
@objective gmp Max Mom(1, mu)

# measure constraint: mu <= given Lebesgue Measure ----> mu + mu_s = mu_Leb  
con = @constraint gmp mu + mu_s == mu_Leb

# moment SDP
set_approximation_degree(gmp, 2*d)
optimize!(gmp)

# status
println("termination_status: $(termination_status(gmp))")

# Polynomial Indicator function from the dual solution
P = MomentOpt.approx_crefs(gmp)[index(con)].dual

# moments of Uniform probability distribution over [-1,1] for uncertain variable w1~mu_w1=U[-1,1]
u=1; l=-1;
yw1 = Array{Float64}(undef, 2*d+1,1); yw1[1]=1; for i in 1:2*d; yw1[i+1]= (1/(u-l))*( u^(i+1) - l^(i+1) )/(i+1); end 

# vector of coefficients and monomials of polynomial indicator function P 
Co=coefficients(P) 
Mon = monomials([x1,w1],0:2*d); m=length(Mon)
# vector of powers of monomials of polynomial indicator function 
deg = Array{Int64}(undef, m,2)
# Chance Constrained Set = {x1: Pcc(x1) >= Delta  } where Pcc(x1)=integral P(x1,w1) d(mu_w1)
Pcc= 0; for i in 1:m; deg[i,1:2]=exponents(Mon[i])'; global Pcc+=Co[i]*x1^(deg[i,1])*yw1[deg[i,2]+1]; end

println("Chance Constrained Set={all x1: probability( g(x1,w1)>=0 ) >= 1-Delta}")
println("Chance Constrained Set={all x1: Pcc(x1)>= 1-Delta}")
println("Pcc(x1)=$(Pcc)")

# save the Co as text to plot the Pd using Matlab
open("Sol_Co.txt", "w") do io; writedlm(io, Co, ','); end;