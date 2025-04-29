# % MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
# Lecture 7: Chance Constrained/Chance Optimization
# Ashkan Jasour, Research Scientist, MIT 2020
#
# jasour.mit.edu  rarnop.mit.edu
# This code obtains the optimal values of the SDP in (3.7)
#  which approximates the chance optimization problem in (1.2)
#  e.q, Find x such that Probability of set { q : p_j(x,q)>0, j=1,...l } becomes maximum 
#       where q: random variable with probability distribution muq. 
#
#  "SEMIDEFINITE PROGRAMMING FOR CHANCE CONSTRAINED OPTIMIZATION OVER SEMIALGEBRAIC SETS"
# SIAM J. OPTIM. Vol. 25, No. 3, pp. 1411–1440
# ASHKAN. M. JASOUR, N. S. AYBAT, AND C. M. LAGOA
############################################################################

using MomentOpt
using DynamicPolynomials
using MosekTools
# using DelimitedFiles # to save the resuults in text file
# using ProxSDP # low-rank SDP solver

# number of decision variables
nx=1
# number of uncertain parameters 
nq=1
# relaxation order: SDP uses 2*d number of the moments of uncertainties
d = 10

# polnomial variables
@polyvar x1 q1

# nonlinear safety constraint
K = @set(((3*x1)/2 + 1)^3/4 - ((3*x1)/2 + 1)^2/4 - ((3*x1)/2 + 1)^4/16 - (9*q1^2)/100 + 29/400>=0)
# decision variables \in X
X = @set(1-x1^2>=0)

# define the optimization with SDP solver
  gmp = GMPModel(with_optimizer(Mosek.Optimizer))
# low-rank SDP solver: ProxSDP
# gmp = GMPModel(with_optimizer(ProxSDP.Optimizer,tol_primal=0.001, max_iter=10^4 ,log_verbose=true))

# measure mu, slack measure mu_s, and measure mu_x
@variable gmp mu Meas([x1,q1], support = K)
@variable gmp mu_s Meas([x1,q1])
@variable gmp mu_x Meas([x1], support = X )


# given moments of uncertain parameter q1 ~ mu_q=uniform[-1,1]
u=1; l=-1; yq1 = Array{Float64}(undef, 2*d+1,1); yq1[1]=1; for i in 1:2*d; yq1[i+1]= (1/(u-l))*( u^(i+1) - l^(i+1) )/(i+1); end 

# moments of mu, mu_s, and mu_x*mu_q
mono_xq = monomials([x1,q1],0:2*d); m_xq=length(mono_xq)
deg = Array{Int64}(undef, m_xq,2)
y=Mom(x1,mu)*zeros(m_xq,1)
y_s=Mom(x1,mu_s)*zeros(m_xq,1)
yxq=Mom(x1,mu_x)*zeros(m_xq,1)
for i in 1:m_xq; 
  deg[i,1:2]=exponents(mono_xq[i])'; 
  global y[i]=Mom(x1^(deg[i,1])*q1^(deg[i,2]),mu); 
  global y_s[i]=Mom(x1^(deg[i,1])*q1^(deg[i,2]),mu_s); 
  global yxq[i]=Mom(x1^(deg[i,1]),mu_x)*yq1[deg[i,2]+1]; 
end

# constrain: mu_x is a probability measure
@constraint gmp Mom(1, mu_x)==1

# constrain: mu <= mu_x*mu_q --->  (mu_x*mu_q - mu) is a measure mu_s
@constraint gmp y_s.==yxq - y

# obj: maximize the volume of mu
@objective gmp Max Mom(1, mu)

# moment SDP relaxation order
set_approximation_degree(gmp, 2*d)

# solve moment SDP
optimize!(gmp)

# Obtained results
optimum = objective_value(gmp);
mmu_x = measure(mu_x);
mmu = measure(mu);
mmu_s = measure(mu_s);

# Obtained moments of mu_x
yy_x = moment_value.(moments(mmu_x));
Lx = monomials(mmu_x); # order of moments

# Obtained moments of mu
yy = moment_value.(moments(mmu));
Ly = monomials(mmu); # order of moments

# Obtained moments of mu_s
yy_s = moment_value.(moments(mmu_s));


# extracted solution of chance optimization
x_sol = atomic(mu_x);

println("optimum: $(optimum)")
println("extracted solution: $(x_sol)")


# extract the solution of dual SDP: gamma>=0 in R and sos polynomial Pd(x1,q1)
# sos cond1: Pd(x1,q1)>=1 on set K
# sos cond2: gamma-int Pd(x1,q1) dmu_q >=0 on set X
# sos cond3: Pd(x1,q1)>=0
gamma=MomentOpt.approx_crefs(gmp)[1].dual[1];
Co= MomentOpt.approx_crefs(gmp)[1].dual[1]*zeros(m_xq,1)
for i=1:m_xq
Co[i]=MomentOpt.approx_crefs(gmp)[i+1].dual[1]
end
#Pd=Co'*Ly;
# save the Co as text to plot the Pd using Matlab
#open("Sol.txt", "w") do io; writedlm(io, Co, ','); end;