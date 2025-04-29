# Risk = probability (g(x)>=0) where x is a random vector and g(x)>=0 is a nonlinear safety constraint.
#############################################################################
# measure-moment SDP formulation:
# max_{mu} vol(mu)                     : obj
#   s.t mu supported on g(x)>=0        : support constraint
#       mu <= given Lebesgue Measure   : measure constraint
# where mu is unknown measure defined on safety set g(x)>=0.
# Risk = Expected value of {Polynomial Indicator function of the set {x:g(x)>=0} } with respect to the probability distribution of random vector x.
# Solution of the dual SDP is the coefficients of the polynomial indicator function.
# Lecture 10: Probabilistic Nonlinear Safety Verification, rarnop.mit.edu 
# Ashkan Jasour, Research Scientist, MIT 2020
# jasour.mit.edu  rarnop.mit.edu
###########################################################################
# MomentOpt julia package
# https://github.com/lanl-ansi/MomentOpt.jl
###########################################################################
# We use Low-Rank Relaxation to solve the SDP
# https://github.com/mariohsouto/ProxSDP.jl
###########################################################################

using MomentOpt
using DynamicPolynomials
using MosekTools
using ProxSDP


# relaxation order: SDP uses 2*d number of the moments of uncertainties to calculate the Risk.
d = 5

# polnomial variables
@polyvar x1 x2

# nonlinear safety constraint
K = @set(-x1^4+0.5*(x1^2-x2^2)+0.1>=0)

# Lebesgue Measure over [-1,1]^2 to calculate the integral
B = @set(1-x1^2>=0 && 1-x2^2>=0)


#SDP solver
#gmp = GMPModel(with_optimizer(Mosek.Optimizer))
 gmp = GMPModel(with_optimizer(ProxSDP.Optimizer,tol_primal=0.001, max_iter=10^4 ,log_verbose=true))

# measure mu, slack measure mu_s, and Lebesgue measure mu_Leb
@variable gmp mu Meas([x1,x2], support = K)
@variable gmp mu_s Meas([x1,x2], support = B)
mu_Leb = lebesgue_measure_box(variable_box(x1 => [-1, 1], x2 => [-1, 1]); normalize = true)

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

# moments of Uniform probability distribution over [0,1] for uncertain variable x1
u=1; l=0;
yx1 = Array{Float64}(undef, 2*d+1,1); yx1[1]=1; for i in 1:2*d; yx1[i+1]= (1/(u-l))*( u^(i+1) - l^(i+1) )/(i+1); end 

# moments of Beta(aB,bB) probability distribution over [0,1] for uncertain variable x2
aB=1.3;bB=3;
yx2 = Array{Float64}(undef, 2*d+1,1);yx2[1]=1;for i in 1:2*d; yx2[i+1]=(aB+i-1)/(aB+bB+i-1)*yx2[i]; end; 

# vector of coefficients and monomials of polynomial indicator function P 
Co=coefficients(P) 
Mon = monomials([x1,x2],0:2*d); m=length(Mon)
# vector of powers of monomials of polynomial indicator function 
deg = Array{Int64}(undef, m,2)

# Risk = Expected value of {Polynomial Indicator function of the set {x:g(x)>=0} } with respect to the probability distribution of random vector x.
#      = Replace the monomials with the moments of probability distribution of random vector x.
Risk=Array{Float64}(undef, 1,1);Risk[1]=0
for i in 1:m
     deg[i,1:2]=exponents(Mon[i])'
     global Risk.+=Co[i]*yx1[deg[i,1]+1]*yx2[deg[i,2]+1]
end

# Result
println("Upper Bound on Risk: $(Risk)")
