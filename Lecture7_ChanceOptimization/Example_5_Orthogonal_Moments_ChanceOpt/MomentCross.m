function MAx= MomentCross(x,Yx,Yq)

global Mind_y
Yxq=x(Yx)'.*Yq; 
MAx= Yxq(Mind_y);