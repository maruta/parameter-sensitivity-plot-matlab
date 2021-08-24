% Display parameter sensitivity plot for coupled spring-mass damper system.

clear
close all
syms s m1 m2 m3 k1 k2 k3 d1 d2 d3 x1 x2 x3 F
vars    = {m1,k1,d1,m2,k2,d2,m3,k3,d3};
varvals = {100,100,100,10,10,10,1,1,1};

sol=solve(...
    m1*s^2*x1==-k1*x1-d1*s*x1-k2*(x1-x2)-d2*s*(x1-x2),...
    m2*s^2*x2==-k2*(x2-x1)-d2*s*(x2-x1)-k3*(x2-x3)-d3*s*(x2-x3),...
    m3*s^2*x3==-k3*(x3-x2)-d3*s*(x3-x2)+F,{x1,x2,x3});

G = simplify(sol.x1/F)

mag = sensitibity_plot(G,logspace(-2,0,100),vars,varvals);

