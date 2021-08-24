% Display parameter sensitivity plot for spring-mass damper system.

clear
close all
syms s m k d
vars    = {m, k, d};
varvals = {1, 1, 1};

G = 1/(m*s^2+d*s+k)

mag = sensitibity_plot(G,logspace(-2,0,100),vars,varvals);
