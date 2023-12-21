%% Script to plot L2norm... and print convergence rates to screen!!
clear all

C=load('../data_out/DM_eigens_8');

%% NOISE FREE...
figure(1)
plot(C(:,1),C(:,2),'k.','linewidth',2);
%axis([-2 2 -40 40])
set(gca,'Fontsize',18)
set(gca,'Fontname','Times')
xlabel('Re \lambda');ylabel('Im \lambda');
print('-depsc','-S400,300',"eigens.eps")



