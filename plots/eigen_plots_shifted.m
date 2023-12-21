%% Script to plot L2norm... and print convergence rates to screen!!
clear all

C4=load('../data_out/DM_eigens_4');C4=C4./60;
C6=load('../data_out/DM_eigens_6');C6=C6./60;
C8=load('../data_out/DM_eigens_8');C8=C8./60;
C10=load('../data_out/DM_eigens_10');C10=C10./60;



plot(C4(:,1),C4(:,2),'k.','linewidth',2,C6(:,1).+0.1,C6(:,2),'k.','linewidth',2,C8(:,1).+0.2,C8(:,2),'k.','linewidth',2,C10(:,1).+0.3,C10(:,2),'k.','linewidth',2);
axis('equal')
set(gca,'Fontsize',18)
set(gca,'Fontname','Times')
xlabel('Re \lambda');ylabel('Im \lambda');
print('-depsc','-S400,300',"eigens_shifted.eps")



