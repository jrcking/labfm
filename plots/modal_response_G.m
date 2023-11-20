%% Script to plot L2norm... and print convergence rates to screen!!
clear all

A=load('../data_out/L2norm');
figure(1)
%%{
x=0:0.01:1;
y=x;%.**2;
%}


plot(A(:,1),A(:,2),'k.',x,y,'k-')

set(gca,'Fontsize',18)
set(gca,'Fontname','Times')
xlabel('k/k_{Ny}');ylabel('k_{eff}/k_{Ny}');
axis([0 1 0 1])
%[hleg1, hobj1] = legend('k=2','k=4','k=6','k=8');
%set(hleg1,'position',[0.7 0.25 0.2 0.25])

%line
%text(3e-4,7e-0,'h/\delta{r}=1','Fontsize',14,'Fontname','Times')

%text(5.0e-2,5e1,'Distribution noise: \epsilon/\delta{r}=0.2','Fontsize',14,'Fontname','Times')



print('-depsc','-S400,300',"modal_response_G.eps")

