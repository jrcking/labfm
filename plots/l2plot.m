%% Script to plot L2norm... and print convergence rates to screen!!
clear all

A=load('../data_out/L2norm');
figure(1)
%{
nn=length(A(:,1));
x=A(3:nn,1);
y=49*1e-16*x.**-2;
%}


loglog(A(:,1),A(:,3),'r*-')%,A(:,1),A(:,4),'ko-')%,A(:,1),A(:,4),'k*-',A(:,1),A(:,5),'k>-')%,x,y,'k--','linewidth',1.5)

%loglog(A(:,1),A(:,2),'k<:',A(:,1),A(:,3),'ko:',A(:,1),A(:,4),'k*-',A(:,1),A(:,5),'k>-',A(:,1),A(:,6),'k+-',A(:,1),A(:,7),'k^-',A(:,1),A(:,8),'kx-')

%loglog(A(:,1),A(:,2),'k<-',A(:,1),A(:,3),'ko-',A(:,1),A(:,4),'k*-',A(:,1),A(:,5),'k>-')

%loglog(A(:,1),A(:,2),'kx-','linewidth',1.5,A(:,1),A(:,3),'ko-','linewidth',1.5,A(:,1),A(:,4),'k^-','linewidth',1.5)%,x,y,'k--','linewidth',1.5)


set(gca,'ytick',[1e-15 1e-12 1e-9 1e-6 1e-3 1 1e3])

set(gca,'Fontsize',18)
set(gca,'Fontname','Times')
xlabel('h/H');ylabel('L_{2} norm');
%title('Error in Laplacian approximation')
%text(2e-7, 5e-10,'f(x)=sin(2\pi{x}/\lambda)cos(2\pi{y}/\lambda)','fontsize',16,'fontname','Times')
%text(1.1e-5, 1e-11,'WendlandC6,O(4),h/dx=3,ss=2','fontsize',12,'fontname','Times')

%[hleg1, hobj1] = legend('k=2','k=4','k=6','k=8');
%set(hleg1,'fontname','times');set(hleg1,'fontsize',12)
%set(hleg1,'position',[0.6 0.25 0.3 0.25])

%text(3e-4,7e-0,'h/\delta{r}=1','Fontsize',14,'Fontname','Times')

%text(5.0e-2,5e1,'Distribution noise: \epsilon/\delta{r}=0.2','Fontsize',14,'Fontname','Times')



print('-depsc','-S400,300',"convergence.eps")
[a,b]=size(A);
for i=1:(a-1)
for j=1:b
C(i,j) = log(A(i,j)/A(i+1,j))/log(A(i,1)/A(i+1,1));
endfor
endfor

C


