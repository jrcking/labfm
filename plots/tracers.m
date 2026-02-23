%%
hold off
% Script to plot some particle paths

Nt=100; % Number of tracers
figure(1),hold on

for i=1:Nt
   A=load(strcat('../fort.',num2str(900+i-1)));
   plot(A(:,1),A(:,2),'.','linewidth',1,'color',[i/Nt 0 1-i/Nt])
   

endfor
