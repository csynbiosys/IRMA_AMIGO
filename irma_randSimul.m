%IRMA simulation script for generating random inputs

clear;clc;

fileID = fopen('results_rand.txt','w');
randx=linspace(0,1,25)
for n=20:length(randx)
    load K
    load y0glu

    tt=linspace(0,3000,20);

    u=rand(1,22)>randx(n)

    fprintf(fileID, '\nu= ');
    fprintf(fileID,'%d ',u);

    for i=1:length(tt)-1
        sol = dde23(@(t,sol,Z) IRMA5b(t,sol,Z,K,u(i),1,tt(i+1)), [100] ,y0,[tt(i) tt(i+1)]);
        out=[sol.x;sol.y]';
        plot(out(:,1),out(:,2),out(:,1),out(:,3),out(:,1),out(:,4),out(:,1),out(:,5),out(:,1),out(:,6));
        legend('cbf1', 'gal4','swi5','gal80','ash1');
        axis([0 3000 0 0.05])
        title('State Evolution for IRMA')
        xlabel('time')
        ylabel('state concentration')
        y0=sol;
        pause(0.05)
    end
    hold off
    filename=['irma',num2str(randx(n)),'.png'];
    saveas(gcf,filename);
    
    auc=gp4_roc(randx(n),out)
    fprintf(fileID,'AUC= ');
    fprintf(fileID,'%f\n ',auc);

end
