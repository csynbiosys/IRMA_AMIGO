%IRMA simulation script


clear;clc;

fileID = fopen('results_TSseries.txt','w');
randx=linspace(0,1,3);%change 3 to user-defined
M_sol = cell(length(randx), 1) ; %cell to store matrices
for n=1:length(randx)
    load K
    load y0glu

    tt=linspace(0,2000,20);

    u=rand(1,22)>randx(n)

    fprintf(fileID, '\n Interpolation %d \nu= ',n);
    fprintf(fileID,'%d ',u);

    for i=1:length(tt)-1
        %fprintf('Da %d a %d \n',tt(i),tt(i+1));
        sol = dde23(@(t,sol,Z) IRMA5b(t,sol,Z,K,u(i),1,tt(i+1)), [100] ,y0,[tt(i) tt(i+1)]);
        out=[sol.x;sol.y]';
        y0=sol;
        plot(out(:,1),out(:,2),out(:,1),out(:,3),out(:,1),out(:,4),out(:,1),out(:,5),out(:,1),out(:,6));
        legend('cbf1', 'gal4','swi5','gal80','ash1');
        axis([0 2000 0 0.05])
        title('State Evolution for IRMA')
        xlabel('time')
        ylabel('state concentration')
        
        pause(0.05)
    end
    M_sol{n}=sol.y'
    hold off
    filename=['irmaTS',num2str(randx(n)),'.png'];
    saveas(gcf,filename);
    
    auc_dat1=gp4_roc(n,out)
    auc_dat=roc_TSseries(M_sol,n)
    fprintf(fileID,'\nORDER \t\t\t AuC');
    %fprintf(fileID,'AUC= ');
    fprintf(fileID,'\n%s \t %f \n',num2str(n),auc_dat1);

    [nrows,ncols] = size(auc_dat);%print cell rows to file
    formatSpec =('\n%s \t %f \n');
    for row = 1:nrows
    fprintf(fileID,formatSpec,auc_dat{row,:});
    end


end
   
