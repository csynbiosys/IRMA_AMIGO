function [mydataAB,aucAB]= myROC(n,TimeseriesA,TimeseriesB,gs,xA,yA,xB,yB,xAB,yAB)
%Combine n timepoints from A and all from B
TimeseriesAxB=[TimeseriesA(1:n,1:5);TimeseriesB];
mydataAxB= gp4grn(TimeseriesAxB,[],[n,11])
headers={'gene1','gene2','gene3','gene4','gene5'};
r_headers(1:5,1)={'gene1','gene2','gene3','gene4','gene5'};
filename1 = [ 'InteractionsA', num2str(n), 'B','.xls' ];
xlswrite(filename1,mydataAxB,'Sheet1','B2');     %Write column header
xlswrite(filename1,r_headers,'Sheet1','A2');    
xlswrite(filename1,headers,'Sheet1','B1');   


AxB=mydataAxB';
AxB=AxB(:);
[xAxB,yAxB,~,aucAxB]=perfcurve(gs,AxB,1);
AzB=['A', num2str(n),'B'];

%plotting
hold on
imgAxB=plot(xAxB,yAxB,'-*black');
imgAxB=plot(xA,yA,'-ro');
imgAxB=plot(xB,yB,'-bo');
imgAxB=plot(xAB,yAB,'-*y');

xlabel('False positive rate');
ylabel('True positive rate');
legend(AzB,'A','B','AB');
title(sprintf('ROC for %s combined vs separate timeseries', AzB));
filename2 = [ 'A', num2str(n), 'B_vs A_B_AB','.png' ];
saveas(imgAxB,filename2);
mydataAB=mydataAxB;
aucAB=aucAxB;
end
