function auc= gp4_roc(n,out)

mydat= gp4grn(out(:,2:6),[],20, 20);

headers={'cbf1','gal4','swi5','gal80','ash1'};
r_headers(1:5,1)={'cbf1','gal4','swi5','gal80','ash1'};

filename1 = [ 'Interactions', num2str(n),'.xls' ];
xlswrite(filename1,mydat,'Sheet1','B2');
xlswrite(filename1,headers,'Sheet1','B1');
xlswrite(filename1,r_headers,'Sheet1','A2');

mydat=mydat';
mydat=mydat(:);


load gs
[x_roc,y_roc,~,auc]=perfcurve(gs,mydat,1)
img1=plot(x_roc,y_roc);
xlabel('False positive rate');
ylabel('True positive rate');
title('ROC');
filename2 = ['ROC',num2str(n),'.png'];
saveas(img1,filename2);

