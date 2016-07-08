function auc_dat= roc_TSseries(M,n)
size_arr=[];
for l =1:n
    size_arr=[size_arr;length(M{l})]; %stores no. of timepoints
end

mydataAB= gp4grn(cell2mat(M),[],[size_arr],35) %changed itermax to 35 since hard computation gives errors
headers={'cbf1','gal4','swi5','gal80','ash1'};
r_headers(1:5,1)={'cbf1','gal4','swi5','gal80','ash1'};  
filename1 = [ 'Interactions_S', num2str(1:n),'.xls' ];
xlswrite(filename1,mydataAB,'Sheet1','B2');     %Write column header
xlswrite(filename1,r_headers,'Sheet1','A2');    
xlswrite(filename1,headers,'Sheet1','B1');   

load gs
AB=mydataAB';
AB=AB(:);
[xAB,yAB,~,aucAB]=perfcurve(gs,AB,1);
auc_dat=[{num2str(1:n),aucAB}]
img1=plot(xAB,yAB);
xlabel('False positive rate');
ylabel('True positive rate');
title(sprintf('ROC for %s ',num2str(1:n)));
filename2 = ['ROC_S',num2str(1:n),'.png'];
saveas(img1,filename2);

end

