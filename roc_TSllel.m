function auc_dat= roc_TSllel(M_sol)


ind=perms(1:size(M_sol,1)); %method to take every permutation of the timeseries
auc_dat=[]; % cell array to store AuC of ROC curves
for i=1:length(ind)
    TS=M_sol(ind(i,:),:); %permutation step
    size_arr=[];
    for l =1:length(TS)
        size_arr=[size_arr;length(TS{l})]; %stores no. of timepoints
    end

    mydataAB= gp4grn(cell2mat(TS),[],[size_arr],35) %changed itermax to 35 since hard computation gives errors
    headers={'cbf1','gal4','swi5','gal80','ash1'};
    r_headers(1:5,1)={'cbf1','gal4','swi5','gal80','ash1'};  
    filename1 = [ 'Interactions', strcat(num2str(ind(i,:))),'.xls' ];
    xlswrite(filename1,mydataAB,'Sheet1','B2');     %Write column header
    xlswrite(filename1,r_headers,'Sheet1','A2');    
    xlswrite(filename1,headers,'Sheet1','B1');   

    load gs
    AB=mydataAB';
    AB=AB(:);
    [xAB,yAB,~,aucAB]=perfcurve(gs,AB,1);
    area_dat={num2str(ind(i,:)),aucAB};
    img1=plot(xAB,yAB);
    xlabel('False positive rate');
    ylabel('True positive rate');
    title('ROC for Timeseries AB' );
    filename2 = ['ROC',strcat(num2str(ind(i,:))),'.png' ];
    saveas(img1,filename2);
    auc_dat=[auc_dat; area_dat];
end

end