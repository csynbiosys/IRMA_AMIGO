function out=irmasimul(u)

load K
load y0glu

tt=linspace(0,2000,20);

for i=1:length(tt)-1
    %fprintf('Da %d a %d \n',tt(i),tt(i+1));
    sol = dde23(@(t,sol,Z) IRMA5b(t,sol,Z,K,u(i),1,tt(i+1)), [100] ,y0,[tt(i) tt(i+1)]);
    out=[sol.x;sol.y]';
    y0=sol;
   % plot(out(:,1),out(:,2),out(:,1),out(:,3),out(:,1),out(:,4),out(:,1),out(:,5),out(:,1),out(:,6));
    %legend('cbf1', 'gal4','swi5','gal80','ash1');
    %axis([0 2000 0 0.05])
    %title('State Evolution for IRMA')
    %xlabel('time')
    %ylabel('state concentration')
    %pause(0.05)
end

end