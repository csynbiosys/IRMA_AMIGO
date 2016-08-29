clear;clc;
load K
load y0glu

tt=linspace(0,2000,20)
sol = dde23(@(t,sol,Z) IRMA5b(t,sol,Z,K,1,1,t), [100] ,y0,[0 3000]);
disp(sol.y)
out=[sol.x;sol.y]';
plot(out(:,1),out(:,2));
axis([0 3000 0 0.05])
pause
hold on
u=repmat([1 0],[1 int8((length(tt)-1)/2)])
for i=1:length(tt)-1
    fprintf('Da %d a %d \n',tt(i),tt(i+1));
    sol = dde23(@(t,sol,Z) IRMA5b(t,sol,Z,K,u(i)), [100] ,y0,[tt(i) tt(i+1)]);out=[sol.x;sol.y]';plot(out(:,1),out(:,2));axis([0 3000 0 0.05])
    y0=sol;
    pause(0.05)
end
hold off
