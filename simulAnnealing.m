clear;clc;

fileID = fopen('res_SimulAnnealing.txt','w');

%define temperature
Tk=1;
Tmin= 0.1; %cooling temperature
alpha=0.9;
u=double(rand(1,22)>rand());
U_res=cell(133,1);
C_res=zeros(133,1);
n=1;
while(Tk>Tmin)
    
    i=1;
    
    while(i<50)
    %generate random starting point

        out1= irmasimul(u);


        % call obj function 
        c1 = 1-gp4_roc(n, out1);
        %create neighbour state:
        %call mutation(x)
        u2= mutation(u);


        out2= irmasimul(u2);
        %create output for state2
        %call obj function
        c2= 1-gp4_roc(n, out2);
        if (c2<c1)
            u=u2;
            c=c2;
            U_res{n}=u;
            C_res(n)=c2;

        else
           % call acceptance probability

            [u,c]= acptProb(c2,c1,Tk,u,u2);
            U_res{n}=u;
            C_res(n)=c;

        end
        
        i=i+1;
    end
    Tk=Tk*alpha;
    fprintf(fileID, 'u= %s',mat2str(U_res{n}));
    fprintf(fileID, '\t n=%d \t\t auc= %f\n',n,(1-C_res(n)));
        
    n=n+1;

end

