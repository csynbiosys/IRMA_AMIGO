                                 %IRMA5.m
%                      Smith predictor/FeedBack Control/Relay
%=========================================================================
%|                   Code written by Gianfranco Fiore                     |
%|                       gianfrancofiore@inwind.it                        | 
%=========================================================================      
%|                     Stand Alone Control Algorithm                      |
%|                            code version 0.0.3                          |
%|                               24/07/2011                               |
%=========================================================================

function dy = IRMA5b(t,y,Z,par,u)
%IRMA
%alfa=teta(2);
%K=[0,0.0404324055000000,1,0.0356000000000000,0.0221812275000000,0.000149286100000000,0.000882735200000000,0.0372000000000000,0.0477845334000000,0.201360986100000,0.00297081400000000,0.00223775600000000,0.200000000000000,0.0937500000000000,0.421662678000000,0.000740000000000000,0.0146950000000000,1.81400000000000,0.0980450000000000,0.167615000000000,0.000610000000000000,0.0181941480000000,1.81400000000000,0.0500000000000000,9,3,9;];
%variabili da minimizzare con il Simulated annealing
%disp(size(y))
K1=par(1);
K2=par(2);
K3=par(3);
K4=par(4);
K5=par(5);
K6=par(6);
K7=par(7);
K8=par(8);
K9=par(9);
K10=par(10);
K11=par(11);
K12=par(12);
K13=par(13);
K14=par(14);
K15=par(15);
K16=par(16);
K17=par(17);
K18=par(18);
K19=par(19);
K20=par(20);
K21=par(21);
K22=par(22);
K23=par(23);
K24=par(24);
K25=par(25);
K26=par(26);
K27=par(27);
K28= 0.0900;
K29= 0.0154;
K30= 0.0154;

%u=u(1)+(t-tlast)*pend(1);
% var=load('galtime.mat');
% galtime=var.galtime;
% u=1/2*(sign(t+galtime*60)-sign(t-galtime*60));

%u=1;

v=1;
deg=1;
alfa3=1;
deg3=1;
alfa4=1;
deg4=1;
alfa6=1;
deg6=1;
dec=1;
% CBF1 mRNA
dy(1,1)= ((K1+ v*K2* (Z(3)/( (K3+Z(3)) *( 1+(y(5)/K4))  ))  -deg* K5*y(1))); 
% GAL4 mRNA
dy(2,1) = (K6+K7*(y(1) /(K8 +y(1)) ) -  K9 *y(2));
% SWI5 mRNA (note that the values of 3 parameters change depending on the medium, galactose
% or glucose)
dy(3,1) =(K11*alfa3+ (K12*(1-u)+u*(K12*K25) )  *(y(2).^4./(  (K14*(1-u)+u*(K14/K27) ).^4 +y(2).^4  .*( 1+(y(4).^4./(   (K13*(1-u)+u*(K13*K26) ).^4)  ))  )    ) -deg3*K15*y(3));
% GAL80 mRNA
dy(4,1) =(K16+K17*(y(3)/(K18+y(3)) )- K19 *y(4) );
% ASH1 mRNA
dy(5,1) =(K21*alfa4+ K22*(y(3)/(K23+y(3)))-K24*deg4*y (5));

%disp(size(dy))
end