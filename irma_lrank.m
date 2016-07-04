% $Header: svn://172.19.32.13/trunk/AMIGO2R2016/Examples/Arabidopsis_circadian/circadian_sobs.m 2398 2015-12-04 07:06:07Z evabalso $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TITLE: The circadian clock in Arabidopsis thaliana
%
%        Type :
%                > help circadian_tutorial
%        for a more detailed description of the model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       INPUT FILE TO GENERATE PSEUDO-EXPERIMENTAL DATA (USEFUL FOR
%        NUMERICAL EXAMPLES)
%        This is the minimum input file to generate pseudo-experimental data 
%        Default values are assigned to non defined inputs.
%
%        Minimum required inputs:
%           > Paths related data
%           > Model:               model_type; n_st; n_par; n_stimulus; 
%                                  st_names; par_names; stimulus_names;  
%                                  eqns; par
%           > Experimental scheme: n_exp; exp_y0{iexp}; t_f{iexp}; 
%                                  u_interp{iexp}; t_con{iexp}; u{iexp}
%                                  n_obs{iexp}; obs_names{iexp}; obs{iexp} 
%
%                (AMIGO_SData)==>> n_s{iexp}; t_s{iexp}; 
%                                  data_type; noise_type; std_dev{iexp}
%                                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%======================
%PATHS RELATED DATA
%======================

inputs.pathd.results_folder='irma_test';         % Folder to keep results (in Results) for a given problem          
inputs.pathd.short_name='irma';                      % To identify figures and reports for a given problem   

%======================
%MODEL RELATED DATA
%======================

inputs.model.input_model_type='blackboxmodel';                % Model introduction: 'charmodelC'|'c_model'|'charmodelM'|'matlabmodel'|'sbmlmodel'|                        
                                                                              % 'blackboxmodel'|'blackboxcost                             
inputs.model.blackboxmodel_file='irmabbmodel';
inputs.model.n_st=5; 
inputs.model.n_par=27;
inputs.model.n_stimulus=1; 
inputs.model.names_type='custom';

inputs.model.st_names=char('cbf1','gal4','swi5','gal80','ash1');     % --Names of the states
inputs.model.par_names=char('K1','K2','K3','K4','K5','K6','K7','K8',...
                            'K9','K10','K11','K12','K13','K14','K15',...
                           'K16','K17','K18','K19','K20','K21','K22',...
                            'K23','K24','K25','K26','K27');                           % Names of the parameters                                %-- Names of the parameters-%omitting
inputs.model.stimulus_names=char('galactose'); 

%{
inputs.model.eqns=...                                      % Equations describing system dynamics. Time derivatives are regarded 'd'st_name''
               char('dy(1,1)= ((K1+ v*K2* (Z(3)/( (K3+Z(3)) *( 1+(y(5)/K4))  ))  -deg* K5*y(1)))',...
                    'dy(2,1) = (K6+K7*(y(1) /(K8 +y(1)) ) -  K9 *y(2))',...
                    'dy(3,1) =(K11*alfa3+ (K12*(1-u)+u*(K12*K25) )  *(y(2).^4./(  (K14*(1-u)+u*(K14/K27) ).^4 +y(2).^4  .*( 1+(y(4).^4./(   (K13*(1-u)+u*(K13*K26) ).^4)  ))  )',...
                    'dy(4,1) =(K16+K17*(y(3)/(K18+y(3)) )- K19 *y(4) )',...
                    'dy(5,1) =(K21*alfa4+ K22*(y(3)/(K23+y(3)))-K24*deg4*y (5))'); 
                
%}
inputs.model.par=[0,0.0404324055000000,1,0.0356000000000000,0.0221812275000000,0.000149286100000000,0.000882735200000000,0.0372000000000000,0.0477845334000000,0.201360986100000,0.00297081400000000,0.00223775600000000,0.200000000000000,0.0937500000000000,0.421662678000000,0.000740000000000000,0.0146950000000000,1.81400000000000,0.0980450000000000,0.167615000000000,0.000610000000000000,0.0181941480000000,1.81400000000000,0.0500000000000000,9,3,9;];

%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================


 inputs.exps.n_exp=2;                                  %Number of experiments                                                                            
 for iexp=1:inputs.exps.n_exp   
     inputs.exps.exp_y0{iexp}=[0.046735005043180 0.013409857422389 0.042059926621195 0.010943944574845 0.020445852394023];  %Initial conditions for each experiment          %Initial conditions for each experiment          
     inputs.exps.t_f{iexp}=3000; %-- change time to 2000                           %Experiments duration

    % OBSEVABLES DEFINITION ---OMITTING ---- 
     inputs.exps.n_obs{iexp}=5; %--- changed 1 to 5                       % Number of observed quantities per experiment  
     inputs.exps.obs_names{iexp}=char('CBF1','GAL4','SWI5','GAL80','ASH1'); % --added observables    % Name of the observed quantities per experiment    
     inputs.exps.obs{iexp}=char('CBF1=cbf1','GAL4=gal4','SWI5=swi5','GAL80=gal80','ASH1=ash1');   % Observation function
 
 
 end 
 
 inputs.exps.u_interp{1}='sustained';                  %Stimuli definition for experiment 1:
                                                       %OPTIONS:u_interp: 'sustained' |'step'|'linear'(default)|'pulse-up'|'pulse-down' 
 inputs.exps.t_con{1}=[0 3000];%-- change time 1000 to 2000                 % Input swithching times: Initial and final time    
 inputs.exps.u{1}=[1];                                 % Values of the inputs 
 
 inputs.exps.u_interp{2}='pulse-down';                 %Stimuli definition for experiment 2
 inputs.exps.n_pulses{2}=5;                            %Number of pulses |-|_|-|_|-|_|-|_|-|_    
 inputs.exps.u_min{2}=0;
 inputs.exps.u_max{2}=1;        %Minimum and maximum value for the input
 inputs.exps.t_con{2}=[0 :300 : 3000];  %--add time                 %Times of switching: Initial time, Intermediate times, Final time

 
% SAMPLING RELATED DATA (OPTIONAL)      

inputs.exps.n_s{1}=25;                                % [] Number of sampling times for each experiment.
inputs.exps.n_s{2}=25;                                %    Optative input. By default "continuous" measurements are assumed.
% inputs.exps.t_s{1}=[1 2 3 ...];                      % [] Sampling times for each experiment, by default equidistant
% inputs.exps.t_s{2}=[0 5 7 ...];                      % [] Sampling times for each experiment, by default equidistant

%==================================
% UNKNOWNS RELATED DATA
%==================================

% GLOBAL UNKNOWNS (SAME VALUE FOR ALL EXPERIMENTS)
% to be considered in the local rank

inputs.PEsol.id_global_theta='all';                   % 'all'|User selected 
%inputs.PEsol.global_theta_guess=[];                   % To calculate the rank for this value of the parameters

% GLOBAL INITIAL CONDITIONS
% inputs.PEsol.id_global_theta_y0=char();              % [] 'all'|User selected| 'none' (default)




