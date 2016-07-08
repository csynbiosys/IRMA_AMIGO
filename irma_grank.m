%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TITLE: The circadian clock in Arabidopsis thaliana
%
%        Type :
%                > help circadian_tutorial
%        for a more detailed description of the model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%        INPUT FILE FOR GLOBAL RANK
%
%        This is the minimum input file for global rank. 
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
%                (AMIGO_GRank)==>> [n_s{iexp}]; [t_s{iexp}]; 
%                                  id_global_theta; [id_global_theta_y0]
%                                  [id_local_theta{iexp}];[id_local_theta_y0{iexp}]global_theta_max; global_theta_min
%                                  [global_theta_y0_max];[global_theta_y0_min]
%                                  [local_theta_max{iexp}];[local_theta_min{iexp}]
%                                  [local_theta_y0_max{iexp}];[local_theta_yo_min{iexp}]
%                                  []:optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%======================
%PATHS RELATED DATA
%======================

inputs.pathd.results_folder='irma_res';         % Folder to keep results (in Results) for a given problem          
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

inputs.model.par=[0,0.0404324055000000,1,0.0356000000000000,0.0221812275000000,0.000149286100000000,0.000882735200000000,0.0372000000000000,0.0477845334000000,0.201360986100000,0.00297081400000000,0.00223775600000000,0.200000000000000,0.0937500000000000,0.421662678000000,0.000740000000000000,0.0146950000000000,1.81400000000000,0.0980450000000000,0.167615000000000,0.000610000000000000,0.0181941480000000,1.81400000000000,0.0500000000000000,9,3,9;];

%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================


 inputs.exps.n_exp=2;                                  %Number of experiments                                                                            
 for iexp=1:inputs.exps.n_exp   
     inputs.exps.exp_y0{iexp}=[0.046735005043180 0.013409857422389 0.042059926621195 0.010943944574845 0.020445852394023];  %Initial conditions for each experiment          %Initial conditions for each experiment          
     inputs.exps.t_f{iexp}=3000;                            %Experiments duration

    % OBSEVABLES DEFINITION 
     inputs.exps.n_obs{iexp}=5;                       % Number of observed quantities per experiment  
     inputs.exps.obs_names{iexp}=char('CBF1','GAL4','SWI5','GAL80','ASH1');     % Name of the observed quantities per experiment    
     inputs.exps.obs{iexp}=char('CBF1=cbf1','GAL4=gal4','SWI5=swi5','GAL80=gal80','ASH1=ash1');   % Observation function
 
 
 end 
 
 inputs.exps.u_interp{1}='sustained';                  %Stimuli definition for experiment 1:
                                                       %OPTIONS:u_interp: 'sustained' |'step'|'linear'(default)|'pulse-up'|'pulse-down' 
 inputs.exps.t_con{1}=[0 3000];               % Input swithching times: Initial and final time    
 inputs.exps.u{1}=[1];                                 % Values of the inputs 
 
 inputs.exps.u_interp{2}='pulse-down';                 %Stimuli definition for experiment 2
 inputs.exps.n_pulses{2}=5;                            %Number of pulses |-|_|-|_|-|_|-|_|-|_    
 inputs.exps.u_min{2}=0;
 inputs.exps.u_max{2}=1;        %Minimum and maximum value for the input
 inputs.exps.t_con{2}=[0 :300 : 3000];                %Times of switching: Initial time, Intermediate times, Final time

 
% SAMPLING RELATED DATA (OPTIONAL)                                                             
 inputs.exps.n_s{1}=15;                                % [] Number of sampling times for each experiment.
 inputs.exps.n_s{2}=25;                                %    Optative input. By default "continuous" measurements are assumed.
% inputs.exps.t_s{1}=[1 2 3 ...];                      % [] Sampling times for each experiment, by default equidistant
% inputs.exps.t_s{2}=[0 5 7 ...];                      % [] Sampling times for each experiment, by default equidistant


%==================================
% UNKNOWNS RELATED DATA
%==================================

% GLOBAL UNKNOWNS (SAME VALUE FOR ALL EXPERIMENTS)
% to be considered in the rank

inputs.PEsol.id_global_theta='all';                     %  'all'|User selected 
inputs.PEsol.global_theta_max=2.*inputs.model.par; %100.*[7.5038 0.6801 1.4992 10.0982 2.3422 7.2482 1.8981 1.2 3.8045 ...
                                 %2.5356 1.4420 4.8600 1.2 9.4440 0.5 0.4364 7.3021 4.5703 1.0];    % Maximum allowed values for the paramters
inputs.PEsol.global_theta_min=0.5.*inputs.model.par;%(1/100).*[7.5038 0.6801 1.4992 10.0982 2.3422 7.2482 1.8981 1.2 3.8045 ...
                                 %2.5356 1.4420 4.8600 1.2 9.4440 0.5 0.4364 7.3021 4.5703 1.0];   % Minimum allowed values for the parameters
inputs.PEsol.global_theta_guess=inputs.model.par;      % [] Initial guess
                            
                             
% GLOBAL INITIAL CONDITIONS
%inputs.PEsol.id_global_theta_y0='none';               % [] 'all'|User selected| 'none' (default)
% inputs.PEsol.global_theta_y0_max=[];                % Maximum allowed values for the initial conditions
% inputs.PEsol.global_theta_y0_min=[];                % Minimum allowed values for the initial conditions
% inputs.PEsol.global_theta_y0_guess=[];              % [] Initial guess

% LOCAL UNKNOWNS (DIFFERENT VALUES FOR DIFFERENT EXPERIMENTS)

%inputs.PEsol.id_local_theta{1}='none';                % [] 'all'|User selected| 'none' (default)
% inputs.PEsol.local_theta_max{iexp}=[];              % Maximum allowed values for the paramters
% inputs.PEsol.local_theta_min{iexp}=[];              % Minimum allowed values for the parameters
% inputs.PEsol.local_theta_guess{iexp}=[];            % [] Initial guess
%inputs.PEsol.id_local_theta_y0{1}='none';             % [] 'all'|User selected| 'none' (default)
% inputs.PEsol.local_theta_y0_max{iexp}=[];           % Maximum allowed values for the initial conditions
% inputs.PEsol.local_theta_y0_min{iexp}=[];           % Minimum allowed values for the initial conditions
% inputs.PEsol.local_theta_y0_guess{iexp}=[];         % [] Initial guess

inputs.rank.gr_samples = 100;