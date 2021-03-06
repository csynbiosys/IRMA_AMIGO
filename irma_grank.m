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

irma_model  %Loads model related data

%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
 inputs.exps.n_exp=4; 

 %Number of experiments                                                                            
 for iexp=1:inputs.exps.n_exp   
     inputs.exps.n_obs{iexp}=5;                      % Number of observed quantities per experiment  
     inputs.exps.obs_names{iexp}=char('CBF1','GAL4','SWI5','GAL80','ASH1');   % Name of the observed quantities per experiment    
     inputs.exps.obs{iexp}=char('CBF1=cbf1','GAL4=gal4','SWI5=swi5','GAL80=gal80','ASH1=ash1');   % Observation function
     inputs.exps.exp_y0{iexp}=[0.046735005043180 0.013409857422389 0.042059926621195 0.010943944574845 0.020445852394023];  %Initial conditions for each experiment          %Initial conditions for each experiment          
     inputs.exps.t_f{iexp}=2000;                          %Experiments duration
     
     
     
 end 


 inputs.exps.u_interp{1}='sustained';                  %Stimuli definition for experiment 1:
                                                       %OPTIONS:u_interp: 'sustained' |'step'|'linear'(default)|'pulse-up'|'pulse-down' 
 inputs.exps.t_con{1}=[0 2000];%-- change time 1000 to 2000                 % Input swithching times: Initial and final time    
 inputs.exps.u{1}=[1]; 
  % Values of the inputs 
 
 inputs.exps.u_interp{2}='pulse-down';                 %Stimuli definition for experiment 2
 inputs.exps.n_pulses{2}=2;                            %Number of pulses |-|_|-|_|-|_|-|_|-|_    
 inputs.exps.u_min{2}=0;
 inputs.exps.u_max{2}=1;        %Minimum and maximum value for the input
 inputs.exps.t_con{2}=[0 500 1000 1500 2000];            %Times of switching: Initial time, Intermediate times, Final time
                  
 inputs.exps.u_interp{3}='pulse-up';                 %Stimuli definition for experiment 2
 inputs.exps.n_pulses{3}=5;                            %Number of pulses |-|_|-|_|-|_|-|_|-|_    
 inputs.exps.u_min{3}=0;
 inputs.exps.u_max{3}=1;        %Minimum and maximum value for the input
 inputs.exps.exp_y0{3}=[0.00922634499300000 0.00679535315300000 0.00704562182800000 0.00812744071900000 0.0204458591000000];%Initial conditions for each experiment                        
 inputs.exps.t_con{3}=[0 200 400 600 800 1000 1200 1400 1600 1800 1999.9999 2000];            %Times of switching: Initial time, Intermediate times, Final time          %Times of switching: Initial time, Intermediate times, Final time

%u= [1 1 1 1 0 0 0 0 1 1 0 1 0 1 0 0 0 1 1]
 inputs.exps.u_interp{4}='pulse-down';                 %Stimuli definition for experiment 2
 inputs.exps.n_pulses{4}=6;                            %Number of pulses |-|_|-|_|-|_|-|_|-|_    
 inputs.exps.u_min{4}=0;
 inputs.exps.u_max{4}=1;        %Minimum and maximum value for the input
 inputs.exps.t_con{4}=[0 421.0526 842.1053 1052.6 1157.9 1263.2 1368.4 1473.7  1578.9 1684.2 1789.5 1999.999 2000];            %Times of switching: Initial time, Intermediate times, Final time          %Times of switching: Initial time, Intermediate times, Final time

 
% SAMPLING RELATED DATA (OPTIONAL)      
   
inputs.exps.n_s{1}=30;                                % [] Number of sampling times for each experiment.
inputs.exps.n_s{2}=30;                                %    Optative input. By default "continuous" measurements are assumed.
inputs.exps.n_s{3}=30
inputs.exps.n_s{4}=30


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
