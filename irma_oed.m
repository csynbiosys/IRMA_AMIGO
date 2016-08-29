irma_model  %Loads model related data
%======================
%PATHS RELATED DATA
%======================

inputs.pathd.results_folder='irma_res_Thesis';         % Folder to keep results (in Results) for a given problem          
inputs.pathd.short_name='irma';                      % To identify figures and reports for a given problem   


%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
 inputs.exps.n_exp=5; 
                                % Total number of experiments (available + experiments to be designed)     
 inputs.exps.exp_type{1}='fixed';                     % Indicates if the the experiment should
 inputs.exps.exp_type{2}='fixed';
 inputs.exps.exp_type{3}='fixed'
 inputs.exps.exp_type{4}='fixed'% or not be optimally designed : 'od' and 'fixed'
 inputs.exps.exp_type{5}='od';   

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
inputs.exps.n_s{5}=30

 inputs.exps.data_type='real';                         % Type of experimental data: 'real'|'pseudo'|'pseudo_pos'(>=0)  
 inputs.exps.noise_type='homo_var';                    % Type of experimental noise: Gaussian with zero mean and 
                                                       %                             Homoscedastic with constant variance: 'homo'
                                                       %                             Homoscedastic with varying variance:'homo_var'
                                                       %                             Heteroscedastic: 'hetero' 
temp=load(strcat('C:\Users\s1565325\Documents\MATLAB\AMIGO2R2016a\Results\irma_13Aug_Thesis\SData_irma_psdata\strreport_irma_psdata.mat'),'results'); 

%EXPERIMENT 1

inputs.exps.exp_data{1}=temp.results.sim.exp_data{1};
inputs.exps.error_data{1}=temp.results.sim.error_data{1};

%EXPERIMENT 2

inputs.exps.exp_data{2}=temp.results.sim.exp_data{2};
inputs.exps.error_data{2}=temp.results.sim.error_data{2};

%EXPERIMENT 3

inputs.exps.exp_data{3}=temp.results.sim.exp_data{3};
inputs.exps.error_data{3}=temp.results.sim.error_data{3};

%EXPERIMENT 4

inputs.exps.exp_data{4}=temp.results.sim.exp_data{4};
inputs.exps.error_data{4}=temp.results.sim.error_data{4};

clear temp;


%  INPUTS FOR THE EXPERIMENT TO BE OPTIMALLY DESIGNED
%
 
 inputs.exps.u_type{5}='od';                           % Type of stimulation: 'fixed' | 'od' (to be designed) 
 inputs.exps.u_interp{5}='pulse-down';                       % Stimuli definition for experiment 3:
                                                       % OPTIONS:u_interp:
                                                       % 'sustained' |'step'|'linear'(default)|'pulse-up'|'pulse-down' 
 inputs.exps.n_pulses{5}=3;                             % Number of pulses _|-|_|-|_ 
 inputs.exps.u_min{5}=0*ones(1,inputs.exps.n_pulses{5});
 inputs.exps.u_max{5}=1*ones(1,inputs.exps.n_pulses{5});% Minimum and maximum value for the input
 inputs.exps.tf_type{5}='fixed';                       % [] Type of experiment duration: 'fixed'(default) | 'od' (to be designed) 
 inputs.exps.t_f{5}=2000;                               % Experiment duration
 inputs.exps.ts_type{5}='fixed';                       % [] Type of sampling times: 'fixed'(default) | 'od' (to be designed) 
 inputs.exps.std_dev{5}=0.1;                           % Standard deviation of the noise for each experiment: Ex: 0.05 <=> 5%


%======================================
% PARAMETERS TO BE CONSIDERED FOR OED
%======================================


% Parameters to be considered for OED 
inputs.PEsol.id_global_theta='all';  %  'all'|User selected  
% Nominal value of the parameters to compute the FIM
inputs.PEsol.global_theta_guess=inputs.model.par; 


%==================================
% COST FUNCTION RELATED DATA
%==================================    

 inputs.exps.noise_type='homo_var';                    % Type of experimental noise: 'homo' |'homo_var'| 'hetero'     
 inputs.OEDsol.OEDcost_type='Eopt';                % FIM based criterium: 'Dopt'|'Eopt'|'Aopt'|'Emod'|'DoverE'
 
 
%==================================
% NUMERICAL METHODS RELATED DATA
%==================================
%
% SIMULATION

 inputs.ivpsol.senssolver='fdsens2';

 inputs.ivpsol.rtol=1.0D-7;                            % [] IVP solver integration tolerances
 inputs.ivpsol.atol=1.0D-7; 
 
 
% OPTIMIZATION
%
inputs.nlpsol.nlpsolver='eSS';                        % [] NLP solver: 
%                                                       % LOCAL: 'local_fmincon'|'local_n2fb'|'local_dn2fb'|'local_dhc'|
%                                                       %        'local_ipopt'|'local_solnp'|'local_nomad'|'local_fsqp'|'local_misqp'|'local_nl2sol'
%                                                       %        'local_lsqnonlin'
%                                                       % MULTISTART:'multi_fmincon'|'multi_n2fb'|'multi_dn2fb'|'multi_dhc'|
%                                                       %            'multi_ipopt'|'multi_solnp'|'multi_nomad'|'multi_fsqp'|'multi_misqp'|'multi_nl2sol'
%                                                       %            'multi_lsqnonlin'
%                                                       % GLOBAL: 'de'|'sres'
%                                                       % HYBRID: 'hyb_de_fmincon'|'hyb_de_n2fb'|'hyb_de_dn2fb'|'hyb_de_dhc'|'hyp_de_ipopt'|
%                                                       %         'hyb_de_solnp'|'hyb_de_nomad'|'hyb_de_fsqp'|'hyb_de_misqp'
%                                                       %         'hyb_sres_fmincon'|'hyb_sres_n2fb'|'hyb_sres_dn2fb'|'hyb_sres_dhc'|
%                                                       %         'hyp_sres_ipopt'|'hyb_sres_solnp'|'hyb_sres_nomad'|'hyb_sres_fsqp'|'hyb_sres_misqp'
%                                                       % METAHEURISTICS:
%                                                       % 'ssm'(DEFAULT)|'fssm'
%                                                       % Note that the corresponding defaults are in files: 
%                                                       % OPT_solvers\DE\de_options.m; OPT_solvers\SRES\sres_options.m; 
%                                                       % OPT_solvers\SSm_**\ssm_options.m 
%                                                       

inputs.nlpsol.eSS.maxeval = 100000;
inputs.nlpsol.eSS.maxtime = 600;
inputs.nlpsol.eSS.local.solver = 'local_dhc';
inputs.nlpsol.eSS.local.finish = 'local_dhc';
                                                       
inputs.nlpsol.multi_starts=500;                        % [] Number of different initial guesses to run local methods in the multistart approach
inputs.nlpsol.multistart.maxeval = 100000;            % Maximum number of function evaluations for the multistart
inputs.nlpsol.multistart.maxtime = 300;                % Maximum allowed time for the optimization
inputs.nlpsol.eSS.local.nl2sol.maxiter             =     500;        % max number of iteration
inputs.nlpsol.eSS.local.nl2sol.maxfeval            =     500;         % max number of function evaluation
% 
% 
% %==================================
% % DISPLAY OF RESULTS
% %==================================
% 
 inputs.plotd.plotlevel='full';                       % [] Display of figures: 'full'|'medium'(default)|'min' |'noplot' 
% inputs.plotd.epssave=0;                              % [] Figures may be saved in .eps (1) or only in .fig format (0) (default)
% inputs.plotd.number_max_states=8;                    % [] Maximum number of states per figure
% inputs.plotd.number_max_obs=8;                       % [] Maximum number of observables per figure
% inputs.plotd.n_t_plot=100;                           % [] Number of times to be used for observables and states plots
% inputs.plotd.contour_rtol=1e-7;                      % [] Integration tolerances for the contour plots. 
% inputs.plotd.contour_atol=1e-7;                      %    ADVISE: These tolerances should be a little bit strict
%  inputs.plotd.nx_contour=100;                          % [] Number of points for plotting the contours x and y direction
%  inputs.plotd.ny_contour=100;                          %    ADVISE: >=50
% inputs.plotd.number_max_hist=8;                      % [] Maximum number of unknowns histograms per figure (multistart)
% 
%  
