
%======================
% PATHS RELATED DATA
%======================

inputs.pathd.results_folder='irma_test';         % Folder to keep results (in Results) for a given problem          
inputs.pathd.short_name='irma';                      % To identify figures and reports for a given problem   

%======================
% MODEL RELATED DATA
%======================

inputs.model.input_model_type='blackboxmodel';                % Model introduction: 'charmodelC'|'c_model'|'charmodelM'|'matlabmodel'|'sbmlmodel'|                        
                                                           %                     'blackboxmodel'|'blackboxcost                             
inputs.model.blackboxmodel_file='irmabbmodel';%-- change model name % File including the simulation of the given model
inputs.model.n_st=5; %-- changed states from 1 to 5
inputs.model.n_par=27; %- ??do we need this? changing par to K
inputs.model.n_stimulus=1; %-- changed from 1 to 2
inputs.model.names_type='custom';

inputs.model.st_names=char('cbf1','gal4','swi5','gal80','ash1');     % --Names of the states
inputs.model.par_names=char('K1','K2','K3','K4','K5','K6','K7','K8',...
                            'K9','K10','K11','K12','K13','K14','K15',...
                           'K16','K17','K18','K19','K20','K21','K22',...
                            'K23','K24','K25','K26','K27');                           % Names of the parameters                                %-- Names of the parameters-%omitting
inputs.model.stimulus_names=char('galactose'); %--added glucose
%-- changing par to K (taken from IRMA5)
inputs.model.par=[0,0.0404324055000000,1,0.0356000000000000,0.0221812275000000,0.000149286100000000,0.000882735200000000,0.0372000000000000,0.0477845334000000,0.201360986100000,0.00297081400000000,0.00223775600000000,0.200000000000000,0.0937500000000000,0.421662678000000,0.000740000000000000,0.0146950000000000,1.81400000000000,0.0980450000000000,0.167615000000000,0.000610000000000000,0.0181941480000000,1.81400000000000,0.0500000000000000,9,3,9;];


%=====================================
%  EXPERIMENTAL SCHEME RELATED DATA
%=====================================
 inputs.exps.n_exp=3;                                 % Total number of experiments (available + experiments to be designed)     
 inputs.exps.exp_type{1}='fixed';                     % Indicates if the the experiment should
 inputs.exps.exp_type{2}='fixed';                     % or not be optimally designed : 'od' and 'fixed'
 inputs.exps.exp_type{3}='od';             



 %inputs.exps.n_exp=2;                                  %Number of experiments                                                                            
 for iexp=1:inputs.exps.n_exp   
              
    % OBSEVABLES DEFINITION 
     inputs.exps.n_obs{iexp}=5; %--- changed 1 to 5                       % Number of observed quantities per experiment  
     inputs.exps.obs_names{iexp}=char('CBF1','GAL4','SWI5','GAL80','ASH1'); % --added observables    % Name of the observed quantities per experiment    
     inputs.exps.obs{iexp}=char('CBF1=cbf1','GAL4=gal4','SWI5=swi5','GAL80=gal80','ASH1=ash1');   % Observation function
     inputs.exps.exp_y0{iexp}=[0.046735005043180 0.013409857422389 0.042059926621195 0.010943944574845 0.020445852394023];  %Initial conditions for each experiment
 
 end 
 
 
 inputs.exps.u_interp{1}='sustained';                  % Stimuli definition for experiment 1:
                                                       % OPTIONS:u_interp: 'sustained' |'step'|'linear'(default)|'pulse-up'|'pulse-down' 
 inputs.exps.t_f{1}=480;                               % Experiment duration
 inputs.exps.t_con{1}=[0 480];                         % Input swithching times: Initial and final time    
 inputs.exps.u{1}=[1];                                 % Values of the inputs 
 inputs.exps.n_s{1}=25;                                % Number of sampling times for each experiment.
 
  inputs.exps.error_data{1}=[                                      % Experimental noise, n_s{iexp}x n_obs{iexp}
		0.00233776  0.000670633  0.00210508  0.000547285  0.0010224
		0.00233776  0.000670633  0.00210508  0.000547285  0.0010224
		0.00233776  0.000670633  0.00210508  0.000547285  0.0010224
		0.00233776  0.000670633  0.00210508  0.000547285  0.0010224
		0.00233776  0.000670633  0.00210508  0.000547285  0.0010224
		0.00233776  0.000670633  0.00210508  0.000547285  0.0010224
		0.00233776  0.000670633  0.00210508  0.000547285  0.0010224
		0.00233776  0.000670633  0.00210508  0.000547285  0.0010224
		0.00233776  0.000670633  0.00210508  0.000547285  0.0010224
		0.00233776  0.000670633  0.00210508  0.000547285  0.0010224
		0.00233776  0.000670633  0.00210508  0.000547285  0.0010224
		0.00233776  0.000670633  0.00210508  0.000547285  0.0010224
		0.00233776  0.000670633  0.00210508  0.000547285  0.0010224
		0.00233776  0.000670633  0.00210508  0.000547285  0.0010224
		0.00233776  0.000670633  0.00210508  0.000547285  0.0010224
		0.00233776  0.000670633  0.00210508  0.000547285  0.0010224
		0.00233776  0.000670633  0.00210508  0.000547285  0.0010224
		0.00233776  0.000670633  0.00210508  0.000547285  0.0010224
		0.00233776  0.000670633  0.00210508  0.000547285  0.0010224
		0.00233776  0.000670633  0.00210508  0.000547285  0.0010224
		0.00233776  0.000670633  0.00210508  0.000547285  0.0010224
		0.00233776  0.000670633  0.00210508  0.000547285  0.0010224
		0.00233776  0.000670633  0.00210508  0.000547285  0.0010224
		0.00233776  0.000670633  0.00210508  0.000547285  0.0010224
		0.00233776  0.000670633  0.00210508  0.000547285  0.0010224
		];
  
 
 
 inputs.exps.u_interp{2}='pulse-down';                 % Stimuli definition for experiment 2
 inputs.exps.n_pulses{2}=2;                            % Number of pulses |-|_|-|_|-|_|-|_|-|_    
 inputs.exps.u_min{2}=0;inputs.exps.u_max{2}=1;        % Minimum and maximum value for the input
 inputs.exps.t_f{2}=480;                               % Experiment duration
 inputs.exps.t_con{2}=[0 :120: 480];                    % Times of switching: Initial time, Intermediate times, Final time
 inputs.exps.n_s{2}=25;                                % Number of sampling times for each experiment.
 
  
inputs.exps.error_data{2}=[
		0.00254788  0.000684793  0.00210302  0.000547198  0.00102229
		0.00254788  0.000684793  0.00210302  0.000547198  0.00102229
		0.00254788  0.000684793  0.00210302  0.000547198  0.00102229
		0.00254788  0.000684793  0.00210302  0.000547198  0.00102229
		0.00254788  0.000684793  0.00210302  0.000547198  0.00102229
		0.00254788  0.000684793  0.00210302  0.000547198  0.00102229
		0.00254788  0.000684793  0.00210302  0.000547198  0.00102229
		0.00254788  0.000684793  0.00210302  0.000547198  0.00102229
		0.00254788  0.000684793  0.00210302  0.000547198  0.00102229
		0.00254788  0.000684793  0.00210302  0.000547198  0.00102229
		0.00254788  0.000684793  0.00210302  0.000547198  0.00102229
		0.00254788  0.000684793  0.00210302  0.000547198  0.00102229
		0.00254788  0.000684793  0.00210302  0.000547198  0.00102229
		0.00254788  0.000684793  0.00210302  0.000547198  0.00102229
		0.00254788  0.000684793  0.00210302  0.000547198  0.00102229
		0.00254788  0.000684793  0.00210302  0.000547198  0.00102229
		0.00254788  0.000684793  0.00210302  0.000547198  0.00102229
		0.00254788  0.000684793  0.00210302  0.000547198  0.00102229
		0.00254788  0.000684793  0.00210302  0.000547198  0.00102229
		0.00254788  0.000684793  0.00210302  0.000547198  0.00102229
		0.00254788  0.000684793  0.00210302  0.000547198  0.00102229
		0.00254788  0.000684793  0.00210302  0.000547198  0.00102229
		0.00254788  0.000684793  0.00210302  0.000547198  0.00102229
		0.00254788  0.000684793  0.00210302  0.000547198  0.00102229
		0.00254788  0.000684793  0.00210302  0.000547198  0.00102229
		];
 
 %{
 inputs.exps.u_interp{1}='sustained';                  %Stimuli definition for experiment 1:
                                                       %OPTIONS:u_interp: 'sustained' |'step'|'linear'(default)|'pulse-up'|'pulse-down' 
 inputs.exps.t_con{1}=[0 2000];%-- change time 1000 to 2000                 % Input swithching times: Initial and final time    
 inputs.exps.u{1}=[1];                                 % Values of the inputs 
 
 inputs.exps.u_interp{2}='pulse-down';                 %Stimuli definition for experiment 2
 inputs.exps.n_pulses{2}=5;                            %Number of pulses |-|_|-|_|-|_|-|_|-|_    
 inputs.exps.u_min{2}=0;
 inputs.exps.u_max{2}=1;        %Minimum and maximum value for the input
 inputs.exps.t_con{2}=[0 :200 : 2000];  %--add time                 %Times of switching: Initial time, Intermediate times, Final time

 
 
 %}
 
 
 
 
 
 
  
%
%  INPUTS FOR THE EXPERIMENT TO BE OPTIMALLY DESIGNED
%
 
 inputs.exps.u_type{3}='od';                           % Type of stimulation: 'fixed' | 'od' (to be designed) 
 inputs.exps.u_interp{3}='step';                       % Stimuli definition for experiment 3:
                                                       % OPTIONS:u_interp:
                                                       % 'sustained' |'step'|'linear'(default)|'pulse-up'|'pulse-down' 
 inputs.exps.n_steps{3}=3;                             % Number of pulses _|-|_|-|_ 
 inputs.exps.u_min{3}=0*ones(1,inputs.exps.n_steps{3});
 inputs.exps.u_max{3}=1*ones(1,inputs.exps.n_steps{3});% Minimum and maximum value for the input
 inputs.exps.tf_type{3}='fixed';                       % [] Type of experiment duration: 'fixed'(default) | 'od' (to be designed) 
 inputs.exps.t_f{3}=480;                               % Experiment duration
 inputs.exps.ts_type{3}='fixed';                       % [] Type of sampling times: 'fixed'(default) | 'od' (to be designed) 
 inputs.exps.n_s{3}=25;
 inputs.exps.std_dev{3}=0.1;                           % Standard deviation of the noise for each experiment: Ex: 0.05 <=> 5%


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

inputs.nlpsol.eSS.maxeval = 10000;
inputs.nlpsol.eSS.maxtime = 60;
inputs.nlpsol.eSS.local.solver = 'fmincon';
inputs.nlpsol.eSS.local.finish = 'fmincon';
                                                       
inputs.nlpsol.multi_starts=500;                        % [] Number of different initial guesses to run local methods in the multistart approach
inputs.nlpsol.multistart.maxeval = 20000;            % Maximum number of function evaluations for the multistart
inputs.nlpsol.multistart.maxtime = 120;                % Maximum allowed time for the optimization
inputs.nlpsol.eSS.local.nl2sol.maxiter             =     300;        % max number of iteration
inputs.nlpsol.eSS.local.nl2sol.maxfeval            =     400;         % max number of function evaluation
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

