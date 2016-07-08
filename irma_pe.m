
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TITLE: AMIGO implementation of IRMA codes from
%|          Code written by Gianfranco Fiore                     |
%|              gianfrancofiore@inwind.it                        | 
%
%        INPUT FILE TO ESTIMATE MODEL UKNOWNS
%
%        This is the minimum input file to simualate with real data.
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
%                (AMIGO_PE)==>>    n_s{iexp}; t_s{iexp}; 
%                                  data_type; noise_type; 
%                                  exp_data{iexp}; [error_data{iexp}]
%                                  id_global_theta; [id_global_theta_y0]
%                                  [id_local_theta{iexp}];[id_local_theta_y0{iexp}]global_theta_max; global_theta_min
%                                  [global_theta_y0_max];[global_theta_y0_min]
%                                  [local_theta_max{iexp}];[local_theta_min{iexp}]
%                                  [local_theta_y0_max{iexp}];[local_theta_yo_min{iexp}]
%                                  [global_theta_guess];[global_theta_y0_guess];
%                                  [local_theta_guess{iexp}];[local_theta_y0_guess{iexp}]
%                                  [PEcost_type];[lsq_type];[llk_type]
%                                  []:optional inputs
%                                  
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


 inputs.exps.n_exp=1;                                  %Number of experiments                                                                            
 for iexp=1:inputs.exps.n_exp   
     inputs.exps.exp_y0{iexp}=[0.046735005043180 0.013409857422389 0.042059926621195 0.010943944574845 0.020445852394023];  %Initial conditions for each experiment          %Initial conditions for each experiment          
     inputs.exps.t_f{iexp}=480;                           %Experiments duration

    % OBSEVABLES DEFINITION 
     inputs.exps.n_obs{iexp}=5;                     % Number of observed quantities per experiment  
     inputs.exps.obs_names{iexp}=char('CBF1','GAL4','SWI5','GAL80','ASH1');  % Name of the observed quantities per experiment    
     inputs.exps.obs{iexp}=char('CBF1=cbf1','GAL4=gal4','SWI5=swi5','GAL80=gal80','ASH1=ash1');   % Observation function
 
 
 end 
 
 inputs.exps.u_interp{1}='sustained';                  %Stimuli definition for experiment 1:
                                                       %OPTIONS:u_interp: 'sustained' |'step'|'linear'(default)|'pulse-up'|'pulse-down' 
 inputs.exps.t_con{1}=[0 480];                % Input swithching times: Initial and final time    
 inputs.exps.u{1}=[1];                                 % Values of the inputs 
 
 inputs.exps.u_interp{2}='pulse-down';                 %Stimuli definition for experiment 2
 inputs.exps.n_pulses{2}=2;                            %Number of pulses |-|_|-|_|-|_|-|_|-|_    
 inputs.exps.u_min{2}=0;
 inputs.exps.u_max{2}=1;        %Minimum and maximum value for the input
 inputs.exps.t_con{2}=[0 :120: 480];                 %Times of switching: Initial time, Intermediate times, Final time

 
%==================================
% EXPERIMENTAL DATA RELATED INFO
%==================================                                                            
                                                            
 inputs.exps.n_s{1}=25;                                % [] Number of sampling times for each experiment.
 inputs.exps.n_s{2}=25;                                %    Optative input. By default "continuous" measurements are assumed.
% inputs.exps.t_s{1}=[1 2 3 ...];                      % [] Sampling times for each experiment, by default equidistant
% inputs.exps.t_s{2}=[0 5 7 ...];                      % [] Sampling times for each experiment, by default equidistant


 inputs.exps.data_type='real';                         % Type of experimental data: 'real'|'pseudo'|'pseudo_pos'(>=0)  
 inputs.exps.noise_type='homo_var';                    % Type of experimental noise: Gaussian with zero mean and 
                                                       %                             Homoscedastic with constant variance: 'homo'
                                                       %                             Homoscedastic with varying variance:'homo_var'
                                                       %                             Heteroscedastic: 'hetero' 

                                                       
                                                       
inputs.exps.exp_data{1}=[
		0.0472798  0.0121953  0.0413647  0.0101386  0.0197962        % Experimental noise, n_s{iexp}x n_obs{iexp}
		0.0460758  0.013747  0.0424616  0.0115098  0.0196396
		0.0505431  0.0135944  0.0437605  0.0106818  0.0209236
		0.0416907  0.0135212  0.0406552  0.0123243  0.0195301
		0.0475498  0.013154  0.0436204  0.0109366  0.0214382
		0.048284  0.0142501  0.0416775  0.0108037  0.0188963
		0.0468362  0.0139757  0.0414207  0.0109913  0.0196974
		0.0498236  0.0125802  0.0427548  0.0110647  0.0225251
		0.0448943  0.0130599  0.0416163  0.0112454  0.0195794
		0.0468346  0.0136828  0.0407792  0.0123813  0.0199668
		0.0498825  0.0131507  0.0456452  0.00987892  0.0207264
		0.0446341  0.0141302  0.0421605  0.0103899  0.0208098
		0.0496661  0.012565  0.0420842  0.01057  0.0234405
		0.0425912  0.0137602  0.0412051  0.0111121  0.021876
		0.0462018  0.0138954  0.0389958  0.0107375  0.0199264
		0.0500876  0.0136433  0.0446877  0.0112525  0.0211266
		0.048026  0.0135757  0.0400105  0.0101835  0.0218996
		0.0467102  0.0124987  0.0422968  0.0115756  0.0203946
		0.0473261  0.0145251  0.0421077  0.0110905  0.017965
		0.0445762  0.0133109  0.0444085  0.0111245  0.0201303
		0.04319  0.0125554  0.04153  0.0102691  0.0191141
		0.0418479  0.0127455  0.0403789  0.0106169  0.0210422
		0.0481699  0.0124431  0.0408568  0.011035  0.0203372
		0.0422404  0.0127826  0.0423477  0.0105249  0.0215181
		0.047665  0.0138198  0.0416832  0.0113897  0.0194077
		];
    
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

		

inputs.exps.exp_data{2}=[
		0.0496739  0.0125419  0.04082  0.0115628  0.0197936
		0.0437744  0.0137822  0.0395302  0.0116137  0.021363
		0.050586  0.0129688  0.0401626  0.0107288  0.0205672
		0.0514153  0.0121795  0.00716841  0.00850717  0.0142983
		0.0148714  0.0082881  0.00433735  0.00805967  0.0117069
		0.00570201  0.00729452  0.0114616  0.00858453  0.0128263
		0.020116  0.00946619  0.0208229  0.010105  0.0170071
		0.0266408  0.011201  0.0341724  0.0103521  0.0166909
		0.0386894  0.0127911  0.00275728  0.00865394  0.0155906
		0.0145271  0.00796503  0.00807907  0.00839552  0.0139355
		0.0102698  0.00817283  0.0154739  0.0083929  0.0163299
		0.0173468  0.008978  0.0215862  0.00976161  0.014935
		0.0285415  0.0106741  0.0310968  0.0101257  0.0191293
		0.0269018  0.0129801  0.00582945  0.00745719  0.0141201
		0.00693953  0.00823404  0.00570441  0.00823932  0.0147786
		0.00373518  0.00796532  0.0161601  0.00872879  0.0146981
		0.0131608  0.00955634  0.0232557  0.00915307  0.0152997
		0.0307817  0.011556  0.00525604  0.00841801  0.014421
		0.0123473  0.00850488  0.00810656  0.0076637  0.0130728
		0.00904589  0.00673359  0.00910096  0.00916868  0.0132771
		0.00936323  0.00654273  0.0105508  0.0078392  0.0148042
		0.0161671  0.008294  0.0258392  0.00857785  0.0169119
		0.035088  0.0127443  0.00488102  0.00872065  0.016166
		0.0128629  0.00923731  0.00295915  0.00770742  0.0127535
		0.00742734  0.00755035  0.007203  0.00777749  0.0129253
		];



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
		0.0025788  0.000684793  0.00210302  0.000547198  0.00102229
		0.00254788  0.000684793  0.00210302  0.000547198  0.00102229
		0.00254788  0.000684793  0.00210302  0.000547198  0.00102229
		0.00254788  0.000684793  0.00210302  0.000547198  0.00102229
		0.00254788  0.000684793  0.00210302  0.000547198  0.00102229
		0.00254788  0.000684793  0.00210302  0.000547198  0.00102229
		0.00254788  0.000684793  0.00210302  0.000547198  0.00102229
		0.00254788  0.000684793  0.00210302  0.000547198  0.00102229
		];


    
    

%==================================
% UNKNOWNS RELATED DATA
%==================================

inputs.PEsol.id_global_theta='all';     		   %'all'|User selected 
inputs.PEsol.global_theta_max=5.*inputs.model.par; % Maximum allowed values for the parameters
inputs.PEsol.global_theta_min=0.2.*inputs.model.par; % Minimum allowed values for the parameters

       
% inputs.PEsol.global_theta_guess=[   5.082824695230836
%    0.745885509439732
%    7.750612254239743
%    3.167859540291811
%    2.476253863567990
%    6.008054327045663
%    4.696372709251819
%    3.346472917415490
%    3.185933721977777]';
% inputs.PEsol.global_theta_max=2.*inputs.PEsol.global_theta_guess;  % Maximum allowed values for the paramters
% inputs.PEsol.global_theta_min=0*ones(1,9); % Minimum allowed values for the paramters






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


%==================================
% COST FUNCTION RELATED DATA
%==================================
         
inputs.PEsol.PEcost_type='lsq';                       % 'lsq' (weighted least squares default) | 'llk' (log likelihood) | 'user_PEcost' 
inputs.PEsol.llk_type='homo_var';                     % [] To be defined for llk function, 'homo' | 'homo_var' | 'hetero' 


%{
%==================================
 %
 % SIMULATION
 inputs.ivpsol.ivpsolver='cvodes';                     % [] IVP solver: 'radau5'(default, fortran)|'rkf45'|'lsodes'|


 inputs.ivpsol.senssolver='cvodes';                    % [] Sensitivities solver: 'cvodes' (C)


 inputs.ivpsol.rtol=1.0D-7;                            % [] IVP solver integration tolerances
 inputs.ivpsol.atol=1.0D-7; 
% 
% %

% % OPTIMIZATION
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
                                                       
%inputs.nlpsol.eSS.log_var = 1:9;
inputs.nlpsol.eSS.maxeval = 10;
inputs.nlpsol.eSS.maxtime = 1000;

inputs.nlpsol.eSS.local.solver = 'fminsearch';
inputs.nlpsol.eSS.local.finish = 'fminsearch';

inputs.nlpsol.multi_starts=500;                        % [] Number of different initial guesses to run local methods in the multistart approach
inputs.nlpsol.multistart.maxeval = 100;              % Maximum number of function evaluations for the multistart
inputs.nlpsol.multistart.maxtime = 100;                % Maximum allowed time for the optimization

inputs.nlpsol.DE.NP = inputs.model.n_par*10;                           % Initial population size (around 10*npar)
inputs.nlpsol.DE.itermax = 100;                      % Maximum number of iteratios in DE
inputs.nlpsol.DE.F = 1; %0.75;  %1                   % F: DE-stepsize F ex [0, 2]
inputs.nlpsol.DE.CR =0.85;                           %CR: crossover probabililty constant ex [0, 1]
inputs.nlpsol.DE.strategy =2;                        % strategy       1 --> DE/best/1/exp                                            
                                                     %                2 --> DE/rand/1/exp           
                                                     %                3 --> DE/rand-to-best/1/exp   
                                                     %                4 --> DE/best/2/exp          
                                                     %                5 --> DE/rand/2/exp           


% 
% %==================================
% % RIdent or GRank DATA
% %==================================
% %
% 
% inputs.rid.conf_ntrials=500;                          % [] Number of trials for the robust confidence computation (default: 500)
% inputs.rank.gr_samples=10000;                         % [] Number of samples for global sensitivities and global rank within LHS (default: 10000)    
% 
% 
% %==================================
% % DISPLAY OF RESULTS
% %==================================
% 
% 
inputs.plotd.plotlevel='full';                        % [] Display of figures: 'full'|'medium'(default)|'min' |'noplot' 
%inputs.plotd.figsave=1;
% inputs.plotd.epssave=0;                              % [] Figures may be saved in .eps (1) or only in .fig format (0) (default)
% inputs.plotd.number_max_states=8;                    % [] Maximum number of states per figure
% inputs.plotd.number_max_obs=8;                       % [] Maximum number of observables per figure
% inputs.plotd.n_t_plot=100;                           % [] Number of times to be used for observables and states plots
% inputs.plotd.number_max_hist=8;                      % [] Maximum number of unknowns histograms per figure (multistart)
%inputs.plotd.nx_contour=100;                          % Number of points for plotting the contours x and y direction
%inputs.plotd.ny_contour=100;                          % ADVISE: >50

%}