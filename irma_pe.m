
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
% PATHS RELATED DATA
%======================

inputs.pathd.results_folder='irma_res';         % Folder to keep results (in Results) for a given problem          
inputs.pathd.short_name='irma';                      % To identify figures and reports for a given problem   

%======================
% MODEL RELATED DATA
%======================

inputs.model.input_model_type='blackboxmodel';                % Model introduction: 'charmodelC'|'c_model'|'charmodelM'|'matlabmodel'|'sbmlmodel'|                        
                                                           %                     'blackboxmodel'|'blackboxcost                             
inputs.model.blackboxmodel_file='irmabbmodel'; % File including the simulation of the given model
inputs.model.n_st=5; 
inputs.model.n_par=25; 
inputs.model.n_stimulus=1; 
inputs.model.names_type='custom';

inputs.model.st_names=char('cbf1','gal4','swi5','gal80','ash1');     % --Names of the states
inputs.model.par_names=char('K1','K2','K3','K4','K5','K6','K7','K8',...
                            'K9','K11','K12','K13','K14','K15',...
                           'K16','K17','K18','K19','K21','K22',...
                            'K23','K24','K25','K26','K27');                           % Names of the parameters                                %-- Names of the parameters-%omitting
inputs.model.stimulus_names=char('galactose'); 
inputs.model.par=[0,0.0404324055000000,1,0.0356000000000000,0.0221812275000000,0.000149286100000000,0.000882735200000000,0.0372000000000000,0.0477845334000000,0.00297081400000000,0.00223775600000000,0.200000000000000,0.0937500000000000,0.421662678000000,0.000740000000000000,0.0146950000000000,1.81400000000000,0.0980450000000000,0.000610000000000000,0.0181941480000000,1.81400000000000,0.0500000000000000,9,3,9;];

%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================

 inputs.exps.n_exp=2; 

 %Number of experiments                                                                            
 for iexp=1:inputs.exps.n_exp   
     inputs.exps.exp_y0{iexp}=[0.046735005043180 0.013409857422389 0.042059926621195 0.010943944574845 0.020445852394023];  %Initial conditions for each experiment          %Initial conditions for each experiment          
     inputs.exps.t_f{iexp}=2000;                          %Experiments duration

     % OBSEVABLES DEFINITION
     inputs.exps.n_obs{iexp}=5;                      % Number of observed quantities per experiment  
     inputs.exps.obs_names{iexp}=char('CBF1','GAL4','SWI5','GAL80','ASH1');   % Name of the observed quantities per experiment    
     inputs.exps.obs{iexp}=char('CBF1=cbf1','GAL4=gal4','SWI5=swi5','GAL80=gal80','ASH1=ash1');   % Observation function
 
 end 
 
 inputs.exps.u_interp{1}='sustained';                  %Stimuli definition for experiment 1:
                                                       %OPTIONS:u_interp: 'sustained' |'step'|'linear'(default)|'pulse-up'|'pulse-down' 
 inputs.exps.t_con{1}=[0 2000];%-- change time 1000 to 2000                 % Input swithching times: Initial and final time    
 inputs.exps.u{1}=[1];                                 % Values of the inputs 
 
 inputs.exps.u_interp{2}='pulse-down';                 %Stimuli definition for experiment 2
 inputs.exps.n_pulses{2}=5;                            %Number of pulses |-|_|-|_|-|_|-|_|-|_    
 inputs.exps.u_min{2}=0;
 inputs.exps.u_max{2}=1;        %Minimum and maximum value for the input
 inputs.exps.t_con{2}=[0 :200 : 2000];            %Times of switching: Initial time, Intermediate times, Final time
                  

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
		0.0479914  0.0141036  0.0402424  0.0101766  0.021305
		0.0510189  0.0138973  0.0422218  0.0101657  0.0195381
		0.0414557  0.0132065  0.0395127  0.0112109  0.0205481
		0.0487489  0.0136068  0.0397129  0.0108469  0.0198892
		0.0474796  0.012882  0.0420738  0.010836  0.0207554
		0.0436743  0.0140052  0.0452839  0.0117205  0.019832
		0.0457226  0.012641  0.0404408  0.0111036  0.0209468
		0.0475331  0.0126932  0.0428411  0.0110522  0.0212017
		0.0550953  0.0128671  0.041585  0.0118127  0.0221959
		0.0532047  0.0114357  0.0444111  0.0105037  0.0202474
		0.0435787  0.0143743  0.0397679  0.0113251  0.0182598
		0.053825  0.0136279  0.0421282  0.0114009  0.0195876
		0.0484282  0.0129037  0.0432224  0.0108106  0.0218307
		0.0465857  0.0143287  0.0443758  0.011062  0.0193498
		0.0484032  0.0122623  0.0453094  0.010306  0.0214282
		0.0462541  0.0133413  0.0422405  0.0103158  0.0205727
		0.0464429  0.013248  0.0389208  0.0110013  0.0219146
		0.0502141  0.0136239  0.0404976  0.0113392  0.0184412
		0.0500257  0.0136196  0.0398257  0.0123587  0.0202438
		0.0500447  0.0128299  0.047006  0.010579  0.0192111
		0.0483022  0.0133897  0.0407642  0.0110465  0.0234187
		0.0439114  0.0132993  0.043634  0.0108988  0.0212895
		0.0484091  0.0138307  0.0416548  0.0098862  0.0218556
		0.0505426  0.0141429  0.0439297  0.0107037  0.0193641
		0.0478755  0.0141536  0.0404501  0.0099619  0.0199668
		];
    
inputs.exps.error_data{1}=[
		0.00233679  0.000670508  0.00210441  0.0005472  0.0010223
		0.00233679  0.000670508  0.00210441  0.0005472  0.0010223
		0.00233679  0.000670508  0.00210441  0.0005472  0.0010223
		0.00233679  0.000670508  0.00210441  0.0005472  0.0010223
		0.00233679  0.000670508  0.00210441  0.0005472  0.0010223
		0.00233679  0.000670508  0.00210441  0.0005472  0.0010223
		0.00233679  0.000670508  0.00210441  0.0005472  0.0010223
		0.00233679  0.000670508  0.00210441  0.0005472  0.0010223
		0.00233679  0.000670508  0.00210441  0.0005472  0.0010223
		0.00233679  0.000670508  0.00210441  0.0005472  0.0010223
		0.00233679  0.000670508  0.00210441  0.0005472  0.0010223
		0.00233679  0.000670508  0.00210441  0.0005472  0.0010223
		0.00233679  0.000670508  0.00210441  0.0005472  0.0010223
		0.00233679  0.000670508  0.00210441  0.0005472  0.0010223
		0.00233679  0.000670508  0.00210441  0.0005472  0.0010223
		0.00233679  0.000670508  0.00210441  0.0005472  0.0010223
		0.00233679  0.000670508  0.00210441  0.0005472  0.0010223
		0.00233679  0.000670508  0.00210441  0.0005472  0.0010223
		0.00233679  0.000670508  0.00210441  0.0005472  0.0010223
		0.00233679  0.000670508  0.00210441  0.0005472  0.0010223
		0.00233679  0.000670508  0.00210441  0.0005472  0.0010223
		0.00233679  0.000670508  0.00210441  0.0005472  0.0010223
		0.00233679  0.000670508  0.00210441  0.0005472  0.0010223
		0.00233679  0.000670508  0.00210441  0.0005472  0.0010223
		0.00233679  0.000670508  0.00210441  0.0005472  0.0010223
		];
		

inputs.exps.exp_data{2}=[
		0.0460599  0.0137624  0.04301  0.0110442  0.0199913
		0.049455  0.0133963  0.041785  0.0103805  0.0202865
		0.0460457  0.0133864  0.0424541  0.0114632  0.020728
		0.0512899  0.0130178  0.00604625  0.00832245  0.0139781
		0.0256884  0.0133947  0.0088605  0.00820158  0.0140714
		0.0116959  0.00836325  0.0189894  0.00936259  0.015595
		0.00740836  0.00652839  0.0162101  0.00897103  0.0141174
		0.0160568  0.0103493  0.0243704  0.00910407  0.016031
		0.0264987  0.00994428  0.00634178  0.00804146  0.0129979
		0.0170914  0.00923965  0.00820881  0.00804757  0.0130911
		0.0103314  0.00719943  0.0190612  0.00864365  0.0149353
		0.00974278  0.00659831  0.0131627  0.00969512  0.0151972
		0.0209122  0.00832796  0.0271592  0.00903258  0.013723
		0.0262237  0.0118562  0.00843464  0.00786317  0.0131973
		0.0128575  0.00988949  0.00690315  0.00773786  0.014879
		0.00980696  0.00736145  0.0154137  0.00823056  0.0142279
		0.0148442  0.00692539  0.0175132  0.00879969  0.016335
		0.0137628  0.0085122  0.0064222  0.00832649  0.0155401
		0.026182  0.0103366  0.00709476  0.00896482  0.0136025
		0.0101412  0.00878107  0.00715369  0.00799119  0.0137948
		0.00706048  0.00611413  0.0170033  0.00824223  0.013694
		0.0140413  0.00673781  0.0237348  0.00997991  0.0156779
		0.0177294  0.00820871  0.00802806  0.00887094  0.0159306
		0.017294  0.0106771  0.0066052  0.00800186  0.013719
		0.0136669  0.00811592  0.00836076  0.00730328  0.0136504
		];





inputs.exps.error_data{2}=[
		0.00247759  0.000677948  0.00210339  0.000547199  0.0010223
		0.00247759  0.000677948  0.00210339  0.000547199  0.0010223
		0.00247759  0.000677948  0.00210339  0.000547199  0.0010223
		0.00247759  0.000677948  0.00210339  0.000547199  0.0010223
		0.00247759  0.000677948  0.00210339  0.000547199  0.0010223
		0.00247759  0.000677948  0.00210339  0.000547199  0.0010223
		0.00247759  0.000677948  0.00210339  0.000547199  0.0010223
		0.00247759  0.000677948  0.00210339  0.000547199  0.0010223
		0.00247759  0.000677948  0.00210339  0.000547199  0.0010223
		0.00247759  0.000677948  0.00210339  0.000547199  0.0010223
		0.00247759  0.000677948  0.00210339  0.000547199  0.0010223
		0.00247759  0.000677948  0.00210339  0.000547199  0.0010223
		0.00247759  0.000677948  0.00210339  0.000547199  0.0010223
		0.00247759  0.000677948  0.00210339  0.000547199  0.0010223
		0.00247759  0.000677948  0.00210339  0.000547199  0.0010223
		0.00247759  0.000677948  0.00210339  0.000547199  0.0010223
		0.00247759  0.000677948  0.00210339  0.000547199  0.0010223
		0.00247759  0.000677948  0.00210339  0.000547199  0.0010223
		0.00247759  0.000677948  0.00210339  0.000547199  0.0010223
		0.00247759  0.000677948  0.00210339  0.000547199  0.0010223
		0.00247759  0.000677948  0.00210339  0.000547199  0.0010223
		0.00247759  0.000677948  0.00210339  0.000547199  0.0010223
		0.00247759  0.000677948  0.00210339  0.000547199  0.0010223
		0.00247759  0.000677948  0.00210339  0.000547199  0.0010223
		0.00247759  0.000677948  0.00210339  0.000547199  0.0010223
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
%inputs.PEsol.llk_type='homo_var';                     % [] To be defined for llk function, 'homo' | 'homo_var' | 'hetero' 



%==================================
 
 inputs.ivpsol.senssolver='fdsens2';

 inputs.ivpsol.rtol=1.0D-7;                            % [] IVP solver integration tolerances
 inputs.ivpsol.atol=1.0D-7; 

%{
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