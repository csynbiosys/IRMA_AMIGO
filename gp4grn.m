function vertices = gp4grn(timeseries, steady, n_timepoints, itermax, delay, zom, max_reg, B, b)
%GP4GRN  GRN inference method based on the use of ODE and Gaussian processes
%   
%   usage: vertices = gp4grn(timeseries, steady, n_timepoints, itermax, delay, zom, max_reg, B, b)
%
%   where:
%       vertices        is the connectivity matrix (element [vertices]_{ij}
%                       is the probability that gene i regulates gene j)       
%       timeseries      is a K by N matrix of combined timeseries data (N
%                       is the number of genes and K is the number of 
%                       timepoints in M timeseries experiments), give [] in
%                       case of no timeseries data
%       steady          is a L by N matrix of combined steady-state data (N
%                       is the number of genes and L is the number of
%                       steady-state experiments), give [] in case of no
%                       steady-state data
%       n_timepoints    number of timepoints in each of the timeseries 
%                       measurements (1 by M vector, default value is the
%                       number of rows in the timeseries matrix)
%       itermax         is the number of maximum iterations in the
%                       optimization (default is 50)
%       delay           is the delay in the model for timeseries data
%                       (default is 0, which should be a safe choice) 
%       zom             use zero-order model (no=0 and yes=1, default is 1))
%       max_reg         maximum number of regulators per gene (default is 3)
%       B               covariance matrix of the linear regression 
%                       coefficients, (basal and degradation, 2 by 2, 
%                       default is [0.01 0;0 0.01], you should really 
%                       change this)
%       b               mean vector of the linear regression coefficients,
%                       (basal and degradation, 2 by 1, default is [0; 0])
%
%   The GPML Matlab code is written by Carl Edward Rasmussen.
%
%   Author: Tarmo Äijö <tarmo.aijo@tut.fi>
%
%   Updated: 24.3.2009

% default parameter values
if nargin < 2, error('not enough parameters'), end
if nargin < 9, b = [0; 0]; end
if nargin < 8, B = [0.01 0;0 0.01]; end
if nargin < 7, max_reg = 3; end
if nargin < 6, zom = 1; end
if nargin < 5, delay = 0; end
if nargin < 4, itermax = 50; end
if nargin < 3, n_timepoints = size(timeseries,1); end

if isempty(timeseries) & isempty(steady)
    error('timeseries and steady are empty matrices')
elseif ~isempty(timeseries) & ~isempty(steady) % timeseries and steady
    if size(timeseries, 2) ~= size(steady, 2)
        error('timeseries and steady matrices should have same number of columns')
    end
    n_variables = size(timeseries, 2);
    use_data = 2;
elseif isempty(timeseries) % only steady
    n_variables = size(steady, 2);
    use_data = 1;
else % only timeseries
    n_variables = size(timeseries, 2);
    use_data = 0;
end

% calculate the connectivity matrix gene by gene (this loop can be easily
% distributed)
vertices = zeros(n_variables);
for variable_i=1:n_variables
    tic
    fprintf('response gene: %d of %d\n',variable_i,n_variables)
    vertices = vertices + gp4grn_caller(variable_i,timeseries,steady,...
        use_data,n_timepoints,itermax,delay,zom,max_reg,n_variables,B,b);
    fprintf('time consumed: %.2f seconds\n\n',toc)
end
vertices = vertices';
end

function result = gp4grn_caller(variable,timeseries,steady,use_data,n_timepoints,itermax,delay,zom,max_reg,n_variables,B,b)
%GP4GRN_CALLER  Convert the data suitable for the model and calls GP4GRN_RES function
%       
%   Author: Tarmo Äijö <tarmo.aijo@tut.fi>

indices = 1:n_variables ~= variable; % indices of explanatory variables
%% regression (GP)
if use_data == 0 % only timeseries data
    t = [0:size(timeseries,1)-1]'; % time
    % approximate the derivatives
    diff_approx = approx_diff(timeseries(:,variable), t); 
    % response data (remove first 1+delay samples from each of the separate
    % timeseries) 
    diff_approx = trunc(diff_approx', n_timepoints, -delay)'; 
    % response data (degradation term) (remove first 1+delay samples from
    % each of the separate timeseries)
    x = trunc(timeseries(:,variable)', n_timepoints, -delay)';
    % explanatory data (remove last 1+delay samples from each of the
    % separate timeseries)
    Y = trunc(timeseries(:,indices)', n_timepoints, delay)'; 
    result=gp4grn_res(diff_approx, x, Y, itermax, zom, indices, variable, max_reg, B, b);
elseif use_data == 1 % only steady-state data
    % response data, derivatives are ~0
    diff_approx = zeros(size(steady,1), 1); 
    x = steady(:,variable); % response data (degradation term)
    Y = steady(:,indices); % explanatory data
    result=gp4grn_res(diff_approx, x, Y, itermax, zom, indices, variable, max_reg, B, b);
elseif use_data == 2 % timeseries and steady-state data
    t = [0:size(timeseries,1)-1]'; % time
    % approximate the derivatives
    diff_approx = approx_diff(timeseries(:,variable), t); 
    % response data (remove first 1+delay samples from each of the 
    % separate timeseries)
    diff_approx = [trunc(diff_approx', n_timepoints, -delay)'; zeros(size(steady,1), 1)];
    % response data (degradation term) (remove first 1+delay samples from
    % each of the separate timeseries)
    x = [trunc(timeseries(:,variable)', n_timepoints, -delay)'; steady(:,variable)]; 
    Y = [trunc(timeseries(:,indices)', n_timepoints, delay)'; steady(:,indices)];
    % explanatory data (remove last 1+delay samples from each of the
    % separate timeseries)
    result=gp4grn_res(diff_approx, x, Y, itermax, zom, indices, variable, max_reg, B, b);
end
end

function result = gp4grn_res(diff_approx, x, Y, itermax, zom, indices_orig, variable_i, max_reg, B, b)
%GP4GRN_RES  Calculate probabilities for the interactions to a given gene
%       
%   Author: Tarmo Äijö <tarmo.aijo@tut.fi>

% the linear part of the regression
H = [ones(size(x)) x]'; % training set for beta

% truncated powerset of the transcription factor candidates set
subps = powerset((1:size(Y,2))',max_reg);

% calculate negative log-marginal likelihood of the zero-order model
% (no regulators)
zom_nllh = [];
if zom == 1
    loghyper = -1;
    loghyper = minimize(loghyper, 'gpr', -itermax, 'covNoise', ...
        zeros(size(diff_approx)), diff_approx, H, B, b);
    zom_nllh = gpr(loghyper, 'covNoise', zeros(size(diff_approx)), ...
        diff_approx, H, B, b);
end

% calculate negative log-likelihoods of different explanatory models
neg_log_likelihoods = zeros(size(subps));
lSubps = length(subps);
lstrSubps = num2str(length(num2str(length(subps))));
strTmp = repmat('\b',1,4+2*str2double(lstrSubps)+19);
fprintf(['explanatory model: %' lstrSubps 'd of %' lstrSubps 'd'],0,lSubps)
for j=1:lSubps
    fprintf([strTmp 'explanatory model: %' lstrSubps 'd of %' lstrSubps 'd'],...
        j,lSubps)
    % you may also use other covariance functions, please see function
    % GPREGRESSION_CALLER
    neg_log_likelihoods(j) = gpregression_caller(Y(:,subps{j}), diff_approx,...
        H, B, b, 'covMatern3iso', itermax);
end
fprintf('\n')
% add the negative log-marginal likelihood of the zero-order model
neg_log_likelihoods = [zom_nllh; neg_log_likelihoods]; 

% scaling, because of the finite wordlength
neg_log_likelihoods = neg_log_likelihoods+abs(min(neg_log_likelihoods));

% apply the bayesian rule for different explanatory models
ss_l = cellfun('size', subps, 1); % cardinalities of the subsets
priors = ones(length(ss_l)+zom,1);
priors = priors/sum(priors); % uniform prior
% calculate the models posterior probabilities with priors and negative
% log-likelihoods
probs = priors.*exp(-neg_log_likelihoods)/(sum(priors.*...
    exp(-neg_log_likelihoods)));

% sum up the interactions probabilities (discard the zero-order model)
probs = probs(1+zom:end);
[probs I] = sort(probs, 1, 'descend');
subps = subps(I);
result = zeros(length(indices_orig));
for subps_i=1:length(subps)
    mask = zeros(1,sum(indices_orig));
    mask(subps{subps_i}) = probs(subps_i);
    result(variable_i,indices_orig) =...
        result(variable_i,indices_orig) + mask;
end
end

function subsets = powerset(set, lub)
%POWERSET  A Matlab function to generate the powerset of a set
%
%   usage: subsets = powerset(set)
%
%   where:
%       set         is a vector which holds the members of a set
%       lub         maximum cardinality of a subset
%       subsets     is a cell vector which holds all subsets of a set
%       
%   Author: Tarmo Äijö <tarmo.aijo@tut.fi>

n_elements = 0;
for i=1:lub
    n_elements = n_elements + nchoosek(length(set), i);
end
subsets = cell(n_elements,1); % initialize the power set
card_set = length(set);
k = 1;
for i=1:lub
    subset_indices = nchoosek([1:card_set], i);
    l = 1;
    for j=k:k+size(subset_indices, 1) - 1
        subsets{j} = set([subset_indices(l,:)]);
        l = l + 1;
    end
    k = j + 1;
end
end

function diff_approx = approx_diff(x, t)
%APPROX_DIFF  Approximate the derivative
%       
%   Author: Tarmo Äijö <tarmo.aijo@tut.fi>

diff_approx = zeros(size(x));

for idx=1:length(x)-1
    diff_approx(idx) = (x(idx+1)-x(idx))/(t(idx+1)-t(idx));
end
end

function Y_trunc = trunc(Y, n_timepoints, delay)
%TRUNC  Truncate matrices based on the information in the variables 
%       n_timepoints and delay
%       
%   Author: Tarmo Äijö <tarmo.aijo@tut.fi>

indices = logical(ones(1, size(Y, 2)));
if delay < 0
    for idx=1:abs(delay)
        indices(idx) = 0;
    end
end
if length(n_timepoints) == 1
    for idx=n_timepoints:n_timepoints:size(Y,2)
        for delay_idx=0:abs(delay)
            if delay >= 0
                indices(idx-delay_idx) = 0;
            elseif delay <= 0
                indices(idx+delay_idx) = 0;
            end     
        end
    end
else
    c_n_tp = cumsum(n_timepoints);
        for idx=c_n_tp
            for delay_idx=0:abs(delay)
                if delay >= 0
                    indices(idx-delay_idx) = 0;
                elseif delay <= 0
                    indices(idx+delay_idx) = 0;
                end     
            end
        end
end
Y_trunc = Y(:,indices);
end

function nlml = gpregression_caller(x, y, H, B, b, covfunction, imax)
%GPREGRESSION_CALLER  A Matlab wrapper for the GPML regression code
%   The GPML Matlab code is written by Carl Edward Rasmussen.
%
%   usage: nlml = gpregression_caller(x, y, covfunction, imax)
%
%   where:
%       x           is a n by D matrix of training inputs
%       y           is a column vector (of size n) of targets
%       covfunction is the covariance function
%       imax        is the number of maximum iterations of the optimization
%       nlml        is the negative log marginal likelihood
%       
%   Possible values for the covfunction are
%       - 'covSEiso'
%       - 'covRQiso'
%       - 'covMatern3iso'
%       - 'covMatern5iso'
%       - 'covPeriodic'
%       - 'covLINone'
%       - 'covNNone'
%       - 'covConst'
%   
%   See also minimize, gpr
%
%   Author: Tarmo Äijö <tarmo.aijo@tut.fi>

covfunc = {'covSum', {covfunction,'covNoise'}};

if(strcmp(covfunction, 'covConst') || strcmp(covfunction,'covLINone'))
    loghyper = [-1; -1];
elseif(strcmp(covfunction,'covMatern3iso') || ...
       strcmp(covfunction, 'covMatern5iso') || ...
       strcmp(covfunction, 'covSEiso') || ...
       strcmp(covfunction, 'covPeriodic') || ...
       strcmp(covfunction, 'covNNone'))
    loghyper = [-1; -1; -1];
elseif(strcmp(covfunction, 'covRQiso'))
    loghyper = [-1; -1; -1; -1];
else
    error('Invalid covariance function');
end

loghyper = minimize(loghyper, 'gpr', -imax, covfunc, x, y, H, B, b);

nlml = gpr(loghyper, covfunc, x, y, H, B, b);

end

function [X, fX, i] = minimize(X, f, length, varargin)
% Minimize a differentiable multivariate function. 
%
% Usage: [X, fX, i] = minimize(X, f, length, P1, P2, P3, ... )
%
% where the starting point is given by "X" (D by 1), and the function named in
% the string "f", must return a function value and a vector of partial
% derivatives of f wrt X, the "length" gives the length of the run: if it is
% positive, it gives the maximum number of line searches, if negative its
% absolute gives the maximum allowed number of function evaluations. You can
% (optionally) give "length" a second component, which will indicate the
% reduction in function value to be expected in the first line-search (defaults
% to 1.0). The parameters P1, P2, P3, ... are passed on to the function f.
%
% The function returns when either its length is up, or if no further progress
% can be made (ie, we are at a (local) minimum, or so close that due to
% numerical problems, we cannot get any closer). NOTE: If the function
% terminates within a few iterations, it could be an indication that the
% function values and derivatives are not consistent (ie, there may be a bug in
% the implementation of your "f" function). The function returns the found
% solution "X", a vector of function values "fX" indicating the progress made
% and "i" the number of iterations (line searches or function evaluations,
% depending on the sign of "length") used.
%
% The Polack-Ribiere flavour of conjugate gradients is used to compute search
% directions, and a line search using quadratic and cubic polynomial
% approximations and the Wolfe-Powell stopping criteria is used together with
% the slope ratio method for guessing initial step sizes. Additionally a bunch
% of checks are made to make sure that exploration is taking place and that
% extrapolation will not be unboundedly large.
%
% See also: checkgrad 
%
% Copyright (C) 2001 - 2006 by Carl Edward Rasmussen (2006-09-08).

INT = 0.1;    % don't reevaluate within 0.1 of the limit of the current bracket
EXT = 3.0;                  % extrapolate maximum 3 times the current step-size
MAX = 20;                         % max 20 function evaluations per line search
RATIO = 10;                                       % maximum allowed slope ratio
SIG = 0.1; RHO = SIG/2; % SIG and RHO are the constants controlling the Wolfe-
% Powell conditions. SIG is the maximum allowed absolute ratio between
% previous and new slopes (derivatives in the search direction), thus setting
% SIG to low (positive) values forces higher precision in the line-searches.
% RHO is the minimum allowed fraction of the expected (from the slope at the
% initial point in the linesearch). Constants must satisfy 0 < RHO < SIG < 1.
% Tuning of SIG (depending on the nature of the function to be optimized) may
% speed up the minimization; it is probably not worth playing much with RHO.

% The code falls naturally into 3 parts, after the initial line search is
% started in the direction of steepest descent. 1) we first enter a while loop
% which uses point 1 (p1) and (p2) to compute an extrapolation (p3), until we
% have extrapolated far enough (Wolfe-Powell conditions). 2) if necessary, we
% enter the second loop which takes p2, p3 and p4 chooses the subinterval
% containing a (local) minimum, and interpolates it, unil an acceptable point
% is found (Wolfe-Powell conditions). Note, that points are always maintained
% in order p0 <= p1 <= p2 < p3 < p4. 3) compute a new search direction using
% conjugate gradients (Polack-Ribiere flavour), or revert to steepest if there
% was a problem in the previous line-search. Return the best value so far, if
% two consecutive line-searches fail, or whenever we run out of function
% evaluations or line-searches. During extrapolation, the "f" function may fail
% either with an error or returning Nan or Inf, and minimize should handle this
% gracefully.

if max(size(length)) == 2, red=length(2); length=length(1); else red=1; end
if length>0, S='Linesearch'; else S='Function evaluation'; end 

i = 0;                                            % zero the run length counter
ls_failed = 0;                             % no previous line search has failed
[f0 df0] = feval(f, X, varargin{:});          % get function value and gradient
fX = f0;
i = i + (length<0);                                            % count epochs?!
s = -df0; d0 = -s'*s;           % initial search direction (steepest) and slope
x3 = red/(1-d0);                                  % initial step is red/(|s|+1)

while i < abs(length)                                      % while not finished
  i = i + (length>0);                                      % count iterations?!

  X0 = X; F0 = f0; dF0 = df0;                   % make a copy of current values
  if length>0, M = MAX; else M = min(MAX, -length-i); end

  while 1                             % keep extrapolating as long as necessary
    x2 = 0; f2 = f0; d2 = d0; f3 = f0; df3 = df0;
    success = 0;
    while ~success && M > 0
      try
        M = M - 1; i = i + (length<0);                         % count epochs?!
        [f3 df3] = feval(f, X+x3*s, varargin{:});
        if isnan(f3) || isinf(f3) || any(isnan(df3)+isinf(df3)), error(''), end
        success = 1;
      catch                                % catch any error which occured in f
        x3 = (x2+x3)/2;                                  % bisect and try again
      end
    end
    if f3 < F0, X0 = X+x3*s; F0 = f3; dF0 = df3; end         % keep best values
    d3 = df3'*s;                                                    % new slope
    if d3 > SIG*d0 || f3 > f0+x3*RHO*d0 || M == 0  % are we done extrapolating?
      break
    end
    x1 = x2; f1 = f2; d1 = d2;                        % move point 2 to point 1
    x2 = x3; f2 = f3; d2 = d3;                        % move point 3 to point 2
    A = 6*(f1-f2)+3*(d2+d1)*(x2-x1);                 % make cubic extrapolation
    B = 3*(f2-f1)-(2*d1+d2)*(x2-x1);
    x3 = x1-d1*(x2-x1)^2/(B+sqrt(B*B-A*d1*(x2-x1))); % num. error possible, ok!
    if ~isreal(x3) || isnan(x3) || isinf(x3) || x3 < 0 % num prob | wrong sign?
      x3 = x2*EXT;                                 % extrapolate maximum amount
    elseif x3 > x2*EXT                  % new point beyond extrapolation limit?
      x3 = x2*EXT;                                 % extrapolate maximum amount
    elseif x3 < x2+INT*(x2-x1)         % new point too close to previous point?
      x3 = x2+INT*(x2-x1);
    end
  end                                                       % end extrapolation

  while (abs(d3) > -SIG*d0 || f3 > f0+x3*RHO*d0) && M > 0  % keep interpolating
    if d3 > 0 || f3 > f0+x3*RHO*d0                         % choose subinterval
      x4 = x3; f4 = f3; d4 = d3;                      % move point 3 to point 4
    else
      x2 = x3; f2 = f3; d2 = d3;                      % move point 3 to point 2
    end
    if f4 > f0           
      x3 = x2-(0.5*d2*(x4-x2)^2)/(f4-f2-d2*(x4-x2));  % quadratic interpolation
    else
      A = 6*(f2-f4)/(x4-x2)+3*(d4+d2);                    % cubic interpolation
      B = 3*(f4-f2)-(2*d2+d4)*(x4-x2);
      x3 = x2+(sqrt(B*B-A*d2*(x4-x2)^2)-B)/A;        % num. error possible, ok!
    end
    if isnan(x3) || isinf(x3)
      x3 = (x2+x4)/2;               % if we had a numerical problem then bisect
    end
    x3 = max(min(x3, x4-INT*(x4-x2)),x2+INT*(x4-x2));  % don't accept too close
    [f3 df3] = feval(f, X+x3*s, varargin{:});
    if f3 < F0, X0 = X+x3*s; F0 = f3; dF0 = df3; end         % keep best values
    M = M - 1; i = i + (length<0);                             % count epochs?!
    d3 = df3'*s;                                                    % new slope
  end                                                       % end interpolation

  if abs(d3) < -SIG*d0 && f3 < f0+x3*RHO*d0          % if line search succeeded
    X = X+x3*s; f0 = f3; fX = [fX' f0]';                     % update variables
    %fprintf('%s %6i;  Value %4.6e\r', S, i, f0); % commented by Tarmo Äijö
    s = (df3'*df3-df0'*df3)/(df0'*df0)*s - df3;   % Polack-Ribiere CG direction
    df0 = df3;                                               % swap derivatives
    d3 = d0; d0 = df0'*s;
    if d0 > 0                                      % new slope must be negative
      s = -df0; d0 = -s'*s;                  % otherwise use steepest direction
    end
    x3 = x3 * min(RATIO, d3/(d0-realmin));          % slope ratio but max RATIO
    ls_failed = 0;                              % this line search did not fail
  else
    X = X0; f0 = F0; df0 = dF0;                     % restore best point so far
    if ls_failed || i > abs(length)         % line search failed twice in a row
      break;                             % or we ran out of time, so we give up
    end
    s = -df0; d0 = -s'*s;                                        % try steepest
    x3 = 1/(1-d0);                     
    ls_failed = 1;                                    % this line search failed
  end
end
%fprintf('\n'); % commented by Tarmo Äijö 
end

function [out1, out2] = gpr(logtheta, covfunc, x, y, H, B, b, xstar, Hstar) % modified by Tarmo Äijö
% function [out1, out2] = gpr(logtheta, covfunc, x, y, xstar);
% gpr - Gaussian process regression, with a named covariance function. Two
% modes are possible: training and prediction: if no test data are given, the
% function returns minus the log likelihood and its partial derivatives with
% respect to the hyperparameters; this mode is used to fit the hyperparameters.
% If test data are given, then (marginal) Gaussian predictions are computed,
% whose mean and variance are returned. Note that in cases where the covariance
% function has noise contributions, the variance returned in S2 is for noisy
% test targets; if you want the variance of the noise-free latent function, you
% must substract the noise variance.
%
% usage: [nlml dnlml] = gpr(logtheta, covfunc, x, y)
%    or: [mu S2]  = gpr(logtheta, covfunc, x, y, xstar)
%
% where:
%
%   logtheta is a (column) vector of log hyperparameters
%   covfunc  is the covariance function
%   x        is a n by D matrix of training inputs
%   y        is a (column) vector (of size n) of targets
%   xstar    is a nn by D matrix of test inputs
%   nlml     is the returned value of the negative log marginal likelihood
%   dnlml    is a (column) vector of partial derivatives of the negative
%                 log marginal likelihood wrt each log hyperparameter
%   mu       is a (column) vector (of size nn) of prediced means
%   S2       is a (column) vector (of size nn) of predicted variances
%
% For more help on covariance functions, see "help covFunctions".
%
% (C) copyright 2006 by Carl Edward Rasmussen (2006-03-20).

if ischar(covfunc), covfunc = cellstr(covfunc); end % convert to cell if needed
[n, D] = size(x);
if eval(feval(covfunc{:})) ~= size(logtheta, 1)
  error('Error: Number of parameters do not agree with covariance function')
end

K_y = feval(covfunc{:}, logtheta, x);    % compute training set covariance matrix (modified by Tarmo Äijö)
K = K_y + H'*B*H; % added by Tarmo Äijö

L = chol(K)';                        % cholesky factorization of the covariance
alpha = solve_chol(L',(H'*b-y)); % modified  by Tarmo Äijö

if nargin == 7 % if no test cases, compute the negative log marginal likelihood (modified by Tarmo Äijö)
  out1 = 0.5*(H'*b-y)'*alpha + sum(log(diag(L))) + 0.5*n*log(2*pi); % modified by Tarmo Äijö

  if nargout == 2               % ... and if requested, its partial derivatives
    out2 = zeros(size(logtheta));       % set the size of the derivative vector
    W = L'\(L\eye(n))-alpha*alpha';                % precompute for convenience
    for i = 1:length(out2)
      out2(i) = sum(sum(W.*feval(covfunc{:}, logtheta, x, i)))/2;
    end
  end

else                    % ... otherwise compute (marginal) test predictions ...

  [Kss, Kstar] = feval(covfunc{:}, logtheta, x, xstar);     %  test covariances

  K_y_inv = inv(K_y); % added by Tarmo Äijö
  beta = inv(inv(B)+H*K_y_inv*H')*(H*K_y_inv*y+inv(B)*b); % added by Tarmo Äijö
  R = Hstar-H*K_y_inv*Kstar; % added by Tarmo Äijö
  out1 = Hstar'*beta+Kstar'*K_y_inv*(y-H'*beta); % predicted means (modified by Tarmo Äijö)
  
  if nargout == 2
    v = L\Kstar;
    out2 = Kss-sum(v.*v)'+diag(R'*inv(inv(B)+H*K_y_inv*H')*R); % modified by Tarmo Äijö
  end  

end
end

function C = sq_dist(a, b, Q)
% sq_dist - a function to compute a matrix of all pairwise squared distances
% between two sets of vectors, stored in the columns of the two matrices, a
% (of size D by n) and b (of size D by m). If only a single argument is given
% or the second matrix is empty, the missing matrix is taken to be identical
% to the first.
%
% Special functionality: If an optional third matrix argument Q is given, it
% must be of size n by m, and in this case a vector of the traces of the
% product of Q' and the coordinatewise squared distances is returned.
%
% NOTE: The program code is written in the C language for efficiency and is
% contained in the file sq_dist.c, and should be compiled using matlabs mex
% facility. However, this file also contains a (less efficient) matlab
% implementation, supplied only as a help to people unfamiliar with mex. If
% the C code has been properly compiled and is avaiable, it automatically
% takes precendence over the matlab code in this file.
%
% Usage: C = sq_dist(a, b)
%    or: C = sq_dist(a)  or equiv.: C = sq_dist(a, [])
%    or: c = sq_dist(a, b, Q)
% where the b matrix may be empty.
%
% where a is of size D by n, b is of size D by m (or empty), C and Q are of
% size n by m and c is of size D by 1.
%
% Copyright (c) 2003, 2004, 2005 and 2006 Carl Edward Rasmussen.
% 2006-03-09.

if nargin < 1 | nargin > 3 | nargout > 1
  error('Wrong number of arguments.');
end

if nargin == 1 | isempty(b)                   % input arguments are taken to be
  b = a;                                   % identical if b is missing or empty
end 

[D, n] = size(a); 
[d, m] = size(b);
if d ~= D
  error('Error: column lengths must agree.');
end

if nargin < 3
  C = zeros(n,m);
  for d = 1:D
    C = C + (repmat(b(d,:), n, 1) - repmat(a(d,:)', 1, m)).^2;
  end
  % C = repmat(sum(a.*a)',1,m)+repmat(sum(b.*b),n,1)-2*a'*b could be used to 
  % replace the 3 lines above; it would be faster, but numerically less stable.
else
  if [n m] == size(Q)
    C = zeros(D,1);
    for d = 1:D
      C(d) = sum(sum((repmat(b(d,:), n, 1) - repmat(a(d,:)', 1, m)).^2.*Q));
    end
  else
    error('Third argument has wrong size.');
  end
end
end

function x = solve_chol(A, B)
% solve_chol - solve linear equations from the Cholesky factorization.
% Solve A*X = B for X, where A is square, symmetric, positive definite. The
% input to the function is R the Cholesky decomposition of A and the matrix B.
% Example: X = solve_chol(chol(A),B);
%
% NOTE: The program code is written in the C language for efficiency and is
% contained in the file solve_chol.c, and should be compiled using matlabs mex
% facility. However, this file also contains a (less efficient) matlab
% implementation, supplied only as a help to people unfamiliar with mex. If
% the C code has been properly compiled and is avaiable, it automatically
% takes precendence over the matlab code in this file.
%
% Copyright (c) 2004, 2005, 2006 by Carl Edward Rasmussen. 2006-02-08.

if nargin ~= 2 | nargout > 1
  error('Wrong number of arguments.');
end

if size(A,1) ~= size(A,2) | size(A,1) ~= size(B,1)
  error('Wrong sizes of matrix arguments.');
end

x = A\(A'\B);
end

function [A, B] = covSum(covfunc, logtheta, x, z)
% covSum - compose a covariance function as the sum of other covariance
% functions. This function doesn't actually compute very much on its own, it
% merely does some bookkeeping, and calls other covariance functions to do the
% actual work.
%
% For more help on design of covariance functions, try "help covFunctions".
%
% (C) Copyright 2006 by Carl Edward Rasmussen, 2006-03-20.

for i = 1:length(covfunc)                   % iterate over covariance functions
  f = covfunc(i);
  if iscell(f{:}), f = f{:}; end          % dereference cell array if necessary
  j(i) = cellstr(feval(f{:}));
end

if nargin == 1,                                   % report number of parameters
  A = char(j(1)); for i=2:length(covfunc), A = [A, '+', char(j(i))]; end
  return
end

[n, D] = size(x);

v = [];              % v vector indicates to which covariance parameters belong
for i = 1:length(covfunc), v = [v repmat(i, 1, eval(char(j(i))))]; end

switch nargin
case 3                                              % compute covariance matrix
  A = zeros(n, n);                       % allocate space for covariance matrix
  for i = 1:length(covfunc)                  % iteration over summand functions
    f = covfunc(i);
    if iscell(f{:}), f = f{:}; end        % dereference cell array if necessary
    A = A + feval(f{:}, logtheta(v==i), x);            % accumulate covariances
  end

case 4                      % compute derivative matrix or test set covariances
  if nargout == 2                                % compute test set cavariances
    A = zeros(size(z,1),1); B = zeros(size(x,1),size(z,1));    % allocate space
    for i = 1:length(covfunc)
      f = covfunc(i);
      if iscell(f{:}), f = f{:}; end      % dereference cell array if necessary
      [AA BB] = feval(f{:}, logtheta(v==i), x, z);   % compute test covariances
      A = A + AA; B = B + BB;                                  % and accumulate
    end
  else                                            % compute derivative matrices
    i = v(z);                                       % which covariance function
    j = sum(v(1:z)==i);                    % which parameter in that covariance
    f = covfunc(i);
    if iscell(f{:}), f = f{:}; end        % dereference cell array if necessary
    A = feval(f{:}, logtheta(v==i), x, j);                 % compute derivative
  end

end
end

function [A, B] = covSEiso(loghyper, x, z)
% Squared Exponential covariance function with isotropic distance measure. The 
% covariance function is parameterized as:
%
% k(x^p,x^q) = sf2 * exp(-(x^p - x^q)'*inv(P)*(x^p - x^q)/2) 
%
% where the P matrix is ell^2 times the unit matrix and sf2 is the signal
% variance. The hyperparameters are:
%
% loghyper = [ log(ell)
%              log(sqrt(sf2)) ]
%
% For more help on design of covariance functions, try "help covFunctions".
%
% (C) Copyright 2006 by Carl Edward Rasmussen (2007-06-25)

if nargin == 0, A = '2'; return; end              % report number of parameters

[n D] = size(x);
ell = exp(loghyper(1));                           % characteristic length scale
sf2 = exp(2*loghyper(2));                                     % signal variance

if nargin == 2
  A = sf2*exp(-sq_dist(x'/ell)/2);
elseif nargout == 2                              % compute test set covariances
  A = sf2*ones(size(z,1),1);
  B = sf2*exp(-sq_dist(x'/ell,z'/ell)/2);
else                                                % compute derivative matrix
  if z == 1                                                   % first parameter
    A = sf2*exp(-sq_dist(x'/ell)/2).*sq_dist(x'/ell);  
  else                                                       % second parameter
    A = 2*sf2*exp(-sq_dist(x'/ell)/2);
  end
end
end

function [A, B] = covNoise(logtheta, x, z)

% Independent covariance function, ie "white noise", with specified variance.
% The covariance function is specified as:
%
% k(x^p,x^q) = s2 * \delta(p,q)
%
% where s2 is the noise variance and \delta(p,q) is a Kronecker delta function
% which is 1 iff p=q and zero otherwise. The hyperparameter is
%
% logtheta = [ log(sqrt(s2)) ]
%
% For more help on design of covariance functions, try "help covFunctions".
%
% (C) Copyright 2006 by Carl Edward Rasmussen, 2006-03-24.

if nargin == 0, A = '1'; return; end              % report number of parameters

s2 = exp(2*logtheta);                                          % noise variance

if nargin == 2                                      % compute covariance matrix
  A = s2*eye(size(x,1));
elseif nargout == 2                              % compute test set covariances
  A = s2;
  B = 0;                               % zeros cross covariance by independence
else                                                % compute derivative matrix
  A = 2*s2*eye(size(x,1));
end
end

function [A, B] = covMatern3iso(loghyper, x, z)
% Matern covariance function with nu = 3/2 and isotropic distance measure. The
% covariance function is:
%
% k(x^p,x^q) = s2f * (1 + sqrt(3)*d(x^p,x^q)) * exp(-sqrt(3)*d(x^p,x^q))
%
% where d(x^p,x^q) is the distance sqrt((x^p-x^q)'*inv(P)*(x^p-x^q)), P is ell
% times the unit matrix and sf2 is the signal variance. The hyperparameters
% are:
%
% loghyper = [ log(ell)
%              log(sqrt(sf2)) ]
%
% For more help on design of covariance functions, try "help covFunctions".
%
% (C) Copyright 2006 by Carl Edward Rasmussen (2006-03-24)

if nargin == 0, A = '2'; return; end

persistent K;
[n, D] = size(x);
ell = exp(loghyper(1));
sf2 = exp(2*loghyper(2));

x = sqrt(3)*x/ell;

if nargin == 2                                      % compute covariance matrix
  A = sqrt(sq_dist(x'));
  K = sf2*exp(-A).*(1+A);
  A = K;
elseif nargout == 2                              % compute test set covariances
  z = sqrt(3)*z/ell;
  A = sf2;
  B = sqrt(sq_dist(x',z'));
  B = sf2*exp(-B).*(1+B);
else                                              % compute derivative matrices
  if z == 1
    A = sf2*sq_dist(x').*exp(-sqrt(sq_dist(x')));
  else
    % check for correct dimension of the previously calculated kernel matrix
    if any(size(K)~=n)  
      K = sqrt(sq_dist(x'));
      K = sf2*exp(-K).*(1+K);
    end
    A = 2*K;
    clear K;
  end
end
end

function [A, B] = covMatern5iso(loghyper, x, z)
% Matern covariance function with nu = 5/2 and isotropic distance measure. The
% covariance function is:
%
% k(x^p,x^q) = s2f * (1 + sqrt(5)*d + 5*d/3) * exp(-sqrt(5)*d)
%
% where d is the distance sqrt((x^p-x^q)'*inv(P)*(x^p-x^q)), P is ell times
% the unit matrix and sf2 is the signal variance. The hyperparameters are:
%
% loghyper = [ log(ell)
%              log(sqrt(sf2)) ]
%
% For more help on design of covariance functions, try "help covFunctions".
%
% (C) Copyright 2006 by Carl Edward Rasmussen (2006-03-24)

if nargin == 0, A = '2'; return; end

persistent K;
[n, D] = size(x);
ell = exp(loghyper(1));
sf2 = exp(2*loghyper(2));

x = sqrt(5)*x/ell;

if nargin == 2                                      % compute covariance matrix
  A = sq_dist(x');
  K = sf2*exp(-sqrt(A)).*(1+sqrt(A)+A/3);
  A = K;
elseif nargout == 2                              % compute test set covariances
  z = sqrt(5)*z/ell;
  A = sf2;
  B = sq_dist(x',z');
  B = sf2*exp(-sqrt(B)).*(1+sqrt(B)+B/3);
else                                              % compute derivative matrices
  if z == 1
    A = sq_dist(x');
    A = sf2*(A+sqrt(A).^3).*exp(-sqrt(A))/3;
  else
    % check for correct dimension of the previously calculated kernel matrix
    if any(size(K)~=n)  
      K = sq_dist(x');
      K = sf2*exp(-sqrt(K)).*(1+sqrt(K)+K/3);
    end
    A = 2*K;
    clear K;
  end
end
end

function [A, B] = covNNone(loghyper, x, z)
% Neural network covariance function with a single parameter for the distance
% measure. The covariance function is parameterized as:
%
% k(x^p,x^q) = sf2 * asin(x^p'*P*x^q / sqrt[(1+x^p'*P*x^p)*(1+x^q'*P*x^q)])
%
% where the x^p and x^q vectors on the right hand side have an added extra bias
% entry with unit value. P is ell^-2 times the unit matrix and sf2 controls the
% signal variance. The hyperparameters are:
%
% loghyper = [ log(ell)
%              log(sqrt(sf2) ]
%
% For more help on design of covariance functions, try "help covFunctions".
%
% (C) Copyright 2006 by Carl Edward Rasmussen (2006-03-24)

if nargin == 0, A = '2'; return; end              % report number of parameters

persistent Q K;                 
[n D] = size(x);
ell = exp(loghyper(1)); em2 = ell^(-2);
sf2 = exp(2*loghyper(2));
x = x/ell;

if nargin == 2                                             % compute covariance
  Q = x*x';
  K = (em2+Q)./(sqrt(1+em2+diag(Q))*sqrt(1+em2+diag(Q)'));
  A = sf2*asin(K);                 
elseif nargout == 2                              % compute test set covariances
  z = z/ell; 
  A = sf2*asin((em2+sum(z.*z,2))./(1+em2+sum(z.*z,2)));
  B = sf2*asin((em2+x*z')./sqrt((1+em2+sum(x.*x,2))*(1+em2+sum(z.*z,2)')));
else                                                % compute derivative matrix
  % check for correct dimension of the previously calculated kernel matrix
  if any(size(Q)~=n)  
    Q = x*x';
  end
  % check for correct dimension of the previously calculated kernel matrix
  if any(size(K)~=n)  
    K = (em2+Q)./(sqrt(1+em2+diag(Q))*sqrt(1+em2+diag(Q)'));
  end
  if z == 1                                                   % first parameter
    v = (em2+sum(x.*x,2))./(1+em2+diag(Q));
    A = -2*sf2*((em2+Q)./(sqrt(1+em2+diag(Q))*sqrt(1+em2+diag(Q)'))- ...
                            K.*(repmat(v,1,n)+repmat(v',n,1))/2)./sqrt(1-K.^2);
    clear Q;
  else                                                       % second parameter
    A = 2*sf2*asin(K);
    clear K;
  end
end
end

function [A, B] = covRQiso(loghyper, x, z)
% Rational Quadratic covariance function with isotropic distance measure. The
% covariance function is parameterized as:
%
% k(x^p,x^q) = sf2 * [1 + (x^p - x^q)'*inv(P)*(x^p - x^q)/(2*alpha)]^(-alpha)
%
% where the P matrix is ell^2 times the unit matrix, sf2 is the signal
% variance and alpha is the shape parameter for the RQ covariance. The
% hyperparameters are:
%
% loghyper = [ log(ell)
%              log(sqrt(sf2))
%              log(alpha) ]
%
% For more help on design of covariance functions, try "help covFunctions".
%
% (C) Copyright 2006 by Carl Edward Rasmussen (2006-09-08)

if nargin == 0, A = '3'; return; end

[n, D] = size(x);

persistent K;
ell = exp(loghyper(1));
sf2 = exp(2*loghyper(2));
alpha = exp(loghyper(3));

if nargin == 2                                      % compute covariance matrix
  K = (1+0.5*sq_dist(x'/ell)/alpha);
  A = sf2*(K.^(-alpha));
elseif nargout == 2                              % compute test set covariances
  A = sf2*ones(size(z,1),1);
  B = sf2*((1+0.5*sq_dist(x'/ell,z'/ell)/alpha).^(-alpha));
else                                              % compute derivative matrices
  % check for correct dimension of the previously calculated kernel matrix
  if any(size(K)~=n)  
    K = (1+0.5*sq_dist(x'/ell)/alpha);
  end
  if z == 1                                           % length scale parameters
    A = sf2*K.^(-alpha-1).*sq_dist(x'/ell);
  elseif z == 2                                           % magnitude parameter
    A = 2*sf2*(K.^(-alpha));
  else
    A = sf2*K.^(-alpha).*(0.5*sq_dist(x'/ell)./K - alpha*log(K));
    clear K;
  end
end
end