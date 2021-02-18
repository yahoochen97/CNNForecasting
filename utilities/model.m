function [meanfunc, covfunc, likfunc, inffunc, prior] = model()

%
%   Define the mean/cov/lik/inf function of gp models, and the prior distribution
%   

    % feature vector includes days before election, polling porportion and sample size
    % define mask over feature vector
    % mask of days before election
    meanmask = [true, false, false];
    % mask of polling porportion and sample size
    covmask = [false, true, true];
    
    % define masked mean function of linear trends
    %   - mmc: masled linear intercept mean
    mc = {@meanConst};
    mmc = {@meanMask, meanmask, mc};
    meanfunc = mmc;
    
    % define masked covariance function of non-linear trends
    %    - cmd: masked diagonal in-poll covariance
    %    - cmm: masked Matérn covariance
    %    - cms: masked covariance of linear trend slope
    %    - ci: constant covariance of linear trend intercept
    cm = {@covMaterniso, 3};
    mb = {@logsqrtbinom};
    cdb = {@covDiag, mb};
    cmd = {@covMask, {covmask, cdb}};
    cmm = {@covMask, {meanmask, cm}};
    cs = {@covLINiso};
    ci = {@covConst};
    cms = {@covMask, {meanmask, cs}};
    covfunc = {@covSum, {cmm, cmd, cms, ci}};
    
    % define Gaussian likelihood and inference function
    likfunc = @likGauss;
    inffunc = @infExact;

    % define prior distribution on linear trend
    mu_ml = 0.0; sigma_ml = 0.002;
    
    % try difference heuristic for intercept std
    % heuristic 1: (distance to 0.5) / 2
    % heuristic 2: 5%
    mu_mc = 0.5; sigma_mc = 0.1;
    
    % heuristic 2
    % sigma_mc = 0.05;
    prior.slope = [mu_ml, sigma_ml];
    prior.intercept = [mu_mc, sigma_mc];
end
