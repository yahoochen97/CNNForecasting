function [allRaces,fts,s2s] = gpm(hyperparameter, xs, ys, raceinfos, plot_path, parms)

%  Obtain posterior belief of voter preferences on election day using gp model
%  input:
%    - hyperparameter: hyperparameter of mean/cov/lik
%    - xs: cell array of polling data in form of [tau, polling proportion, sample size]
%    - ys: cell array of [polling proportion]
%    - raceinfos: election metadata [year, state, candidatename, actual vote share, pvi, experienced, isRep]
%    - plot_path: path of plots if parms.plot==1
%    - parms: function keyword arguments
%       - tau: forecasting horizon
%       - j: current index of searching sequence
%       - type: prior model, 'GP' or 'LM'
%       - plot:
%           - 1: generate plots of underlying voter preference for each election race
%           - 0: just obtain posterior belief of voter preference on election day
%       - test_year: validation year in LOYO process or test year in testing process
%       - coefs: the coefs of prior linear model for the linear trend intercept
%
% Caller of this function should specify the keyword arguments parms.
%

    % define model
    [meanfunc, covfunc, likfunc, inffunc, prior] = model();
    mu_ml = prior.slope(1);
    sigma_mc = prior.intercept(2);
    
    
    % AllRaces is a struct with key year+state.
    % Value of allRaces struct is an array of [model posterior mean, actual
    % vote, ...] of all candidates for the year/state race
    allRaces = struct;
    n = numel(xs);
    
    % fts/s2s are the array of posterior mean/var on election day of all races
    fts = zeros(n,1);
    s2s = zeros(n,1);
    
    % iterate every race
    for i = 1:n
        % obtain metadata
        year = raceinfos{i}{1};
        state = raceinfos{i}{2}{1};
        candidateName = raceinfos{i}{3};
        trueVote = raceinfos{i}{4};
        pvi = raceinfos{i}{5};
        experienced = raceinfos{i}{6};
        party = raceinfos{i}{7};
        
        % key for allRaces struct
        fn = char(state+""+year);
        fn = fn(~isspace(fn));
        
        % compute prior mean on intercept
        mu_b = computePrior(pvi, experienced, party, parms);
        
        % heuristic 1:  (distance to 0.5) / 2
        % sigma_mc = abs(mu_b-0.5)/2;
        % sigma_mc = 0.05;
    
        if numel(xs{i})==0
            % if there is no data avaiable
            % use prior
            if ~isfield(allRaces, fn)
                allRaces.(fn) = [mu_b, trueVote];
            else
                allRaces.(fn) = [allRaces.(fn), mu_b , trueVote];
            end
            fts(i) = mu_b;
            s2s(i) = sigma_mc^2 + exp(2*hyperparameter.cov(2));
        else
            % if there is data avaiable   
            % define hyp for one specific race
            % hyp is exactly the same of hyperparameter except potentially
            % different prior mean of linear trend intercept
            hyp.mean = [mu_b];
            hyp.cov = hyperparameter.cov;
            hyp.lik = hyperparameter.lik;
            hyp.cov(4) = log(sigma_mc);
            
            % obtain gp posterior
            if parms.plot==0
                % obtain gp posterior on election day only
                % test position is [day 0, 0 polling porportion, 1 samplesize]
                xstar = [0,0,1];
                [~, ~, fts(i), s2s(i)] = gp(hyp, inffunc, meanfunc, covfunc, likfunc, xs{i}, ys{i}, xstar);
            else 
                % plot upon required
                % define 200 test location from the earliest poll day to
                % election day
                nz = 200;
                
                XMIN = floor(xs{i}(1,1)/parms.BIN)*parms.BIN;
                xstar = [linspace(XMIN,0,nz).',zeros(1,nz)',ones(1,nz)'];
                [~, ~, fmu, fs2] = gp(hyp, inffunc, meanfunc, covfunc, likfunc, xs{i}, ys{i}, xstar);
                fts(i) = fmu(end);
                s2s(i) = s2s(end);
                parms.prior = [mu_b, sigma_mc];
                fig = plot_posterior(fmu, fs2, xs{i}(:,1), ys{i}, xstar(:,1), parms);
                plot_title = year + " " + state + " " + candidateName;
%                 title(plot_title);
                
                % save plot to files
                yearFolder = fullfile(plot_path, num2str(year));
                stateFolder = fullfile(yearFolder, state);
                if ~exist(plot_path, 'dir')
                    mkdir(plot_path);
                end
                if ~exist(yearFolder, 'dir')
                    mkdir(yearFolder);
                end
                if ~exist(stateFolder, 'dir')
                    mkdir(stateFolder);
                end
                filename = fullfile(stateFolder, plot_title + num2str(parms.j)+".pdf");
                set(fig, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
                set(fig, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
                print(fig, filename, '-dpdf','-r300');
                close;
            end   
             
            % update allRaces struct
            if ~isfield(allRaces, fn)
                allRaces.(fn) = [fts(i), trueVote];
            else
                allRaces.(fn) = [allRaces.(fn), fts(i), trueVote];
            end
        end
    end
end