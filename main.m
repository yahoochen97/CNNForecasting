function [varout]=main(TYPE, mode, tau, plot)

%  main function of obtaining gp posteriors for each election race
%  input:
%    - TYPE: 'GP' or 'LM', specifying the prior model
%    - mode: 
%      -: 1 for doing LOYO on year 1992 to 2016
%      -: 2 for testing 2018 races or 2016 in R&R
%      -: 3 for forecasting 2020 races
%    - tau: forecasting horizon
%    - plot: whether to plot gp
%

    % define output
    varout = "";
    
    % load gmpl, utilities and data
    addpath("gpml-matlab-v3.6-2015-07-07");
    addpath("utilities");
    addpath("data");
    startup;

    % define horizons
    % 8/6/4/3/2/1/0 weeks
    taus = [0,7, 14,21, 28,42, 56];

    % define search space
    search_size = 100;
    if strcmp(TYPE, "GP")==1
        p = sobolset(3, 'Skip', 1);
    else 
        search_size = 20;
        p = linspace(0,0.1,search_size);
    end
    
    if mode>=2
        % no search needed
        % forecasting year 2018 or 2020 or 2016 in R&R
        ts = readData("results/"+TYPE+"_opthyp.csv");
        ts = ts.opt_idx;

        for i=1:numel(taus)
            j = ts(i);
            if strcmp(TYPE, "GP")==1
                ls = p(j,1)*(56-7)+7; % 7-56
                os = p(j,2)/20; % 0%-5%
                lik = p(j,3)/20; % 0%-5%
                myrun(taus(i),TYPE, ls, os, lik, j, mode, plot);
            else 
                % linear model does not have ls/os
                ls = 0; os = 0;
                lik = p(j); % 0%-10%
                myrun(taus(i),TYPE, ls, os, lik, j, mode, plot);
            end
        end
    else
        % performing a sobel sequence search
%         for i=1:numel(taus)
            for j=1:search_size
                if strcmp(TYPE, "GP")==1
                    ls = p(j,1)*(56-7)+7; % 7-56
                    os = p(j,2)/20; % 0%-5%
                    lik = p(j,3)/20; % 0%-5%
                    myrun(tau,TYPE, ls, os, lik, j, mode, plot);
                else 
                    % linear model does not have ls/os
                    ls = 0; os = 0;
                    lik = p(j); % 0%-10%
                    myrun(tau,TYPE, ls, os, lik, j, mode, plot);
                end
            end
%         end
    end
end

function myrun(tau,type, ls, os, lik, j, mode, plot)
%  load 1992-2016 data
%  feature of data includes: 
%       - cycle: election year
%       - state: state that holds election in the election year
%       - Candidateidentier: include election year, state and name of candidate
%       - Percentage_of_Vote_won_x: actual vote share, 0 for future elections
%       - samplesize: number of respondences in the poll
%       - daysLeft: days prior to election when the poll was conducted
%       - numberSupport: number of respondences that supports the candidate
%       - pollster: name of the polling company/institution
%       - Republican: 1 if the candidate is republican, 0 otherwise
%       - Democrat: 1 if the candidate is democratic, 0 otherwise
%       - pvi: Cook partisan voting index of the election state
%       - experienced: 1 if the candidate has served in any political office, 0 otherwise

    CNNdata = readData("data/CNNdata1992-2016.csv");
    
    if mode==2
        % test on 2018 data
        % load 2018 data
        
        % 2016 in R&R
        % do not need to load 2018/2020 data

        CNNdata2018 = readData("data/CNNData2018.csv");
        CNNdata2018(:, ["candidate_name"]) = [];
        CNNdata = vertcat(CNNdata, CNNdata2018);
    % elseif mode==3
        % forecast 2020 races
        CNNdata2020 = readData("data/CNNdata2020.csv");
        CNNdata2020(:, ["candidate_name"]) = [];
        result2020 = readData("data/2020results.csv");
        for i=1:size(CNNdata2020)
            v = result2020.Percentage_of_Vote_won_x(strcmp(result2020.Candidateidentifier,CNNdata2020.Candidateidentifier(i)));
            CNNdata2020.Percentage_of_Vote_won_x(i) = v;
        end
        
        % CNNdata2020.Percentage_of_Vote_won_x = zeros(size(CNNdata2020,1),1);
        CNNdata = vertcat(CNNdata, CNNdata2020);
    end

    % split data into cell arrays
    % obtain features (days before election) and regressor (polling proportions)
    years = unique(CNNdata.cycle);
    states = unique(CNNdata.state);
    [xs, ys, raceinfos] = buildTrainCellArrays(CNNdata, years, states);

    % shrink data according to the forecasting horizons
    counter = size(xs,1);
    for i=1:counter
        idx = xs{i}(:,1) <= -tau;
        xs{i} = xs{i}(idx,:);
        ys{i} = ys{i}(idx);
    end

    % define hyp (hyperparameter) struct
    %   - hyp.mean:
    %       - prior mean of the linear trend slope
    %       - prior mean of the linear trend intercept
    %   - hyp.cov:
    %       - log of length scale in Matérn kernel
    %       - log of outpur scale in Matérn kernel
    %       - 1 / log of prior std of the linear trend slope
    %       - log of prior std of the linear trend intercept
    %   - hyp.lik: log of observation noise std
    
    hyp.cov(1) = log(ls);
    hyp.cov(2) = log(os);
    
    % model() defines prior distribution of linear trends
    [~,~,~,~, prior] = model();
    sigma_ml = prior.slope(2);
    sigma_mc = prior.intercept(2);
    hyp.cov(3) = log(1/sigma_ml);
    hyp.cov(4) = log(sigma_mc);
    hyp.lik = log(lik);

    % parms struct specifies keyword arguments
    %   - tau: forecasting horizon
    %   - j: current index of searching sequence
    %   - type: prior model, 'GP' or 'LM'
    %   - plot:
    %       - 1: generate plots of underlying voter preference for each election race
    %       - 0: just obtain posterior belief of voter preference on election day
    %   - test_year: validation year in LOYO process or test year in testing process
    %   - coefs: the coefs of prior linear model for the linear trend intercept
    parms.tau = tau;
    parms.j = j;
    parms.type = type;
    parms.plot = plot;
    
    % plot days bin
    parms.BIN = 30;
    
    if mode==2
        % tesing 2018 races
        parms.test_year = 2018;
        
        % tesing 2016 races in R&R
        % parms.test_year = 2016;
        winner_all = readData("data/winners_all.csv");
        parms.test_year = 2020;
        disp(parms.test_year);

        % add incumbency
        for i=1:numel(raceinfos)
           year = raceinfos{i}{1};
           state = raceinfos{i}{2};
           c = cell2mat(raceinfos{i}{3});
           incumbent = cell2mat(winner_all.incumbency(winner_all.year==(year-6) & strcmp(winner_all.state,state)));
           if strcmp(c,incumbent)
              raceinfos{i}{8} = 1; 
           else
              raceinfos{i}{8} = 0;
           end
        end
        
        % precompute coefs of prior linear model of the linear trend intercept
        parms.coefs = priorModel(CNNdata, parms.test_year);
        plot_path = "plots/" + type + "MargLinTre"+num2str(parms.test_year)+"_"+num2str(tau);
        plot_path = "RR/plots/" + type + "_"+num2str(parms.test_year)+"_"+num2str(tau);
        if strcmp(type, "GP")==1
            [allRaces, fts, s2s] = gpm(hyp, xs, ys, raceinfos, plot_path, parms);
        else
            [allRaces, fts, s2s] = lm(hyp, xs, ys, raceinfos, plot_path, parms);
        end
        % post train process
        % compute corr, rmse, accuracy, coverage rate and nlz
        % write results to csv files
        posttrain(raceinfos,fts,s2s,allRaces,hyp, tau, parms);   
    elseif mode==3
        % forecasting 2020 races
        parms.test_year = 2020;
        % precompute coefs of prior linear model of the linear trend intercept
        parms.coefs = priorModel(CNNdata, parms.test_year);
        plot_path = "plots/" + type + "MargLinTre"+num2str(parms.test_year)+"_"+num2str(tau);
        if strcmp(type, "GP")==1
            [allRaces, fts, s2s] = gpm(hyp, xs, ys, raceinfos, plot_path, parms);
        else
            [allRaces, fts, s2s] = lm(hyp, xs, ys, raceinfos, plot_path, parms);
        end
        % post train process
        % compute corr, rmse, accuracy, coverage rate and nlz
        % write results to csv files
        posttrain(raceinfos,fts,s2s,allRaces,hyp, tau, parms);   
    else
        % Leave-one-year-out process
        % iterate over 1992-2016
        for y=1:numel(years)
            parms.test_year = years(y);
            disp(parms.test_year);
            % precompute coefs of prior linear model of the linear trend intercept
            parms.coefs = priorModel(CNNdata, parms.test_year);
            plot_path = "plots/" + type + "MargLinTre"+num2str(parms.test_year)+"_"+num2str(tau);
            if strcmp(type, "GP")==1
                [allRaces, fts, s2s] = gpm(hyp, xs, ys, raceinfos, plot_path, parms);
            else
                [allRaces,fts,s2s] = lm(hyp, xs, ys, raceinfos, plot_path, parms);
            end        
            % post train process
            % compute corr, rmse, accuracy, coverage rate and nlz
            % write results to csv files
            posttrain(raceinfos,fts,s2s,allRaces,hyp, tau, parms);
        end
    end
end
