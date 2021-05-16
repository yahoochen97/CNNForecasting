% load everything
addpath("gpml-matlab-v3.6-2015-07-07");
addpath("utilities");
addpath("data");
startup;

CNNdata = readData("data/CNNdata1992-2016.csv");
CNNdata2018 = readData("data/CNNData2018.csv");
CNNdata2018(:, ["candidate_name"]) = [];
CNNdata = vertcat(CNNdata, CNNdata2018);

CNNdata2020 = readData("data/CNNData2020.csv");
CNNdata2020(:, ["candidate_name"]) = [];
result2020 = readData("data/2020results.csv");
for i=1:size(CNNdata2020)
    v = result2020.Percentage_of_Vote_won_x(strcmp(result2020.Candidateidentifier,CNNdata2020.Candidateidentifier(i)));
    CNNdata2020.Percentage_of_Vote_won_x(i) = v;
end

CNNdata = vertcat(CNNdata, CNNdata2020);


% % split data into cell arrays
% % obtain features (days before election) and regressor (polling proportions)
years = unique(CNNdata.cycle);
states = unique(CNNdata.state);
[xs, ys, raceinfos, ~] = buildTrainCellArrays(CNNdata, years, states);

% % number of polls versus final margin
y = [];
s = cell(numel(raceinfos),1);
numpoll = [];
v =cell(numel(raceinfos),1);

for i=1:numel(raceinfos)
    year = raceinfos{i}{1};
state = raceinfos{i}{2}{1};
trueVote = raceinfos{i}{4}/100;
numpoll = [numpoll;size(xs{i},1)];
v{i}=  trueVote;
y =  [y;year];
s{i}=  state;
end
result = table(y,s,numpoll,v);
writetable(result,strcat('RR/results/numpoll.csv'));


js = [];
for j=1:numel(states)
    if states{j}(end)=="S"
        js = [js, j];
    end
end
states(js)=[]; 

a = [];
b = [];
for i=1:numel(years)
    data = readData("results/stan_LOOGP_" + int2str(years(i)) + "day0_37.csv");
    for j=1:numel(states)
        if numel(data(strcmp(data.state, states{j}),:))==0, continue; end
       [m, idx] = max(data(strcmp(data.state, states{j}),:).median);
       v = data(strcmp(data.state, states{j}),:).vote(idx);
        b = [b,m-v];
        a = [a, max(result(result.y==years(i) & strcmp(result.s, states{j}),:).numpoll)];
    end
end

c = [];
d = [];
for i=1:numel(years)
    data = readData("results/stan_LOOGP_" + int2str(years(i)) + "day0_37.csv");
    for j=1:numel(states)
        if numel(data(strcmp(data.state, states{j}),:))==0, continue; end
        c = [c; data(strcmp(data.state, states{j}),:).median];
        d = [d; data(strcmp(data.state, states{j}),:).vote];
    end
end


fig = figure(1);
scatter(a,abs(b),8,"filled");
xlabel("# of polls"); ylabel("Final margin");
filename = fullfile("/Users/yahoo/Desktop/pollvsmargin.pdf");
set(fig, 'PaperPosition', [0 0 10 10]); %Position plot at left hand corner with width 5 and height 5.
set(fig, 'PaperSize', [10 10]); %Set the paper to have width 5 and height 5.
print(fig, filename, '-dpdf','-r300'); 

exit;

% % different prior heuristics for intercept std
p = sobolset(3, 'Skip', 1);
ts = readData("results/GP_opthyp.csv");
ts = ts.opt_idx;
j = ts(1);

ls = p(j,1)*(56-7)+7; % 7-56
os = p(j,2)/20; % 0%-5%
lik = p(j,3)/20; % 0%-5%

hyp.cov(1) = log(ls);
hyp.cov(2) = log(os);

% model() defines prior distribution of linear trends
[~,~,~,~, prior] = model();
sigma_ml = prior.slope(2);
sigma_mc = prior.intercept(2);
hyp.cov(3) = log(1/sigma_ml);
hyp.cov(4) = log(sigma_mc);
hyp.lik = log(lik);

parms.tau = 0;
parms.j = j;
parms.type = 'GP';
parms.plot = 1;

% plot days bin
parms.BIN = 30;

parms.test_year = 2016;

% % precompute coefs of prior linear model of the linear trend intercept
parms.coefs = priorModel(CNNdata, parms.test_year);
plot_path = "RR/plots/" + parms.type + "_"+num2str(parms.test_year)+"_0";

% define model
[meanfunc, covfunc, likfunc, inffunc, prior] = model();
mu_ml = prior.slope(1);
sigma_mc = prior.intercept(2);

sigma_mc = 0.1;
i = 808; % 2016 Alaska Metcalfe
% i = 842; % 2016 Missouri Blunt

for i=820:842

year = raceinfos{i}{1};
state = raceinfos{i}{2}{1};
candidateName = raceinfos{i}{3};
trueVote = raceinfos{i}{4};
pvi = raceinfos{i}{5};
experienced = raceinfos{i}{6};
party = raceinfos{i}{7};

mu_b = computePrior(pvi, experienced, party, parms);
hyp.mean = [mu_b];
hyp.cov(4) = log(sigma_mc);

nz = 200;
XMIN = floor(xs{i}(1,1)/parms.BIN)*parms.BIN;
% XMIN = -60;
xstar = [linspace(XMIN,0,nz).',zeros(1,nz)',ones(1,nz)'];
[~, ~, fmu1, fs21] = gp(hyp, inffunc, meanfunc, covfunc, likfunc, xs{i}, ys{i}, xstar);
parms.prior1 = [mu_b, sigma_mc];

sigma_mc = 0.05;
hyp.cov(4) = log(sigma_mc);
[~, ~, fmu2, fs22] = gp(hyp, inffunc, meanfunc, covfunc, likfunc, xs{i}, ys{i}, xstar);
parms.prior2 = [mu_b, sigma_mc];


sigma_mc = abs(mu_b-0.5)/2;
hyp.cov(4) = log(sigma_mc);
[~, ~, fmu3, fs23] = gp(hyp, inffunc, meanfunc, covfunc, likfunc, xs{i}, ys{i}, xstar);
parms.prior = [mu_b, sigma_mc];


fig = figure('visible', 'off');
hold on;
plot(xs{i}(:,1), ys{i}, 'k.','MarkerSize',8); 
plot(xstar(:,1), fmu1); 
plot(xstar(:,1), fmu2);
plot(xstar(:,1), fmu3);

% scale range of polls to [0,1]
xlim([XMIN ,0]);
ylim([0.15,0.5]);


BIN = parms.BIN;
Nx = (abs(min(xs{i}(:,1)))/BIN);
Nx = abs(XMIN)/BIN;
XTICK = -BIN*[Nx:-1:0];
XTICKLABELS = cell(numel(XTICK),1);
for i=1:numel(XTICK)
XTICKLABELS{i} = num2str(-XTICK(i));
end

YT = yticks;
YTLABELS = cell(numel(YT),1);
for i=1:numel(YT)
YTLABELS{i} = num2str(100*YT(i))+"%";
end

set(gca, 'box', 'off', ...
'tickdir', 'out', ...
'xtick', XTICK, ...
'xticklabels', XTICKLABELS, ...
'ytick', YT, ...
'yticklabels', YTLABELS);


legend('Polling data' , 'Posterior mean: Main text', 'Posterior mean: H1', 'Posterior mean: H2', 'Location', 'Best');
xlabel("Days to election"); ylabel("Latent voter preference");


plot_title = year + " " + state + " " + candidateName;
filename = fullfile(plot_path+ '/'+plot_title + num2str(parms.j)+" mean.pdf");
set(fig, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
set(fig, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
print(fig, filename, '-dpdf','-r300');
close;
end
% 
% 
% % add incumbency
% 
% CNNdata = readData("data/CNNdata1992-2016.csv");
% years = unique(CNNdata.cycle);
% states = unique(CNNdata.state);
% YEAR = [];
% STATE = [];
% WINNER = [];
% for i = 1:numel(years)
%    for j = 1:numel(states)
%       [~, ~, candidateNames, v, pvi, experienced, parties] = getRaceCandidateData(CNNdata, years(i), states(j));
%       if isempty(v), continue; end
%       YEAR = [YEAR, years(i)];
%       STATE = [STATE, states(j)];
%       [~, k] = max(v);
%       WINNER = [WINNER, candidateNames(k)];
%    end
% end
% 
% WINNER_TABLE = table(YEAR',STATE',WINNER','VariableNames',{'year','state','incumbency'});
% writetable(WINNER_TABLE, "RR/winners.csv");
% 
% winner_all = readData("RR/winners_all.csv");
% 
% for i=1:numel(raceinfos)
%    year = raceinfos{i}{1};
%    state = raceinfos{i}{2};
%    c = cell2mat(raceinfos{i}{3});
%    incumbent = cell2mat(winner_all.incumbency(winner_all.year==(year-6) & strcmp(winner_all.state,state)));
%    if strcmp(c,incumbent)
%       raceinfos{i}{8} = 1; 
%    else
%       raceinfos{i}{8} = 0;
%    end
% end

% % presidential/senatoral vote share correlations
% president = readData("data/1976-2020-president.csv");
% senate = readData("data/1976-2020-senate.csv");
% 
% years = unique(president.year);
% states = unique(president.state);
% president_dem = zeros(length(states),length(years));
% president_rep = zeros(length(states),length(years));
% 
% 
% for i=1:length(years)
%     year = years(i);
%     for j=1:length(states)
%         state = cell2mat(states(j));
%         president_dem(j,i) = voteshare(president, 'DEMOCRAT', year, state);
%         president_rep(j,i) = voteshare(president, 'REPUBLICAN', year, state);
%     end 
% end
% 
% president_rep_cor = zeros(length(states),length(states));
% for i=1:length(states)
%     for j=1:length(states)
%         president_rep_cor(i,j) = corr(president_rep(i,:)', president_rep(j,:)');
%     end 
% end
% 
% president_dem_cor = zeros(length(states),length(states));
% for i=1:length(states)
%     for j=1:length(states)
%         president_dem_cor(i,j) = corr(president_dem(i,:)', president_dem(j,:)');
%     end 
% end
% writematrix(president_dem_cor, "RR/results/president_dem_cor.csv");
% 
% p = symrcm(president_dem_cor);
% 
% % imagesc(president_rep_cor);        % draw image and scale colormap to values range
% % colorbar;          % show color scale
% writematrix(president_rep_cor, "RR/results/president_rep_cor.csv");
% 
% 
% years = unique(senate.year);
% states = unique(senate.state);
% senate_dem = zeros(length(states),length(years));
% senate_rep = zeros(length(states),length(years));
% 
% for i=1:length(years)
%     year = years(i);
%     for j=1:length(states)
%         state = cell2mat(states(j));
%         senate_dem(j,i) = voteshare(senate, 'DEMOCRAT', year, state);
%         senate_rep(j,i) = voteshare(senate, 'REPUBLICAN', year, state);
%     end 
% end
% 
% senate_dem(senate_dem==0) = mean(senate_dem(senate_dem~=0));
% senate_rep(senate_rep==0) = mean(senate_rep(senate_rep~=0));
% 
% senate_rep_cor = zeros(length(states),length(states));
% for i=1:length(states)
%     for j=1:length(states)
%         senate_rep_cor(i,j) = corr(nonzeros(senate_rep(i,:)'), nonzeros(senate_rep(j,:)'));
%     end 
% end
% writematrix(senate_rep_cor, "RR/results/senate_rep_cor.csv");
% 
% senate_dem_cor = zeros(length(states),length(states));
% for i=1:length(states)
%     for j=1:length(states)
%         senate_dem_cor(i,j) = corr(nonzeros(senate_dem(i,:)'), nonzeros(senate_dem(j,:)'));
%     end 
% end
% writematrix(senate_dem_cor, "RR/results/senate_dem_cor.csv");
% 
% 
% function v = voteshare(data, party, year, state)
%     tmp = data(data.year==year & strcmp(data.state,state) & strcmp(data.party_simplified,party),:);
%     if size(tmp,1) == 0
%        v = 0;
%        return
%     end
%     v = sum(tmp.candidatevotes ./ tmp.totalvotes);
% end