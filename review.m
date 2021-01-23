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


% number of polls versus final margin
% fig = figure(1);
% a = CNNdata.samplesize(CNNdata.samplesize>1);
% b = CNNdata.numberSupport(CNNdata.samplesize>1) ./ CNNdata.samplesize(CNNdata.samplesize>1) - CNNdata.Percentage_of_Vote_won_x(CNNdata.samplesize>1)/100;
% scatter(a,abs(b),8,"filled");
% xlabel("# of polls"); ylabel("Final margin");
% filename = fullfile("/Users/yahoo/Desktop/pollvsmargin.pdf");
% set(fig, 'PaperPosition', [0 0 10 10]); %Position plot at left hand corner with width 5 and height 5.
% set(fig, 'PaperSize', [10 10]); %Set the paper to have width 5 and height 5.
% print(fig, filename, '-dpdf','-r300');


% split data into cell arrays
% obtain features (days before election) and regressor (polling proportions)
years = unique(CNNdata.cycle);
states = unique(CNNdata.state);
[xs, ys, raceinfos] = buildTrainCellArrays(CNNdata, years, states);

% different prior heuristics for intercept std
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

% precompute coefs of prior linear model of the linear trend intercept
parms.coefs = priorModel(CNNdata, parms.test_year);
plot_path = "RR/plots/" + parms.type + "_"+num2str(parms.test_year)+"_0";

% define model
[meanfunc, covfunc, likfunc, inffunc, prior] = model();
mu_ml = prior.slope(1);
sigma_mc = prior.intercept(2);

sigma_mc = 0.1;
i = 808; % 2016 Alaska Metcalfe
% i = 842; % 2016 Missouri Blunt

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
plot(xs{i}(:,1), ys{i}, 'k.'); 
plot(xstar(:,1), fmu1); 
plot(xstar(:,1), fmu2);
plot(xstar(:,1), fmu3);

% scale range of polls to [0,1]
% ylim([0,1]);
BIN = parms.BIN;
Nx = (abs(min(xs{i}(:,1)))/BIN);
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

xlim([min(XTICK) ,0]);

set(gca, 'box', 'off', ...
'tickdir', 'out', ...
'xtick', XTICK, ...
'xticklabels', XTICKLABELS, ...
'ytick', YT, ...
'yticklabels', YTLABELS);


legend('Polling data' , 'Posterior mean: our model', 'Posterior mean: h1', 'Posterior mean: h2', 'Location', 'Best');
xlabel("Days to election"); ylabel("Latent voter preference");


plot_title = year + " " + state + " " + candidateName;
filename = fullfile(plot_path+ '/'+plot_title + num2str(parms.j)+".pdf");
set(fig, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
set(fig, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
print(fig, filename, '-dpdf','-r300');
close;
