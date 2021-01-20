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
fig = figure(1);
a = CNNdata.samplesize(CNNdata.samplesize>1);
b = CNNdata.numberSupport(CNNdata.samplesize>1) ./ CNNdata.samplesize(CNNdata.samplesize>1) - CNNdata.Percentage_of_Vote_won_x(CNNdata.samplesize>1)/100;
scatter(a,b,8,"filled");
xlabel("# of polls"); ylabel("Final margin");
filename = fullfile("/Users/yahoo/Desktop/pollvsmargin.pdf");
set(fig, 'PaperPosition', [0 0 10 10]); %Position plot at left hand corner with width 5 and height 5.
set(fig, 'PaperSize', [10 10]); %Set the paper to have width 5 and height 5.
print(fig, filename, '-dpdf','-r300');
 