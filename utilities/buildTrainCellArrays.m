function [xs, ys, raceinfos, g3] = buildTrainCellArrays(data, years, states)
%
% Obtain polling/fundemental features and election race metadata of given years/states
% Results are in the form of cell arrays. 
% input: 
%   - data: dataset that contains experience,
%           cycle,state,Candidateidentifier,Percentage_of_Vote_won_x,
%           samplesize,daysLeft,numberSupport,pollster,Republican,Democrat
%   - years: election years of to-be-collected features and metadata 
%   - states: election states of to-be-collected features and metadata
% 
% output:
%   - xs: cell array of polling/fundemental features.
%          Each element contains days before election, polling proportion and samplesize
%   - ys: cell array of pollings. Each element contains polling proportion.
%   - raceinfos: cell array of election race metadata. Each element
%   contains year, state, candidate name, actual vote share, pvi,
%   experience and party (-1 if republican, 1 if democratic, 0 if third party)

    % predefine cell array
    % assuming number of candidates does not exceed 1,000
    xs = cell(1000,1);
    ys = cell(1000,1);
    raceinfos = cell(1000,1);
    counter = 1;
    g3 = 0;
    for i = 1:numel(years)
       for j = 1:numel(states)
          [x, y, candidateNames, v, pvi, experienced, parties] = getRaceCandidateData(data, years(i), states(j));
          if isempty(x), continue; end
%            if numel(candidateNames)>2, g3 = g3+1; disp(years(i));disp(states(i)); end
          g3 = g3 + 1;
          for k = 1:numel(x)
             xs{counter} = x{k};
             ys{counter} = y{k};
             raceinfos{counter} = {years(i), states(j), candidateNames(k), v(k), pvi(k), experienced(k), parties(k)};
             counter = counter + 1;
          end
       end
    end

    counter = counter - 1;
    xs = xs(1:counter);
    ys = ys(1:counter);
    raceinfos = raceinfos(1:counter);
end