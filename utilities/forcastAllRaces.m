function allRaces = forcastAllRaces(besthyp, xs, ys, raceinfos, plot_path, parms)
    % iterate cycle/state race
    [meanfunc, covfunc, likfunc, inffunc, prior] = model(parms);
    im = {@infPrior, inffunc, prior};
    allRaces = struct;
    nz = 200;
    for i = 1:size(xs,1)
        if isempty(xs{i}), continue; end
        republican = xs{i}(1,5);
        xstar = [linspace(xs{i}(1,1)-10,0,nz).',zeros(1,nz)',ones(1,nz)',...
            parms.nfirm*ones(1,nz)',republican*ones(1,nz)'];
        im{3}.mean{2}{2} = parms.a(i);
        hyp = full2one(besthyp, i, parms.ncandidates, parms.nfirm);
        [~, ~, fmu, fs2] = gp(hyp, im, meanfunc, covfunc, likfunc, xs{i}, ys{i}, xstar);
        
        fig = plot_posterior(fmu, fs2, xs{i}(:,1), ys{i}, xstar(:,1), i);
        predPoll = fmu(end);
        year = raceinfos{i}{1};
        state = raceinfos{i}{2}{1};
        candidateName = raceinfos{i}{3};
        trueVote = raceinfos{i}{4};
        plot_title = year + " " + state + " " + candidateName;
        title(plot_title);
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
        filename = fullfile(stateFolder, plot_title + ".jpg");
        saveas(fig, filename);
        close;
        disp(plot_title + " predicted winning rate: " + predPoll);
        disp(plot_title + " actual votes won: " + trueVote + newline);

        fn = char(state+""+year);
        fn = fn(~isspace(fn));
        if ~isfield(allRaces, fn)
            allRaces.(fn) = [predPoll, trueVote];
        else
            allRaces.(fn) = [allRaces.(fn), predPoll, trueVote];
        end
    end
end