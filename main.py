import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import us
import csv


def plot_win(gp_path, stan_path):
    stan_data = pd.read_csv(stan_path)
    gp_data = pd.read_csv(gp_path)

    dem_win = []
    rep_win = []
    dem_metadata = []
    rep_metadata = []
    states = []

    for state in stan_data.state.unique():
        data = stan_data[stan_data.state==state]
        if len(data.candidate.unique())>2:
            continue
        states.append(state)
        for candidate in data.candidate.unique():
            win = data[(data.state==state) & (data.candidate==candidate)].win.values[0]
            party = gp_data[(gp_data.cycle==2020) & (gp_data.state==state) & (gp_data.candidate==candidate)].party.values[0]
            if party==1: 
                rep_win.append(win)
                rep_metadata.append((state,candidate))
            else:
                dem_win.append(win)
                dem_metadata.append((state,candidate))


    idx = np.argsort(rep_win)
    rep_win= np.array(rep_win)[idx]
    dem_win = np.array(dem_win)[idx]
    states = np.array(states)[idx]
    print(states)

    ind = np.arange(len(dem_win))  
    width = 0.5
    p1 = plt.barh(ind, dem_win, width, color='blue', alpha=0.7)
    p2 = plt.barh(ind, rep_win, width, color='red', left=dem_win, alpha=0.7)

    to_abbr = us.states.mapping('name', 'abbr')
    labels = []
    for state in states:
        if state == "Georgia Special":
            labels.append('GA S')
        else:
            labels.append(to_abbr[state])

    plt.xlabel('Winning Probability')
    plt.title('Winning Probability by party at day 42')
    plt.yticks(ind, labels)
    plt.xticks(np.linspace(0, 1, num=11))
    # plt.legend((p1[0], p2[0]), ('Dem', 'Rep'))
    ax = plt.axes()
    ax.xaxis.set_ticks_position('none') 
    ax.yaxis.set_ticks_position('none') 
    plt.savefig("plots/wp42.pdf")
    plt.show()
    plt.close()


def plot_vote(gp_path, stan_path):
    stan_data = pd.read_csv(stan_path)
    gp_data = pd.read_csv(gp_path)
    stan_data = stan_data.rename(columns={"median":"m"})

    dem_vote = []
    rep_vote = []
    dem_metadata = []
    rep_metadata = []
    states = []

    for state in stan_data.state.unique():
        data = stan_data[stan_data.state==state]
        if len(data.candidate.unique())>2:
            continue
        states.append(state)
        for candidate in data.candidate.unique():
            vote = data[(data.state==state) & (data.candidate==candidate)].m.values[0]
            party = gp_data[(gp_data.cycle==2020) & (gp_data.state==state) & (gp_data.candidate==candidate)].party.values[0]
            if party==1: 
                rep_vote.append(vote)
                rep_metadata.append((state,candidate))
            else:
                dem_vote.append(vote)
                dem_metadata.append((state,candidate))

    idx = np.argsort(rep_vote)
    dem_vote = np.array(dem_vote)[idx]
    rep_vote = np.array(rep_vote)[idx]
    states = np.array(states)[idx]

    ind = np.arange(len(dem_vote))  
    width = 0.5
    p1 = plt.barh(ind, dem_vote, width, color='blue', alpha=0.8)
    p2 = plt.barh(ind, rep_vote, width, color='red', left=1-np.array(rep_vote), alpha=0.8)

    to_abbr = us.states.mapping('name', 'abbr')
    labels = []
    for state in states:
        if state == "Georgia Special":
            labels.append('GA S')
        else:
            labels.append(to_abbr[state])

    plt.ylabel('Vote')
    plt.title('Forecasted Vote by party at day 42')
    plt.yticks(ind, labels)
    plt.xticks(np.linspace(0, 1, num=11))
    # plt.legend((p1[0], p2[0]), ('Dem', 'Rep'))
    plt.savefig("plots/vote42.png")
    plt.show()
    plt.close()


def plot_vote_percentage(gp_path, stan_path):
    stan_data = pd.read_csv(stan_path)
    gp_data = pd.read_csv(gp_path)
    stan_data = stan_data.rename(columns={"median":"m"})

    dem_win = []
    rep_win = []
    dem_vote = []
    rep_vote = []
    dem_u = []
    rep_u = []
    dem_l = []
    rep_l = []
    dem_metadata = []
    rep_metadata = []
    states = []

    for state in stan_data.state.unique():
        data = stan_data[stan_data.state==state]
        if len(data.candidate.unique())>2:
            continue
        states.append(state)
        for candidate in data.candidate.unique():
            win = data[(data.state==state) & (data.candidate==candidate)].win.values[0]
            vote = data[(data.state==state) & (data.candidate==candidate)].m.values[0]
            party = gp_data[(gp_data.cycle==2020) & (gp_data.state==state) & (gp_data.candidate==candidate)].party.values[0]
            u = data[(data.state==state) & (data.candidate==candidate)].upper95.values[0]
            l = data[(data.state==state) & (data.candidate==candidate)].lower95.values[0]
            if party==1: 
                rep_vote.append(vote)
                rep_metadata.append((state,candidate))
                rep_u.append(u)
                rep_l.append(l)
                rep_win.append(win)
            else:
                dem_vote.append(vote)
                dem_metadata.append((state,candidate))
                dem_u.append(u)
                dem_l.append(l)
                dem_win.append(win)

    idx = np.argsort(dem_win)
    dem_win = np.array(dem_win)[idx]
    rep_win = np.array(rep_win)[idx]
    dem_vote = np.array(dem_vote)[idx]
    rep_vote = np.array(rep_vote)[idx]
    dem_u = np.array(dem_u)[idx]
    rep_u= np.array(rep_u)[idx]
    dem_l = np.array(dem_l)[idx]
    rep_l = np.array(rep_l)[idx]
    states = np.array(states)[idx]

    ind = np.arange(len(dem_vote))  
    p1 = plt.errorbar(ind, dem_vote, dem_ci, marker='x', mfc='blue',
         mec='blue', ms=2, mew=2)
    p2 = plt.errorbar(ind, rep_vote, rep_ci, marker='o', mfc='red',
         mec='red', ms=2, mew=2)

    # p1 = plt.scatter(ind, dem_vote, color='blue', marker='x', alpha=0.8)
    # p2 = plt.scatter(ind, rep_vote, color='red', marker='o', alpha=0.8)

    to_abbr = us.states.mapping('name', 'abbr')
    labels = []
    for state in states:
        if state == "Georgia Special":
            labels.append('GA S')
        else:
            labels.append(to_abbr[state])

    plt.ylabel('Vote')
    plt.title('Forecasted Vote by party at day 42')
    plt.xticks(ind, labels)
    plt.yticks(np.linspace(0, 1, num=11))
    plt.legend((p1, p2), ('Dem', 'Rep'))
    plt.savefig("plots/vote42.png")
    plt.show()
    plt.close()


def utility(data, party):
    to_abbr = us.states.mapping('name', 'abbr')
    abbrs = []
    for state in data.state.unique():
        if state == "MinnesotaS":
            abbrs.append('MNS')
        elif state == "MississippiS":
            abbrs.append('MSS')
        else:
            abbrs.append(to_abbr[state])

    abbrs = np.sort(abbrs)

    abbrs  = abbrs[::-1][:len(abbrs )]
    
    s2x = {state:i for i,state in enumerate(abbrs)}

    data = data[data.party==party]
    y = data.m*100
    l = data.lower95*100
    u = data.upper95*100
    v = data.vote*100
    
    x = []
    for state in data.state:
        if state == "MinnesotaS":
            x.append(s2x['MNS'])
        elif state == "MississippiS":
            x.append(s2x['MSS'])
        else:
            x.append(s2x[to_abbr[state]])

    s = []
    for state in data.state:
        if state == "MinnesotaS":
            s.append('MNS')
        elif state == "MississippiS":
            s.append('MSS')
        else:
            s.append(to_abbr[state])


    idx = np.argsort(s)[::-1][:len(s)]
    s = np.array(s)[idx]
    x = np.array(x)[idx]
    y = np.array(y)[idx]
    l = np.array(l)[idx]
    u = np.array(u)[idx]
    v = np.array(v)[idx]

    lower_error = y - l
    upper_error = u - y
    asymmetric_error = [lower_error, upper_error]

    return x, y, asymmetric_error, v, s


def plot_2018():
    stan_path = "results/stan_LOOGP_2018day7_59.csv"
    data = pd.read_csv(stan_path)
    data = data.rename(columns={"median":"m"})
    x1, y1, asymmetric_error1, v1, s1 = utility(data, party=1)
    x2, y2, asymmetric_error2, v2, s2 = utility(data, party=-1)
   
    fig= plt.figure(figsize=(11,11))
    shift1 = np.array([0.1 for _ in range(len(x1))])
    shift1[-3] = -0.1
    p1 = plt.errorbar(y1, x1-shift1, xerr=asymmetric_error1, fmt='.', elinewidth=2,  label="DEM", color="blue", alpha=0.5)
    p2 = plt.errorbar(y2, x2+0.1, xerr=asymmetric_error2, fmt='.', elinewidth=2, label="REP", color="red", alpha=0.5)
    p3 = plt.scatter(v1, x1-shift1, marker='*',label="DEM VOTE", color="blue")
    p4 = plt.scatter(v2, x2+0.1, marker='*',label="REP VOTE", color="red")

    plt.axvline(x=50, color='grey', linestyle='--')
    plt.axvline(x=25, color='grey', linestyle='--')
    plt.axvline(x=75, color='grey', linestyle='--')

    s1 = np.unique(s1)
    s1 = np.sort(s1)
    s1 = s1[::-1][:len(s1)]

    for i in range(len(s1)):
        plt.axhline(y=i, color='grey', linestyle='-', alpha=0.1)

    plt.yticks(range(len(s1)), s1)
    plt.xticks([0,25,50,75,100])


    STATE_COLORS = []
    for s in s1:
        if s in ['AZ', 'NV']:
            STATE_COLORS.append('red')
        else:
            STATE_COLORS.append('black')

    for ticklabel, tickcolor in zip(plt.gca().get_yticklabels(), STATE_COLORS):
        ticklabel.set_color(tickcolor)

    plt.xlabel('Posterior Vote (%)')
    plt.title('Posterior credible intervals of vote share for major party candidates')
    plt.legend((p1, p2, p3, p4), ('DEM 95% CI', 'REP 95% CI',"DEM VOTE","REP VOTE"), loc='best')
    plt.savefig('plots/2018day7.pdf', dpi=300)
    plt.show()

def plot_2016():
    stan_path = "results/stan_LOOGP_2016day7_59.csv"
    data = pd.read_csv(stan_path)
    data = data.rename(columns={"median":"m"})
    x1, y1, asymmetric_error1, v1, s1 = utility(data, party=1)
    x2, y2, asymmetric_error2, v2, s2 = utility(data, party=-1)
   
    fig= plt.figure(figsize=(11,11))
    shift1 = np.array([0.1 for _ in range(len(x1))])
    # shift1[-3] = -0.1
    p1 = plt.errorbar(y1, x1-shift1, xerr=asymmetric_error1, fmt='.', elinewidth=2,  label="DEM", color="blue", alpha=0.5)
    p2 = plt.errorbar(y2, x2+0.1, xerr=asymmetric_error2, fmt='.', elinewidth=2, label="REP", color="red", alpha=0.5)
    p3 = plt.scatter(v1, x1-shift1, marker='*',label="DEM VOTE", color="blue")
    p4 = plt.scatter(v2, x2+0.1, marker='*',label="REP VOTE", color="red")

    plt.axvline(x=50, color='grey', linestyle='--')
    plt.axvline(x=25, color='grey', linestyle='--')
    plt.axvline(x=75, color='grey', linestyle='--')

    s1 = np.unique(s1)
    s1 = np.sort(s1)
    s1 = s1[::-1][:len(s1)]

    for i in range(len(s1)):
        plt.axhline(y=i, color='grey', linestyle='-', alpha=0.1)

    plt.yticks(range(len(s1)), s1)
    plt.xticks([0,25,50,75,100])


    STATE_COLORS = []
    for s in s1:
        if s in ['MO', 'WI','PA']:
            STATE_COLORS.append('red')
        else:
            STATE_COLORS.append('black')

    for ticklabel, tickcolor in zip(plt.gca().get_yticklabels(), STATE_COLORS):
        ticklabel.set_color(tickcolor)

    plt.xlabel('Posterior Vote (%)')
    plt.title('Posterior credible intervals of vote share for major party candidates')
    plt.legend((p1, p2, p3, p4), ('DEM 95% CI', 'REP 95% CI',"DEM VOTE","REP VOTE"), loc='best')
    plt.savefig('plots/2016day7.pdf', dpi=300)
    plt.show()


def plot_2020():
    stan_path = "data/2020day0_result.csv"
    data = pd.read_csv(stan_path)
    data = data[data.state!="Arkansas"]
    data = data.reset_index(drop=True)
    data.loc[data.candidate=="Harrington", "party"] = 1
    data = data.rename(columns={"median":"m"})
    x1, y1, asymmetric_error1, v1, s1 = utility(data, party=1)
    x2, y2, asymmetric_error2, v2, s2 = utility(data, party=-1)
   
    fig= plt.figure(figsize=(11,11))
    shift1 = np.array([0.1 for _ in range(len(x1))])
    p1 = plt.errorbar(y1, x1-shift1, xerr=asymmetric_error1, fmt='.', elinewidth=2,  label="DEM", color="blue", alpha=0.5)
    p2 = plt.errorbar(y2, x2+0.1, xerr=asymmetric_error2, fmt='.', elinewidth=2, label="REP", color="red", alpha=0.5)
    p3 = plt.scatter(v1, x1-shift1, marker='*',label="DEM VOTE", color="blue")
    p4 = plt.scatter(v2, x2+0.1, marker='*',label="REP VOTE", color="red")

    plt.axvline(x=50, color='grey', linestyle='--')
    plt.axvline(x=25, color='grey', linestyle='--')
    plt.axvline(x=75, color='grey', linestyle='--')

    s1 = np.unique(s1)
    s1 = np.sort(s1)
    s1 = s1[::-1][:len(s1)]

    for i in range(len(s1)):
        plt.axhline(y=i, color='grey', linestyle='-', alpha=0.1)

    plt.yticks(range(len(s1)), s1)
    plt.xticks([0,25,50,75,100])


    STATE_COLORS = []
    for s in s1:
        if s in ['ME', 'NC']:
            STATE_COLORS.append('red')
        else:
            STATE_COLORS.append('black')

    for ticklabel, tickcolor in zip(plt.gca().get_yticklabels(), STATE_COLORS):
        ticklabel.set_color(tickcolor)

    plt.xlabel('Posterior Vote (%)')
    plt.title('Posterior credible intervals of vote share for major party candidates')
    plt.legend((p1, p2, p3, p4), ('DEM 95% CI', 'REP 95% CI',"DEM VOTE","REP VOTE"), loc=7)
    plt.savefig('plots/2020day0_result.pdf', dpi=300)
    plt.show()


def economist_2020():
    stan_path = "data/Economist_day_0.csv"
    data = pd.read_csv(stan_path, index_col=False)

    data = data[~data.state.isin(["GA-S3","LA-S2"])]
    data = data.reset_index(drop=True)
    for i in range(len(data)):
        tmp = data.state.loc[i]
        data.loc[i,'state'] = tmp[0:2]


    data = data.rename(columns={"pmean":"m"})

    PARTY_DICT = {"democratic": 1, "republican":-1}
    data['party'] = data['party'].map(PARTY_DICT)

    result2020 = pd.read_csv("data/2020results.csv", index_col=False)
    result2020 = result2020[~result2020.Candidateidentifier.isin(
        ["2020GACollins", "2020GALoeffler","2020GATarver","2020GAWarnock"])]
    result2020 = result2020.reset_index(drop=True)

    PS_V = {}
    for i in range(len(result2020)):
        s = result2020.Candidateidentifier[i][4:6]
        p = result2020.Party[i]
        v = result2020.loc[i,"Percentage.of.Vote.won.x"]
        if s in PS_V:
            PS_V[s][p] = v
        else:
            PS_V[s] = {p: v}

    for s in PS_V:
        dem_v = PS_V[s][1]
        rep_v = PS_V[s][-1]
        tmp = dem_v + rep_v
        dem_v = dem_v / tmp
        rep_v = rep_v / tmp
        PS_V[s][1] = dem_v
        PS_V[s][-1] = rep_v


    vote = [0 for _ in range(len(data))]

    for i in range(len(data)):
        p = data.party[i]
        s = data.state[i]
        v = PS_V[s][p]
        vote[i] = v

    data['vote'] = vote


    data1 = data[data.party==1]
    y1 = data1.m*100
    l1 = data1.lower95*100
    u1 = data1.upper95*100
    v1 = data1.vote*100

    lower_error1 = y1 - l1
    upper_error1 = u1 - y1
    asymmetric_error1 = [lower_error1, upper_error1]

    data2 = data[data.party==-1]
    y2 = data2.m*100
    l2 = data2.lower95*100
    u2 = data2.upper95*100
    v2 = data2.vote*100

    # print(np.sum((v2>u2) | (v2<l2)) + np.sum((v1>u1) | (v1<l1)) )

    lower_error2 = y2 - l2
    upper_error2 = u2 - y2
    asymmetric_error2 = [lower_error2, upper_error2]

    x1 = np.array([i for i in range(len(y1))])
    x2 = np.array([i for i in range(len(y2))])
     
    s1 = data1.state
   
    fig= plt.figure(figsize=(11,11))
    shift1 = np.array([0.1 for _ in range(len(x1))])
    p1 = plt.errorbar(y1, x1-0.1, xerr=asymmetric_error1, fmt='.', elinewidth=2,  label="DEM", color="blue", alpha=0.5)
    p2 = plt.errorbar(y2, x2+0.1, xerr=asymmetric_error2, fmt='.', elinewidth=2, label="REP", color="red", alpha=0.5)
    p3 = plt.scatter(v1, x1-0.1, marker='*',label="DEM VOTE", color="blue")
    p4 = plt.scatter(v2, x2+0.1, marker='*',label="REP VOTE", color="red")

    plt.axvline(x=50, color='grey', linestyle='--')
    plt.axvline(x=25, color='grey', linestyle='--')
    plt.axvline(x=75, color='grey', linestyle='--')


    for i in range(len(s1)):
        plt.axhline(y=i, color='grey', linestyle='-', alpha=0.1)

    plt.yticks(range(len(s1)), s1)
    plt.xticks([0,25,50,75,100])


    STATE_COLORS = []
    for s in s1:
        if s in ['ME', 'IA','GA']:
            STATE_COLORS.append('red')
        else:
            STATE_COLORS.append('black')

    for ticklabel, tickcolor in zip(plt.gca().get_yticklabels(), STATE_COLORS):
        ticklabel.set_color(tickcolor)

    plt.xlabel('Posterior Vote (%)')
    plt.title('Posterior credible intervals of vote share for major party candidates')
    plt.legend((p1, p2, p3, p4), ('DEM 95% CI', 'REP 95% CI',"DEM VOTE","REP VOTE"), loc=7)
    plt.savefig('plots/economist_result.pdf', dpi=300)
    plt.show()


def fivethirtyeight(TYPES, stan_path, result_path):

    coverage = np.zeros((len(TYPES)+2,))
    acc = np.zeros((len(TYPES)+2,))
    rmse = np.zeros((len(TYPES)+2,))
    data = pd.read_csv(stan_path, index_col=False)
    data = data.rename(columns={"pmean":"m"})
    data = data[data.state!="Arkansas"]
    data = data.reset_index(drop=True)

    result2020 = pd.read_csv(result_path, index_col=False)
    result2020= result2020[~result2020.Candidateidentifier.isin(
        ["2020GACollins", "2020GALoeffler","2020GATarver","2020GAWarnock","2020MESavage"])]
    result2020 = result2020.reset_index(drop=True)

    result2020 = result2020.rename(columns={"Percentage.of.Vote.won.x":"vote","Party":"party"})

    result2020_state = []
    for i in range(len(result2020)):
        result2020_state.append(result2020.Candidateidentifier[i][4:6])

    result2020["state"] = result2020_state

    result2020 = result2020[~result2020.state.isin(["LA"])]
    result2020 = result2020.reset_index(drop=True)

    for i, TYPE in enumerate(TYPES):
        data538 = pd.read_csv("data/538_"+TYPE+"_day_0.csv", index_col=False)
        data538 = data538.rename(columns={"pmean":"m"})
        PARTY_DICT = {"D": 1, "R":-1, "I":0, "O":2}
        data538['party'] = data538['party'].map(PARTY_DICT)
        to_abbr = us.states.mapping('name', 'abbr')

        for state in data.state.unique():
            state_data = data[data.state==state]
            s = to_abbr[state]
            parties = state_data.party.values
            votes = result2020.loc[result2020.state==s, ["vote","party"]].sort_values(by=["party"]).vote.values
            state_538data = data538[(data538.state==s) & (data538.party.isin(parties))].sort_values(by=["party"])
            # tmp = np.sum(state_538data.m.values)
            tmp = 1
            upper90 = state_538data.upper90.values/tmp
            lower10 = state_538data.lower10.values/tmp
            m = state_538data.m.values

            norm_votes = votes/np.sum(votes)

            if np.argmax(votes) == np.argmax(m):
                acc[1+i] += 1

            data_m = state_data.sort_values(by=["party"]).m.values
            if np.argmax(votes) == np.argmax(data_m):
                acc[0] += 1

            for j,vote in enumerate(votes):
                rmse[1+i] += ((vote-m[j])/100)**2
                rmse[0] += (norm_votes[j]-state_data.sort_values(by=["party"]).m.values[j])**2
                if vote>upper90[j] or vote<lower10[j]:
                    coverage[1+i] += 1
                if norm_votes[j]>state_data.sort_values(by=["party"]).upper90.values[j] or norm_votes[j]<state_data.sort_values(by=["party"]).lower10.values[j]:
                    coverage[0] += 1


    dataEcon = pd.read_csv("data/Economist_day_0.csv", index_col=False)
    dataEcon = dataEcon.rename(columns={"pmean":"m"})
    PARTY_DICT = {"democratic": 1, "republican":-1}
    dataEcon['party'] = dataEcon['party'].map(PARTY_DICT)
    to_abbr = us.states.mapping('name', 'abbr')

    dataEcon['state'] = dataEcon['state'].astype(str).str[0:2]

    for state in data.state.unique():
        state_data = data[data.state==state]
        s = to_abbr[state]
        parties = state_data.party.values
        votes = result2020.loc[result2020.state==s, ["vote","party"]].sort_values(by=["party"]).vote.values
        state_Econdata = dataEcon[(dataEcon.state==s) & (dataEcon.party.isin(parties))].sort_values(by=["party"])
        # tmp = np.sum(state_538data.m.values)
        tmp = 1
        upper90 = state_Econdata.upper95.values/tmp
        lower10 = state_Econdata.lower95.values/tmp
        m = state_Econdata.m.values
        if len(m)==1:
            continue
            m = [1,0]
            upper90 = [upper90, 0]
            lower10 = [lower10, 0]

        votes = votes/np.sum(votes)
        if np.argmax(votes) == np.argmax(m):
            acc[-1] += 1

        for j,vote in enumerate(votes):
            rmse[-1] += ((vote-m[j]))**2
            if vote>upper90[j] or vote<lower10[j]:
                coverage[-1] += 1
                print(s)


    coverage[0] /= len(TYPES)
    rmse[0] /= len(TYPES)
    rmse /= len(data)
    rmse = np.sqrt(rmse)
    coverage = 1 - coverage / len(data)
    acc[0] = (acc[0] / 3) / 32
    acc[1:] = acc[1:] / len(data.state.unique())
    print(acc)
    print(coverage)
    print(rmse)
    
    with open('results/coverage80.csv','w') as file:
        file.write(',GP+DR,')
        for TYPE in TYPES:
            file.write("538 "+TYPE.upper()+',')
        file.write("Economist,")
        file.write('\n')
        file.write("80% CI,")
        for c in coverage:
            file.write(str(c)+",")
        file.write('\n')
        file.write("RMSE,")
        for c in rmse:
            file.write(str(c)+",")
        file.write('\n')
        file.write("Acc,")
        for c in acc:
            file.write(str(c)+",")
        file.write('\n')


def main():
    # gp_path = "results/LOOGP_2018day7_59.csv"
    # stan_path = "results/stan_LOOGP_2018day7_59.csv"
    # plot_win(gp_path, stan_path)
    # plot_2016()
    # plot_2020()
    # economist_2020()

    TYPES = ["classic", "deluxe", "lite"]

    fivethirtyeight(TYPES, "data/2020day0_80.csv", "data/2020results.csv")

if __name__ == "__main__":
    main()
