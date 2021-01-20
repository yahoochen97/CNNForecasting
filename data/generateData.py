import pandas as pd
import numpy as np
import us
from datetime import datetime

# 
# Generate up-to-date polling data for 2018 and 2020 races
# Required files:
#   - 2018results.csv: actual vote shares for 2018 candidates [candidate, vote share]
#   - 92-16nopoll.csv: null polls of no poll races from 1992 to 2016
#   - 2018nopolls.csv: null polls of no poll races in 2018
#   - 2020nopolls.csv: null polls of no poll races in 2020
#   - experiences.csv: binary experience of 2018/2020 candidates 
#   - primary.csv: primary results of 2020 races
#   - pvi2018.csv: cook partisan voting index of 2018 cycle
#   - pvi2020.csv: cook partisan voting index of 2020 cycle
#
# Output:
#   - CNNdata2018.csv
#   - CNNdata2020.csv
#

# read updated polling data
url = "https://projects.fivethirtyeight.com/polls-page/senate_polls.csv"
data = pd.read_csv(url, index_col=None)
# data = pd.read_csv("senate_polls.csv", index_col=None)

# drop useless columns and rename
data = data.drop(columns=["display_name"])
data = data.rename(columns={"answer":"name","candidate_party": "party", "sample_size":"samplesize", 
                            "fte_grade": "grade"})

# drop null sample size rows                           
data.drop(data[data.samplesize.isnull()].index, inplace=True)

# convert samplesize to integer
data.samplesize = data.samplesize.astype(int)

# mutate null election day
data.loc[data.election_date.isnull(), 'election_date'] = data.election_date[4]

# define state name to abbr mapping
toabbr = us.states.mapping('name', 'abbr')

# create emtpy feature lists
Candidateidentifier = ["" for i in range(data.shape[0])]
daysLeft = [0 for i in range(data.shape[0])]
nsup = [0 for i in range(data.shape[0])]
Republican = [0 for i in range(data.shape[0])]
Democrat = [0 for i in range(data.shape[0])]

# iterate observations and compute features
#   - Candidateidentifier: cycle+state_abbr+candidate_last_name
#   - daysLeft:  poll end date - election date (negative)
#   - nsup: number of supporters
#   - Republican: binary indicator of whether candidate is republican
#   - Democrat: binary indicator of whether candidate is demoratic
for i in range(data.shape[0]):
    Candidateidentifier[i] = str(data.cycle.iloc[i]) + toabbr[data.state.iloc[i]] + data.name.iloc[i]
    daysLeft[i] = -pd.date_range(data.start_date.iloc[i], data.election_date.iloc[i]).shape[0]
    nsup[i] = data.samplesize.iloc[i]*data.pct.iloc[i]//100
    if data.party.iloc[i] == "REP":
        Republican[i] = 1
    if data.party.iloc[i] == "DEM":
        Democrat[i] = 1
    
# drop useless columns
data["Candidateidentifier"] =  Candidateidentifier
data["daysLeft"] = daysLeft
data["numberSupport"] = nsup
data["Democrat"] = Democrat
data["Republican"] = Republican
data = data.drop(columns=["start_date", "end_date", "election_date", "pct", "name", "party"])

# filter polls over one year
data = data[data.daysLeft>=-365]
data = data.reset_index(drop=True)

# select 2018 data
data2018 = data[data.cycle==2018]

# the warning came from here
data2018.loc[(data2018.Candidateidentifier=="2018MSHyde-Smith"), "Candidateidentifier"] = "2018MSSmith"

# Remove duplicate polls with different types
# One polling id may have multiple question ids, corresponding to different types.
# Types of polls includes
#   - lv: likely voters
#   - rv: register voters
#   - v: voters
#   - a: adults
#
# Rules:
#   1. if the question has one candidate, remove it since the other candidate may not win the primary
#   2. if one polling id has multiple questions including likely voters, only keep likely voters.
#   3. if one polling id has multiple questions not including likely voters, remove adult voters.
#


# obtain 2018 actual vote shares
results =  pd.read_csv("./2018results.csv", index_col=None)


# build vote share dict
votes = {}
for i in range(results.shape[0]):
    if  results.iloc[i,1] > 5:
        votes[results.iloc[i,0]] = results.iloc[i,1]
        
# remove candidates not participating elections
c = list(votes.keys())
data2018 = data2018[data2018.Candidateidentifier.isin(c)]
data2018 = data2018.reset_index(drop=True)

# keep if flags == 0
flags = (data2018.cycle!=2018)
for pid in data2018.poll_id.unique():
    tmp = data2018[data2018.poll_id==pid]

    for qid in tmp.question_id.unique():
        if tmp[tmp.question_id==qid].candidate_name.unique().shape[0]==1:
            flags[(data2018.poll_id==pid) & (data2018.question_id==qid)] = 1

    if tmp.population.unique().shape[0]>=2:
        flags[(data2018.poll_id==pid) & (data2018.population!='lv')] = 1
        # if 'a' in tmp.population.unique():
        #     flags[(data2018.poll_id==pid) & (data2018.population=='a')] = 1
        # else:
        #     flags[(data2018.poll_id==pid) & (data2018.population!='lv')] = 1

# flags[data2018.Candidateidentifier.isin(['2018MSSmith', '2018MSEspy'])] = 0

data2018 = data2018[flags==0]

for pid in data2018.poll_id.unique():
    tmp = data2018[data2018.poll_id==pid]

    for qid in tmp.question_id.unique():
        data2018.loc[data2018.question_id==qid,"samplesize"] = sum(data2018[data2018.question_id==qid].numberSupport.values)

# Select non-partisan polls and selection columns.
# Partisan polls are sponsored by parties.

data2018 = data2018[data2018.partisan.isnull()]
data2018 = data2018[['cycle', 'state', 'pollster',
        'samplesize', 'candidate_name','Candidateidentifier', 'daysLeft',
       'numberSupport', 'Democrat', 'Republican']]
data2018 = data2018.drop_duplicates()
data2018 = data2018.reset_index(drop=True)


# add vote share colum to 2018 data
percentage = [0 for i in range(data2018.shape[0])]
for i in range(data2018.shape[0]):
    percentage[i] = votes[data2018.Candidateidentifier.iloc[i]]
    
data2018["Percentage.of.Vote.won.x"] = percentage

# create experience dict
data_experience = pd.read_csv("./experiences.csv", index_col=None)
c2e = {}
for c in data_experience.candidate.unique():
    c2e[c] = data_experience[data_experience.candidate==c].experienced.values[0]
    
data_pvi = pd.read_csv("pvi2018.csv", index_col=None)
s2pvi = {}
for s in data_pvi.state.unique():
    s2pvi[s] = data_pvi[data_pvi.state==s].pvi.values[0]

        
# add experience column to 2018 data
n = data2018.shape[0]
experienceds = np.zeros((n,))
pvis = np.zeros((n,))
for i in range(n):
    c = data2018.candidate_name[i]
    experienceds[i] = c2e[c]
    state = data2018.state[i]
    pvis[i] = s2pvi[state]
    
data2018['pvi'] = pvis
data2018['experienced'] = experienceds

# Code Angus King and Bernard Sanders as dem
data2018.loc[data2018.Candidateidentifier.isin(['2018VTSanders', '2018MEKing']),'Democrat'] = 1

# Indicate 2018 special election
data2018.loc[data2018.Candidateidentifier.isin(['2018MSSmith', '2018MSEspy']),'state'] = 'MississippiS'
data2018.loc[data2018.Candidateidentifier.isin(['2018MNSmith', '2018MNHousley']),'state'] = 'MinnesotaS'

# Concatenate no poll races
data2018nopoll = pd.read_csv("./2018nopolls.csv",na_filter = False)
data2018 = pd.concat([data2018,data2018nopoll])

# Write results to disk
data2018.to_csv("./CNNdata2018.csv",  index = False)


# drop null party sponsor data
data.drop(data[data.population.isnull()].index, inplace=True)

# select 2020 subset
data2020 = data[data.cycle==2020]

# print(data2020[data2020.state=='Arkansas'][['partisan',
#         'samplesize', 'candidate_name','daysLeft',
#        'numberSupport']])

# filter 2020 data with 2020 primary results
primaries = pd.read_csv("./primary.csv", index_col=None)
nominees = primaries.nominee.unique()
data2020 = data2020.loc[data2020.candidate_name.isin(nominees)]

# Remove duplicate polls with different types
# keep if flags == 0
flags = (data2020.cycle!=2020)

for pid in data2020.poll_id.unique():
    tmp = data2020[data2020.poll_id==pid]
    for qid in tmp.question_id.unique():
        if tmp[tmp.question_id==qid].candidate_name.unique().shape[0]==1:
            flags[(data2020.poll_id==pid) & (data2020.question_id==qid)] = 1

    if tmp.population.unique().shape[0]>=2:
        flags[(data2020.poll_id==pid) & (data2020.population!='lv')] = 1

data2020 = data2020[flags==0]

# Filter out party-sponsored polls
data2020 = data2020[~data2020.partisan.isin(['REP','DEM'])]

# Add Linn to Savage in Maine 2020

for pid in data2020.poll_id.loc[data2020.candidate_name=="Max Linn"]:
    data2020.loc[((data2020.candidate_name=="Lisa Savage") & (data2020.poll_id==pid)),"numberSupport"] = data2020.numberSupport.loc[((data2020.candidate_name=="Lisa Savage") & (data2020.poll_id==pid))].values[0] + data2020.numberSupport.loc[(data2020.candidate_name=="Max Linn") & (data2020.poll_id==pid)].values[0]

data2020 = data2020[data2020.candidate_name!="Max Linn"]

# drop third party
data2020 = data2020[(data2020.Democrat==1) | (data2020.Republican==1)]
# data2020 = data2020[(data2020.Democrat==1) | (data2020.Republican==1) | (data2020.candidate_name=="Lisa Savage")]

# Select feature columns
data2020 = data2020[['cycle', 'state', 'pollster',
        'samplesize', 'candidate_name','Candidateidentifier', 'daysLeft',
       'numberSupport', 'Democrat', 'Republican']]
data2020 = data2020.drop_duplicates()
data2020 = data2020.reset_index(drop=True)

# add experience column
data_experience = pd.read_csv("./experiences.csv", index_col=None)
c2e = {}
for i in range(data_experience.shape[0]):
    c2e[data_experience.candidate[i]] = data_experience.experienced[i]
    
# add pvi column
data_pvi = pd.read_csv("./pvi2020.csv", index_col=None)
tofull = us.states.mapping('abbr','name')
s2pvi = {}
for i in range(data_pvi.shape[0]):
    abbr = data_pvi.DISTRICT[i][0:2]
    state = tofull[abbr]
    if state in s2pvi:
        s2pvi[state].append(data_pvi.Raw_PVI[i])
    else:
        s2pvi[state] = [data_pvi.Raw_PVI[i]]
        
for state in s2pvi:
    s2pvi[state] = np.around(np.mean(s2pvi[state]),2)
    
n = data2020.shape[0]
experienceds = np.zeros((n,))
pvis = np.zeros((n,))
for i in range(n):
    c = data2020.candidate_name[i]
    experienceds[i] = c2e[c]
    state = data2020.state[i]
    pvis[i] = s2pvi[state]
    
data2020['pvi'] = pvis
data2020['experienced'] = experienceds

# Indicate Georgia Special Election
data2020.loc[data2020.candidate_name.isin(['Kelly Loeffler','Doug Collins','Matt Lieberman','Raphael Warnock','Ed Tarver']),'state']='GeorgiaS'

# Concatenate no poll races
# Write 2020 data to disk
data2020nopoll = pd.read_csv("./2020nopolls.csv",na_filter = False)
data2020 = pd.concat([data2020, data2020nopoll])
data2020.to_csv("./CNNdata2020.csv",  index = False)

