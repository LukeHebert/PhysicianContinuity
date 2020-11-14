'''
Author: Luke Hebert
Date begun: Oct 20th, 2020
Description:
  -takes N input .CSV files (1 per year) and reads patient continuity data
  -statistically analyzes data by year using ANOVA (& checks assumptions)
    -evaluates each year's data for normality using QQ plot & Shapiro Wilks test
    -compares the variance between distributions
  -plots that data in multiple formats
'''

'''NOTES+++++++++++++++++++++++++++++++++++++++++++++++
*change the distribution plots so that all Y axes are the same (more fair &
    actually more impressive)
-add y label to jitter
-go re-learn & explain everything
-update GitHub
*consider another statistical test if there does seem to be a gradient difference
*consider making additional QQ plots to test against constant distribution (flat horizontal line)
'''


import sys, os, statistics
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

#the folder where the multiple CSVs are located
#one CSV for each year of collected data, named with the following format:
#YYYYMMDD_YYYYMMDD
in_folder = sys.argv[1]

#get the input file paths and put them in a list
in_paths = []
for root, dirs, files in os.walk(in_folder):
    in_paths.extend(files)
    break

#input files should be named in the format YYYYMMDD_YYYYMMDD
#so that they can be sorted by year
in_paths.sort()

#the slash string is used to help manipulate file or directory pathways
slash = '\\' if os.name == 'nt' else '/'

#dictionary with one key per input file (i.e. year)
#the values are lists of sets of this format:
#(age, resident MMCI, attending MMCI)
ara_dict = {}
year_labels = []
for path in in_paths:
    year = path.split(slash)[-1].replace('.csv','').replace('.CSV','')
    year_labels.append(year)
    ara_dict[year] = []
    with open(in_folder + slash + path, 'r') as in_file:
        for i,line in enumerate(in_file):
            line_list = line.replace('\r','').replace('\n','').split(',')
            if i == 0:
                age_index = line_list.index('Age')
                resident_index = line_list.index('MMCI resident')
                attending_index = line_list.index('MMCI attending')
            else:
                age = line_list[age_index]
                resident = line_list[resident_index]
                attending = line_list[attending_index]
                ara_dict[year].append((int(age),float(resident),float(attending))) #"subset"

#a list of colors to help with plot making. Must be changed if using more than 6 input years
colors = ['#723529','#D85A23','#F49E26','#F7D314','#A5CD0C','#4FB107']

#these "super"-lists are for the ANOVA later
res_lists, att_lists, age_lists = [], [], []

#check for normal distributions using QQ plots & Shapiro Wilks test
#QQ plot helpful video for understanding: https://youtu.be/okjYjClSjOg
for year,setlist in ara_dict.items():
    residents = [subset[1] for subset in setlist]
    Xs,Ys = stats.probplot(residents, dist='norm', fit=False, plot = plt)
    max_x = max(Xs)#used for text placement
    max_y = max(Ys)#used for text placement
    min_x = min(Xs)#used for text placement
    min_y = min(Ys)#used for text placement
    plt.title('QQ plot of Resident MMCI values\n (against normal distribution)')
    #calculate & place the Shapiro Wilk stats
    w_res, p_res = stats.shapiro(residents)
    shapiro_loc = [(0.01*(max_x-min_x))+min_x, (0.99*(max_y-min_y))+min_y]#location of shapiro stats text
    plt.text(shapiro_loc[0],shapiro_loc[1],'Shapiro W ' + str(round(w_res,2)) + ' & P ' + str('{:0.2e}'.format(p_res)), color='gray')
    #calculate and place the variance
    variance = round(statistics.pvariance(residents),4)
    var_loc = [(0.01*(max_x-min_x))+min_x, (0.9*(max_y-min_y))+min_y]#location of variance text
    plt.text(var_loc[0], var_loc[1], 'Variance: ' + str(variance), color='gray')
    plt.savefig(in_folder + slash + year + '_resQQ.png', dpi=800)
    plt.close()

    attendings = [subset[2] for subset in setlist]
    Xs,Ys = stats.probplot(attendings, dist='norm', fit=False, plot = plt)
    max_x = max(Xs)#used for text placement
    max_y = max(Ys)#used for text placement
    min_x = min(Xs)#used for text placement
    min_y = min(Ys)#used for text placement
    plt.title('QQ plot of Attending MMCI values\n (against normal distribution)')
    #calculate and place the Shapiro Wilk statistics
    w_att, p_att = stats.shapiro(attendings)
    shapiro_loc = [(0.01*(max_x-min_x))+min_x, (0.99*(max_y-min_y))+min_y]#location of shapiro stats text
    plt.text(shapiro_loc[0],shapiro_loc[1],'Shapiro W ' + str(round(w_att,2)) + ' & P ' + str('{:0.2e}'.format(p_att)), color='gray')
    #calculate and place the variance
    variance = round(statistics.variance(attendings),4)
    var_loc = [(0.01*(max_x-min_x))+min_x, (0.9*(max_y-min_y))+min_y]#location of variance text
    plt.text(var_loc[0], var_loc[1], 'Variance: ' + str(variance), color='gray')
    plt.savefig(in_folder + slash + year + '_attQQ.png', dpi=800)
    plt.close()

    #this part is for the ANOVA later
    ages = [subset[0] for subset in setlist]
    res_lists.append(residents)
    att_lists.append(attendings)
    age_lists.append(ages)

#apply ANOVA
#asterisk lets us pass dynamic number of arguments
anoF_res, anoP_res = stats.f_oneway(*res_lists)
anoF_att, anoP_att = stats.f_oneway(*att_lists)

with open(in_folder + slash + 'ANOVA_stats.txt', 'w') as out_file:
    out_file.write('One-way ANOVA F value for resident MMCI data: ' + str(anoF_res))
    out_file.write('\nOne-way ANOVA P value for resident MMCI data: ' + str(anoP_res))
    out_file.write('\nOne-way ANOVA F value for attending MMCI data: ' + str(anoF_att))
    out_file.write('\nOne-way ANOVA P value for attending MMCI data: ' + str(anoP_att))

#plot comparative distributions
##plot resident data
##this max subsection is just for establishing MMCI axis range
    #so that we can keep the same axis range for all subplots
axis_max = 0
for sublist in res_lists:
    hist_ys,hist_xs,temp = plt.hist(sublist)
    max_y = max(hist_ys)
    if max_y>axis_max:
        axis_max = int(max_y)
axis_max += 5
plt.close()
years_num = len(in_paths)
fig, axs = plt.subplots(years_num)
fig.suptitle('Yearly Resident MMCI Distributions')
plt.ylim(0,axis_max)
for i,sublist in enumerate(res_lists):
    axs[i].hist(sublist, bins='fd', color=colors[i].upper())#create the historgram
    axs[i].axvline(sum(sublist)/len(sublist), alpha=0.5, linestyle='--')#vertical line for arithmetic mean
    #this part also relies on input files being formatted YYYYMMDD_YYYYMMDD
    dates = year_labels[i].split('_')
    axs[i].set_ylabel(dates[0][:4] + '-' + dates[1][:4], fontsize=8)
#making a second y axis label right in the middle of the figure requires us to know which is the middle subplot
middle_index = int(round(len(att_lists)/2, 0)-1)
axs_right = axs[middle_index].twinx()#make another y axis that shares an x axis for this figure
axs_right.get_yaxis().set_ticks([1])#workaround for preventing right-axis label/figure overlap
axs_right.tick_params(colors='white')#workaround for preventing right-axis label/figure overlap
axs_right.set_ylabel('One Way ANOVA P-value: ' + str('{:0.2e}'.format(anoP_res)) + '\n ', fontsize=8, rotation=270)
plt.xlabel('MMCI values')
plt.savefig(in_folder + slash + 'resident_dists.png', dpi=800)
plt.close()

##plot attending data
##this max subsection is just for establishing MMCI axis range
    #so that we can keep the same axis range for all subplots
axis_max = 0
for sublist in att_lists:
    hist_ys,hist_xs,temp = plt.hist(sublist)
    max_y = max(hist_ys)
    if max_y>axis_max:
        axis_max = int(max_y)
axis_max += 5
plt.close()
fig, axs = plt.subplots(years_num)
fig.suptitle('Yearly Attending MMCI Distributions')
plt.ylim(0,axis_max)
for i,sublist in enumerate(att_lists):
    axs[i].hist(sublist, bins='fd', color=colors[i].upper())
    axs[i].axvline(sum(sublist)/len(sublist), alpha=0.5, linestyle='--')
    #this part also relies on input files being formatted YYYYMMDD_YYYYMMDD
    dates = year_labels[i].split('_')
    axs[i].set_ylabel(dates[0][:4] + '-' + dates[1][:4], fontsize=8)
#making a second y axis label right in the middle of the figure requires us to know which is the middle subplot
middle_index = int(round(len(att_lists)/2, 0)-1)
axs_right = axs[middle_index].twinx()#make another y axis that shares an x axis for this figure
axs_right.get_yaxis().set_ticks([1])#workaround for preventing right-axis label/figure overlap
axs_right.tick_params(colors='white')#workaround for preventing right-axis label/figure overlap
axs_right.set_ylabel('One Way ANOVA P-value: ' + str('{:0.2e}'.format(anoP_att)) + '\n ', fontsize=8, rotation=270)
plt.xlabel('MMCI values')
plt.savefig(in_folder + slash + 'attending_dists.png', dpi=800)
plt.close()

#make swarmplots
##plot resident data
fig, axs = plt.subplots(ncols=years_num, sharey=True)
fig.suptitle('Yearly Resident MMCI Values')
for i,sublist in enumerate(res_lists):
    color_list, younger, older = [],[],[]
    for j,age in enumerate(age_lists[i]):
        if age<24:
            color_list.append('blue')
            younger.append(res_lists[i][j])
        else:
            color_list.append('orange')
            older.append(res_lists[i][j])
    axs[i].scatter(x=[np.random.random() for each in sublist], y=sublist, c=color_list, s=1)
    axs[i].axhline((sum(sublist)/len(sublist)), alpha=0.5, linestyle='--')
    axs[i].axhline((sum(older)/len(older)), alpha=0.5, linestyle='-', color='orange')
    axs[i].axhline((sum(younger)/len(younger)), alpha=0.5, linestyle='-', color='blue')
    axs[i].tick_params(axis='x',colors='white')
    #this part also relies on input files being formatted YYYYMMDD_YYYYMMDD
    dates = year_labels[i].split('_')
    axs[i].set_xlabel(dates[0][:4] + '-' + dates[1][:4], fontsize=8)
axs[0].set_ylabel('MMCI')
plt.savefig(in_folder + slash + 'resident_jitter.png', dpi=800)
plt.close()

##plot attending data
fig, axs = plt.subplots(ncols=years_num, sharey=True)
fig.suptitle('Yearly Attending MMCI Values')
for i,sublist in enumerate(att_lists):
    color_list, younger, older = [],[],[]
    for j,age in enumerate(age_lists[i]):
        if age<24:
            color_list.append('blue')
            younger.append(att_lists[i][j])
        else:
            color_list.append('orange')
            older.append(att_lists[i][j])
    axs[i].scatter(x=[np.random.random() for each in sublist], y=sublist, c=color_list, s=1)
    axs[i].axhline((sum(sublist)/len(sublist)), alpha=0.5, linestyle='--')
    axs[i].axhline((sum(older)/len(older)), alpha=0.5, linestyle='-', color='orange')
    axs[i].axhline((sum(younger)/len(younger)), alpha=0.5, linestyle='-', color='blue')
    axs[i].tick_params(axis='x',colors='white')
    #this part also relies on input files being formatted YYYYMMDD_YYYYMMDD
    dates = year_labels[i].split('_')
    axs[i].set_xlabel(dates[0][:4] + '-' + dates[1][:4], fontsize=8)
axs[0].set_ylabel('MMCI')
plt.savefig(in_folder + slash + 'attending_jitter.png', dpi=800)
plt.close()
