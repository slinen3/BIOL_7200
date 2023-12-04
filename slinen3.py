#!/usr/bin/env python3

"""
COVID-19 Hospitilization Rate Affect Script

This script reads COVID-19 hospitalization rate data from a CSV file downloaded from the CDC, performs data filtering and aggregation,
calculates average rates before and after the year 2020, conducts t-tests, and generates
cumulative and weekly prevalence plots based on gender, age, and race categories.


Parameters:
-----------
input_file : str
    The path to the CSV file containing COVID-19 data.

Returns:
--------
None

Usage:
------
./script.py <path_to_data_file.csv>

Note:
-----
Overall, COVID-19 did increase flu hospitilization rates cumulatively each year, however it did not affect average weekly rates yearly. 


"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
import numpy as np

input_file = sys.argv[1]

# specify col names
column_names = ["CATCHMENT", "NETWORK", "YEAR", "MMWR-YEAR", "MMWR-WEEK", "AGE CATEGORY", "SEX CATEGORY", "RACE CATEGORY", "CUMULATIVE RATE", "WEEKLY RATE"]

# read in data 
df = pd.read_csv(input_file, sep=',', skiprows=3, names=column_names)

# Select only the columns after "YEAR"
df = df.loc[:, 'MMWR-YEAR':]
df = df.dropna(subset=['CUMULATIVE RATE', 'WEEKLY RATE'])

#for pplot normalization later
scaling_factor = 0.1

# filter data for OVERALL changes (race, sex, age)
overall_df = df[(df['SEX CATEGORY'] == 'Overall') & (df['AGE CATEGORY'] == 'Overall') & (df['RACE CATEGORY'] == 'Overall')]
max_vals_overall = overall_df.groupby('MMWR-YEAR').max()
max_vals_overall['CUMULATIVE RATE']/=10

# filter data for female changes overall
female_df = df[(df['SEX CATEGORY'] == 'Female')]
max_vals_female = female_df.groupby('MMWR-YEAR').max()
max_vals_female['CUMULATIVE RATE']/=10

# filter data for male changes overall
male_df = df[(df['SEX CATEGORY'] == 'Male')]
max_vals_male = male_df.groupby('MMWR-YEAR').max()
max_vals_male['CUMULATIVE RATE']/=10

# filter data for AGE changes 
under18_df = df[(df['SEX CATEGORY'] == 'Overall') & (df['AGE CATEGORY'].isin(['0-4 yr', '< 18', '5-17 yr', '5-11 yr', '12-17 yr']))]
max_vals_under18 = under18_df.groupby('MMWR-YEAR').max()
max_vals_under18['CUMULATIVE RATE']/=10

over18_df = df[(df['SEX CATEGORY'] == 'Overall') & (df['AGE CATEGORY'].isin(['18-49 yr', '>= 18', '50-64 yr', '65+ yr', '65-74 yr', '75-84 yr', '85+']))]
max_vals_over18 = over18_df.groupby('MMWR-YEAR').max()
max_vals_under18['CUMULATIVE RATE']/=10

# filter data for RACE changes 
white_df = df[(df['SEX CATEGORY'] == 'Overall') & (df['AGE CATEGORY'] == 'Overall') & (df['RACE CATEGORY'] == 'White')]
max_vals_white = white_df.groupby('MMWR-YEAR').max()
max_vals_white['CUMULATIVE RATE']/=10

black_df = df[(df['SEX CATEGORY'] == 'Overall') & (df['AGE CATEGORY'] == 'Overall') & (df['RACE CATEGORY'] == 'Black')]
max_vals_black = black_df.groupby('MMWR-YEAR').max()
max_vals_black['CUMULATIVE RATE'] /= 10

hispanic_df = df[(df['SEX CATEGORY'] == 'Overall') & (df['AGE CATEGORY'] == 'Overall') & (df['RACE CATEGORY'] == 'Hispanic/Latino')]
max_vals_hispanic = hispanic_df.groupby('MMWR-YEAR').max()
max_vals_hispanic['CUMULATIVE RATE'] /= 10

asian_df = df[(df['SEX CATEGORY'] == 'Overall') & (df['AGE CATEGORY'] == 'Overall') & (df['RACE CATEGORY'] == 'Asian/Pacific Islander')]
max_vals_asian = asian_df.groupby('MMWR-YEAR').max()
max_vals_asian['CUMULATIVE RATE'] /= 10

amer_indian_df = df[(df['SEX CATEGORY'] == 'Overall') & (df['AGE CATEGORY'] == 'Overall') & (df['RACE CATEGORY'] == 'American Indian/Alaska Native')]
max_vals_amer_indian = amer_indian_df.groupby('MMWR-YEAR').max()
max_vals_amer_indian['CUMULATIVE RATE'] /= 10


# calc avg rates before and after 2020
rates_before_2020_cumulative = (df[df['MMWR-YEAR'] < 2020]['CUMULATIVE RATE'])/10
rates_after_2020_cumulative = (df[df['MMWR-YEAR'] >= 2020]['CUMULATIVE RATE'])/10

average_before_2020_cumulative = rates_before_2020_cumulative.mean()
average_after_2020_cumulative = rates_after_2020_cumulative.mean()

# calc p-value
p_value_cumulative = ttest_ind(rates_before_2020_cumulative, rates_after_2020_cumulative).pvalue
print(p_value_cumulative)

# calc std error mean
sem_before_2020_cumulative = rates_before_2020_cumulative.sem()
sem_after_2020_cumulative = rates_after_2020_cumulative.sem()

# plot cumulative data
fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(16, 8))

# race (cumulative)
ax1 = axes[0, 0].twinx()
ax1.plot(max_vals_white.index, max_vals_white['CUMULATIVE RATE'], linestyle='dotted', color='lightcoral', alpha=0.5, label='White')
ax1.plot(max_vals_black.index, max_vals_black['CUMULATIVE RATE'], linestyle='dotted', color='navajowhite', alpha=0.5, label='Black')
ax1.plot(max_vals_hispanic.index, max_vals_hispanic['CUMULATIVE RATE'], linestyle='dotted', color='thistle', alpha=0.5, label='Hispanic/Latino')
ax1.plot(max_vals_asian.index, max_vals_asian['CUMULATIVE RATE'], linestyle='dotted', color='cornflowerblue', alpha=0.5, label='Asian/Pacific Islander')
ax1.plot(max_vals_amer_indian.index, max_vals_amer_indian['CUMULATIVE RATE'], linestyle='dotted', color='lightgreen', alpha=0.5, label='American Indian/Alaska Native')
legend1 = ax1.legend(loc='upper left', fontsize='small')  # Adjust fontsize here
ax1.tick_params(axis='y', labelright=False)
ax1.yaxis.set_ticks([])
axes[0, 0].set_xlabel('Year')
axes[0, 0].set_ylabel('Max Cumulative Rate')
axes[0, 0].set_title('Cumulative Race Prevalence')


# overlay overall data plot 
axes[0, 0].plot(max_vals_overall.index, max_vals_overall['CUMULATIVE RATE'], label='Overall', color='k')
axes[0, 0].axvline(x=2020, color='red', linestyle='--', label='COVID-19')
axes[0, 0].errorbar([2019.5, 2020.5], [average_before_2020_cumulative, average_after_2020_cumulative], yerr=[sem_before_2020_cumulative, sem_after_2020_cumulative], fmt='o-', color='k', label=f'Overall (p={p_value_cumulative:.3f})')

axes[0, 0].plot([2019.5, 2020.5], [average_before_2020_cumulative, average_after_2020_cumulative], linestyle=':', color='black')


# gender (cumulative)
ax2 = axes[1, 0].twinx()
ax2.plot(max_vals_female.index, max_vals_female['CUMULATIVE RATE'], linestyle='dotted', color='lightcoral', alpha=0.5, label='Female')
ax2.plot(max_vals_male.index, max_vals_male['CUMULATIVE RATE'], linestyle='dotted', color='navajowhite', alpha=0.5, label='Male')
ax2.legend(loc='upper left')
ax2.tick_params(axis='y', labelright=False)
ax2.yaxis.set_ticks([])
axes[1, 0].set_xlabel('Year')
axes[1, 0].set_ylabel('Max Cumulative Rate')
axes[1, 0].set_title('Cumulative Gender Prevalence')

# overlay overall data plot 
axes[1, 0].plot(max_vals_overall.index, max_vals_overall['CUMULATIVE RATE'], label='Overall', color='k')
axes[1, 0].axvline(x=2020, color='red', linestyle='--', label='COVID-19')
axes[1, 0].errorbar([2019.5, 2020.5], [average_before_2020_cumulative, average_after_2020_cumulative], yerr=[sem_before_2020_cumulative, sem_after_2020_cumulative], fmt='o-', color='k', label=f'Overall (p={p_value_cumulative:.3f})')

axes[1, 0].plot([2019.5, 2020.5], [average_before_2020_cumulative, average_after_2020_cumulative], linestyle=':', color='black')

# age (cumulative)
ax3 = axes[2, 0].twinx()
ax3.plot(max_vals_under18.index, max_vals_under18['CUMULATIVE RATE'], linestyle='dotted', color='lightcoral', alpha=0.5, label='Under 18')
ax3.plot(max_vals_over18.index, max_vals_over18['CUMULATIVE RATE'], linestyle='dotted', color='navajowhite', alpha=0.5, label='18 and Over')
ax3.legend(loc='upper left')
ax3.tick_params(axis='y', labelright=False)
ax3.yaxis.set_ticks([])
axes[2, 0].set_xlabel('Year')
axes[2, 0].set_ylabel('Max Cumulative Rate')
axes[2, 0].set_title('Cumulative Age Prevalence')

# overlay overall data plot 
axes[2, 0].plot(max_vals_overall.index, max_vals_overall['CUMULATIVE RATE'], label='Overall', color='k')
axes[2, 0].axvline(x=2020, color='red', linestyle='--', label='COVID-19')
axes[2, 0].errorbar([2019.5, 2020.5], [average_before_2020_cumulative, average_after_2020_cumulative], yerr=[sem_before_2020_cumulative, sem_after_2020_cumulative], fmt='o-', color='k', label=f'Overall (p={p_value_cumulative:.3f})')

axes[2, 0].plot([2019.5, 2020.5], [average_before_2020_cumulative, average_after_2020_cumulative], linestyle=':', color='black')

# plot p-value of means 
for i in range(3):
	if p_value_cumulative < 0.05:
		axes[i, 0].text(2020.5, average_after_2020_cumulative, f'***p={p_value_cumulative:.10f}', verticalalignment='bottom', horizontalalignment='left', color='red', fontsize=10)
	else:
		axes[i, 0].text(2020.5, average_after_2020_cumulative, f'p={p_value_cumulative:.10f}', verticalalignment='bottom', horizontalalignment='left', color='red', fontsize=10)

## WEEKLY DATA PLOTTING
# calc avg rates before and after 2020 (weekly)
rates_before_2020_weekly = df[df['MMWR-YEAR'] < 2020]['WEEKLY RATE']
rates_after_2020_weekly = df[df['MMWR-YEAR'] >= 2020]['WEEKLY RATE']

average_before_2020_weekly = rates_before_2020_weekly.mean()
average_after_2020_weekly = rates_after_2020_weekly.mean()

# calc p-val (weekly)
p_value_weekly = ttest_ind(rates_before_2020_weekly, rates_after_2020_weekly).pvalue
print(p_value_weekly)
# Calculate standard error of the mean for Weekly rates
sem_before_2020_weekly = rates_before_2020_weekly.sem()
sem_after_2020_weekly = rates_after_2020_weekly.sem()

# race (weekly)
ax4 = axes[0, 1].twinx()
ax4.plot(max_vals_white.index, max_vals_white['WEEKLY RATE'], linestyle='dotted', color='lightcoral', alpha=0.5, label='White')
ax4.plot(max_vals_black.index, max_vals_black['WEEKLY RATE'], linestyle='dotted', color='navajowhite', alpha=0.5, label='Black')
ax4.plot(max_vals_hispanic.index, max_vals_hispanic['WEEKLY RATE'], linestyle='dotted', color='thistle', alpha=0.5, label='Hispanic/Latino')
ax4.plot(max_vals_asian.index, max_vals_asian['WEEKLY RATE'], linestyle='dotted', color='cornflowerblue', alpha=0.5, label='Asian/Pacific Islander')
ax4.plot(max_vals_amer_indian.index, max_vals_amer_indian['WEEKLY RATE'], linestyle='dotted', color='lightgreen', alpha=0.5, label='American Indian/Alaska Native')
legend4 = ax4.legend(loc='upper left', fontsize='x-small')  # Adjust fontsize here
ax4.tick_params(axis='y', labelright=False)
ax4.yaxis.set_ticks([])
axes[0, 1].set_xlabel('Year')
axes[0, 1].set_ylabel('Max Weekly Rate')
axes[0, 1].set_title('Weekly Race Prevalence')

# overlay overall data plot (weekly)
axes[0, 1].plot(max_vals_overall.index, max_vals_overall['WEEKLY RATE'], label='Overall', color='k')
axes[0, 1].axvline(x=2020, color='red', linestyle='--', label='COVID-19')
axes[0, 1].errorbar([2019.5, 2020.5], [average_before_2020_weekly, average_after_2020_weekly], yerr=[sem_before_2020_weekly, sem_after_2020_weekly], fmt='o-', color='k', label=f'Overall (p={p_value_weekly:.3f})')

axes[0, 1].plot([2019.5,2020.5], [average_before_2020_weekly, average_after_2020_weekly], linestyle=':', color='black')

# gender (weekly)
ax5 = axes[1, 1].twinx()
ax5.plot(max_vals_female.index, max_vals_female['WEEKLY RATE'], linestyle='dotted', color='lightcoral', alpha=0.5, label='Female')
ax5.plot(max_vals_male.index, max_vals_male['WEEKLY RATE'], linestyle='dotted', color='navajowhite', alpha=0.5, label='Male')
ax5.legend(loc='upper left')
ax5.tick_params(axis='y', labelright=False)
ax5.yaxis.set_ticks([])
axes[1, 1].set_xlabel('Year')
axes[1, 1].set_ylabel('Max Weekly Rate')
axes[1, 1].set_title('Weekly Gender Prevalence')

# overlay overall data plot (weekly)
axes[1, 1].plot(max_vals_overall.index, max_vals_overall['WEEKLY RATE'], label='Overall', color='k')
axes[1, 1].axvline(x=2020, color='red', linestyle='--', label='COVID-19')
axes[1, 1].errorbar([2019.5, 2020.5], [average_before_2020_weekly, average_after_2020_weekly], yerr=[sem_before_2020_weekly, sem_after_2020_weekly], fmt='o-', color='k', label=f'Overall (p={p_value_weekly:.3f})')

axes[1, 1].plot([2019.5, 2020.5], [average_before_2020_weekly, average_after_2020_weekly], linestyle=':', color='black')


# age (weekly)
ax6 = axes[2, 1].twinx()
ax6.plot(max_vals_under18.index, max_vals_under18['WEEKLY RATE'], linestyle='dotted', color='lightcoral', alpha=0.5, label='Under 18')
ax6.plot(max_vals_over18.index, max_vals_over18['WEEKLY RATE'], linestyle='dotted', color='navajowhite', alpha=0.5, label='18 and Over')
ax6.legend(loc='upper left')
ax6.tick_params(axis='y', labelright=False)
ax6.yaxis.set_ticks([])
axes[2, 1].set_xlabel('Year')
axes[2, 1].set_ylabel('Max Weekly Rate')
axes[2, 1].set_title('Weekly Age Prevalence')

# overlay overall data plot (weekly)
axes[2, 1].plot(max_vals_overall.index, max_vals_overall['WEEKLY RATE'], label='Overall', color='k')
axes[2, 1].axvline(x=2020, color='red', linestyle='--', label='COVID-19')
axes[2, 1].errorbar([2019.5, 2020.5], [average_before_2020_weekly, average_after_2020_weekly], yerr=[sem_before_2020_weekly, sem_after_2020_weekly], fmt='o-', color='k', label=f'Overall (p={p_value_weekly:.3f})')

axes[2, 1].plot([2019.5, 2020.5], [average_before_2020_weekly, average_after_2020_weekly], linestyle=':', color='black')

# plot p-value of means 
for i in range(3):
	if p_value_weekly < 0.05:
		axes[i, 1].text(2020.5, average_after_2020_weekly, f'***p={p_value_weekly:.10f}', verticalalignment='bottom', horizontalalignment='left', color='red', fontsize=10)
	else:
		axes[i, 1].text(2020.5, average_after_2020_weekly, f'p={p_value_weekly:.10f}', verticalalignment='bottom', horizontalalignment='left', color='red', fontsize=10)



plt.tight_layout()
plt.show()





 