# Plot timeseries of regional mean OSTIA SST, multi years overlaid
# Arguments:
#  $1  region id (11=NWShelf, 16=Northern North Sea)

# exec(open('./marine_heatwave_ostiaTS_maremon.py').read())

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import xarray as xr
from datetime import datetime,timedelta
from dateutil.relativedelta import relativedelta
import sys
import seaborn as sns

import netCDF4
import cf_units as unit
import iris
import pandas as pd
import datetime as dt
import numpy as np
import matplotlib

def to_annual(dseries, stat='mean'):
    # dseries = dseries.to_frame()
    dseries['day_of_year'] = dseries.index.map(lambda x: x.strftime("%j"))
    if stat == 'mean':
        annual_stat = dseries.groupby('day_of_year').mean()
    elif stat == 'q10':
        annual_stat = dseries.groupby('day_of_year').quantile(0.1)
    elif stat == 'q90':
        annual_stat = dseries.groupby('day_of_year').quantile(0.9)
    else:
        print('stat {} not supported, only mean, q10, q90'.format('stat'))
    annual_stat = annual_stat.query('day_of_year != "060"')
    print(annual_stat.index)
    annual_stat.index = pd.date_range("2023-01-01", "2023-12-31")
    return annual_stat

# plt.rcParams['text.usetex'] = True
font = {#'family': 'serif',
        'weight' : 'normal',
        'size'   : 12}
matplotlib.rc('font', **font)

region = 'WAKELIN_SHELF'

DIR = "/data/users/frjr/Heatwave2023/RegionalMeans.OSTIA"
DIR_CLIM = ""

for region in ['NWS', 'UKV', 'North_Sea']:
    fsize = (10, 4)
    fig, ax = plt.subplots(figsize=fsize)

    #### Load climatology (Segolene) #####
    for region_wak in ('ostia_cdr_{}_sst'.format(region), 'ostia_nrt_{}_sst'.format(region)):  # range(1, 14):
        # region = region_wak + 13

        # region_dict = {'WAKELIN_3': 16, # Northern North Sea
        #                'WAKELIN_SHELF': 11,
        #                'WAKELIN_2': 15, # Central North Sea
        #                'WAKELIN_8': 21  # Irish shelf
        #                }

        filename = '/data/users/ofrd-mopa/maremon/datasets/database_ofrd-mopa.csv'
        # f = netCDF4.Dataset(filename)
        # time = f.variables['time'][:]
        # model_seconds = unit.num2date(time[:], 'seconds since 1981-01-01', calendar='standard').data
        # cubes = iris.load_cube('/home/h02/frjr/Heatwave2023/OSTIA/regmeans_OSTIA.nc', 'Regional mean sea surface temperature')
        # out_series = pd.Series(cubes[:, region, 0, 0].data, index=model_seconds.astype(dt.datetime))
        a = pd.read_csv(filename)
        a = a.set_index('date')
        a.index = pd.to_datetime(a.index, format='%Y-%m-%d')
        out_series = pd.DataFrame(a[region_wak])
        out_series.to_pickle(region_wak+".pkl")
        # library(reticulate)
        # out_series <- py_read_pickle("ostia_nrt_NWS_sst.pkl")
        # out_series <- py_load_object("ostia_nrt_NWS_sst.pkl")
        out_series.to_csv(region_wak+".csv")






        # mask = out_series.index.values < dt.datetime(2013, 1, 1, 0, 0, 0)
        cropped_series = out_series.query("date >= '1982-01-01' and date < '2013-01-01'")
        mean_1982_2012 = to_annual(cropped_series, 'mean')
        q10 = to_annual(cropped_series, 'q10')
        q90 = to_annual(cropped_series, 'q90')
        cropped_series = out_series.query("date < '2003-01-01'")
        mean_1982_2002 = to_annual(cropped_series, 'mean')
        # mask = out_series.index.values > dt.datetime(2003, 1, 1, 0, 0, 0)
        # mask1 = out_series.index.values < dt.datetime(2023, 1, 1, 0, 0, 0)
        cropped_series = out_series.query("date >= '2003-01-01' and date < '2023-01-01'")
        mean_2003_2022 = to_annual(cropped_series, 'mean')
        a = q90 - mean_1982_2012
        b = mean_1982_2012 - q10
        a = a.rolling(31, min_periods=1, center=True).mean()
        b = b.rolling(31, min_periods=1, center=True).mean()
        q90 = mean_1982_2012 + a
        q10 = mean_1982_2012 - b
        # a = a.rolling(10, min_periods=1, center=True).mean()
        # b = b.rolling(10, min_periods=1, center=True).mean()

        yr_st = 1982
        yr_end = 2025
        #clim = 'CLIM_1982-2022'
        # clim = 'CLIM_2000-2019'
        clim = 'CLIM_1982-2016'

        # days per month  - we ignore Feb 29th
        dpm = [31,28,31,30,31,30,31,31,30,31,30,31]

        # set months to plot (1=January)
        #mn_st, mn_end = 3, 9
        mn_st, mn_end = 1, 12
        date_start = datetime(2021,mn_st,1)
        date_end = date_start + relativedelta(months=mn_end-mn_st+1)
        dt_list = [date_start]
        while dt_list[-1] < date_end:
            dt_list.append(dt_list[-1]+timedelta(days=1))
        dt_list = dt_list[:-1]
        day0 = int(np.sum(dpm[0:mn_st-1]))
        day1 = day0 + len(dt_list)
        # print(f'{day0=} {day1=}')

        # remove Year from the x-axis label so we can overplot multiple years
        month_day_fmt = mdates.DateFormatter('%b') # B for full month name, b for abbreviated

        # set up axes

        # gather the data
        reg_means = {}
        for year in list(range(yr_st,yr_end+1)):  # + [clim]:
        #   print(f'processing year {year}')
        #     if year==2007:
        #         file = f'{DIR}/regmeans_{year}.nc'
        #         nc = xr.open_dataset(file)
            # mask = out_series.index.values >= dt.datetime(year, 1, 1, 0, 0, 0)
            # mask1 = out_series.index.values < dt.datetime(year+1, 1, 1, 0, 0, 0)
            reg_means[year] = out_series.query("date >= '{}-01-01' and date < '{}-01-01'".format(year, year+1))

        # set colour for each year
        normal_year_colour = 'dimgrey' #'#dbdbdb' # grey
        colors = {}
        for i in range(yr_st,yr_end+1):
            colors[i] =  normal_year_colour
        # colors[2023] = 'orange' # red #ff0000
        # colors[2022] = '' # blue
        # colors[2024] = 'r'  # red #ff0000
        colors[2024] = 'b' #'#ff7f0e' # orange
        colors[2025] = 'r'  # '#ff7f0e' # orange
        # colors[2017] = '#2ca02c' # green
        # colors[2008] = '#d62728' # brown
        colors[clim] = 'k' # black #000000

        # plot the data for each year
        for year in list(range(1982,2026)):
            color = colors[year]
            lw = 0.5
            label = None
            if color != normal_year_colour:
                lw = 1.3
                label = f'{year}'
            if year == 2025:
                lw = 2.0
            # if year == 2024:
            #     lw = 2.0
            if year == clim:
                lw = 2.0
                label = clim[5:] + ' climatology'
            if year == 2021:
                label = f'{yr_st} to 2024'
            ndays = len(reg_means[year]) - day0
            values = reg_means[year]
            ax.plot(dt_list[:ndays],values[day0:day1],color=color,label=label,lw=lw)

            # if year == 2023: # add June 2023-last 20-year anomaly
            #     print('day0, day1', day0, day1)
            #     print('plotting removed background trend')
            #     print('val', values[151:212])
            #     print('2003-2022 val', np.array(mean_2003_2022[151:212].values)[:, 0])
            #     print('1982-2002 val', np.array(mean_1982_2002[151:212].values)[:, 0])
            #     print(dt_list[151:212])
                # ax.plot(dt_list[151:190], values[151:190]-np.array(mean_2003_2022[151:190].values)[:, 0]+np.array(mean_1982_2002[151:190].values)[:,0],
                #         color=color, linestyle='--', label='June 2023-last 20-year anomaly', lw=lw)

        ndays = len(reg_means[2007]) - day0

        if 'cdr' in region_wak:
            ax.plot(dt_list[:ndays], mean_1982_2012.values, 'k', label='1982-2012 clim')
            # ax.plot(dt_list[:ndays], mean_1982_2002.values, 'k-', label='1982-2002')
            # ax.plot(dt_list[:ndays], mean_2003_2022.values, 'k:', label='2003-2022')
            ax.fill_between(dt_list[:ndays], q10.values[:, 0] + q10.values[:, 0] - mean_1982_2012.values[:, 0],
                                                     q90.values[:, 0] + q90.values[:, 0] - mean_1982_2012.values[:, 0],
                                                     color='grey', alpha=0.7, label = 'cat I')
            ax.fill_between(dt_list[:ndays], q10.values[:, 0] + 2*q10.values[:, 0] - 2*mean_1982_2012.values[:, 0],
                                                     q90.values[:, 0] + 2*q90.values[:, 0] - 2*mean_1982_2012.values[:, 0],
                                                     color='lightgrey', alpha=0.5, label = 'cat II')
            ax.fill_between(dt_list[:ndays], q10.values[:, 0] + 3*q10.values[:, 0] - 3*mean_1982_2012.values[:, 0],
                                                     q90.values[:, 0] + 3*q90.values[:, 0] - 3*mean_1982_2012.values[:, 0],
                                                     color='lightgrey', alpha=0.3, label = 'cat III')
            ax.fill_between(dt_list[:ndays], q10.values[:, 0], q90.values[:, 0], color='k', alpha=0.4,
                            label='10$^{th}$-90$^{th}$')

        # print('dt_list: ', dt_list)
    ax.legend(loc='upper left', fontsize=10)
    ax.set_ylabel('SST (deg C)')
    # ax.set_title(f'mean OSTIA SST {nc.region_names.values[region].decode()}')
    ax.xaxis.set_major_formatter(month_day_fmt)
    ax.margins(x=0)
    ax.grid(True, linestyle=':')
    # # 10/06 - 18/06
    # plt.axvline(x=dt_list[156], color='k', linestyle='--')
    # # plt.axvline(x=dt_list[163], color='k', linestyle='--')
    # plt.axvline(x=dt_list[168], color='k', linestyle='--')
    #
    #
    # plt.axvline(x=dt_list[204], color='k', linestyle='--')
    # # plt.axvline(x=dt_list[209], color='k', linestyle='--')  # 29/07
    # plt.axvline(x=dt_list[214], color='k', linestyle='--')
    #
    # plt.axvline(x=dt_list[256], color='k', linestyle='--')  #
    # # plt.axvline(x=dt_list[261], color='k', linestyle='--')  # 19/09
    # plt.axvline(x=dt_list[266], color='k', linestyle='--')  #
    # remove top and right spines from plot
    # sns.despine()

    plt.tight_layout()
    #plt.show()
    # save the plot
    # if region_wak < 0:
    #     region_wak = 'SHELF'
    # plt.show()
    plt.savefig('/home/users/simon.brown/extremes/heatwaves/mhw/OSTIA_region_WAKELIN_{}_maremon.png'.format(region_wak), dpi=300)
    plt.close()


# Monthly mean record calculation:
# a_month = out_tseries.resample('M').mean()
# aa['month'] = aa.index.month
# a_may = a_month.query("month == 5")
# print('Maximum May temperature:', a_may['ostia_cdr_NWS_sst'].values.max())
# # second maximum
# print('Second maximum May temperature:', a_may.query("date < '2024-05-31'")['ostia_cdr_NWS_sst'].values.max())