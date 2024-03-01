'''
## Burn_day_calculation.py

This code computes prescribed fire burn days 

## Authors: Julia Oliveto (joliveto@lanl.gov)
            Alex Jonko (ajonko@lanl.gov)
            Teresa Beaty (teresa.beaty@nnmc.edu)

## Reference: Jonko, A., J. Oliveto, T. Beaty, A. Atchley, M. A. Battaglia, 
              M.B. Dickinson, M.R. Gallagher, A. Gilbert, D. Godwin, J.A. Kupfer, 
              J.K. Hiers, C. Hoffman, M. North, J. Restaino, C. Sieg, and N. Skowronski 
              "How will future climate change impact prescribed fire across the 
              continguous United States?" npj Climate and Atmospheric Science (submitted).

© 2024. Triad National Security, LLC. All rights reserved. This program was produced under 
U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), 
which is operated by Triad National Security, LLC for the U.S. Department of Energy/National 
Nuclear Security Administration. All rights in the program are reserved by Triad National 
Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, 
irrevocable worldwide license in this material to reproduce, prepare derivative works, 
distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
'''

### Load python modules ###
import netCDF4 as nc
import numpy as np
import os
import pandas as pd
import geopandas as gpd
import rasterio 
import rasterio.mask
from shapely.geometry import *
from haversine import haversine
import itertools
from tqdm import tqdm
import netCDF4
import dateutil.parser

print('=== RUNNING BURN_DAY_CALCULATION.PY ===')

climate_timeframe  = ['2006_2010','2011-2015','2051-2055','2056-2069']
model_name         = ['bcc-csm1-1-m','bcc-csm1-1','BNU-ESM','CanESM2','CSIRO-Mk3-6-0','CNRM-CM5','GFDL-ESM2G','GFDL-ESM2M','HadGEM2-CC365','HadGEM2-ES365','inmcm4','IPSL-CM5A-LR','IPSL-CM5A-MR','IPSL-CM5B-LR','MIROC5','MIROC-ESM','MIROC-ESM-CHEM','MRI-CGCM3']

### Specify Data Loactions ###
ecoregion_shape_pf = './ECO_L2_Cont_States.shp'
lffms              = './LF2020_FBFM40_200_CONUS_coarsedata.nc'
lffms_csv          = './LF16_F40_200.csv'
burnRx             = './Rx_windows_github_Rev2.csv'
outpath            = #SPECIFY AN OUTPUT DIRECTORY
climate_pf         = #SPECIFY MACA INPUT FILE LOCATION. DATA SHOULD HAVE LEAP YEARS REMOVED

### Open ecoregion shapefile ###
eco_shapefile      = gpd.read_file(ecoregion_shape_pf)
eco_shapefile      = eco_shapefile.to_crs(4326) #switch coordinate system to LANDFIRE

### Read Rx window data ###
df = pd.read_csv(burnRx)
df['latitude']    = df['lat'].astype(float)
df['longitude']   = df['lon'].astype(float)
df['Fuels']       = df['fuels']
df['maxT']        = df['maxT'].astype(float)
df['minT']        = df['minT'].astype(float)
df['maxRH']       = df['maxRH'].astype(float)
df['minRH']       = df['minRH'].astype(float)
df['maxWS']       = df['maxWS'].astype(float)
df['minWS']       = df['minWS'].astype(float)
df['start1']      = df['start1'].astype(int)
df['end1']        = df['end1'].astype(int)
df['start2']      = df['start2'].astype(float)
df['end2']        = df['end2'].astype(float)

### Create full geodataframe ###
df_fuelmap = df
df_fuelmap['geometry'] = df_fuelmap.apply(lambda x: Point((float(x.longitude), float(x.latitude))), axis=1) 
point = gpd.GeoDataFrame(df_fuelmap, geometry='geometry',crs=4326) 

### Intersect Rx sites with ecoshapefile ###
matches = gpd.overlay(point,eco_shapefile,how='intersection') 

for timeframe in climate_timeframe:
    for mmn in range(len(model_name)):
        climate_model_only = model_name[mmn] 
        climate_model = climate_model_only+'_rcp45_'+timeframe+'.nc'

        print('Model: ',climate_model)
        maca_ta_min = climate_pf+'macav2metdata_tasmin_'+climate_model
        maca_ta_max = climate_pf+'macav2metdata_tasmax_'+climate_model
        maca_rh_min = climate_pf+'macav2metdata_rhsmin_'+climate_model
        maca_ws     = climate_pf+'macav2metdata_ws_'+climate_model

        print('Opening LANDFIRE fuel data')
        lf_fuel_src   = rasterio.open('netcdf:'+lffms, mode='r')
        lf_fuel_csv   = pd.read_csv(lffms_csv)

        print('Opening climate model data')
        src_maca_ta_min   = netCDF4.Dataset(maca_ta_min, 'r')
        src_maca_ta_max   = netCDF4.Dataset(maca_ta_max, 'r') 
        src_maca_rh_min   = netCDF4.Dataset(maca_rh_min, 'r')
        src_maca_ws       = netCDF4.Dataset(maca_ws    ,'r')

### Create array where burn days will be saved ###
        final_raster = np.ones((5,lf_fuel_src.height, lf_fuel_src.width),dtype='int')*9999
        cols, rows = np.meshgrid(np.arange(lf_fuel_src.width), np.arange(lf_fuel_src.height))
        xs, ys = rasterio.transform.xy(lf_fuel_src.transform, rows, cols)
        lons = np.array(xs)
        lats = np.array(ys)
        lon = lons[0,:] 
        lat = lats[:,0] 
       
        cnt = 0 
        all_times = src_maca_ta_min.variables['time']
        for ttime in range(int(timeframe.split('_')[0]), int(timeframe.split('_')[1])+1):
            print('Running Year: ',ttime)
            fill_raster = np.ones((365,lf_fuel_src.height, lf_fuel_src.width),dtype='int')*9999
            sdt = dateutil.parser.parse(str(ttime)+"-01-01")
            edt = dateutil.parser.parse(str(ttime)+"-12-31")
            st_idx = netCDF4.date2index(sdt, all_times)
            et_idx = netCDF4.date2index(edt, all_times)

            maca_data_ta_min = src_maca_ta_min.variables['air_temperature'][st_idx:et_idx+1,:,:]
            maca_data_ta_max = src_maca_ta_max.variables['air_temperature'][st_idx:et_idx+1,:,:]
            maca_data_rh_min = src_maca_rh_min.variables['relative_humidity'][st_idx:et_idx+1,:,:]
            maca_data_ws     = src_maca_ws.variables['wind_speed'][st_idx:et_idx+1,:,:]
            maca_data_ws     = np.flip(maca_data_ws,axis=1)
            maca_data_rh_min = np.flip(maca_data_rh_min,axis=1)
            maca_data_ta_max = np.flip(maca_data_ta_max,axis=1)
            maca_data_ta_min = np.flip(maca_data_ta_min,axis=1)
         
### Select all cells that match ecoregion ###
            for ec in tqdm(eco_shapefile['NA_L2CODE'].unique()):
                
                eco_choice = eco_shapefile[eco_shapefile['NA_L2CODE']==ec]   
                match_sites = matches.loc[matches['NA_L2CODE']==ec] 
                print('Ecoregion: ',ec)

                if not match_sites.empty:
                    lf_fuel, lf_fuel_transform = rasterio.mask.mask(lf_fuel_src,eco_choice['geometry'], crop=False)
                    print(np.shape(maca_data_ta_min), np.shape(lf_fuel))

                    match_sites['Fuels'] = [f for f in match_sites['Fuels'].str.split(",").replace(" ", "")]
                    fms = match_sites['Fuels']
                    print(fms)
                    fms = set(itertools.chain.from_iterable(fms))
                    fms = list(map(lambda x: x.replace(" ", ""), list(fms)))
                    fm_val = lf_fuel_csv['VALUE'][lf_fuel_csv['FBFM40'].str.contains('|'.join(fms))].to_numpy() 
                    lf_fuel = lf_fuel.reshape((lf_fuel.shape[1], lf_fuel.shape[2]))
                    lf_points = np.argwhere(np.isin(lf_fuel,fm_val)).T

                    if np.size(lf_points[0]):
                        print('Cycling through ', np.size(lf_points[0]),' data points')
                        for x,y in tqdm(zip(lf_points[0], lf_points[1])):

                            lon_c, lat_c = rasterio.transform.xy(lf_fuel_transform, x,y)
                            fuel_label = lf_fuel_csv['FBFM40'][lf_fuel_csv['VALUE']==lf_fuel[x,y]].to_list()
                            matches_similarfuel = pd.DataFrame()
                            for ms in match_sites.index:
                                match_list = match_sites.loc[ms]['Fuels']
                                match_list = [k.replace(" ", "") for k in match_list]
                                if (fuel_label[0] in match_list):
                                    matches_similarfuel = matches_similarfuel.append(match_sites.loc[ms])

### For each grid cell within the ecoregion, find the closest Rx location based on lat/lon coordinates using the haversine function ###
                            prox_idx = 999999
                            xxxx = []
                            for xx in matches_similarfuel.index:
                                prox_idx1 = prox_idx
                                prox_idx = min(prox_idx1,haversine( (lat_c,lon_c) , ( matches_similarfuel.loc[xx]['latitude'], matches_similarfuel.loc[xx]['longitude']) ))
                                if prox_idx< prox_idx1:
                                    xxx = xx
                            xxxx.append(xxx)
                            for xx in matches_similarfuel.index:
                                if ((matches_similarfuel.loc[xxx]['latitude']==matches_similarfuel.loc[xx]['latitude']) &
                                    (matches_similarfuel.loc[xxx]['longitude']==matches_similarfuel.loc[xx]['longitude']) &
                                    (xx!=xxx)):
                                    xxxx.append(xx)
                            
                            for xloc in xxxx:
                                prox_match = matches_similarfuel.loc[xloc]
                                match_sitename = prox_match['site']                                

### Determine burn days based on prescription ranges for closest prescription within ecoregion ###
                                maxtemp    = prox_match['maxT']      
                                mintemp    = prox_match['minT']  
                                maxrhs     = prox_match['maxRH']       
                                minrhs     = prox_match['minRH']    
                                maxwindspd = prox_match['maxWS']   
                                minwindspd = prox_match['minWS']
                                startday1  = prox_match['start1']
                                startday2  = prox_match['start2']
                                endday1    = prox_match['end1']
                                endday2    = prox_match['end2']

### Several prescriptions have different ranges depending on season, they have two start and end days to discriminate between the seasons ###
                                if np.isnan(startday2):
                                    startday2 = startday1 
                                    endday2   = endday1
                                else:
                                    startday2 = startday2.astype(int)
                                    endday2   = endday2.astype(int)

                                daysindx = np.unique(np.concatenate((np.arange(startday1,endday1),np.arange(startday2,endday2)), axis=0)).astype(int)
                                nr_days_pos_all = np.nonzero((maca_data_ta_min[daysindx,x,y]>mintemp) &
                                                            ( maca_data_ta_max[daysindx,x,y]<maxtemp) &  
                                                            ( maca_data_rh_min[daysindx,x,y]>minrhs) &  
                                                            (  maca_data_ws[daysindx,x,y]>minwindspd) & 
                                                            (  maca_data_ws[daysindx,x,y]<maxwindspd) )

                                b = np.abs(lon - lon_c).argmin().astype(int); a = (np.abs(lat - lat_c)).argmin().astype(int)
                    
                                if np.count_nonzero(fill_raster[:,a,b]==9999)==len(fill_raster[:,a,b]):
                                    fill_raster[:,a,b] = 0
                                if np.size(nr_days_pos_all[0]):
                                    fill_raster[nr_days_pos_all[0],a,b] = 1
                            
                else:
                    print('No grid cells found in ecoregion ',ec)
                    
                cnt+=1
                print('Percent Finished: ', int((cnt/len(eco_shapefile['NA_L2CODE'].unique()))*100),'%')
            pre_final = np.sum(fill_raster, axis=0)
            pre_final[pre_final>365] = 9999
            final_raster[int(ttime%int(timeframe.split('_')[0])),:,:] = pre_final

        del maca_data_ta_min
        del maca_data_ta_max
        del maca_data_rh_min
        del maca_data_ws
        del lf_fuel

        src_maca_ta_min.close()
        src_maca_ta_max.close()
        src_maca_rh_min.close()
        src_maca_ws.close()
        lf_fuel_src.close()

### Write output file ###
        fname = outpath+climate_model_only+'_'+timeframe+'_burn_days.nc'
        if os.path.isfile(fname):
            os.remove(fname)

        print('Saving burn day file for '+climate_model_only,':')
        print(fname)
        outputfile        = nc.Dataset(fname, "w")
        lat_dim           = outputfile.createDimension('lat', len(lat))
        lon_dim           = outputfile.createDimension('lon', len(lon))
        time_dim          = outputfile.createDimension('time', 5)
        lat_new           = outputfile.createVariable('lat', np.float32, ('lat',))
        lat_new[:]        = lat
        lon_new           = outputfile.createVariable('lon', np.float32, ('lon',))
        lon_new[:]        = lon
        time_new          = outputfile.createVariable('time', np.int64, ('time',))
        time_new[:]       = np.array(np.arange(int(timeframe.split('_')[0]), int(timeframe.split('_')[1])+1))
        days              = final_raster
        avg_var           = outputfile.createVariable('Days',np.int64,('time','lat','lon'), fill_value=9999) 
        avg_var[:,:,:]    = days
        outputfile.close()

        del final_raster

