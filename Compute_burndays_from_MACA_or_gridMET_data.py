# Load packages
import netCDF4 as nc
import numpy as np
import os
import pandas as pd
import geopandas as gpd
import rasterio 
import rasterio.mask
from haversine import haversine
import itertools

#Climate model controls
climate_flag = 2 #if 1-> use MACA data, if 2 -> use GridMET data
#c = ['future']
c = ['present']

#maca_model_name = ['CanESM2', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'GFDL-ESM2G', 'GFDL-ESM2M', 'HadGEM2-ES365', 'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'MIROC5', 'MIROC-ESM', 'MRI-CGCM3']
maca_model_name = ['GridMET']

#Load ecoregion shapefile
ecoregion_shape_pf = 'ECO_L2_Cont_States.shp'

#Load LANDFIRE fuels data
lffms = 'LF2020_FBFM40_200_CONUS_coarsedata.nc'
lffms_csv = 'LF16_F40_200.csv'

#Load prescription data
burnRx = 'Rx_windows_github.xlsx'

#output path for saving files at end of run
outpath = ''

#open ecoregion shapefile
eco_shapefile = gpd.read_file(ecoregion_shape_pf)
eco_shapefile = eco_shapefile.to_crs(4326) #switch coordinate system to LANDFIRE

#read in prescription data
df_raw = pd.read_excel(burnRx)
df = df_raw
df['latitude']                = df['latitude'].astype(float)
df['longitude']               = df['longitude'].astype(float)
df['Fuels']                   = df['Fuels (Scott&Burgan)']
df['max temp (K)']                = df['max temp (K)'].astype(float)
df['min temp (K)']            = df['min temp (K)'].astype(float)
df['min RH (%)']              = df['min RH (%)'].astype(float)
df['max WS (m/s)']                  = df['max WS (m/s)'].astype(float)
df['min WS (m/s)']            = df['min WS (m/s)'].astype(float)

df_fuelmap = df

#create full geodataframe
df_fuelmap['geometry'] = df_fuelmap.apply(lambda x: point((float(x.longitude), float(x.latitude))), axis=1) 
point = gpd.GeoDataFrame(df_fuelmap, geometry='geometry',crs=4326) 

#intersect sites with ecoregion shapefile
matches = gpd.overlay(point,eco_shapefile,how='intersection')

#compute burn days for specified time frame(s)
for climate_timeframe in c:
    if (c=='future')&(climate_flag==2):
        print('No gridmet future data, exiting now....')
        exit()

#Loop over MACA models (or use gridMET)
    for mmn in range(len(maca_model_name)):
        if climate_timeframe=='present':
            timeframe = ['2006','2015']
        else:
            timeframe = ['2051','2060']

        if climate_flag==1:
            climate_pf = '/Volumes/Project_Data/RxFire/climate_models/MACA/'
            climate_model_only = maca_model_name[mmn] #change the 0 to another index in maca_model_name to change models
            climate_model = climate_model_only+'_'+timeframe[0]+'_'+timeframe[1]+'_CONUS_10yr_avg.nc'

        if climate_flag==2:
            climate_pf = '/Volumes/Project_Data/RxFire/climate_models/GRIDMET/'
            climate_model_only = 'GridMet' #change the 0 to another index in maca_model_name to change models
            climate_model = climate_model_only+'_10y_avg_orig.nc'

        #MACA path names
        print(climate_model)
        maca_ta_min = climate_pf+'tasmin_'+climate_model
        maca_ta_max = climate_pf+'tasmax_'+climate_model
        maca_rh_min = climate_pf+'rhsmin_'+climate_model
        maca_ws     = climate_pf+'ws_'+climate_model

        print('Opening LF data')
        lf_fuel_src   = rasterio.open('netcdf:'+lffms, mode='r')
        lf_fuel_csv   = pd.read_csv(lffms_csv)

        print('Opening MACA data')
        src_maca_ta_min   = rasterio.open('netcdf:'+maca_ta_min, mode='r') 
        src_maca_ta_max   = rasterio.open('netcdf:'+maca_ta_max, mode='r') 
        src_maca_rh_min   = rasterio.open('netcdf:'+maca_rh_min, mode='r')
        src_maca_ws   = rasterio.open('netcdf:'+maca_ws, mode='r')

        #create array to hold burn day results
        fill_raster = np.ones((365,lf_fuel_src.height, lf_fuel_src.width),dtype='int')*9999
        cols, rows = np.meshgrid(np.arange(lf_fuel_src.width), np.arange(lf_fuel_src.height))
        xs, ys = rasterio.transform.xy(lf_fuel_src.transform, rows, cols)
        lons = np.array(xs)
        lats = np.array(ys)
        lon = lons[0,:] #coordinate array 1D
        lat = lats[:,0] #coordinate array 1D
       
        cnt = 0 
        for ec in eco_shapefile['NA_L2CODE'].unique():
            
            eco_choice = eco_shapefile[eco_shapefile['NA_L2CODE']==ec]            
            
            #select all sites based on ecoregion
            match_sites = matches.loc[matches['NA_L2CODE']==ec] 
            print('Ecoregion: ',ec)

            #Loop over ecoregions that have data
            if not match_sites.empty:
                lf_fuel, lf_fuel_transform                     = rasterio.mask.mask(lf_fuel_src,eco_choice['geometry'])
                maca_data_ta_min, maca_data_ta_transform_min   = rasterio.mask.mask(src_maca_ta_min,eco_choice['geometry'])
                maca_data_ta_max, maca_data_ta_transform_max   = rasterio.mask.mask(src_maca_ta_max,eco_choice['geometry'])
                maca_data_rh_min, maca_data_rh_transform_min   = rasterio.mask.mask(src_maca_rh_min,eco_choice['geometry'])
                maca_data_ws, maca_data_ws_transform           = rasterio.mask.mask(src_maca_ws,eco_choice['geometry']    )

                match_sites['Fuels'] = [f for f in match_sites['Fuels'].str.split(",").replace(" ", "")]
                fms = match_sites['Fuels']
                print(fms)
                fms = set(itertools.chain.from_iterable(fms))
                fms = list(map(lambda x: x.replace(" ", ""), list(fms)))
                fm_val = lf_fuel_csv['VALUE'][lf_fuel_csv['FBFM40'].str.contains('|'.join(fms))].to_numpy() 

                lf_fuel = lf_fuel.reshape((lf_fuel.shape[1], lf_fuel.shape[2]))
                lf_points = np.argwhere(np.isin(lf_fuel,fm_val)).T

                if np.size(lf_points[0]):
                    for x,y in zip(lf_points[0], lf_points[1]):
                        #get coordinates of x,y climate block
                        lon_c, lat_c = rasterio.transform.xy(lf_fuel_transform, x,y)
                        fuel_label = lf_fuel_csv['FBFM40'][lf_fuel_csv['VALUE']==lf_fuel[x,y]].to_list()
                        matches_similarfuel = pd.DataFrame()
                        for ms in match_sites.index:
                            match_list = match_sites.loc[ms]['Fuels']
                            match_list = [k.replace(" ", "") for k in match_list]
                            if (fuel_label[0] in match_list):
                                matches_similarfuel = matches_similarfuel.append(match_sites.loc[ms])

                        prox_idx = 99999999
                        #given a climate box, find the closest site with prescription data based on lat/lons
                        for xx in matches_similarfuel.index:
                            prox_idx1 = prox_idx
                            prox_idx = min(prox_idx1,haversine( (lat_c,lon_c) , ( matches_similarfuel.loc[xx]['latitude'], matches_similarfuel.loc[xx]['longitude']) ))
                            if prox_idx< prox_idx1:
                                xxx = xx

                        prox_match = matches_similarfuel.loc[xxx]
                        match_sitename = prox_match['site name']

                        #if fuel is in block, do climate analysis
                        #determine burn days based on site info/find the points in the climate data that match up with the LF points
                        maxtemp    = prox_match['max temp (K)']    
                        mintemp    = prox_match['min temp (K)']      
                        minrhs     = prox_match['min RH (%)']    
                        maxwindspd = prox_match['max WS (m/s)']  
                        minwindspd = prox_match['min WS (m/s)']

                        #create index of days
                        daysindx = (np.arange(1,365)).astype(int)
                        
                        #find days in maca data that fall within their respective prescriptions
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


                #delete rasters for faster processing
                del lf_fuel
                del maca_data_ta_min, maca_data_ta_transform_min
                del maca_data_ta_max, maca_data_ta_transform_max
                del maca_data_rh_min, maca_data_rh_transform_min
                del maca_data_ws, maca_data_ws_transform
            
            else:
                print('No sites found in ecoregion ',ec)
                
            cnt+=1
            print('Percent Finished: ', int((cnt/len(eco_shapefile['NA_L2CODE'].unique()))*100),'%')
            
            
        #close files
        src_maca_ta_min.close()
        src_maca_ta_max.close()
        src_maca_rh_min.close()
        src_maca_ws.close()
        lf_fuel_src.close()

        fname = outpath+climate_timeframe+'/'+climate_model_only+'_'+climate_timeframe+'_burndays.nc'
        if os.path.isfile(fname):
            os.remove(fname)

        #save burn days
        print('saving burn days '+climate_model_only)
        burndays_netfile = nc.Dataset(fname, "w")
        lat_dim = burndays_netfile.createDimension('lat', len(lat))
        lon_dim = burndays_netfile.createDimension('lon', len(lon))
        time_dim = burndays_netfile.createDimension('time', 365)
        lat_new = burndays_netfile.createVariable('lat', np.float32, ('lat',))
        lat_new[:] = lat
        lon_new = burndays_netfile.createVariable('lon', np.float32, ('lon',))
        lon_new[:] = lon
        time_new = burndays_netfile.createVariable('time', np.int64, ('time',))
        time_new[:] = np.array(np.arange(365))
        days=fill_raster
        avg_var = burndays_netfile.createVariable('Days',np.int64,('time','lat','lon'), fill_value=9999) 
        avg_var[:,:,:] = days
        print('closing new NC file')
        burndays_netfile.close()

        del fill_raster
