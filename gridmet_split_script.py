# Libraries
import pickle
import time
import geopandas as gpd
import geopandas.geodataframe
import grd2shp_xagg
import xagg as xa
import xarray as xr
import os
import datetime

# Function definitions
def get_gridmet_datasets(variable, start_date, end_date, polygon_for_bbox = None, lon_min = None, lat_min = None, lon_max = None, lat_max = None):

    """
    :param str/list variable: data variable short name or list of data variables in short name. Must be one of the following: ['tmmx', 'tmmn', 'pr', 'srad', 'vs', 'rmax', 'rmin', 'sph']
    :param str start_date: Start date of data collection yyyy-mm-dd
    :param str end_date: End date of data collection yyyy-mm-dd
    :param gpd.GeoDataFrame or str polygon_for_bbox: either a geodataframe of the polygons you are aggregating to or a path to a shapefile (or other geo file) that can be read into a geodataframe
    :param str lon_min: bbox of aoi longitude min. Not used if polygon_for_bbox is given (i.e. polygon_bbox != None)
    :param str lat_min: bbox of aoi latitude min. Not used if polygon_for_bbox is given (i.e. polygon_bbox != None)
    :param str lon_max: bbox of aoi longitude max. Not used if polygon_for_bbox is given (i.e. polygon_bbox != None)
    :param str lat_max: bbox of aoi latitude max. Not used if polygon_for_bbox is given (i.e. polygon_bbox != None)
    :return: dictionary of xarray stored by variable name
    """

    ## check/define bounds for data slicing
    if polygon_for_bbox is not None:
        if isinstance(polygon_for_bbox, geopandas.geodataframe.GeoDataFrame):
            print('polygon is geodataframe')
            print(polygon_for_bbox.total_bounds)
            pass
        elif os.path.isfile(polygon_for_bbox):
            print('polygon is path to shapefile')
            polygon_for_bbox = gpd.read_file(polygon_for_bbox)
            print(polygon_for_bbox.total_bounds)
        else:
             raise TypeError('polygon_for_bbox should str path to shapefile or geodataframe')

        if polygon_for_bbox.crs != 'EPSG:4326':
            polygon_for_bbox = polygon_for_bbox.to_crs('EPSG:4326')
            print('geodataframe crs changed to EPSG:4326')

        lon_min, lat_min, lon_max, lat_max = polygon_for_bbox.total_bounds

    ## if only 1 var (short) is provided, change to list
    if not isinstance(variable, list):
        variable = [variable]

    ## Initiate list of empty datasets
    xarray_dict = dict()

    ## Loop through variables and pull data + place in dict
    for var in variable:
        ## Pulling data
        # source: http://thredds.northwestknowledge.net:8080/thredds/reacch_climate_MET_aggregated_catalog.html
        url = f'http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_met_{var}_1979_CurrentYear_CONUS.nc'
        # call data from url
        start = time.perf_counter()
        ds = xr.open_dataset(url + '#fillmismatch')

        ## Subset to timeframe and bbox of interest:
        ds_subset = ds.sel(day = slice(start_date, end_date),
                           lon = slice(lon_min, lon_max),
                           lat = slice(lat_max, lat_min))
        end = time.perf_counter()
        print(f'finish in {round(end - start, 2)} seconds')

        # Append to dict of xr.datasets
        xarray_dict[var] = ds_subset

    return xarray_dict

def create_weightmap(xarray_dict, polygon, output_data_folder, weightmap_var = None):

    """
    :param dict xarray_dict: dictionary of gridmet data. the output of get_gridmet_datasets()
    :param gpd.GeoDataFrame polygon: geodataframe polygon or multipolygon that is used for regridding
    :param str output_data_folder: path to folder for function output
    :param str weightmap_var: variable used to make weightfile + naming output e.i. weightmap_var = 'tmin'. If none, first var of dict  will be used.
    :return: output path to weightmap pickle file
    """
    if isinstance(polygon, geopandas.geodataframe.GeoDataFrame):
        print('polygon is geodataframe')
        pass
    elif os.path.isfile(polygon):
        print('polygon is path to shapefile')
        polygon = gpd.read_file(polygon)
    else:
        raise TypeError('polygon should str path to shapefile or geodataframe')

    ## Name output weightmap_file
    if weightmap_var is None:
        # If no var is specified, first var of xarray_dict will be used for generating weight file, and name of output will be generic: 'grd2shp_weights.pickle'
        weightmap_file = os.path.join(output_data_folder, 'grd2shp_weights.pickle')
        weightmap_var = list(xarray_dict.keys())[0]
    else:
        weightmap_file = os.path.join(output_data_folder, f'grd2shp_weights_{weightmap_var}.pickle')


    ## Produce weightmap
    start = time.perf_counter()
    weightmap = xa.pixel_overlaps(xarray_dict[weightmap_var], polygon)
    end = time.perf_counter()
    print(f'finished weightmap calculation in {round(end - start, 2)} second(s)')

    with open(weightmap_file, 'wb') as file:
        pickle.dump(weightmap, file)

    return weightmap_file

def g2shp_regridding(xarray_dict, polygon, weightmap_file, out_file):

    """
    :param dict xarray_dict: dictionary of gridmet data. the output of get_gridmet_datasets()
    :param gpd.GeoDataFrame polygon: geodataframe polygon or multipolygon that is used for regridding
    :paran str weightmap_file: path to weight file for redridding. Output of create_weightmap()
    :param str output_data_folder: path to folder for function output
    :return: regriddied g2s object, which is saved as .nc to output_data_folder
    """

    ## params for g2s.initialise()
    # Lsit of xarray.datasets
    vars_grd_list = list(xarray_dict.values())
    
    combined_ds = xr.merge(vars_grd_list)
    with open(weightmap_file, "rb") as file:
        weightmap = pickle.load(file) 

    # Run regridding
    start = time.perf_counter()
    ds_agg = xa.aggregate(combined_ds, weightmap)
    end = time.perf_counter()
    print('finished agg in {} second(s)'.format(round(end - start, 2)))

    ## Saving output
    ds_agg.to_netcdf(out_file)
    return ds_agg

# Variables and run functions
if __name__ =='__main__':

    ## Variable definitions
    ### official list of variables needed for drb-inland-salinity model
    data_vars_shrt_all = ['tmmx', 'tmmn', 'pr', 'srad', 'vs','rmax','rmin','sph']
    data_vars_shrt_all = data_vars_shrt_all[:2]
    ### drb catchment polygons
    gdf_nhru02_path = './data/drb_prms_basins_fixed_from_nhru02.shp'
    gdf = gpd.read_file(gdf_nhru02_path)
    ### date range
    start_date = '1979-01-01'
    end_date = '1985-01-01'
    ### output folder
    output_path = './data/'

    xarray_dict = get_gridmet_datasets(variable = data_vars_shrt_all,
                                       start_date= start_date, end_date = end_date,
                                       polygon_for_bbox = gdf)
                                      # lat_min=round(lat_min,2), lon_min=round(lon_min,2),
                                      # lat_max=round(lat_max,2), lon_max=round(lon_max,2))

    weight_map_path = create_weightmap(xarray_dict = xarray_dict,
                     polygon=gdf,
                     output_data_folder = output_path,
                     weightmap_var = 'tmmx')

    ### Subset for streamlined testing
    subset = {key: xarray_dict[key] for key in ['tmmx']}

    start = time.perf_counter()
    g2shp_regridding(xarray_dict= xarray_dict,
                     polygon=gdf,
                     weightmap_file= weight_map_path,
                     out_file = "test_out.nc")
    end = time.perf_counter()
    print("total regrid time: ", end-start)
