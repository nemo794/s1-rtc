from zipfile import ZipFile
import os
import math
import ruamel.yaml
from datetime import datetime
from maap.maap import MAAP
import subprocess
import glob

import sys
import argparse
import subprocess
import numpy as np

from maap.Result import Result
from maap.dps.DpsHelper import DpsHelper
from urllib.parse import urlparse

import eof


def build_rtc_runcnfg(user_rncfg, input_s1_zip, orbit_file, dem_file, scratch_dir, output_dir):
    
    # print("user_rncfg: ", user_rncfg)
    # print("input_s1_zip: ", input_s1_zip)
    # print("orbit_file: ", orbit_file)
    # print("dem_file: ", dem_file)
    # print("scratch_dir: ", scratch_dir)
    # print("output_dir: ", output_dir)
    
    # Construct path to the RTC default runconfig
    path = np.__file__      # Ex: '/opt/conda/lib/python3.10/site-packages/numpy/__init__.py'
    path = path.split("site-packages")
    path = os.path.join(path[0], "site-packages", "rtc/defaults/rtc_s1.yaml")
    
    # Read defaults yaml into memory
    # with open("RTC/src/rtc/defaults/rtc_s1.yaml", 'r') as defaults_file:
    # with open("/opt/conda/lib/python3.10/site-packages/rtc/defaults/rtc_s1.yaml", 'r') as defaults_file:
    with open(path, 'r') as defaults_file:
        defaults = ruamel.yaml.safe_load(defaults_file)


    # Construct the output filenames
    s1_filename = os.path.basename(input_s1_zip).lstrip(".zip")
    orbit_num = s1_filename.split("_")[-3]
    data_start = s1_filename.split("_")[-5]
    now_dt_string = datetime.now().strftime("%Y%m%dT%H%M%SZ")     # processing datetime
    proc_str = f"{orbit_num}_{data_start}_{now_dt_string}"
    filenaming_convention = "MAAP_S1RTC_{burst_id}_%s" % proc_str
        
    # Update the desired fields
    defaults['runconfig']['groups']['input_file_group']['safe_file_path'] = [input_s1_zip]
    defaults['runconfig']['groups']['input_file_group']['orbit_file_path'] = [orbit_file]
    # defaults['runconfig']['groups']['input_file_group']['burst_id'] = ['t071_151227_iw2']  # Descending over Mt. Wilson
    # defaults['runconfig']['groups']['input_file_group']['burst_id'] = ['t064_135520_iw1']  # Ascending over Mt. Wilson
    defaults['runconfig']['groups']['dynamic_ancillary_file_group']['dem_file'] = dem_file
    
    defaults['runconfig']['groups']['product_group']['product_path'] = output_dir
    defaults['runconfig']['groups']['product_group']['scratch_path'] = scratch_dir
    defaults['runconfig']['groups']['product_group']['output_dir'] = output_dir
    defaults['runconfig']['groups']['product_group']['product_id'] = filenaming_convention
    defaults['runconfig']['groups']['product_group']['save_bursts'] = True
    defaults['runconfig']['groups']['product_group']['save_mosaics'] = False
    defaults['runconfig']['groups']['product_group']['save_secondary_layers_as_hdf5'] = False
    defaults['runconfig']['groups']['product_group']['save_metadata'] = False
    
    defaults['runconfig']['groups']['processing']['polarization'] = 'cross-pol'
    defaults['runconfig']['groups']['processing']['geocoding']['save_incidence_angle'] = False
    defaults['runconfig']['groups']['processing']['geocoding']['save_local_inc_angle'] = False
    defaults['runconfig']['groups']['processing']['geocoding']['save_projection_angle'] = False
    defaults['runconfig']['groups']['processing']['geocoding']['save_rtc_anf_psi'] = False
    defaults['runconfig']['groups']['processing']['geocoding']['save_range_slope'] = False
    defaults['runconfig']['groups']['processing']['geocoding']['save_nlooks'] = False
    defaults['runconfig']['groups']['processing']['geocoding']['save_rtc_anf'] = False
    defaults['runconfig']['groups']['processing']['geocoding']['save_dem'] = False
    defaults['runconfig']['groups']['processing']['geocoding']['save_layover_shadow_mask'] = False

    defaults['runconfig']['groups']['processing']['geocoding']['bursts_geogrid']['output_epsg'] = None
    defaults['runconfig']['groups']['processing']['geocoding']['bursts_geogrid']['x_posting'] = 30
    defaults['runconfig']['groups']['processing']['geocoding']['bursts_geogrid']['y_posting'] = 30
    defaults['runconfig']['groups']['processing']['geocoding']['bursts_geogrid']['x_snap'] = 30
    defaults['runconfig']['groups']['processing']['geocoding']['bursts_geogrid']['y_snap'] = 30

    defaults['runconfig']['groups']['processing']['mosaicking']['mosaic_geogrid']['output_epsg'] = None
    defaults['runconfig']['groups']['processing']['mosaicking']['mosaic_geogrid']['x_posting'] = None
    defaults['runconfig']['groups']['processing']['mosaicking']['mosaic_geogrid']['y_posting'] = None
    defaults['runconfig']['groups']['processing']['mosaicking']['mosaic_geogrid']['x_snap'] = None
    defaults['runconfig']['groups']['processing']['mosaicking']['mosaic_geogrid']['y_snap'] = None

    with open(user_rncfg, 'w') as out_file:
        ruamel.yaml.dump(defaults, out_file, default_flow_style=False)
    

def get_s1_approx_bbox(input_s1_zip, buffer_distance_in_km=30):
    
    # extract lat/lon string
    with ZipFile(input_s1_zip) as z:
        for zf in z.infolist():
            if zf.filename.endswith(".kml"):
                kml_file = zf.filename
                break
        else:
            raise FileNotFoundError(f"No kml file found in: {safe_file}")
        
        with z.open(kml_file, 'r') as kml:
            kml_str = kml.read().decode('utf-8')

    # Example:
    # <coordinates>-119.081970,35.028236 -116.312889,35.430531 -115.988983,33.809345 -118.702644,33.405102</coordinates>

    tag_start = "<coordinates>"
    tag_stop = "</coordinates>"
    lonlat_str = ''.join(kml_str.split(tag_start)[1].split(tag_stop)[0])
    
    lonlat_list = lonlat_str.split(" ")
    lonlat_pairs = []
    for pair in lonlat_list:
        lon, lat = pair.split(",")
        lonlat_pairs.append((float(lon), float(lat)))
        
    # Initialize min and max longitude and latitude variables with the first pair
    min_longitude, min_latitude = lonlat_pairs[0]
    max_longitude, max_latitude = min_longitude, min_latitude

    # Iterate over the remaining pairs to find the min and max longitude and latitude
    for longitude, latitude in lonlat_pairs[1:]:
        min_longitude = min(min_longitude, longitude)
        max_longitude = max(max_longitude, longitude)
        min_latitude = min(min_latitude, latitude)
        max_latitude = max(max_latitude, latitude)
        
    # Add buffer to the min and max longitude and latitude
    # S1 swath is 250km.

    min_longitude -= buffer_distance_in_km / (111.32 * math.cos(math.radians(min_latitude)))
    max_longitude += buffer_distance_in_km / (111.32 * math.cos(math.radians(max_latitude)))
    min_latitude -= buffer_distance_in_km / 111.32
    max_latitude += buffer_distance_in_km / 111.32

    # find N/S/E/W
    north = max_latitude
    south = min_latitude
    west = min_longitude
    east = max_longitude
    
    return north, south, east, west


def get_query_results(bbox, dateRange):
    """
    Query for Sentinel-1 IW_SLC .zip files from ASF.
    
    Uses the maap.py api so that we use MAAP's EarthData credentials
    
    Arguments
    ---------
    
    
    Returns
    -------
    results : 
        List of results returned from a maap.searchGranule(..) call.
    """

    # invoke the MAAP search client
    # As of May 24, 2023, set maap_host="api.ops.maap-project.org". 
    # Dev team thinks there are still issues with using "api.maap-project.org"
    maap = MAAP(maap_host="api.ops.maap-project.org")
    # maap = MAAP()

    # Set the maximum number of results
    MAX_RESULTS = 500

    from time import time
    start = time()
    results = []
    for short_name in ("SENTINEL-1A_SLC", "SENTINEL-1B_SLC"):
        results += maap.searchGranule(cmr_host="cmr.earthdata.nasa.gov",
                                        # concept_id=concept_id,
                                        short_name=short_name,
                                        bounding_box=bbox,
                                        limit=MAX_RESULTS,
                                        temporal=dateRange
                                        )

    return results


def filter_direction(results, direction='both', verbose=False):
    """
    Filter the MAAP searchGranule query results for only ascending
    or descending granules.
    
    Parameters
    ----------
    results : 
        List of results returned from a maap.searchGranule(..) call.
    direction : str
        The direction of the satellite orbit to filter for.
        Input must be one of: One of: 'asc', 'dsc', 'both'.
        If 'both', then no filtering will be applied.
    
    Returns
    -------
    filtered_results :
        a copy of `results`, but filtered to include only `direction` granules.
    """
    
    if direction not in ('asc', 'dsc', 'both'):
        raise ValueError(f"{direction=}, must be one of 'asc', 'dsc', 'both')")

    if direction == 'both':
        return results
    elif direction == 'asc':
        dir = 'ASCENDING'
    else:
        dir = 'DESCENDING'

    # Collect only the ascending tracks
    filtered_results = []
    for result in results:
        tmp = result['Granule']['AdditionalAttributes']['AdditionalAttribute']
        for item in tmp:
            if item['Name'] == 'ASCENDING_DESCENDING':
                if item['Values']['Value'] == dir:
                    if verbose:
                        gran_id = result['Granule']['DataGranule']['ProducerGranuleId']
                        print("Keeping GranuleID: ", gran_id)
                    filtered_results.append(result)

    return filtered_results


def download_s1_zip(maap_item, input_dir):
    
    os.makedirs(input_dir, exist_ok=True)
    
    input_file = maap_item.getData(input_dir)
    
    return input_file


def granule2bursts(input_s1_zip, output_dir="./output"):
    """
    Downloads a granule and corresponding DEM and EOF files,
    and then processes through OPERA RTC.
    
    Parameters
    ----------
    input_s1_zip : str
        Path to the Sentinel-1 safe file to process
    """

    # extract the granule name from the path
    granule_name = os.path.basename(input_s1_zip)
    granule_name = granule_name.split(".")[0]  # remove the ".zip" extension
    
    # Make new directory for processing and outputs for this granule
    output_dir = os.path.join(output_dir, granule_name)
    os.makedirs(output_dir, exist_ok=True)
    
    # Make a scratch_dir for storing temp files
    scratch_dir = os.path.join("scratch", granule_name)
    os.makedirs(scratch_dir, exist_ok=True)
    
    # download the EOF file

    # New way to access S1 EOF files, as of Oct 31, 2023.
    eof_https_link = get_eof_https_link(input_s1_zip)
    orbit_file = download_file_from_https(eof_https_link)

    # # Use subprocess.run() instead of subprocess.Popen() so that we wait for the command to finish.
    # cmd = ["eof", "--force-asf", "--sentinel-file", input_s1_zip, "--save-dir", scratch_dir]
    # result = subprocess.run(cmd, capture_output=True, text=True, check=True)

    # !eof --sentinel-file {input_s1_zip} --save-dir {scratch_dir}
    # orbit_file = glob.glob(os.path.join(scratch_dir, "*.EOF"))[0]
    print("orbit_file: ", orbit_file)
    
    # Get DEM for this granule    
    north, south, east, west = get_s1_approx_bbox(input_s1_zip)
    # "...four latitude and longitude values in the order of [W,S,E,N]"
    granule_bbox = f"{west} {south} {east} {north}"
    dem_file = os.path.join(scratch_dir, "dem.tif")

    # Annoyingly, rasterio cannot find the PROJ_DATA directory when running in the ADE.
    # The workaround in a Jupyter notebook on MAAP is this:
    #     !PROJ_DATA=/opt/conda/share/proj sardem --bbox {granule_bbox} --data-source COP -o {dem_file} --output-format GTiff
    
    # When running from a python module, I can't figure out how to set the environment
    # variable on MAAP and then also have it "stick" for running sardem via subprocess.run().
    
    # So, use the old school os.system(..) calls.
    # NOTE: This path might not work for DPS jobs. It was found via the ADE Terminal: echo $PROJ_DATA
    os.environ['PROJ_DATA'] = '/opt/conda/share/proj'
    os.system(f"sardem --bbox {granule_bbox} --data-source COP -o {dem_file} --output-format GTiff")
    
    # build RTC runcfg
    user_rncfg = os.path.join(output_dir, "rtc_rncfg.yaml")
    build_rtc_runcnfg(user_rncfg, input_s1_zip, orbit_file, dem_file, scratch_dir, output_dir)
    
    # Run RTC
    log_file = os.path.join(output_dir, "rtc_log.txt")

    # Use rtc_s1.py to have each burst processed in parallel, and ISCE run single-threaded
    # !PROJ_DATA=/opt/conda/share/proj python RTC/app/rtc_s1.py {user_rncfg} --log-file {log_file}
    os.system(f"rtc_s1.py {user_rncfg} --log-file {log_file}")
    
    # Use rtc_s1_single_job.py to have each burst processed in serial, but ISCE run multi-threaded
    # !PROJ_DATA=/opt/conda/share/proj python RTC/app/rtc_s1_single_job.py {user_rncfg} --log-file {log_file}
    # os.system(f"rtc_s1_single_job.py {user_rncfg} --log-file {log_file}")
    

def main(bbox, dateRange, direction):
    """
    Run the RTC pipeline as a MAAP DPS job.
    
    Arguments
    ---------
    bbox : str    
        Bounding Box for AOI
           "...four latitude and longitude values in the order of [W,S,E,N]"
        Examples:
            Altadena, CA: '-118.253784,34.138953,-118.078002,34.265594'
            Mt. Wilson: '-118.06817369989554,34.22169196645805,-118.05801589904763,34.2282243912438'
    dateRange : str
        Date range to query
        Example: '2020-08-31T00:00:00Z,2020-09-16T23:59:59Z'
    direction : str
        The direction of the satellite orbit to filter for.
        Input must be one of: One of: 'asc', 'dsc', 'both'.
        If 'both', then no filtering will be applied.
    """
    # Search for granules
    results = utils.get_query_results(bbox, dateRange)
    
    results = filter_direction(results, direction)

    maap = MAAP(maap_host='api.dit.maap-project.org')

    for item in results:
        maap.submitJob(identifier="test-job",
                       algo_id="s1-rtc",
                       version="delay10",
                       username="anonymous",
                       queue="maap-dps-worker-8gb",
                       input_file = item.getDownloadUrl()
                      )
    granule2bursts(input_s1_zip, output_dir)


def download_file_from_https(url_of_file):
    
    maap = MAAP()
    filename = os.path.basename(urlparse(url_of_file).path)
    filepath = os.path.join("input", filename)

    proxy = Result({})
    proxy._cmrFileUrl = maap._SEARCH_GRANULE_URL
    proxy._apiHeader = maap._get_api_header()

    # clear pgt value to simulate a DPS context
    # Update 10/24/2023: This header value is needed to authorize your request, so do not update to ''
    # proxy._apiHeader['proxy-ticket'] = ''

    # Update 10/24/2023: update 10/24/2023: update 10/24/2023: update 10/24/2023: update 10/24/2023: update 10/24/2023: update 10/24/2023: update 10/24/2023: update 10/24/2023: update 10/24/2023: update to just `proxy._dps = maap._DPS`
    # proxy._dps = DpsHelper(proxy._apiHeader, "")
    proxy._dps = maap._DPS

    proxy._getHttpData(url_of_file, False, filepath)

    return filepath


def get_eof_https_link(sentinel_file: str, orbit_type: str = "precise") -> str:
    """Downloads and saves EOF files for specific dates

    Args:
        sentinel_file (str): path to Sentinel-1 filename to download one .EOF for
        orbit_type (str): precise or restituted

    Returns:
        str: https link at ASF to the EOF file corresponding to `sentinel_file`.
    """
    sent = eof.products.Sentinel(sentinel_file)
    orbit_dts, missions = [sent.start_time], [sent.mission]

    # First make sure all are datetimes if given string
    orbit_dts = [dt for dt in orbit_dts]

    asfclient = eof.asf_client.ASFClient()
    urls = asfclient.get_download_urls(orbit_dts, missions, orbit_type=orbit_type)

    return urls[0]


if __name__ == "__main__":
    '''
    Command line script to process a Sentinel-1 RTC.

    Example cmd line call:
    python s1_rtc.py 
        --in_file /projects/zip_dir/S1A_IW_SLC__1SDV_20200827T015033_20200827T015101_034086_03F51E_AE5F.zip
        -o /projects/my-private-bucket/0.0.1-test/
    '''
    # Set up and aws permissions to public bucket (copied from iscesat2_boreal product)
    os.environ['AWS_NO_SIGN_REQUEST'] = 'YES'

    parser = argparse.ArgumentParser()
    
    msg = "The input filename with path for a Sentinel-1 SLC .zip file to be processed through RTC."
    parser.add_argument("-i", "--in_file", type=str, help=msg)
    
    msg = "Path for the output directory to place the S1 bursts that have been processed via RTC"
    parser.add_argument("-o", "--output_dir", type=str, help=msg)
    
    args = parser.parse_args()

    in_file = download_file_from_https(args.in_file)
    
    if not os.path.isfile(in_file):
        raise FileNotFoundError(f"Input File does not exist {args.in_file}")
    
    granule2bursts(in_file, args.output_dir)
    
