# ------------------------------------
# PYTHON WORKFLOW TO CALL SMART
# ------------------------------------

# NEEDED MODULES
import imas,os,yaml,datetime
import numpy as np
from smart.actor import smart as smart_actor
#from imas_rt_mapping import imas_rt_mapper

# INPUT/OUTPUT CONFIGURATION
file = open('input/scenario.yaml', 'r')
config = yaml.load(file,Loader=yaml.CLoader)
file.close()
shot                = config['shot']
run_in              = config['run_in']
input_user_or_path  = config['input_user_or_path']
input_database      = config['input_database']
if config['output_user_or_path'] == 'default':
    output_user_or_path = os.getenv('USER')
else:
    output_user_or_path = config['output_user_or_path']
output_database     = config['output_database']
run_out             = config['run_out']
time_slice          = config['time_slice']
debug               = config['debug']
use_pellets_ids     = config['use_pellets_ids']

# DISPLAY SIMULATION INFORMATION
print('---------------------------------')
print('shot                = ',shot)
print('run_in              = ',run_in)
print('run_out             = ',run_out)
print('input_user_or_path  = ',input_user_or_path)
print('input_database      = ',input_database)
print('output_user_or_path = ',output_user_or_path)
print('output_database     = ',output_database)
print('time_slice          = ',time_slice)
print('---------------------------------')

# OPEN INPUT DATAFILE TO GET DATA FROM IMAS SCENARIO DATABASE
print('=> Open input datafile')
input = imas.DBEntry(imas.imasdef.MDSPLUS_BACKEND,input_database,shot,run_in,input_user_or_path)
input.open()

# READ INPUT IDSS FROM LOCAL DATABASE
print('=> Read input IDSs')
input_equilibrium = input.get_slice('equilibrium',time_slice,1)
input_core_profiles = input.get_slice('core_profiles',time_slice,1)
input.close()

# PELLETS WRITTEN ON THE FLY (TO BE LATER FILLED VIA WAVEFORM-COOKER OR TAKEN FROM PCSSP)
input_pellets = imas.pellets()
if use_pellets_ids == 1:
    input_pellets.ids_properties.homogeneous_time = 1
    input_pellets.ids_properties.provider = os.getenv('USER')
    input_pellets.ids_properties.creation_date = datetime.datetime.now().strftime("%y-%m-%d")
    input_pellets.time.resize(1)
    input_pellets.time[0] = 0.
    input_pellets.time_slice.resize(1)
    input_pellets.time_slice[0].pellet.resize(1)
    input_pellets.time_slice[0].pellet[0].shape.size = np.array([92.e-9])
    input_pellets.time_slice[0].pellet[0].species.resize(1)
    input_pellets.time_slice[0].pellet[0].species[0].a = 2.5 # (2.5 for 50:50 DT)
    input_pellets.time_slice[0].pellet[0].velocity_initial = 0.3e5

#real_time_data = imas_rt_mapper(input_pellets)
#exit()

# IF LOCAL DATABASE DOES NOT EXIST: CREATE IT
local_database = os.getenv("HOME") + "/public/imasdb/" + output_database + "/3/0"
if os.path.isdir(local_database) == False:
    print("-- Create local database " + local_database)
    os.makedirs(local_database)

# CREATE OUTPUT DATAFILE
print('=> Create output datafile')
output = imas.DBEntry(imas.imasdef.MDSPLUS_BACKEND,output_database,shot,run_out,output_user_or_path)
output.create()

# INITIALIZE THE ACTOR
smart = smart_actor()
code_parameters = smart.get_code_parameters()
code_parameters.parameters_path = 'input/smart.xml'
runtime_settings = smart.get_runtime_settings()
if debug == 1:
    from smart.common.runtime_settings import DebugMode
    runtime_settings.debug_mode = DebugMode.STANDALONE
smart.initialize(code_parameters=code_parameters,runtime_settings=runtime_settings)

# EXECUTE SMART
print('=> Execute SMART')
try:
    output_core_profiles = smart(input_equilibrium, input_core_profiles, input_pellets)
except Exception as error_message:
    print('ERROR in run_smart',str(error_message))
    exit(1)
#input_core_profiles = output_core_profiles

# FINALIZE THE ACTOR
smart.finalize()

# SAVE IDS INTO OUTPUT FILE
print('=> Save IDSs to local database')
output.put(output_core_profiles)
output.put(input_equilibrium)
if use_pellets_ids == 1:
    output.put(input_pellets)
    
output.close()
print('Done exporting.')
