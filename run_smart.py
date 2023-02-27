# ------------------------------------
# PYTHON WORKFLOW TO CALL SMART
# ------------------------------------

# NEEDED MODULES
import imas,os,yaml
from smart.actor import smart as smart_actor

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
smart.initialize(code_parameters=code_parameters)

# EXECUTE SMART
print('=> Execute SMART')
try:
    output_core_profiles = smart(input_equilibrium, input_core_profiles)
except Exception as error_message:
    print('ERROR in run_smart',str(error_message))
    exit(1)
#input_core_profiles = output_core_profiles

# FINALIZE THE ACTOR
smart.finalize()

# SAVE IDS INTO OUTPUT FILE
print('=> Append IDS slice to local database')
output.put(output_core_profiles)
    
output.close()
print('Done exporting.')
