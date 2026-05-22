# ------------------------------------
# PYTHON WORKFLOW TO CALL SMART
# ------------------------------------

import os, copy, datetime, yaml
import numpy as np
import imas
from smart.actor import smart as smart_actor
from time_compute import time_compute

# INPUT/OUTPUT CONFIGURATION
with open('input/scenario.yaml', 'r') as f:
    config = yaml.load(f, Loader=yaml.CLoader)

pulse               = config['pulse']
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
time_no_more_pellet = config['time_no_more_pellet']
dt_required         = config['dt_required']
ntimes              = config['ntimes']
debug               = config['debug']
use_pellets_ids     = config['use_pellets_ids']
scenario_backend    = config.get('scenario_backend', 'hdf5')

# DISPLAY SIMULATION INFORMATION
print('---------------------------------')
print('pulse               = ', pulse)
print('run_in              = ', run_in)
print('run_out             = ', run_out)
print('input_user_or_path  = ', input_user_or_path)
print('input_database      = ', input_database)
print('output_user_or_path = ', output_user_or_path)
print('output_database     = ', output_database)
print('time_slice          = ', time_slice)
print('ntimes              = ', ntimes)
print('---------------------------------')

# Build URIs for input (read) and output (write).
# 'version=3' refers to the imasdb on-disk layout, NOT the DD version.
input_uri = (
    f"imas:{scenario_backend}?user={input_user_or_path}"
    f";database={input_database};version=3;pulse={pulse};run={run_in}"
)
output_uri = (
    f"imas:hdf5?user={output_user_or_path}"
    f";database={output_database};version=3;pulse={pulse};run={run_out}"
)

# If the legacy imasdb directory does not exist yet, create it.
local_database = os.path.join(os.getenv('HOME'), 'public', 'imasdb',
                              output_database, '3', '0')
if not os.path.isdir(local_database):
    print('-- Create local database ' + local_database)
    os.makedirs(local_database)

print('=> Open input datafile')
print('=> Create output datafile')
with imas.DBEntry(input_uri, 'r') as input_entry, \
     imas.DBEntry(output_uri, 'w') as output_entry:

    # PELLETS WRITTEN ON THE FLY (TO BE LATER FILLED VIA WAVEFORM-COOKER OR
    # TAKEN FROM PCSSP)
    input_pellets = imas.IDSFactory().pellets()
    if use_pellets_ids == 1:
        input_pellets.ids_properties.homogeneous_time = (
            imas.ids_defs.IDS_TIME_MODE_HOMOGENEOUS
        )
        input_pellets.ids_properties.provider = os.getenv('USER')
        input_pellets.ids_properties.creation_date = (
            datetime.datetime.now().strftime('%y-%m-%d')
        )
        input_pellets.time.resize(1)
        input_pellets.time[0] = 0.
        input_pellets.time_slice.resize(1)
        input_pellets.time_slice[0].pellet.resize(1)
        input_pellets.time_slice[0].pellet[0].shape.type.index = 2
        input_pellets.time_slice[0].pellet[0].shape.size.resize(2)
        input_pellets.time_slice[0].pellet[0].shape.size[0] = 5.0 / 2.0 * 1.0e-3
        input_pellets.time_slice[0].pellet[0].shape.size[1] = (
            33.0 / (np.pi * 2.5**2) * 1.0e-3
        )
        input_pellets.time_slice[0].pellet[0].species.resize(1)
        input_pellets.time_slice[0].pellet[0].species[0].a = 2.5   # 2.5 for 50:50 DT
        input_pellets.time_slice[0].pellet[0].velocity_initial = 0.3e3

    # READ FULL TIME VECTOR OF EQUILIBRIUM IDS TO GET THE TIME BASE
    time_array, it = time_compute(input_entry, time_slice, ntimes, dt_required)

    # INITIALIZE THE ACTOR
    smart = smart_actor()
    code_parameters = smart.get_code_parameters()
    code_parameters.parameters_path = 'input/smart.xml'
    runtime_settings = smart.get_runtime_settings()
    if debug == 1:
        from smart.common.runtime_settings import DebugMode
        runtime_settings.debug_mode = DebugMode.STANDALONE
    smart.initialize(code_parameters=code_parameters,
                     runtime_settings=runtime_settings)

    # TIME LOOP
    FirstTime = True
    for itime in range(it, it + ntimes):

        if len(time_array) > 1:
            print('Time = %5.2f' % time_array[itime],
                  's, itime = ', itime, '/', it + ntimes - 1)
            time = time_array[itime]
        else:
            time = time_array[0]

        # READ IDSes FROM INPUT SCENARIO
        print('=> Read input IDSes')
        input_equilibrium = input_entry.get_slice(
            'equilibrium', time, imas.ids_defs.CLOSEST_INTERP)
        if FirstTime:
            input_core_profiles = input_entry.get_slice(
                'core_profiles', time, imas.ids_defs.CLOSEST_INTERP)
            FirstTime = False

        # Stop pellet injection after a while
        print('time >= time_no_more_pellet', time, time_no_more_pellet)
        if time >= time_no_more_pellet:
            input_pellets.time_slice[0].pellet[0].shape.size = np.array([1.e-37])

        # EXECUTE SMART
        print('=> Execute SMART')
        try:
            output_core_profiles = smart(input_equilibrium,
                                         input_core_profiles, input_pellets)
            input_core_profiles = copy.deepcopy(output_core_profiles)
        except Exception as error_message:
            print('ERROR in run_smart', str(error_message))
            exit(1)

        output_entry.put_slice(input_equilibrium)
        output_entry.put_slice(input_pellets)
        output_entry.put_slice(output_core_profiles)
        print('Output time = %5.2f s' % (output_core_profiles.time[0]))

    # FINALIZE THE ACTOR
    smart.finalize()

    # SAVE IDS INTO OUTPUT FILE
    print('=> Save IDSs to local database')
    if use_pellets_ids == 1:
        output_entry.put(input_pellets)

print('Done exporting.')
