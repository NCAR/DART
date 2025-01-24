import os
import sys
import numpy as np
import onnx
import onnxruntime as ort

from datetime import datetime
from datetime import timedelta

# The directory of your input and output data
model_path = '/glade/work/chennuo/code/pangu/pangu_weather_6.onnx'
model_6 = onnx.load(model_path)

old_date = sys.argv[1]
ninst = sys.argv[2]
input_data_dir = sys.argv[3]
output_data_dir = sys.argv[3]
input_sfc_fname = 'input_surface'+ninst+'.postassim-'+old_date+'.npy'
input_upper_fname = 'input_upper'+ninst+'.postassim-'+old_date+'.npy'

old_date = datetime.strptime(old_date, '%Y-%m-%d-%H')
new_date = old_date + timedelta(hours=6)
new_date = datetime.strftime(new_date, '%Y-%m-%d-%H')

output_sfc_fname = 'output_surface'+ninst+'.forecast-'+new_date+'.npy'
output_upper_fname = 'output_upper'+ninst+'.forecast-'+new_date+'.npy'

# Set the behavier of onnxruntime
options = ort.SessionOptions()
options.enable_cpu_mem_arena=False
options.enable_mem_pattern = False
options.enable_mem_reuse = False
# Increase the number for faster inference and more memory consumption
options.intra_op_num_threads = 1

# Set the behavier of cuda provider
cuda_provider_options = {'arena_extend_strategy':'kSameAsRequested',}

# Initialize onnxruntime session for Pangu-Weather Models
ort_session_6 = ort.InferenceSession(model_path, sess_options=options, providers=['CPUExecutionProvider'])

# Load the upper-air numpy arrays
input = np.load(os.path.join(input_data_dir, input_upper_fname)).astype(np.float32)
# Load the surface numpy arrays
input_surface = np.load(os.path.join(input_data_dir, input_sfc_fname)).astype(np.float32)

# Run the inference session
output, output_surface = ort_session_6.run(None, {'input':input, 'input_surface':input_surface})

# Save the results
np.save(os.path.join(output_data_dir, output_upper_fname), output)
np.save(os.path.join(output_data_dir, output_sfc_fname), output_surface)
