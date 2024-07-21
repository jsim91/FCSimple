import flowkit as fk
import numpy as np
import pandas as pd
import sys
import os

in_file = sys.argv[1]
out_file = sys.argv[2]
in_t = float(sys.argv[3])
in_w = float(sys.argv[4])
in_m = float(sys.argv[5])
in_a = float(sys.argv[6])

data = pd.read_csv(in_file)

try:
    os.remove(in_file)
except OSError:
    pass

hyperlog_xform = fk.transforms.HyperlogTransform(
    'hyperlog',
    param_t=in_t,
    param_w=in_w,
    param_m=in_m,
    param_a=in_a
)

sample_from_df = fk.Sample(fcs_path_or_data = data, sample_id='my_sample_from_dataframe', channel_labels=['x'], subsample=0)
sample_from_df.apply_transform(hyperlog_xform)
df_events = sample_from_df.as_dataframe(source='xform')
df_events.to_csv(out_file + '/__tmp_exprs__.csv', index = False)
