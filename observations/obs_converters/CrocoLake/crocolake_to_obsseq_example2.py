#!/usr/bin/env python3
import datetime
from convert_crocolake_obs import ObsSequence

crocolake_path = "/path/to/CrocoLake/"

LAT0 = 5
LAT1 = 60
LON0 = -100
LON1 = -30

selected_variables = [
    "DB_NAME",
    "JULD",
    "LATITUDE",
    "LONGITUDE",
    "PRES",
    "TEMP",
    "PRES_QC",
    "TEMP_QC",
    "PRES_ERROR",
    "TEMP_ERROR",
    "PSAL",
    "PSAL_QC",
    "PSAL_ERROR"
]

year0 = 2010
month0 = 5
wmo_id = 1900256
for j in range(10):

    day0 = 1+j
    day1 = day0+1
    date0 = datetime.datetime(year0, month0, day0, 0, 0, 0)
    date1 = datetime.datetime(year0, month0, day1, 0, 0, 0)

    print(f"Converting obs between {date0} and {date1}")

    and_filters = (
        ("LATITUDE",'>',LAT0),  ("LATITUDE",'<',LAT1),
        ("LONGITUDE",'>',LON0), ("LONGITUDE",'<',LON1),
        ("PRES",'>',-1e30), ("PRES",'<',1e30),
        ("JULD",">",date0), ("JULD","<",date1),
        ("PLATFORM_NUMBER","==",str(wmo_id))
    )

    db_filters = [
        list(and_filters) + [("PSAL", ">", -1e30), ("PSAL", "<", 1e30)],
        list(and_filters) + [("TEMP", ">", -1e30), ("TEMP", "<", 1e30)],
    ]

    obs_seq_out = f"obs_seq_CL.ARGO{wmo_id}.{year0}{month0:02d}{day0:02d}.out"
    obsSeq = ObsSequence(
        crocolake_path,
        selected_variables,
        db_filters,
        obs_seq_out=obs_seq_out,
        loose=True
    )
    obsSeq.write_obs_seq()
