#!/usr/bin/env python3
##########################################################################
import os
from datetime import datetime
import dask.dataframe as dd
import gsw
import numpy as np
import pandas as pd
##########################################################################


class ObsSequence:
    """class ObsSequence: methods to load observations from CrocoLake and store
    them in obs_seq format

    """

    # ------------------------------------------------------------------ #
    # Constructors/Destructors                                           #
    # ------------------------------------------------------------------ #

    def __init__(self, crocolake_path, selected_vars, db_filters,
                 fill_na_qc=None, fill_na_error=None, loose=False,
                 obs_seq_out="obs_seq.out"):
        """Constructor

        Arguments:
        crocolake_path (str):  path to desired crocolake database
        selected_vars (list):  list of variables to be extracted from the database
        db_filters (list):     list of db_filters to be applied to the database
        fill_na_qc (int):      replace value for NA in QC flags
        fill_na_error (float): replace value for NA in error variables
        obs_seq_out (str):     obs_seq file name
        loose (bool):          if True, store observation values also when
                               their QC and error are not present (default: False)

        """
        self.crocolake_path = crocolake_path
        self.selected_vars = selected_vars
        self.db_filters = db_filters
        self.obs_seq_out = obs_seq_out
        self.loose = loose

        self.obs_condition = lambda row, var: (
            pd.notna(row[var]) and pd.notna(row[var + "_ERROR"]) and pd.notna(row[var + "_QC"])
        )
        if self.loose:
            self.obs_condition = lambda row, var: pd.notna(row[var])

        self._load_data(fill_na_qc,fill_na_error)


    # ------------------------------------------------------------------ #
    # Methods                                                            #
    # ------------------------------------------------------------------ #

#------------------------------------------------------------------------------#
    def _load_data(self,fill_na_qc,fill_na_error):
        """Load data into dataframe

        Arguments:
         fill_na_qc (int):      replace value for NA in QC flags
         fill_na_error (float): replace value for NA in error variables
        """

        ddf = dd.read_parquet(
            self.crocolake_path,
            engine="pyarrow",
            columns=self.selected_vars,
            filters=self.db_filters
        )

        ddf = ddf.dropna(subset=["JULD"])

        if fill_na_qc is not None:
            fill_cols = [
                varname + "_QC"
                for varname in self.selected_vars
                if varname+"_QC" in ddf.columns
            ]
            ddf[fill_cols] = ddf[fill_cols].fillna(fill_na_qc)
        if fill_na_error is not None:
            fill_cols = [
                varname + "_ERROR"
                for varname in self.selected_vars
                if varname+"_ERROR" in ddf.columns
            ]
            ddf[fill_cols] = ddf[fill_cols].fillna(fill_na_error)

        self.ddf = ddf
        return

#------------------------------------------------------------------------------#
    def _get_obs_type_vars(self):
        """Build dictionary of available variable names for obs_type"""

        obs_type_vars = {}
        obs_type_vars["TEMP"] = "TEMPERATURE"
        obs_type_vars["PSAL"] = "SALINITY"
        obs_type_vars["DOXY"] = "OXYGEN"
        obs_type_vars["TOT_ALKALINITY"] = "ALKALINITY"
        obs_type_vars["TCO2"] = "INORGANIC_CARBON"
        obs_type_vars["NITRATE"] = "NITRATE"
        obs_type_vars["SILICATE"] = "SILICATE"
        obs_type_vars["PHOSPHATE"] = "PHOSPHATE"
        return obs_type_vars

#------------------------------------------------------------------------------#
    def _get_obs_type_source(self,varname):
        """Build dictionary of available source names for obs_type"""
        if varname=="ARGO":
            return "ARGO"
        elif varname=="GLODAP":
            return "BOTTLE"
        elif varname=="SprayGliders":
            return "GLIDER"
        else:
            return ""

#------------------------------------------------------------------------------#
    def _timestamp_to_days_seconds(self,timestamp_str):
        """Convert YYYYMMDD HH:MM:SS timestamp to number of days, number of
        seconds since 1601-01-01

        Arguments:
        timestamp_str (str): timestamp in string format

        Returns:
        days (int): number of days since 1601-01-01
        seconds (int): number of seconds since (1601-01-01 + days)
        """

        timestamp = datetime.strptime(timestamp_str, "%Y-%m-%d %H:%M:%S")
        reference_date = datetime(1601, 1, 1)
        time_difference = timestamp - reference_date
        days = time_difference.days
        seconds_in_day = time_difference.seconds

        return days, seconds_in_day

#------------------------------------------------------------------------------#
    def write_obs_seq(self):
        """Write observations to obs_seq file"""

        ddf = self.ddf

        if not len(ddf)>0:
            print("No observations found. Interrupting execution.")
            return

        # Rename time column
        if "JULD" in ddf.columns:
            ddf = ddf.rename(columns={"JULD": "TIME"})

        ddf = ddf.rename(columns={"LONGITUDE": "LONGITUDE_DEG", "LATITUDE": "LATITUDE_DEG"})

        # Generate vertical coordinate in meters (positive downwards, gsw returns positive upwards)
        ddf['VERTICAL'] = -gsw.conversions.z_from_p(ddf['PRES'], ddf['LATITUDE_DEG'], 0)

        # Convert lat and lon to radians and lon to 0:360 range (CrocoLake is in -180:180)
        # not elegant but pyarrow backend does not support modulo operator
        ddf['LONGITUDE_DEG'] = ddf['LONGITUDE_DEG'].astype("float64")
        ddf['LONGITUDE_DEG'] = ddf['LONGITUDE_DEG'] % 360
        ddf['LONGITUDE_DEG'] = ddf['LONGITUDE_DEG'].astype("float64[pyarrow]")

        ddf['LONGITUDE'] = np.deg2rad(ddf['LONGITUDE_DEG'])
        ddf['LATITUDE'] = np.deg2rad(ddf['LATITUDE_DEG'])

        # Convert salinity from PSUs (g/kg) to kg/kg
        if 'PSAL' in ddf.columns:
            ddf['PSAL'] = ddf['PSAL']*1e-3

        # Sort dataframe by TIME columns
        ddf = ddf.repartition(npartitions=1)
        ddf = ddf.sort_values(by='TIME')
        ddf = ddf.set_index("TIME")
        ddf = ddf.reset_index(drop=False).rename(columns={'index': 'obs_num'})

        # Repartition for memory-friendly sizes
        ddf = ddf.repartition(partition_size="300MB")

        obs_tmp_file = self.obs_seq_out + ".tmp"

        # Find unique values in the "DB_NAME" column
        obs_type = ddf['DB_NAME'].unique().compute()

        print("Found the following data sources: " + str(obs_type))

        # Build dictionaries for obs_type sources and measurements
        obs_type_sources = {varname: self._get_obs_type_source(varname) for varname in obs_type}
        obs_type_vars = self._get_obs_type_vars()
        obs_type_vars = {
            varname: obs_type_vars[varname]
            for varname in obs_type_vars
            if varname in ddf.columns.to_list()
        }

        # Start building header
        header = [
            ' obs_sequence',
            'obs_type_definitions'
        ]

        # Append observations kinds
        obs_type_nb = (len(obs_type_sources)+1)*len(obs_type_vars)
        header.append(f'          {obs_type_nb}')
        idx=10
        obs_dict = {}
        for varname in obs_type_vars.values():
            idx+=1
            kind_name = varname
            obs_dict[kind_name] = idx
            header.append(f'          {idx:02d} {kind_name}')
            for source in obs_type_sources.values():
                if not source=="":
                    idx+=1
                    kind_name = source+"_"+varname
                    obs_dict[kind_name] = idx
                    header.append(f'          {idx:02d} {kind_name}')

        cumul_obs_count = 0 # initializing value outside with statement
        with open(obs_tmp_file, 'w', encoding='ascii') as f:

            #------------------------------------------------------------------------------#
            def generate_linked_list_pattern(obs_ind,indmax):
                """Create a list of strings with the linked list pattern for indmax lines.

                Arguments:
                obs_ind : current observation number
                indmax (int) :    total number of observations to convert

                Returns:
                result : nested list for observation obs_ind
                """
                ind0 = obs_ind
                result = []
                col1 = ind0 if ind0 > 0 else -1
                col2 = ind0 + 2
                col3 = -1
                result = f"{col1:<12}{col2:<11}{col3}"
                if (ind0+1)==indmax:
                    result = f"{indmax-1:<12}{'-1':<11}{'-1'}"

                return result

            #------------------------------------------------------------------------------#
            def row_to_obs(df,cols_obs,indmax):
                """Convert dataframe row into obs_seq format

                Arguments:
                df (row) :        current row of observations dataframe
                cols_obs (list) : observations names to convert
                indmax (int) :    total number of observations to convert

                Returns:
                obs (list) : rows to print to file
                """
                # df technically is just a row
                obs = []
                # for var in params["CROCOLAKE_BGC"]:
                ref_count_obs = df["obs_count"] # nb of obs to include for this row
                obs_num_row = 0
                for var in cols_obs:
                    #skip depending on loose condition or not
                    if self.obs_condition(df,var):
                        obs_num_row += 1 #current number of the obs within this row
                        if obs_num_row > ref_count_obs:
                            raise ValueError(
                                "Retrieving more observations than expected. "
                                "There might be a bug."
                            )

                        # current number of the obs within the whole db
                        obs_num = df["cumul_obs_count"] + obs_num_row
                        obs.append('OBS        ' + str(obs_num))
                        obs.append(df[var])
                        if pd.notna(df[var+"_QC"]):
                            obs.append(df[var+"_QC"])  # qc of current obs value
                        else:
                            obs.append(0.)  # qc of current obs value

                        obs_ind = obs_num - 1 #index starts from 0
                        linked_list = generate_linked_list_pattern(obs_ind,indmax)
                        obs.append(linked_list)  # linked list info

                        obs.append('obdef')
                        obs.append('loc3d')
                        #location x, y, z, vert-type (meters = 3)
                        obs.append(
                            '   '.join(map(str, df[ ['LONGITUDE','LATITUDE','VERTICAL'] ])) + '   ' + str(3)
                        )

                        obs.append('kind') # this is type of observation
                        obs_type_name = obs_type_vars[var]
                        obs_src = obs_type_sources[ df["DB_NAME"] ]
                        if not obs_src == "":
                            obs_type_name = obs_src  + "_" + obs_type_name
                        obs_type = obs_dict[obs_type_name]
                        obs.append(obs_type)  # observation type

                        days, seconds = self._timestamp_to_days_seconds(
                            df["TIME"].strftime('%Y-%m-%d %H:%M:%S')
                        )
                        obs.append(' '.join(map(str, [seconds, days])))  # seconds, days

                        if pd.notna(df[var+"_ERROR"]):
                            obs.append(df[var+"_ERROR"])  # obs error variance
                        else:
                            obs.append(-888888.)

                return obs

            #------------------------------------------------------------------------------#
            def write_row(row,cols_obs,indmax):
                """Convert observation and write it to file

                Arguments:
                row (row) :       current row of observations dataframe
                cols_obs (list) : observations names to convert
                indmax (int) :    total number of observations to convert
                """

                ob_write = row_to_obs(row,cols_obs,indmax)
                for line in ob_write:
                    f.write(str(line) + '\n')

            #------------------------------------------------------------------------------#
            def build_cols_obs(column_names):
                """Build list of names of observations to convert

                Arguments:
                column_names (list) : dataframe's column names
                """
                j = 0
                while j < len(column_names):
                    var = column_names[j]
                    if (
                        var in ["PLATFORM_NUMBER","TIME","DATA_MODE","DB_NAME","PRES","VERTICAL"]
                        or ("QC" in var)
                        or ("ERROR" in var)
                        or ("LATITUDE" in var)
                        or ("LONGITUDE" in var)
                    ):
                        column_names.pop(j)
                    else:
                        j += 1
                return column_names

            #------------------------------------------------------------------------------#
            def count_obs(row, column_names):
                """Count total number of observations to store to obs_seq file"""
                obs_count = 0
                for var in column_names:
                    if self.obs_condition(row,var):
                        obs_count += 1
                return obs_count

            ind0 = 0
            indmax = -1
            cols_obs = build_cols_obs(ddf.columns.to_list())

            print("Writing the following observations: " + str(cols_obs))

            for p in range(ddf.npartitions):
                # converting dask dataframe to pandas df for writing operations
                df = ddf.partitions[p].compute()

                # count obs for each row; count cumulative obs in previous rows
                df["obs_count"] = df.apply(
                    count_obs,
                    axis=1,
                    args=(cols_obs,)
                )
                # shift and fill_value because cumsum() includes current row, but I
                # want to count up to previous row
                df["cumul_obs_count"]=df["obs_count"].cumsum().shift(1,fill_value=0) + cumul_obs_count
                if p == (ddf.npartitions-1):
                    indmax = df["cumul_obs_count"].iloc[-1] + df["obs_count"].iloc[-1]

                if p%10==0 or p == (ddf.npartitions):
                    print()
                    print("Writing partition "+str(p+1)+" of " + str(ddf.npartitions) + ".")

                # add obs per each row
                df.apply(
                    write_row,
                    axis=1,
                    args=(cols_obs,indmax)
                )
                ind0 += len(df)
                cumul_obs_count = df["cumul_obs_count"].iloc[-1] + df["obs_count"].iloc[-1]
                print("cumul_obs_count: " + str(cumul_obs_count) )
                print("df cumul_obs_count: " + str(df["cumul_obs_count"].iloc[-1]) )
                print("df obs_count: " + str(df["obs_count"].iloc[-1]) )
            print("Done.")

        # Finish building header (cannot do it earlier because we need to know the number of observations)
        num_obs = cumul_obs_count
        print("Total number of observations: " + str(num_obs))

        # Add number of observations and qc
        num_copies = 1
        num_qcs = 1
        header.append(f"num_copies: {num_copies:>7} num_qc: {num_qcs:>15}")
        header.append(f'num_obs: {num_obs:>10} max_num_obs: {num_obs:>10}')
        header.append('CROCOLAKE observation')
        header.append('CROCOLAKE QC')
        first = 1
        header.append(f'first: {first:>12} last: {num_obs:>12}')

        # Write obs_seq.out file
        with open(self.obs_seq_out, 'w', encoding='ascii') as f:

            # Write header
            for line in header:
                f.write(str(line) + '\n')

            # Write observations
            with open(obs_tmp_file, 'r', encoding='ascii') as ftmp:
                for line in ftmp:
                    f.write(line)

        # Remove temporary file
        os.remove(obs_tmp_file)

        print("Observations written to " + self.obs_seq_out)

        return
