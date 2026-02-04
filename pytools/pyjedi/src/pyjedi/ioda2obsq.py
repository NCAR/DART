import argparse
import numpy as np
import pandas as pd
from ._utils import _parseYamlConfig
from ._utils import _ioda2iodaDF
from ._utils import _iodaDF2obsqDF
from ._utils import _buildObsSeqFromObsqDF

def main():
    # Parse arguments
    argParser = argparse.ArgumentParser(description="Convert IODA obs file to an obs_seq file.")
    
    # Required arguments
    argParser.add_argument("yamlConfigFile", help="YAML configuration file")
    argParser.add_argument("iodaInFile", help="IODA input file")
    argParser.add_argument("obsqOutFile", help="obs_seq output file")

    # Optional arguements
    argParser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose output')

    # Parse the arguments
    args = argParser.parse_args()
    configFile = args.yamlConfigFile
    inFile = args.iodaInFile
    outFile = args.obsqOutFile
    verbose = args.verbose

    print("INFO: Converting IODA file to obsq file:")
    print("INFO:    YAML configuration file: ", configFile)
    print("INFO:    IODA input file: ", inFile)
    print("INFO:    obs_seq output file: ", outFile)
    print("")

    # Parse the YAML config
    (iodaVarsConfig, obsCategoryConfig) = _parseYamlConfig(configFile)
    if (verbose):
        print("DEBUG: IODA variable configuration: ", iodaVarsConfig)
        print("DEBUG: Observation category configuration: ", obsCategoryConfig)
        print("")

    # Conversion steps
    #   1. read ioda file into a pandas dataframe (which will have a different layout
    #      compared to what the obsq writer requires)
    #   2. convert the ioda dataframe layout to the obsq dataframe layout
    #      Do this one variable at a time. When there are multiple variables, iterate
    #      through the variables, convert the current variable to an obsq dataframe
    #      and then concatenate that dataframe to the main dataframe
    #   3. build a pyDARTdiags ObsSequence object from the main dataframe
    #   4. call the ObsSequence.write_obs_seq method to create the output file
    (iodaDF, epochDT) = _ioda2iodaDF(inFile, obsCategoryConfig)
    if (verbose):
        print("DEBUG: IODA dataframe columns: ", iodaDF.columns.tolist())
        print("DEBUG: IODA dataframe shape: ", iodaDF.shape)
        print("")

    obsqDF = pd.DataFrame()
    for iodaVarConfig in iodaVarsConfig:
        iodaVarName = iodaVarConfig['name']
        iodaVarType = iodaVarConfig['type']
        print("INFO: IODA input variable: ", iodaVarName, " (", iodaVarType, ")")
        obsqDF = pd.concat([obsqDF, _iodaDF2obsqDF(iodaDF, epochDT, iodaVarName,
                                                   iodaVarType, obsCategoryConfig)],
                           ignore_index = True)
    print("")

    obs_seq = _buildObsSeqFromObsqDF(obsqDF)
    obs_seq.write_obs_seq(outFile)

    print("INFO: Converted", len(obsqDF), "observations")

if __name__ == "__main__":
    main()
