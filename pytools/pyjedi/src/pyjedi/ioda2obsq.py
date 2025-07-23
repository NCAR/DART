import numpy as np
from ._utils import _ioda2iodaDF
from ._utils import _iodaDF2obsqDF
from ._utils import _buildObsSeqFromObsqDF

def main():
    inFile = 'sondes_obs_2018041500_m.200.nc4'
    outFile = 'obs_seq.radiosonde.ioda2obsq.out'

    print("INFO: Trying simple radiosonde example")
    print("INFO:    input file: ", inFile)
    print("INFO:    output file: ", outFile)

    # Try simple radiosonde example, put in generic code with configuration later.
    (iodaDF, epochDT) = _ioda2iodaDF(inFile)
    obsqDF = _iodaDF2obsqDF(iodaDF,epochDT)
    obs_seq = _buildObsSeqFromObsqDF(obsqDF)
    obs_seq.write_obs_seq(outFile)

    print("INFO: Converted ", len(obsqDF), " observations")

if __name__ == "__main__":
    main()
