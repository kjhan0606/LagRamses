#!/bin/bash
#SBATCH --job-name=rar-ddo154-hi
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/ddo154_hi_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/ddo154_hi_%j.err
set -Eeuo pipefail
umask 077
module purge
base=/gpfs/kjhan/RAR_QUMOND
out="$base/observations/ddo154_littlethings_${SLURM_JOB_ID:?}"
url=https://bulk.cv.nrao.edu/littlethings/ddo154/HI/DDO154_NA_X0_P_R.FITS
test ! -e "$out"
mkdir -p "$base/observations"
mkdir "$out"
cp -p "${BASH_SOURCE[0]}" "$out/job.sbatch"
curl --fail --location --connect-timeout 20 --max-time 180 --retry 2 -o "$out/hi_moment0.fits.part" "$url"
test -s "$out/hi_moment0.fits.part"
mv "$out/hi_moment0.fits.part" "$out/hi_moment0.fits"
printf '%s\n' "$url" > "$out/source_url.txt"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
env -u LD_LIBRARY_PATH /home/kjhan/miniconda3/bin/python3 - "$out" <<'PY'
import json,sys
from pathlib import Path
import astropy
from astropy.io import fits
import numpy as np
out=Path(sys.argv[1])
with fits.open(out/'hi_moment0.fits',checksum=True) as h:
    h.verify('exception')
    data=np.asarray(h[0].data)
    header=h[0].header
    if data.ndim < 2 or not np.isfinite(data).any():
        raise ValueError('empty or invalid FITS image')
    keys=['OBJECT','BUNIT','BTYPE','BMAJ','BMIN','BPA','CTYPE1','CTYPE2',
          'CUNIT1','CUNIT2','CRVAL1','CRVAL2','CDELT1','CDELT2',
          'RESTFRQ','RESTFREQ','TELESCOP','DATE-OBS','ORIGIN']
    result={'source':'NRAO LITTLE THINGS DDO154 natural-weight image',
            'scope':'FITS integrity and header inventory only; flux calibration, primary-beam correction and radial gas model unverified',
            'astropy_version':astropy.__version__,'shape':list(data.shape),
            'finite_count':int(np.count_nonzero(np.isfinite(data))),
            'header':{k:header.get(k) for k in keys if k in header}}
    (out/'header_inventory.json').write_text(json.dumps(result,indent=2,default=str)+'\n')
    (out/'primary_header.txt').write_text(header.tostring(sep='\n',endcard=True,padding=False)+'\n')
    print(json.dumps(result,indent=2,default=str))
PY
sha256sum "$out/job.sbatch" "$out/source_url.txt" "$out/hi_moment0.fits" "$out/header_inventory.json" "$out/primary_header.txt" > "$out/artifacts.sha256"
sha256sum -c "$out/artifacts.sha256" --status
printf '%s\n' 'OBSERVATION_HEADER_ONLY_NOT_A_VALIDATED_GAS_MODEL' > "$out/COMPLETE"
