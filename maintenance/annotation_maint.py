import allel
import zarr
import numpy as np
import stdpopsim as stp
import logging
import warnings
import urllib.request
import os

logger = logging.getLogger(__name__)
# make root directory for zarr annotations
annot_path = "annotations"
os.mkdir(annot_path)
# loop through species and download
for spc in stp.all_species():
    if spc.annotations:
        address = spc.annotations[0].url
        genome_version = os.path.basename(address).split(".")[1]
        logger.info(f"Downloading GFF file {spc.id}")
        tmp_path = f"{spc.id}.tmp.gff.gz"
        try:
            x, y = urllib.request.urlretrieve(address, tmp_path)
        except FileNotFoundError:
            warnings.warn("can't connnect to url")
        logger.info(f"creating zarr arrays {spc.id}")
        # create zarr store and zarr root
        spc_path = os.path.join(annot_path, spc.id + "." + genome_version + ".zip")
        store = zarr.ZipStore(spc_path)
        root = zarr.group(store=store, overwrite=True)
        x = allel.gff3_to_dataframe(tmp_path)
        for col_name in x.columns:
            if x[col_name].dtype == "O":
                tmp = root.array(col_name, np.array(x[col_name], dtype=str))
            else:
                tmp = root.array(col_name, np.array(x[col_name]))
        # cleanup
        os.unlink(tmp_path)
        store.close()
