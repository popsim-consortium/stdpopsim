"""
Infrastructure for defining chromosome information for different species.
"""
import os.path
import tempfile
import tarfile

import appdirs
import requests
import msprime

recombination_maps = {
        "HapmapII_GRCh37": "http://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/20110106_recombination_hotspots/HapmapII_GRCh37_RecombinationHotspots.tar.gz"  # NOQA

}

def download_map(map_name, map_dir):
    """
    Downloads the specified map and extracts it into the specified directory.
    """
    os.makedirs(map_dir, exist_ok=True)

    url = recombination_maps[map_name]
    r = requests.get(url, stream=True)
    with tempfile.TemporaryDirectory() as tempdir:
        download_file = os.path.join(tempdir, "download")
        print("Downloading to", download_file)
        with open(download_file, 'wb') as f:
            for chunk in r.iter_content(chunk_size=1024):
                f.write(chunk)
        print("Extracting")
        os.chdir(map_dir)
        with tarfile.open(download_file, 'r') as tf:
            tf.extractall()


class Chromosome(object):

    def __init__(self, name, length, mean_recombination_rate, mean_mutation_rate):
        self.name = name
        self.length = length
        self.mean_recombination_rate = mean_recombination_rate
        self.mean_mutation_rate = mean_mutation_rate


    def _get_recombination_map(self, map_name, chr_name):
        print("get", map_name, chr_name)

        cache_dir = appdirs.user_cache_dir("stdpopsim", "popgensims")
        map_dir = os.path.join(cache_dir, map_name)

        if not os.path.exists(map_dir):
            download_map(map_name, map_dir)
        # TODO need to make the filename pattern here a parameter somehow.
        filename = os.path.join(
            map_dir, "genetic_map_GRCh37_{}.txt".format(self.name))
        return msprime.RecombinationMap.read_hapmap(filename)





