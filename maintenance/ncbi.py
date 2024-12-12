"""
Utilites for working with NCBI

"""

import logging
import urllib

from Bio import Entrez

logger = logging.getLogger("ncbi")

Entrez.email = "stdpopsim@popsimconsortium.org"

# TODO This is very badly factored, we're making the same call 3 or 4 times.
# Fix up!


def get_species_name(uid):
    esummary_handle = Entrez.esummary(db="assembly", id=uid, report="full")
    esummary_record = Entrez.read(esummary_handle)
    summary = esummary_record["DocumentSummarySet"]["DocumentSummary"][0]
    return summary["SpeciesName"]


def get_species_data(uid):
    esummary_handle = Entrez.esummary(db="assembly", id=uid, report="full")
    esummary_record = Entrez.read(esummary_handle)
    summary = esummary_record["DocumentSummarySet"]["DocumentSummary"][0]
    species_data = {
        "scientific_name": summary["SpeciesName"],
        "display_name": summary["Organism"],
    }
    return species_data


def get_assembly_summary(id):
    """Get esummary for an entrez id"""
    esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
    esummary_record = Entrez.read(esummary_handle)
    return esummary_record


def get_genome_data(uid):
    """
    Get chromosome data for the specified NCBI UID.
    """
    logger.info(f"Getting genome data for id: {uid}")
    summary = get_assembly_summary(uid)
    sum_dict = summary["DocumentSummarySet"]["DocumentSummary"][0]
    data = {
        "assembly_accession": sum_dict["AssemblyAccession"],
        "assembly_name": sum_dict["AssemblyName"],
    }
    # get ftp link
    url = summary["DocumentSummarySet"]["DocumentSummary"][0]["FtpPath_Stats_rpt"]
    chromosomes = {}
    chrom_list = []
    with urllib.request.urlopen(url) as f:
        stats = f.read().decode("utf-8")
        for line in stats.split("\n"):
            if line and line[0] != "#":
                tokens = line.strip().split("\t")
                if (
                    (tokens[2] == "Chromosome" or tokens[2] == "Mitochondrion")
                    and tokens[3] == "all"
                    and tokens[4] == "total-length"
                ):
                    chromosomes[tokens[1]] = {
                        "length": tokens[5],
                        "synonyms": [],  # not listed in NCBI
                    }
                    chrom_list.append(tokens[1])
    data["chromosomes"] = {name: chromosomes[name] for name in chrom_list}
    return data
