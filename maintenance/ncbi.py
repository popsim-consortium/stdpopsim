"""
Utilites for working with NCBI

"""
import logging
from Bio import Entrez
import urllib
from datetime import datetime

logger = logging.getLogger("ncbi")


def get_species_data(id):
    species_data = {
        "scientific_name": id,
        "display_name": id,
    }
    return species_data


def get_assembly_summary(id):
    """Get esummary for an entrez id"""
    esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
    esummary_record = Entrez.read(esummary_handle)
    return esummary_record


def get_genome_data(species_name):
    """
    Searches NCBI for assemblies associated with species name
    finds most recently submitted assembly and sucks in that
    chromosome length information
    """

    from Bio import Entrez

    # provide your own mail here
    Entrez.email = "stdpopsim@popsimconsortium.org"
    handle = Entrez.esearch(db="assembly", term=species_name, retmax="200")
    record = Entrez.read(handle)
    ids = record["IdList"]
    logger.info(f"Searching NCBI for {species_name}")
    logger.info(f"found {len(ids)} ids")

    # find most recent assembly
    recent = datetime.strptime("1900/01/01 11:11", "%Y/%m/%d %H:%M")
    recent_id = "666"
    for id in ids:
        # get summary
        summary = get_assembly_summary(id)
        dt = summary["DocumentSummarySet"]["DocumentSummary"][0]["SubmissionDate"]
        if datetime.strptime(dt, "%Y/%m/%d %H:%M") > recent:
            recent = datetime.strptime(dt, "%Y/%m/%d %H:%M")
            recent_id = id
    logger.info(f"most recent id: {recent_id}")
    # pull the summary we want
    summary = get_assembly_summary(recent_id)
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
                    tokens[2] == "Chromosome"
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
