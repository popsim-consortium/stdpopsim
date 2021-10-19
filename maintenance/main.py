"""
Command line interface for managing routine maintenance tasks.

This interface is for stdpopsim developers only.
"""
import contextlib
import shutil
import string
import pathlib
import logging

import click
import black
import daiquiri

import stdpopsim
from . import ensembl
from . import ncbi
from . import annotation_maint

logger = logging.getLogger("maint")

species_template = string.Template(
    """
import stdpopsim

from . import genome_data

# [The following are notes for implementers and should be deleted
#  once the recombination rates have been inserted]
# This is the per-chromosome recombination rate, typically the mean
# rate along the chromosome.
# Values in this dictionary are set to -1 by default, so you have
# to update each one. These should be derived from the most reliable
# data and how they were obtained should be documented here.
# The appropriate citation must be added to the list of
# recombination_rate_citations in the Genome object.

_recombination_rate = $chromosome_rate_dict

# [The following are notes for implementers and should be deleted
#  once the mutation rates have been inserted]
# This is the per-chromosome mutation rate, typically the mean
# rate along the chromosome. If per chromosome rates are not available,
# the same value should be used for each chromosome. In this case,
# please use a variable to store this value, rather than repeating
# the same numerical constant, e.g.
# _mutation_rate = {
#    1: _overall_rate,
#    2: _overall_rate,
#    ...
# Values in this dictionary are set to -1 by default, so you have
# to update each one. These should be derived from the most reliable
# data and how they were obtained should be documented here.
# The appropriate citation must be added to the list of
# mutation_rate_citations in the Genome object.

_mutation_rate = $chromosome_rate_dict

_genome = stdpopsim.Genome.from_data(
    genome_data.data,
    recombination_rate=_recombination_rate,
    mutation_rate=_mutation_rate,
    # [ Implementers: please insert citations for the papers you are basing
    # the estimates for recombination and mutation rates. The assembly
    # citation is optional and can be deleted if not needed.]
    citations=[
        stdpopsim.Citation(
            author="", year=-1, doi="", reasons={stdpopsim.CiteReason.ASSEMBLY}),
        stdpopsim.Citation(
            author="", year=-1, doi="", reasons={stdpopsim.CiteReason.REC_RATE}),
        stdpopsim.Citation(
            author="", year=-1, doi="", reasons={stdpopsim.CiteReason.MUT_RATE})
    ]
)

_species = stdpopsim.Species(
    id="$sps_id",
    ensembl_id="$ensembl_id",
    name="$scientific_name",
    common_name="$common_name",
    genome=_genome,
    # [Implementers: you must provide an estimate of the generation_time.
    # Please also add a citation for this.]
    generation_time=0,
    # [Implementers: you must provide an estimate of the population size.
    # TODO: give a definition of what this should be.
    # Please also add a citation for this below..]
    population_size=0,
    citations=[
        stdpopsim.Citation(
            author="", year=-1, doi="", reasons={stdpopsim.CiteReason.POP_SIZE}),
        stdpopsim.Citation(
            author="", year=-1, doi="", reasons={stdpopsim.CiteReason.GEN_TIME})
    ],
)

stdpopsim.register_species(_species)
"""
)

species_test_template = string.Template(
    """
import pytest

import stdpopsim
from tests import test_species


class TestSpeciesData(test_species.SpeciesTestBase):

    species = stdpopsim.get_species("$sps_id")

    def test_ensembl_id(self):
        assert self.species.ensembl_id == "$ensembl_id"

    def test_name(self):
        assert self.species.name == "$scientific_name"

    def test_common_name(self):
        assert self.species.common_name == "$common_name"

    # QC Tests. These tests are performed by another contributor
    # independently referring to the citations provided in the
    # species definition, filling in the appropriate values
    # and deleting the pytest "skip" annotations.
    @pytest.mark.skip("Population size QC not done yet")
    def test_qc_population_size(self):
        assert self.species.population_size == -1

    @pytest.mark.skip("Generation time QC not done yet")
    def test_qc_generation_time(self):
        assert self.species.generation_time == -1

class TestGenomeData(test_species.GenomeTestBase):

    genome = stdpopsim.get_species("$sps_id").genome

    @pytest.mark.skip("Recombination rate QC not done yet")
    @pytest.mark.parametrize(
        ["name", "rate"],
        $chromosome_rate_dict.items())
    def test_recombination_rate(self, name, rate):
        assert rate == pytest.approx(self.genome.get_chromosome(name).recombination_rate)

    @pytest.mark.skip("Mutation rate QC not done yet")
    @pytest.mark.parametrize(
        ["name", "rate"],
        $chromosome_rate_dict.items())
    def test_mutation_rate(self, name, rate):
        assert rate == pytest.approx(self.genome.get_chromosome(name).mutation_rate)
"""
)


def black_format(code):
    return black.format_file_contents(code, fast=False, mode=black.FileMode())


def ensembl_stdpopsim_id(ensembl_id):
    tmp = ensembl_id.split("_")[:2]
    sps_id = "".join([x[0:3].capitalize() for x in tmp])
    if len(sps_id) != 6:
        raise ValueError(f"Cannot extract six character id from {ensembl_id}")
    return sps_id


def ncbi_stdpopsim_id(ncbi_id):
    tmp = ncbi_id.split(" ")[:2]
    sps_id = "".join([x[0:3].capitalize() for x in tmp])
    if len(sps_id) != 6:
        raise ValueError(f"Cannot extract six character id from {ncbi_id}")
    return sps_id


def catalog_path(sps_id):
    return pathlib.Path(f"stdpopsim/catalog/{sps_id}")


def write_catalog_stub(*, path, sps_id, ensembl_id, species_data, genome_data):
    """
    Writes stub files to the catalog for a new species.
    """
    with open(path / "__init__.py", "w") as f:
        print('"""', file=f)
        print(f"Catalog definitions for {sps_id} (Ensembl ID='{ensembl_id}')", file=f)
        print('"""', file=f)
        print("from . import species  # noqa: F401", file=f)

    scientific_name = species_data["scientific_name"]
    common_name = species_data["display_name"]
    logger.info(f"{sps_id}: name={scientific_name}, common_name={common_name}")
    chr_names = genome_data["chromosomes"].keys()
    chromosome_rate_template = {name: -1 for name in chr_names}
    species_code = species_template.substitute(
        ensembl_id=ensembl_id,
        sps_id=sps_id,
        scientific_name=scientific_name,
        common_name=common_name,
        chromosome_rate_dict=chromosome_rate_template,
    )
    path = path / "species.py"
    logger.info(f"Writing species definition stub to {path}")
    with open(path, "w") as f:
        f.write(black_format(species_code))

    test_code = species_test_template.substitute(
        ensembl_id=ensembl_id,
        sps_id=sps_id,
        scientific_name=scientific_name,
        common_name=common_name,
        chromosome_rate_dict=chromosome_rate_template,
    )
    test_path = pathlib.Path("tests") / f"test_{sps_id}.py"
    logger.info(f"Writing species test stub to {test_path}")
    with open(test_path, "w") as f:
        f.write(black_format(test_code))


class DataWriter:
    """
    Writes data obtained from upstream sources into the stdpopsim
    package hierarchy.
    """

    def __init__(self):
        self.ensembl_client = ensembl.EnsemblRestClient()

    @contextlib.contextmanager
    def write(self, path):
        with open(path, "w") as f:
            print("# File autogenerated from Ensembl REST API. Do not edit.", file=f)
            yield f

    def add_species(self, ensembl_id, force=False):
        sps_id = ensembl_stdpopsim_id(ensembl_id)
        logger.info(f"Adding new species {sps_id} for Ensembl ID {ensembl_id}")
        root = catalog_path(sps_id)
        if force:
            shutil.rmtree(root, ignore_errors=True)
        root.mkdir()
        genome_data = self.write_genome_data(ensembl_id)
        species_data = self.ensembl_client.get_species_data(ensembl_id)
        write_catalog_stub(
            path=root,
            sps_id=sps_id,
            ensembl_id=ensembl_id,
            species_data=species_data,
            genome_data=genome_data,
        )

    def add_species_ncbi(self, ncbi_id, force=False):
        species_name = ncbi.get_species_name(ncbi_id)
        tmp = species_name.split(" ")[:2]
        sps_id = "".join([x[0:3].capitalize() for x in tmp])
        if len(sps_id) != 6:
            raise ValueError(f"Cannot extract six character id from {species_name}")
        logger.info(f"Adding new species {sps_id} for NCBI ID {ncbi_id}")
        root = catalog_path(sps_id)
        if force:
            shutil.rmtree(root, ignore_errors=True)
        root.mkdir()
        genome_data = self.write_genome_data_ncbi(ncbi_id, sps_id)
        species_data = ncbi.get_species_data(ncbi_id)
        write_catalog_stub(
            path=root,
            sps_id=sps_id,
            ensembl_id=ncbi_id,
            species_data=species_data,
            genome_data=genome_data,
        )

    def write_genome_data(self, ensembl_id):
        sps_id = ensembl_stdpopsim_id(ensembl_id)
        path = catalog_path(sps_id)
        if not path.exists():
            raise ValueError(
                f"Directory {id} corresponding to {ensembl_id} does" + "not exist"
            )
        logger.info(f"Writing genome data for {sps_id} {ensembl_id}")
        path = path / "genome_data.py"
        data = self.ensembl_client.get_genome_data(ensembl_id)
        code = f"data = {data}"

        # Format the code with Black so we don't get noisy diffs
        with self.write(path) as f:
            f.write(black_format(code))
        return data

    def write_genome_data_ncbi(self, ncbi_id, sps_id):
        path = catalog_path(sps_id)
        if not path.exists():
            raise ValueError(
                f"Directory {id} corresponding to {ncbi_id} does" + "not exist"
            )
        logger.info(f"Writing genome data for {sps_id} {ncbi_id}")
        path = path / "genome_data.py"
        data = ncbi.get_genome_data(ncbi_id)
        code = f"data = {data}"

        # Format the code with Black so we don't get noisy diffs
        with self.write(path) as f:
            f.write(black_format(code))
        return data

    def write_ensembl_release(self):
        release = self.ensembl_client.get_release()
        logger.info(f"Using Ensembl release {release}")
        path = pathlib.Path("stdpopsim/catalog/ensembl_info.py")
        code = f"release = {release}"
        with self.write(path) as f:
            f.write(black_format(code))


################
# User interface
################


@click.group()
@click.option("--debug", is_flag=True)
def cli(debug):
    log_level = logging.INFO
    if debug:
        log_level = logging.DEBUG
    daiquiri.setup(level=log_level)

    # Black's logging is very noisy
    black_logger = logging.getLogger("blib2to3")
    black_logger.setLevel(logging.CRITICAL)


@cli.command()
def list_species():
    """
    List species in stdpopsim with their Ensembl IDs.
    """
    click.echo("ID       Ensembl ID")
    for species in stdpopsim.all_species():
        click.echo(f"{species.id}   {species.ensembl_id}")


@cli.command()
@click.argument("species", nargs=-1)
def update_genome_data(species):
    """
    Update the species genome data from Ensembl. Ensembl IDs for
    species can optionally be provided, e.g.

    update-genome-data homo_sapiens

    will update the genome data for humans. Multiple species can
    be specified. By default all species are updated.
    """
    # TODO make this work for NCBI as well
    if len(species) == 0:
        embl_ids = [s.ensembl_id for s in stdpopsim.all_species()]
    else:
        embl_ids = [s.lower() for s in species]
    writer = DataWriter()
    for eid in embl_ids:
        writer.write_genome_data(eid)
    writer.write_ensembl_release()


@cli.command()
@click.argument("ensembl-id")
@click.option("--force", is_flag=True)
def add_species(ensembl_id, force):
    """
    Add a new species to the catalog using its ensembl ID.
    """
    writer = DataWriter()
    writer.add_species(ensembl_id.lower(), force=force)
    writer.write_ensembl_release()


# TODO refactor this so that it's an option to add-species. By default
# we assume the repository is Ensembl.
@cli.command()
@click.argument("NCBI-id")
@click.option("--force", is_flag=True)
def add_species_ncbi(ncbi_id, force):
    """
    Add a new species to the catalog using its NCBI UID. UIDs can be
    found by searching NCBI for the species in question and looking
    at the assembly page.
    """
    writer = DataWriter()
    writer.add_species_ncbi(ncbi_id, force=force)


@cli.command()
def process_annotations():
    """
    downloads annotations in catalog and preps them for
    upload to aws
    """
    click.echo("downloading annotations")
    annotation_maint.download_process_annotations()


def main():
    cli()
