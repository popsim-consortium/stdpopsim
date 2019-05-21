"""
The command line interface for stdpopsim. Allows provides standard simulations
at the command line.
"""
import argparse
import json
import logging
import platform
import sys
import msprime
import tskit

import stdpopsim
from stdpopsim import homo_sapiens


logger = logging.getLogger(__name__)
log_format = '%(asctime)s %(levelname)s %(message)s'


def exit(message):
    """
    Exit with the specified error message, setting error status.
    """
    sys.exit("{}: {}".format(sys.argv[0], message))


def setup_logging(args):
    log_level = "WARN"
    if args.verbosity > 0:
        log_level = "INFO"
    if args.verbosity > 1:
        log_level = "DEBUG"
    logging.basicConfig(level=log_level, format=log_format)


def get_environment():
    """
    Returns a dictionary describing the environment in which stdpopsim
    is currently running.
    """
    env = {
        "os": {
            "system": platform.system(),
            "node": platform.node(),
            "release": platform.release(),
            "version": platform.version(),
            "machine": platform.machine(),
        },
        "python": {
            "implementation": platform.python_implementation(),
            "version": platform.python_version(),
        },
        "libraries": {
            "msprime": {"version": msprime.__version__},
            "tskit": {"version": tskit.__version__},
        }
    }
    return env


def get_provenance_dict():
    """
    Returns a dictionary encoding an execution of stdpopsim conforming to the
    tskit provenance schema.
    """
    document = {
        "schema_version": "1.0.0",
        "software": {
            "name": "stdpopsim",
            "version": stdpopsim.__version__
        },
        "parameters": {
            "command": sys.argv[0],
            "args": sys.argv[1:]
        },
        "environment": get_environment()
    }
    return document


def write_output(ts, args):
    """
    Adds provenance information to the specified tree sequence (ensuring that the
    output is reproducible) and write the resulting tree sequence to output.
    """
    tables = ts.dump_tables()
    logger.debug("Updating provenance")
    provenance = get_provenance_dict()
    tables.provenances.add_row(json.dumps(provenance))
    ts = tables.tree_sequence()
    logger.info(f"Writing to {args.output}")
    ts.dump(args.output)


def write_citations(chromosome, model, args):
    """
    Write out citation information so that the user knows what papers to cite
    for the simulation engine, the model and the mutation/recombination rate
    information.
    """
    if not args.quiet:
        # TODO say this better
        print(
            "If you use this simulation in published work, please cite the following "
            "papers:")
        print("******************")
        print("Simulation engine:")
        print("******************")
        print(
            "\tmsprime: Kelleher et al. 2016: "
            "https://doi.org/10.1371/journal.pcbi.1004842")
        print("******************")
        print("Genetic map:")
        print("******************")
        print("\tTODO")
        # TODO need some way to get a GeneticMap instance from the chromosome. We'll also
        # want to be able to output mutation map, and perhaps other information too, so
        # we want to keep some flexibility for this in mind.
        print("Simulation model:")
        print(f"\t{model.name}: {model.author} {model.year} {model.doi}")


# TODO This setup is quite repetive and we can clearly abstract some of this
# logic about sampling into the model object. What we want is for the model
# to contain the names of the populations plus some documtation descriptions
# from which we automatically generate the CLI and this runner inferface
# can then be entirely generic.

def run_human_gutenkunst_three_pop_ooa(args):
    if args.num_yri_samples + args.num_ceu_samples + args.num_chb_samples < 2:
        exit("Must specify at least 2 samples from YRI, CEU or CHB populations")

    chromosome = homo_sapiens.chromosome_factory(
            args.chromosome, genetic_map=args.genetic_map,
            length_multiplier=args.length_multiplier)
    logger.info(
        f"Running GutenkunstThreePopOutOfAfrica with YRI={args.num_yri_samples}"
        f" CEU={args.num_ceu_samples} CHB={args.num_chb_samples}")
    samples = (
        [msprime.Sample(population=0, time=0)] * args.num_yri_samples +
        [msprime.Sample(population=1, time=0)] * args.num_ceu_samples +
        [msprime.Sample(population=2, time=0)] * args.num_chb_samples)
    model = homo_sapiens.GutenkunstThreePopOutOfAfrica()
    ts = msprime.simulate(
        samples=samples,
        recombination_map=chromosome.recombination_map,
        mutation_rate=chromosome.mutation_rate,
        **model.asdict())
    write_output(ts, args)
    write_citations(chromosome, model, args)


def run_tennessen_two_pop_ooa(args):
    if args.num_european_samples + args.num_african_samples < 2:
        exit("Must specify at least 2 samples")

    chromosome = homo_sapiens.chromosome_factory(
            args.chromosome, genetic_map=args.genetic_map,
            length_multiplier=args.length_multiplier)
    model = homo_sapiens.TennessenTwoPopOutOfAfrica()
    logger.info(
        f"Running TennessenTwoPopOutOfAfrica with European={args.num_european_samples}"
        f" African={args.num_african_samples}")
    samples = (
        [msprime.Sample(population=0, time=0)] * args.num_european_samples +
        [msprime.Sample(population=1, time=0)] * args.num_african_samples)
    ts = msprime.simulate(
        samples=samples,
        recombination_map=chromosome.recombination_map,
        mutation_rate=chromosome.mutation_rate,
        **model.asdict())
    write_output(ts, args)
    write_citations(chromosome, model, args)


def run_browning_america(args):
    total_samples = (
        args.num_european_samples + args.num_african_samples
        + args.num_asian_samples + args.num_american_samples)
    if total_samples < 2:
        exit("Must specify at least 2 samples")

    chromosome = homo_sapiens.chromosome_factory(
            args.chromosome, genetic_map=args.genetic_map,
            length_multiplier=args.length_multiplier)
    logger.info(
        f"Running BrowningAmerica with European={args.num_european_samples} "
        f"African={args.num_african_samples}, Asian={args.num_asian_samples} "
        f"African={args.num_american_samples}")

    samples = (
        [msprime.Sample(population=0, time=0)] * args.num_european_samples +
        [msprime.Sample(population=1, time=0)] * args.num_african_samples +
        [msprime.Sample(population=2, time=0)] * args.num_asian_samples +
        [msprime.Sample(population=3, time=0)] * args.num_american_samples)

    model = homo_sapiens.BrowningAmerica()
    ts = msprime.simulate(
        samples=samples,
        recombination_map=chromosome.recombination_map,
        mutation_rate=chromosome.mutation_rate,
        **model.asdict())
    write_output(ts, args)
    write_citations(chromosome, model, args)


def add_output_argument(parser):
    parser.add_argument(
        "output",
        help="The file to write simulated tree sequence to")


def stdpopsim_cli_parser():

    # TODO the CLI defined by this hierarchical and clumsy, but it's the best
    # I could figure out. It can definitely be improved!
    top_parser = argparse.ArgumentParser(
        description="Run simulations defined by stdpopsim from the command line")
    top_parser.add_argument(
        "-V", "--version", action='version',
        version='%(prog)s {}'.format(stdpopsim.__version__))
    top_parser.add_argument(
        "-v", "--verbosity", action='count', default=0,
        help="Increase the verbosity")
    top_parser.add_argument(
        "-q", "--quiet", action='store_true',
        help="Do not write out citation information")
    subparsers = top_parser.add_subparsers(dest="subcommand")
    subparsers.required = True

    species_parser = subparsers.add_parser(
        "homo-sapiens",
        help="Run simulations of human history.")
    species_parser.add_argument(
        "-g", "--genetic-map", default=None,
        # TODO use the genetic map registry
        choices=["HapmapII_GRCh37", "Decode_2010_sex_averaged"],
        help="Specify a particular genetic map. Use a flat map by default.")
    species_parser.add_argument(
        "-c", "--chromosome", default="chr22",
        help="Simulate a specific chromosome")
    species_parser.add_argument(
        "-l", "--length-multiplier", default=1, type=float,
        help="Simulate a chromsome of length l times the named chromosome")
    subsubparsers = species_parser.add_subparsers(dest="subcommand")
    subsubparsers.required = True

    parser = subsubparsers.add_parser(
        "GutenkunstThreePopOutOfAfrica",
        help="Runs the three population out-of-Africa model")
    parser.set_defaults(runner=run_human_gutenkunst_three_pop_ooa)
    add_output_argument(parser)
    parser.add_argument(
        "--num-yri-samples", default=0, type=int,
        help="The number of samples to take from the YRI population")
    parser.add_argument(
        "--num-ceu-samples", default=0, type=int,
        help="The number of samples to take from the CEU population")
    parser.add_argument(
        "--num-chb-samples", default=0, type=int,
        help="The number of samples to take from the CHB population")

    parser = subsubparsers.add_parser(
        "TennessenTwoPopOutOfAfrica",
        help="Runs the TennessenTwoPopOutOfAfrica model.")
    add_output_argument(parser)
    parser.add_argument(
        "--num-african-samples", default=0, type=int,
        help="The number of samples to take from the African population")
    parser.add_argument(
        "--num-european-samples", default=0, type=int,
        help="The number of samples to take from the European population")
    parser.set_defaults(runner=run_tennessen_two_pop_ooa)

    parser = subsubparsers.add_parser(
        "BrowningAmerica",
        help="Runs the BrowningAmerica model.")
    add_output_argument(parser)
    parser.add_argument(
        "--num-african-samples", default=0, type=int,
        help="The number of samples to take from the African population")
    parser.add_argument(
        "--num-european-samples", default=0, type=int,
        help="The number of samples to take from the European population")
    parser.add_argument(
        "--num-asian-samples", default=0, type=int,
        help="The number of samples to take from the Asian population")
    parser.add_argument(
        "--num-american-samples", default=0, type=int,
        help="The number of samples to take from Admixed American population")
    parser.set_defaults(runner=run_browning_america)

    # Add stubs for discussion
    species_parser = subparsers.add_parser(
        "arabadopsis-thaliana",
        help="Run simulations of Arabadopsis thaliana.")

    species_parser = subparsers.add_parser(
        "e-coli",
        help="Run simulations of E-coli.")

    return top_parser


def stdpopsim_main(arg_list=None):
    parser = stdpopsim_cli_parser()
    args = parser.parse_args(arg_list)
    setup_logging(args)
    if not args.quiet:
        print("*****************************")
        print("**        WARNING          **")
        print("*****************************")
        print("This inferface is highly experimental and *will* change. ")
        print("It is provided as a work-in-progress to get some feedback on ")
        print("how to better structure the interfaces")
        print("Disable this message with the -q option")
        print("*****************************")
        print("**      END WARNING        **")
        print("*****************************")
    args.runner(args)
