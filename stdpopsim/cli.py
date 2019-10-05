"""
The command line interface for stdpopsim. Provides standard simulations
at the command line and methods to manage resources used by stdpopsim.
"""
import argparse
import json
import logging
import platform
import sys
import textwrap

import msprime
import tskit
import humanize
import daiquiri

import stdpopsim

# resource is from the standard library, but it's not available on
# Windows. We break from the usual import grouping conventions here
# to avoid lots of pep8 complaints about mixing imports and code.
_resource_module_available = False
try:
    import resource
    _resource_module_available = True
except ImportError:
    pass


logger = logging.getLogger(__name__)


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
    daiquiri.setup(level=log_level)


def get_models_help(species_id, model_id):
    """
    Generate help text for the specified species. If model_id is None, generate
    help for all models. Otherwise, it must be a string with a valid model ID.
    """
    species = stdpopsim.get_species(species_id)
    if model_id is None:
        models_text = f"\nAll simulation models for {species.name}\n\n"
        models = [model.id for model in species.models]
    else:
        models = [model_id]
        models_text = f"\nModel description\n\n"

    # TODO improve this text formatting.
    indent = " " * 4
    wrapper = textwrap.TextWrapper(initial_indent=indent, subsequent_indent=indent)
    for model_id in models:
        model = species.get_model(model_id)
        models_text += f"{model.id}: {model.name}\n"
        models_text += wrapper.fill(textwrap.dedent(model.description))
        models_text += "\n\n"

        models_text += indent + "Populations:\n"

        for population in model.populations:
            if population.allow_samples:
                models_text += indent * 2
                models_text += f"{population.name}: {population.description}\n"
        models_text += "\n"

    return models_text


class HelpModels(argparse.Action):
    """
    Action used to produce model help text.
    """
    def __call__(self, parser, namespace, values, option_string=None):
        help_text = get_models_help(namespace.species, values)
        print(help_text, file=sys.stderr)
        parser.exit()


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


def write_citations(contig, model):
    """
    Write out citation information so that the user knows what papers to cite
    for the simulation engine, the model and the mutation/recombination rate
    information.
    """
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
    print("Simulation model:", model.name)
    for citation in model.citations:
        print("\t", citation, sep="")


def add_output_argument(parser):
    parser.add_argument(
        "output",
        help="The file to write simulated tree sequence to")


def summarise_usage():
    # Don't report usage on Windows as the resource module is not available.
    #  We could do this using the psutil external library, if demand exists.
    if _resource_module_available:
        rusage = resource.getrusage(resource.RUSAGE_SELF)
        user_time = humanize.naturaldelta(rusage.ru_utime)
        sys_time = rusage.ru_stime
        max_mem = rusage.ru_maxrss
        if sys.platform != 'darwin':
            max_mem *= 1024  # Linux and other OSs (e.g. freeBSD) report maxrss in kb
        max_mem_str = humanize.naturalsize(max_mem, binary=True)
        logger.info("rusage: user={}; sys={:.2f}s; max_rss={}".format(
            user_time, sys_time, max_mem_str))


def add_simulate_species_parser(parser, species):
    header = (
        f"Run simulations for {species.name} using up-to-date genome information, "
        "genetic maps and simulation models from the literature."
    )

    description_text = textwrap.fill(header)

    species_parser = parser.add_parser(
        f"{species.id}",
        description=description_text,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help=f"Run simulations for {species.name}.")
    species_parser.set_defaults(species=species.id)
    species_parser.set_defaults(genetic_map=None)
    species_parser.set_defaults(chromosome=None)
    species_parser.add_argument(
        "--help-models", action=HelpModels, nargs="?",
        help=(
            "Print descriptions of simulation models and exit. If a model ID "
            "is provided as an argument show help for this model; otherwise "
            "show help for all available models"))
    species_parser.add_argument(
        "-q", "--quiet", action='store_true',
        help="Do not write out citation information")

    # Set metavar="" to prevent help text from writing out the explicit list
    # of options, which can be too long and ugly.
    choices = [gm.name for gm in species.genetic_maps]
    if len(species.genetic_maps) > 0:
        species_parser.add_argument(
            "-g", "--genetic-map",
            choices=choices, metavar="", default=None,
            help=(
                "Specify a particular genetic map. If no genetic map is specified "
                "use a flat map by default. Available maps: "
                f"{', '.join(choices)}. "))

    if len(species.genome.chromosomes) > 1:
        choices = [chrom.id for chrom in species.genome.chromosomes]
        species_parser.add_argument(
            "-c", "--chromosome", choices=choices, metavar="", default=choices[0],
            help=(
                f"Simulate a specific chromosome. "
                f"Options: {', '.join(choices)}. "
                f"Default={choices[0]}."))
    species_parser.add_argument(
        "-l", "--length-multiplier", default=1, type=float,
        help="Simulate a chromsome of length l times the named chromosome")
    species_parser.add_argument(
        "-s", "--seed", default=None, type=int,
        help=(
            "The random seed to use for simulations. If not specified a "
            "high-quality random seed will be generated automatically. "
            "For msprime, seeds must be > 0 and < 2^32."))

    model_help = (
        "Specify a simulation model. If no model is specified, a single population"
        "constant size model is used. Available models:"
        f"{', '.join(model.id for model in species.models)}"
        ". Please see --help-models for details of these models.")
    species_parser.add_argument(
        "-m", "--model", default=None, metavar="",
        choices=[model.id for model in species.models],
        help=model_help)

    species_parser.add_argument(
        "samples", type=int, nargs="+",
        help=(
            "The number of samples to draw from each population. At least "
            "two samples must be specified. The number of arguments that "
            "will be accepted depends on the simulation model that is "
            "specified: for a model that has n populations, we can specify "
            "the number of samples to draw from each of these populations."
            "We do not need to provide sample numbers of each of the "
            "populations; those that are omitted are set to zero."))
    species_parser.add_argument(
        "output",
        help="Where to write the output tree sequence file.")

    def run_simulation(args):
        if args.model is None:
            model = stdpopsim.PiecewiseConstantSize(species.population_size)
            model.citations = species.population_size_citations
        else:
            model = species.get_model(args.model)
        if len(args.samples) > model.num_sampling_populations:
            exit(
                f"Cannot sample from more than {model.num_sampling_populations} "
                "populations")
        samples = model.get_samples(*args.samples)

        contig = species.get_contig(
            args.chromosome, genetic_map=args.genetic_map,
            length_multiplier=args.length_multiplier)
        logger.info(
            f"Running simulation model {model.name} for {species.name} on "
            f"{contig} with {len(samples)} samples")
        ts = model.simulate(contig, samples, args.seed)
        summarise_usage()
        write_output(ts, args)
        if not args.quiet:
            write_citations(contig, model)

    species_parser.set_defaults(runner=run_simulation)


def run_download_genetic_maps(args):
    species_names = [args.species]
    if args.species is None:
        species_names = [species.id for species in stdpopsim.all_species()]
    for species_id in species_names:
        species = stdpopsim.get_species(species_id)
        if len(args.genetic_maps) == 0:
            genetic_maps = [gmap.name for gmap in species.genetic_maps]
        else:
            genetic_maps = args.genetic_maps
        for genetic_map_id in genetic_maps:
            genetic_map = species.get_genetic_map(genetic_map_id)
            genetic_map.download()


def stdpopsim_cli_parser():

    # TODO the CLI defined by this hierarchical and clumsy, but it's the best
    # I could figure out. It can definitely be improved!
    top_parser = argparse.ArgumentParser(
        description="Command line interface for stdpopsim.")
    top_parser.add_argument(
        "-V", "--version", action='version',
        version='%(prog)s {}'.format(stdpopsim.__version__))
    top_parser.add_argument(
        "-v", "--verbosity", action='count', default=0,
        help="Increase the verbosity")
    top_parser.add_argument(
        "-c", "--cache-dir", type=str, default=None,
        help=(
            "Set the cache directory to the specified value. "
            "Note that this can also be set using the environment variable "
            "STDPOPSIM_CACHE. If both the environment variable and this "
            "option are set, the option takes precedence. "
            f"Default: {stdpopsim.get_cache_dir()}"))

    subparsers = top_parser.add_subparsers(dest="subcommand")
    subparsers.required = True

    for species in stdpopsim.all_species():
        add_simulate_species_parser(subparsers, species)

    download_maps_parser = subparsers.add_parser(
        "download-genetic-maps",
        help="Download genetic maps",
        description=(
            "Download genetic maps and store them in the cache directory. "
            "Maps are downloaded regardless of whether they are already "
            "in the cache or not. Please use the --cache-dir option to "
            "download maps to a specific directory. "))
    download_maps_parser.add_argument(
        "species", nargs="?",
        help=(
            "Download genetic maps for this species. If not specified "
            "download all known genetic maps."))
    download_maps_parser.add_argument(
        "genetic_maps", type=str, nargs="*",
        help=(
            "If specified, download these genetic maps. If no maps "
            "are provided, download all maps for this species."))

    download_maps_parser.set_defaults(runner=run_download_genetic_maps)

    return top_parser


# This function only exists to make mocking out the actual running of
# the program easier.
def run(args):
    args.runner(args)


def stdpopsim_main(arg_list=None):
    parser = stdpopsim_cli_parser()
    args = parser.parse_args(arg_list)
    setup_logging(args)
    if args.cache_dir is not None:
        stdpopsim.set_cache_dir(args.cache_dir)
    run(args)
