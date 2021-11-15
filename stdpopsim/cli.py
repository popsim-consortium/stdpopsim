"""
The command line interface for stdpopsim. Provides standard simulations
at the command line and methods to manage resources used by stdpopsim.
"""
import argparse
import json
import logging
import warnings
import platform
import sys
import textwrap
import tempfile
import pathlib
import shutil
import os
import re
import inspect

import msprime
import tskit
import humanize

import stdpopsim
import numpy as np

# resource is from the standard library, but it's not available on
# Windows. We break from the usual import grouping conventions here
# to avoid lots of pep8 complaints about mixing imports and code.
_resource_module_available = False
try:
    import resource

    _resource_module_available = True
except ImportError:
    pass


IS_WINDOWS = sys.platform.startswith("win")

logger = logging.getLogger(__name__)


def exit(message):
    """
    Exit with the specified error message, setting error status.
    """
    sys.exit(f"{sys.argv[0]}: {message}")


class CLIFormatter(logging.Formatter):
    # Log levels
    # ERROR: only output errors. This is the level when --quiet is specified.
    # WARN: Output warnings or any other messages that the user should be aware of
    # INFO: Write out log messages for things the user might be interested in.
    # DEBUG: Write out debug messages for the user/developer.
    def format(self, record):
        if record.name == "py.warnings":
            # trim the ugly warnings.warn message
            match = re.search(
                r"Warning:\s*(.*?)\s*warnings.warn\(", record.args[0], re.DOTALL
            )
            record.args = (match.group(1),)
            self._style = logging.PercentStyle("WARNING: %(message)s")
        else:
            if record.levelno == logging.WARNING:
                self._style = logging.PercentStyle("%(message)s")
            else:
                self._style = logging.PercentStyle("%(levelname)s: %(message)s")
        return super().format(record)


def setup_logging(args):
    log_level = "ERROR"
    if args.verbose == 1:
        log_level = "WARN"
    elif args.verbose == 2:
        log_level = "INFO"
    elif args.verbose > 2:
        log_level = "DEBUG"
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(CLIFormatter())
    logging.basicConfig(handlers=[handler], level=log_level)
    logging.captureWarnings(True)


def get_species_wrapper(species_id):
    try:
        return stdpopsim.get_species(species_id)
    except ValueError as ve:
        exit(str(ve))


def get_model_wrapper(species, model_id):
    try:
        return species.get_demographic_model(model_id)
    except ValueError as ve:
        exit(str(ve))


def get_genetic_map_wrapper(species, genetic_map_id):
    try:
        return species.get_genetic_map(genetic_map_id)
    except ValueError as ve:
        exit(str(ve))


def get_annotation_wrapper(species, annotation_id):
    try:
        return species.get_annotations(annotation_id)
    except ValueError as ve:
        exit(str(ve))


def get_dfe_wrapper(species, dfe_id):
    try:
        return species.get_dfe(dfe_id)
    except ValueError as ve:
        exit(str(ve))


def get_models_help(species_id, model_id):
    """
    Generate help text for the specified species. If model_id is None, generate
    help for all models. Otherwise, it must be a string with a valid model ID.
    """
    species = stdpopsim.get_species(species_id)
    if model_id is None:
        models_text = f"\nAll simulation models for {species.name}\n\n"
        models = [model.id for model in species.demographic_models]
    else:
        models = [model_id]
        models_text = "\nModel description\n\n"

    # TODO improve this text formatting.
    indent = " " * 4
    wrapper = textwrap.TextWrapper(initial_indent=indent, subsequent_indent=indent)
    for model_id in models:
        model = get_model_wrapper(species, model_id)
        models_text += f"{model.id}: {model.description}\n"
        models_text += wrapper.fill(textwrap.dedent(model.long_description))
        models_text += "\n\n"

        models_text += indent + "Populations:\n"

        for population in model.populations:
            if population.allow_samples:
                models_text += indent * 2
                models_text += f"{population.id}: {population.description}\n"
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


def get_genetic_maps_help(species_id, genetic_map_id):
    """
    Generate help text for the given genetic map. If map_id is None, generate
    help for all genetic maps. Otherwise, it must be a string with a valid map
    ID.
    """
    species = stdpopsim.get_species(species_id)
    if genetic_map_id is None:
        maps_text = f"\nAll genetic maps for {species.name}\n\n"
        maps = [genetic_map.id for genetic_map in species.genetic_maps]
    else:
        maps = [genetic_map_id]
        maps_text = "\nGenetic map description\n\n"

    indent = " " * 4
    wrapper = textwrap.TextWrapper(initial_indent=indent, subsequent_indent=indent)
    for map_id in maps:
        gmap = get_genetic_map_wrapper(species, map_id)
        maps_text += f"{gmap.id}\n"
        maps_text += wrapper.fill(textwrap.dedent(gmap.long_description))
        maps_text += "\n\n"

    return maps_text


class HelpGeneticMaps(argparse.Action):
    """
    Action used to produce genetic map help text.
    """

    def __call__(self, parser, namespace, values, option_string=None):
        help_text = get_genetic_maps_help(namespace.species, values)
        print(help_text, file=sys.stderr)
        parser.exit()


def get_dfes_help(species_id, dfe_id):
    """
    Generate help text for the given DFE. If dfe_id is None, generate
    help for all DFEs. Otherwise, it must be a string with a valid DFE
    ID.
    """
    species = stdpopsim.get_species(species_id)
    if dfe_id is None:
        dfes_text = f"\nAll DFEs for {species.name}\n\n"
        dfes = [dfe.id for dfe in species.dfes]
    else:
        dfes = [dfe_id]
        dfes_text = "\nDFE description\n\n"

    indent = " " * 4
    wrapper = textwrap.TextWrapper(initial_indent=indent, subsequent_indent=indent)
    for dfe_id in dfes:
        gdfe = get_dfe_wrapper(species, dfe_id)
        dfes_text += f"{gdfe.id}\n"
        dfes_text += wrapper.fill(textwrap.dedent(gdfe.long_description))
        dfes_text += "\n\n"

    return dfes_text


class HelpDFEs(argparse.Action):
    """
    Action used to produce DFE help text.
    """

    def __call__(self, parser, namespace, values, option_string=None):
        help_text = get_dfes_help(namespace.species, values)
        print(help_text, file=sys.stderr)
        parser.exit()


def get_species_help(species_id):
    """
    Generate help text for the given species with some of the species attributes
    that are not covered by the other helps.
    """
    species = stdpopsim.get_species(species_id)
    species_text = f"\nDefault population parameters for {species.name}:\n"

    species_text += f"Generation time: {species.generation_time}\n"
    species_text += f"Population size: {species.population_size}\n"
    species_text += f"Mutation rate: {species.genome.mean_mutation_rate:.4g}\n"
    species_text += (
        f"Recombination rate: {species.genome.mean_recombination_rate:.4g}\n"
    )
    species_text += (
        f"Gene conversion rate: {species.genome.mean_gene_conversion_rate:.4g}\n"
    )
    species_text += (
        f"Gene conversion tract length: "
        f"{species.genome.range_gene_conversion_lengths}\n"
    )

    return species_text


class HelpAnnotations(argparse.Action):
    """
    Action used to produce genetic map help text.
    """

    def __call__(self, parser, namespace, values, option_string=None):
        help_text = get_annotations_help(namespace.species, values)
        print(help_text, file=sys.stderr)
        parser.exit()


def get_annotations_help(species_id, annotation_id):
    """
    Generate help text for the given annotation. If annotation_id is None, generate
    help for all annotations. Otherwise, it must be a string with a valid
    annotation ID.
    """
    species = stdpopsim.get_species(species_id)
    if annotation_id is None:
        annotation_text = f"\nAll annotations for {species.name}\n\n"
        annotations = [annotation.id for annotation in species.annotations]
    else:
        annotations = [annotation_id]
        annotation_text = "\nAnnotation  description\n\n"

    indent = " " * 4
    wrapper = textwrap.TextWrapper(initial_indent=indent, subsequent_indent=indent)
    for an_id in annotations:
        annot = get_annotation_wrapper(species, an_id)
        annotation_text += f"{annot.id}\n"
        annotation_text += wrapper.fill(textwrap.dedent(annot.description))
        annotation_text += "\n\n"

    return annotation_text


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
        },
    }
    return env


def get_provenance_dict():
    """
    Returns a dictionary encoding an execution of stdpopsim conforming to the
    tskit provenance schema.
    """
    document = {
        "schema_version": "1.0.0",
        "software": {"name": "stdpopsim", "version": stdpopsim.__version__},
        "parameters": {"command": sys.argv[0], "args": sys.argv[1:]},
        "environment": get_environment(),
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
    if args.output is None:
        # There's no way to get tskit to write directly to stdout, so we write
        # to a tempfile first.
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpfile = pathlib.Path(tmpdir) / "tmp.trees"
            ts.dump(tmpfile)
            with open(tmpfile, "rb") as f:
                shutil.copyfileobj(f, sys.stdout.buffer)
    else:
        logger.debug(f"Writing to {args.output}")
        ts.dump(args.output)


def get_citations(engine, model, contig, species, dfe):
    """
    Return a list of all the citations.
    """
    citations = [stdpopsim.citations._stdpopsim_citation]
    citations.extend(engine.citations)
    citations.extend(species.citations)
    citations.extend(species.genome.citations)
    if contig.genetic_map is not None:
        citations.extend(contig.genetic_map.citations)
    citations.extend(model.citations)
    if dfe is not None:
        citations.extend(dfe.citations)
    return stdpopsim.Citation.merge(citations)


def write_bibtex(engine, model, contig, species, bibtex_file, dfe):
    """
    Write bibtex for available citations to a file."""
    citations = get_citations(engine, model, contig, species, dfe)
    for citation in citations:
        bibtex_file.write(citation.fetch_bibtex())
    bibtex_file.close()


def write_citations(engine, model, contig, species, dfe):
    """
    Write out citation information so that the user knows what papers to cite
    for the simulation engine, the model and the mutation/recombination rate
    information.
    """
    cite_str = ["If you use this simulation in published work, please cite:"]
    for citation in get_citations(engine, model, contig, species, dfe):
        if (
            stdpopsim.citations.CiteReason.MUT_RATE in citation.reasons
            and model.mutation_rate is not None
        ):
            continue
        cite_str.append("\n")
        cite_str.append(citation.displaystr())
    logger.warning("".join(cite_str))


def summarise_usage():
    # Don't report usage on Windows as the resource module is not available.
    #  We could do this using the psutil external library, if demand exists.
    if _resource_module_available:
        rusage = resource.getrusage(resource.RUSAGE_SELF)
        user_time = humanize.naturaldelta(rusage.ru_utime)
        sys_time = rusage.ru_stime
        max_mem = rusage.ru_maxrss
        if sys.platform != "darwin":
            max_mem *= 1024  # Linux and other OSs (e.g. freeBSD) report maxrss in kb
        max_mem_str = humanize.naturalsize(max_mem, binary=True)
        logger.info(
            "rusage: user={}; sys={:.2f}s; max_rss={}".format(
                user_time, sys_time, max_mem_str
            )
        )


def add_simulate_species_parser(parser, species):
    header = (
        f"Run simulations for {species.name} using up-to-date genome information, "
        "genetic maps and simulation models from the literature. "
        "NOTE: By default, the tskit '.trees' binary file is written to stdout,"
        "so you should either redirect this to a file or use the '--output' "
        "option to specify a filename."
    )

    description_text = textwrap.fill(header) + "\n" + get_species_help(species.id)

    species_parser = parser.add_parser(
        f"{species.id}",
        description=description_text,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help=f"Run simulations for {species.name}.",
    )
    species_parser.set_defaults(species=species.id)
    species_parser.set_defaults(genetic_map=None, dfe_annotation=None)
    species_parser.add_argument(
        "--help-models",
        action=HelpModels,
        nargs="?",
        help=(
            "Print descriptions of simulation models and exit. If a model ID "
            "is provided as an argument show help for this model; otherwise "
            "show help for all available models"
        ),
    )
    species_parser.add_argument(
        "-b",
        "--bibtex-file",
        type=argparse.FileType("w"),
        help="Write citations to a given bib file. This will overwrite the file.",
        default=None,
        action="store",
    )

    # Set metavar="" to prevent help text from writing out the explicit list
    # of options, which can be too long and ugly.
    choices = [gm.id for gm in species.genetic_maps]
    if len(species.genetic_maps) > 0:
        species_parser.add_argument(
            "--help-genetic-maps",
            action=HelpGeneticMaps,
            nargs="?",
            help=(
                "Print list of genetic maps and exit. If a genetic map ID is "
                "given as an argument, show help for this map. Otherwise show "
                "help for all available genetic maps"
            ),
        )

    species_parser.add_argument(
        "-D",
        "--dry-run",
        action="store_true",
        default=False,
        help="Do not run actual simulation",
    )

    if len(species.genetic_maps) > 0:
        species_parser.add_argument(
            "-g",
            "--genetic-map",
            choices=choices,
            metavar="",
            default=None,
            help=(
                "Specify a particular genetic map. By default, a chromosome-specific "
                "uniform recombination rate is used. These default rates are listed in "
                "the catalog: "
                "<https://popsim-consortium.github.io/stdpopsim-docs/"
                "stable/catalog.html> "
                "Available maps: "
                f"{', '.join(choices)}. "
            ),
        )

    if len(species.dfes) > 0:
        species_parser.add_argument(
            "--help-dfes",
            action=HelpDFEs,
            nargs="?",
            help=(
                "Print list of DFEs and exit. If a DFE ID is "
                "given as an argument, show help for this DFE. Otherwise show "
                "help for all available DFEs"
            ),
        )

    # To avoid listing too much stuff out in the help, we only list
    # the actual IDs. We make all synonyms available as choices though.
    choices = []
    all_choices = []
    for chrom in species.genome.chromosomes:
        choices.append(chrom.id)
        all_choices.extend([chrom.id] + chrom.synonyms)
    species_parser.add_argument(
        "-c",
        "--chromosome",
        choices=all_choices,
        metavar="",
        default=None,
        help=(
            f"Simulate a specific chromosome. If no chromosome is given, "
            f"simulate a generic contig of given length, specified using --length. "
            f"Options: {', '.join(choices)}. "
            f"Default=None."
        ),
    )

    species_parser.add_argument(
        "-L",
        "--length",
        default=None,
        type=float,
        help="Simulate a default contig of given length.",
    )
    species_parser.add_argument(
        "-i",
        "--inclusion-mask",
        default=None,
        type=str,
        help="Path to inclusion mask specified in bed format.",
    )
    species_parser.add_argument(
        "-e",
        "--exclusion-mask",
        default=None,
        type=str,
        help="Path to exclusion mask specified in bed format.",
    )
    species_parser.add_argument(
        "-l",
        "--length-multiplier",
        default=1,
        type=float,
        help="Simulate a sequence of length l times the named chromosome's length, "
        "using the named chromosome's mutation and recombination rates.",
    )
    species_parser.add_argument(
        "-s",
        "--seed",
        default=None,
        type=int,
        help=(
            "The random seed to use for simulations. If not specified a "
            "high-quality random seed will be generated automatically. "
            "For msprime, seeds must be > 0 and < 2^32."
        ),
    )

    model_help = (
        "Specify a simulation model. If no model is specified, a single population "
        "constant size model is used. Available models: "
        f"{', '.join(model.id for model in species.demographic_models)}"
        ". Please see --help-models for details of these models."
    )
    species_parser.add_argument(
        "-d",
        "--demographic-model",
        default=None,
        metavar="",
        choices=[model.id for model in species.demographic_models],
        help=model_help,
    )
    dfe_help = (
        "Specify a Distribution of Fitness Effects (DFE) model. "
        "If no DFE is specified, all mutations are neutral. "
        "Available DFE models:"
        f"{', '.join(dfe.id for dfe in species.dfes)}"
        ". Please see --help-dfes for details of these DFE models."
    )
    species_parser.add_argument(
        "--dfe",
        default=None,
        type=str,
        metavar="",
        choices=[dfe.id for dfe in species.dfes],
        help=dfe_help,
    )
    interval_help = (
        "Specify the interval where selection (given a DFE) is simulated. "
        "Anything outside of the interval is simulated as neutral "
        "If no interval is specified, "
        "selection is simulated across the entire region (contig). "
        "Interval (for now) is writen as 'left,right' (separated by a comma, "
        "with no space, for instance: --dfe-interval 1000,2000"
        # say something here about multiple dfes?
    )
    species_parser.add_argument(
        "--dfe-interval",
        default=None,
        type=str,
        metavar="",
        help=interval_help,
    )
    bed_help = (
        "A bed file specifing the intervals where selection (given a DFE) is simulated. "
        "Non-overlaping intervals belonging to the same chromosome are required. "
        "If no interval is specified, "
        "selection is simulated across the entire contig. "
        "See also --dfe-interval and --dfe-annotation."
    )
    species_parser.add_argument(
        "--dfe-bed-file",
        default=None,
        type=str,
        metavar="",
        help=bed_help,
    )
    annot_choices = [gm.id for gm in species.annotations]
    if len(species.annotations) > 0:
        species_parser.add_argument(
            "--help-annotations",
            action=HelpAnnotations,
            nargs="?",
            help=(
                "Print list of annotations and exit. If a annotation ID is "
                "given as an argument, show help for this annotation. Otherwise show "
                "help for all available annotations"
            ),
        )
    if len(species.annotations) > 0:
        species_parser.add_argument(
            "--dfe-annotation",
            choices=annot_choices,
            metavar="",
            default=None,
            help=(
                "Specify the intervals over which a given DFE acts "
                "using a genomic annotation track "
                "Available maps: "
                f"{', '.join(annot_choices)}. "
            ),
        )

    species_parser.add_argument(
        "-o",
        "--output",
        help=(
            "Where to write the output tree sequence file. Defaults to "
            "stdout if not specified"
        ),
    )

    species_parser.add_argument(
        "samples",
        type=int,
        nargs="+",
        help=(
            "The number of samples to draw from each population. At least "
            "two samples must be specified. The number of arguments that "
            "will be accepted depends on the simulation model that is "
            "specified: for a model that has n populations, we can specify "
            "the number of samples to draw from each of these populations."
            "We do not need to provide sample numbers of each of the "
            "populations; those that are omitted are set to zero."
        ),
    )

    def run_simulation(args):
        if args.demographic_model is None:
            model = stdpopsim.PiecewiseConstantSize(species.population_size)
            model.generation_time = species.generation_time
            for citation in species.citations:
                reasons = {
                    stdpopsim.CiteReason.POP_SIZE,
                    stdpopsim.CiteReason.GEN_TIME,
                }
                if len(citation.reasons & reasons) > 0:
                    model.citations.append(citation)
            qc_complete = True
        else:
            model = get_model_wrapper(species, args.demographic_model)
            qc_complete = model.qc_model is not None
        if len(args.samples) > model.num_sampling_populations:
            exit(
                f"Cannot sample from more than {model.num_sampling_populations} "
                "populations"
            )
        samples = model.get_samples(*args.samples)
        contig = species.get_contig(
            args.chromosome,
            genetic_map=args.genetic_map,
            length_multiplier=args.length_multiplier,
            length=args.length,
            inclusion_mask=args.inclusion_mask,
            exclusion_mask=args.exclusion_mask,
            mutation_rate=model.mutation_rate,
        )
        engine = stdpopsim.get_engine(args.engine)
        logger.info(
            f"Running simulation model {model.id} for {species.id} on "
            f"{contig} with {len(samples)} samples using {engine.id}."
        )
        # DFE assignment
        dfe = None
        intervals_summary_str = None
        if args.dfe is None:
            if args.dfe_interval is not None:
                exit(
                    "A DFE interval has been assigned without a DFE. "
                    "Please specify a DFE."
                )
            if args.dfe_annotation is not None:
                exit(
                    "A DFE annotation has been assigned without a DFE. "
                    "Please specify a DFE."
                )
            if args.dfe_bed_file is not None:
                exit(
                    "A DFE bed file has been assigned without a DFE. "
                    "Please specify a DFE."
                )
        else:
            if args.dfe_interval is not None:
                if args.dfe_annotation is not None:
                    exit(
                        "A DFE annotation and a DFE interval have been "
                        "selected. Please only use one."
                    )
                if args.dfe_bed_file is not None:
                    exit(
                        "A DFE bed file and a DFE interval have been "
                        "selected. Please only use one."
                    )
                left, right = args.dfe_interval.split(",")
                intervals = np.array([[int(left), int(right)]])
                intervals_summary_str = f"[{left}, {right})"
            if args.dfe_annotation is not None:
                if args.dfe_bed_file is not None:
                    exit(
                        "A DFE bed file and a DFE annotation have been "
                        "selected. Please only use one."
                    )
                annot = species.get_annotations(args.dfe_annotation)
                intervals = annot.get_chromosome_annotations(args.chromosome)
                intervals_summary_str = f"{annot.id} elements on {args.chromosome}"
            if args.dfe_bed_file is not None:
                intervals = np.loadtxt(args.dfe_bed_file, usecols=[1, 2], dtype="int")
                left = np.min(intervals)
                right = np.max(intervals)
                intervals_summary_str = f"[{left}, {right})"
            if intervals_summary_str is None:
                # case where no intervals specified but we have a DFE
                intervals = np.array(
                    [[0, int(contig.recombination_map.sequence_length)]]
                )
                intervals_summary_str = f"[{intervals[0][0]}, {intervals[0][1]})"

            dfe = species.get_dfe(args.dfe)
            contig.add_dfe(
                intervals=intervals,
                DFE=dfe,
            )
            logger.info(
                f"Applying selection under the DFE model {dfe.id} "
                f"in intervals {intervals_summary_str}."
            )
        write_simulation_summary(
            engine=engine,
            model=model,
            contig=contig,
            samples=samples,
            dfe=args.dfe,
            dfe_interval=intervals_summary_str,
            seed=args.seed,
        )
        if not qc_complete:
            warnings.warn(
                stdpopsim.QCMissingWarning(
                    f"{model.id} has not been QCed. Use at your own risk! "
                    "Demographic models that have not undergone stdpopsim's "
                    "Quality Control procedure may contain implementation "
                    "errors, leading to differences between simulations "
                    "and the model described in the original publication. "
                    "More information about the QC process can be found in "
                    "the developer documentation. "
                    "https://popsim-consortium.github.io/stdpopsim-docs/"
                    "latest/development.html#demographic-model-review-process"
                )
            )

        # extract simulate() parameters from CLI args
        accepted_params = inspect.signature(engine.simulate).parameters.keys()
        kwargs = {k: v for k, v in vars(args).items() if k in accepted_params}
        kwargs.update(demographic_model=model, contig=contig, samples=samples)
        ts = engine.simulate(**kwargs)

        summarise_usage()
        if ts is not None:
            write_output(ts, args)
        # Non-QCed models shouldn't be used in publications, so we skip the
        # "If you use this simulation in published work..." citation request.
        if qc_complete:
            write_citations(engine, model, contig, species, dfe)
        if args.bibtex_file is not None:
            write_bibtex(engine, model, contig, species, args.bibtex_file, dfe)

    species_parser.set_defaults(runner=run_simulation)


def write_simulation_summary(
    engine, model, contig, samples, dfe, dfe_interval, seed=None
):
    indent = " " * 4
    # Header
    dry_run_text = "Simulation information:\n"
    # Get information about the engine
    dry_run_text += f"{indent}Engine: {engine.id} ({engine.get_version()})\n"
    # Get information about model
    dry_run_text += f"{indent}Model id: {model.id}\n"
    dry_run_text += f"{indent}Model description: {model.description}\n"
    # Add seed information if extant
    if seed is not None:
        dry_run_text += f"{indent}Seed: {seed}\n"
    # Get information about the number of samples
    dry_run_text += f"{indent}Population: number_samples (sampling_time_generations):\n"
    sample_counts = [0] * model.num_sampling_populations
    for s in samples:
        sample_counts[s.population] += s.num_samples
    for p in range(0, model.num_sampling_populations):
        pop_name = model.populations[p].id
        sample_time = model.populations[p].default_sampling_time
        dry_run_text += f"{indent * 2}{pop_name}: "
        dry_run_text += f"{sample_counts[p]} ({sample_time})\n"
    # Get information about relevant contig
    gmap = "None" if contig.genetic_map is None else contig.genetic_map.id
    mean_recomb_rate = contig.recombination_map.mean_rate
    # use the model mutation rate, if provided with the demographic model
    if model.mutation_rate is not None:
        mut_rate = model.mutation_rate
        logging.info(
            f"using mutation rate from demographic model ({model.mutation_rate})"
        )
    else:
        mut_rate = contig.mutation_rate
        logging.info(
            f"using mutation rate from species contig ({contig.mutation_rate})"
        )
    contig_len = contig.recombination_map.sequence_length
    dry_run_text += "Contig Description:\n"
    dry_run_text += f"{indent}Contig length: {contig_len}\n"
    dry_run_text += f"{indent}Mean recombination rate: {mean_recomb_rate}\n"
    dry_run_text += f"{indent}Mean mutation rate: {mut_rate}\n"
    dry_run_text += f"{indent}Genetic map: {gmap}\n"
    if dfe is not None:
        dry_run_text += f"{indent}DFE: {dfe}\n"
        dry_run_text += f"{indent}DFE applied to: {dfe_interval}\n"
    logger.warning(dry_run_text)


def run_download_genetic_maps(args):
    species_names = [args.species]
    if args.species is None:
        species_names = [species.id for species in stdpopsim.all_species()]
    for species_id in species_names:
        species = get_species_wrapper(species_id)
        if len(args.genetic_maps) == 0:
            genetic_maps = [gmap.id for gmap in species.genetic_maps]
        else:
            genetic_maps = args.genetic_maps
        for genetic_map_id in genetic_maps:
            genetic_map = get_genetic_map_wrapper(species, genetic_map_id)
            logger.warning(f"Downloading map {species_id}/{genetic_map_id}")
            genetic_map.download()


def stdpopsim_cli_parser():
    class QuietAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            namespace.verbose = 0

    top_parser = argparse.ArgumentParser(
        description="Command line interface for stdpopsim."
    )
    top_parser.add_argument(
        "-V",
        "--version",
        action="version",
        version="%(prog)s {}".format(stdpopsim.__version__),
    )
    logging_group = top_parser.add_mutually_exclusive_group()
    logging_group.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=1,
        help="Increase logging verbosity (can use be used multiple times).",
    )
    logging_group.add_argument(
        "-q",
        "--quiet",
        nargs=0,
        help="Do not write any non-essential messages",
        action=QuietAction,
    )

    top_parser.add_argument(
        "-c",
        "--cache-dir",
        type=str,
        default=None,
        help=(
            "Set the cache directory to the specified value. "
            "Note that this can also be set using the environment variable "
            "STDPOPSIM_CACHE. If both the environment variable and this "
            "option are set, the option takes precedence. "
            f"Default: {stdpopsim.get_cache_dir()}"
        ),
    )

    top_parser.add_argument(
        "-e",
        "--engine",
        default=stdpopsim.get_default_engine().id,
        choices=[e.id for e in stdpopsim.all_engines()],
        help="Specify a simulation engine.",
    )

    supported_models = stdpopsim.get_engine("msprime").supported_models
    msprime_parser = top_parser.add_argument_group("msprime specific parameters")
    msprime_parser.add_argument(
        "--msprime-model",
        default=supported_models[0],
        choices=supported_models,
        help="Specify the simulation model used by msprime. "
        "See msprime API documentation for details.",
    )

    def time_or_model(
        arg,
        _arg_is_time=[
            True,
        ],
        parser=top_parser,
    ):
        if _arg_is_time[0]:
            try:
                arg = float(arg)
            except ValueError:
                parser.error(f"`{arg}' is not a number")
        else:
            if arg not in supported_models:
                parser.error(f"`{arg}' is not a supported model")
        _arg_is_time[0] = not _arg_is_time[0]
        return arg

    msprime_parser.add_argument(
        "--msprime-change-model",
        metavar=("T", "MODEL"),
        type=time_or_model,
        default=[],
        action="append",
        nargs=2,
        help="Change to the specified simulation MODEL at generation T. "
        "This option may provided multiple times.",
    )

    # SLiM is not available for windows.
    if not IS_WINDOWS:

        def slim_exec(path):
            # Hack to set the SLIM environment variable at parse time,
            # before get_version() can be called.
            os.environ["SLIM"] = path
            return path

        slim_parser = top_parser.add_argument_group("SLiM specific parameters")
        slim_parser.add_argument(
            "--slim-path",
            metavar="PATH",
            type=slim_exec,
            default=None,
            help="Full path to `slim' executable.",
        )
        slim_parser.add_argument(
            "--slim-script",
            action="store_true",
            default=False,
            help="Write script to stdout and exit without running SLiM.",
        )
        slim_parser.add_argument(
            "--slim-scaling-factor",
            metavar="Q",
            default=1,
            type=float,
            help="Rescale model parameters by Q to speed up simulation. "
            "See SLiM manual: `5.5 Rescaling population sizes to "
            "improve simulation performance`. "
            "[default=%(default)s].",
        )
        slim_parser.add_argument(
            "--slim-burn-in",
            metavar="X",
            default=10,
            type=float,
            help="Length of the burn-in phase, in units of N generations "
            "[default=%(default)s].",
        )

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
            "download maps to a specific directory. "
        ),
    )
    download_maps_parser.add_argument(
        "species",
        nargs="?",
        help=(
            "Download genetic maps for this species. If not specified "
            "download all known genetic maps."
        ),
    )
    download_maps_parser.add_argument(
        "genetic_maps",
        type=str,
        nargs="*",
        help=(
            "If specified, download these genetic maps. If no maps "
            "are provided, download all maps for this species."
        ),
    )

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
