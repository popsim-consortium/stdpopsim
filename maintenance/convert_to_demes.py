#!/usr/bin/env python3
# Convert all stdpopsim models into demes YAML files.

import copy
import sys
import pathlib
import textwrap

import demes
import msprime
import stdpopsim


def change_units(graph: demes.Graph, time_units: str, generation_time: float):
    # TODO: put a function like this in demes.

    # stdpopsim models are always in generations
    assert graph.time_units == "generations"

    # return a copy instead of modifying the original
    graph = copy.deepcopy(graph)

    if time_units == "generations":
        assert generation_time == 1
        return graph

    for deme in graph.demes:
        deme.start_time *= generation_time
        for epoch in deme.epochs:
            epoch.start_time *= generation_time
            epoch.end_time *= generation_time
    for migration in graph.migrations:
        migration.start_time *= generation_time
        migration.end_time *= generation_time
    for pulse in graph.pulses:
        pulse.time *= generation_time
    graph.time_units = time_units
    graph.generation_time = generation_time

    # Check for stupid mistakes.
    graph2 = demes.Graph.fromdict(graph.asdict())
    graph2.assert_close(graph)

    return graph


# Convert the models to the given time_units and generation_time.
# This table was filled by inspecting the stdpopsim model definitions.
time_conversion = {
    # AnaPla
    "MallardBlackDuck_2L19": ("years", 4),
    # AnoGam
    "GAS_1A17": ("generations", 1),
    # AraTha
    "SouthMiddleAtlas_1D17": ("generations", 1),
    "African2Epoch_1H18": ("generations", 1),
    "African3Epoch_1H18": ("generations", 1),
    # BosTau
    "HolsteinFriesian_1M13": ("generations", 1),
    # DroMel
    "African3Epoch_1S16": ("generations", 1),
    "OutOfAfrica_2L06": ("generations", 1),
    # HomSap
    "OutOfAfricaExtendedNeandertalAdmixturePulse_3I21": (
        "thousands of years",
        25 / 1000,
    ),
    "OutOfAfrica_3G09": ("years", 25),
    "OutOfAfrica_2T12": ("years", 25),
    "Africa_1T12": ("years", 25),
    "AmericanAdmixture_4B11": ("generations", 1),
    "OutOfAfricaArchaicAdmixture_5R19": ("years", 29),
    "Zigzag_1S14": ("generations", 1),
    "AncientEurasia_9K19": ("years", 25),
    "PapuansOutOfAfrica_10J19": ("generations", 1),
    "AshkSub_7G19": ("generations", 1),
    "OutOfAfrica_4J17": ("years", 29),
    "Africa_1B08": ("generations", 1),
    # PanTro
    "BonoboGhost_4K19": ("thousands of years", 25 / 1000),
    # PonAbe
    "TwoSpecies_2L11": ("years", 20),
}

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"usage: {sys.argv[0]} output_folder/")
        exit(1)

    output_folder = pathlib.Path(sys.argv[1])
    output_folder.mkdir(exist_ok=True)

    for species in stdpopsim.all_species():
        species_folder = output_folder / species.id
        species_folder.mkdir(exist_ok=True)
        for model in species.demographic_models:
            graph = msprime.Demography.to_demes(model.model)

            # Change the time units.
            time_units, generation_time = time_conversion[model.id]
            graph = change_units(graph, time_units, generation_time)

            # Add description.
            graph.description = (
                model.description
                + "\n"
                + " ".join(
                    textwrap.wrap(textwrap.dedent(model.long_description))
                ).strip()
            )

            # Add citations.
            for citation in model.citations:
                graph.doi.append(str(citation))

            # Add metadata.
            if model.mutation_rate is not None:
                graph.metadata["mutation_rate"] = model.mutation_rate

            demes.dump(graph, species_folder / f"{model.id}.yaml")
