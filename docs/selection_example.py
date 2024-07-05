import random
import logging
import stdpopsim

logger = logging.getLogger(__name__)

# TODO: This example has been updated to reflect changes in the extended events
# API (see PRs #1306 and #1341) but it should be run and checked for
# correctness at some point


def adaptive_introgression(seed):
    """
    Adaptive introgression in HomSap/PapuansOutOfAfrica_10J19.

    A neutral mutation is drawn in Denisovans, transmitted to Papuans via a
    migration pulse, and is then positively selected in the Papuan population.
    The time of mutation introduction, the time of selection onset, and the
    selection coefficient, are each random variables.
    """
    species = stdpopsim.get_species("HomSap")
    model = species.get_demographic_model("PapuansOutOfAfrica_10J19")
    contig = species.get_contig("chr1", length_multiplier=0.001)
    samples = {"YRI": 50, "Papuan": 50, "DenA": 1, "NeaA": 1}

    # We need some demographic model parameters to set bounds on the timing
    # of random variables and extended_events (below).
    # These values were copied from the PapuansOutOfAfrica_10J19 model
    # implementation, but can also be found in the catalog documentation.
    T_Den_Nea_split = 15090
    T_DenA_Den1_split = 9750
    T_Den1_Papuan_mig = 29.8e3 / model.generation_time
    # The drawn mutation is transmitted via Den1.
    T_Den_split = T_DenA_Den1_split
    T_mig = T_Den1_Papuan_mig

    # Draw random variables.
    rng = random.Random(seed)
    t_delta = 1000 / model.generation_time
    # Time of mutation introduction. As a lower bound, we use the split time of
    # the Denisovan lineages, plus a small offset (t_delta). The offset avoids
    # an edge case, where T_mut could be rounded to T_Den_split (see note below
    # about avoiding start_time < end_time).
    T_mut = rng.uniform(T_Den_split + t_delta, T_Den_Nea_split)
    # Time of selection onset. We use a non-zero lower bound, so that
    # weaker selection has time to act.
    T_sel = rng.uniform(t_delta, T_mig)
    # Selection coefficient.
    s = rng.uniform(0.001, 0.1)
    logger.info(f"Parameters: T_mut={T_mut:.3f}, T_sel={T_sel:.3f}, s={s:.3g}")

    # Place the drawn mutation in the middle of the contig.
    locus_id = "introgressed_locus"
    coordinate = round(contig.recombination_map.sequence_length / 2)
    contig.add_single_site(id=locus_id, coordinate=coordinate)

    # Thinking forwards in time, we define a number of extended events that
    # correspond to drawing the mutation, conditioning on the new allele not
    # being lost, and its selection. The extended events will be time-sorted
    # by the SLiM engine, so the ordering here is only to aid clarity.
    # Like msprime.DemographicEvents, all the times used here have units of
    # generations before present.
    extended_events = [
        # Draw mutation in DenA.
        stdpopsim.DrawMutation(
            time=T_mut,
            single_site_id=locus_id,
            population="DenA",
        ),
        # Because the drawn mutation is neutral at the time of introduction,
        # it's likely to be lost due to drift. To avoid this, we condition on
        # the mutation having AF > 0 in DenA. If this condition is false at any
        # time between start_time and end_time, the simulation will be
        # returned to the point where the mutation was introduced.
        # Conditioning should start one generation after T_mut (not at T_mut!),
        # to avoid checking for the mutation before SLiM can introduce it.
        stdpopsim.ConditionOnAlleleFrequency(
            # Note: if T_mut ~= T_Den_split, then we end up with:
            #       GenerationAfter(T_mut) < T_Den_split,
            #       which will give an error due to "start_time < end_time".
            start_time=stdpopsim.GenerationAfter(T_mut),
            end_time=T_Den_split,
            single_site_id=locus_id,
            population="DenA",
            op=">",
            allele_frequency=0,
        ),
        # Denisovans split into DenA and Den1 at time T_Den_split,
        # so now we condition on having AF > 0 in Den1.
        stdpopsim.ConditionOnAlleleFrequency(
            start_time=stdpopsim.GenerationAfter(T_Den_split),
            end_time=T_mig,
            single_site_id=locus_id,
            population="Den1",
            op=">",
            allele_frequency=0,
        ),
        # The Den1 lineage has migrants entering the Papaun lineage at T_mig,
        # so condition on AF > 0 in Papuans.
        stdpopsim.ConditionOnAlleleFrequency(
            start_time=stdpopsim.GenerationAfter(T_mig),
            end_time=0,
            single_site_id=locus_id,
            population="Papuan",
            op=">",
            allele_frequency=0,
        ),
        # The mutation is positively selected in Papuans at T_sel.
        # Note that this will have no effect, unless/until a mutation with the
        # specified mutation_type_id is found in the population.
        stdpopsim.ChangeMutationFitness(
            start_time=T_sel,
            end_time=0,
            single_site_id=locus_id,
            population="Papuan",
            selection_coeff=s,
            dominance_coeff=0.5,
        ),
        # Condition on AF > 0.05 in Papuans at the end of the simulation.
        stdpopsim.ConditionOnAlleleFrequency(
            start_time=0,
            end_time=0,
            single_site_id=locus_id,
            population="Papuan",
            op=">",
            allele_frequency=0.05,
        ),
    ]

    # Simulate.
    engine = stdpopsim.get_engine("slim")
    ts = engine.simulate(
        model,
        contig,
        samples,
        seed=rng.randrange(1, 2**32),
        extended_events=extended_events,
        slim_scaling_factor=10,
        slim_burn_in=0.1,
        # Set slim_script=True to print the script instead of running it.
        # slim_script=True,
    )

    return ts, T_mut, T_sel, s


if __name__ == "__main__":
    import sys
    import stdpopsim.cli
    import collections

    if len(sys.argv) == 2:
        seed = int(sys.argv[1])
    else:
        seed = 1234

    # Setup logging at verbosity level 2.
    args = collections.namedtuple("_", ["quiet", "verbose"])(False, 2)
    stdpopsim.cli.setup_logging(args)

    ts, T_mut, T_sel, s = adaptive_introgression(seed)
