import random
import logging
import numpy as np
import stdpopsim

logger = logging.getLogger(__name__)


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
    samples = model.get_samples(
        100, 0, 0, 100, 2, 2  # YRI, CEU, CHB, Papuan, DenA, NeaA
    )

    # One mutation type, which we'll use for the positively selected mutation.
    # Neutral mutations will be added by the SLiM engine as usual, after the
    # SLiM phase of the simulation has completed.
    positive = stdpopsim.ext.MutationType(convert_to_substitution=False)
    mutation_types = [positive]
    mut_id = len(mutation_types) - 1

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
    coordinate = round(contig.recombination_map.get_length() / 2)

    pop = {p.id: i for i, p in enumerate(model.populations)}

    # Thinking forwards in time, we define a number of extended events that
    # correspond to drawing the mutation, conditioning on the new allele not
    # being lost, and its selection. The extended events will be time-sorted
    # by the SLiM engine, so the ordering here is only to aid clarity.
    # Like msprime.DemographicEvents, all the times used here have units of
    # generations before present.
    extended_events = [
        # Draw mutation in DenA.
        stdpopsim.ext.DrawMutation(
            time=T_mut,
            mutation_type_id=mut_id,
            population_id=pop["DenA"],
            coordinate=coordinate,
            # Save state before the mutation is introduced.
            save=True,
        ),
        # Because the drawn mutation is neutral at the time of introduction,
        # it's likely to be lost due to drift. To avoid this, we condition on
        # the mutation having AF > 0 in DenA. If this condition is false at any
        # time between start_time and end_time, the simulation will be
        # restored to the most recent save point.
        # Conditioning should start one generation after T_mut (not at T_mut!),
        # to avoid checking for the mutation before SLiM can introduce it.
        stdpopsim.ext.ConditionOnAlleleFrequency(
            # Note: if T_mut ~= T_Den_split, then we end up with:
            #       GenerationAfter(T_mut) < T_Den_split,
            #       which will give an error due to "start_time < end_time".
            start_time=stdpopsim.ext.GenerationAfter(T_mut),
            end_time=T_Den_split,
            mutation_type_id=mut_id,
            population_id=pop["DenA"],
            op=">",
            allele_frequency=0,
        ),
        # Denisovans split into DenA and Den1 at time T_Den_split,
        # so now we condition on having AF > 0 in Den1.
        stdpopsim.ext.ConditionOnAlleleFrequency(
            start_time=stdpopsim.ext.GenerationAfter(T_Den_split),
            end_time=T_mig,
            mutation_type_id=mut_id,
            population_id=pop["Den1"],
            op=">",
            allele_frequency=0,
            # Update save point at start_time (if the condition is met).
            save=True,
        ),
        # The Den1 lineage has migrants entering the Papaun lineage at T_mig,
        # so condition on AF > 0 in Papuans.
        stdpopsim.ext.ConditionOnAlleleFrequency(
            start_time=stdpopsim.ext.GenerationAfter(T_mig),
            end_time=0,
            mutation_type_id=mut_id,
            population_id=pop["Papuan"],
            op=">",
            allele_frequency=0,
            # Update save point at start_time (if the condition is met).
            # If the Den1 migrants didn't carry the mutation, we will
            # instead restore to the previous save point.
            save=True,
        ),
        # The mutation is positively selected in Papuans at T_sel.
        # Note that this will have no effect, unless/until a mutation with the
        # specified mutation_type_id is found in the population.
        stdpopsim.ext.ChangeMutationFitness(
            start_time=T_sel,
            end_time=0,
            mutation_type_id=mut_id,
            population_id=pop["Papuan"],
            selection_coeff=s,
            dominance_coeff=0.5,
        ),
        # Condition on AF > 0.05 in Papuans at the end of the simulation.
        stdpopsim.ext.ConditionOnAlleleFrequency(
            start_time=0,
            end_time=0,
            mutation_type_id=mut_id,
            population_id=pop["Papuan"],
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
        seed=rng.randrange(1, 2 ** 32),
        extended_events=extended_events,
        slim_scaling_factor=10,
        slim_burn_in=0.1,
        # Set slim_script=True to print the script instead of running it.
        # slim_script=True,
    )

    return ts, T_mut, T_sel, s


def KimDFE():
    """
    Return neutral and negative MutationType()s representing a human DFE.
    Kim et al. (2018), p.23, http://doi.org/10.1371/journal.pgen.1007741
    """
    neutral = stdpopsim.ext.MutationType()
    gamma_shape = 0.186  # shape
    gamma_mean = -0.01314833  # expected value
    h = 0.5 / (1 - 7071.07 * gamma_mean)  # dominance coefficient
    negative = stdpopsim.ext.MutationType(
        dominance_coeff=h,
        distribution_type="g",
        distribution_args=[gamma_mean, gamma_shape],
    )
    # note neutral mutations have 0 proportion because they are not simulated by SLiM
    return {"mutation_types": [neutral, negative], "proportions": [0.0, 0.7]}


def OutOfAfrica_3G09_with_DFE(seed):
    """
    The Gutenkunst et al. HomSap/OutOfAfrica_3G09 model, simulated with a DFE.
    """
    species = stdpopsim.get_species("HomSap")
    model = species.get_demographic_model("OutOfAfrica_3G09")
    contig = species.get_contig("chr1", length_multiplier=0.001)
    samples = model.get_samples(100, 100, 100)  # YRI, CEU, CHB

    # neutral and deleterious mutations occur across the whole contig
    contig.add_genomic_element_type(
        intervals=np.array([[0, int(contig.length)]]), **KimDFE()
    )

    # Simulate.
    engine = stdpopsim.get_engine("slim")
    ts = engine.simulate(
        model,
        contig,
        samples,
        seed=seed,
        slim_scaling_factor=10,
        slim_burn_in=10,
        # Set slim_script=True to print the script instead of running it.
        # slim_script=True,
    )
    return ts


def gene_with_noncoding_OutOfAfrica_3G09(seed):
    """
    Simulating a 1kb gene flanked by 1kb neutral regions. Within genes,
    30% of the total influx of mutations are neutral and 70% are deleterious,
    with the DFE from Kim et al. The HomSap/OutOfAfrica_3G09 model was simulated.
    """
    species = stdpopsim.get_species("HomSap")
    model = species.get_demographic_model("OutOfAfrica_3G09")
    contig = species.get_contig(length=3000)
    samples = model.get_samples(100, 100, 100)  # YRI, CEU, CHB

    # within the gene, KimDFE is used, outside genomic elements
    # neutral muts are added with msprime
    gene_interval = np.array([[1000, 2000]])
    contig.add_genomic_element_type(intervals=gene_interval, **KimDFE())

    # Simulate.
    engine = stdpopsim.get_engine("slim")
    ts = engine.simulate(
        model,
        contig,
        samples,
        seed=seed,
        slim_scaling_factor=10,
        slim_burn_in=10,
        # Set slim_script=True to print the script instead of running it.
        # slim_script=True,
    )
    return ts


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

    ts = OutOfAfrica_3G09_with_DFE(seed)

    ts = gene_with_noncoding_OutOfAfrica_3G09(seed)
