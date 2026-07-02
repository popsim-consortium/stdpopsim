import stdpopsim
import numpy as np


traits = [
    stdpopsim.Trait(id="add1", type="additive"),
    stdpopsim.Trait(id="add2", type="additive"),
    stdpopsim.Trait(id="mult", type="multiplicative"),
]
tm = stdpopsim.TraitsModel(
    traits=traits,
)


tm.add_environment(
    id="env1",
    trait_ids=["add1", "add2"],
    distribution_type="mvn",
    distribution_args=[np.zeros(2), np.eye(2)],
    time_intervals=[(0, 10)],
    population_list=[0]
)


tm.add_environment(
    id="env2",
    trait_ids=["add1", "add2"],
    distribution_type="mvn",
    distribution_args=[2*np.zeros(2), 2*np.eye(2)]
)


tm.add_environment(
    id="env3",
    trait_ids=["add1", "add2"],
    distribution_type="mvn",
    distribution_args=[2*np.zeros(2), 2*np.eye(2)],
    time_intervals=[(10, float('inf'))]
)


tm.add_fitness_function(
    id="fit1",
    trait_ids=["add1", "add2"],
    function_type="gaussian",
    function_args=[np.zeros(2), np.eye(2)],
    time_intervals=[(0, 10)],
    population_list=[0]
)

tm.add_fitness_function(
    id="fit2",
    trait_ids=["add1"],
    function_type="gaussian",
    function_args=[np.array([0.3]), np.array([[2.0]])],
)

tm.add_fitness_function(
    id="fit3",
    trait_ids=["add1"],
    function_type="gaussian",
    function_args=[np.array([0.3]), np.array([[2.0]])],
    time_intervals=[(10, float('inf'))]
)

model = stdpopsim.PiecewiseConstantSize(10)
engine = stdpopsim.get_engine("slim")
species = stdpopsim.get_species("AnaPla")
contig = species.get_contig("chr1", left=1e7, right=1e7 + 1e4)

mt1 = stdpopsim.MutationType(
    trait_ids=['add1', 'add2'],
    distribution_type='mvn',
    distribution_args=[np.zeros(2), np.eye(2)],
)

mt2 = stdpopsim.MutationType(
    trait_ids=['mult'],
    distribution_type='e',
    distribution_args=[1.],
)

dme = stdpopsim.DistributionOfMutationEffects(
    mutation_types=[mt1, mt2],
    proportions=[0.3, 0.7]
)
contig.add_dme(intervals=np.array([[1e7, 1e7+1e4]]), DME=dme)

ts = engine.simulate(
    model, contig, samples={"pop_0": 3}, traits_model=tm, seed=7,
    slim_script=True
)
