"""
Tests for classes that hold information about genomic region to be simulated.
"""

import numpy as np
import pytest
import msprime
import stdpopsim


class TestContig(object):

    example_dfe = stdpopsim.DFE(
        id="abc",
        description="example DFE",
        long_description="test test beep boop beep",
        proportions=[0.3, 0.7],
        mutation_types=[stdpopsim.MutationType() for _ in range(2)],
    )

    def verify_dfe_breakpoints(self, contig):
        breaks, dfe_labels = contig.dfe_breakpoints()
        for j, intervals in enumerate(contig.interval_list):
            for left, right in intervals:
                assert left in breaks
                assert right in breaks
                k = np.searchsorted(breaks, left)
                assert dfe_labels[k] == j

    def test_default_gc(self):
        contig = stdpopsim.Contig.basic_contig(length=42)
        assert contig.gene_conversion_fraction is None
        assert contig.gene_conversion_length is None
        assert contig.bacterial_recombination is False

    def test_gc_errors(self):
        # should not error
        contig = stdpopsim.Contig.basic_contig(
            length=100, bacterial_recombination=True, gene_conversion_length=14
        )
        assert contig.bacterial_recombination is True
        contig = stdpopsim.Contig.basic_contig(
            length=100, gene_conversion_fraction=1.0, gene_conversion_length=14
        )
        assert contig.gene_conversion_fraction == 1.0
        assert contig.gene_conversion_length == 14
        # can't set gc fraction for bacterial recomb
        with pytest.raises(ValueError, match="Cannot set.*bacterial recomb"):
            _ = stdpopsim.Contig.basic_contig(
                length=100,
                bacterial_recombination=True,
                gene_conversion_fraction=0.4,
                gene_conversion_length=20,
            )
        # must set gc length for bacterial recomb
        with pytest.raises(ValueError, match="Must set.*bacterial recomb"):
            _ = stdpopsim.Contig.basic_contig(length=100, bacterial_recombination=True)
        # must set gc length and fraction together otherwise
        with pytest.raises(ValueError, match="without setting"):
            _ = stdpopsim.Contig.basic_contig(length=100, gene_conversion_fraction=0.4)
        with pytest.raises(ValueError, match="without setting"):
            _ = stdpopsim.Contig.basic_contig(length=100, gene_conversion_length=14)
        # gc_frac must be between 0 and 1 (inclusive)
        for bad_frac in [-0.1, 1.1]:
            with pytest.raises(ValueError, match="must be between 0 and 1"):
                _ = stdpopsim.Contig.basic_contig(
                    length=100,
                    gene_conversion_fraction=bad_frac,
                    gene_conversion_length=14,
                )
        # gc_length must be at least 1
        for bad_len in [-0.1, -10]:
            with pytest.raises(ValueError, match="must be greater than 1"):
                _ = stdpopsim.Contig.basic_contig(
                    length=100,
                    gene_conversion_fraction=0.8,
                    gene_conversion_length=bad_len,
                )

    def test_default_dfe(self):
        contig = stdpopsim.Contig.basic_contig(length=100)
        assert len(contig.dfe_list) == 1
        assert len(contig.interval_list) == 1
        assert contig.dfe_list[0].id == "neutral"
        assert np.all(contig.interval_list[0] == np.array([[0, contig.length]]))

    def test_add_dfe_errors(self):
        contig = stdpopsim.Contig.basic_contig(length=100)
        # bad intervals
        with pytest.raises(ValueError):
            contig.add_dfe(np.array([10, 20]), self.example_dfe)
        with pytest.raises(ValueError):
            contig.add_dfe("abc", self.example_dfe)

    def test_dfe_breakpoints(self):
        contig = stdpopsim.Contig.basic_contig(length=100)
        contig.clear_dfes()
        mt = stdpopsim.MutationType()
        for j in range(3):
            dfe = stdpopsim.DFE(
                id=str(j),
                description="test",
                long_description="test test",
                mutation_types=[mt],
            )
            contig.add_dfe(
                np.array(
                    [[(j + 1) * 5, (j + 1) * 10], [(j + 1) * 20, (j + 1) * 20 + 1]],
                    dtype="int",
                ),
                dfe,
            )
        self.verify_dfe_breakpoints(contig)

    def test_add_dfe(self):
        for clear in (True, False):
            contig = stdpopsim.Contig.basic_contig(length=100)
            if clear:
                contig.clear_dfes()
            props = [0.3, 0.7]
            mt = [stdpopsim.MutationType() for _ in props]
            dfes = [
                stdpopsim.DFE(
                    id=str(j),
                    description="test",
                    long_description="test test",
                    proportions=props,
                    mutation_types=mt,
                )
                for j in range(3)
            ]
            contig.add_dfe(
                intervals=np.array([[10, 30], [50, 100]]),
                DFE=dfes[0],
            )
            contig.add_dfe(intervals=np.array([[30, 40]]), DFE=dfes[1])
            contig.add_dfe(intervals=np.array([[20, 60]]), DFE=dfes[2])
            assert len(contig.dfe_list) == 4 - clear
            assert len(contig.mutation_types()) == 7 - clear
            if clear:
                dfe_ids = []
                true_ints = []
            else:
                true_ints = [np.array([[0, 10]])]
                dfe_ids = ["neutral"]
            dfe_ids += [dfe.id for dfe in dfes]
            true_ints += [
                np.array([[10, 20], [60, 100]]),
                np.empty((0, 2)),
                np.array([[20, 60]]),
            ]
            for d, i in zip(contig.dfe_list, dfe_ids):
                assert d.id == i
            for a1, a2 in zip(contig.interval_list, true_ints):
                assert np.all(a1.shape == a2.shape)
                assert np.all(a1 == a2)
            self.verify_dfe_breakpoints(contig)

    def test_add_dfe_interval_formats(self):
        L = 50818
        for intervals in (
            [[0, L]],
            [[0, int(0.2 * L)]],
            ([0, int(0.2 * L)], [int(0.6 * L), L]),
        ):
            contig = stdpopsim.Contig.basic_contig(length=L)
            contig.add_dfe(intervals=intervals, DFE=self.example_dfe)
            np.testing.assert_array_equal(intervals, contig.interval_list[1])
        with pytest.raises(ValueError, match="must be a numpy array"):
            stdpopsim.utils._check_intervals_validity([[0, 1]])

    def test_is_neutral(self):
        for neutral in (True, False):
            for dist in ("f", "e"):
                contig = stdpopsim.Contig.basic_contig(length=100)
                contig.clear_dfes()
                props = [0.3, 0.7]
                if neutral:
                    s = 0
                else:
                    s = 0.1
                mt = [
                    stdpopsim.MutationType(
                        distribution_type=dist, distribution_args=[s]
                    )
                    for _ in props
                ]
                dfes = [
                    stdpopsim.DFE(
                        id=str(j),
                        description="test",
                        long_description="test test",
                        proportions=props,
                        mutation_types=mt,
                    )
                    for j in range(2)
                ]
                contig.add_dfe(
                    intervals=np.array([[10, 30], [50, 100]]),
                    DFE=dfes[0],
                )
                contig.add_dfe(intervals=np.array([[30, 40]]), DFE=dfes[1])
                # exponential with mean zero doesn't count as neutral!
                assert contig.is_neutral is (neutral and dist == "f")

    def test_chromosome_segment(self):
        species = stdpopsim.get_species("AnaPla")
        chrom = species.genome.chromosomes[1]
        length = chrom.length
        for interval in [(0, 1341), (5020, 12850), (4249, length), (0, length)]:
            contig = species.get_contig(chrom.id, left=interval[0], right=interval[1])
            chromosome, left, right = contig.coordinates
            assert chromosome == chrom.id
            assert left == interval[0]
            assert right == interval[1]
            assert contig.origin == f"{chromosome}:{left}-{right}"
            assert contig.length == length

    def test_not_simulated_outside_region(self):
        # test that when left, right are specified
        # we legit don't simulate anything outside that region
        species = stdpopsim.get_species("AraTha")
        model = stdpopsim.PiecewiseConstantSize(100)
        samples = {"pop_0": 50}

        left, right = 100000, 900000
        contig = species.get_contig("1", left=left, right=right)

        engine = stdpopsim.get_engine("msprime")
        ts = engine.simulate(
            model,
            contig,
            samples,
            seed=236,
        )

        assert ts.sequence_length > right
        assert ts.num_sites > 0
        assert left in ts.breakpoints()
        assert right in ts.breakpoints()
        assert left <= min(ts.sites_position)
        assert max(ts.sites_position) < right
        assert ts.num_trees > 2
        for t in ts.trees(root_threshold=2):
            tl, tr = t.interval
            if tl > right or tr < left:
                assert t.num_roots == 0

    def test_chromosome_segment_with_genetic_map(self):
        chr_id = "chr2"
        left = 1000001
        right = 3000000
        species = stdpopsim.get_species("HomSap")
        genmap = "HapMapII_GRCh38"
        contig = species.get_contig(chr_id, left=left, right=right, genetic_map=genmap)
        chrom = species.get_contig(chr_id, genetic_map=genmap)
        assert contig.recombination_map.get_rate(
            left
        ) == chrom.recombination_map.get_rate(left)
        assert contig.recombination_map.get_rate(
            right - 1
        ) == chrom.recombination_map.get_rate(right - 1)

    def test_chromosome_segment_with_mask(self):
        chr_id = "chr2"
        left = 1000001
        right = 3000000
        length = right - left
        offset = int(length * 0.1)
        mask = [
            (left - 2 * offset, left - offset),
            (left - offset, left + offset),
            (left + offset, left + 2 * offset),
            (right - offset, right + offset),
            (right + offset, right + 2 * offset),
        ]
        clipped_mask = [
            (left, left + offset),
            (left + offset, left + 2 * offset),
            (right - offset, right),
        ]
        species = stdpopsim.get_species("HomSap")
        contig_inclusion = species.get_contig(
            chromosome=chr_id,
            left=left,
            right=right,
            inclusion_mask=mask,
        )
        contig_exclusion = species.get_contig(
            chromosome=chr_id,
            left=left,
            right=right,
            exclusion_mask=mask,
        )
        np.testing.assert_array_equal(
            contig_inclusion.inclusion_mask,
            clipped_mask,
        )
        np.testing.assert_array_equal(
            contig_exclusion.exclusion_mask,
            clipped_mask,
        )

    def test_generic_contig_origin(self):
        length = 42
        contig = stdpopsim.Contig.basic_contig(length=length)
        chromosome, left, right = contig.coordinates
        assert chromosome is None
        assert left == 0
        assert right == length
        assert contig.origin is None

    def test_add_single_site_coordinate_system(self):
        chrom = "chr2"
        species = stdpopsim.get_species("HomSap")
        interval = [100000, 200000]
        original_sweep_coord = 100100
        bad_sweep_coord = original_sweep_coord - interval[0]
        contig = species.get_contig(chrom, left=interval[0], right=interval[1])
        # Correctly specified coordinate system
        contig.add_single_site(
            id="sweep",
            coordinate=original_sweep_coord,
        )
        # Coordinate system mismatch. Shifted coordinate falls outside of
        # `interval` so the added DFE has empty intervals: throw a warning
        with pytest.warns(UserWarning, match="No intervals remain"):
            contig.add_single_site(
                id="sweep",
                coordinate=bad_sweep_coord,
            )

    def test_original_coordinates(self):
        species = stdpopsim.get_species("AnaPla")
        contig = species.get_contig("2", left=10000, right=392342)
        c, x, y = contig.coordinates
        with pytest.warns(
            stdpopsim.DeprecatedFeatureWarning, match="no longer shifted"
        ):
            oc, ox, oy = contig.original_coordinates
        assert c == oc
        assert x == ox
        assert y == oy

    def test_add_dfe_coordinate_system(self):
        chrom = "chr2"
        species = stdpopsim.get_species("HomSap")
        contig_interval = [200000, 300000]
        original_dfe_interval = np.array([[190000, 210000], [290000, 310000]])
        bad_dfe_interval = np.array([[0, 10000], [90000, 100000]])
        contig = species.get_contig(
            chrom, left=contig_interval[0], right=contig_interval[1]
        )
        # Correctly specified coordinate system
        contig.add_dfe(
            intervals=original_dfe_interval,
            DFE=self.example_dfe,
        )
        # Coordinate system mismatch. Shifted DFE intervals all fall outside of
        # `contig.origin` so the added DFE is empty: throw a warning
        with pytest.warns(UserWarning, match="No intervals remain"):
            contig.add_dfe(
                intervals=bad_dfe_interval,
                DFE=self.example_dfe,
            )

    def test_dfe_breakpoints_coordinate_system(self):
        chrom = "chr2"
        species = stdpopsim.get_species("HomSap")
        contig_interval = [200000, 300000]
        dfe_interval = [[190000, 210000], [290000, 310000]]
        contig = species.get_contig(
            chrom, left=contig_interval[0], right=contig_interval[1]
        )
        contig.add_dfe(
            intervals=dfe_interval,
            DFE=self.example_dfe,
        )
        breakpoints, _ = contig.dfe_breakpoints()
        for x, y in zip(breakpoints, [0, 200000, 210000, 290000, 300000]):
            assert x == y

        with pytest.warns(
            stdpopsim.DeprecatedFeatureWarning, match="relative_coordinates argument"
        ):
            ob, _ = contig.dfe_breakpoints(relative_coordinates=True)
        for x, y in zip(breakpoints, ob):
            assert x == y

    @pytest.mark.filterwarnings("ignore::stdpopsim.DeprecatedFeatureWarning")
    def test_chromosome_segment_fails_with_length_multiplier(self):
        chrom = "chr2"
        species = stdpopsim.get_species("HomSap")
        with pytest.raises(ValueError, match="specifying left or right"):
            species.get_contig(
                chrom,
                left=200000,
                right=300000,
                length_multiplier=8.2,
            )

    def test_no_chromosome_segment_with_generic_contig(self):
        species = stdpopsim.get_species("HomSap")
        interval = [200000, 300000]
        with pytest.raises(ValueError, match="coordinates with generic contig"):
            species.get_contig(
                left=interval[0],
                right=interval[1],
            )

    def test_chromosome_segment_exceeds_length(self):
        chr_id = "chr22"
        species = stdpopsim.get_species("HomSap")
        for interval in [[50e6, 80e6], [80e6, 100e6], [-1, 10e6]]:
            with pytest.raises(ValueError, match="the length of"):
                species.get_contig(
                    chromosome=chr_id,
                    left=interval[0],
                    right=interval[1],
                )

    def test_basic_contig_is_diploid(self):
        contig = stdpopsim.Contig.basic_contig(length=1000)
        assert contig.ploidy == 2


class TestGeneConversion(object):
    def test_mean_gene_conversion(self):
        dro_mel = stdpopsim.get_species("DroMel")
        contig = dro_mel.get_contig(length=1000, use_species_gene_conversion=False)
        assert contig.gene_conversion_fraction is None
        assert contig.gene_conversion_length is None
        contig = dro_mel.get_contig(length=1000, use_species_gene_conversion=True)
        mean_gc = dro_mel.genome.mean_gene_conversion_fraction
        assert np.isclose(mean_gc, contig.gene_conversion_fraction)

    def test_no_mean_gene_conversion(self):
        # this will need to be changed if we add GC rates to AnaPla
        ana_pla = stdpopsim.get_species("AnaPla")
        contig = ana_pla.get_contig(length=1000, use_species_gene_conversion=True)
        mean_gc = ana_pla.genome.mean_gene_conversion_fraction
        assert mean_gc == 0.0
        assert contig.gene_conversion_fraction is None

    def test_modifiying_gene_conversion(self):
        esc_col = stdpopsim.get_species("EscCol")
        assert esc_col.genome.mean_gene_conversion_fraction == 0
        contig = esc_col.get_contig(chromosome="Chromosome")
        assert contig.bacterial_recombination is True
        assert contig.gene_conversion_fraction is None
        assert contig.gene_conversion_length > 0

        rm = msprime.RateMap(position=[0.0, contig.length], rate=[0.1])
        contig.recombination_map = rm
        contig.gene_conversion_length = 1
        assert contig.recombination_map == rm
        assert contig.gene_conversion_length == 1

    def test_use_species_gene_conversion_errors(self):
        esc_col = stdpopsim.get_species("EscCol")
        with pytest.raises(ValueError, match="with bacterial recombination"):
            _ = esc_col.get_contig(
                chromosome="Chromosome", use_species_gene_conversion=True
            )

    def test_use_species_gene_conversion_unnamed_contig_undefined_gc(self):
        species = stdpopsim.get_species("AnaPla")
        contig = species.get_contig("1")
        assert contig.gene_conversion_fraction is None
        assert contig.gene_conversion_length is None
        assert contig.bacterial_recombination is False
        contig = species.get_contig(length=42, use_species_gene_conversion=True)
        assert contig.gene_conversion_fraction is None
        assert contig.gene_conversion_length is None
        assert contig.bacterial_recombination is False

    def test_get_species_contig_rates(self):
        species = stdpopsim.get_species("AnaPla")
        chrom = species.genome.chromosomes[0]
        gc_frac = 0.7
        gc_len = 123
        # this is *not* the recommended workflow; just doing it for testing
        # see test_modifiying_gene_conversion for the recommended way to modify it
        chrom.gene_conversion_fraction = gc_frac
        chrom.gene_conversion_length = gc_len
        assert species.genome.chromosomes[0].gene_conversion_fraction == gc_frac
        assert species.genome.chromosomes[0].gene_conversion_length == gc_len
        # without use_species_gene_conversion
        contig = species.get_contig(chrom.id)
        assert contig.gene_conversion_fraction is None
        assert contig.gene_conversion_length is None
        assert contig.recombination_map.mean_rate == chrom.recombination_rate
        # now, with use_species_gene_conversion
        contig = species.get_contig(chrom.id, use_species_gene_conversion=True)
        assert contig.gene_conversion_fraction == gc_frac
        assert contig.gene_conversion_length == gc_len
        assert np.isclose(
            contig.recombination_map.mean_rate, chrom.recombination_rate / (1 - gc_frac)
        )
        species = stdpopsim.get_species("HomSap")
        with pytest.raises(
            ValueError, match="Cannot use recombination rate with genetic map"
        ):
            _ = species.get_contig(
                "chr22", genetic_map="HapMapII_GRCh38", recombination_rate=1e-8
            )
