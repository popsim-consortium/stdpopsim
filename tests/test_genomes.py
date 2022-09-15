"""
Tests for classes that hold information about genomic region to be simulated.
"""
import numpy as np
import pytest
import copy
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
        assert contig.gene_conversion_rate is None
        assert contig.gene_conversion_length is None

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
        chrom = "chr2"
        interval = [100, 130]
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig(chrom, left=interval[0], right=interval[1])
        chromosome, left, right = contig.original_coordinates
        assert chromosome == chrom
        assert left == interval[0]
        assert right == interval[1]
        assert contig.length == interval[1] - interval[0]
        assert contig.origin == f"{chromosome}:{left}-{right}"

    def test_chromosome_segment_with_genetic_map(self):
        chr_id = "chr2"
        left = 1000001
        right = 3000000
        species = stdpopsim.get_species("HomSap")
        genmap = "HapMapII_GRCh38"
        contig = species.get_contig(chr_id, left=left, right=right, genetic_map=genmap)
        chrom = species.get_contig(chr_id, genetic_map=genmap)
        assert contig.recombination_map.get_rate(0) == chrom.recombination_map.get_rate(
            left
        )
        assert contig.recombination_map.get_rate(
            contig.length - 1
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
        clipped_mask = [(0, offset), (offset, 2 * offset), (length - offset, length)]
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
        chromosome, left, right = contig.original_coordinates
        assert chromosome is None
        assert left == 0
        assert right == length
        assert contig.origin is None

    def test_add_single_site_coordinate_system(self):
        chrom = "chr2"
        species = stdpopsim.get_species("HomSap")
        interval = [100000, 200000]
        original_sweep_coord = 100100
        shifted_sweep_coord = original_sweep_coord - interval[0]
        contig = species.get_contig(chrom, left=interval[0], right=interval[1])
        # Correctly specified coordinate system
        contig.add_single_site(
            id="sweep",
            coordinate=original_sweep_coord,
        )
        assert contig.interval_list[1][0, 0] == shifted_sweep_coord
        # Coordinate system mismatch. Shifted coordinate falls outside of
        # `interval` so the added DFE has empty intervals: throw a warning
        with pytest.warns(UserWarning, match="No intervals remain"):
            contig.add_single_site(
                id="sweep",
                coordinate=shifted_sweep_coord,
            )

    def test_add_dfe_coordinate_system(self):
        chrom = "chr2"
        species = stdpopsim.get_species("HomSap")
        contig_interval = [200000, 300000]
        original_dfe_interval = [[190000, 210000], [290000, 310000]]
        shifted_dfe_interval = [[0, 10000], [90000, 100000]]
        contig = species.get_contig(
            chrom, left=contig_interval[0], right=contig_interval[1]
        )
        # Correctly specified coordinate system
        contig.add_dfe(
            intervals=original_dfe_interval,
            DFE=self.example_dfe,
        )
        np.testing.assert_array_equal(contig.interval_list[1], shifted_dfe_interval)
        # Coordinate system mismatch. Shifted DFE intervals all fall outside of
        # `contig.origin` so the added DFE is empty: throw a warning
        with pytest.warns(UserWarning, match="No intervals remain"):
            contig.add_dfe(
                intervals=shifted_dfe_interval,
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
        shifted_breakpoints, _ = contig.dfe_breakpoints(relative_coordinates=True)
        original_breakpoints, _ = contig.dfe_breakpoints(relative_coordinates=False)
        for x, y in zip(shifted_breakpoints, original_breakpoints):
            assert x + contig_interval[0] == y

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
        for interval in [[50e6, 80e6], [80e6, 100e6]]:
            with pytest.raises(ValueError, match="falls outside chromosome"):
                species.get_contig(
                    chromosome=chr_id,
                    left=interval[0],
                    right=interval[1],
                )


class TestGeneConversion(object):
    species = stdpopsim.get_species("EscCol")

    def test_modifiying_gene_conversion(self):
        genome = self.species.get_contig(
            chromosome="Chromosome", gene_conversion_rate=0.1, gene_conversion_length=1
        )
        assert genome.gene_conversion_rate == 0.1
        assert genome.gene_conversion_length == 1

    def test_use_species_gene_conversion(self):
        genome = self.species.get_contig(
            chromosome="Chromosome", use_species_gene_conversion=True
        )
        assert genome.gene_conversion_rate == 8.9e-11
        assert genome.gene_conversion_length == 345

    def test_unnamed_contig_with_gc(self):
        genome = self.species.get_contig(
            length=42, gene_conversion_rate=0.1, gene_conversion_length=1
        )
        assert genome.gene_conversion_rate == 0.1
        assert genome.gene_conversion_length == 1

    def test_use_species_gene_conversion_unnamed_contig(self):
        genome = self.species.get_contig(length=42, use_species_gene_conversion=True)
        assert genome.gene_conversion_rate == 8.9e-11
        assert genome.gene_conversion_length == 345

    def test_use_species_gene_conversion_unnamed_contig_undefined_gc(self):
        mod_species = copy.deepcopy(self.species)
        mod_species.genome.get_chromosome("Chromosome").gene_conversion_rate = None
        mod_species.genome.get_chromosome("Chromosome").gene_conversion_length = None
        genome = mod_species.get_contig(length=42, use_species_gene_conversion=True)
        assert genome.gene_conversion_rate is None
        assert genome.gene_conversion_length is None
