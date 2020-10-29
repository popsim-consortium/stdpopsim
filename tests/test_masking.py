"""
Tests for the genetic maps management.
"""
import unittest

import stdpopsim.utils
import msprime
import numpy as np
import os
import sys

IS_WINDOWS = sys.platform.startswith("win")


class TestMasking(unittest.TestCase):
    def test_load_intervals(self):
        intervals_in = {
            "chr1": [(10, 10000)],
            "chr22": [(100, 1000), (2000, 5000), (6000, 9000)],
        }
        with open("temp_file.bed", "w+") as fout:
            for c, i in intervals_in.items():
                for (l, r) in i:
                    fout.write(f"{c}\t{l}\t{r}\n")
        intervals_chr1 = stdpopsim.utils.read_bed("temp_file.bed", "chr1")
        intervals_chr22 = stdpopsim.utils.read_bed("temp_file.bed", "chr22")
        os.remove("temp_file.bed")
        for interval in intervals_chr1:
            self.assertTrue(tuple(interval) in intervals_in["chr1"])
        for interval in intervals_chr22:
            self.assertTrue(tuple(interval) in intervals_in["chr22"])
        self.assertTrue(len(intervals_chr1) == len(intervals_in["chr1"]))
        self.assertTrue(len(intervals_chr22) == len(intervals_in["chr22"]))

    def test_length_interval_invalid(self):
        species = stdpopsim.get_species("HomSap")
        with self.assertRaises(ValueError):
            contig = species.get_contig(
                "chr22", length_multiplier=0.1, inclusion_mask="test.bed"
            )
            print(contig.recombination_map.get_sequence_length())

    def test_read_bed(self):
        def intervals_to_keep_test_func(mask_fpath, chrom):
            intervals = []
            with open(mask_fpath, "r") as fin:
                for line in fin:
                    if line.split()[0] == chrom:
                        intervals.append((int(line.split()[1]), int(line.split()[2])))
            intervals = sorted(intervals)
            return intervals

        intervals_in = {"chr1": [(0, 10), (100, 1000), (2000, 5000), (6000, 10000)]}
        with open("temp_file.bed", "w+") as fout:
            for c, i in intervals_in.items():
                for (l, r) in i:
                    fout.write(f"{c}\t{l}\t{r}\n")
        intervals_test = intervals_to_keep_test_func("temp_file.bed", "chr1")
        intervals_utils = stdpopsim.utils.read_bed("temp_file.bed", "chr1")
        os.remove("temp_file.bed")
        for i1, i2 in zip(intervals_utils, intervals_test):
            self.assertTrue(i1[0] == i2[0])
            self.assertTrue(i1[1] == i2[1])

    def test_mask_from_intervals(self):
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr22", exclusion_mask=[(0, 16.5e6)])
        self.assertTrue(contig.inclusion_mask is None)
        self.assertTrue(contig.exclusion_mask[0][0] == 0)
        self.assertTrue(contig.exclusion_mask[0][1] == 16.5e6)
        contig = species.get_contig("chr22", inclusion_mask=[(16.5e6, 45e6)])
        self.assertTrue(contig.inclusion_mask[0][0] == 16.5e6)
        self.assertTrue(contig.inclusion_mask[0][1] == 45e6)
        self.assertTrue(contig.exclusion_mask is None)

    def test_multiple_masks(self):
        species = stdpopsim.get_species("HomSap")
        with self.assertRaises(ValueError):
            species.get_contig(
                "chr22", exclusion_mask=[(0, 16.5e6)], inclusion_mask=[(16.5e6, 50e6)]
            )

    def test_read_masks_from_bed(self):
        intervals_in = {"chr1": [(0, 10), (100, 1000), (2000, 5000), (6000, 10000)]}
        with open("temp_file.bed", "w+") as fout:
            for c, i in intervals_in.items():
                for (l, r) in i:
                    fout.write(f"{c}\t{l}\t{r}\n")
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr1", inclusion_mask="temp_file.bed")
        self.assertTrue(contig.exclusion_mask is None)
        self.assertTrue(len(contig.inclusion_mask) == 4)
        contig = species.get_contig("chr1", exclusion_mask="temp_file.bed")
        self.assertTrue(contig.inclusion_mask is None)
        self.assertTrue(contig.exclusion_mask[1][0] == 100)
        os.remove("temp_file.bed")

    def test_mask_tree_sequence(self):
        intervals = np.array([[0, 10], [100, 200], [500, 1000]])
        ts = msprime.simulate(Ne=1000, samples=[(0, 0), (0, 0)], length=1000)
        ts = msprime.mutate(ts, rate=1e-4)

        exclude = False
        ts_mask = stdpopsim.utils.mask_tree_sequence(ts, intervals, exclude)
        # check we get the right number of trees (simulated with no recomb)
        self.assertTrue(ts_mask.num_trees == len(intervals) * 2 - 1)

        exclude = True
        intervals = np.array([[0, 100], [200, 1000]])
        ts_mask = stdpopsim.utils.mask_tree_sequence(ts, intervals, exclude)
        # check we get the right number of trees (simulated with no recomb)
        self.assertTrue(ts_mask.num_trees == len(intervals) * 2 - 1)
        # check that positions of mutations in ts_mask are within mask
        positions = np.array([m.position for m in ts_mask.mutations()])
        self.assertTrue(np.all(np.logical_and(100 <= positions, positions < 200)))


@unittest.skipIf(IS_WINDOWS, "SLiM not available on windows")
class TestSimulate(unittest.TestCase):
    def test_simulate_with_mask(self):
        engines = ["msprime", "slim"]
        species = stdpopsim.get_species("HomSap")
        L = 1000
        contig = species.get_contig(length=L)
        contig.mutation_rate = 1e-3
        contig.recombination_map = msprime.RecombinationMap.uniform_map(L, 0)
        samples = [msprime.Sample(0, 0), msprime.Sample(0, 0)]
        model = stdpopsim.PiecewiseConstantSize(100)
        for engine_name in engines:
            engine = stdpopsim.get_engine(engine_name)

            # test engine with exclusion mask
            contig.inclusion_mask = None
            contig.exclusion_mask = np.array([[0, L // 2]])
            ts = engine.simulate(
                demographic_model=model, contig=contig, samples=samples
            )
            # check that positions of mutations are within mask
            positions = np.array([m.position for m in ts.mutations()])
            self.assertTrue(np.all(positions >= L // 2))

            # test engine with inclusion mask
            contig.exclusion_mask = None
            contig.inclusion_mask = np.array([[0, L // 2]])
            ts = engine.simulate(
                demographic_model=model, contig=contig, samples=samples
            )
            # check that positions of mutations are within mask
            positions = np.array([m.position for m in ts.mutations()])
            self.assertTrue(np.all(positions < L // 2))
