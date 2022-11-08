"""
Tests for the utils module.
"""
import functools
import os
import pathlib
import tarfile
import tempfile

import numpy as np
import scipy.stats
import pytest

from stdpopsim import utils
from stdpopsim import Chromosome, Genome


class TestValidDemographicModelId:
    """
    Tests for the is_valid_demographic_model_id function.
    """

    def test_empty_string(self):
        assert not (utils.is_valid_demographic_model_id(""))

    def test_contains_spaces(self):
        bad_ids = ["CamelCase 4X19", " CamelCase4X19", "CamelCase4X19"]
        for bad_id in bad_ids:
            assert not (utils.is_valid_demographic_model_id(bad_id))

    def test_too_many_sections(self):
        bad_ids = [
            "_CamelCase_4X19",
            "CamelCase_4X19_",
            "CamelCase__4X19",
            "CamelCase_4X19_1234_5678",
        ]
        for bad_id in bad_ids:
            assert not (utils.is_valid_demographic_model_id(bad_id))

    def test_bad_name(self):
        bad_ids = ["camelCase_4X19", "1CamelCase_4X19", "Camel-Case_4X19"]
        for bad_id in bad_ids:
            assert not (utils.is_valid_demographic_model_id(bad_id))

    def test_bad_populations(self):
        bad_ids = ["CamelCase_0X19", "CamelCase_0.1X19", "CamelCase_OneX19"]
        for bad_id in bad_ids:
            assert not (utils.is_valid_demographic_model_id(bad_id))

    def test_bad_author(self):
        bad_ids = ["CamelCase_1019", "CamelCase_1-19", "CamelCase_1;19"]
        for bad_id in bad_ids:
            assert not (utils.is_valid_demographic_model_id(bad_id))

    def test_bad_year(self):
        bad_ids = [
            "CamelCase_10X1",
            "CamelCase_1X",
            "CamelCase_1X2001",
            "CamelCase_1Xfive",
        ]
        for bad_id in bad_ids:
            assert not (utils.is_valid_demographic_model_id(bad_id))

    def test_good_ids(self):
        good_ids = [
            "CamelCase_10X10",
            "CamelCase_1Y00",
            "CamelCase_100A05",
            "Camel_1X10",
            "C_1Y00",
            "C01234_1Y00",
        ]
        for good_id in good_ids:
            assert utils.is_valid_demographic_model_id(good_id)


class TestValidGeneticMapId:
    """
    Tests for the is_valid_demographic_model_id function.
    """

    def test_empty_string(self):
        assert not (utils.is_valid_genetic_map_id(""))

    def test_contains_spaces(self):
        bad_ids = ["CamelCase ABCD", " CamelCase_ABCD", "CamelCase_ABCD "]
        for bad_id in bad_ids:
            assert not (utils.is_valid_genetic_map_id(bad_id))

    def test_bad_name(self):
        bad_ids = ["camelCase_ABCD", "1CamelCase_ABCD", "Camel-Case_ABCD"]
        for bad_id in bad_ids:
            assert not (utils.is_valid_genetic_map_id(bad_id))

    def test_bad_assembly(self):
        bad_ids = [
            "CamelCase_AB CD",
            "CamelCase_",
            "CamelCase_AB-CD",
            "CamelCase_AB.CD",
        ]
        for bad_id in bad_ids:
            assert not (utils.is_valid_genetic_map_id(bad_id))

    def test_good_ids(self):
        good_ids = ["CamelCase_ABCD", "CamelCase_a", "CamelCase_1000", "C_ABCD", "C_X"]
        for good_id in good_ids:
            assert utils.is_valid_genetic_map_id(good_id)


class TestValidSpeciesId:
    """
    Tests for the is_valid_demographic_model_id function.
    """

    def test_empty_string(self):
        assert not (utils.is_valid_species_id(""))

    def test_contains_spaces(self):
        bad_ids = ["Cam Cas", " CamCas", "CamCas "]
        for bad_id in bad_ids:
            assert not (utils.is_valid_species_id(bad_id))

    def test_bad_ids(self):
        bad_ids = ["camCas", "camcas", "CAMCAS", "C1mCas", "Camcas", "CamelCase"]
        for bad_id in bad_ids:
            assert not (utils.is_valid_species_id(bad_id))

    def test_good_ids(self):
        bad_ids = ["CamCas", "HomSap"]
        for bad_id in bad_ids:
            assert utils.is_valid_species_id(bad_id)


class TestValidSpeciesName:
    """
    Tests for the is_valid_demographic_model_id function.
    """

    def test_empty_string(self):
        assert not (utils.is_valid_species_name(""))

    def test_extra_spaces(self):
        bad_names = ["Pan  pan", " Pan pan", "Pan pan "]
        for bad_name in bad_names:
            assert not (utils.is_valid_species_name(bad_name))

    def test_bad_names(self):
        bad_names = [
            "PAN PAN",
            "Pan Pan",
            "Panpan",
            "PanPan",
            "pan pan",
            "Pan p0n",
            "Pan, pan",
        ]
        for bad_name in bad_names:
            assert not (utils.is_valid_species_name(bad_name))

    def test_good_names(self):
        good_names = ["Pan pan", "Homo sapiens"]
        for good_name in good_names:
            assert utils.is_valid_species_name(good_name)

    @pytest.mark.skip("Implement more flexible species name")
    def test_three_or_more(self):
        bad_names = ["Pan pan pan"]
        for bad_name in bad_names:
            assert utils.is_valid_species_name(bad_name)


class TestValidSpeciesCommonName:
    """
    Tests for the is_valid_demographic_model_id function.
    """

    def test_empty_string(self):
        assert not (utils.is_valid_species_common_name(""))

    def test_bad_common_names(self):
        # TODO add more examples when we define restrictions (#330)
        bad_common_names = ["0", "human"]
        for bad_common_name in bad_common_names:
            assert not (utils.is_valid_species_common_name(bad_common_name))

    def test_good_common_names(self):
        # TODO add more examples when we define restrictions (#330)
        good_names = ["Human", "Stuff and things"]
        for good_name in good_names:
            assert utils.is_valid_species_common_name(good_name)


class TestDownload:
    def test_download_local_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = pathlib.Path(tmpdir)
            input_filename = tmpdir / "input.foo"
            output_filename = tmpdir / "output.foo"
            with open(input_filename, "w") as f:
                print("foofoofoo", file=f)
            url = input_filename.resolve().as_uri()
            utils.download(url, output_filename)
            assert output_filename.exists()
            with open(output_filename) as f:
                contents = f.read()
            assert contents.strip() == "foofoofoo"

    def test_download_nonexistent(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = pathlib.Path(tmpdir)
            input_filename = tmpdir / "nonexistant"
            output_filename = tmpdir / "output.foo"
            for url in (
                input_filename.resolve().as_uri(),  # local file
                "http://example.com/nonexistant",  # remote file
            ):
                with pytest.raises(OSError, match="nonexistant"):
                    utils.download(url, output_filename)
                assert not output_filename.exists()


class TestSha256:
    # These tests are not intended to comprehensively test sha256
    # calculations, as we assume the Python standard library tests
    # have that covered. Instead, these tests are to ensure that
    # we've used the standard library sha256 functionality correctly.
    def test_sha256_small(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            filename = pathlib.Path(tmpdir) / "test.foo"
            # write binary file, to avoid windows/mac/unix line ending differences
            with open(filename, "wb") as f:
                f.write(b"foo-bar-baz")
            sha256sum = utils.sha256(filename)
            assert (
                sha256sum
                ==
                # Calculated with `sha256sum` from GNU coreutils.
                "269dce1a5bb90188b2d9cf542a7c30e410c7d8251e34a97bfea56062df51ae23"
            )

    def test_sha256_big(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            filename = pathlib.Path(tmpdir) / "test.foo"
            # write binary file, to avoid windows/mac/unix line ending differences
            with open(filename, "wb") as f:
                for i in range(1024 * 1024):
                    f.write(str(i).encode())
            sha256sum = utils.sha256(filename)
            assert (
                sha256sum
                ==
                # Calculated with `sha256sum` from GNU coreutils.
                "995e0fde646f7dc98423af9a862be96014574bfa76be1186b484f796c4e58533"
            )


class TestCd:
    def test_cd_context_manager(self):
        # On Mac, the path we enter with "cd" may differ from
        # the path we get with cwd() due to symlinks. So we
        # resolve all paths here to ignore symlink-only differences.
        old_cwd = pathlib.Path.cwd().resolve()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = pathlib.Path(tmpdir).resolve()
            assert old_cwd != tmpdir
            with utils.cd(tmpdir):
                assert pathlib.Path.cwd().resolve() == tmpdir
            assert pathlib.Path.cwd().resolve() == old_cwd


def rm_f(filename):
    try:
        os.unlink(filename)
    except OSError:
        pass


class TestUntar:
    def test_untar(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = pathlib.Path(tmpdir)
            test_files = ["foo", "bar", "baz"]
            sha_list = []
            for name in test_files:
                filename = tmpdir / name
                with open(filename, "wb") as f:
                    f.write(name.encode())
                # Record checksums for later.
                sha_list.append(utils.sha256(filename))

            tar = tmpdir / "test.tgz"
            with utils.cd(tmpdir):
                with tarfile.open(tar, mode="w:gz") as tf:
                    for name in test_files:
                        tf.add(name)

            dest = tmpdir / "dest"
            dest.mkdir()
            utils.untar(tar, dest)
            for name, sha in zip(test_files, sha_list):
                filename = dest / name
                assert filename.exists()
                # Check that extracted files have the same checksums as
                # the files we put in the tar.
                assert utils.sha256(filename) == sha

    # Security related tests for utils.untar().
    # For an extended discussion of tar-related security issues see:
    # https://bugs.python.org/issue21109

    def test_tar_path_traversal_attack(self):
        # Test for vulnerability to path-traversal attacks.
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = pathlib.Path(tmpdir)
            dest = tmpdir / "dest"
            dest.mkdir()
            for link_dest in ("../nonexistant", "/nonexistant"):
                tar = tmpdir / "symlink-path-traversal.tgz"
                filename = tmpdir / "link"
                filename.symlink_to(link_dest)
                with utils.cd(tmpdir):
                    with tarfile.open(tar, mode="w:gz") as tf:
                        tf.add("link")
                with pytest.raises(ValueError):
                    utils.untar(tar, dest)
                rm_f(filename)
                rm_f(tar)

            for name in ("../nonexistant", "/nonexistant"):
                tar = tmpdir / "path-traversal.tgz"
                filename = tmpdir / "test-thing"
                with open(filename, "w") as f:
                    print("blah", file=f)
                with utils.cd(tmpdir):
                    with tarfile.open(tar, mode="w:gz") as tf:

                        def filt(info):
                            info.name = name  # path the file will be extracted to
                            return info

                        tf.add("test-thing", filter=filt)
                with pytest.raises(ValueError):
                    utils.untar(tar, dest)
                rm_f(filename)
                rm_f(tar)

    def test_bad_tar_members(self):
        # Pretend we downloaded a tarball containing a FIFO or device file.
        # There is no reasonable use for these types of files in stdpopsim,
        # so their presence likely indicates a maliciously crafted tarball.
        # Creating a character or block special device file requires root
        # privileges, so we instead modify the ``type`` field of each file
        # in the tar.
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = pathlib.Path(tmpdir)
            dest = tmpdir / "dest"
            dest.mkdir()
            for filename, type_ in [
                ("fifo", tarfile.FIFOTYPE),
                ("char-device", tarfile.CHRTYPE),
                ("block-device", tarfile.BLKTYPE),
            ]:
                tar = tmpdir / "irregular-type.tgz"
                filename = tmpdir / "irregular"
                with open(filename, "w") as f:
                    print("blah", file=f)
                with utils.cd(tmpdir):
                    with tarfile.open(tar, mode="w:gz") as tf:

                        def filt(info):
                            info.type = type_  # lie about the type
                            return info

                        tf.add("irregular", filter=filt)
                with pytest.raises(ValueError):
                    utils.untar(tar, dest)
                rm_f(filename)
                rm_f(tar)


class TestSynonyms:
    def test_add_digit_autosomes(self):
        chroms = [
            Chromosome(
                id="1",
                length=100,
                synonyms=[],
                mutation_rate=1e-8,
                recombination_rate=1e-8,
            )
        ]
        genome = Genome(chroms)
        utils.append_common_synonyms(genome)
        assert "chr1" in genome.chromosomes[0].synonyms

    def test_add_drosophila_like_synonyms(self):
        chroms = [
            Chromosome(
                id="1a",
                length=100,
                synonyms=[],
                mutation_rate=1e-8,
                recombination_rate=1e-8,
            )
        ]
        genome = Genome(chroms)
        utils.append_common_synonyms(genome)
        assert "chr1a" in genome.chromosomes[0].synonyms

    def test_add_sex_chrom_synonyms(self):
        chroms = [
            Chromosome(
                id=s,
                length=100,
                synonyms=[],
                mutation_rate=1e-8,
                recombination_rate=1e-8,
            )
            for s in ["X", "Y", "Z", "W"]
        ]
        genome = Genome(chroms)
        utils.append_common_synonyms(genome)
        for chrom in genome.chromosomes:
            assert "chr" + chrom.id in chrom.synonyms

    def test_add_only_unique_synonyms(self):
        chroms = [
            Chromosome(
                id="1",
                length=100,
                synonyms=["chr1"],
                mutation_rate=1e-8,
                recombination_rate=1e-8,
            )
        ]
        genome = Genome(chroms)
        utils.append_common_synonyms(genome)
        assert "chr1" in genome.chromosomes[0].synonyms
        assert len(genome.chromosomes[0].synonyms) == 1


class TestIntervalUtilities:
    def test_intervals_array_shape(self):
        # shape (n, x) where x >= 2
        bad_arrays = (np.array([10, 20]), np.array([[[10, 20]]]))
        for bad_array in bad_arrays:
            with pytest.raises(ValueError):
                utils._check_intervals_array_shape(intervals=bad_array)

    def test_invalid_intervals(self):
        # right <= left, overlapping, uncastable float->int
        invalid_arrays = (
            np.array([[20, 10]]),
            np.array([[10, 20], [15, 20]]),
            np.array([[10, 20], [13, 41]]),
        )
        for invalid_array in invalid_arrays:
            with pytest.raises(ValueError):
                utils._check_intervals_validity(intervals=invalid_array)
        # interval outside [start, end)
        with pytest.raises(ValueError):
            utils._check_intervals_validity(intervals=np.array([[10, 50]]), end=30)
        with pytest.raises(ValueError):
            utils._check_intervals_validity(intervals=np.array([[10.2, 50.3]]))

    def intervals_to_set(self, intervals):
        # note: we may want to use intervaltree for more complicated operations
        # in the future
        utils._check_intervals_validity(intervals)
        return functools.reduce(
            lambda a, b: a | b, [set(range(a, b)) for a, b in intervals], set()
        )

    def test_mask_intervals(self):
        empty_array = np.array([]).reshape((0, 2))
        for a, b in [
            (empty_array, [[1, 2]]),
            (
                [[0, 5], [5, 10]],
                [[0, 10]],
            ),
            (
                [[3, 4], [5, 12]],
                [[0, 10]],
            ),
            (
                [[0, 3], [3, 4], [5, 9]],
                [[1, 10]],
            ),
            (
                [[k, k + 1] for k in range(12)],
                [[0, 10]],
            ),
            (
                [[1, 4], [6, 9], [11, 13], [14, 17], [19, 21]],
                [[0, 3], [5, 9], [11, 13], [15, 17], [20, 22]],
            ),
            (
                [[1, 2], [6, 7], [11, 12], [14, 15], [19, 20], [20, 21], [22, 23]],
                [[0, 3], [5, 7], [11, 12], [14, 16], [19, 23]],
            ),
        ]:
            for u, v in [(a, b), (b, a), (a, empty_array), (empty_array, a)]:
                u = np.array(u, dtype="int32")
                v = np.array(v, dtype="int32")
                umv = utils.mask_intervals(u, v)
                assert u.dtype == umv.dtype
                x = self.intervals_to_set(umv)
                y = self.intervals_to_set(u) - self.intervals_to_set(v)
                assert x == y

    def test_mask_intervals_errors(self):
        with pytest.raises(ValueError):
            utils.mask_intervals(intervals=np.array([[50, 10]]), mask=np.array([[]]))
        with pytest.raises(ValueError):
            utils.mask_intervals(intervals=np.array([[]]), mask=np.array([[10, 5]]))


class TestGammaPdf:
    def test_gamma_pdf_basic(self):
        np.testing.assert_allclose(
            utils.gamma_pdf(x=1, a=1), scipy.stats.gamma.pdf(x=1, a=1)
        )
        np.testing.assert_allclose(
            utils.gamma_pdf(x=range(5), a=1), scipy.stats.gamma.pdf(x=range(5), a=1)
        )
        np.testing.assert_allclose(
            utils.gamma_pdf(x=1, a=range(1, 5)),
            scipy.stats.gamma.pdf(x=1, a=range(1, 5)),
        )

    @pytest.mark.parametrize(
        "x,a,loc,scale",
        [
            (1, 1, 0, 1),
            # Note that some of these values are (deliberately) invalid,
            # which should produce NaN in the corresponding output.
            (np.arange(5), np.linspace(0, 1, 5), 0, 1),
            (np.arange(5), np.linspace(0, 1, 5), 1, 2),
            (np.arange(5), np.linspace(0, 1, 5), np.arange(5), np.arange(5)),
            (np.eye(5), 1, 0, 1),
        ],
    )
    def test_gamma_pdf(self, x, a, loc, scale):
        np.testing.assert_allclose(
            utils.gamma_pdf(x=x, a=a, loc=loc, scale=scale),
            scipy.stats.gamma.pdf(x=x, a=a, loc=loc, scale=scale),
        )
