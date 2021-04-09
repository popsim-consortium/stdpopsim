"""
Tests for the utils module.
"""
import unittest
import os
import pathlib
import tarfile
import tempfile
import numpy as np

from stdpopsim import utils


class TestValidDemographicModelId(unittest.TestCase):
    """
    Tests for the is_valid_demographic_model_id function.
    """

    def test_empty_string(self):
        self.assertFalse(utils.is_valid_demographic_model_id(""))

    def test_contains_spaces(self):
        bad_ids = ["CamelCase 4X19", " CamelCase4X19", "CamelCase4X19"]
        for bad_id in bad_ids:
            self.assertFalse(utils.is_valid_demographic_model_id(bad_id))

    def test_too_many_sections(self):
        bad_ids = [
            "_CamelCase_4X19",
            "CamelCase_4X19_",
            "CamelCase__4X19",
            "CamelCase_4X19_1234_5678",
        ]
        for bad_id in bad_ids:
            self.assertFalse(utils.is_valid_demographic_model_id(bad_id))

    def test_bad_name(self):
        bad_ids = ["camelCase_4X19", "1CamelCase_4X19", "Camel-Case_4X19"]
        for bad_id in bad_ids:
            self.assertFalse(utils.is_valid_demographic_model_id(bad_id))

    def test_bad_populations(self):
        bad_ids = ["CamelCase_0X19", "CamelCase_0.1X19", "CamelCase_OneX19"]
        for bad_id in bad_ids:
            self.assertFalse(utils.is_valid_demographic_model_id(bad_id))

    def test_bad_author(self):
        bad_ids = ["CamelCase_1019", "CamelCase_1-19", "CamelCase_1;19"]
        for bad_id in bad_ids:
            self.assertFalse(utils.is_valid_demographic_model_id(bad_id))

    def test_bad_year(self):
        bad_ids = [
            "CamelCase_10X1",
            "CamelCase_1X",
            "CamelCase_1X2001",
            "CamelCase_1Xfive",
        ]
        for bad_id in bad_ids:
            self.assertFalse(utils.is_valid_demographic_model_id(bad_id))

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
            self.assertTrue(utils.is_valid_demographic_model_id(good_id))


class TestValidGeneticMapId(unittest.TestCase):
    """
    Tests for the is_valid_demographic_model_id function.
    """

    def test_empty_string(self):
        self.assertFalse(utils.is_valid_genetic_map_id(""))

    def test_contains_spaces(self):
        bad_ids = ["CamelCase ABCD", " CamelCase_ABCD", "CamelCase_ABCD "]
        for bad_id in bad_ids:
            self.assertFalse(utils.is_valid_genetic_map_id(bad_id))

    def test_bad_name(self):
        bad_ids = ["camelCase_ABCD", "1CamelCase_ABCD", "Camel-Case_ABCD"]
        for bad_id in bad_ids:
            self.assertFalse(utils.is_valid_genetic_map_id(bad_id))

    def test_bad_assembly(self):
        bad_ids = [
            "CamelCase_AB CD",
            "CamelCase_",
            "CamelCase_AB-CD",
            "CamelCase_AB.CD",
        ]
        for bad_id in bad_ids:
            self.assertFalse(utils.is_valid_genetic_map_id(bad_id))

    def test_good_ids(self):
        good_ids = ["CamelCase_ABCD", "CamelCase_a", "CamelCase_1000", "C_ABCD", "C_X"]
        for good_id in good_ids:
            self.assertTrue(utils.is_valid_genetic_map_id(good_id))


class TestValidSpeciesId(unittest.TestCase):
    """
    Tests for the is_valid_demographic_model_id function.
    """

    def test_empty_string(self):
        self.assertFalse(utils.is_valid_species_id(""))

    def test_contains_spaces(self):
        bad_ids = ["Cam Cas", " CamCas", "CamCas "]
        for bad_id in bad_ids:
            self.assertFalse(utils.is_valid_species_id(bad_id))

    def test_bad_ids(self):
        bad_ids = ["camCas", "camcas", "CAMCAS", "C1mCas", "Camcas", "CamelCase"]
        for bad_id in bad_ids:
            self.assertFalse(utils.is_valid_species_id(bad_id))

    def test_good_ids(self):
        bad_ids = ["CamCas", "HomSap"]
        for bad_id in bad_ids:
            self.assertTrue(utils.is_valid_species_id(bad_id))


class TestValidSpeciesName(unittest.TestCase):
    """
    Tests for the is_valid_demographic_model_id function.
    """

    def test_empty_string(self):
        self.assertFalse(utils.is_valid_species_name(""))

    def test_extra_spaces(self):
        bad_names = ["Pan  pan", " Pan pan", "Pan pan "]
        for bad_name in bad_names:
            self.assertFalse(utils.is_valid_species_name(bad_name))

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
            self.assertFalse(utils.is_valid_species_name(bad_name))

    def test_good_names(self):
        good_names = ["Pan pan", "Homo sapiens"]
        for good_name in good_names:
            self.assertTrue(utils.is_valid_species_name(good_name))

    @unittest.skip("Implement more flexible species name")
    def test_three_or_more(self):
        bad_names = ["Pan pan pan"]
        for bad_name in bad_names:
            self.assertTrue(utils.is_valid_species_name(bad_name))


class TestValidSpeciesCommonName(unittest.TestCase):
    """
    Tests for the is_valid_demographic_model_id function.
    """

    def test_empty_string(self):
        self.assertFalse(utils.is_valid_species_common_name(""))

    def test_bad_common_names(self):
        # TODO add more examples when we define restrictions (#330)
        bad_common_names = ["0", "human"]
        for bad_common_name in bad_common_names:
            self.assertFalse(utils.is_valid_species_common_name(bad_common_name))

    def test_good_common_names(self):
        # TODO add more examples when we define restrictions (#330)
        good_names = ["Human", "Stuff and things"]
        for good_name in good_names:
            self.assertTrue(utils.is_valid_species_common_name(good_name))


class TestDownload(unittest.TestCase):
    def test_download_local_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = pathlib.Path(tmpdir)
            input_filename = tmpdir / "input.foo"
            output_filename = tmpdir / "output.foo"
            with open(input_filename, "w") as f:
                print("foofoofoo", file=f)
            url = input_filename.resolve().as_uri()
            utils.download(url, output_filename)
            self.assertTrue(output_filename.exists())
            with open(output_filename) as f:
                contents = f.read()
            self.assertEqual(contents.strip(), "foofoofoo")

    def test_download_nonexistent(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = pathlib.Path(tmpdir)
            input_filename = tmpdir / "nonexistant"
            output_filename = tmpdir / "output.foo"
            for url in (
                input_filename.resolve().as_uri(),  # local file
                "http://example.com/nonexistant",  # remote file
            ):
                with self.assertRaises(OSError):
                    utils.download(url, output_filename)
                self.assertFalse(output_filename.exists())


class TestSha256(unittest.TestCase):
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
            self.assertEqual(
                sha256sum,
                # Calculated with `sha256sum` from GNU coreutils.
                "269dce1a5bb90188b2d9cf542a7c30e410c7d8251e34a97bfea56062df51ae23",
            )

    def test_sha256_big(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            filename = pathlib.Path(tmpdir) / "test.foo"
            # write binary file, to avoid windows/mac/unix line ending differences
            with open(filename, "wb") as f:
                for i in range(1024 * 1024):
                    f.write(str(i).encode())
            sha256sum = utils.sha256(filename)
            self.assertEqual(
                sha256sum,
                # Calculated with `sha256sum` from GNU coreutils.
                "995e0fde646f7dc98423af9a862be96014574bfa76be1186b484f796c4e58533",
            )


class TestCd(unittest.TestCase):
    def test_cd_context_manager(self):
        # On Mac, the path we enter with "cd" may differ from
        # the path we get with cwd() due to symlinks. So we
        # resolve all paths here to ignore symlink-only differences.
        old_cwd = pathlib.Path.cwd().resolve()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = pathlib.Path(tmpdir).resolve()
            self.assertNotEqual(old_cwd, tmpdir)
            with utils.cd(tmpdir):
                self.assertEqual(pathlib.Path.cwd().resolve(), tmpdir)
            self.assertEqual(pathlib.Path.cwd().resolve(), old_cwd)


def rm_f(filename):
    try:
        os.unlink(filename)
    except OSError:
        pass


class TestUntar(unittest.TestCase):
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
                self.assertTrue(filename.exists())
                # Check that extracted files have the same checksums as
                # the files we put in the tar.
                self.assertEqual(utils.sha256(filename), sha)

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
                with self.assertRaises(ValueError):
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
                with self.assertRaises(ValueError):
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
                with self.assertRaises(ValueError):
                    utils.untar(tar, dest)
                rm_f(filename)
                rm_f(tar)


class TestIntervalUtilities(unittest.TestCase):
    def test_intervals_array_shape(self):
        # shape (n, x) where x >= 2
        bad_arrays = (np.array([10, 20]), np.array([[[10, 20]]]))
        for bad_array in bad_arrays:
            with self.assertRaises(ValueError):
                utils.check_intervals_array_shape(intervals=bad_array)
            with self.assertRaises(ValueError):
                utils.build_intervals_array(intervals=bad_array)
        good_arrays = (
            np.array([[10, 20], [30, 40]]),
            np.array([[10, 20, 1], [30, 40, 2]]),
        )
        for good_array in good_arrays:
            utils.check_intervals_array_shape(intervals=good_array)
            built_array = utils.build_intervals_array(intervals=good_array)
            self.assertTrue((built_array == good_array).all())

    def test_invalid_intervals(self):
        # right <= left, overlapping, uncastable float->int
        invalid_arrays = (
            np.array([[20, 10]]),
            np.array([[10, 20], [15, 20]]),
            np.array([[10, 20], [13, 41]]),
        )
        for invalid_array in invalid_arrays:
            with self.assertRaises(ValueError):
                utils.build_intervals_array(intervals=invalid_array)
            with self.assertRaises(ValueError):
                utils.check_intervals_validity(intervals=invalid_array)
        # interval outside [start, end)
        with self.assertRaises(ValueError):
            utils.build_intervals_array(intervals=[[10, 50]], start=20)
        with self.assertRaises(ValueError):
            utils.check_intervals_validity(intervals=np.array([[10, 50]]), end=30)

    def test_intervals_casting(self):
        castable_intervals = (
            [[10, 20], [20, 40]],
            [[1, 2, 2], [2, 3, 2]],
            np.array([[10, 20]]),
        )
        for intervals in castable_intervals:
            casted = utils.build_intervals_array(intervals)
            self.assertTrue(isinstance(casted, np.ndarray) and casted.dtype == np.int64)
            self.assertTrue(casted.shape[0] == len(intervals))
            self.assertTrue(casted.shape[1] >= 2)

    def test_interval_sorting(self):
        unsorted_intervals = (
            np.array([[20, 30], [10, 20]]),
            np.array([[10, 20], [20, 30], [2, 3]]),
        )
        for intervals in unsorted_intervals:
            casted = utils.build_intervals_array(intervals)
            self.assertFalse(np.all(np.diff(intervals[:, 0]) >= 0))
            self.assertTrue(np.all(np.diff(casted[:, 0]) >= 0))
