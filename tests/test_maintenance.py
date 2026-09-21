"""
Tests for the maintenance utilities.
"""

from unittest import mock
import importlib.util
import json
import sys
import urllib
import urllib.request

import stdpopsim
import maintenance as maint
from maintenance import main


class TestCatalogStub:
    def test_gene_conversion(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        path = tmp_path / "stub_species"
        path.mkdir()
        (tmp_path / "tests").mkdir()
        genome_data = {
            "assembly_accession": "test_accession",
            "assembly_name": "test_assembly",
            "assembly_source": "ensembl",
            "assembly_build_version": "test_version",
            "chromosomes": {
                "1": {"length": 1000, "synonyms": []},
                "X": {"length": 500, "synonyms": ["chrX"]},
            },
        }
        (path / "genome_data.py").write_text(f"data = {genome_data!r}\n")
        main.write_catalog_stub(
            path=path,
            sps_id="TesSpe",
            ensembl_id="test_species",
            species_data={
                "scientific_name": "Test species",
                "display_name": "Test species",
            },
            genome_data=genome_data,
        )
        # Both the catalog entry and the independent QC test stub must parse.
        for generated in [path / "species.py", tmp_path / "tests/test_TesSpe.py"]:
            compile(generated.read_text(), str(generated), "exec")

        spec = importlib.util.spec_from_file_location(
            "stub_species", path / "__init__.py"
        )
        package = importlib.util.module_from_spec(spec)
        with (
            mock.patch.dict(sys.modules),
            mock.patch("stdpopsim.register_species") as register,
            mock.patch(
                "stdpopsim.Genome.from_data", wraps=stdpopsim.Genome.from_data
            ) as from_data,
        ):
            sys.modules[spec.name] = package
            spec.loader.exec_module(package)
            species = package.species
            expected = {"1": None, "X": None}
            assert species._gene_conversion_fraction == expected
            assert species._gene_conversion_length == expected
            assert from_data.call_args.kwargs["gene_conversion_fraction"] == expected
            assert from_data.call_args.kwargs["gene_conversion_length"] == expected
            register.assert_called_once_with(species._species)
            for chrom in species._genome.chromosomes:
                assert chrom.gene_conversion_fraction is None
                assert chrom.gene_conversion_length is None


class MockedResponse:
    def __init__(self, value={}):
        self.value = value

    def read(self):
        return json.dumps(self.value)


class TestEnsemblClient:
    """
    Tests for the Ensembl rest client.
    """

    def test_defaults(self):
        client = maint.EnsemblRestClient()
        assert client.server == "http://rest.ensembl.org"
        assert client.max_requests_per_second == 15

    def test_request_params(self):
        client = maint.EnsemblRestClient("http://example.org")
        request = client._make_request("a/b", params={"a": "b"})
        assert request.full_url == "http://example.org/a/b?a=b"
        request = client._make_request("", params={"a": 0, "b": 1})
        assert request.full_url == "http://example.org?a=0&b=1"
        request = client._make_request("a b", params={"a": "a b"})
        assert request.full_url == "http://example.org/a%20b?a=a+b"

    def test_basic_example(self):
        test_server = "http://example.com"
        client = maint.EnsemblRestClient(test_server)
        assert client.server == test_server

        returned_response = MockedResponse()
        with mock.patch(
            "urllib.request.urlopen", autospec=True, return_value=returned_response
        ) as mocked_open:
            value = client.get("test_endpoint")
            assert value == {}
            mocked_open.assert_called_once()
            request_sent = mocked_open.call_args[0][0]
            assert isinstance(request_sent, urllib.request.Request)
            assert request_sent.full_url == test_server + "/test_endpoint"
            assert request_sent.headers == {"Content-type": "Application/json"}

    def test_extra_headers(self):
        test_server = "http://example.com"
        client = maint.EnsemblRestClient(test_server)

        returned_response = MockedResponse()
        with mock.patch(
            "urllib.request.urlopen", autospec=True, return_value=returned_response
        ) as mocked_open:
            extra_headers = {"A": "sdf", "B": "xyz"}
            value = client.get("/test_endpoint", extra_headers)
            # We don't modify the input
            assert len(extra_headers) == 2
            assert value == {}
            mocked_open.assert_called_once()
            request_sent = mocked_open.call_args[0][0]
            assert isinstance(request_sent, urllib.request.Request)
            assert request_sent.full_url == test_server + "/test_endpoint"
            extra_headers["Content-type"] = "Application/json"
            assert request_sent.headers == extra_headers

    def test_rate_limit(self):
        client = maint.EnsemblRestClient(max_requests_per_second=1)
        with mock.patch("time.sleep", autospec=True) as mocked_sleep:
            client._sleep_if_needed()
            assert mocked_sleep.call_count == 0
            client._sleep_if_needed()
            assert mocked_sleep.call_count == 1
            client._sleep_if_needed()
            assert mocked_sleep.call_count == 1
            client._sleep_if_needed()
            assert mocked_sleep.call_count == 2

        client = maint.EnsemblRestClient(max_requests_per_second=3)
        with mock.patch("time.sleep", autospec=True) as mocked_sleep:
            client._sleep_if_needed()
            assert mocked_sleep.call_count == 0
            client._sleep_if_needed()
            assert mocked_sleep.call_count == 0
            client._sleep_if_needed()
            assert mocked_sleep.call_count == 0
            client._sleep_if_needed()
            assert mocked_sleep.call_count == 1
