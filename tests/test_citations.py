"""
Test cases for methods related to citations.
"""

from unittest import mock
import urllib.error

import stdpopsim
import pytest


class mocked_response:
    """
    This exists solely to mock out the response from the
    urllib.request call so that we can test the bibtex
    retrieval. The enter and exit methods are for calling
    it as a context manager."""

    def __init__(self, code, text):
        self.code = code
        self.text = text

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    def read(self):
        return self.text

    def close(self):
        pass


class TestFetchBibtex:
    """
    Test the fetching of bibtex files."""

    def test_get_bibtex_success(self):
        # Tests a success
        citation = stdpopsim.Citation(doi="DOI", author="Authors", year="2000")
        with mock.patch(
            "urllib.request.urlopen",
            return_value=mocked_response(code=200, text=b"test"),
            autospec=True,
        ):
            with mock.patch("urllib.request.Request", autospec=True):
                bib = citation.fetch_bibtex()
                assert bib == "test"

    def test_get_bibtex_bad_connection(self):
        # Tests an invalid URL
        # Asserts that it raises a value error.
        citation = stdpopsim.Citation(doi="DOI", author="Authors", year="2000")
        with pytest.raises(ValueError):
            citation.fetch_bibtex()

    def test_get_bibtex_404(self):
        # Tests a 404.  The previous URL here, http://doi.org/url-does-not-exist,
        # started just not responding.
        citation = stdpopsim.Citation(
            doi="http://github.com/popsim-consortium/this_is_not_a_repo",
            author="Authors",
            year="2000",
        )
        with pytest.raises(urllib.error.HTTPError):
            citation.fetch_bibtex()


class TestEverythingHappensForAReason:
    """
    Test that citations have a reason.
    """

    def test_reasons_for_engine_citations(self):
        for engine in stdpopsim.all_engines():
            for citation in engine.citations:
                assert len(citation.reasons) > 0, (
                    f"No reason given for '{citation.author}' citation "
                    f"in engine {engine.id}"
                )

    def test_reason_for_genetic_map_citations(self):
        for genetic_map in stdpopsim.all_genetic_maps():
            for citation in genetic_map.citations:
                assert len(citation.reasons) > 0, (
                    f"No reason given for '{citation.author}' citation "
                    f"in genetic map "
                    f"{genetic_map.species.id}/{genetic_map.id}"
                )

    def test_reason_for_model_citations(self):
        for model in stdpopsim.all_demographic_models():
            for citation in model.citations:
                assert len(citation.reasons) > 0, (
                    f"No reason given for '{citation.author}' citation "
                    f"in model {model.id}"
                )

    def test_reason_for_dfe_citations(self):
        for dfe in stdpopsim.all_dfes():
            for citation in dfe.citations:
                assert len(citation.reasons) > 0, (
                    f"No reason given for '{citation.author}' citation "
                    f"in dfe {dfe.id}"
                )
