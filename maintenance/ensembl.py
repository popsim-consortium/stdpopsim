"""
Utilites for working with the ensembl Rest API.
"""
import json
import time
import logging

import urllib.parse
import urllib.request


logger = logging.getLogger("ensembl")


class EnsemblRestClient:
    """
    A client for the Ensembl REST API. Based on the example code
    given at
    https://github.com/Ensembl/ensembl-rest/wiki/Example-Python-Client


    The ``max_requests_per_second`` parameter is the number of requests per second.
    """

    def __init__(self, server="http://rest.ensembl.org", max_requests_per_second=15):
        self.server = server
        self.max_requests_per_second = max_requests_per_second
        self.num_requests = 0
        self.last_request_time = 0

    def _make_request(self, endpoint, headers=None, params=None):
        if headers is None:
            headers = {}
        else:
            headers = dict(headers)
        if "Content-type" not in headers:
            headers["Content-type"] = "Application/json"

        quoted_endpoint = urllib.parse.quote(endpoint)
        url = urllib.parse.urljoin(self.server, quoted_endpoint)
        parsed_url = urllib.parse.urlparse(url)
        new_url = list(parsed_url)
        if params is not None:
            new_url[4] = urllib.parse.urlencode(params, doseq=True)
        return urllib.request.Request(urllib.parse.urlunparse(new_url), headers=headers)

    def _sleep_if_needed(self):
        # check if we need to rate limit ourselves
        if self.num_requests >= self.max_requests_per_second:
            delta = time.time() - self.last_request_time
            if delta < 1:
                logger.info("Rate limiting REST API")
                time.sleep(1 - delta)
            self.num_requests = 0
        else:
            self.num_requests += 1
        self.last_request_time = time.time()

    def get(self, endpoint, headers=None, params=None):
        """
        Runs a HTTP get request to the specified endpoint with the
        specified headers and parameters, and returns a dictionary
        consisting of the response.
        """
        self._sleep_if_needed()
        request = self._make_request(endpoint, headers, params)
        logger.info("making request to %s", request.full_url)
        response = urllib.request.urlopen(request)
        content = response.read()
        logger.debug("Response: %s", content)
        data = json.loads(content)
        return data

    def get_release(self):
        """
        Returns the current ensembl release number.
        """
        output = self.get(endpoint="/info/data")
        releases = output["releases"]
        # Docs say: may return more than one release
        # (unfrequent non-standard Ensembl configuration).
        assert len(releases) == 1
        return releases[0]

    def get_species_data(self, ensembl_id):
        """
        Returns species information for the specified ensembl_id.
        """
        output = self.get(endpoint=f"/info/genomes/{ensembl_id}")
        return output

    def get_genome_data(self, ensembl_id):
        """
        Returns the genome data for the specified Ensembl species
        ID (e.g., "homo_sapiens"). This is a subset of what
        is returned by the Ensembl API, to restrict to the information
        that we're interested in for stdpopsim.
        """
        output = self.get(
            endpoint=f"/info/assembly/{ensembl_id}", params={"synonyms": "1"}
        )
        data = {
            "assembly_accession": output["assembly_accession"],
            "assembly_name": output["assembly_name"],
        }
        chromosomes = {}
        for region in output["top_level_region"]:
            if region["coord_system"] in ("chromosome", "primary_assembly"):
                synonyms = []
                for record in region.get("synonyms", []):
                    # We're only interested in UCSC synonyms
                    if record["dbname"] == "UCSC":
                        synonyms.append(record["name"])
                chromosomes[region["name"]] = {
                    "length": region["length"],
                    "synonyms": synonyms,
                }
        # Make sure the chromosomes are sorted by name
        data["chromosomes"] = {
            name: chromosomes[name]
            for name in output["karyotype"]
            if name in chromosomes
        }
        return data
