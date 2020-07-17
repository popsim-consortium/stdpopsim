"""
Tests for the maintenance utilities.
"""
import unittest
from unittest import mock
import json
import urllib
import urllib.request


import maintenance as maint


class MockedResponse:
    def __init__(self, value={}):
        self.value = value

    def read(self):
        return json.dumps(self.value)


class TestEnsemblClient(unittest.TestCase):
    """
    Tests for the Ensembl rest client.
    """

    def test_defaults(self):
        client = maint.EnsemblRestClient()
        self.assertEqual(client.server, "http://rest.ensembl.org")
        self.assertEqual(client.max_requests_per_second, 15)

    def test_request_params(self):
        client = maint.EnsemblRestClient("http://example.org")
        request = client._make_request("a/b", params={"a": "b"})
        self.assertEqual(request.full_url, "http://example.org/a/b?a=b")
        request = client._make_request("", params={"a": 0, "b": 1})
        self.assertEqual(request.full_url, "http://example.org?a=0&b=1")
        request = client._make_request("a b", params={"a": "a b"})
        self.assertEqual(request.full_url, "http://example.org/a%20b?a=a+b")

    def test_basic_example(self):
        test_server = "http://example.com"
        client = maint.EnsemblRestClient(test_server)
        self.assertEqual(client.server, test_server)

        returned_response = MockedResponse()
        with mock.patch(
            "urllib.request.urlopen", autospec=True, return_value=returned_response
        ) as mocked_open:
            value = client.get("test_endpoint")
            self.assertEqual(value, {})
            mocked_open.assert_called_once()
            request_sent = mocked_open.call_args[0][0]
            self.assertIsInstance(request_sent, urllib.request.Request)
            self.assertEqual(request_sent.full_url, test_server + "/test_endpoint")
            self.assertEqual(request_sent.headers, {"Content-type": "Application/json"})

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
            self.assertEqual(len(extra_headers), 2)
            self.assertEqual(value, {})
            mocked_open.assert_called_once()
            request_sent = mocked_open.call_args[0][0]
            self.assertIsInstance(request_sent, urllib.request.Request)
            self.assertEqual(request_sent.full_url, test_server + "/test_endpoint")
            extra_headers["Content-type"] = "Application/json"
            self.assertEqual(request_sent.headers, extra_headers)

    def test_rate_limit(self):
        client = maint.EnsemblRestClient(max_requests_per_second=1)
        with mock.patch("time.sleep", autospec=True) as mocked_sleep:
            client._sleep_if_needed()
            self.assertEqual(mocked_sleep.call_count, 0)
            client._sleep_if_needed()
            self.assertEqual(mocked_sleep.call_count, 1)
            client._sleep_if_needed()
            self.assertEqual(mocked_sleep.call_count, 1)
            client._sleep_if_needed()
            self.assertEqual(mocked_sleep.call_count, 2)

        client = maint.EnsemblRestClient(max_requests_per_second=3)
        with mock.patch("time.sleep", autospec=True) as mocked_sleep:
            client._sleep_if_needed()
            self.assertEqual(mocked_sleep.call_count, 0)
            client._sleep_if_needed()
            self.assertEqual(mocked_sleep.call_count, 0)
            client._sleep_if_needed()
            self.assertEqual(mocked_sleep.call_count, 0)
            client._sleep_if_needed()
            self.assertEqual(mocked_sleep.call_count, 1)
