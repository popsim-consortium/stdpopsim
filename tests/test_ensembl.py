import stdpopsim


# Make sure we don't update the release without realising it.
def test_version():
    release = stdpopsim.catalog.ensembl_info.release
    assert release == 111
