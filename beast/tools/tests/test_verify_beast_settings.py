import pytest
import warnings
from beast.physicsmodel.stars import isochrone
from beast.tools import verify_beast_settings


class settings_mock:
    """A simple mock beast_settings"""

    obsfile = "photometry.fits"
    astfile = "ast.fits"
    logt = [6.0, 10.13, 0.05]
    age_prior_model = {"name": "flat_log"}
    mass_prior_model = {"name": "kroupa"}
    z = [0.00452, 0.00254, 0.00143]
    met_prior_model = {"name": "flat"}
    avs = [0.0, 4.0, 0.1]
    av_prior_model = {"name": "lognormal", "max_pos": 0.76, "sigma": 0.73, "N": 10.0}
    rvs = [3.0, 6.0, 1.0]
    rv_prior_model = {"name": "flat"}
    fAs = [0.0, 1.0, 0.25]
    fA_prior_model = {"name": "flat"}
    distances = [50.0, 70.0, 5.0]
    distance_prior_model = {"name": "flat"}
    oiso = isochrone.MISTWeb()


class settings_mock_nofA(settings_mock):
    """Mock beast_settings w/ fAs=None"""

    fAs = None
    fA_prior_model = None


class settings_mock_allowwarn(settings_mock_nofA):
    """Mock beast_settings w/ fAs=None and allow_warnings = True"""

    allow_warnings = True


def test_verifyparams_nowarning():
    """Test: verify_beast_settings for case of no warnings or exceptions."""
    # with pytest.warns(None) as record:
    with warnings.catch_warnings(record=True) as record:
        verify_beast_settings.verify_input_format(settings_mock())
    assert len(record) == 0


def test_verifyparams_nofA():
    """Test: verify_beast_settings for case with fA=None (allowed)."""
    with warnings.catch_warnings(record=True) as record:
        verify_beast_settings.verify_input_format(settings_mock_nofA())
    assert len(record) == 0


@pytest.mark.skip(reason="failing, test needs updating as fAs==None allowed")
def test_verifyparams_allowwarn():
    """Test: verify_beast_settings for case of warning with no exception."""
    with pytest.warns(UserWarning, match="fAs is not defined."):
        verify_beast_settings.verify_input_format(settings_mock_allowwarn())


class settings_mock_RV(settings_mock):
    """Mock beast_settings w/ single-valued R_V"""

    rvs = [3.1, 3.1, 1.0]


class settings_mock_noallowRV(settings_mock):
    """Mock beast_settings w/ single-valued R_V"""

    rvs = [3.1, 3.1, 1.0]
    allow_warnings = False


class settings_mock_allowwarnRV(settings_mock_RV):
    """Mock beast_settings w/ single-valued R_V and allow_warnings = True"""

    allow_warnings = True


def test_verifyparams_errorRV():
    """Test: verify_beast_settings for case of warning raising exception."""
    with pytest.raises(UserWarning) as exc:
        verify_beast_settings.verify_input_format(settings_mock_RV())
    assert exc.value.args[0] == "Note: rvs grid is single-valued."


def test_verifyparams_noallowRV():
    """Test: verify_beast_settings when warn raising except w/ allow_warnings=False"""
    with pytest.raises(UserWarning) as exc:
        verify_beast_settings.verify_input_format(settings_mock_noallowRV())
    assert exc.value.args[0] == "Note: rvs grid is single-valued."


def test_verifyparams_allowwarnRV():
    """Test: verify_beast_settings for case of warning with no exception."""
    with pytest.warns(UserWarning, match="Note: rvs grid is single-valued."):
        verify_beast_settings.verify_input_format(settings_mock_allowwarnRV())


class settings_mock_dist_explicit(settings_mock):
    """Mock beast_settings w/ an explicit, non-uniform distance grid (>3 values)"""

    distances = [55.0, 58.0, 60.0, 61.0, 62.0, 63.0, 65.0, 70.0]


class settings_mock_dist_single(settings_mock):
    """Mock beast_settings w/ a single distance"""

    distances = [60.0]


class settings_mock_dist_bare(settings_mock):
    """Mock beast_settings w/ a single distance given as a bare number"""

    distances = 60.0


class settings_mock_dist_int(settings_mock):
    """Mock beast_settings w/ an integer distance grid"""

    distances = [50, 70, 5]


class settings_mock_dist_two(settings_mock):
    """Mock beast_settings w/ an unusable two-element distance list"""

    distances = [50.0, 70.0]


class settings_mock_dist_dup(settings_mock):
    """Mock beast_settings w/ duplicate values in an explicit distance grid"""

    distances = [55.0, 60.0, 60.0, 65.0, 70.0]


class settings_mock_dist_negative(settings_mock):
    """Mock beast_settings w/ an unphysical negative distance"""

    distances = [-10.0, 55.0, 60.0, 65.0]


class settings_mock_dist_str(settings_mock):
    """Mock beast_settings w/ a non-numeric entry in the distance grid"""

    distances = [55.0, 60.0, "65.0", 70.0]


class settings_mock_dist_none(settings_mock):
    """Mock beast_settings w/ distances undefined"""

    distances = None


class settings_mock_dist_none_allowwarn(settings_mock_dist_none):
    """Mock beast_settings w/ distances undefined and allow_warnings = True"""

    allow_warnings = True


def test_verifyparams_dist_explicit():
    """Test: an explicit >3-value distance grid passes verification."""
    with warnings.catch_warnings(record=True) as record:
        verify_beast_settings.verify_input_format(settings_mock_dist_explicit())
    assert len(record) == 0


def test_verifyparams_dist_single():
    """Test: a single-element distance list passes verification."""
    with warnings.catch_warnings(record=True) as record:
        verify_beast_settings.verify_input_format(settings_mock_dist_single())
    assert len(record) == 0


def test_verifyparams_dist_bare():
    """Test: a bare numeric distance passes verification."""
    with warnings.catch_warnings(record=True) as record:
        verify_beast_settings.verify_input_format(settings_mock_dist_bare())
    assert len(record) == 0


def test_verifyparams_dist_int():
    """Test: an integer-valued distance grid passes verification."""
    with warnings.catch_warnings(record=True) as record:
        verify_beast_settings.verify_input_format(settings_mock_dist_int())
    assert len(record) == 0


def test_verifyparams_dist_two():
    """Test: a two-element distance list is rejected."""
    with pytest.raises(ValueError, match="distances needs to be a single number"):
        verify_beast_settings.verify_input_format(settings_mock_dist_two())


def test_verifyparams_dist_duplicates():
    """Test: duplicate values in an explicit distance grid are rejected."""
    with pytest.raises(ValueError, match="distances explicit grid has duplicate"):
        verify_beast_settings.verify_input_format(settings_mock_dist_dup())


def test_verifyparams_dist_negative():
    """Test: a negative distance is rejected."""
    with pytest.raises(ValueError, match="distances min value not physical."):
        verify_beast_settings.verify_input_format(settings_mock_dist_negative())


def test_verifyparams_dist_str():
    """Test: a non-numeric distance entry is rejected."""
    with pytest.raises(TypeError, match="distances is not in the right format"):
        verify_beast_settings.verify_input_format(settings_mock_dist_str())


def test_verifyparams_dist_none():
    """Test: undefined distances raises, as the warning is turned into an error."""
    with pytest.raises(UserWarning) as exc:
        verify_beast_settings.verify_input_format(settings_mock_dist_none())
    assert exc.value.args[0] == "distances is not defined."


def test_verifyparams_dist_none_allowwarn():
    """Test: undefined distances only warns when allow_warnings is set."""
    with pytest.warns(UserWarning, match="distances is not defined."):
        verify_beast_settings.verify_input_format(settings_mock_dist_none_allowwarn())
