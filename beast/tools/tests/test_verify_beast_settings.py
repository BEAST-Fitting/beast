import pytest
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


class settings_mock_nofA(settings_mock):
    """Mock beast_settings w/ fAs=None"""

    fAs = None


class settings_mock_allowwarn(settings_mock_nofA):
    """Mock beast_settings w/ fAs=None and allow_warnings = True"""

    allow_warnings = True


def test_verifyparams_nowarning():
    """Test: verify_beast_settings for case of no warnings or exceptions."""
    with pytest.warns(None) as record:
        verify_beast_settings.verify_input_format(settings_mock())
    assert len(record) == 0


def test_verifyparams_error():
    """Test: verify_beast_settings for case of warning raising exception."""
    with pytest.raises(UserWarning) as exc:
        verify_beast_settings.verify_input_format(settings_mock_nofA())
    assert exc.value.args[0] == "fAs is not defined."


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
