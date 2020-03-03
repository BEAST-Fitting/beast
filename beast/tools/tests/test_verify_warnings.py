import pytest
from beast.tools import verify_params

class datamodel_mock:
    """A simple mock datamodel"""
    obsfile = 'photometry.fits'
    astfile = 'ast.fits'
    logt = [6.0, 10.13, 0.05]
    age_prior_model = {'name': 'flat_log'}
    mass_prior_model = {'name': 'kroupa'}
    z = [0.00452, 0.00254, 0.00143]
    met_prior_model = {'name': 'flat'}
    avs = [0.0, 4.0, 0.1]
    av_prior_model = {'name': 'lognormal',
                      'max_pos': 0.76,
                      'sigma': 0.73,
                      'N': 10.}
    rvs = [3.0, 6.0, 1.0]
    rv_prior_model = {'name': 'flat'}
    fAs = [0.0, 1.0, 0.25]
    fA_prior_model = {'name': 'flat'}

class datamodel_mock_nofA(datamodel_mock):
    """Mock datamodel w/ fAs=None"""
    fAs = None

class datamodel_mock_allowwarn(datamodel_mock_nofA):
    """Mock datamodel w/ fAs=None and allow_warnings = True"""
    allow_warnings = True

def test_verifyparams_nowarning():
    """Test: verify_params for case of no warnings or exceptions."""
    with pytest.warns(None) as record:
        verify_params.verify_input_format(datamodel_mock())
    assert len(record) == 0

def test_verifyparams_error():
    """Test: verify_params for case of warning raising exception."""
    with pytest.raises(UserWarning) as exc:
        verify_params.verify_input_format(datamodel_mock_nofA())
    assert exc.value.args[0] == "fAs is not defined."

def test_verifyparams_allowwarn():
    """Test: verify_params for case of warning with no exception."""
    with pytest.warns(UserWarning, match="fAs is not defined."):
        verify_params.verify_input_format(datamodel_mock_allowwarn())
