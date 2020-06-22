from beast.tools.tests.beast_settings_mockup import beast_settings_mockup
from astropy.tests.helper import remote_data


@remote_data
def test_beast_settings():
    """
    Test that a given text file creates the expected beast_examples class.
    """

    # grab the settings
    settings = beast_settings_mockup()

    # assert it's the correct class
    assert isinstance(
        settings, beast_settings.beast_settings
    ), "Did not produce the correct class"
