import numpy as np
import pytest

from beast.physicsmodel.prior_weights_dust import PriorWeightsDust

def test_dust_prior_weights():
    """
    Test the dust prior weights
    """

    av_vals = [0., 2., 4., 6., 8., 10.]
    rv_vals = [2., 3., 4., 5., 6.]
    fA_vals = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    av_model = {"name":"flat"}
    rv_model = {"name":"flat"}
    fA_model = {"name":"flat"}

    prior_weights = PriorWeightsDust(av_vals, av_model, rv_vals, rv_model, fA_vals, fA_model)

    # test if priors are set correctly when creating an object of the class PriorWeightsDust
    np.testing.assert_allclose(prior_weights.av_priors,np.full((len(av_vals)), 1.0), err_msg="A problem occurred with setting the A(V) priors when creating a PriorWeightsDust object. They should be set to flat by default.")
    np.testing.assert_allclose(prior_weights.rv_priors,np.full((len(rv_vals)), 1.0), err_msg="A problem occurred with setting the R(V) priors when creating a PriorWeightsDust object. They should be set to flat by default.")
    np.testing.assert_allclose(prior_weights.fA_priors,np.full((len(fA_vals)), 1.0), err_msg="A problem occurred with setting the f_A priors when creating a PriorWeightsDust object. They should be set to flat by default.")


    # Test A(V) priors
    # test default: model=flat
    prior_weights.set_av_weights()
    np.testing.assert_allclose(prior_weights.av_priors,np.full((len(av_vals)), 1.0), err_msg="A problem occurred while setting the A(V) priors. They should be flat by default.")

    # test other models (with random input values)
    expected_lognormal_av = [0., 2.00684868, 1.57828555, 1.09756182, 0.7677122, 0.54959176]
    prior_weights.set_av_weights(model={"name" : "lognormal", "max_pos" : 2.0, "sigma" : 1.0})
    np.testing.assert_allclose(prior_weights.av_priors, expected_lognormal_av, err_msg="A problem occurred while setting the A(V) priors to lognormal.")

    expected_two_lognormal_av = [0.0, 5.962474, 3.07642414e-02, 4.46347150e-03, 1.60916384e-03, 6.89125670e-04]
    prior_weights.set_av_weights(model={"name" : "two_lognormal", "max_pos1" : 0.2, "max_pos2": 2.0, "sigma1" : 1.0, "sigma2" : 0.2, "N1_to_N2" : 1.0 / 5.0})
    np.testing.assert_allclose(prior_weights.av_priors, expected_two_lognormal_av, err_msg="A problem occurred while setting the A(V) priors to two_lognormal.")

    expected_exponential_av = [5.18802018, 7.02122180e-01, 9.50219041e-02, 1.28598163e-02, 1.74038688e-03, 2.35535752e-04]
    prior_weights.set_av_weights(model={"name" : "exponential", "a" : 1.0})
    np.testing.assert_allclose(prior_weights.av_priors, expected_exponential_av, err_msg="A problem occurred while setting the A(V) priors to exponential.")

    # test error output when unsupported model is given
    with pytest.raises(NotImplementedError) as exc_av:
        prior_weights.set_av_weights(model={"name" : "exp"})
    assert str(exc_av.value) == "**Error in setting the A(V) dust prior weights!****model exp not supported**"


    # Test R(V) priors
    prior_weights.set_rv_weights()
    np.testing.assert_allclose(prior_weights.rv_priors,np.full((len(rv_vals)), 1.0), err_msg="A problem occurred while setting the R(V) priors. They should be flat by default.")

    expected_lognormal_rv = [0.53984107, 2.48794534, 1.4922695, 0.4032949, 0.0766492]
    prior_weights.set_rv_weights(model={"name" : "lognormal", "max_pos" : 3.1, "sigma" : 0.25})
    np.testing.assert_allclose(prior_weights.rv_priors, expected_lognormal_rv, err_msg="A problem occurred while setting the R(V) priors to lognormal.")

    expected_two_lognormal_rv = [5.19299866e-04, 1.85412702, 1.32008619, 1.29608363, 5.29183855e-01]
    prior_weights.set_rv_weights(model={"name" : "two_lognormal", "max_pos1" : 3.1, "max_pos2": 4.5, "sigma1" : 0.1, "sigma2" : 0.2, "N1_to_N2" : 2.0 / 5.0})
    np.testing.assert_allclose(prior_weights.rv_priors, expected_two_lognormal_rv, err_msg="A problem occurred while setting the R(V) priors to two_lognormal.")

    with pytest.raises(NotImplementedError) as exc_rv:
        prior_weights.set_rv_weights(model={"name" : "exponential"})
    assert str(exc_rv.value) == "**Error in setting the R(V) dust prior weights!****model exponential not supported**"


    # Test f_A priors
    prior_weights.set_fA_weights()
    np.testing.assert_allclose(prior_weights.fA_priors,np.full((len(fA_vals)), 1.0), err_msg="A problem occurred while setting the f_A priors. They should be flat by default.")

    expected_lognormal_fA = [0.0, 1.01292523e-41, 2.01507540e-10, 8.71092362e-02, 5.46004145, 4.52849316e-01]
    prior_weights.set_fA_weights(model={"name" : "lognormal", "max_pos" : 0.8, "sigma" : 0.1})
    np.testing.assert_allclose(prior_weights.fA_priors, expected_lognormal_fA, err_msg="A problem occurred while setting the f_A priors to lognormal.")

    expected_two_lognormal_fA = [0.0, 8.76235139e-10, 7.80598377e-03, 1.12556575, 3.16704186, 1.69958640]
    prior_weights.set_fA_weights(model={"name" : "two_lognormal", "max_pos1" : 0.1, "max_pos2": 0.8, "sigma1" : 0.1, "sigma2" : 0.2, "N1_to_N2" : 2.0 / 5.0})
    np.testing.assert_allclose(prior_weights.fA_priors, expected_two_lognormal_fA, err_msg="A problem occurred while setting the f_A priors to two_lognormal.")

    with pytest.raises(NotImplementedError) as exc_fA:
        prior_weights.set_fA_weights(model={"name" : "exponential"})
    assert str(exc_fA.value) == "**Error in setting the f_A dust prior weights!****model exponential not supported**"
