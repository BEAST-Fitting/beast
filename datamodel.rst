
..  _datamodel:

#####################################
Specifying Parameters in datamodel.py
#####################################



* project : string
      The name of the output results directory.
      ::
  
          project = 'beast_example_phat'
          
  * filters : list of strings
      Full filter names in the BEAST filter database.
      ::
      
          filters = ['HST_WFC3_F275W','HST_WFC3_F336W','HST_ACS_WFC_F475W', 
              'HST_ACS_WFC_F814W', 'HST_WFC3_F110W','HST_WFC3_F160W']
      
  * basefilters : list of strings
      Short versions of filter names.
      ::
      
          basefilters = ['F275W','F336W','F475W', 'F814W','F110W','F160W']
               
  * obs_colnames : list of strings
      Names of columns for filters matching the observed catalog. The input data MUST be in fluxes, NOT in magnitudes and the fluxes MUST be in normalized Vega units.
      ::
      
          obs_colnames = [ f.lower() + '_rate' for f in basefilters ]
          
  * obsfile : string
      Pathname to the observed catalog.
      ::
      
          obsfile = 'data/b15_4band_det_27_A.fits'
          
The BEAST generates artificial star test (AST) input files based on additional
input parameters from datamodel.py. The parameters which need to be 
specified (and example values) include:

  * ast_models_selected_per_age : integer
      Number of models to pick per age (Default = 70).
      ::
      
          ast_models_selected_per_age = 70  

  * ast_bands_above_maglimit : integer 
      Number of filters that must be above the magnitude limit for an AST to be included in the list (Default = 3).
      ::
      
          ast_bands_above_maglimit = 3  
                             

  * ast_realization_per_model : integer
      Number of Realizations of each included AST model to be put into the list. (Default = 20).
      ::
      
          ast_realization_per_model = 20
                             

  * ast_maglimit : float (single value or array with one value per filter)
      (1) option 1: [number] to change the number of mags fainter than
      the 90th percentile faintest star in the photometry catalog to be used for
      the mag cut (Default = 1).
      
      (2) option 2: [space-separated list of numbers] to set custom faint end limits
      (one value for each band).
      ::
      
          ast_maglimit = [1.] 

  * ast_with_positions :  (bool,optional)
      If True, the AST list is produced with X,Y positions. If False, the AST list is produced with only magnitudes.
      ::
      
          ast_with_positions = True
                         
  * ast_pixel_distribution : float (optional; used if ast_with_positions is True)
      Minimum pixel separation between AST position and catalog star used to determine the AST spatial distribution.
      ::
      
          ast_pixel_distribution = 10.0 

  * ast_reference_image : string (optional, but required if ast_with_positions is True and no X,Y information is present in the photometry catalog)	
      Name of the reference image used by DOLPHOT when running the measured photometry.	            
      ::
      
          ast_reference_image = None
          
Next, you must specify the parameters for the AST noise model. These
parameters (and example values) include:

  * astfile : string
      Pathname to the AST files (single camera ASTs).
      ::
      
          astfile = 'data/fake_stars_b15_27_all.hd5'

  * noisefile : string
      A name for the noise model.
      ::
      
          noisefile = project + '/' + project + '_noisemodel.hd5'

  * distance modulus to the galaxy : float
      ::
      
          distanceModulus = 24.47 * units.mag
          
  * ast_colnames : list of strings 
      Names of columns for filters in the AST catalog (matches basefilters; not necessary/recommended to change).
      ::
      
          ast_colnames = np.array(basefilters)
          
  * absflux calibration covariance matrix
      Currently specific to HST filters (not necessary/recommended to change).
      ::
      
          absflux_a_matrix = absflux_covmat.hst_frac_matrix(filters)        
          
Next, you must define the parameters for the grid of stellar models
used by the BEAST. These parameters (and example values) include:

  * log10(Age) : list of floats ([min,max,step]) to generate the isochrones in years.
      ::
      
          logt = [6.0, 10.13, 1.0]

  * Mass : not sampled, instead the isochrone supplied mass spacing is used.


  * Metallicity : list of floats
      ::
      
          z = [0.03, 0.019, 0.008, 0.004]

  * Isochrone Model Grid : current choices: Padova or MIST
      
      PadovaWeb() -- `modeltype` param for iso sets from ezpadova
      (choices: parsec12s_r14, parsec12s, 2010, 2008, 2002)
      
      MISTWeb() -- `rotation` param (choices: vvcrit0.0=default, vvcrit0.4)

      Default: PARSEC+CALIBRI
      ::
      
          oiso = isochrone.PadovaWeb(modeltype='parsec12s', filterPMS=True)
      Alternative: MIST -- v1, no rotation
      ::
      
          oiso = isochrone.MISTWeb()

  * Stellar Atmospheres library definition
      Options include Kurucz, `Tlusty`_, `BTSettl`_, Munari, Elodie and BaSel. You can also generates an object from the union of multiple individual libraries.
      ::
      
          osl = stellib.Tlusty() + stellib.Kurucz()
          
Finally, you must specify the parameters for the dust extinction grid. 
These parameters (and example values) include:

  * A(V): dust column in magnitudes ([min,max,step]), and prior model.
      ::
      
          avs = [0.0, 10.055, 1.0]
          av_prior_model = {'name': 'flat'}

  * R(V): dust average grain size ([min,max,step]), and prior model.
      ::
      
          rvs = [2.0,6.0,1.0]
          rv_prior_model = {'name': 'flat'}

  * fA: mixture factor between "MW" and "SMCBar" extinction curves ([min,max,step]), and prior model.
      ::
      
          fAs = [0.0,1.0, 0.25]
          fA_prior_model = {'name': 'flat'}
