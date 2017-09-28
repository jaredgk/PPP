### Impute missing GT w/ beagle

# argparse call for estimate file
# phase_parser.add_argument('--estimate-file', help = 'Defines the estimated genotype frequency filename. Required for the beagle algorithm', default = 'estimated_gt', action = parser_confirm_no_file())


# Impute calls

# Assign the algorithm
algorithm_call_args = ['java', '-jar', 'bin/beagle.jar']

# Assign the arguments for the algorithm
likelihood_call_args = ['gtgl=' + phase_args.vcfname, 'out=' + phase_args.estimate_file]

logging.info('beagle estimate parameters assigned')

# beagle estimated genotype frequency subprocess call
likelihood_call = subprocess.Popen(algorithm_call_args + likelihood_call_args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
likelihood_out, likelihood_err = likelihood_call.communicate()

# Confirms that beagle finished without error
if likelihood_err:
    logging.error('Error creating the estimated genotype frequency file. Please check input file.')
    sys.exit('Error creating the estimated genotype frequency file. Please check input file.')

logging.info('beagle estimate file created')
