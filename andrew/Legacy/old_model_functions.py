def run (passed_arguments = []):
    '''

    Parameters
    ----------

    Returns
    -------

    Raises
    ------

    '''

    # Grab VCF arguments from command line
    model_args = model_creator_parser(passed_arguments)

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(model_args, 'model_creator')

    # Loop each model specified to confirm parameters are valid
    for current_model in model_args.model:

        # Check that a tree has been assigned to the model
        if current_model not in model_args.model_tree:
            print 'No tree assigned to: %s' % current_model
            sys.exit()

        # Check that pops have been assigned to the model
        if current_model not in model_args.pops:
            print 'No populations assigned to: %s' % current_model
            sys.exit()

        # Check if the number of populations has been assigned
        if model_args.assign_npop:
            # Check if the number of populations has been assigned for the current model
            if model_args.assign_npop[current_model]:
                if len(model_args.pops[current_model]) != model_args.assign_npop[current_model]:
                    print  len(model_args.pops[current_model]),  model_args.assign_npop[current_model]
                    print 'Number of populations (--npop) assigned to %s does not match the named popultions' % current_model
                    sys.exit()

        # Loop each population in the model
        for current_pop in model_args.pops[current_model]:

            # Check that the population has been assigned to at least one model
            if current_pop not in set(itertools.chain.from_iterable(model_args.pops.values())):
                print 'Population (%s) not assigned to any model' % current_pop
                sys.exit()

            # Check that inds have been assigned to the population
            if current_pop not in model_args.inds:
                print 'No individuals assigned to: %s' % current_pop
                sys.exit()

            # Check that inds have been assigned to the population
            if current_pop not in model_args.inds:
                print 'No individuals assigned to: %s' % current_pop
                sys.exit()

            # Check if the number of inds has been assigned
            if model_args.assign_nind:
                # Check if the number of inds has been assigned for the current pop
                if model_args.assign_nind[current_pop]:
                    if len(model_args.inds[current_pop]) != model_args.assign_nind[current_pop]:
                        print 'Number of individuals (--nind) assigned to %s does not match the named individuals' % current_pop
                        sys.exit()

    model_output = open(model_args.out, 'w')
    # Loop each model specified to generate the output
    for current_model in model_args.model:
        model_output.write('MODEL:\t%s\n' % current_model)
        model_output.write('  TREE:\t%s\n' % model_args.model_tree[current_model])
        model_output.write('  NPOP:\t%s\n' % len(model_args.pops[current_model]))

        # Loop each population in the model
        for current_pop in model_args.pops[current_model]:
            model_output.write('  POP:\t%s\n' % current_pop)
            model_output.write('    NIND:\t%s\n' % len(model_args.inds[current_pop]))
            model_output.write('    INDS:\t%s\n' % ', '.join(model_args.inds[current_pop]))

    model_output.close()
    
def read_model_file (filename):

    # Used to store each model within the file
    models_in_file = {}

    # Used to store the current model
    current_model = None
    # Used to store the expected number of populations in a model.
    model_npop = 0
    # Used to store the current population name
    pop_name = ''
    # Used to store the expected number of individuals in a population.
    pop_nind = 0

    with open(filename, 'r') as model_file:
        for model_line in model_file:
            # Remove whitespace
            model_line = model_line.strip()

            if 'MODEL:' in model_line:

                # Check if previous model is stored
                if current_model != None:

                    # Make sure the population counts match
                    if current_model.npop != model_npop:
                        print ('Number of populations found in %s does not match NPOP' % current_model.name)
                        sys.exit()

                    # Store model
                    models_in_file[current_model.name] = current_model

                    # Clear variables before next model
                    current_model = None
                    model_npop = 0

                # Create the model
                current_model = Model(model_line.split(':')[1].strip())

            elif 'TREE:' in model_line:
                # Assign the tree to the model
                current_model.assign_tree(model_line.split(':')[1].strip())

            elif 'NPOP:' in model_line:
                model_npop = int(model_line.split(':')[1].strip())

            elif 'POP:' in model_line:
                pop_name = model_line.split(':')[1].strip()

            elif 'NIND:' in model_line:
                pop_nind = int(model_line.split(':')[1].strip())

            elif 'INDS:' in model_line:
                # Get the list of individuals
                ind_list = [ind.strip() for ind in model_line.split(':')[1].split(',')]
                # Make sure the individual counts match
                if pop_nind != len(ind_list):
                    print ('Number of individuals specified in INDS does not match NIND')
                    sys.exit()

                # Assign the current population
                current_model.assign_pop(pop_name, ind_list)

                # Clear variables before next population/model
                pop_name = ''
                pop_nind = None

    # Make sure the population counts match
    if current_model.npop != model_npop:
        print ('Number of populations found in %s does not match NPOP' % current_model.name)
        sys.exit()

    # Store model
    models_in_file[current_model.name] = current_model

    return models_in_file
