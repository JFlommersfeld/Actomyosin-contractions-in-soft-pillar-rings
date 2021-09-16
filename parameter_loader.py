import sys

full_model_parameter_names = ['x_catch', 'x_slip', 'k_off^catch', 'k_off^slip',	'k_on', 'k_on_fil', 'a/kBT', 'N_h',
                              'N_max', 'h*eta_am^eff', 'xi*rho_a(t=0)^2', 'rho_a^max/rho_a(t=0)', 'R0']
density_model_parameter_names = ['h*eta_am^eff', 'xi*NM*rho_a(t=0)^2', 'rho_a^max/rho_a(t=0)', 'R0']


def load_parameters(model_type, parameter_file):
    """
    Loads in all parameter values given in a parameter file.

        Parameters:
            model_type (str):
                sets the type of model, which determines the exact parameter set that is needed. Possible values for
                the parameter model_type are: 'full model', 'density model' and 'constant force model'.
            parameter_file (str):
                the path to the parameter file.

        Returns:
            list: the values of the parameters in the order as specified in the <model_type>_parameter_names lists.
    """
    if model_type == 'full model':
        parameters = [0.]*len(full_model_parameter_names)
        parameters_found = [0]*len(full_model_parameter_names)
        parameter_names = full_model_parameter_names
    elif model_type == 'density model':
        parameters = [0.]*len(density_model_parameter_names)
        parameters_found = [0]*len(density_model_parameter_names)
        parameter_names = density_model_parameter_names
    else:
        print("ERROR: the parameter <model_type> has to be given one of the three values: "
              "'full model' or 'density model'.")
        sys.exit(1)

    f = open(parameter_file)
    for line in f.readlines():
        line_split = line.split()
        try:
            idx = parameter_names.index(line_split[0])
            parameters_found[idx] = 1
            if line_split[0] == 'N_h' or line_split[0] == 'N_max':
                parameters[idx] = int(line_split[2])
            else:
                parameters[idx] = float(line_split[2])
        except ValueError:
            print("WARNING: Parameter {} cannot be interpreted for the model type "
                  "'{}'!".format(line_split[0], model_type))

    f.close()
    if 0 in parameters_found:
        print("ERROR: Not all necessary parameters for the model type '{}' were defined in the given "
              "file!".format(model_type))
        sys.exit(1)

    return parameters
