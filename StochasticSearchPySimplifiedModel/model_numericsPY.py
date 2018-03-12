from model_numerics import NumericOdeSolution
nm = NumericOdeSolution()
file = './OutputParameters/parameters-16-Feb-2018-22:55:09.yml'
nm.load_parameters(file)
nm.ode_int_solution()
#nm.parameters_sampling()
nm.solution_plot()
