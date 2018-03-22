from model_numerics_per_week import NumericOdeSolutionPerWeek
# nm = NumericOdeSolution()
nm = NumericOdeSolutionPerWeek()
# file = './OutputParameters/parameters-16-Feb-2018-22:55:09.yml'
# nm.load_parameters(file)
nm.parameters_sampling_per_week()
nm.ode_int_solution_per_week()
nm.solution_plot()
r01, r02, r0 = nm.compute_r_zero()
print (r01, r02, r0)
