###  SCRIPT from Adrien Tran Lu Y
### 2021
### Dadi version 2.1 with local research using dual annealing from  https://gitlab.mbb.univ-montp2.fr/khalid/dadi/-/tree/master~
###

GPU_e=False  ### set true is using GPU instead of CPU 

import sys
import os
# set le directory de travail
dir=os.getcwd()
os.chdir(dir)

##### size of the grid
pts_1=20
pts_2=25
pts_3=30



modeldemo=sys.argv[1] # argv is the name of the model called ie 'IM , IM2NG... etc... " 
Fs_file=sys.argv[2] #  name of the FS_file
iteration= sys.argv[3] # integer that give the interation of the run. (i.e 1 , 2...)

if modeldemo=="SI":
	#params= nu1,nu2,Ts
	upper_bound = [100,100,10]
	lower_bound = [0.00001,0.00001, 0.00001]
	p0=[1,1,1]
elif modeldemo=="SIG":
	#params= nu1,nu2,b1,b2,Ts
	upper_bound = [100,100,1000,1000,10]
	lower_bound = [0.00001,0.00001, 0.00001, 0.00001, 0.00001]
	p0=[1,1,1,1,1]
elif modeldemo=="SI2N":
	#params= nu1, nu2, hrf, Ts, Q
	upper_bound = [100,100,1,10,0.95]
	lower_bound = [0.00001,0.00001, 0.01, 0.00001, 0.01]
	p0=[1,1,0.5,1,0.2]
elif modeldemo=="SI2NG":
	#params= nu1, nu2,b1 ,b2, hrf, Ts, Q
	upper_bound = [100,100,100,100,1,10,0.95]
	lower_bound = [0.00001,0.00001,0.00001,0.00001, 0.01, 0.00001, 0.01]
	p0=[1,1,1,1,0.5,1,0.2]
elif modeldemo=="SC":
	#params= nu1,nu2,Ts,Tsc,m12,m21
	upper_bound = [100,100,10,10,10,10]
	lower_bound = [0.00001,0.00001, 0.00001, 0.00001, 0.00001, 0.00001]
	p0 = [2,2,1,1,1,1]
elif modeldemo=="SCG":
	#params= nu1,nu2,b1,b2,Ts,Tsc,m12,m21
	upper_bound = [100,100,100,100,10,10,10,10]
	lower_bound = [0.00001,0.00001,0.1,0.1,0.00001,0.00001,0.00001,0.00001]
	p0 = [2,2,1,1,1,1,1,1]
elif modeldemo=="SC2N":
	#params= nu1, nu2, hrf,Tam,Ts ,m12, m21, Q
	upper_bound = [100,100,1,10,10,10,10,0.99]
	lower_bound = [0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001]
	p0 = [2,2,0.5,1,1,1,1,0.2]
elif modeldemo=="SC2NG":
	#params= nu1, nu2,b1,b2, hrf,Tam,Ts ,m12, m21, Q
	upper_bound = [100,100,100,100,1,10,10,10,10,0.99]
	lower_bound = [0.00001,0.00001,0.1,0.1,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001]
	p0 = [2,2,1,1,0.5,1,1,1,1,0.2]
elif modeldemo=="SC2m":
	#params=nu1,nu2,Ts,Tsc, m12, m21, me12, me21, P
	upper_bound = [100,100,100,100,10,10,10,10,0.99]
	lower_bound = [0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001]
	p0 = [2,2,2,2,1,1,1,1,0.2]
elif modeldemo=="SC2mG":
	#params=nu1, nu2, b1, b2,Ts,Tsc m12, m21, me12, me21, P
	upper_bound = [100,100,100,100,100,100,10,10,10,10,0.99]
	lower_bound = [0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001]
	p0 = [2,2,1,1,2,2,1,1,1,1,0.2]
elif modeldemo=="SC2N2m":
	#params= nu1, nu2, hrf,Ts,Tsc, m12, m21, me12, me21, P, Q
	upper_bound = [100,100,1,100,100,10,10,10,10,0.99,0.99]
	lower_bound = [0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001]
	p0 = [2,2,0.2,2,2,1,1,1,1,0.2,0.2]
elif modeldemo=="SC2N2mG":
	#params= nu1, nu2, b1, b2, hrf,Ts,Tsc,m12, m21, me12, me21, P,Q
	upper_bound = [100,100,100,100,1,100,100,10,10,10,10,0.99,0.99]
	lower_bound = [0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001]
	p0 = [2,2,1,1,0.2,1,1,1,1,1,1,0.2,0.2]
elif modeldemo=="AM":
	#params= nu1,nu2,Tam,Ts,m12,m21
	upper_bound = [100,100,10,10,10,10]
	lower_bound = [0.00001,0.00001, 0.00001, 0.00001, 0.00001, 0.00001]
	p0 = [2,2,1,1,1,1]
elif modeldemo=="AMG":
	#params= nu1,nu2,b1,b2,Tam,Ts,m12,m21
	upper_bound = [100,100,100,100,10,10,10,10]
	lower_bound = [0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001]
	p0 = [2,2,1,1,1,1,1,1]
elif modeldemo=="AM2N":
	#params= nu1, nu2, hrf,Tam,Ts ,m12, m21, Q
	upper_bound = [100,100,1,10,10,10,10,0.99]
	lower_bound = [0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001]
	p0 = [2,2,0.5,1,1,1,1,0.2]
elif modeldemo=="AM2NG":
	#params= nu1, nu2,b1,b2, hrf,Tam,Ts ,m12, m21, Q
	upper_bound = [100,100,100,100,1,10,10,10,10,0.99]
	lower_bound = [0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001]
	p0 = [2,2,1,1,0.5,1,1,1,1,0.2]
elif modeldemo=="AM2m":
	#params=nu1, nu2,Tam,Ts, m12, m21, me12, me21, P
	upper_bound = [100,100,100,100,10,10,10,10,0.99]
	lower_bound = [0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001]
	p0 = [2,2,2,2,1,1,1,1,0.2]
elif modeldemo=="AM2mG":
	#params=nu1, nu2, b1, b2,Tam,Ts, m12, m21, me12, me21, P
	upper_bound = [100,100,100,100,100,100,10,10,10,10,0.99]
	lower_bound = [0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001]
	p0 = [2,2,1,1,2,2,1,1,1,1,0.2]
elif modeldemo=="AM2N2m":
	#params= nu1, nu2, hrf,Tam,Ts, m12, m21, me12, me21, P,Q
	upper_bound = [100,100,1,100,100,10,10,10,10,0.99,0.99]
	lower_bound = [0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001]
	p0 = [2,2,0.2,2,2,1,1,1,1,0.2,0.2]
elif modeldemo=="AM2N2mG":
	#params= nu1, nu2, b1, b2, hrf,Tam,Ts,m12, m21, me12, me21, P,Q
	upper_bound = [100,100,100,100,1,100,100,10,10,10,10,0.99, 0.99]
	lower_bound = [0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001]
	p0 = [2,2,1,1,0.2,2,2,1,1,1,1,0.2, 0.2]
elif modeldemo=="IM":
	#params= nu1,nu2,Tsm,m12,m21
	upper_bound = [100,100,10,10,10]
	lower_bound = [0.00001,0.00001, 0.00001, 0.00001, 0.00001]
	p0 = [2,2,1,1,1]
elif modeldemo=="IMG":
	#params= nu1,nu2,b1,b2,Tsm,m12,m21
	upper_bound = [100,100,100,100,10,10,10]
	lower_bound = [0.00001,0.00001,0.00001, 0.00001, 0.00001, 0.00001, 0.00001]
	p0 = [2,2,1,1,1,1,1]
elif modeldemo=="IM2N":
	#params= nu1, nu2, hrf, Tsm,m12, m21, Q
	upper_bound = [100,100,1,10,10,10,0.99]
	lower_bound = [0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001]
	p0 = [2,2,0.5,1,1,1,0.2]
elif modeldemo=="IM2NG":
	#params= nu1, nu2,b1,b2, hrf, Tsm ,m12, m21, Q
	upper_bound = [100,100,100,100,1,10,10,10,0.99]
	lower_bound = [0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001]
	p0 = [2,2,1,1,0.5,1,1,1,0.2]
elif modeldemo=="IM2m":
	#params=nu1, nu2,Tsm, m12, m21, me12, me21, P
	upper_bound = [100,100,10,10,10,10,10,0.99]
	lower_bound = [0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001]
	p0 = [2,2,1,1,1,1,1,0.2]
elif modeldemo=="IM2mG":
	#params=nu1, nu2, b1, b2, Tsm, m12, m21, me12, me21, P
	upper_bound = [100,100,100,100,100,10,10,10,10,0.99]
	lower_bound = [0.00001,0.00001,0.01,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001]
	p0 = [2,2,1,2,2,1,1,1,1,0.2]
elif modeldemo=="IM2N2m":
	#params= nu1, nu2, hrf, Tsm, m12, m21, me12, me21, P,Q
	upper_bound = [100,100,1,10,10,10,10,10,0.99, 0.99]
	lower_bound = [0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001]
	p0 = [2,2,0.2,2,2,1,1,1,0.2, 0.2]
elif modeldemo=="IM2N2mG":
	#params= nu1, nu2, b1, b2, hrf,Tsm,m12, m21, me12, me21,, P, Q
	upper_bound = [100,100,100,100,1,10,10,10,10,10,0.99,0.99]
	lower_bound = [0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001]
	p0 = [2,2,1,1,0.2,2,1,1,1,1,0.2,0.2]

else:
	print ("error, check code")

print(upper_bound)
print(lower_bound)
print(p0)
print(" ")
from numpy import array
import timeit
import model_no_NA_change
import dadi 

# rewritten model from rougueux et al 2017 with IM2N2M and IM2N2mG in addition and considering only growth phase (G) during migration phase (except for SI)
# Load the data
data = dadi.Spectrum_mod.Spectrum.from_file(Fs_file)
ns = data.sample_sizes

####### load all model
method_to_call=getattr(model_no_NA_change,modeldemo)
print("model called "+modeldemo)

func=method_to_call

# Make the extrapolating version of our demographic model function.
func_ex = dadi.Numerics.make_extrap_log_func(func)

# Perturb our parameters before optimization. This does so by taking each
# parameter a up to a factor of two up or down.
p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)
# Do the optimization. By default we assume that theta is a free parameter,
# since it's trivial to find given the other parameters. If you want to fix
# theta, add a multinom=False to the call.
# The maxiter argument restricts how long the optimizer will run. For real 
# runs, you will want to set this value higher (at least 10), to encourage
# better convergence. You will also want to run optimization several times
# using multiple sets of intial parameters, to be confident you've actually
# found the true maximum likelihood parameters.
print('Beginning optimization ************************************************')


# These are the grid point settings will use for extrapolation.
pts_l = [pts_1,pts_2,pts_3]

########  Benchmark optimize_dual_anneal ######################

maxiterGlobal=5 #maximum global search iterations - in each iteration, it explores twice the number of parameters
accept=1 #parameter for acceptance distribution (lower values makes the probability of acceptance smaller)
visit=1.01 #parameter for visiting distribution (higher values makes the algorithm jump to a more distant region)
Tini=50 #initial temperature
no_local_search=False #If set to True, a Generalized Simulated Annealing will be performed with no local search strategy applied
local_method='L-BFGS-B' #local search method
maxiterLocal=5 #maximum local search iterations
verbose=True  #### set false if you want to remove iteration display in the file during research of the optimal Likelihood

dadi.cuda_enabled(GPU_e)
popt = dadi.MyInference.optimize_dual_anneal(p0=p0, data=data, model_func=func_ex, pts=pts_l,
                                                     lower_bound=lower_bound, upper_bound=upper_bound,
                                                     no_local_search=no_local_search, local_method=local_method, local_maxiter=maxiterLocal,
                                                     maxiter=maxiterGlobal, Tini=Tini, accept=accept, visit=visit, verbose=verbose, full_output=True)   

print('Finishing  optimization ************************************************')



params=popt.x
model = func_ex(params, ns, pts_l)
theta = dadi.Inference.optimal_sfs_scaling(model, data)
ll_opt = dadi.Inference.ll_multinom(model, data)
AIC = 2*len(params)-2*ll_opt
# Print results
#print('Model log-likelihood:', ll_model)
print( 'Optimized parameters', repr(popt))#LogL =>Fmax et xmax => xmin
print( 'Optimized log-likelihood:', ll_opt)
print( 'theta:', theta)
print( 'AIC', AIC)	


Output_path= os.path.join(os.getcwd(),modeldemo)


if not os.path.exists(Output_path):
	os.mkdir(Output_path)


os.chdir(Output_path)

#write_result
line=("\n"+ str(modeldemo) + "\n" + "Model ouptut optimized parameters :" + repr(popt) + "\n" + "Optimized log-likelihood:" + repr(ll_opt) + "\n"+ "Theta:" + repr(theta) +"\n" + "AIC:" +repr(AIC))

aa = open((modeldemo+"_"+iteration+"_.txt"),"w") 
aa.write(line)
aa.close


# Plot a comparison of the resulting fs with the data.
import pylab
pylab.figure(figsize=(8,6))
dadi.Plotting.plot_2d_comp_multinom(model, data, vmin=1, resid_range=3,
                                    pop_ids =('Leg1','Leg2'), show=False)
# Save the figure
pylab.savefig(modeldemo+"_"+iteration+"_.png", dpi=600)
pylab.savefig(modeldemo+"_"+iteration+"_.svg", dpi=600)
pylab.savefig(modeldemo+"_"+iteration+"_.pdf", dpi=600)


