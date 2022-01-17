"""
Simple example of extracting a frequency spectrum from a SNP data file.
"""
import dadi
import pylab
import matplotlib.pyplot as pyplot
import sys

VCF=sys.argv[1]    #### VCF_file
popmap=sys.argv[2] ### Popmap file containing individual and associated populaiton (tab delimited, similar to stacks popmap)
size=sys.argv[3]   #### The size of the spectrum (symmetric spectrum) 
size=int(size)
nom=sys.argv[4]    #### Output name of the files (plot and spectrum)
pop1=sys.argv[5]   #### population 1
pop2=sys.argv[6]   #### population 2

dd =dadi.Misc.make_data_dict_vcf(VCF,popmap)
fs = dadi.Spectrum.from_data_dict(dd, [pop1,pop2], [size,size], polarized=False)
spectrum_01 = ()
fs.mask[0,1]=True  # mask singleton pop1   comment these lines to keep singleton
fs.mask[1,0]=True  # mask singleton pop2
fs.to_file("spectrum_"+str(size)+"_"+str(size)+nom)
a=dadi.Plotting.plot_single_2d_sfs(fs,vmin=1) 
pyplot.savefig('Plot_FS_dadi_'+str(size)+'_'+str(size)+nom+'.png')
pyplot.savefig('Plot_FS_dadi_'+str(size)+'_'+str(size)+nom+'.svg')
print("FST_measure_by_dadi",fs.Fst())
