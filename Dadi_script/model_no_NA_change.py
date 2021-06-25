
import numpy
import dadi

############ Adrien Tran Lu Y ############
############


def IM(params, ns, pts):
    """
    ns = ns
    params = (s,nu1,nu2,T,m12,m21)
    Isolation-with-migration model
    s: Size of pop 1 after split. (Pop 2 has size 1-s.)
    nu1: Final size of pop 1.
    nu2: Final size of pop 2.
    T: Time in the past of split (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1,nu2,Ts,m12,m21 = params

    xx = dadi.Numerics.default_grid(pts)
    
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def IMa(params, ns, pts):
    """
    ns = ns
    params = (s,nu1,nu2,T,m12,m21)

    Isolation-with-migration model with exponential pop growth.

    s: Size of pop 1 after split. (Pop 2 has size 1-s.)
    nu1: Final size of pop 1.
    nu2: Final size of pop 2.
    T: Time in the past of split (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    s,nu1,nu2,T,m12,m21 = params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    nu1_func = lambda t: s * (nu1/s)**(t/T)
    nu2_func = lambda t: (1-s) * (nu2/(1-s))**(t/T)
    phi = dadi.Integration.two_pops(phi, xx, T, nu1_func, nu2_func,
                               m12=m12, m21=m21)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs
	
def IMG(params, ns, pts):
    """
    ns = ns
    params = (s,nu1,nu2,T,m12,m21)
    Isolation-with-migration model with exponential growth
    s: Size of pop 1 after split. (Pop 2 has size 1-s.)
    nu1: Final size of pop 1.
    nu2: Final size of pop 2.
    b1 : growth parameter of pop 1
    b2 : growth parameter of pop 2
    T: Time in the past of split (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1, nu2,b1, b2, Ts, m12, m21= params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)

    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    nu1_func = lambda t: nu1 * b1**(t/Ts)
    nu2_func = lambda t: nu2 * b2**(t/Ts)
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1_func, nu2_func,
                               m12=m12, m21=m21)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs	
		
def IM2N(params, ns, pts):
    nu1, nu2, hrf,Ts,m12, m21, Q = params
    """
    Model with split and migration, heterogenous effective population size (with 2 classes of loci shared by the two populations = Hill-Robertson effects)
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and present (in units of 2*Na generations).
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)

    phinr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to m12 and m21 
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    # oriented
    fsnrO = dadi.Spectrum.from_phi(phinr, ns, (xx,xx))
    # mis-oriented

    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split to hrf*nu1 and hrf*nu2 and the migration rates to m12 and m21 
    philr = dadi.Integration.two_pops(philr, xx, Ts, nu1*hrf, nu2*hrf, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    # oriented 
    fslrO = dadi.Spectrum.from_phi(philr, ns, (xx,xx))
    # mis-oriented

    #### Sum the two spectra in proportion O (and Q) #only use O for unfolded spectrum
    fs= ((1-Q)*fsnrO + Q*fslrO)
    return fs		
	
def IM2NG(params, ns, pts):
    nu1, nu2, b1, b2, hrf,Ts, m12, m21, Q = params

    """
    Model with split with migration, heterogenous effective population size (with 2 classes of loci shared by the two populations = Hill-Robertson effects)
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the splitÂ and present (in units of 2*Na generations).
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """ 
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)
    
    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)

    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phinr = phi
    # We start the population size change after the split independantly in each population and set the migration rates to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ###
    ## Finally, calculate the spectrum.
    # oriented
    fsnrO = dadi.Spectrum.from_phi(phinr, ns, (xx,xx))
    # mis-oriented
	
	
    #### Spectrum of low-recombining regions

    # Now do the divergence event
    philr = phi
    # We start the population size change after the split independantly in each population (bnu{1,2}_func) & integrate the hrf for low-recombining regions and set the migration rates to m12 and m21
    bnu1hrf_func = lambda t: (nu1 * b1**(t/Ts)) * hrf
    bnu2hrf_func = lambda t: (nu2 * b2**(t/Ts)) * hrf
    philr = dadi.Integration.two_pops(philr, xx, Ts, bnu1hrf_func, bnu2hrf_func, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    # oriented
    fslrO = dadi.Spectrum.from_phi(philr, ns, (xx,xx))
    # mis-oriented

    #### Sum the two spectra in proportion O (and Q)
    fs= ((1-Q)*fsnrO + Q*fslrO)
    return fs	
	
def IM2m(params, ns, pts): 
    nu1, nu2,Ts, m12, m21, me12, me21, P = params
    """
    Model with migration during the divergence with two type of migration.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and present (in units of 2*Na generations).
    P: The porportion of the genome evolving neutrally
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    ### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We start the population size change after the split, independently in each population, and set the migration rates to m12 and m21

    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    ###
    ## Finally, calculate the spectrum.
    # oriented
    fsNO = dadi.Spectrum.from_phi(phiN, ns, (xx,xx))
    # mis-oriented

    ### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We start the population size change after the split, independently in each population, and set the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=me12, m21=me21)
    ###
    ## Finally, calculate the spectrum.
    # oriented
    fsIO = dadi.Spectrum.from_phi(phiI, ns, (xx,xx))
    # mis-oriented

    ### Sum the two spectra in proportion O (and P)
    fs = (P*fsNO + (1-P)*fsIO)
    return fs

def IM2mG(params, ns, pts): 
    nu1, nu2, b1, b2, Ts, m12, m21, me12, me21, P = params
    """
    Model with migration during the divergence with two type of migration and population growth
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and present (in units of 2*Na generations).
    P: The porportion of the genome evolving neutrally
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    ### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)

    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phiN= phi
    # We start the population size change after the split, independently in each population, and set the migration rates to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ###
    ## Finally, calculate the spectrum.
    # oriented
    fsNO = dadi.Spectrum.from_phi(phiN, ns, (xx,xx))
    # mis-oriented

    ### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population


    phiI= phi
    # We start the population size change after the split, independently in each population, and set the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, bnu1_func, bnu2_func, m12=me12, m21=me21)
    ###
    ## Finally, calculate the spectrum.
    # oriented
    fsIO = dadi.Spectrum.from_phi(phiI, ns, (xx,xx))
    # mis-oriented

    ### Sum the two spectra in proportion O (and P)
    fs = (P*fsNO + (1-P)*fsIO) 
    return fs

def IM2N2m(params, ns, pts): 
    nu1, nu2, hrf,Ts,m12, m21, me12, me21, P, Q = params
    """
    Model of semi permeability with split and migration with 2 migration rates
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and present(in units of 2*Na generations).
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    P: The proportion of the genome evolving neutrally
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    #### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and m21
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    # oriented
    fsNO = dadi.Spectrum.from_phi(phiN, ns, (xx,xx))
    # mis-oriented

    #### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=me12, m21=me21)
    ###
    ## calculate the spectrum.
    # oriented
    fsIO = dadi.Spectrum.from_phi(phiI, ns, (xx,xx))
    # mis-oriented

    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We keep the population sizes after the split and isolation to nu1 and nu2 and set the migration rates to zero
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    # oriented
    fsnrO = dadi.Spectrum.from_phi(phinr, ns, (xx,xx))
    # mis-oriented
    
    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We keep the population sizes after the split and isolation to nu1 and nu2 and set the migration rate to m12 and m21
    philr = dadi.Integration.two_pops(philr, xx, Ts, nu1*hrf, nu2*hrf, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    # oriented
    fslrO = dadi.Spectrum.from_phi(philr, ns, (xx,xx))
    # mis-oriente

    ### Sum the spectra 
    fs = (Q*fslrO+(1-Q)*fsnrO+P*fsNO+(1-P)*fsIO)
    return fs

def IM2N2mG(params, ns, pts):
    nu1, nu2, b1, b2, hrf,Ts ,m12, m21, me12, me21, P, Q = params
    """
    Model of semi permeability with split and migration with 2 migration rates
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and present(in units of 2*Na generations).
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    P: The proportion of the genome evolving neutrally
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)
    #### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)

    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phiN= phi
    # We start the population sizes change independently in each populations after the split to bnu1 and bnu2 and set the migration rates to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    # oriented
    fsNO = dadi.Spectrum.from_phi(phiN, ns, (xx,xx))
    # mis-oriented

    phiI = phi

    # We keep the population sizes change after the split to bnu1 and bnu2 and set the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, bnu1_func, bnu2_func, m12=me12, m21=me21)
    ###
    ## calculate the spectrum.
    # oriented
    fsIO = dadi.Spectrum.from_phi(phiI, ns, (xx,xx))
    # mis-oriented

    phinr = phi
    # We keep the population sizes change independently in each populations after the isolation to bnu1 and bnu2 and set the migration rates to m12 & m21
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    # oriented
    fsnrO = dadi.Spectrum.from_phi(phinr, ns, (xx,xx))


    philr = phi
    # Now do the divergence event
    # We keep the population sizes after the split and isolation to nu1 and nu2 and set the migration rate to m12 and m21
    bnu1hrf_func = lambda t: (nu1 * b1**(t/Ts)) * hrf
    bnu2hrf_func = lambda t: (nu2 * b2**(t/Ts)) * hrf
    philr = dadi.Integration.two_pops(philr, xx, Ts, bnu1hrf_func, bnu2hrf_func, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    # oriented
    fslrO = dadi.Spectrum.from_phi(philr, ns, (xx,xx))
    # mis-oriented

    ### Sum the spectra 
    fs = (Q*fslrO+(1-Q)*fsnrO+P*fsNO+(1-P)*fsIO)
    return fs

	
	
def SI(params, ns, pts):
    """
    ns = ns
    params = (s,nu1,nu2,Ts,Tam,m12,m21)

.
	s: Size of pop 1 after split. (Pop 2 has size 1-s.)
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    Tam: The scale time between the ancient migration and present.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    nu1,nu2,Ts = params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)

    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2,
                               m12=0, m21=0)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs
		
def SIa(params, ns, pts):
    """
    ns = ns
    params = (s,nu1,nu2,Ts,Tam,m12,m21)

.
	s: Size of pop 1 after split. (Pop 2 has size 1-s.)
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    Tam: The scale time between the ancient migration and present.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    s,nu1,nu2,Ts = params

    xx = dadiNumerics.default_grid(pts)

    phi = dadiPhiManip.phi_1D(xx)
    phi = dadiPhiManip.phi_1D_to_2D(xx, phi)
	
    nu1_func = lambda t: s * (nu1/s)**(t/T)
    nu2_func = lambda t: (1-s) * (nu2/(1-s))**(t/T)
    phi = dadiIntegration.two_pops(phi, xx, Ts, nu1_func, nu2_func,
                               m12=0, m21=0)
    fs = dadiSpectrum.from_phi(phi, ns, (xx,xx))
    return fs
	
def SIG(params, ns, pts):
    """
    ns = ns
    params = (s,nu1,nu2,T,m12,m21)

    Isolation-with-migration model with exponential pop growth.

    s: Size of pop 1 after split. (Pop 2 has size 1-s.)
    nu1: Final size of pop 1.
    nu2: Final size of pop 2.
    T: Time in the past of split (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1, nu2, b1, b2, Ts = params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)

    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phi = dadi.Integration.two_pops(phi, xx, Ts, bnu1_func, bnu2_func,
                               m12=0, m21=0)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs	
		
def SI2N(params, ns, pts):
    nu1, nu2, hrf ,Ts, Q = params
    """
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the splitÂ and present (in units of 2*Na generations).
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)

    phinr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to m12 and m21 
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1, nu2, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    # oriented
    fsnrO = dadi.Spectrum.from_phi(phinr, ns, (xx,xx))
    # mis-oriented

    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)

    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split to hrf*nu1 and hrf*nu2 and the migration rates to m12 and m21 
    philr = dadi.Integration.two_pops(philr, xx, Ts, nu1*hrf, nu2*hrf, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    # oriented
    fslrO = dadi.Spectrum.from_phi(philr, ns, (xx,xx))
    # mis-oriented

    #### Sum the two spectra in proportion O (and Q)
    fs= ((1-Q)*fsnrO + Q*fslrO)
    return fs		
	
def SI2NG(params, ns, pts):
    nu1, nu2, b1, b2, hrf, Ts, Q = params

    """
    Model with split, ancient migration, heterogenous effective population size (with 2 classes of loci shared by the two populations = Hill-Robertson effects)
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the splitÂ and present (in units of 2*Na generations).
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """ 
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)
    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phinr = phi
    # We start the population size change after the split independantly in each population and set the migration rates to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, bnu1_func, bnu2_func, m12=0, m21=0)
    ###
    ## Finally, calculate the spectrum.
    # oriented
    fsnrO = dadi.Spectrum.from_phi(phinr, ns, (xx,xx))
    # mis-oriented

    philr = phi
    # Now do the divergence even
    # We start the population size change after the split independantly in each population (bnu{1,2}_func) & integrate the hrf for low-recombining regions and set the migration rates to m12 and m21
    bnu1hrf_func = lambda t: (nu1 * b1**(t/Ts)) * hrf
    bnu2hrf_func = lambda t: (nu2 * b2**(t/Ts)) * hrf
    philr = dadi.Integration.two_pops(philr, xx, Ts, bnu1hrf_func, bnu2hrf_func, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    # oriented
    fslrO = dadi.Spectrum.from_phi(philr, ns, (xx,xx))
    #### Sum the two spectra in proportion O (and Q)
    fs= ((1-Q)*fsnrO + Q*fslrO) 
    return fs	
	
	
def AM(params, ns, pts):
    """
    ns = ns
    params = (s,nu1,nu2,Ts,Tam,m12,m21)

.
	s: Size of pop 1 after split. (Pop 2 has size 1-s.)
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    Tam: The scale time between the ancient migration and present.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    nu1,nu2,Tam,Ts,m12,m21 = params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)


    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, Tam, nu1, nu2,
                               m12=m12, m21=m21)
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2,
                               m12=0, m21=0)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def AMa(params, ns, pts):
    """
    ns = ns
    params = (s,nu1,nu2,Ts,Tam,m12,m21)

.
	s: Size of pop 1 after split. (Pop 2 has size 1-s.)
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    Tam: The scale time between the ancient migration and present.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    s,nu1,nu2,Ts,Tam,m12,m21 = params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    nu1_func = lambda t: s * (nu1/s)**(t/T)
    nu2_func = lambda t: (1-s) * (nu2/(1-s))**(t/T)
    phi = dadi.Integration.two_pops(phi, xx, Tam, nu1_func, nu2_func,
                               m12=m12, m21=m21)
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1_func, nu2_func,
                               m12=0, m21=0)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs
	
def AMG(params, ns, pts):
    """
    ns = ns
    params = (s,nu1,nu2,T,m12,m21)

    Isolation-with-migration model with exponential pop growth.

    s: Size of pop 1 after split. (Pop 2 has size 1-s.)
    nu1: Final size of pop 1.
    nu2: Final size of pop 2.
    T: Time in the past of split (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1, nu2, b1, b2, Ts,Tam, m12, m21= params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)


    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    nu1_func = lambda t: nu1 * b1**(t/Tam)
    nu2_func = lambda t: nu2 * b2**(t/Tam)
    phi = dadi.Integration.two_pops(phi, xx, Tam, nu1_func, nu2_func,
                               m12=m12, m21=m21)
    nu1_b = nu1 * b1
    nu2_b = nu2 * b2   
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1_b, nu2_b,
                               m12=0, m21=0)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs	
		
def AM2N(params, ns, pts):
    nu1, nu2, hrf, Ts, Tam,m12, m21, Q = params

    """
    Model with split, ancient migration, heterogenous effective population size (with 2 classes of loci shared by the two populations = Hill-Robertson effects)

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the splitÂ and present (in units of 2*Na generations).
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phinr = phi
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to m12 and m21 
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    phinr = dadi.Integration.two_pops(phinr, xx, Tam, nu1, nu2, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    # oriented
    fsnrO = dadi.Spectrum.from_phi(phinr, ns, (xx,xx))

    philr = phi
    # Now do the divergence even
    # We set the population sizes after the split to hrf*nu1 and hrf*nu2 and the migration rates to m12 and m21 
    philr = dadi.Integration.two_pops(philr, xx, Ts, nu1*hrf, nu2*hrf, m12=m12, m21=m21)
    philr = dadi.Integration.two_pops(philr, xx, Ts, nu1*hrf, nu2*hrf, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    # oriented
    fslrO = dadi.Spectrum.from_phi(philr, ns, (xx,xx))

    #### Sum the two spectra in proportion O (and Q)
    fs= ((1-Q)*fsnrO + Q*fslrO)
    return fs		
	
def AM2NG(params, ns, pts):
    nu1, nu2, b1, b2, hrf, Ts, Tam, m12, m21, Q = params

    """
    Model with split, ancient migration, heterogenous effective population size (with 2 classes of loci shared by the two populations = Hill-Robertson effects)
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the splitÂ and present (in units of 2*Na generations).
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """ 
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)
    
    phi = dadi.PhiManip.phi_1D(xx)


    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phinr=phi
    # We start the population size change after the split independantly in each population and set the migration rates to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/Tam)
    bnu2_func = lambda t: nu2 * b2**(t/Tam)
    phinr = dadi.Integration.two_pops(phinr, xx, Tam, bnu1_func, bnu2_func, m12=m12, m21=m21)
    nu1_b = nu1 * b1
    nu2_b = nu2 * b2   
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1_b, nu2_b,m12=0, m21=0)

    ###
    ## Finally, calculate the spectrum.
    # oriented
    fsnrO = dadi.Spectrum.from_phi(phinr, ns, (xx,xx))

    philr = phi

    # Now do the divergence event
    # We start the population size change after the split independantly in each population (bnu{1,2}_func) & integrate the hrf for low-recombining regions and set the migration rates to m12 and m21
    bnu1hrf_func = lambda t: (nu1 * b1**(t/Tam)) * hrf
    bnu2hrf_func = lambda t: (nu2 * b2**(t/Tam)) * hrf
    philr = dadi.Integration.two_pops(philr, xx, Tam, bnu1hrf_func, bnu2hrf_func, m12=m12, m21=m21)
    nu1_b = nu1 * b1*hrf
    nu2_b = nu2 * b2*hrf
    philr = dadi.Integration.two_pops(philr, xx, Ts, nu1_b, nu2_b,m12=0, m21=0)
    ###
    ## calculate the spectrum.
    # oriented
    fslrO = dadi.Spectrum.from_phi(philr, ns, (xx,xx))
    # mis-oriented

    #### Sum the two spectra in proportion O (and Q)
    fs= ((1-Q)*fsnrO + Q*fslrO) 
    return fs	
	
def AM2m(params, ns, pts): 
    nu1, nu2,Ts,Tam, m12, m21, me12, me21, P = params
    """
    Model with migration during the divergence with two type of migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and present (in units of 2*Na generations).
    P: The porportion of the genome evolving neutrally
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)


    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phiN=phi    # We start the population size change after the split, independently in each population, and set the migration rates to m12 and m21

    phiN = dadi.Integration.two_pops(phiN, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=0, m21=0)
    ###
    ## Finally, calculate the spectrum.
    # oriented
    fsNO = dadi.Spectrum.from_phi(phiN, ns, (xx,xx))


    ### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population

    # Now do the divergence event
    phiI = phi
    # We start the population size change after the split, independently in each population, and set the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Tam, nu1, nu2, m12=me12, m21=me21)
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=0, m21=0)
    ###
    ## Finally, calculate the spectrum.
    # oriented
    fsIO = dadi.Spectrum.from_phi(phiI, ns, (xx,xx))
    # mis-oriented

    ### Sum the two spectra in proportion O (and P)
    fs = (P*fsNO + (1-P)*fsIO)
    return fs

def AM2mG(params, ns, pts): 
    nu1, nu2, b1, b2, Ts, Tam, m12, m21, me12, me21, P = params
    """
    Model with migration during the divergence with two type of migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and present (in units of 2*Na generations).
    P: The porportion of the genome evolving neutrally
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    ### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)


    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phiN=phi
    # We start the population size change after the split, independently in each population, and set the migration rates to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/Tam)
    bnu2_func = lambda t: nu2 * b2**(t/Tam)
    phiN = dadi.Integration.two_pops(phiN, xx, Tam, bnu1_func, bnu2_func, m12=m12, m21=m21)
    nu1_b = nu1 * b1
    nu2_b = nu2 * b2   
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1_b, nu2_b,m12=0, m21=0)
    ###
    ## Finally, calculate the spectrum.
    # oriented
    fsNO = dadi.Spectrum.from_phi(phiN, ns, (xx,xx))



    phiI = phi
    # We start the population size change after the split, independently in each population, and set the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Tam, bnu1_func, bnu2_func, m12=me12, m21=me21)
    nu1_b = nu1 * b1
    nu2_b = nu2 * b2   
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1_b, nu2_b,m12=0, m21=0)
    ###
    ## Finally, calculate the spectrum.
    # oriented
    fsIO = dadi.Spectrum.from_phi(phiI, ns, (xx,xx))
    # mis-oriented

    ### Sum the two spectra in proportion O (and P)
    fs = (P*fsNO + (1-P)*fsIO)
    return fs

def AM2N2m(params, ns, pts): 
    nu1, nu2, hrf, Ts, Tam, m12, m21, me12, me21, P, Q = params
    """
    Model of semi permeability with split, complete isolation, followed by secondary contact with 2 migration rates

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    P: The proportion of the genome evolving neutrally
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)


    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phiN=phi
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and m21
    phiN = dadi.Integration.two_pops(phiN, xx, Tam, nu1, nu2, m12=m12, m21=m21) 
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    # oriented
    fsNO = dadi.Spectrum.from_phi(phiN, ns, (xx,xx))


    #### Calculate the genomic island spectrum

    phiI = phi
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Tam, nu1, nu2, m12=me12, m21=me21)
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    # oriented
    fsIO = dadi.Spectrum.from_phi(phiI, ns, (xx,xx))

    phinr = phi

    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We keep the population sizes after the split and isolation to nu1 and nu2 and set the migration rates to zero
    phinr = dadi.Integration.two_pops(phinr, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1, nu2, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    # oriented
    fsnrO = dadi.Spectrum.from_phi(phinr, ns, (xx,xx))
    # mis-oriented

    philr = phi
    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We keep the population sizes after the split and isolation to nu1 and nu2 and set the migration rate to m12 and m21
    philr = dadi.Integration.two_pops(philr, xx, Tam, nu1*hrf, nu2*hrf, m12=m12, m21=m21)
    philr = dadi.Integration.two_pops(philr, xx, Ts, nu1*hrf, nu2*hrf, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    # oriented
    fslrO = dadi.Spectrum.from_phi(philr, ns, (xx,xx))


    ### Sum the spectra 
    fs = (Q*fslrO+(1-Q)*fsnrO+P*fsNO+(1-P)*fsIO)
    return fs

def AM2N2mG(params, ns, pts):
    nu1, nu2, b1, b2, hrf, Ts, Tam,m12, m21, me12, me21, P, Q = params
    """
    Model of semi permeability with split, complete isolation, followed by secondary contact with 2 migration rates

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    P: The proportion of the genome evolving neutrally
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)
    #### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)


    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phiN=phi
    # We start the population sizes change independently in each populations after the split to bnu1 and bnu2 and set the migration rates to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/Tam)
    bnu2_func = lambda t: nu2 * b2**(t/Tam)
    phiN = dadi.Integration.two_pops(phiN, xx, Tam, bnu1_func, bnu2_func, m12=m12, m21=m21)
    nu1_b = nu1 * b1
    nu2_b = nu2 * b2   
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1_b, nu2_b,m12=0, m21=0)
    ###
    ## calculate the spectrum.
    # oriented
    fsNO = dadi.Spectrum.from_phi(phiN, ns, (xx,xx))
    # mis-oriented
    phiI = phi
    # Now do the divergence event

    # We keep the population sizes change after the split to bnu1 and bnu2 and set the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Tam, bnu1_func, bnu2_func, m12=me12, m21=me21)
    nu1_b = nu1 * b1
    nu2_b = nu2 * b2   
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1_b, nu2_b,m12=0, m21=0)
    ###
    ## calculate the spectrum.
    # oriented
    fsIO = dadi.Spectrum.from_phi(phiI, ns, (xx,xx))
    # mis-oriented

    phinr = phi

    # Now do the divergence event
    # We keep the population sizes change independently in each populations after the isolation to bnu1 and bnu2 and set the migration rates to m12 & m21
    phinr = dadi.Integration.two_pops(phinr, xx, Tam, bnu1_func, bnu2_func, m12=m12, m21=m21)
    nu1_b = nu1 * b1
    nu2_b = nu2 * b2   
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1_b, nu2_b,m12=0, m21=0)
    ###
    ## calculate the spectrum.
    # oriented
    fsnrO = dadi.Spectrum.from_phi(phinr, ns, (xx,xx))
    # mis-oriented^
    philr = phi

    # Now do the divergence event
    # We keep the population sizes after the split and isolation to nu1 and nu2 and set the migration rate to m12 and m21
    bnu1hrf_func = lambda t: (nu1 * b1**(t/Tam)) * hrf
    bnu2hrf_func = lambda t: (nu2 * b2**(t/Tam)) * hrf
    philr = dadi.Integration.two_pops(philr, xx, Tam, bnu1hrf_func, bnu2hrf_func, m12=m12, m21=m21)
    nu1_b = nu1 * b1*hrf
    nu2_b = nu2 * b2*hrf   
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1_b, nu2_b,m12=0, m21=0)
    ###
    ## calculate the spectrum.
    # oriented
    fslrO = dadi.Spectrum.from_phi(philr, ns, (xx,xx))
    # mis-oriented


    ### Sum the spectra 
    fs = (Q*fslrO+(1-Q)*fsnrO+P*fsNO+(1-P)*fsIO)
    return fs

	
def SC(params, ns, pts):
    """
    ns = ns
    params = (s,nu1,nu2,Ts,Tsc,m12,m21)

.
	s: Size of pop 1 after split. (Pop 2 has size 1-s.)
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    Tam: The scale time between the ancient migration and present.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    nu1,nu2,Ts,Tsc,m12,m21 = params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)


    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2,
                               m12=0, m21=0)
    phi = dadi.Integration.two_pops(phi, xx, Tsc, nu1, nu2,
                               m12=m12, m21=m21)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs	
	
def SCa(params, ns, pts):
    """
    ns = ns
    params = (s,nu1,nu2,Ts,Tsc,m12,m21)

.
	s: Size of pop 1 after split. (Pop 2 has size 1-s.)
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    Tam: The scale time between the ancient migration and present.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    nu1,nu2,Ts,Tsc,m12,m21 = params

    xx = dadi.Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    nu1_func = lambda t: s * (nu1/s)**(t/Ts)
    nu2_func = lambda t: (1-s) * (nu2/(1-s))**(t/Ts)
    phi = Integration.two_pops(phi, xx, Ts, nu1_func, nu2_func,
                               m12=0, m21=0)
    phi = Integration.two_pops(phi, xx, Tsc, nu1_func, nu2_func,
                               m12=m12, m21=m21)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs	
	
def SCG(params, ns, pts):
    """
    ns = ns
    params = (s,nu1,nu2,T,m12,m21)

    Isolation-with-migration model with exponential pop growth.

    s: Size of pop 1 after split. (Pop 2 has size 1-s.)
    nu1: Final size of pop 1.
    nu2: Final size of pop 2.
    T: Time in the past of split (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1, nu2, b1, b2, Ts,Tsc, m12, m21= params
    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)

    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    nu1_func = lambda t: nu1 * b1**(t/Tsc)
    nu2_func = lambda t: nu2 * b2**(t/Tsc)
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2,
                               m12=0, m21=0)
    phi = dadi.Integration.two_pops(phi, xx, Tsc, nu1_func, nu2_func,
                               m12=m12, m21=m21)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs	
		
def SC2N(params, ns, pts):
    nu1, nu2, hrf, Ts, Tsc,m12, m21, Q = params

    """
    Model with split, ancient migration, heterogenous effective population size (with 2 classes of loci shared by the two populations = Hill-Robertson effects)

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the splitÂ and present (in units of 2*Na generations).
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)


    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phinr=phi
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to m12 and m21 
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1, nu2, m12=0, m21=0)
    phinr = dadi.Integration.two_pops(phinr, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    # oriented
    fsnrO = dadi.Spectrum.from_phi(phinr, ns, (xx,xx))
    # mis-oriented



    # Now do the divergence event
    philr = phi
    # We set the population sizes after the split to hrf*nu1 and hrf*nu2 and the migration rates to m12 and m21 
    philr = dadi.Integration.two_pops(philr, xx, Ts, nu1*hrf, nu2*hrf, m12=0, m21=0)
    philr = dadi.Integration.two_pops(philr, xx, Tsc, nu1*hrf, nu2*hrf, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    # oriented
    fslrO = dadi.Spectrum.from_phi(philr, ns, (xx,xx))
    # mis-oriented

    #### Sum the two spectra in proportion O (and Q)
    fs= ((1-Q)*fsnrO + Q*fslrO)
    return fs		
	
def SC2NG(params, ns, pts):
    nu1, nu2, b1, b2, hrf, Ts, Tsc, m12, m21, Q = params

    """
    Model with split, ancient migration, heterogenous effective population size (with 2 classes of loci shared by the two populations = Hill-Robertson effects)
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the splitÂ and present (in units of 2*Na generations).
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """ 
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)
    
    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)

    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phinr=phi
    # We start the population size change after the split independantly in each population and set the migration rates to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/Tsc)
    bnu2_func = lambda t: nu2 * b2**(t/Tsc)
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1, nu2, m12=0, m21=0)
    phinr = dadi.Integration.two_pops(phinr, xx, Tsc, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ###
    ## Finally, calculate the spectrum.
    # oriented
    fsnrO = dadi.Spectrum.from_phi(phinr, ns, (xx,xx))
    # mis-oriented
    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population

    philr = phi
    # We start the population size change after the split independantly in each population (bnu{1,2}_func) & integrate the hrf for low-recombining regions and set the migration rates to m12 and m21
    bnu1hrf_func = lambda t: (nu1 * b1**(t/Tsc)) * hrf
    bnu2hrf_func = lambda t: (nu2 * b2**(t/Tsc)) * hrf
    philr = dadi.Integration.two_pops(philr, xx, Ts,  nu1*hrf, nu2*hrf, m12=0, m21=0)
    philr = dadi.Integration.two_pops(philr, xx, Tsc, bnu1hrf_func, bnu2hrf_func, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    # oriented
    fslrO = dadi.Spectrum.from_phi(philr, ns, (xx,xx))
    # mis-oriented

    #### Sum the two spectra in proportion O (and Q)
    fs= ((1-Q)*fsnrO + Q*fslrO)
    return fs	
	
def SC2m(params, ns, pts): 
    nu1, nu2,Ts,Tsc, m12, m21, me12, me21, P = params
    """
    Model with migration during the divergence with two type of migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and present (in units of 2*Na generations).
    P: The porportion of the genome evolving neutrally
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)


    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phiN=phi
    # We start the population size change after the split, independently in each population, and set the migration rates to m12 and m21

    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=0,m21=0)
    phiN = dadi.Integration.two_pops(phiN, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    ###
    ## Finally, calculate the spectrum.
    # oriented
    fsNO = dadi.Spectrum.from_phi(phiN, ns, (xx,xx))
    # mis-oriented

    # Now do the divergence event
    phiI = phi
    # We start the population size change after the split, independently in each population, and set the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=0, m21=0)
    phiI = dadi.Integration.two_pops(phiI, xx, Tsc, nu1, nu2,  m12=me12, m21=me21)
    ###
    ## Finally, calculate the spectrum.
    # oriented
    fsIO = dadi.Spectrum.from_phi(phiI, ns, (xx,xx))
    # mis-oriented

    ### Sum the two spectra in proportion O (and P)
    fs = (P*fsNO + (1-P)*fsIO)
    return fs

def SC2mG(params, ns, pts): 
    nu1, nu2, b1, b2, Ts, Tsc, m12, m21, me12, me21, P = params
    """
    Model with migration during the divergence with two type of migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and present (in units of 2*Na generations).
    P: The porportion of the genome evolving neutrally
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)


    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phiN=phi
    # We start the population size change after the split, independently in each population, and set the migration rates to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/Tsc)
    bnu2_func = lambda t: nu2 * b2**(t/Tsc)
    phiN = dadi.Integration.two_pops(phiN, xx, Ts,  nu1, nu2, m12=0, m21=0)
    phiN = dadi.Integration.two_pops(phiN, xx, Tsc, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ###
    ## Finally, calculate the spectrum.
    # oriented
    fsNO = dadi.Spectrum.from_phi(phiN, ns, (xx,xx))


    ### Calculate the genomic island spectrum

    # Now do the divergence event
    phiI = phi
    # We start the population size change after the split, independently in each population, and set the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Ts,  nu1, nu2, m12=0, m21=0)
    phiI = dadi.Integration.two_pops(phiI, xx, Tsc, bnu1_func, bnu2_func, m12=me12, m21=me21)
    ###
    ## Finally, calculate the spectrum.
    # oriented
    fsIO = dadi.Spectrum.from_phi(phiI, ns, (xx,xx))


    ### Sum the two spectra in proportion O (and P)
    fs = (P*fsNO + (1-P)*fsIO)
    return fs

def SC2N2m(params, ns, pts): 
    nu1, nu2, hrf, Ts, Tsc, m12, m21, me12, me21, P, Q = params
    """
    Model of semi permeability with split, complete isolation, followed by secondary contact with 2 migration rates

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    P: The proportion of the genome evolving neutrally
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    #### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phiN = phi
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and m21
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=0, m21=0) 
    phiN = dadi.Integration.two_pops(phiN, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    # oriented
    fsNO = dadi.Spectrum.from_phi(phiN, ns, (xx,xx))

    #### Calculate the genomic island spectrum
    phiI = phi
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=0, m21=0)
    phiI = dadi.Integration.two_pops(phiI, xx, Tsc, nu1, nu2, m12=me12, m21=me21)
    ###
    ## calculate the spectrum.
    # oriented
    fsIO = dadi.Spectrum.from_phi(phiI, ns, (xx,xx))

    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population

    # We keep the population sizes after the split and isolation to nu1 and nu2 and set the migration rates to zero
    phinr = phi
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1, nu2, m12=0, m21=0)
    phinr = dadi.Integration.two_pops(phinr, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    # oriented
    fsnrO = dadi.Spectrum.from_phi(phinr, ns, (xx,xx))
    # mis-oriented
    

    # Now do the divergence event
    philr = phi
    # We keep the population sizes after the split and isolation to nu1 and nu2 and set the migration rate to m12 and m21
    philr = dadi.Integration.two_pops(philr, xx, Ts, nu1*hrf, nu2*hrf, m12=0, m21=0)
    philr = dadi.Integration.two_pops(philr, xx, Tsc, nu1*hrf, nu2*hrf, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    # oriented
    fslrO = dadi.Spectrum.from_phi(philr, ns, (xx,xx))
    # mis-oriented
    

    ### Sum the spectra 
    fs = (Q*fslrO+(1-Q)*fsnrO+P*fsNO+(1-P)*fsIO)
    return fs

def SC2N2mG(params, ns, pts):
    nu1, nu2, b1, b2, hrf, Ts, Tsc,m12, m21, me12, me21, P, Q = params
    """
    Model of semi permeability with split, complete isolation, followed by secondary contact with 2 migration rates

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    P: The proportion of the genome evolving neutrally
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)
    #### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)

    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phiN=phi
    # We start the population sizes change independently in each populations after the split to bnu1 and bnu2 and set the migration rates to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/Tsc)
    bnu2_func = lambda t: nu2 * b2**(t/Tsc)
    phiN = dadi.Integration.two_pops(phiN, xx, Ts,  nu1, nu2,  m12=0, m21=0)
    phiN = dadi.Integration.two_pops(phiN, xx, Tsc, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    # oriented
    fsNO = dadi.Spectrum.from_phi(phiN, ns, (xx,xx))
    # mis-oriented


    #### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population

    phiI = phi

    # We keep the population sizes change after the split to bnu1 and bnu2 and set the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Ts,  nu1, nu2, m12=0, m21=0)
    phiI = dadi.Integration.two_pops(phiI, xx, Tsc, bnu1_func, bnu2_func, m12=me12, m21=me21)
    ###
    ## calculate the spectrum.
    # oriented
    fsIO = dadi.Spectrum.from_phi(phiI, ns, (xx,xx))
    # mis-oriented

    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population

    # Now do the divergence event
    phinr = phi
    # We keep the population sizes change independently in each populations after the isolation to bnu1 and bnu2 and set the migration rates to m12 & m21
    phinr = dadi.Integration.two_pops(phinr, xx, Ts,  nu1, nu2,m12=0, m21=0)
    phinr = dadi.Integration.two_pops(phinr, xx, Tsc, bnu1_func, bnu2_func,m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    # oriented
    fsnrO = dadi.Spectrum.from_phi(phinr, ns, (xx,xx))
    # mis-oriented


    # Now do the divergence event
    philr = phi
    # We keep the population sizes after the split and isolation to nu1 and nu2 and set the migration rate to m12 and m21
    bnu1hrf_func = lambda t: (nu1 * b1**(t/Tsc)) * hrf
    bnu2hrf_func = lambda t: (nu2 * b2**(t/Tsc)) * hrf
    philr = dadi.Integration.two_pops(philr, xx, Ts,  nu1*hrf, nu2*hrf, m12=0, m21=0)
    philr = dadi.Integration.two_pops(philr, xx, Tsc, bnu1hrf_func, bnu2hrf_func, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    # oriented
    fslrO = dadi.Spectrum.from_phi(philr, ns, (xx,xx))
    # mis-oriented

    ### Sum the spectra 
    fs = (Q*fslrO)+((1-Q)*fsnrO)+(P*fsNO)+((1-P)*fsIO)
    return fs		
	
	
	
	
	
	
	
	
