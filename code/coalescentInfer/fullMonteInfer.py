import numpy as np
import os
import re
import pandas as pd
import glob
import scipy.stats as stats
import pickle
import mpmath
import warnings
from numba import njit, set_num_threads
from numba import jit
from numba import prange
from bisect import bisect_left

mpmath.mp.dps = 1000

@njit()
def logPriorAlpha(alpha, mean=0, std=0.1):
    # half normal distribution over alpha
    return np.log(1/(std*np.sqrt(np.pi/2))*np.exp(-(alpha-mean)**2/(2*std**2)))

@njit()
def logPriorNu(nu, mean=0, std=0.01):
    # half normal distribution over nu
    return np.log(1/(std*np.sqrt(np.pi/2))*np.exp(-(nu-mean)**2/(2*std**2)))

@njit()
def logPriorTauAlpha(tau_alpha, mean=0, std=1):
    # half normal distribution over tau_alpha
    return np.log(1/(std*np.sqrt(np.pi/2))*np.exp(-(tau_alpha-mean)**2/(2*std**2)))

@njit()
def logPriorTauNu(tau_nu, mean=0, std=1):
    # half normal distribution over tau_nu 
    return np.log(1/(std*np.sqrt(np.pi/2))*np.exp(-(tau_nu-mean)**2/(2*std**2)))

@njit(parallel=True)
def logHyperPrior(alphas, nus, tau_alphas, tau_nus):
    
    logHyperPrior = np.zeros((alphas.shape[0], nus.shape[0], tau_alphas.shape[0], tau_nus.shape[0]))
    for alpha_index, alpha in enumerate(alphas):
        for nu_index, nu in enumerate(nus):
            for tau_alpha_index, tau_alpha in enumerate(tau_alphas):
                for tau_nu_index, tau_nu in enumerate(tau_nus):
                    logHyperPrior[alpha_index, nu_index, tau_alpha_index, tau_nu_index] = logPriorAlpha(alpha)+logPriorNu(nu)+logPriorTauAlpha(tau_alpha)+logPriorTauNu(tau_nu)
    return logHyperPrior

@njit(parallel=True)
def logPopPriorAlphas(alphas, tau_alphas, times, age):
    # Makes the log priors for a given hyperprior value of alpha, tau_alpha and age of donor
    
    logPopPriorAlpha = np.zeros((alphas.shape[0], tau_alphas.shape[0], times.shape[0]))
    for alpha_index, alpha in enumerate(alphas):
            for tau_alpha_index, tau_alpha in enumerate(tau_alphas):
                logPopPriorAlpha[alpha_index, tau_alpha_index, :] = np.log(age/(times*tau_alpha*np.sqrt(2*np.pi))) -(np.log(times/age)-np.log(alpha))**2/(2*tau_alpha**2)
    return logPopPriorAlpha

@njit(parallel=True)
def logPopPriorThetas(nus, tau_nus, thetas, bases):
    # Makes the log priors for a given hyperprior value of nu, tau_nu and number of bases 
    
    logPopPriorTheta = np.zeros((nus.shape[0], tau_nus.shape[0], thetas.shape[0]))
    for nu_index, nu in enumerate(nus):
            for tau_nu_index, tau_nu in enumerate(tau_nus):
                logPopPriorTheta[nu_index, tau_nu_index, :] = np.log(bases/(thetas*tau_nu*np.sqrt(2*np.pi))) - (np.log(thetas/bases)-np.log(nu))**2/(2*tau_nu**2)
    return logPopPriorTheta


def hyperPriorDiff(hyperPriors, nus):
    diffs = np.diff(nus)
    diffs = np.append(diffs,diffs[-1])
    bigDiffs = np.tile(diffs, (hyperPriors.shape[0], hyperPriors.shape[2], hyperPriors.shape[3], 1))
    bigDiffs = np.swapaxes(np.swapaxes(bigDiffs, 2,3), 1,2)
    hyperPriors = bigDiffs* hyperPriors
    return hyperPriors




def logPopMargSum_thetsmal(donorLikely, alphaPopPrior, thetaPopPrior, thetas):
    # This is going to sum over time and mutation rate for a given donor. The output is the logarithm of this sum, with a shape equal to the hyperpriors.
    
    
    donorMarger = np.zeros((alphaPopPrior.shape[0], thetaPopPrior.shape[0], alphaPopPrior.shape[1], thetaPopPrior.shape[1]))
    diffThetaPrior = popPriorDiffReady(thetaPopPrior, thetas)
    
    # Tile this alpha array to include a theta value
    donorAlpha_ThetasAdded = np.tile(alphaPopPrior, (thetaPopPrior.shape[2],1,1,1))
    donorAlpha_ThetasAdded = np.swapaxes(np.swapaxes(np.swapaxes(donorAlpha_ThetasAdded, 0,1), 1,2), 2,3)
    
    # Tile the likelihood to have alpha and tau_alpha
    donorLikelyAlphas = np.tile(donorLikely, (alphaPopPrior.shape[0],alphaPopPrior.shape[1],1,1))
    
    # An array where [i,j,k] has i=alpha, j=tau_alpha, k=theta
    need_nus = (donorAlpha_ThetasAdded * donorLikelyAlphas).sum(axis=2)
    
    # Sum over one axis without constructing dataframe!
    for theta_index in range(need_nus.shape[2]):
    
        # Tile to include nu and tau_nu
        noNeedNus = np.tile(need_nus[:,:,theta_index], (thetaPopPrior.shape[0], thetaPopPrior.shape[1],1,1))
        noNeedNus = np.swapaxes(np.swapaxes(np.swapaxes(noNeedNus, 0,2), 1,2), 2,3)
    
        # Tile to include alpha and tau_alpha
        donorTheta_AlphasAdded = np.tile(thetaPopPrior[:,:,theta_index], (alphaPopPrior.shape[0], alphaPopPrior.shape[1],1,1))
        donorTheta_AlphasAdded = np.swapaxes(donorTheta_AlphasAdded, 1, 2)
    
        # Array [i,j,k,l] where i=alpha, j=nu, k=tau_alpha, l=tau_nu
        donorMarger += (noNeedNus * donorTheta_AlphasAdded)
        
        
    return np.log(donorMarger)

def hyperMargSum_small(nonDonorMarger, hyperPriors, donorAlphaPrior, donorThetaPrior, donorLikelihood, nus):
    
    donor_marger = np.zeros(donorLikelihood.shape)
    hyperDiffPriors = hyperPriorDiff(hyperPriors, nus)
    sumer = hyperDiffPriors*nonDonorMarger

    alphaGone = np.zeros((hyperPriors.shape[1], hyperPriors.shape[3], donorLikelihood.shape[0]))
    
    for alpha_index in range(hyperPriors.shape[0]):
    
        # This is going to look complicated but we are just tiling the dataframe to add in axes for nu and tau_nu
        donorAlpha_NusAdded = np.tile(donorAlphaPrior[alpha_index], (donorThetaPrior.shape[0],donorThetaPrior.shape[1],1,1) )
        donorAlpha_NusAdded = np.swapaxes(donorAlpha_NusAdded, 1,2)

        # In this case we are adding the axis for time (specific to donor)
        donorTiled = np.tile(sumer[alpha_index], (donorLikelihood.shape[0],1,1,1) )
        donorTiled = np.swapaxes(np.swapaxes(np.swapaxes(donorTiled, 0,1), 1,2), 2,3)

        # Now we are going to sum over tau_alpha - so this is nu, tau_nu, time
        alphaGone += (donorTiled * donorAlpha_NusAdded).sum(axis=1)
    
    # Then we add in an axis for theta
    alphaGone = np.tile(alphaGone, (donorLikelihood.shape[1],1,1,1))
    alphaGone = np.swapaxes(np.swapaxes(np.swapaxes(alphaGone, 0,1), 1,2), 2,3)
    
    # This is adding in a time axis for the theta prior
    donorTheta_timeAdded = np.tile(donorThetaPrior, (donorLikelihood.shape[0],1,1,1))
    donorTheta_timeAdded = np.swapaxes(np.swapaxes(donorTheta_timeAdded, 0,1), 1,2)
    
    final_sum = (donorTheta_timeAdded*alphaGone).sum(axis=0).sum(axis=0)
    
    donor_marger = final_sum * donorLikelihood
    
    return donor_marger

def popPriorDiffReady(thetaPrior, thetas):
    diffs = np.diff(thetas)
    diffs = np.append(diffs,diffs[-1])
    
    bigDiffs = np.tile(diffs, (thetaPrior.shape[0], thetaPrior.shape[1], 1))
            
    return thetaPrior*bigDiffs

def createLogDonorLikely_newest(alphaPrior, thetaPrior, donorLikely, thetas):
    #This is going to produce the donor specific hyper log likelhoods integrated over time and theta.
    
    for donor in donorLikely:
        print(donor)
        donorHyperLikely = logPopMargSum_thetsmal(donorLikely[donor], alphaPrior[donor], thetaPrior[donor], thetas)
        print('Saving...')
        with open(f'{save_folder}/logDonorHyper_{donor}.pkl', 'wb') as handle:
            pickle.dump(donorHyperLikely, handle, protocol=pickle.HIGHEST_PROTOCOL)
   
    return

def createDonorMarginal_newest(hyperPriors, alphaPrior, thetaPrior, donorLikely, donor_id, nus):
    #This is going to produce the donor specific marginal
    
    donorMarginal = np.zeros(donorLikely[donor_id].shape)
    
    nonDonorMarger = np.zeros((hyperPriors.shape[0], hyperPriors.shape[1], hyperPriors.shape[2], hyperPriors.shape[3]))
    for donor in donorLikely:
        if donor == donor_id:
            pass
        else:
            with open(f'{save_folder}/logDonorHyper_{donor}.pkl', 'rb') as handle:
                donorLogHyperLikely = pickle.load(handle)
            nonDonorMarger += donorLogHyperLikely
    
    nonDonorMarger = np.exp(nonDonorMarger)
    alpha_donor, theta_donor, likely_donor = alphaPrior[donor_id], thetaPrior[donor_id], donorLikely[donor_id]
    print('hyperMargSum')
    donorMarginal = hyperMargSum_small(nonDonorMarger, hyperPriors, alpha_donor, theta_donor, likely_donor, nus)
    return donorMarginal

def createHyperMarginal_newest(hyperPriors, donorLikely):
    #This is going to produce the marginal for all hyperpriors
    
    hyperMarge = np.zeros(hyperPriors.shape)
    for donor in donorLikely:
        with open(f'{save_folder}/logDonorHyper_{donor}.pkl', 'rb') as handle:
            donorLogHyperLikely = pickle.load(handle)
        hyperMarge += donorLogHyperLikely
    
    hyperMarge = np.exp(hyperMarge)*hyperPriors

    return hyperMarge


@njit(parallel=True)
def makethetaDiffs(times, thetas):

    diffs = np.diff(thetas)
    diffs = np.append(diffs,diffs[-1])
    bigDiffs = np.zeros((len(times), len(thetas)))
    for time_index, time in enumerate(times):
        bigDiffs[time_index, :] = diffs
    return bigDiffs

def take_closest(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
        return after
    else:
        return before
    
if __name__ == '__main__':
    
    import argparse
    from time import time
    import warnings

    warnings.simplefilter("ignore")
    parser = argparse.ArgumentParser(
        description="Compute likelihood for given theta")
    parser.add_argument('-thresh', default=None, type=float,
        help="Threshold of heteroplasmy under consideration")
    parser.add_argument('-n', default=None, type=int,
        help="Number of lineages")
    parser.add_argument('-crypticDF', type=str,
        help="Input cryptic dataframe")
    parser.add_argument('-cellDict', type=str,
        help="Input cell dictionary")
    parser.add_argument('-baseDict', type=str,
        help="Input base dictionary")
    parser.add_argument('-ageDict', type=str,
        help="Input age dictionary")
    parser.add_argument('-logFold',  default='/rds/general/user/ag5818/home/monteCarloInfer/', type=str,
        help="folder with the dictionaries for log likelihoods")
    parser.add_argument('-saveFolder', default='/rds/general/user/ag5818/ephemeral/', type=str,
        help="savePath")
    
    args = parser.parse_args()
    thresh = args.thresh
    n = args.n
    crypticFile = args.crypticDF
    cellFile = args.cellDict
    baseFile = args.baseDict
    ageFile = args.ageDict
    logLikeFolder = args.logFold
    save_folder = args.saveFolder

    set_num_threads(8)
    exp = np.vectorize(mpmath.exp)
    
    with open(f'{logLikeFolder}times.pkl', 'rb') as handle:
        times = pickle.load(handle)
    with open(f'{logLikeFolder}thetas.pkl', 'rb') as handle:
        thetas = pickle.load(handle)
    with open(f'{logLikeFolder}disc_muts_full{n}.pkl', 'rb') as handle:
        disc_muts_full = pickle.load(handle)
    with open(f'{logLikeFolder}logLikelySingleCell{n}.pkl', 'rb') as handle:
        logLikely = pickle.load(handle)
    
    theta_diffs = makethetaDiffs(times, thetas)

    lowest = take_closest(disc_muts_full, thresh)
    disc_muts = disc_muts_full[np.where(disc_muts_full == lowest)[0][0]:]
    
    
    cryps = pd.read_pickle(crypticFile)
    cryps['discHF'] = cryps.apply(lambda x: take_closest(disc_muts_full, x['HF']), axis=1)
    cryps = cryps[(cryps.discHF>=lowest)]
    
    with open(baseFile, 'rb') as handle:
        av_bases = pickle.load(handle)
    with open(cellFile, 'rb') as handle:
        cells = pickle.load(handle)
    with open(ageFile, 'rb') as handle:
        ageDict = pickle.load(handle)
    
    donors = cells.keys()

    alphas = np.linspace(0,0.1,101)[1:]
    nus = np.logspace(-6,-3, 151)
    tau_alphas = np.linspace(0,1.5,151)[1:]
    tau_nus = np.linspace(0,2.0,201)[1:]

    with open(f'{save_folder}/alphas_range.pkl', 'wb') as handle:
        pickle.dump(alphas, handle, protocol=pickle.HIGHEST_PROTOCOL)
    with open(f'{save_folder}/nus_range.pkl', 'wb') as handle:
        pickle.dump(nus, handle, protocol=pickle.HIGHEST_PROTOCOL)
    with open(f'{save_folder}/tau_alphas_range.pkl', 'wb') as handle:
        pickle.dump(tau_alphas, handle, protocol=pickle.HIGHEST_PROTOCOL)
    with open(f'{save_folder}/tau_nus_range.pkl', 'wb') as handle:
        pickle.dump(tau_nus, handle, protocol=pickle.HIGHEST_PROTOCOL)

    logHypes = logHyperPrior(alphas, nus, tau_alphas, tau_nus)
    hyperPriors =  np.exp(logHypes)

    max_count = cryps.groupby(['sample_id','discHF']).count().HF.max()

    short_monte = logLikely[:max_count+1, np.where(disc_muts_full == lowest)[0][0]:,:,:]

    logAlphaDonorPriors={}
    logThetaDonorPriors={}
    alphaDonorPriors={}
    thetaDonorPriors={}

    for donor_id in cells.keys():
        logAlphaDonorPriors[donor_id] = logPopPriorAlphas(alphas, tau_alphas, times, ageDict[donor_id])
        logThetaDonorPriors[donor_id] = logPopPriorThetas(nus, tau_nus, thetas, av_bases[donor_id])
        alphaDonorPriors[donor_id] = np.exp(logAlphaDonorPriors[donor_id])
        thetaDonorPriors[donor_id] = np.exp(logThetaDonorPriors[donor_id])

    # This is the log of the population level likelihood which we will use as input into the hierarchical model.
    full_log = {}
    for donor_id in cells.keys():
        full_log[donor_id] = np.zeros((times.shape[0], thetas.shape[0]))
        total_cells = cells[donor_id]
        tuples = []
        cell_num = []
        for index, hf in enumerate(disc_muts):
            het_counter = cryps[(cryps.donor_id == donor_id)&(cryps.discHF == hf)].groupby(['sample_id']).count().groupby('POS').count()
            tuples.extend([(0, index)])
            cell_num.extend([int(total_cells-het_counter.sum().discHF)])
            for row in het_counter.iterrows():
                tuples.extend([(row[0],index)])
                cell_num.extend([row[1].discHF])        
        for index, count in enumerate(tuples):
            full_log[donor_id] = full_log[donor_id] + cell_num[index] * (short_monte[count[0], count[1], :, :])
    
    full = {}
    for donor_id in cells.keys():
        full[donor_id] = exp(full_log[donor_id])
        full[donor_id] = (full[donor_id]/full[donor_id].sum()).astype(float)

    createLogDonorLikely_newest(alphaDonorPriors, thetaDonorPriors, full, thetas)

    print('Hypermarginal Construction')
    hyperMarginal = createHyperMarginal_newest(hyperPriors, full)
    
    print('Saving...')
    with open(f'{save_folder}/posteriorHyperMarg.pkl', 'wb') as handle:
        pickle.dump(hyperMarginal, handle, protocol=pickle.HIGHEST_PROTOCOL)

    print('Donor Marginal Construction')
    for donor_id in cells.keys():
        print(donor_id)
        donorMarginal = createDonorMarginal_newest(hyperPriors, alphaDonorPriors, thetaDonorPriors, full, donor_id, nus)
        print('Saving...')
        with open(f'{save_folder}/posteriorMarg_{donor_id}.pkl', 'wb') as handle:
            pickle.dump(donorMarginal, handle, protocol=pickle.HIGHEST_PROTOCOL)

    print('Success!')
