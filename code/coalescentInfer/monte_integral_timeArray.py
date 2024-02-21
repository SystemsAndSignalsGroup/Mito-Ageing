import numpy as np # type: ignore
from numba import njit, set_num_threads, get_num_threads, config, threading_layer, prange # type: ignore
import pickle
import os # type: ignore

@njit(parallel=True, fastmath=True)
def length_level(level, lengths, total_l, tot_levels, W, samples):
    
    if level <(tot_levels-1):
        current_l = np.random.exponential(2/((tot_levels - level)*((tot_levels - level) - 1)), samples)
        continue_index = np.where(total_l+current_l<W)
        stop_index = np.where(total_l+current_l>=W)
        lengths[level][continue_index[0]] = current_l[continue_index[0]]
        lengths[level][stop_index[0]] = W - total_l[stop_index[0]]
        total_l = total_l+ lengths[level]
    elif level == (tot_levels-1):
        lengths[level] = W-total_l
    return total_l, lengths

@njit(fastmath=True)
def coalescent_fast(tot_levels, W, samples):
    ## This produces an array of size (levels-1, samples) where the level counts back from tot levels in the first axis and each column in an independent sample
    ## ie: lengths[0,:] are the samples from T_{tot_levels}, lengths[-2,:] are samples from T_2, lengths[-1,:] are the root lengths
    lengths = np.zeros((tot_levels,samples))
    total_l = np.zeros(samples)

    for level, length in enumerate(lengths):
         total_l, lengths = length_level(level, lengths, total_l, tot_levels, W, samples)
    return lengths

@njit()
def children_struct_single(tot_levels):
    
    children = np.zeros((tot_levels, tot_levels+1))
    simer = np.zeros(tot_levels, dtype = np.int64)
    simer = simer+1
    children[0,1] = tot_levels
    for level in range(1, tot_levels):
        indexes =  np.array(list(range(tot_levels-level+1)), dtype = np.int64)
        random = np.random.choice(indexes, 2, replace=False)
        mask = np.zeros(simer.shape[0], dtype=np.int64) == 0
        mask[random] = False
        simer = np.append(simer[mask], simer[random].sum())
        countsim = np.sort(simer)
        counter_childs = countsim[np.append(0, np.where(np.diff(countsim)!= 0)[0]+1)]
        counter_count = np.zeros(counter_childs.shape[0], dtype=np.int64)
        for index, i in enumerate(counter_childs):
            counter_count[index] = np.where(countsim == i)[0].shape[0]
        children[level][counter_childs] = counter_count

    children = children[:,1:].T
    return children

@njit(parallel=True, fastmath=True)
def length_samples(tot_levels, W, coal_lengths, samples):

    lengths = np.zeros((tot_levels, samples))
    for samp in prange(samples):
        lengths[:,samp] = (children_struct_single(tot_levels) * coal_lengths[:, samp]).sum(axis=1)
    return lengths

@njit(parallel = True, fastmath=True)
def likely_theta_speed(factorial, mutant_number, total_length_samples, thetas):
    ########
    # Returns the unnormalised log likelihood for a set of theta and time values and number of mutations.
    # The returned logLikely is accessed so that logLikely[i,j] is the (unnormalised) log likelihood 
    # for the probability a cell has mutant_number mutatations on i leaves for parameter values
    # times[j], thetas[k]

    samples = total_length_samples.shape[-1]
    logLikely = np.zeros((total_length_samples.shape[0], thetas.shape[0]))
    for theta_index, theta in enumerate(thetas):
        logLikely[:, theta_index] += np.log(((theta*total_length_samples)**mutant_number*np.exp(-theta*total_length_samples)/(factorial*samples)).sum(axis=1))

    return logLikely

if __name__ == '__main__':

    import argparse
    import warnings
    
    warnings.simplefilter("ignore")

    parser = argparse.ArgumentParser(
        description="Compute probability of having m mutations for sample size n for range of theta and coalescent time")
    parser.add_argument('-threads', default=8, type=int,
        help="Number of threads you want to parallelise across (Min 8 recommended)")
    parser.add_argument('-tot_levels', default=200, type=int,
        help="Total levels interested in")
    parser.add_argument('-time_index', default=None, type=int,
        help="Time index you want")
    parser.add_argument('-save_path', default=None, type=str,
        help="The file path for both the times and thetas you want to consider, as well as where you will save it")
    
    
    
    args = parser.parse_args()
    threads = args.threads
    tot_levels = args.tot_levels
    time_index = args.time_index
    save_path = args.save_path


    config.THREADING_LAYER = 'threadsafe'
    set_num_threads(threads)
    print("Threading layer chosen: %s" % threading_layer(), flush = True)
    
    with open(f'./times.pkl', 'rb') as handle:
        times = pickle.load(handle)
    with open(f'./thetas.pkl', 'rb') as handle:
        thetas = pickle.load(handle)
               
    time = times[time_index]
    samples = 1000000
    
    print(f'Coalescent times starting', flush=True)
    coal_lengths = coalescent_fast(tot_levels, time, samples)
    print(f'Coalescent times done!', flush=True)

    print(f'length_samples starting, threads = {get_num_threads()}', flush=True)
    
    all_lengths = length_samples(tot_levels, time, coal_lengths, samples = samples)
        
    print(all_lengths.shape)
    print('length_samples done', flush=True)
        
    for mutant_num in range(9):
        
        try: 
            os.mkdir(f'{save_path}/{tot_levels}M={mutant_num}')
        except OSError as error: 
            print(error) 

        fact = np.math.factorial(mutant_num)
            
        print(f'likely_theta_speed starting {mutant_num}', flush=True)
        loger = likely_theta_speed(fact, mutant_num, all_lengths, thetas)
        print('likely_theta_speed done', flush=True)

        with open(f'{save_path}/{tot_levels}M={mutant_num}/{time_index}.pkl', 'wb') as handle:
            pickle.dump(loger, handle)