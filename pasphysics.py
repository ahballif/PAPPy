from scipy.optimize import curve_fit
from scipy.special import erf
from numpy import float_power as power
from numpy import exp, pi, sqrt, trapz, sum, array, mean, std, linspace, transpose, asarray
import pasdatalib as pdl
from numpy.random import poisson

# Current as of July 5, 2021 at 9:30 PM - Addison Ballif

# This file contains most of the physics calculations for the PAPPy program. 
# The variables where the data stored are stored in the mainApplication class in the
# main file, but the actual calculation code is here.  

# -------- BACKGROUND MODELS -----------
# To add a background model, add an option as a string in the background_models list,
# and then add a case to the if statement tree inside the find_background_shape function. You can then 
# return the background function. Data available to you for background calculation is the energy (x-data)
# and the counts (y-data), as well as the lower and upper bounds of the curve part of the data. Data inside these bounds are the curve, and 
# data outside these bounds is considered background that we can use to find the height of the background. 
# 
# Later we will add it so that we return the function and the uncertainty list.  

background_models = ['linear' ,'erf', 'linear monte-carlo', 'none']

#functions used for different models
def gaussian(x, mu, sigma, k, C):
    # simple gaussian
    return k*exp(-(x-mu)**2 / (2 * sigma**2)) + C


def linear(x_list, left_height, right_height, lower_bound, upper_bound):
    # a line. 
    output = []
    for x in x_list:
        if x < lower_bound:
            output.append(left_height)
        elif x > upper_bound:
            output.append(right_height)
        else:
            m = (right_height-left_height)/(upper_bound-lower_bound)
            b = left_height - m*lower_bound
            output.append(m*x + b)
    return output

def comp_erf(x, mu, sigma, k, C):
    # complimentary error function, but using parameters
    return k*erf((mu-x)/(sigma*sqrt(2))) + C


# ------------ THE MAIN BACKGROUND FUNCTION --------------
def find_background_shape(bins, energy, counts, lower, upper, type):
    #finds the background function
    if type == 'linear' or type == 'linear monte-carlo':
        #use a linear model
        lower_idx = lower - min(bins)
        upper_idx = upper - min(bins)
        
        bl = mean(counts[:lower_idx-1])
        bu = mean(counts[upper_idx+1:])

        background_line = linear(bins, bl, bu, lower, upper)

        #uncertainties
        bl_un = sqrt(sum(counts[:lower_idx]))/(lower_idx + 1) #pay attention to the difference between bins and indices
        bu_un = sqrt(sum(counts[upper_idx:]))/(len(bins) - upper_idx) #channels/bins start at 1, indices start at 0

        background_un = sqrt((1 + (lower-array(bins))/(upper-lower))**2 *bl_un**2 + ((array(bins)-lower)/(upper-lower))**2 *bu_un**2)

        return background_line, background_un

    elif type == 'erf':
        # a version where we fit a gaussian to find a standard deviation and mean, and then we fit the s curve to find the height and scale. 
        peakCenterGuess = bins[find_nearest_idx(energy, 511)]
        p0 = [peakCenterGuess, 3, max(counts), 0.1*max(counts)]  
        


        #first fit gaussian to find std and mu. 
        params, cov = curve_fit(gaussian, bins, counts, p0=p0)

        fit_sigma = params[1]
        fit_mu = params[0]


        #cut out the peak for the next curve fit
        lower_idx = lower - min(bins)
        upper_idx = upper - min(bins)

        energy_t = bins[:lower_idx-1] + bins[upper_idx+1:]
        counts_t = counts[:lower_idx-1] + counts[upper_idx+1:]


        #fitting to find height and scale. lambda is used to lock the mu and sigma parameters. 
        params2, cov2 = curve_fit(lambda x, k, C: comp_erf(x, fit_mu, fit_sigma, k, C), energy_t, counts_t )

        #returning 0 for uncertainty because we are going to use monte-carlo method for the erf function
        return comp_erf(bins, fit_mu, fit_sigma, params2[0], params2[1]), array(bins)*0

    elif type == 'none':
        return array(bins)*0, array(bins)*0
    else:
        print('Error - invalid background model')


# --------- OTHER CALCULATION FUNCTIONS ----------
# Someday I want to add deconvolution to the preprocessing. I spent alot of time trying to deconvolve the data. I was not
# able to get something meaningful. The idea was that we could find the point spread function by sticking a known frequency spike in the data. 
# Because of blur, this spike would be measured as a gaussian. We could then use the width of this gaussian to find the width of the psf, which
# would be used for deconvolution.  

def deconvolve_data(bins, counts, psf):
    #this is an idea for the future
    pass


# ------------- BOUND GUESS FUNCTIONS --------------------
# I thing I wanted to add was functions that make guesses for different bounds. Because the parameter values are arbitrary, it is okay
# to use these. I also thought it was useful for the main curve bound guess. The user can still manually adjust the bounds.  

def make_s_parameter_bound_guess(bins, counts):
    # uses a gaussian fit to find the z-value and the bounds for a 0.5 s-parameter. 
    p0 = [5110, 3, max(counts), 0.1*max(counts)]

    #first fit gaussian to find std and mu. 
    params, cov = curve_fit(gaussian, bins, counts, p0=p0)
    
    #standard deviations to get a s-parameter of 0.5
    z_score = 0.725 #should be 0.6745, altered to adjust for the fact that it's not actually gaussian. 

    mu = params[0]
    sigma = params[1]

    lower = mu-sigma*z_score
    upper = mu+sigma*z_score

    return lower, upper

#the w parameter guess functions are a new iteration. To see the old version, go to the end of the file. 

def make_left_w_parameter_bound_guess(bins, counts, curveBounds):
    # uses a gaussian fit to find the z-value and the bounds for a 0.5 s-parameter. 
    p0 = [5110, 3, max(counts), 0.1*max(counts)]

    #first fit gaussian to find std and mu. 
    params, cov = curve_fit(gaussian, bins, counts, p0=p0)
    
    #standard deviations to get a s-parameter of 0.5
    lower_z_score = 5 #should be 0.6745, altered to adjust for the fact that it's not actually gaussian. 
    upper_z_score = 2

    mu = params[0]
    sigma = params[1]

    lower = mu-lower_z_score*sigma
    upper = mu-upper_z_score*sigma

    return lower, upper

def make_right_w_parameter_bound_guess(bins, counts, curveBounds):
    # uses a gaussian fit to find the z-value and the bounds for a 0.5 s-parameter. 
    p0 = [5110, 3, max(counts), 0.1*max(counts)]

    #first fit gaussian to find std and mu. 
    params, cov = curve_fit(gaussian, bins, counts, p0=p0)
    
    #standard deviations to get a s-parameter of 0.5
    lower_z_score = 2 #should be 0.6745, altered to adjust for the fact that it's not actually gaussian. 
    upper_z_score = 5

    mu = params[0]
    sigma = params[1]

    lower = mu+lower_z_score*sigma
    upper = mu+upper_z_score*sigma

    return lower, upper

def make_curve_bound_guess(bins, energy, counts):
    # uses a gaussian fit to find the bounds for z-score of 5
    p0 = [511, 3, max(counts), 0.1*max(counts)]

    #first fit gaussian to find std and mu. 
    params, cov = curve_fit(gaussian, energy, counts, p0=p0)

    #number of standard deviations
    z_score = 6 

    mu = params[0]
    sigma = params[1]

    lower = mu-sigma*z_score
    upper = mu+sigma*z_score

    lowerBin = bins[find_nearest_idx(energy, lower)]
    upperBin = bins[find_nearest_idx(energy, upper)]

    return lowerBin, upperBin

def find_nearest_idx(array, value):
    array = asarray(array)
    idx = (abs(array - value)).argmin()
    return idx


# ------------------ PARAMETER CALCULATION ------------------

def calculate_parameter(bins, counts, background, background_un, curveBounds, paramBounds, background_type, **kwargs):
    #this function will return a parameter given an upper and lower bound in the energy spectrum  
    #the uncertainty for the parameter depends on the type of background. These will be the same times included in the
    # background list above. 

    #curveBounds and paramBounds are each a tuple with two arguments (lower, upper)

    #through kwargs there is an option to add a progress bar for monte-carlo calculations. Because this is used through
    # kwargs, it is optional. This function should be able to used without a GUI, if needed. 
    # There is also a progress_value option to monitor the actual progress as a fraction of 1 

    if background_type == 'erf':
        #using the gaussian error function model 
        #using monte-carlo method and recalculation of the background for each iteration. 
        n_passes = kwargs.get('number_of_iterations', 500) #number of monte-carlo iterations
        progress_bar_update = 20 #the higher it is, the less smooth the progress bar will be, but the faster the calculation speed. 

        data = []
        parameters = []

        #get wobbled data
        for channel_lambda in counts:
            #wobble the original data according to a poisson uncertainty distribution
            data.append(poisson(channel_lambda, n_passes))
        data = transpose(data).tolist() #flip rows and columns to get an array of wobbled data sets
        
        #reset progress bar
        if not kwargs.get('progress_bar', None) == None:
            #update the progress_bar
            progress = kwargs.get('progress_bar', None)
            progress.updateProgress()

        # ----- MONTE CARLO ------
        for each_counts in data:

            #update the progress bar
            if not kwargs.get('progress_bar', None) == None and data.index(each_counts)%progress_bar_update == 0:
                #update the progress_bar
                progress = kwargs.get('progress_bar', None)
                progress.updateProgress(progress_value = float(data.index(each_counts)+1)/float(len(data)))

            #calculate background
            background, _ = find_background_shape(bins, each_counts, curveBounds[0], curveBounds[1], 'erf')

            subtracted = array(each_counts) - array(background)

            #crop to get bounds
            energy_curve, counts_curve = pdl.cropWindow2(bins, subtracted, curveBounds[0], curveBounds[1])
            energy_param, counts_param = pdl.cropWindow2(bins, subtracted, paramBounds[0], paramBounds[1])

            #next integrate to find total and area under the s
            total_area = sum(counts_curve)
            param_area = sum(counts_param)

            if kwargs.get('trapz', False):
                #optional trapezoidal integration 
                    #x scale doesn't matter because our channels have uniform spacing, 
                    # and we are only looking at a ratio of areas. 
                total_area = trapz(counts_curve)
                param_area = trapz(counts_param)
                
            #divide to get s-parameter
            parameter = param_area/total_area
            parameters.append(parameter)

        #now report mean and std of the parameters distribution

        return mean(parameters), std(parameters)
        

    elif background_type == 'linear monte-carlo':
        #using monte-carlo for linear uncertainty 
        n_passes = kwargs.get('number_of_iterations', 500) #number of monte-carlo iterations
        progress_bar_update = 20 #the higher it is, the less smooth the progress bar will be, but the faster the calculation speed. 

        data = []
        parameters = []

        #get wobbled data
        for channel_lambda in counts:
            #wobble the original data according to a poisson uncertainty distribution
            data.append(poisson(channel_lambda, n_passes))
        data = transpose(data).tolist() #flip rows and columns to get an array of wobbled data sets
        
        #reset progress bar
        if not kwargs.get('progress_bar', None) == None:
            #update the progress_bar
            progress = kwargs.get('progress_bar', None)
            progress.updateProgress()

        # ----- MONTE CARLO ------
        for each_counts in data:

            #update the progress bar
            if not kwargs.get('progress_bar', None) == None and data.index(each_counts)%progress_bar_update == 0:
                #update the progress_bar
                progress = kwargs.get('progress_bar', None)
                progress.updateProgress(progress_value = float(data.index(each_counts)+1)/float(len(data)))

            #calculate background
            background, _ = find_background_shape(bins, each_counts, curveBounds[0], curveBounds[1], 'linear')

            subtracted = array(each_counts) - array(background)

            #crop to get bounds
            energy_curve, counts_curve = pdl.cropWindow2(bins, subtracted, curveBounds[0], curveBounds[1])
            energy_param, counts_param = pdl.cropWindow2(bins, subtracted, paramBounds[0], paramBounds[1])

            #next integrate to find total and area in the parameter bounds
            total_area = sum(counts_curve)
            param_area = sum(counts_param)

            if kwargs.get('trapz', False):
                #optional trapezoidal integration 
                    #x scale doesn't matter because our channels have uniform spacing, 
                    # and we are only looking at a ratio of areas. 
                total_area = trapz(counts_curve)
                param_area = trapz(counts_param)
                
            #divide to get s-parameter
            parameter = param_area/total_area
            parameters.append(parameter)

        #now report mean and std of the parameters distribution

        return mean(parameters), std(parameters)

    else:
        #assuming that background is independant and that each channel is independant. 

        #first subtract background
        subtracted = array(counts) - array(background)
        subtracted_un = sqrt(array(counts) + array(background_un)**2)

        #crop to get bounds
        energy_curve, counts_curve = pdl.cropWindow2(bins, subtracted, curveBounds[0], curveBounds[1])
        _,                un_curve = pdl.cropWindow2(bins, subtracted_un, curveBounds[0], curveBounds[1])
        energy_param, counts_param = pdl.cropWindow2(bins, subtracted, paramBounds[0], paramBounds[1])
        _,                un_param = pdl.cropWindow2(bins, subtracted_un, paramBounds[0], paramBounds[1])

        #next integrate to find total and area under the s
        total_area = sum(counts_curve)
        param_area = sum(counts_param)

        total_un = sqrt(sum(array(un_curve)**2))
        param_un = sqrt(sum(array(un_param)**2))

        if kwargs.get('trapz', False):
            #optional trapezoidal integration
            total_area = trapz(counts_curve, x=energy_curve)
            param_area = trapz(counts_param, x=energy_param)

            total_un /= sqrt(2)
            param_un /= sqrt(2)

        

        #divide to get s-parameter
        parameter = param_area/total_area
        parameter_un = sqrt((param_un/total_area)**2 + (param_area*total_un/total_area**2)**2)

    return parameter, parameter_un




# a little something to help people know what to do. 
if __name__ == '__main__':
    print('\nThis is a help file to the PAPPy project. To run the program, run the PAPPy.py file. ')
    


## Below are the old w parameter guess functions. I changed the guess calculations because I think 
# the new technique is a better approach. The old technique picks a random region in between the s and lower or upper curve bounds. 
# The new versions use the standard deviation of the gaussian fit. 
# 
# I left the functions commented incase they are prefered by someone doing alot of work with W parameters. It shouldn't 
# matter too much anyway because the user can still set the bounds.  

# def make_left_w_parameter_bound_guess(bins, counts, curveBounds):
#     # uses a gaussian fit to find the z-value and the bounds for a 0.5 s-parameter. 
#     p0 = [5110, 3, max(counts), 0.1*max(counts)]

#     #first fit gaussian to find std and mu. 
#     params, cov = curve_fit(gaussian, bins, counts, p0=p0)
    
#     #standard deviations to get a s-parameter of 0.5
#     z_score = 0.725 #should be 0.6745, altered to adjust for the fact that it's not actually gaussian. 

#     mu = params[0]
#     sigma = params[1]

#     #pick 3/12 and 7/12 from lower bound and lower_s_bound

#     lower_s_bound = mu-sigma*z_score

#     lower = (lower_s_bound-curveBounds[0])*3/12 + curveBounds[0]
#     upper = (lower_s_bound-curveBounds[0])*7/12 + curveBounds[0]

#     return lower, upper

# def make_right_w_parameter_bound_guess(bins, counts, curveBounds):
#     # uses a gaussian fit to find the z-value and the bounds for a 0.5 s-parameter. 
#     p0 = [5110, 3, max(counts), 0.1*max(counts)]

#     #first fit gaussian to find std and mu. 
#     params, cov = curve_fit(gaussian, bins, counts, p0=p0)
    
#     #standard deviations to get a s-parameter of 0.5
#     z_score = 0.725 #should be 0.6745, altered to adjust for the fact that it's not actually gaussian. 

#     mu = params[0]
#     sigma = params[1]

#     #pick 5/12 and 9/12 from lower bound and lower_s_bound

#     upper_s_bound = mu+sigma*z_score

#     lower = (curveBounds[1] - upper_s_bound)*5/12 + upper_s_bound
#     upper = (curveBounds[1] - upper_s_bound)*9/12 + upper_s_bound

#     return lower, upper
