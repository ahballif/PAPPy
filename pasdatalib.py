
# Current as of July 5, 2021 at 9:30 PM - Addison Ballif

# I library of functions and tools that I (Addison Ballif) am 
# creating to do simple data manipulation for Positron Annihilation
# Spectroscopy (PAS). 
# 
# The file here is slightly different that the one used in my jupyter files. This 
# version is specific to the PAPPy project.  

def importFile(filepath):
    start_reading = False

    bins = []
    energy = []
    counts = []

    with open(filepath) as csv:
        for line in csv:
            data = line.strip().split(',')
            if start_reading: #actually read the data
                bins.append(int(data[0]))
                energy.append(float(data[1]))
                counts.append(int(data[2]))

            if data[0] == 'Channel':
                #wait until the misc data has been passed and the real information has started
                start_reading = True
            elif not start_reading:
                #here we can get the misc data from the file
                if 'Elapsed Computational' in data[0]:
                    number_of_samples = int(data[1])

    return bins, energy, counts


# --------- CROPPING FUNCTIONS --------------

# This version is different because it includes the option to use an upper and lower bound rather than center and pm. 
def cropWindow(bins, energy, counts, lower, upper):
    #lower and upper bounds are in channels
    lower_idx = int(lower - min(bins)) #convert bins to indices
    upper_idx = int(upper - min(bins))

    energy_t = energy[lower_idx:upper_idx]
    bins_t = bins[lower_idx:upper_idx]
    counts_t = counts[lower_idx:upper_idx]

    return bins_t, energy_t, counts_t 

# the cropping function that doesn't take in energy. This one is used for parameter calculations. 
def cropWindow2(bins, counts, lower, upper):
    #same as cropWindow, but it doesn't include energy. 
    #lower and upper bounds are in channels
    lower_idx = int(lower - min(bins)) #convert bins to indices
    upper_idx = int(upper - min(bins))

    bins_t = bins[lower_idx:upper_idx]
    counts_t = counts[lower_idx:upper_idx]

    return bins_t, counts_t 


###### Conversion functions #########
def convertTokeV(value, bins, energy):
    #converts channels to keV. 
    try:
        return float(value) *(max(energy)-min(energy))/(max(bins)-min(bins))
    except ValueError:
        return value #failed because bins and energy lists are empty

def convertToChannels(value, bins, energy):
    #converts keV to channels
    try:
        return int(float(value) *(max(bins)-min(bins))/(max(energy)-min(energy)))
    except ValueError:
        return value #failed because bins and energy lists are empty



####### EXPORT PARAMETERS ###########
# a function that exports the parameters in csv file format
def exportParameters(filepath, headers, *args):
    try:
        with open(filepath, 'w') as file:
            lines = []

            #make the header line
            headerline = ''
            for header in headers:
                headerline += header+','
            lines.append(headerline+'\n')

            #add the data rows
            for i in range(len(args[0])):
                line = ''
                for arg in args:
                    line += str(arg[i]) +','
                lines.append(line+'\n')

            file.writelines(lines)
        print(f'Save successfully to "{filepath}"')
    except FileNotFoundError:
        print(f'Invalid Filepath for CSV Export. ')        





# a little something to help people know what to do. 
if __name__ == '__main__':
    print('\nThis is a help file to the PAPPy project. To run the program, run the PAPPy.py file. ')
    