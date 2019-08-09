#!/usr/bin/python

'''
Main handler for pySS
'''

# imports
import os, sys, re, csv
import pytraj as pt
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
# from vmd import *


def pow_curve(r):
    '''
    Rough model for H bond distance-energy relationship. Takes r in
    Angstroms and returns energy in kcal/mol.
    '''
    # energy in kJ/mol
    kJ = (-1.2*(10**10))/((r*100)**3.78)
    return kJ/4.184
# end function pow_curve

def lint(r,acc_elem,don_elem):
    '''
    Linear model with intercept from Wendler 2010.
    (https://pubs.acs.org/doi/abs/10.1021/jp103470e). Parameterized based
    on donor and acceptor heteroatom elements. Will implement later
    if pow_curve isn't good enough.
    '''

    # parameterizations for different heteroatoms as acceptors
    pacc = {
        "O": 0
    }

    # parameterizations for different heteroatoms as acceptors
    pdon = {
        "O": 0
    }

    return x

# end function lint

def rmse(predictions, targets=None):
    '''
    Root mean square error. If optional arg targets is not provided,
    defaults to arithmetic mean of predictions.
    '''
    if targets == None:
        targets = np.ones_like(predictions) * np.mean(predictions)

    return np.sqrt(((predictions - targets) ** 2).mean())

# end function rmse

def plot_correl(x, y, func="pearson"):
    '''
    Calculates Pearson coefficient and p value between normalized ndarrays
    x and y. Plots using matplotlib and returns tuple like scipy.stats.pearsonr.
    '''

    if func == "pearson":
        x = norm_dev(x)
        y = norm_dev(y)
        pf = stats.pearsonr(x, y)
    elif func == "spearman":
        x = x - np.mean(x)
        y = y - np.mean(y)
        pf = stats.spearmanr(x, y)
    else:
        raise Exception

    # plot it!

    # Pearson coefficient plot
    plt.plot(x, y, "r*")
    plt.plot(x, pf[0]*x, "k-")
    plt.show()
    return pf

# end function plot_correl

def norm_dev(array):
    '''
    Given numpy array, calculates normalized deviation from the mean. Deviation
    is normalized to 1 standard deviation.
    '''

    return (array - np.mean(array)) / np.std(array)

# end function norm_dev

def reject_outliers(array, m=2):

    outliers = np.where(abs(array - np.mean(array)) >= m * np.std(array))[0]
    array[outliers] = np.mean(array)
    return array

# end function reject_outliers

def is_calc_root(path):
    '''
    Determines whether the directory at path is a root directory, e.g. whether
    it has *_wat.prmtop and *_wat.inpcrd files and a in/ directory.
    '''

    # check if it is a directory
    if not os.path.isdir(path):
        return False
    # has paths to *_wat.prmtop, *_wat.inpcrd, and ./in/
    req_ct = 0
    for f in os.listdir(path):
        if "_wat.prmtop" in f:
            req_ct += 1
        elif "_wat.inpcrd" in f:
            req_ct += 1
        elif f == "in":
            req_ct += 1
    if req_ct >= 3:
        return True
    else:
        return False

#end method is_calc_root

class OutReader(object):

    def __init__(self):
        '''
        Class responsible for reading, storing, and re-writing an AMBER .out
        file.
        '''

        self.cycles = {}

    #end method __init__

    def load_out(self,read_fp):
        '''
        Reads AMBER.sander .out file at file path self.read_fp and loads
        NSTEP information into the self.cycles dict of dicts.

        Arguments:
            str read_fp - relative file path of .out file to read
        '''

        def parse_raw(raw):
            '''
            Parses list of words in NSTEP entry into a dictionary. Tries float
            conversion of each element after '=', and loads as dict value.

            Arguments:
                list raw - input from below. see below
            Returns:
                dict
            '''

            # empty dict
            d = {}
            # placeholder
            key = ""
            # 1 means key is being read, 2 is '=', 3 is reading value
            i = 1
            for word in raw:
                if word == "=": #equals sign. get ready for value
                    i = 3
                    continue
                elif i == 1: # key
                    key += ' ' + word
                elif i == 3: # value is being read
                    d[key.strip()] = float(word) # try float conversion
                    key = "" # reset key
                    i = 1 # get ready to read next key
                else: # what happened here
                    raise Exception()
            # NSTEP value should be an integer
            d['NSTEP'] = int(d['NSTEP'])
            return d

        #end function parse_raw

        # clear out the cycles dict
        self.cycles = {}
        # list of all words in NSTEP entry. includes multiple lines
        raw_entry = []
        # switch for loading into cycles dict or not
        loading = False

        with open(read_fp,'r') as f:
            for full_line in f:
                # strip whitespace
                full_line = full_line.strip()
                # found start of NSTEP info
                if full_line[:5] == "NSTEP":
                    raw_entry = [] # start new entry
                    loading = True # switch goes on
                elif loading and (full_line[:5] == "-----"):
                    # line at end of entry. wrap it up
                    entry = parse_raw(raw_entry)
                    # key is int type value of NSTEP
                    self.cycles[entry['NSTEP']] = entry
                    loading = False
                if loading: # still reading the same NSTEP entry
                    raw_entry.extend(full_line.split())
        # okay, it has some error handling
        if not raw_entry:
            raise Exception(
                "No cycle entries were successfully read in file '{}'".format(read_fp))
        return self.cycles

    #end method load_out

    def write(self,write_fp,append=False):
        '''
        Writes information stored in self.cycles dict of dicts to a tab
        delimited .txt file at write_fp.

        Arguments:
            str write_fp - relative file path of tab-delimited .txt to write
            bool append[=False] - whether to append or overwrite this file
        '''

        # make base directory if it does not exist
        basedir = os.path.dirname(write_fp)
        if not os.path.exists(basedir):
            print("Write directory at path '{}' did not ".format(basedir) +
                  "exist at runtime. Making directory.")
            os.makedirs(basedir)

        # mode passed to open() function
        if append:
            mode = 'a'
        else:
            mode = 'w'

        # list of lists of cells separated horizontally by \t and vertically
        # by \n
        text_lst = []
        # same thing but str type
        text = ""

        # get header line and put in str type text var
        header_lst = self.cycles.values()[0].keys()
        text_lst.append(header_lst)

        for nstep, entry in self.cycles.items():
            row = []
            for name in header_lst:
                try:
                    row.append(str(entry[name]))
                except KeyError:
                    row.append("")
            text_lst.append(row)
        # add to the str type var
        for row in text_lst:
            text += "\t".join(row) + "\n"
        #print(text_lst)
        # write text to file at write_fp
        print("Writing to file at '{}'...".format(write_fp))
        with open(write_fp,mode) as f:
            f.write(text)
    #end method write

    def get_EPtot(self):
        '''
        Returns EPtot values for all NSTEPS as numpy array. kJ/mol
        '''

        EPtot_lst = []
        for cycle,attrs in self.cycles.items():
            EPtot_lst.append(attrs['EPtot'])
        return np.array(EPtot_lst)

    #end method get_Etot

#end class OutReader

class Struct(object):

    def __init__(self):
        '''
        Handling class for analyses of EndoG dimer bound to ssDNA. Looks at
        H-bond between Cys110 and the 5HM hydroxyl.
        '''

        pass

    #end method __init__

    def loadtraj(self,prmtop,mdcrd):
        '''
        Loads a trajectory and prmtop files into Struct using pytraj.
        '''

        print("Loading prmtop file '{}' with mdcrd file '{}'...".format(prmtop,mdcrd))
        self.traj = pt.iterload(mdcrd,prmtop)
        #print(self.traj)
        return self.traj

    #end method loadmol

    def load_EndoG_bonds(self, calc_eelec=False):
        '''
        Loads relevant EndoG H-bonds as Bond objects. These include:
            -Cys110H-5HMO (CH_O5)
            -Cys110S-5HMH (CS_HO5)
            -Cys110H-PO4 (CH_OP2)
        for each monomer = 6 total H bonds per dimer structure.

        Arguments:
            bool calc_eelec[=False] - whether to calculate electronic portion
                of LIE interaction energy
        '''

        # dict of distances in Ang considered to be maximum for H-bond
        max_d = {
            'SG_O5': 3.9,
            'SG_OP2': 3.9,
            'SG_HO5': 2.6,
            'HG_O5': 2.6,
            'HG_OP2': 2.2
        }

        # cpptraj type masks for atom selections
        masks = {
            'SG_O5_A': ":51@SG :446@O5",
            'SG_OP2_A': ":51@SG :446@OP2",
            'SG_HO5_A': ":51@SG :446@HO5",
            'HG_O5_A': ":51@HG :446@O5",
            'HG_OP2_A': ":51@HG :446@OP2",
            'SG_O5_B': ":272@SG :451@O5",
            'SG_OP2_B': ":272@SG :451@OP2",
            'SG_HO5_B': ":272@SG :451@HO5",
            'HG_O5_B': ":272@HG :451@O5",
            'HG_OP2_B': ":272@HG :451@OP2"
        }

        # energy calculation masks. Uses the masks if eelec_masks throws
        # KeyError
        eelec_masks = {
            'HG_O5_A': ":51@HG :446@O5",
            'HG_OP2_A': ":51@HG :446@OP2",
            'HG_O5_B': ":272@HG :451@O5",
            'HG_OP2_B': ":272@HG :451@OP2"
        }

        hb_dist = {}
        int_freq = {}
        eelec = {}
        for chain in ('A','B'):
            for key in ('SG_O5','SG_OP2','SG_HO5','HG_O5','HG_OP2'):
                keychain = "{}_{}".format(key,chain) # key with chain identifier
                print("Loading EndoG bond {}...".format(keychain))
                # load interaction distances over all frames in self.traj
                hb_dist[keychain] = pt.distance(self.traj, masks[keychain])
                # frequencies at which each interaction is below the max_d
                int_freq[keychain] = np.mean(hb_dist[keychain] < max_d[key])
                # calculate electronic portion of LIE
                if calc_eelec:
                    try: # in kcal/mol
                        eelec[keychain] = pt.analysis.energy_analysis.lie(
                            self.traj, mask=eelec_masks.get(keychain, masks[keychain]),
                            dtype='ndarray')[0] / 4.184
                    except KeyError: # skip if key is not in eelec_masks
                        continue

        # load into attrs
        self.int_freq = int_freq
        self.hb_dist = hb_dist
        if calc_eelec:
            eelec['HG_O5_sum'] = eelec['HG_O5_A'] + eelec['HG_O5_B']
            eelec['HG_OP2_sum'] = eelec['HG_OP2_A'] + eelec['HG_OP2_B']
            #print((eelec['HG_O5_sum'],eelec['HG_O5_A'],eelec['HG_O5_B']))
            self.eelec = eelec

    # end function load_EndoG_bonds

    def get_val(self, d_name, keys):
        '''
        Given str self.d_name dictionary, retrieves values for each str key in
        list keys. Supports * single character wildcard, which will return the
        sum of the values for all the matching keys. Appends None to return
        list if KeyError.

        Arguments:
            str d_name - name of dictionary to retrieve from self
            list keys - list of str type keys to retrieve from dict d_name
        Returns:
            list of values for each key in keys
        '''
        assert False
        assert type(d_name) == str
        assert type(keys) == list

        # safely get dictionary
        d = getattr(self, d_name, None)
        if not d or type(d) != dict:
            raise Exception("Could not find dictionary Struct.{}".format(d_name))

        # loop through keys
        return_list = []
        for k in keys:
            # has wildcard *
            if "*" in k:
                for existing_k in d.keys():
                    k_mod = k.replace("*",".")
                    match = re.search(k_mod, existing_k)
                    print(match.group())
            else: # no wildcard. safely get value
                return_list.append(d.get(k, None))

    # end method get_val

    def load_OutReader(self,fp):
        '''
        Loads OutReader object reading from .out file fp. Also runs get_EPtot
        for the OutReader.
        '''

        self.out_reader = OutReader()
        self.out_reader.load_out(fp)
        self.EPtot = self.out_reader.get_EPtot()
        return self.out_reader

    #end method load_OutReader

    def write_low_EPtot(self,fps,write_files=True):
        '''
        Gets the len(fps) frames with the lowest EPtot and writes them to files
        at str file paths fps.

        Arguments:
            list fps - list of file paths to write to
            bool write_files[=True] - whether to actually write frames to fps
                file names
        '''

        # get EPtot and the n=len(fps) frames with lowest EPtot
        try:
            EPtot = self.EPtot
        except AttributeError as e:
            print("Please run Struct.load_OutReader first")
            raise e
        self.lowest_EP_idx = EPtot.argsort()[:len(fps)]

        # check to make sure number of frames in self.traj is same as
        # the number of EPtot values over time
        if len(self.traj) != len(EPtot):
            raise Exception(
                "mdcrd file has {} frames, but .out ".format(len(self.traj)) +
                "file has {} frames.".format(len(EPtot)))

        # write lowest EPtot frames from trajectory using pytraj
        # writes to .rst extension as a NetCDF restart file
        if write_files:
            pt.write_traj('{}.rst'.format(fps[0]), self.traj, "ncrestart",
                          frame_indices=self.lowest_EP_idx, overwrite=True,
                          options="keepext")

        return self.lowest_EP_idx

    #end method write_low_EPtot

#end class Struct

def convert_lmin_frames(nframes,calc_dir_parent="",plot_EPtot=False):
    '''
    Starting at a directory calc_dir_parent, looks for Amber calculation
    directories (defined by function is_calc_root). Finds every matching pair
    of .mdcrd and .out files in each calculation directory, and finds the
    nframes frames with the lowest EPtot values in the .mdcrd file. Then writes
    each frame as a rst7 (ASCII) formatted coordinate file at
    {calc_dir_parent}/{calc_dir}/md/{trial}/*_md_f*.rst. These files are
    formatted for input into cpptraj NetCDF file conversion, functionality that
    is built into the batch_run.py script as of v0.3.8. Returns dictionary where
    keys are string file paths for found .mdcrd files and values are arrays of
    ints corresponding to the indices of the lowest EPtot frames found.

    Arguments:
        int nframes - number of frames to convert. Soft cap at 100
        str calc_dir_parent[=""] - file path pointing to parent dir of
            calculation directories. By default, set to the directory containing
            the pySS.py script.
        bool plot_EPtot[=False] - determines whether to show plot of EPtot for
            every mdcrd file found.
    '''
    # assertions
    assert type(nframes) == int
    assert 1 < nframes < 101

    # script parent dir abs file path. set default if necessary
    if not calc_dir_parent:
        script_dir_abs = os.path.dirname(os.path.realpath(__file__)) + "/"
    else:
        script_dir_abs = calc_dir_parent

    # return dict
    lmin_dict = {}

    # loop through pySS.py parent dir, looking for calc dirs containing
    # mdcrd files from md runs
    for d in os.listdir(script_dir_abs):
        prmtop = ""
        calc_dir = script_dir_abs + d + "/"
        # skip if not a calculation dir
        if not is_calc_root(calc_dir):
            continue
        # look for .prmtop
        for pt in os.listdir(calc_dir):
            if "_wat.prmtop" in pt and pt[0] != ".":
                prmtop = script_dir_abs + d + "/" + pt
        # loop through trial directories
        for td in os.listdir(calc_dir + "md/"):
            trial_dir = calc_dir + "md/" + td + "/"
            if not os.path.isdir(trial_dir):
                continue
            # loop through mdcrd files in trial directories
            for mdcrd in os.listdir(trial_dir):
                if mdcrd[-6:] != ".mdcrd":
                    continue
                print("-----------------")
                print("Found mdcrd file at '{}{}'".format(
                    trial_dir,mdcrd))
                # look for .out files matching same name
                out_file = trial_dir + mdcrd[:-6] + ".out"
                if not os.path.isfile(out_file):
                    print("No .out file found at '{}'".format(out_file))
                    continue
                else:
                    print("Found .out file at '{}'".format(out_file))
                # everything's okay, write the file
                s = Struct()
                s.loadtraj(prmtop,trial_dir + mdcrd)
                s.load_OutReader(out_file)
                s.write_low_EPtot([trial_dir + mdcrd[:-6]] * nframes)
                lmin_dict[trial_dir + mdcrd] = s.lowest_EP_idx
                # plot EPtot using matplotlib
                if plot_EPtot:
                    plt.scatter(range(len(s.EPtot)),s.EPtot,marker="*",
                                c="blue")
                    plt.scatter(s.lowest_EP_idx,s.EPtot[s.lowest_EP_idx],
                                marker="o",c="orange")
                    plt.title("EPtot vs. frames for trajectory file: \n" +
                              trial_dir + mdcrd)
                    plt.show()
    return lmin_dict

#end function convert_lmin_frames

def EndoG_HB_analysis(calc_dir_parent=""):
    '''
    Handler for EndoG MD simulation trajectory analysis.
    '''

    # script parent dir abs file path. set default if necessary
    if not calc_dir_parent:
        script_dir_abs = os.path.dirname(os.path.realpath(__file__)) + "/"
    else:
        script_dir_abs = calc_dir_parent

    # loop through pySS.py parent dir, looking for calc dirs containing
    # mdcrd files from md runs
    for d in os.listdir(script_dir_abs):
        prmtop = ""
        calc_dir = script_dir_abs + d + "/"
        # skip if not a calculation dir
        if not is_calc_root(calc_dir):
            continue
        # look for .prmtop
        for pt in os.listdir(calc_dir):
            if "_wat.prmtop" in pt and pt[0] != ".":
                prmtop = script_dir_abs + d + "/" + pt
        # loop through trial directories
        for td in os.listdir(calc_dir + "md/"):
            trial_dir = calc_dir + "md/" + td + "/"
            if not os.path.isdir(trial_dir):
                continue
            # loop through mdcrd files in trial directories
            for mdcrd in os.listdir(trial_dir):
                if mdcrd[-6:] != ".mdcrd":
                    continue
                print("-----------------")
                print("Found mdcrd file at '{}{}'".format(
                    trial_dir,mdcrd))
                # look for .out files matching same name
                out_file = trial_dir + mdcrd[:-6] + ".out"
                if not os.path.isfile(out_file):
                    print("No .out file found at '{}'".format(out_file))
                    continue
                else:
                    print("Found .out file at '{}'".format(out_file))
                # everything's there, do the analysis

                # load trajectory file
                s = Struct()
                s.loadtraj(prmtop,trial_dir + mdcrd)

                # load bond distances, energies, and frequencies
                s.load_EndoG_bonds(False)

                # get EPtot from .out file
                s.load_OutReader('./test.out')
                lowest_EP_idx = s.write_low_EPtot([''] * 10, False)
                s.EPtot = reject_outliers(s.EPtot) / 4.184 # to kcal/mol

                # write to CSV file
                with open('{}_{}_{}.csv'.format(d,td,mdcrd), 'a', newline='') as csvfile:
                    writer = csv.writer(csvfile, delimiter=' ', quotechar='|',
                                        quoting=csv.QUOTE_MINIMAL)
                    writer.writerow(['Structure' + s.hb_dist.keys()])
                    writer.writerow(['Spam', 'Lovely Spam', 'Wonderful Spam'])

                print(s.int_freq)

                # correlate deviances in EPtot to deviation in LIE energy
                # for a specific interaction
                '''
                for key in s.eelec.keys():
                    print(key)
                    print(plot_correl(s.eelec[key], s.EPtot, 'spearman'))
                    try:
                        print(plot_correl(s.eelec[key], s.hb_dist[key], 'pearson'))
                    except KeyError:
                        continue
                    print(s.eelec[key])

                # all points as lines
                #plt.plot(s.hb_dist['SG_O5'],"r-")
                #plt.plot(s.hb_dist['HG_OP2'],"g-")
                #plt.plot(norm_dev(s.eelec['HG_OP2_A']),"b-")
                #plt.plot(norm_dev(s.hb_dist['HG_OP2_A']),"r-")
                #plt.plot(norm_dev(s.EPtot),"k-")
                #plt.plot(norm_dev(s.eelec['HG_OP2_A'])+norm_dev(s.EPtot),"r-")

                # lowest EPtot frames as dots
                #plt.plot(lowest_EP_idx, s.hb_dist['SG_O5'][lowest_EP_idx],"ro")
                #plt.plot(lowest_EP_idx, s.hb_dist['HG_OP2'][lowest_EP_idx],"go")
                #plt.plot(lowest_EP_idx,s.eelec['HG_O5'][lowest_EP_idx],"bo")
                #plt.plot(lowest_EP_idx,s.eelec['HG_OP2'][lowest_EP_idx],"ro")
                #plt.show()
                '''

    # plot and print




#end function EndoG_HB_analysis

def main_handler():
    '''
    Main handling function for pySS.py.
    '''

    # Write 5 lowest EPtot frames of trajectories as rst7 files
    #d = convert_lmin_frames(5, plot_EPtot=True)

    # EndoG H-bond analysis
    EndoG_HB_analysis()


    '''
    the old way using VMD
    #               Get EPtot from out file

    # instantiate .out file reader
    out_reader = OutReader()
    out_reader.load_out("./dev/1_md.out")

    # check if the nsteps in the .out file are as long as the number of
    # frames in the mdcrd file, and load all Etot values into a np array
    if len(out_reader.cycles) == len(be): # equals frames in mdcrd
        EPtot = out_reader.get_EPtot()
    else:
        raise Exception(
            "Output file number of NSTEPs does not equal frames in mdcrd")

    #               Get n frames with lowest EPtot values

    n = 5
    lowest_EP_idx = EPtot.argsort()[:n]

    #              Plot bond energies and Etot

    # scaling factors
    scf_EPtot = 400.
    scf_be = 1.

    # plot Etot and bond energy as normalized deviance from mean
    be_mean = np.multiply(scf_be,be/np.mean(be)) - scf_be
    EPtot_mean = np.multiply(scf_EPtot,EPtot/np.mean(EPtot)) - scf_EPtot

    # all EPtots and bond energies
    plt.scatter(range(len(EPtot)),EPtot_mean,marker="x",c="blue")
    plt.scatter(range(len(EPtot)),be_mean,marker="*",c="red")

    # frames with lowest Etot
    plt.scatter(lowest_EP_idx,EPtot_mean[lowest_EP_idx],marker="o",c="purple")
    plt.scatter(lowest_EP_idx,be_mean[lowest_EP_idx],marker="o",c="orange")

    plt.show()
    '''
#end function main_handler

if __name__ == '__main__':
    main_handler()
#end if
