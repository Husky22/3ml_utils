import numpy as np

#optional 3ML imports
from threeML import DataList, clone_model

import h5py

import matplotlib.pyplot as plt

def compute_ppc(analysis, result, n_sims, file_name):
    """ 
    Compute a posterior predictive check from a 3ML DispersionLike
    Plugin. The resulting posterior data simulations are stored
    in an HDF5 file which can be read by the PPC class

    :param analysis: 3ML bayesian analysis object 
    :param result: 3ML analysis result
    :param n_sims: the number of posterior simulations to create
    :param file_name: the filename to save to
    :returns: None
    :rtype: 

    """

    with h5py.File(file_name, 'w', libver='latest') as database:

        # first we collect the real data data and save it so that we will not have to
        # look it up in the future

        data_names = []

        database.attrs['n_sims'] = n_sims

        for data in analysis.data_list.values():

            data_names.append(data.name)
            grp = database.create_group(data.name)
            grp.attrs['exposure'] = data.exposure
            grp.create_dataset('ebounds', data=data.response.ebounds, compression='lzf')
            grp.create_dataset('obs_counts', data=data.observed_counts, compression='lzf')
            grp.create_dataset('bkg_counts', data=data.background_counts, compression='lzf')
            grp.create_dataset('mask', data=data.mask, compression='lzf')

        # select random draws from the posterior

        choices = np.random.choice(len(result.samples.T), replace=False, size=n_sims)

        # for each posterior sample

        for j, choice in enumerate(choices):

            # get the parameters of the choice

            params = result.samples.T[choice]

            # set the analysis free parameters to the value of the posterior
            for i, (k, v) in enumerate(analysis.likelihood_model.free_parameters.items()):
                v.value = params[i]

            # create simulated data sets with these free parameters
            sim_dl = DataList(*[data.get_simulated_dataset() for data in analysis.data_list.values()])

            # set the model of the simulated data to the model of the simulation
            for i, data in enumerate(sim_dl.values()):

                # clone the model for saftey's sake
                # and set the model. For now we do nothing with this

                data.set_model(clone_model(analysis.likelihood_model))

                # store the PPC data in the file
                grp = database[data_names[i]]
                grp.create_dataset('ppc_counts_%d' % j, data=data.observed_counts, compression='lzf')
                grp.create_dataset('ppc_background_counts_%d' % j, data=data.background_counts, compression='lzf')
            #sim_dls.append(sim_dl)


class PPC(object):

    def __init__(self, filename):
        """
        Reads a PPC HDF5 created by compute_ppc. This applies to DispersionLike
        data types only. Each detector is read from the file and an associated
        detector attribute is added to the class allowing the user to access the
        observed and PPC information of the detector
        

        :param filename: the file name to read
        :returns: 
        :rtype: 

        """

        # open the file

        with h5py.File(filename, 'r') as f:

            n_sims = f.attrs['n_sims']

            dets = f.keys()

            # go thru all the detectors and grab
            # their data

            for d in dets:

                ppc_counts = []
                ppc_bkg = []
                obs_counts = f[d]['obs_counts'].value
                background_counts = f[d]['bkg_counts'].value
                mask = f[d]['mask'].value

                ebounds = f[d]['ebounds'].value

                exposure = f[d].attrs['exposure']

                # scroll thru the PPCS and build up PPC matrix

                for n in range(n_sims):
                    ppc_counts.append(f[d]['ppc_counts_%d' % n].value.tolist())
                    ppc_bkg.append(f[d]['ppc_background_counts_%d' % n].value.tolist())

                # build a detector object and attach it to the class
                det_obj = PPCDetector(d, obs_counts, background_counts, mask, ebounds, exposure, np.array(ppc_counts),
                                      np.array(ppc_bkg))

                setattr(self, d, det_obj)

            self._n_sims = n_sims
            self._dets = dets
            self._filename = filename

    @property
    def n_sims(self):

        return self._n_sims


class PPCDetector(object):

    def __init__(self, name, obs_counts, obs_background, mask, ebounds, exposure, ppc_counts, ppc_background):
        """
        This is simply a container object that stores the observed and PPC information of each detector for examination

        :param name: 
        :param obs_counts: 
        :param obs_background: 
        :param mask: 
        :param ebounds: 
        :param exposure: 
        :param ppc_counts: 
        :param ppc_background: 
        :returns: 
        :rtype: 

        """

        self._exposure = exposure
        self._obs_counts = obs_counts
        self._obs_background = obs_background
        self._mask = mask
        self._ebounds = ebounds
        self._channel_width = ebounds[1:] - ebounds[:-1]
        self._ppc_counts = ppc_counts
        self._ppc_background = ppc_background
        self._name = name

    @property
    def name(self):
        return self._name

    @property
    def obs_counts(self):
        return self._obs_counts

    @property
    def obs_background(self):
        return self._obs_background

    @property
    def mask(self):
        return self._mask

    @property
    def ebounds(self):
        return self._ebounds

    @property
    def channel_width(self):
        return self._channel_width

    @property
    def exposure(self):
        return self._exposure

    @property
    def ppc_counts(self):
        return self._ppc_counts

    @property
    def ppc_background(self):
        return self._ppc_background

    def plot(self, ax=None, levels=[95, 75, 55], colors=['#ABB2B9', '#566573', '#17202A'], lc='#FFD100', lw=.9):

        assert len(levels) == len(colors), 'Number of levels and number of colors MUST be the same'

        if ax is None:

            fig, ax = plt.subplots()

        else:

            fig = ax.get_figure()


        # compute all the percentiles
            
        ppc_low = []
        ppc_high = []

        for level in levels:

            tmp_low = np.percentile(self._ppc_counts / self._channel_width / self._exposure, 50. - level / 2., axis=0)
            tmp_high = np.percentile(self._ppc_counts / self._channel_width / self._exposure, 50. + level / 2., axis=0)

            ppc_low.append(tmp_low)
            ppc_high.append(tmp_high)

        true_rate = self._obs_counts / self._channel_width / self._exposure

        #colors = [light,mid,dark]

        for j, (lo, hi) in enumerate(zip(ppc_low, ppc_high)):

            for i in range(len(self._ebounds) - 1):
                if self._mask[i]:

                    ax.fill_between([self._ebounds[i], self._ebounds[i + 1]], lo[i], hi[i], color=colors[j])

        n_chan = len(self._ebounds) - 1

        for i in range(len(self._ebounds) - 1):
            if self._mask[i]:

                ax.hlines(true_rate[i], self._ebounds[i], self._ebounds[i + 1], color=lc, lw=lw)

                if i < n_chan - 1:
                    if self._mask[i + 1]:

                        ax.vlines(self._ebounds[i + 1], true_rate[i], true_rate[i + 1], color=lc, lw=lw)

        ax.set_xscale('log')
        ax.set_yscale('log')

        ax.set_ylabel(r'Rate [cnt s$^{-1}$ keV$^{-1}$]')
        ax.set_xlabel(r'Energy [keV]')

        return fig
