import threeML
from threeML import threeML_config
from threeML.io.plotting.cmap_cycle import cmap_intervals
from threeML.io.plotting.data_residual_plot import ResidualPlot

import numpy as np
from threeML_utils.colors import Colors



NO_REBIN = 1E-99


def scale_colour(self, colour, scalefactor):  # pragma: no cover
        if isinstance(colour, np.ndarray):
            r, g, b = colour[:3] * 255.0
        else:
            hexx = colour.strip('#')
            if scalefactor < 0 or len(hexx) != 6:
                return hexx
            r, g, b = int(hexx[:2], 16), int(hexx[2:4], 16), int(hexx[4:], 16)
        r = self._clamp(int(r * scalefactor))
        g = self._clamp(int(g * scalefactor))
        b = self._clamp(int(b * scalefactor))
        return "#%02x%02x%02x" % (r, g, b)

def display_posterior_model_counts(bayesian_analysis, result = None ,thin=100, shade=True, q_level=68, gradient=0.6, data=(), **kwargs):
    """

    Display the fitted model count spectrum of one or more Spectrum plugins

    NOTE: all parameters passed as keyword arguments that are not in the list below, will be passed as keyword arguments
    to the plt.subplots() constructor. So for example, you can specify the size of the figure using figsize = (20,10)

    :param args: one or more instances of Spectrum plugin
    :param min_rate: (optional) rebin to keep this minimum rate in each channel (if possible). If one number is
    provided, the same minimum rate is used for each dataset, otherwise a list can be provided with the minimum rate
    for each dataset
    :param data_cmap: (str) (optional) the color map used to extract automatically the colors for the data
    :param model_cmap: (str) (optional) the color map used to extract automatically the colors for the models
    :param data_colors: (optional) a tuple or list with the color for each dataset
    :param model_colors: (optional) a tuple or list with the color for each folded model
    :param data_color: (optional) color for all datasets
    :param model_color: (optional) color for all folded models
    :param show_legend: (optional) if True (default), shows a legend
    :param step: (optional) if True (default), show the folded model as steps, if False, the folded model is plotted
    :param model_subplot: (optional) axe(s) to plot to for overplotting
    with linear interpolation between each bin
    :return: figure instance


    """

    # If the user supplies a subset of the data, we will use that

    if not data:

        data_keys = bayesian_analysis.data_list.keys()

    else:

        data_keys = data

    # Now we want to make sure that we only grab OGIP plugins

    new_data_keys = []

    for key in data_keys:

        # Make sure it is a valid key
        if key in bayesian_analysis.data_list.keys():

            if isinstance(bayesian_analysis.data_list[key], threeML.plugins.SpectrumLike.SpectrumLike):

                new_data_keys.append(key)

            else:

                custom_warnings.warn("Dataset %s is not of the SpectrumLike kind. Cannot be plotted by "
                                     "display_spectrum_model_counts" % key)

    if not new_data_keys:
        RuntimeError(
            'There were no valid SpectrumLike data requested for plotting. Please use the detector names in the data list'
        )

    data_keys = new_data_keys

    # default settings

    # Default is to show the model with steps
    step = True

    data_cmap = threeML_config['ogip']['data plot cmap']    # plt.cm.rainbow
    model_cmap = threeML_config['ogip']['model plot cmap']    # plt.cm.nipy_spectral_r

    # Legend is on by default
    show_legend = False

    show_residuals = False

    # Default colors

    data_colors = cmap_intervals(len(data_keys), data_cmap)
    model_colors = cmap_intervals(len(data_keys), model_cmap)

    # Now override defaults according to the optional keywords, if present

    model_kwargs = dict(alpha=0.05, zorder=-5000)
    if 'model_kwargs' in kwargs:

        model_kwargs_tmp = kwargs.pop('model_kwargs')

        for k,v in model_kwargs_tmp.items():

            model_kwargs[k] = v

        

    data_kwargs = dict(
        alpha=1.,
        fmt=threeML_config['residual plot']['error marker'],
        markersize=threeML_config['residual plot']['error marker size'],
        linestyle='',
        elinewidth=threeML_config['residual plot']['error line width'],
        capsize=0,
        zorder=-1)

            
    if 'data_kwargs' in kwargs:

        data_kwargs_tmp = kwargs.pop('data_kwargs')

        for k,v in data_kwargs_tmp.items():

            data_kwargs[k] = v
        


    if 'show_data' in kwargs:

        show_data = bool(kwargs.pop('show_data'))

    else:

        show_data = True

    if 'show_legend' in kwargs:
        show_legend = bool(kwargs.pop('show_legend'))

    if 'show_residuals' in kwargs:
        show_residuals = bool(kwargs.pop('show_residuals'))

    if 'step' in kwargs:
        step = bool(kwargs.pop('step'))

    if 'min_rate' in kwargs:
        min_rate = kwargs.pop('min_rate')

        # If min_rate is a floating point, use the same for all datasets, otherwise use the provided ones

        try:

            min_rate = float(min_rate)

            min_rates = [min_rate] * len(data_keys)

        except TypeError:

            min_rates = list(min_rate)

            assert len(min_rates) >= len(
                data_keys), "If you provide different minimum rates for each data set, you need" \
                            "to provide an iterable of the same length of the number of datasets"

    else:

        # This is the default (no rebinning)

        min_rates = [NO_REBIN] * len(data_keys)

    if 'data_cmap' in kwargs:
        data_cmap = plt.get_cmap(kwargs.pop('data_cmap'))
        data_colors = cmap_intervals(len(data_keys), data_cmap)

    if 'model_cmap' in kwargs:
        model_cmap = kwargs.pop('model_cmap')
        model_colors = cmap_intervals(len(data_keys), model_cmap)

    if 'data_colors' in kwargs:
        data_colors = kwargs.pop('data_colors')

        assert len(data_colors) >= len(data_keys), "You need to provide at least a number of data colors equal to the " \
                                                   "number of datasets"

    elif 'data_color' in kwargs:

        data_colors = [kwargs.pop('data_color')] * len(data_keys)

    if 'model_colors' in kwargs:
        model_colors = kwargs.pop('model_colors')

        assert len(model_colors) >= len(
            data_keys), "You need to provide at least a number of model colors equal to the " \
                        "number of datasets"

    ratio_residuals = False
    if 'ratio_residuals' in kwargs:
        ratio_residuals = bool(kwargs['ratio_residuals'])

    elif 'model_color' in kwargs:

        model_colors = [kwargs.pop('model_color')] * len(data_keys)

    if 'model_labels' in kwargs:
        model_labels = kwargs.pop('model_labels')

        assert len(model_labels) == len(data_keys), 'you must have the same number of model labels as data sets'

    else:

        model_labels = ['%s Model' % bayesian_analysis.data_list[key]._name for key in data_keys]

    #fig, (ax, ax1) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [2, 1]}, **kwargs)

    residual_plot = ResidualPlot(show_residuals=show_residuals, **kwargs)

    if show_residuals:

        axes = [residual_plot.data_axis, residual_plot.residual_axis]

    else:

        axes = residual_plot.data_axis



    # extract the samples
    if result is None:
        samples = bayesian_analysis.results.samples.T[::thin]

    else:

        samples = result.samples.T[::thin]
        
     


    if shade:

        color_config = Colors()

        shade_y = []
        shade_x = []
    
    for params in samples:

        for i, (k, v) in enumerate(bayesian_analysis.likelihood_model.free_parameters.items()):

            v.value = params[i]

        # first with no data
        if shade:
            per_det_y = []
            per_det_x = []

        for key, data_color, model_color, min_rate, model_label in zip(data_keys, data_colors, model_colors,
                                                                       min_rates, model_labels):



            # NOTE: we use the original (unmasked) vectors because we need to rebin ourselves the data later on

            data = bayesian_analysis.data_list[key]    # type: threeML.plugins.SpectrumLike.SpectrumLike

            if not shade:

                data.display_model(
                    data_color=data_color,
                    model_color=model_color,
                    min_rate=min_rate,
                    step=step,
                    show_residuals=False,
                    show_data=False,
                    show_legend=show_legend,
                    ratio_residuals=ratio_residuals,
                    model_label=None,
                    model_subplot=axes,
                    data_kwargs=data_kwargs,
                    model_kwargs=model_kwargs)
            else:

                # this is private for now
                rebinned_quantities = data._construct_counts_arrays(min_rate, ratio_residuals)

                if step:

                    pass

                else:

                    y = (rebinned_quantities['expected_model_rate'] / rebinned_quantities['chan_width'])[data.mask]


                    x =  np.mean([rebinned_quantities['energy_min'],
                         rebinned_quantities['energy_max']], axis=0)[data.mask]

                    per_det_y.append(y)
                    per_det_x.append(x)

        if shade:
            shade_y.append(per_det_y)
            shade_x.append(per_det_x)

    if shade:

        # convert to per detector
        shade_y = np.array(shade_y).T
        shade_x = np.array(shade_x).T
        model_kwargs.pop('zorder')
        for key, data_color, model_color, min_rate, model_label, x, y, in zip(data_keys, data_colors, model_colors, min_rates,
                                                                              model_labels, shade_x, shade_y):


        
            
            # we have to do a little reshaping because... life

            y = np.array([yy.tolist() for yy in y ])

            q_levels = np.atleast_1d(q_level)
            q_levels.sort()

            scale = 1.
            zorder = -100
            for level in q_levels:

                color = color_config.format(model_color)
                color_scale = color_config.scale_colour(color,scale)
            
                # first we need to get the quantiles along the energy axis
                low = np.percentile(y, 50 - level*0.5, axis=0)
                high = np.percentile(y, 50 + level*0.5, axis=0)
                
                residual_plot.data_axis.fill_between(x[0], low, high, color=color_scale,zorder=zorder, **model_kwargs)

                scale *= gradient

                zorder-=1
        

        
                
    for key, data_color, model_color, min_rate, model_label in zip(data_keys, data_colors, model_colors, min_rates,
                                                                   model_labels):

        # NOTE: we use the original (unmasked) vectors because we need to rebin ourselves the data later on

        data = bayesian_analysis.data_list[key]    # type: threeML.plugins.SpectrumLike.SpectrumLike

        data.display_model(
            data_color=data_color,
            model_color=model_color,
            min_rate=min_rate,
            step=step,
            show_residuals=False,
            show_data=show_data,
            show_legend=show_legend,
            ratio_residuals=ratio_residuals,
            model_label=model_label,
            model_subplot=axes,
            model_kwargs=dict(alpha=0.),    # no model
            data_kwargs=data_kwargs)

    return residual_plot.figure
