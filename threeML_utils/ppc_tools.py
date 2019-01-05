import numpy as np

#optional 3ML imports
from threeML import DataList, clone_model

import h5py


def compute_ppc(analysis, result, n_sims, file_name=None):
    """FIXME! briefly describe function

    :param analysis: 
    :param result: 
    :param n_sims: 
    :param file_name: 
    :returns: 
    :rtype: 

    """
    
    totals = []
    sim_totals = []
    sim_dls = []

    if file_name is not None:

        database = h5py.File(file_name,'w',libver='latest')
        data_names =[]

        database.attrs['n_sims'] = n_sims
        
        for data in analysis.data_list.values():

            data_names.append(data.name)
            grp = database.create_group(data.name)
            grp.attrs['exposure'] = data.exposure
            grp.create_dataset('ebounds', data = data.response.ebounds, compression = 'lzf')
            grp.create_dataset('obs_counts', data=data.observed_counts, compression='lzf')
            grp.create_dataset('bkg_counts', data=data.background_counts, compression='lzf')
            grp.create_dataset('mask', data=data.mask, compression='lzf')
            
    # select random draws from the posterior
    choices = np.random.choice(len(result.samples.T), replace=False, size=n_sims)
    
    
    # for each posterior
    for j, choice in enumerate(choices):

        params = result.samples.T[choice]
        
        # set the analysis free parameters to the value of the posterior
        for i, (k,v) in enumerate(analysis.likelihood_model.free_parameters.items()):
            v.value = params[i]



            
        # create simulated data sets with these free parameters
        sim_dl = DataList(*[ data.get_simulated_dataset() for data in analysis.data_list.values()])

        # set the model of the simulated data to the model of the simulation
        for i, data in enumerate(sim_dl.values()):

            data.set_model(clone_model(analysis.likelihood_model))

            if file_name is not None:

                grp = database[data_names[i]]
                grp.create_dataset('ppc_counts_%d'%j, data=data.observed_counts, compression='lzf')


                #grp.create_dataset('ppc_bkg_%d'%j, data=data.background_counts, compression='lzf')
                grp.create_dataset('ppc_background_counts_%d'%j, data=data.background_counts, compression='lzf')
                
            
        sim_dls.append(sim_dl)

    if file_name is not None:

        database.close()

        
    return sim_dls


# light="#DCBCBC"
# light_highlight="#C79999"
# mid="#B97C7C"
# mid_highlight="#A25050"
# dark="#8F2727"
# dark_highlight="#7C0000"
# green="#00FF00"
# yellow='#EFFF2D'

# def compute_graphical_ppc(analysis, result, n_sims, lle=False):
    
    
#     n_dets = len(analysis.data_list.values())
    
#     names = [d.name for d in analysis.data_list.values()]
    
#     sim_dls = compute_ppc(analysis, result, n_sims)
    
#     n_rows = int(np.ceil(n_dets/2.))
    
#     row_count = 0
#     col_count = 0
    
    
#     fig, axes = plt.subplots(n_rows,2)
    
    
#     for det in range(n_dets):
        
#         if names[det][0] =='b':
#             name = 'BGO'+names[det][1]
            

#             tmp_rate = 0
#         elif names[det] == 'LLE':

#             name = 'LLE'
#             tmp_rate = 0
        
#         else:
            
#             name = 'NaI'+names[det][1]
            
#             tmp_rate = 5
            
        
        
#         ax = axes[row_count,col_count]
        
        
    
            
            
#         min_rate = 1E-99

#         mean_energy = []
#         rate = []
#         lo_edge = []
#         hi_edge = []


#         this_det = analysis.data_list.values()[det]

#         minus = 0
#         flag = True
#         while(flag):

#             try:
        
#                 ca = this_det._construct_counts_arrays(tmp_rate+minus,False)
#                 print('success')
#                 flag = False
#             except:
#                 print(minus)
#                 minus-=1

                
            
#         data_rebinner = ca['rebinner']
#         data_rate = ca['new_rate']/ca['new_chan_width']
#         data_lo_edge = ca['new_energy_min']
#         data_hi_edge = ca['new_energy_max']

        
#         for x in sim_dls:
    
#             aa = x.values()[det]


#             ca = aa._construct_counts_arrays(min_rate,False, custom_rebinner=data_rebinner)

#             mean_energy.append(ca['mean_energy'])
#             rate.append(ca['new_rate']/ca['new_chan_width'])
#             lo_edge.append(ca['new_energy_min'])
#             hi_edge.append(ca['new_energy_max'])
            
#         rate = np.array(rate)


#         level=95

#         low_95 = np.percentile(rate, 50 - level*0.5, axis=0)
#         high_95 = np.percentile(rate, 50 + level*0.5, axis=0)




#         for lo_rate, hi_rate, lo_e, hi_e in zip(low_95, high_95, lo_edge[0], hi_edge[0]):

#             ax.fill_between([lo_e, hi_e],lo_rate,hi_rate,facecolor=dark,edgecolor=dark_highlight)


#         level=68

#         low = np.percentile(rate, 50 - level*0.5, axis=0)
#         high = np.percentile(rate, 50 + level*0.5, axis=0)




#         for lo_rate, hi_rate, me, lo_e, hi_e in zip(low, high, mean_energy[0], lo_edge[0], hi_edge[0]):

#             ax.fill_between([lo_e, hi_e], lo_rate, hi_rate,facecolor=mid,edgecolor=mid_highlight)


#         for dr, lo_e, hi_e, lob, hib in zip(data_rate, data_lo_edge, data_hi_edge, low_95, high_95):    

#             if (dr>= lob) and (hib>= dr):
#                 color = green
#             else:
#                 color = yellow
            

#             ax.fill_between([lo_e, hi_e], dr-dr*.05,dr+dr*.05 ,facecolor=color,edgecolor='k' ,alpha=0.7)
        
        
        
#         ax.set_xscale('log')
#         ax.set_yscale('log')
        
#         ax.text(0.75, 0.95, name, transform=ax.transAxes,
#           fontsize=6, fontweight='bold', va='top')
        

        
#         if row_count == n_rows -1 or ((row_count == n_rows -2) and (n_dets%2 ==1) and (col_count==1)):
#             ax.set_xlabel("Energy\n(keV)")
            
#         if col_count ==0:
#             ax.set_ylabel("Net rate\n(counts s$^{-1}$ keV$^{-1}$)")
            
            
#         if col_count ==1:
#             col_count =0
#             row_count+=1
        
#         else:
#             col_count+=1
        
#     if n_dets%2 ==1:
        
#         ax = axes[row_count,col_count]
        
#         ax.set_visible(False)
        
        
#     return fig

