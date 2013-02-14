import numpy as np
import figure
import grid
import eztables
import output
import sensetest

def gen_log_likelihoods(fakein, Nstars, ntests, err, filters, outdir, outnames):
    for i in range(len(filters)):
        sensetest.fit_model_seds_pytables(fakein, Nstars, ntests, stellar_filename,err=err,filters = filters[i], outname = outdir+outnames[i]+'.hf5')

def gen_summary_tables(gext_grid, keys, outdir, outnames):
    for i in range(len(outnames)):
        summary_table = output.expectation_values_pytables(gext_grid,keys=keys,filepath = outdir+outnames[i]+'.hf5')
        summary_table.write(outdir+'summary_'+outnames[i] + '.fits',clobber=True,append=False)
        del summary_table

def plot_keys(keys, outdir, outnames):
    summary_table = eztables.Table()
    for i in range(len(outnames)):
        summary_table.read(outdir+'summary_'+outnames[i]+'.fits')
        for j, key in enumerate(keys):
            figure.figure(j)
            rec_vals = summary_table.data[key+'_recovered'] #recovered values
            true_vals = summary_table.data[key] #true values
            rec_vals = rec_vals - true_vals #offset
            uniq_vals = np.unique(true_vals) #unique true values
            avg_offs = np.zeros(uniq_vals.size) #Mean of recovered params for given input param
            for k in range(avg_offs.size):
                sel = np.where(true_vals==uniq_vals[k])
                avg_offs[k] = rec_vals[sel].mean()
            figure.cplot(uniq_vals, avg_offs,label=(outnames[i]).replace('_',' '),marker='+',markersize=10)
    for j, key in enumerate(keys):
        figure.figure(j)
        figure.xlabel(key.replace('_',''))
        figure.ylabel(key.replace('_','')+' (out - in)')
        figure.legend(loc='best')
        figure.cplot(figure.xlim(),[0,0],linestyle='dotted',color='black')
    figure.show()
            
if __name__ == '__main__':
    #Pick Nstars fake stars within delta_logParam of logParam
    Nstars = 300
    ntests = 2
    stellar_filename = 'fbump_only.fits'

    gext = grid.FileSEDGrid(stellar_filename)
    ten_pc_to_andromeda = 1.609e-10
    gext.seds *=  ten_pc_to_andromeda
    logM = 1.0
    logT= 3.6
    delt_logM = 0.1
    delt_logT = 0.1
    subsel = sensetest.getFakeInds(gext,logM,logT,delt_logM,delt_logT)
    fakein = subsel[np.random.randint(low=0,high=subsel.size,size=Nstars)]
    del gext, subsel
    
    err = 0.05
    keys = ['Av','Rv','logT','logM','logA','logL']
    filters = [np.arange(6),[0,1],[2,3]] #All six PHAT filters, F275W and F336W, F475W and F814W
    outdir = 'Tests/'
    outnames = ['all_phat','F275_336','F475_814']

    #gen_log_likelihoods(fakein, Nstars, ntests, err, filters, outdir, outnames)
    
    #gext = grid.FileSEDGrid(stellar_filename) 
    #gext_grid = gext.grid
    #del gext
    #gen_summary_tables(gext_grid, keys, outdir, outnames) #generate and save summary tables
    #del gext_grid
    
    plot_keys(keys,outdir,outnames)
    

