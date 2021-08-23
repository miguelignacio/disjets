(1) Use the Process.ipynb notebook to read in the data files from the DESY batch and produce pickle files with only the events and information we need.  The processed datasets are about 10 GB each and there are about a dozen of them.

(2) Run the unfolding using unfold_fullstats.py.  Before you run this, you need to (once) do

mkdir models log_files storage_files storage_plots

You can then run it like

`python unfold_fullstats.py <mc> <syst> <GPU>`
    
where <mc> is Rapgap or Django, <syst> is nominal, syst_0, syst_1, syst_5, sys_7, or sys_11, and <GPU> is the GPU you want to use (e.g. 0, 1, 2, 3 if the machine you are using has four GPUs).  All of the NNs are saved in models and timing information is saved in the log files.  A number of useful diagnostic plots are saved in storage_plots and the final results (as numpy array histogram bin contents) are stored in storage_files.
    
(3) Run the bootstrapping for the data stats using unfold_fullstats_boot.py.  You can run it like
    
python unfold_fullstats_boot.py Rapgap nominal 0 <boostrap>
    
where <bootstrap> is a number greater than 0 that sets the seed for the bootstrapping.
    
(4) Run the Uncerts.ipynb notebook to produce the error bands and uncertainty breakdown plots for the final results.
    
(5) Run finalplots.ipynb notebook to produce the final normalized differntial cross section plots.
