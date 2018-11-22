"# PPRaman" 

This is essentially the development file for the package. To use and download as a package you must:
1. Initialise the repository from github
2. Open R studio session and go to Build -> Build Binary Package, this build the package and will create a zip file containing it in the directory you're working in
3. Unzip the file, and copy into the folder your packages are stored -- probably something like ".../R/Win-library/3.5/" and then a list of all your packages. Drop PPRaman into there.
4. Load a new R session and to load "PPRaman" use the command "library("PPRaman"), this should also load/download other packages PPRaman imports such as hyperSpec etc.
5. The repository PP_Scripts contains pre-written scripts for processing spectra -- A quick start that mimicks the matlab code is "PreProcessing_Plots_PCA_Seq.R", however this relies heavily on the file name, so adhoc runs wont work. There will be a script for adhoc tasks in PP_Scripts soon.
6. Alternatively, to run preprocessing without generating plots/labels (e.g for building model) using the function process_all_spec(), without arguments - the result is a list with the first element containing corrected spectra, and the second contains patient and date info.

If developing the package itself, branch from the master -- this branch can be changed to the master in the future if necessary.
