"# PPRaman" 

This is essentially the development file for the package. To use and download as a package you must:
1. Initialise the repository from github
2. Open R studio session and go to Build -> Build Binary Package, this build the package and will create a zip file containing it in the directory you're working in
3. Unzip the file, and copy into the folder your packages are stored -- probably something like ".../R/Win-library/3.5/" and then a list of all your packages. Drop PPRaman into there.
4. Load a new R session and to load "PPRaman" use the command "library("PPRaman"), this should also load/download other packages PPRaman imports such as hyperSpec etc.
5. PP_Scripts contains pre-written scripts for processing spectra -- note, often these need to be run in line by line order, rather than highlighting and running the whole script. Using the source() function bypasses this issue, putting the script you're using in the source argument.

If developing the package itself, branch from the master -- this branch can be changed to the master in the future if necessary.
