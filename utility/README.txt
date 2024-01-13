Use this script to automatically download the pesudopotentials from http://www.pseudo-dojo.org/ for use in your code for QE calculations. 
Note, you only need to download them if you do not already have the appropriate pseudopotential files for your calculations.
We however, recommend this potential because it is norm convserving.
We have provided the file "pseudo.info" where you can specify some options. If you want to download all the latest and available PBE pseudopotentials from the website, just remove this file to use the default. It uses chrome brower. We recommend this option. Note, the downloading is done automatically. You do not have to do anything.
 

To use utility, you can run it in two ways:
1. You can install it simple as `pip install .`
Afterwards, simply run "qepotential" in your terminal. Code is optimized to use Google Chrome


OR 

2. python download_qe_pseudo.py
