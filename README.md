# g2c-data
raw genus 2 curve data for upload to the LMFDB
created by Andrew Booker, Andrew Sutherland John Voight, and Dan Yasaki

The main file is ucd_1000000_lmfdb.txt, which contains data for genus 2 curves of discriminant up to 1,000,000 that was originally loaded into the LMFDB at the March 2015 LMFDB workshop held at University College Dublin UCD

The sage script loadscript.sage loads the data into the LMFDB (assuming you have setup an ssh tunnel), it uses loadcurves.sage to do this.

The auxiliary files ucd_1000000_*_lmfdb.txt contain data related to various objects related to the L-functions of these genus 2 curves (e.g. Hilbert modular forms, Weil-restrictions of elliptic curves over quadratic fields, etc...)

IMPORTANT: this repository does not include the file Ldata.txt which contains precomputed data for each of the 64,803 L-functions associated to these curves (even compressed it is over 9 gigabytes, much too large to put in github).  The sage script will work without this file present, but you will then be missing all the L-function data.  Contact Andrew Booker for access to the L-function data file.
