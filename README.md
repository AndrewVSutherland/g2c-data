# g2c-data
raw genus 2 curve data for upload to the LMFDB
created by Andrew Booker, Andrew Sutherland, John Voight, and Dan Yasaki.

The main file is gce_1000000_lmfdb.txt, which contains data for genus 2 curves of discriminant up to 1,000,000 that was originally loaded into the LMFDB at the March 2016 LMFDB workshop held at Bristol.

The sage file loadscript.sage loads the data into the LMFDB (assuming you have setup an ssh tunnel), it uses loadcurves.sage to do this.  Before running loadscript.sage, you need to unzip the files gce_1000000_ldata?.bz2 and concatenate them into a single file gce_1000000_ldata.txt.

The auxiliary files gce_1000000_*_lmfdb.txt contain data related to the L-functions of these genus 2 curves (e.g. Hilbert modular forms, products of elliptic curves, etc...)
