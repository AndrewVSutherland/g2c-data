DIR=/usr/local

plot: plot.c bessel_arb.o fft_acb.o wdata.c
	gcc -O -I$(DIR)/include -I$(DIR)/include/flint -o plot plot.c bessel_arb.o fft_acb.o -L$(DIR)/lib -larb -lflint $(DIR)/lib/libmpfr.a $(DIR)/lib/libsmalljac.a $(DIR)/lib/libff_poly.a $(DIR)/lib/libgmp.a -lm

bessel_arb.o: bessel_arb.c wdata.c
	gcc -O -I$(DIR)/include -I$(DIR)/include/flint -c bessel_arb.c

fft_acb.o: fft_acb.c
	gcc -O -I$(DIR)/include -I$(DIR)/include/flint -c fft_acb.c

clean:
	rm -f bessel_arb.o fft_acb.o plot
