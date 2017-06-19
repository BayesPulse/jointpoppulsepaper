SRC = jt_deconvolution_main.c jt_format_data.c jt_mcmc.c jt_hash.c\
randgen.c jt_birthdeath.c cholesky.c normal.c
OBJ = jt_deconvolution_main.o jt_format_data.o jt_mcmc.o jt_hash.o\
randgen.o jt_birthdeath.o cholesky.o normal.o
CFLAGS =  -O3
LFLAGS =
LIBS = -lm

decon : ${OBJ}
	gcc -o jtdecon ${OBJ} ${LIBS}

${OBJ} : ${SRC}
	gcc -c ${CFLAGS}  ${SRC}

{OBJ} : jt_deconvolution_main.h jt_birthdeath.h jt_hash.h jt_mcmc.h randgen.h

clean:
	rm ${OBJ}
