LIBS= -lm -L/usr/local/lib -lsundials_cvode -lsundials_nvecserial

CPATHS= -L/usr/local/lib -I/usr/local/include

ALL : prog

prog :	prog.c
	gcc $(CPATHS) -O4 -Wall -o prog prog.c evolution_sun_mod.c $(LIBS) 

clean :
	rm prog
