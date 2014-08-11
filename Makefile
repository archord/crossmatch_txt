PERFIX = ../../lib
WCSTOOL = ${PERFIX}/wcstools-3.8.5
ASTROMETRY = ${PERFIX}/astrometry.net-0.38
CFITSIO = ${PERFIX}/cfitsio

#WCSTOOL = /home/xy/program/lib/wcstools-3.8.5/libwcs
#ASTROMETRY = /home/xy/program/lib/astrometry
#CFITSIO = /disk1/opt/basic/prog_os_20110130_from_pc49
POSTGRESQL = /usr


ALL:crossmatch

crossmatch:
	g++ ${IDIR} function.c main.c -o crossmatch  ${LDIR} ${LIBS}

clean:
	rm -rf *.o crossmatch

