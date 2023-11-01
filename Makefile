################################################################################
##                                                                            ##
##    Copyright 2023 Amol Upadhye                                             ##
##                                                                            ##
##    This file is part of Cosmic-Enu.                                        ##
##                                                                            ##
##    Cosmic-Enu is free software: you can redistribute it and/or modify      ##
##    it under the terms of the GNU General Public License as published by    ##
##    the Free Software Foundation, either version 3 of the License, or       ##
##    (at your option) any later version.                                     ##
##                                                                            ##
##    Cosmic-Enu is distributed in the hope that it will be useful,           ##
##    but WITHOUT ANY WARRANTY; without even the implied warranty of          ##
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           ##
##    GNU General Public License for more details.                            ##
##                                                                            ##
##    You should have received a copy of the GNU General Public License       ##
##    along with Cosmic-Enu.  If not, see <http:##www.gnu.org/licenses/>.     ##
##                                                                            ##
################################################################################

CFLAGS=-O3 -fopenmp
PATHS=
LIBS=-lm -lgsl -lgslcblas

Cosmic-Enu: Cosmic-Enu.c AU_D2_D2ehNorm_deciles.h AU_fftgrid.h AU_ncint.h AU_fftm_emu_data.h Makefile
	gcc Cosmic-Enu.c -o Cosmic-Enu $(CFLAGS) $(PATHS) $(LIBS)

clean:
	$(RM) Cosmic-Enu

