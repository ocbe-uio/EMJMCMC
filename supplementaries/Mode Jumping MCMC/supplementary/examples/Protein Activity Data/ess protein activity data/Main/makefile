# This file is part of ESS++.
#      Copyright (c) Marc Chadeau (m.chadeau@imperial.ac.uk)
#                    Leonardo Bottolo (l.bottolo@imperial.ac.uk)
#                    David Hastie (d.hastie@imperial.ac.uk)  
#      2010
#      
# ESS++ is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESS++ is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ESS++.  If not, see <http://www.gnu.org/licenses/>.
#

LIBDIR		= ../Routines/
LIBDIR2		= ../Classes/

CC 		= gcc
CXX 		= g++
IFLAGS		= -I$(LIBDIR)
LFLAGS		= -L$(LIBDIR)
LIBS		= -lgsl -lgslcblas

CXXFLAGS	= 

CPPFLAGS	= -W -Wall --pedantic -O3  $(IFLAGS)
#CPPFLAGS	= -W -Wall -g -pg -O0 --pedantic $(IFLAGS)
#CPPFLAGS	= -W -Wall -O4 -g $(IFLAGS)
# Use the compile option -pg to time C++ code.



all:		ESS .depend 

clean:		
		rm -f *.o
		rm -f .depend
		rm -f ESS
		rm -f *~
		rm -f ../Routines/*.o
		rm -f ../Routines/*~
		rm -f ../Classes/*.o
		rm -f ../Classes/*~
		rm -f ../.depend
dep .depend:
	$(CXX) $(CPPFLAGS) -MM ESS.cc >>.depend

infiles1	= ESS.o  $(LIBDIR)dyn_name.o $(LIBDIR)matrix_handling.o  $(LIBDIR)post_processing.o $(LIBDIR)rand.o $(LIBDIR2)String_Matrices.o $(LIBDIR2)Int_Matrices.o $(LIBDIR2)Double_Matrices.o  $(LIBDIR2)Double_Matrices_cont.o  $(LIBDIR)regression.o $(LIBDIR)cond_post.o $(LIBDIR)xml_file_read.o $(LIBDIR2)Prior_param.o $(LIBDIR)moves.o $(LIBDIR2)Temperatures.o $(LIBDIR2)g_AdMH.o $(LIBDIR2)Move_monitor.o $(LIBDIR2)DR.o $(LIBDIR2)CM.o $(LIBDIR2)Int_Matrices_var_dim.o

ESS:	$(infiles1)
		$(CXX) $(CPPFLAGS) $(LFLAGS) -o ESS $(infiles1) $(LIBS)

