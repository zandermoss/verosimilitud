#
# 'make depend' uses makedepend to automatically generate dependencies 
#			   (dependencies are added to end of Makefile)
# 'make'		build executable file 'mycc'
# 'make clean'  removes all .o and executable files
#

# define the C compiler to use
#CC = gcc
CC = clang++
#CC = g++

# define any compile-time flags
CFLAGS = -Wall -g -fPIC -std=c++11

# define any directories containing header files other than /usr/include
#
INCLUDES = -I ./inc -I/usr/include/python2.7 -I/usr/local/hdf5/include -I/home/pinkpig/src/dlib-18.18

# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib I'd specify
#   their path using -Lpath, something like:
LFLAGS = -L/usr/local/hdf5/lib

# define any libraries to link into executable:
#   if I want to link in libraries (libx.so or libx.a) I use the -llibname 
#   option, something like (this will link in libmylib.so and libm.so:
LIBS = -lpython2.7 -lm -lhdf5 -lhdf5_cpp  

# define the C source files
#SRCS = src/Verosimilitud.cpp src/Tools.cpp
SRCS = src/*

# define the C object files 
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#		 For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
OBJS = $(SRCS:.c=.o)

# define the executable file 
#SHARED = libverosimilitud.so
SHARED = runolv

#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

.PHONY: depend clean

all:	$(SHARED)
		@echo  Caller has been compiled

#$(SHARED): $(OBJS) 
#		$(CC) $(CFLAGS) -shared $(INCLUDES) -o $(SHARED) $(OBJS) $(LFLAGS) $(LIBS)
$(SHARED): $(OBJS) 
		$(CC) $(CFLAGS) $(INCLUDES) -o $(SHARED) $(OBJS) $(LFLAGS) $(LIBS)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
.c.o:
		$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
		$(RM) *.o *~ $(SHARED)

depend: $(SRCS)
		makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
