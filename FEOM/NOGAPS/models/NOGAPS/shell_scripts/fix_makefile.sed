#!/usr/bin/sed -f
#
# This code may (or may not) be part of the NOGAPS distribution,
# So it is not protected by the DART copyright agreement.
#
# DART $Id: fix_makefile.sed 111 2010-06-09 21:55:44Z thoar $
#
# Run this on the NOGAPS makefile, e.g.:
# $ ./fix_makefile.sed < makefile > makefile.lib
# $ make -f makefile.lib all

# Separate the objects into the one holding the main routine and everything else
/^OBJS *=/i\
MAIN_OBJ = forecast.o
/^OBJS *=/,/^$/s/forecast.o//

# Rename for clarity - separate the library from the main application
s/OBJS/LIB_OBJS/g
s/TARGET/TARGET_EXE/g

# Add a target for building the library
/TARGET_EXE=/i\
TARGET_LIB=libnogaps.a
/^\$(TARGET_EXE)/i\
all: $(TARGET_EXE) $(TARGET_LIB)\
 \
$(TARGET_LIB): module $(LIB_OBJS)\
	ar -rc $(TARGET_LIB) $(LIB_OBJS)

# Modify the executable target to link in the library
/^\$(TARGET_EXE)/,/^$|^#/s/\$(LIB_OBJS)/& $(MAIN_OBJ)/


