#
# This Makefile was generated by Code::Blocks IDE.
#

SRCS_f90d1 = \
main.f90 

OBJS_f90d1 = \
main.o 

SRC_DIR_f90d1 = 
OBJS_DIR = obj/Release/
EXE_DIR = bin/Release/

EXE = Generate2DVTKFrom3D
FC = ifort
LD = ifort
IDIR = 
CFLAGS = -integer-size 64  -J$(OBJS_DIR) $(IDIR)
LFLAGS = 
LIBS = 

VPATH = $(SRC_DIR_f90d1):$(OBJS_DIR)
OBJS = $(addprefix $(OBJS_DIR), $(OBJS_f90d1))

all : $(EXE)

$(EXE) : $(OBJS_f90d1)
	@mkdir -p $(EXE_DIR)
	$(LD) -o $(EXE_DIR)$(EXE) $(OBJS) $(LFLAGS) $(LIBS)

$(OBJS_f90d1):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_f90d1)$(@:.o=.f90) -o $(OBJS_DIR)$@

clean :
	rm -f $(OBJS_DIR)*.*
	rm -f $(EXE_DIR)$(EXE)

# Dependencies of files
main.o: \
    main.f90

