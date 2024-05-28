#!/bin/bash

helptest=$(cat <<-END
               EXOPLASIM CONFIGURATION SCRIPT
               
    Usage:
        ./configure.sh -v 3.9
    
    Options:
        -v:   Version of Python that will be used, with major and minor version numbers.
        
        -h:   Output this test
END
)

pyversion=3
while getopts "v:h" opt; do
    case $opt in
        v)
            pyversion=$OPTARG
            ;;
        h)
            echo "$helptext"
            exit 0
            ;;
        \?)
            pyversion=3
            echo "No Python version specified; using default Python 3"
            ;;
    esac
done

rm -f *.x most_* F90_*

# check for FORTRAN-90 compiler
# put your favourite compiler in front if you have more than one
# WARNING: Portland Group Compiler pgf is reported to be buggy

for MOST_F90 in ifort gfortran xlf sunf90 nagfor f90 f95 pgf90 g95 af90 NO_F90
do
   `hash 2>/dev/null $MOST_F90`
   if [ $? = 0 ] ; then break ; fi
done

if [ $MOST_F90 != "NO_F90" ] ; then
   F90_PATH=`which $MOST_F90`
   echo "Found FORTRAN-90 compiler     at: $F90_PATH"
   if [ $MOST_F90 = "sunf90" ] ; then
      MOST_F90_OPTS="-O3"
      DEBUG_F90_OPTS="-g -C -ftrap=common"
      echo > most_precision_options "MOST_PREC=-r8"
      echo > most_precision_optionsx "MOST_PREC=-r"
   elif [ $MOST_F90 = "ifort" ] ; then
      MOST_F90_OPTS="-O3 -cpp"
      DEBUG_F90_OPTS="-g -C -cpp -fpe0 -traceback"
      echo > most_precision_options "MOST_PREC=-r8"
      echo > most_precision_optionsx "MOST_PREC=-r"
      PREC_OPTS="-r8"
   elif [ $MOST_F90 = "nagfor" ] ; then
      MOST_F90_OPTS="-O -kind=byte"
      DEBUG_F90_OPTS="-g -C -fpe0 -traceback -kind=byte"
   elif [ $MOST_F90 = "gfortran" ] ; then
      # get major and minor version number
      ver=$( echo "__GNUC__.__GNUC_MINOR__.__GNUC_PATCHLEVEL__" | gfortran -E -P - | sed 's/ //g' )
      GFMAJ=`echo $ver | head -1 | sed -e 's/.*)//' | awk 'BEGIN {FS="."}{print $1}'`
      GFMIN=`echo $ver | head -1 | sed -e 's/.*)//' | awk 'BEGIN {FS="."}{print $2}'`
      GFVER=`expr 100 '*' $GFMAJ + $GFMIN`
      echo "gfortran version " $GFMAJ.$GFMIN
      # flags for gfortran version >= 4.9 [ -ffpe-summary ]
      if [ "$GFVER" -ge "409" ] ; then
         MOST_F90_OPTS="-O3 -cpp -fcheck=all -ffixed-line-length-132 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=zero"
         DEBUG_F90_OPTS="-ggdb -Og -cpp -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -fcheck=all -finit-real=snan"
         echo > most_precision_options "MOST_PREC=-fdefault-real-8"
         echo > most_precision_optionsx "MOST_PREC=-fdefault-real-"
         PREC_OPTS="-fdefault-real-8"
      # flags for gfortran version 4.5 - 4.8
      elif [ "$GFVER" -ge "405" ] ; then
         MOST_F90_OPTS="-O3 -cpp -fcheck=all -ffixed-line-length-132 -finit-real=zero -ffpe-trap=invalid,zero,overflow"
         DEBUG_F90_OPTS="-ggdb -Og -cpp -ffpe-trap=invalid,zero,overflow -fcheck=all -finit-real=snan"
         echo > most_precision_options "MOST_PREC=-fdefault-real-8"
         echo > most_precision_optionsx "MOST_PREC=-fdefault-real-"
         PREC_OPTS="-fdefault-real-8"
      # flags for gfortran version 4.2, 4.3 and 4.4
      else
         MOST_F90_OPTS="-O3 -cpp -finit-real=zero"
         DEBUG_F90_OPTS="-g -cpp -ffpe-trap=invalid,zero,overflow -fbounds-check"
      fi
   else
      MOST_F90_OPTS="-O"
      DEBUG_F90_OPTS="-g"
   fi
   echo  > most_compiler "MOST_F90=$MOST_F90"
   echo >> most_compiler "MOST_F90_OPTS=$MOST_F90_OPTS"
   echo >> most_compiler "MPIMOD=mpimod_stub"
   echo >> most_compiler "UTILMOD=utilities_stub"
   echo  > most_debug_options "MOST_F90_OPTS=$DEBUG_F90_OPTS"
   
   echo  > plasim/bld/most_snow_build4 "#!/bin/bash"
   echo >> plasim/bld/most_snow_build4 "$MOST_F90 $MOST_F90_OPTS newsnow.f90 -o ../bin/newsnow.x"
   
   echo  > plasim/bld/most_snow_build8 "#!/bin/bash"
   echo >> plasim/bld/most_snow_build8 "$MOST_F90 $MOST_F90_OPTS $PREC_OPTS newsnow.f90 -o ../bin/newsnow.x"
   
   echo  > plasim/bld/most_ice_build4 "#!/bin/bash"
   echo >> plasim/bld/most_ice_build4 "$MOST_F90 $MOST_F90_OPTS buildice.f90 -o ../bin/buildice.x"
   
   echo  > plasim/bld/most_ice_build8 "#!/bin/bash"
   echo >> plasim/bld/most_ice_build8 "$MOST_F90 $MOST_F90_OPTS $PREC_OPTS buildice.f90 -o ../bin/buildice.x"
   
else
   echo "****************************************************"
   echo "* Sorry - didn't find any FORTRAN-90 compiler      *"
   echo "* Please install one (e.g. gfortran or g95)        *"
   echo "* or update list of known compilers in this script *"
   echo "****************************************************"
   exit 1
fi

# We no longer need Xlib
# #check for Xlib
#   if [ -e "/usr/X11/lib64/libX11.so" ] ; then
#    XLIB_PATH="/usr/X11/lib64"
#    XINC_PATH="/usr/X11/lib64/include"
# elif [ -e "/usr/lib64/libX11.so" ] ; then
#    XLIB_PATH="/usr/lib64"
#    XINC_PATH="/usr/lib64/include"
# elif [ -e "/usr/X11/lib" ] ; then
#    XLIB_PATH="/usr/X11/lib"
#    XINC_PATH="/usr/X11/include"
# elif [ -e "/usr/lib/X11" ] ; then
#    XLIB_PATH="/usr/lib/X11"
#    XINC_PATH="/usr/lib/X11/include"
# elif [ -e "/opt/X11" ] ; then
#    XLIB_PATH="/opt/X11"
#    XINC_PATH="/opt/X11/include"
# fi
# if [ -n ${XLIB_PATH} ] ; then
#    echo "Found Xlib (X11)              at: $XLIB_PATH"
#    GUILIB="-L$XLIB_PATH -lX11"
# else
#    echo "**********************************************"
#    echo "* Didn't find Xlib (X11) at standard paths   *"
#    echo "* Hopefully the compiler will find it itself *"
#    echo "**********************************************"
#    GUILIB=" -lX11"
# fi

#check for C compiler
# put your favourite compiler in front if you have more than one

for MOST_CC in gcc cc gxlc suncc cc NO_CC
do
   `hash 2>/dev/null $MOST_CC`
   if [ $? = 0 ] ; then break ; fi
done

if [ $MOST_CC != "NO_CC" ] ; then
   CC_PATH=`which $MOST_CC`
   echo >> most_compiler MOST_CC=$MOST_CC " # " $CC_PATH
   echo >> most_compiler "MOST_CC_OPTS=-O3"
# -I" $XINC_PATH
   echo "Found C compiler              at: $CC_PATH"
   if [ $MOST_CC = "suncc" ] ; then
      echo >> most_debug_options "MOST_CC_OPTS=-g -C -ftrap=common"
   else
      echo >> most_debug_options "MOST_CC_OPTS=-g"
   fi
else
   echo "****************************************************"
   echo "* Sorry - didn't find any C compiler               *"
   echo "* Please install one (e.g. gcc or g++              *"
   echo "* or update list of known compilers in this script *"
   echo "****************************************************"
   exit 2
fi

# check for C++ compiler
# put your favourite compiler in front if you have more than one

for MOST_PP in c++ g++ NO_PP
do
   `hash 2>/dev/null $MOST_PP`
   if [ $? = 0 ] ; then break ; fi
done

if [ $MOST_PP != "NO_PP" ] ; then
   PP_PATH=`which $MOST_PP`
   echo >> most_compiler MOST_PP=$MOST_PP " # " $PP_PATH
   echo "Found C++ compiler            at: $PP_PATH"
else
   echo "****************************************************"
   echo "* Sorry - didn't find any C++ compiler             *"
   echo "* Please install one (e.g. c++ or g++              *"
   echo "* or update list of known compilers in this script *"
   echo "****************************************************"
   exit 2
fi

#check for MPI-FORTRAN-90 compiler
for MPI_F90 in openmpif90 mpif90 mpf90 mpxlf90 NO_F90
do
   `hash 2>/dev/null $MPI_F90`
   if [ $? = 0 ] ; then break ; fi
done

if [ $MPI_F90 != "NO_F90" ] ; then
   F90_PATH=`which $MPI_F90`
   if [ $MPI_F90 = "openmpif90" ] ; then
      MPI_RUN="openmpiexec"
   else
      MPI_RUN="mpiexec"
   fi
   echo  > most_compiler_mpi "MPI_RUN=$MPI_RUN"
   echo >> most_compiler_mpi MOST_F90=$MPI_F90
   echo >> most_compiler_mpi "MOST_F90_OPTS=$MOST_F90_OPTS"
   echo >> most_compiler_mpi "MPIMOD=mpimod"
   echo >> most_compiler_mpi "UTILMOD=utilities"
   echo "Found MPI FORTRAN-90 compiler at: $F90_PATH"
else
   echo "********************************************"
   echo "* No Message Passing Interface (MPI) found *"
   echo "* Models run on single CPU                 *"
   echo "********************************************"
fi


#check for MPI-C compiler
if [ $MPI_F90 != "NO_F90" ] ; then
   for MPI_CC in openmpicc gxlc mpicc mpcc g++ gcc NO_CC
   do
      `hash 2>/dev/null $MPI_CC`
      if [ $? = 0 ] ; then break ; fi
   done
   if [ $MPI_CC != "NO_CC" ] ; then
      CC_PATH=`which $MPI_CC`
      echo >> most_compiler_mpi MOST_CC=$MPI_CC
      echo >> most_compiler_mpi "MOST_CC_OPTS=-O3"
      echo "Found MPI C compiler          at: $CC_PATH"
   else
      echo "********************************************"
      echo "* No MPI C compiler found                  *"
      echo "********************************************"
   fi
fi

#check for MPI-C++ compiler
if [ $MPI_F90 != "NO_F90" ] ; then
   for MPI_PP in openmpicxx mpicxx NO_PP
   do
      `hash 2>/dev/null $MPI_PP`
      if [ $? = 0 ] ; then break ; fi
   done
   if [ $MPI_PP != "NO_PP" ] ; then
      PP_PATH=`which $MPI_PP`
      echo >> most_compiler_mpi MOST_PP=$MPI_PP
      echo >> most_compiler_mpi "MOST_PP_OPTS=-O3 -DMPI"
      echo >> most_debug_options "MOST_PP_OPTS=-g -DMPI"
      echo "Found MPI C++ compiler        at: $PP_PATH"
   else
      echo "********************************************"
      echo "* No MPI C++ compiler found                *"
      echo "********************************************"
   fi
fi

#check for 86_64 operating system
MACHINE=`uname -m`

#check for GNU assembler
for MOST_AS in as NO_AS
do
   `hash 2>/dev/null $MOST_AS`
   if [ $? = 0 ] ; then break ; fi
done

if [ $MOST_AS != "NO_AS" ] ; then
   AS_PATH=`which $MOST_AS`
   echo >> most_compiler     MOST_AS=$MOST_AS
   echo >> most_compiler_mpi MOST_AS=$MOST_AS
   echo "Found GNU assembler           at: $AS_PATH"
fi

# activate fast assembler Legendre module for OS X and Linux 64 bit
if [ $MACHINE = "x86_64" -a $MOST_AS = "as" ] ; then
   # use fast vectorized assembler routines for Legendre transformation
   echo >> most_compiler_mpi LEGMOD=legini
   echo >> most_compiler_mpi LEGFAST=legfast32.o
   echo >> most_compiler     LEGMOD=legini
   echo >> most_compiler     LEGFAST=legfast32.o
   UNAMES=`uname -s`
   # generate correct symbol names for accessing gfortran module variables
   if [ $UNAMES = "Darwin" ] ; then
      echo "Fast Legendre Transformation    : ACTIVE (Darwin)"
   else
      echo "Fast Legendre Transformation    : ACTIVE (Linux)"
   fi
else
   # use FORTRAN routines for Legendre transformation
   echo >> most_compiler_mpi LEGMOD=legsym
   echo >> most_compiler_mpi LEGFAST=
   echo >> most_compiler     LEGMOD=legsym
   echo >> most_compiler     LEGFAST=
fi
# 
# echo >> most_compiler     "GUILIB=$GUILIB"
# echo >> most_compiler_mpi "GUILIB=$GUILIB"
echo >> most_compiler     "GUIMOD=guimod_stub"
echo >> most_compiler_mpi "GUIMOD=guimod_stub"
echo >> most_compiler     "PUMAX=pumax_stub"
echo >> most_compiler_mpi "PUMAX=pumax_stub"
# echo  > makefile "most.x: most.c"
# echo >> makefile "	$MOST_CC -o most.x most.c"
# #-I$XINC_PATH -lm $GUILIB"
# echo >> makefile "clean:"
# echo >> makefile "	rm -f *.o *.x F90* most_*"

#create directories

mkdir -p puma/bld
mkdir -p puma/bin
mkdir -p puma/run
mkdir -p plasim/bld
mkdir -p plasim/bin
mkdir -p plasim/run
mkdir -p plasim/dat/T21_exo
mkdir -p plasim/dat/T42_exo
#mkdir -p cat/bld
#mkdir -p cat/bin
#mkdir -p cat/run
#mkdir -p sam/bld
#mkdir -p sam/bin
#mkdir -p sam/run

#check FORTRAN I/O

export MOST_F90
export MOST_CC
HOSTNAME=`hostname`
export HOSTNAME
MOSTARCH=`uname -a`
export MOSTARCH
make -e -f makecheck
./f90check.x
./cc_check.x
echo >> most_info.txt "FORTRAN Compiler: $MOST_F90"
echo >> most_info.txt "C       Compiler: $MOST_CC"
cat most_info.txt
# make

#On Sunnyvale, the Python/NumPy stack was built with gcc/g++. Even if we would *like* to use 
#Intel compilers, on this machine we can't. This is likely going to be true on most systems,
#so we'll hard-code this, because there probably isn't a way to check this.

#If you're reading this and that's a problem, I'm so sorry. Blame compiler programmers who don't
#talk to each other. Or blame the NumPy devs for having zero cross-version compatibility in f2py.
oldcc=$CC
oldcpp=$CXX
export CC=gcc
export CXX=g++
#Compile pyfft libraries for Python 2
# f2py2 -c -m --f90exec=gfortran --f77exec=gfortran --f90flags="-O3" pyfft2 pyfft.f90 || f2py -c -m --f90exec=gfortran --f77exec=gfortran --f90flags="-O3" pyfft2 pyfft.f90
# f2py2 -c -m --f90exec=gfortran --f77exec=gfortran --f90flags="-O3" pyfft991v2 pyfft991.f90 || f2py -c -m --f90exec=gfortran --f77exec=gfortran --f90flags="-O3" pyfft991v2 pyfft991.f90

#Compile pyfft libraries for Python 3
# f2py$pyversion -c -m --f90exec=gfortran --f77exec=gfortran --f90flags="-O3" pyfft pyfft.f90 && mv pyfft.cpython*.so pyfft.so
# f2py$pyversion -c -m --f90exec=gfortran --f77exec=gfortran --f90flags="-O3" pyfft991 pyfft991.f90 && mv pyfft991.cpython*.so pyfft991.so
echo "python$pyversion -m numpy.f2py -c -m --f90exec=gfortran --f77exec=gfortran --f90flags="-O3" pyfft pyfft.f90 && mv pyfft.cpython*.so pyfft.so"
python$pyversion -m numpy.f2py -c --f90exec=gfortran --f77exec=gfortran --f90flags="-O3" -m pyfft pyfft.f90 && mv pyfft.cpython*.so pyfft.so
python$pyversion -m numpy.f2py -c --f90exec=gfortran --f77exec=gfortran --f90flags="-O3" -m pyfft991 pyfft991.f90 && mv pyfft991.cpython*.so pyfft991.so
    
    
export CC=$oldcc
export CXX=$oldcpp


# echo
# echo "configuration complete - run <most.x>"
# echo
#
#end of configure.sh
