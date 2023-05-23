#!/bin/bash
#
#               EXOPLASIM COMPILATION SCRIPT
#
#  Usage:
#       ./compile.sh -p 8 -r T21 -v 10 -n 16 -O mavx 
#
#
#  Options:
#
#      -p:   Set floating point precision. Argument should be either 4 or 8. Default: 4
#
#      -r:   Set resolution. Can be T21/T42/T63, or alternatively the number of latitudes.
#            Default: T21
#
#      -v:   Set number of vertical levels. Default: 10
#
#      -n:   Number of processing cores (MPI threads) to use. Default: 4
#
#      -s:   Number of light sources in the sky. Default: 1
#
#      -t:   Number of years to use in most_plasim_run script. Default: 10
#
#      -O:   Specify additional compiler optimization flags, like -mavx (leave out prepended
#            hyphen).
#
#      -d:   Compile in debug mode (will produce line-number tracebacks on crash)
#
helptext=$(cat <<-END
               EXOPLASIM COMPILATION SCRIPT

  Usage:
       ./compile.sh -p 8 -r T21 -v 10 -n 16 -O mavx 


  Options:

      -p:   Set floating point precision. Argument should be either 4 or 8. Default: 4

      -r:   Set resolution. Can be T21/T42/T63, or alternatively the number of latitudes.
            Default: T21

      -v:   Set number of vertical levels. Default: 10
      
      -s:   Set number of light sources in the sky. Default: 1

      -n:   Number of processing cores (MPI threads) to use. Default: 4

      -t:   Number of years to use in most_plasim_run script. Default: 10
      
      -f:   Force FFT991mod.f90

      -O:   Specify additional compiler optimization flags, like -mavx (leave out prepended
            hyphen).

      -d:   Compile in debug mode (will produce line-number tracebacks on crash)

      -m:   Compile with Mars routines
      
      -h:   Output this text
END
)


prec=4
resolution="t21"
latitudes=32
longitudes=64
fftopt="fftmod"
levels=10
nlights=1
ncpus=4
debug=0
optimization=""
nopt=0
years=10
nmars=0

while getopts "p:r:v:n:O:t:dhm" opt; do
    case $opt in
        p)
            case $OPTARG in
                4)
                    prec=4 ;;
                8)
                    prec=8 ;;
                single)
                    prec=4 ;;
                double)
                    prec=8 ;;
                \?)
                    echo "INVALID PRECISION. Reverting to single precision."
                    ;;
            esac
            ;;
        r)
            case $OPTARG in #for >T63, recommend 30 levels
                T21)
                    resolution="t21"
                    latitudes=32
                    longitudes=64
                    ;;
                T42)
                    resolution="t42"
                    latitudes=64
                    longitudes=128
                    ;;
                T63)
                    resolution="t63"
                    latitudes=96
                    longitudes=192
                    fftopt="fft991mod"
                    ;;
                T85)
                    resolution="t85"
                    latitudes=128
                    longitudes=256
                    ;;
                T106)
                    resolution="t106"
                    latitudes=160
                    longitudes=320
                    fftopt="fft991mod"
                    ;;
                T127)
                    resolution="t127"
                    latitudes=192
                    longitudes=384
                    ;;
                T170)
                    resolution="t170"
                    latitudes=256
                    longitudes=512
                    ;;
                32)
                    resolution="t21"
                    latitudes=32
                    longitudes=64
                    ;;
                64)
                    resolution="t42"
                    latitudes=64
                    longitudes=128
                    ;;
                96)
                    resolution="t63"
                    latitudes=96
                    longitudes=192
                    fftopt="fft991mod"
                    ;;
                128)
                    resolution="t85"
                    latitudes=128
                    longitudes=256
                    ;;
                160)
                    resolution="t106"
                    latitudes=160
                    longitudes=320
                    fftopt="fft991mod"
                    ;;
                192)
                    resolution="t127"
                    latitudes=192
                    longitudes=384
                    ;;
                256)
                    resolution="t170"
                    latitudes=256
                    longitudes=512
                    ;;
                \?)
                    resolution="t21"
                    latitudes=32
                    longitudes=64
                    echo "INVALID RESOLUTION PASSED! Reverting to T21."
                    ;;
            esac ;;
        f)    
            fftopt="fft991mod"
            ;;
        v)
            levels=$OPTARG
            ;;
        n)
            ncpus=$OPTARG
            ;;
        d)
            debug=1
            ;;
        t)
            years=$OPTARG
            ;;
        O)
            optimization="-"$OPTARG
            nopt=1
            ;;
        h)
            echo "$helptext"
            exit 0
            ;;
        m)
            echo "Compiling for Mars..."
            export PLAMOD=p_mars
            nmars=1
            ;;
        s)
            nlights=$OPTARG
        \?)
            echo "UNRECOGNIZED ARGUMENT PASSED: "$OPTARG
            ;;
    esac
done

echo "PRODUCING: "$optimization" -r"$prec" -o most_plasim_"$resolution"_l"$levels"_s"$nlights"_p"$ncpus".x"
executable="most_plasim_"$resolution"_l"$levels"_s"$nlights"_p"$ncpus".x"

echo "Writing resmod.f90....."

cd plasim/bld/
rm -rf *


#       ! T85L30 on 16/32/64 processors
#       !parameter(NLAT_ATM = 128)
#       !parameter(NLEV_ATM = 30)
#       !!parameter(NPRO_ATM = 16)
#       !parameter(NPRO_ATM = 32)
#       !!parameter(NPRO_ATM = 64)
# 
# 
#       ! T127L30 on 16/32/64 processors
#       !parameter(NLAT_ATM = 192)
#       !parameter(NLEV_ATM = 30)
#       !!parameter(NPRO_ATM = 32)
#       !parameter(NPRO_ATM = 48) ! Does not work so well, 32 more efficient? \
# 
# 
# 
# 
#       ! T170L30 on 16/32/64 processors
#       !parameter(NLAT_ATM = 256)
#       !parameter(NLEV_ATM = 30)
#       !parameter(NPRO_ATM = 32)
#       !parameter(NPRO_ATM = 64) 


echo "      module resmod ! generated by compile.sh ">resmod.f90
echo "      parameter(NLAT_ATM = "$latitudes") ">>resmod.f90
echo "      parameter(NLEV_ATM = "$levels") ">>resmod.f90
echo "      parameter(NPRO_ATM = "$ncpus") ">>resmod.f90
echo "      parameter(NLIGHTS = "$nlights") ">>resmod.f90
echo "      end module resmod ">>resmod.f90
echo " ">>resmod.f90

rm plasim.x
rm ../bin/$executable
rm ../run/$executable
cp -p ../src/* .

if [ "$ncpus" -gt 1 ]
then
    [ ! -e MPI ] && rm -f *.o *.mod *.x
    touch MPI
    cp ../../most_compiler_mpi compilerargs
else
    [ ! -e MPI ] && rm -f *.o *.mod *.x MPI
    cp ../../most_compiler compilerargs
fi

dbgs=""
(($debug)) && dbgs='../../most_debug_options '

cp ../../most_precision_optionsx precisionargsx
if [ "$prec" -gt 4 ]
then
    sed -i.bak '1s/$/'$prec'/' precisionargsx
    rm -rf precisionargsx.bak
else
    echo "">precisionargsx
fi

(($nopt)) && sed -i.bak '3s/$/ '$optimization'/' compilerargs && rm -rf compilerargs.bak

cat compilerargs $dbgs../bld/precisionargsx make_plasim > makefile

echo "Writing makefile..."
echo ""

export OCEANCOUP=cpl_stub
export FFTMOD=$fftopt
#cat makefile

make -e
./most_snow_build$prec
./most_ice_build$prec
cp plasim.x ../bin/$executable
cp ../bin/$executable ../run/
cd ../../

(($nmars)) && cp plasim/dat/T"${resolution[@]:1}"_mars/* plasim/run/ || cp plasim/dat/T"${resolution[@]:1}"/* plasim/run/

namelist=example.nl
(($nmars)) && namelist=mars.nl
(($nmars)) && cp postprocessor/mars.nl plasim/run/

snamelist=snapshot.nl
(($nmars)) && snamelist=mars_snapshot.nl
(($nmars)) && cp postprocessor/mars_snapshot.nl plasim/run/

echo "#!/bin/bash ">plasim/run/most_plasim_run
echo "# run-script generated by compile.sh                       ">>plasim/run/most_plasim_run
echo "EXP=MOST    # Name your experiment here                    ">>plasim/run/most_plasim_run
echo "[ \$# ==1 ] && cd \$1                                      ">>plasim/run/most_plasim_run
echo "rm -f plasim_restart                                       ">>plasim/run/most_plasim_run
echo "rm -f Abort_Message                                        ">>plasim/run/most_plasim_run
echo "YEAR=0                                                     ">>plasim/run/most_plasim_run
echo "YEARS="$years"                                             ">>plasim/run/most_plasim_run
echo "while [ \$YEAR -lt \$YEARS ]                               ">>plasim/run/most_plasim_run
echo "do                                                         ">>plasim/run/most_plasim_run
echo "   YEAR=\`expr \$YEAR + 1\`                                ">>plasim/run/most_plasim_run
echo "   DATANAME=\`printf '%s.%03d' \$EXP \$YEAR\`              ">>plasim/run/most_plasim_run
echo "   SNAPNAME=\`printf '%s_SNAP.%03d' \$EXP \$YEAR\`         ">>plasim/run/most_plasim_run
echo "   DIAGNAME=\`printf '%s_DIAG.%03d' \$EXP \$YEAR\`         ">>plasim/run/most_plasim_run
echo "   RESTNAME=\`printf '%s_REST.%03d' \$EXP \$YEAR\`         ">>plasim/run/most_plasim_run
echo "   SNOWNAME=\`printf '%s_SNOW.%03d' \$EXP \$YEAR\`         ">>plasim/run/most_plasim_run
if [ "$ncpus" -gt 1 ]
then
   MPI_RUN=$(head -n 1 most_compiler_mpi | tr "=" "\n" | tail -1)
   echo "   $MPI_RUN -np $ncpus $executable                      ">>plasim/run/most_plasim_run
else
   echo "   ./$executable                                        ">>plasim/run/most_plasim_run
fi
echo "   [ -e Abort_Message ] && exit 1                          ">>plasim/run/most_plasim_run
echo "   [ -e plasim_output ] && mv plasim_output \$DATANAME     ">>plasim/run/most_plasim_run
echo "   [ -e plasim_snapshot ] && mv plasim_snapshot \$SNAPNAME ">>plasim/run/most_plasim_run
echo "   [[ -e burn7.x && -e "$namelist" && -e \$DATANAME ]] && ./burn7.x -n <"$namelist">burnout \$DATANAME \$DATANAME.nc ">>plasim/run/most_plasim_run
echo "   [[ -e burn7.x && -e "$snamelist" && -e \$SNAPNAME ]] && ./burn7.x -n <"$snamelist">snapout \$SNAPNAME \$SNAPNAME.nc ">>plasim/run/most_plasim_run
echo "   [ -e plasim_diag ] && mv plasim_diag \$DIAGNAME         ">>plasim/run/most_plasim_run
echo "   [ -e plasim_status ] && cp plasim_status plasim_restart ">>plasim/run/most_plasim_run
echo "   [ -e plasim_status ] && mv plasim_status \$RESTNAME     ">>plasim/run/most_plasim_run
echo "   [ -e restart_snow ] && mv restart_snow \$SNOWNAME       ">>plasim/run/most_plasim_run
echo "done                                                       ">>plasim/run/most_plasim_run

