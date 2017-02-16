#/bin/bash
#$ -cwd

echo "WORKDIR:   $SGE_O_WORKDIR"
echo "JOBSCRIPT: $JOB_SCRIPT"

if [ ${0:-3} == '.sh' ]; then
  ROOT=$(dirname $0)
elif [ ! -z "$SGE_O_WORKDIR" ]; then
  ROOT=$SGE_O_WORKDIR
else
  ROOT=$(pwd)
fi

if [ $(basename $ROOT) == 'tools' ]; then
  ROOT=$(dirname $ROOT)
fi

if [ $(basename $ROOT) != 'Siccuracy' ]; then
  echo "Nope, don\' know where I am."
  exit 1
fi

cd $TMPDIR

git clone $ROOT 

. /etc/profile.d/modules.sh
module load intel/2016

mkdir lib

R CMD INSTALL --no-inst --build --library=lib Siccuracy

cd Siccuracy/src

rm -f *.so *.o

for f in *.f95; do
  rm -f ${f%.f95}.f90
  cp $f ${f%.f95}.f90
  ifort -O3 -fpic -c ${f%.f95}.f90 -o ${f%.f95}.o
done

echo FC=ifort > Makevars
echo F77=ifort >> Makevars
export R_MAKEVARS_USER=Makevars
R CMD SHLIB -o Siccuracy.so *.f90
cd ../..

R CMD INSTALL --no-inst --build --library=lib Siccuracy

for f in Siccuracy_*.tar.gz; do
  mv $f $ROOT/../${f/linux-gnu/linux-intel}
done
