#/bin/bash
#$ -cwd
#$ -j y

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
module load R/3.3.2

mkdir lib

R --version
R CMD INSTALL --no-inst --build --library=lib Siccuracy

echo ====
echo Moving on to compiling with ifort
echo ====

cd Siccuracy/src

rm -f *.so *.o

for f in *.f95; do
  rm -f ${f%.f95}.f90
  cp $f ${f%.f95}.f90
  ifort -c -O3 -fpic  ${f%.f95}.f90 -o ${f%.f95}.o
done

for f in *.f95; do
  rm -f ${f%.f95}.f90
  cp $f ${f%.f95}.f90
  ifort -c -O3 -fpic -L/exports/applications/apps/SL7/R/3.3.2/lib64/R/lib -L/usr/local/lib64 ${f%.f95}.f90 -o ${f%.f95}.o -lR 
done

echo FC=ifort > Makevars
echo F77=ifort >> Makevars
export R_MAKEVARS_USER=$(pwd)/Makevar
R CMD SHLIB -o Siccuracy.so *.f90 
# >>  ifort -shared -L/exports/applications/apps/SL7/R/3.3.2/lib64/R/lib -L/usr/local/lib64 -o Siccuracy.so auxil.o converters.o hetero.o impacc.o mutate_chips.o plink.o -L/exports/applications/apps/SL7/R/3.3.2/lib64/R/lib -lR
#ifort -O3 -c -L/exports/applications/apps/SL7/R/3.3.2/lib64/R/lib -L/exports/applications/apps/SL7/intel/parallel_studio_xe_2016/compilers_and_libraries/linux/ipp/lib/intel64  -L/usr/local/lib64 -o Siccuracy.so *.o -lR
cd ../..

R CMD INSTALL --build --library=lib Siccuracy

module unload intel/2016

cat > test_package.R <<EOF
  .libPaths(c(file.path(getwd(), 'lib'), .libPaths()))
  library(testthat)
  library(Siccuracy, lib.loc=file.path(getwd(), 'lib'))
  test_path <- file.path(getwd(), 'Siccuracy', 'tests', 'testthat')
  testthat:::run_tests('Siccuracy', test_path, filter=NULL, reporter='summary')
EOF

Rscript --vanilla test_package.R

for f in Siccuracy_*.tar.gz; do
  mv $f $ROOT/../${f/linux-gnu/linux-intel}
done
