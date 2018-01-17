#!/bin/bash

exe="diff"
#exe=scaling

#
# make_file nodes procs tasks-per-node 
#
function makefile () {
  n=$1
  tpn=$2
  mpi=$3

  cputask=1  # mpi = 1 core per process

  if [ "$mpi" = "pdht" ]; then
    ext=""
    cputask=2                     # pdht = 2 cores per process
    tpn=`expr $tpn / $cputask`    # 32 / 2 = 16
    tasks=`expr $n '*' $tpn`      # nodes * 16
    workers=$tasks
  else
    ext="MPI"
    tpn=`expr $tpn / $cputask`    # 32 / 1 = 32
    tasks=`expr $n '*' $tpn`      # nodes * 32
    workers=`expr $tasks / 2`     # tasks / 2 (1 init + 1 target per worker)
  fi


  sdir="scripts/"
  sfname="$sdir/$exe.$mpi.$workers$runtype.sh"

  runner="srun"
  uenv="env PTL_IGNORE_UMMUNOTIFY=1"

  if [ "$mpi" = "ompi" ]; then
    runner="mpiexec"
  fi 


  echo '#!/bin/bash' > $sfname
  echo '#SBATCH -p haswell' >> $sfname
  echo '#SBATCH -t 30:00' >> $sfname
  echo "#SBATCH -n $tasks" >> $sfname
  echo "#SBATCH -N $n" >> $sfname
  echo "#SBATCH --ntasks-per-node $tpn" >> $sfname
  echo "#SBATCH -c $cputask" >> $sfname
  cat $mpi.skel >> $sfname
  echo ' ' >> $sfname
  echo "$runner $uenv $path/$exe$ext > $HOME/pdht/runs/$mpi/$exe.$workers$runtype 2>&1" >> $sfname
}



if [ "$exe" = "diff" ]; then
  path=$HOME/pdht/bench/diff
  opts=""
  runtype=""
elif [ $exe = "scaling" ]; then
  path=$HOME/pdht/test
  runtype=".10k"
  opts="-quUl -v 10000 -s 1024 -i 5"
  #runtype=50k
  #opts="-quUl -v 50000 -s 1024 -i 5"
fi

#
# main
#
for m in mpich mvapich ompi pdht
do
  for nodes in 1 2 3 4 6 8 12 16 24 32 
  do 
    if [ $nodes -eq 1 ]; then
      # special case 1 2 4 8 processor runs
      for p in 1 2 4 8 16
      do
        makefile $nodes $p $m
      done
    else 
      makefile $nodes 32 $m
    fi
  done
done

