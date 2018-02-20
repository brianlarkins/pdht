#!/bin/bash

#exe="diff"
exe=scaling

#
# make_file nodes procs tasks-per-node 
#
function makefile () {
  n=$1      # total number of srun/mpirun tasks (-n)
  wpn=$2    # number of PDHT workers per node
  mpi=$3    # mpi implementation string

  cputask=2                     # mpi and pdht run 2 cores per process
  tpn=`expr $tpn / $cputask`    # workers per node

  if [ "$mpi" = "pdht" ]; then
    ext=""
    tpn=$wpn                     # 1 tasks per worker
  else
    ext="MPI"
    tpn=`expr $wpn '*' 2`        # 2 tasks per 1 worker (1 init + 1 target)
  fi
 
  tasks=`expr $n '*' $tpn` 
  workers=`expr $n '*' $wpn`

  sdir="scripts/"
  sfname="$sdir/$exe.$mpi.$workers$runtype.sh"

  runner="srun"
  uenv="env PTL_IGNORE_UMMUNOTIFY=1"
  #uenv="env PTL_DISABLE_MEM_REG_CACHE=1"

  mopts=""

  if [ "$mpi" = "ompi" ]; then
    runner="mpirun"
    mopts="--mca pml ob1 --mca btl portals4 --mca coll ^portals4 --bind-to none"
  fi 

  echo "nodes: $n total tasks: $tasks tpn: $tpn cputasks: $cputask workers: $workers mpi: $mpi file: $sfname"

  echo '#!/bin/bash' > $sfname
  echo '#SBATCH -p haswell' >> $sfname
  echo '#SBATCH -t 30:00' >> $sfname
  echo "#SBATCH -n $tasks" >> $sfname
  echo "#SBATCH -N $n" >> $sfname
  echo "#SBATCH --ntasks-per-node $tpn" >> $sfname
  echo "#SBATCH -c $cputask" >> $sfname
  echo "#SBATCH --exclusive" >> $sfname
  cat $mpi.skel >> $sfname
  echo ' ' >> $sfname
  echo "$runner $mopts $uenv $path/$exe$ext $opts > $HOME/pdht/runs/$mpi/$exe.$workers$runtype 2>&1" >> $sfname
}


if [ "$exe" = "diff" ]; then
  path=$HOME/pdht/bench/diff
  opts=""
  runtype=""
elif [ $exe = "scaling" ]; then
  path=$HOME/pdht/test
  runtype=".10k"
  #opts="-quUlL -v 10000 -s 1024 -i 5"
  opts="-quUlL -v 1000 -s 128 -i 1"
  #runtype=50k
  #opts="-quUl -v 50000 -s 1024 -i 5"
fi

maxtasks_per_node=6

#
# main
#
for m in mpich mvapich ompi pdht
do
  mkdir -p $m
  for nodes in 1 2 3 4 6 8 12 16 24 32 
  do 
    if [ $nodes -eq 1 ]; then
      # special case 1 2 4 8 processor runs
      for p in 1 2 4 8 
      do
        makefile $nodes $p $m
      done
    else 
      makefile $nodes $maxtasks_per_node  $m
    fi
  done
done

