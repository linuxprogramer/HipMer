#!/bin/bash

USAGE="PROPER USAGE:

$0 install-path meraculous.config

"

PHYS_MEM_MB=${PHYS_MEM_MB:=$(awk '/MemTotal:/ { t=$2 ; print int(0.80 * t / 1024)}' /proc/meminfo || echo 3000)}
export GASNET_PHYSMEM_MAX=${GASNET_PHYSMEM_MAX:=${PHYS_MEM_MB}MB}
export GASNET_PHYSMEM_NOPROBE=${GASNET_PHYSMEM_NOPROBE:=1}
UPC_PTHREADS=${UPC_PTHREADS:=}
UPC_PTHREADS_OPT=
if [ -n "$UPC_PTHREADS" ] && [ $UPC_PTHREADS -ne 0 ]
then
  UPC_PTHREADS_OPT="-pthreads=${UPC_PTHREADS}"
fi

if [ -z "$1" -o -z "$2" ]
then
  echo "$USAGE" 1>&2
  echo "Please specify the hipmer install path and config file" 1>&2
  exit 1
fi

INSTALL_PATH=$(readlink -f $1 2>/dev/null)
if [ ! -d "${INSTALL_PATH}" ]
then
  echo "$USAGE" 1>&2
  echo "Invalid hipmer install path: ${INSTALL_PATH}" 1>&2
  exit 1
fi

CONFIG=$(readlink -f $2)
if [ ! -f "${CONFIG}" ]
then
  echo "$USAGE
Please specify a config file! '${CONFIG}' is invalid!
" 1>&2
  exit 1
fi
export PATH=${INSTALL_PATH}/bin:${PATH}

HIPMER_VERSION=$(cat ${INSTALL_PATH}/HIPMER_VERSION || echo 'UNKNOWN VERSION')
echo "# Starting '$0 $@' (${HIPMER_VERSION}) at $(date) in $(pwd) on $(hostname) (${SLURM_JOB_ID}${PBS_JOBID}${JOB_ID} pid:$$)"

DRYRUN=${DRYRUN:=0}

LAST_LOG=
bail()
{
  if [ "$DRYRUN" != "0" ]
  then
     echo "# Dryrun: Skipping check $@"
     return
  fi
  echo "

ERROR: at $(date) in $(pwd) on $(hostname):
$@
    ${LAST_LOG}

" 1>&2
  exit 1
}

get_config()
{
  field=$1
  vals=($(grep "^${field}" ${CONFIG} || /bin/true))
  echo ${vals[1]}
}

get_rankpath()
{
  rank="0000000${1}"
  echo "${PER_THREAD_DIR}/${rank: -7}"
}
  

libfiles=
libname=
libinsavg=
libinsstddev=
libreadlen=
libinnie=
librevcomp=
libcontigging=
libonosetid=
libgapclosing=
lib5pwiggle=
lib3pwiggle=
libfilesperpair=
libsplinting=
get_libseq_vars()
{
  fileline=$1
  read -a vals < ${fileline}
  libfiles=($(echo $(echo "${vals[1]}" | tr ',' ' ')))
  libname="${vals[2]}"
  libinsavg=${vals[3]}
  libinsstddev=${vals[4]}
  libreadlen=${vals[5]}
  libinnie=${vals[6]}
  librevcomp=${vals[7]}
  libcontigging=${vals[8]}
  libonosetid=${vals[9]}
  libgapclosing=${vals[10]}
  lib5pwiggle=${vals[11]}
  lib3pwiggle=${vals[12]}
  libfilesperpair=${vals[13]}
  libsplinting=${vals[14]}
  if [ "${vals[1]%,*}" != "${vals[1]}" ]
  then
    [ ${libfilesperpair} -eq 2 ] || bail "Invalid config file: ${libname} has a comma in the filename field but ${libfilesperpair} in FilesPerPair"
    # comma separates paired files, shuffle ordering of files from: a1 b1 c1 a2 b2 c2 to: a1 a2 b1 b2 c1 c2
    declare -a newlibfiles
    numFiles=${#libfiles[*]}
    halfFiles=$((numFiles/2))
    for fnum in $(seq 1 ${halfFiles})
    do
       newlibfiles=(${newlibfiles[@]} ${libfiles[$((fnum-1))]} ${libfiles[$((fnum+halfFiles-1))]})
    done
    libfiles=(${newlibfiles[@]})
  fi
}

# this is the base dir in which all the temporary pipeline files are written
if [ "$USE_SHM" ]
then 
    BASE_DIR="/dev/shm"
    echo "# Using SHM for caching files"
else
    BASE_DIR="."
    echo "# Using base directory '$BASE_DIR'"
fi

MIN_DEPTH_CUTOFF=${MIN_DEPTH_CUTOFF:=$(get_config min_depth_cutoff)}
if [ -z "${MIN_DEPTH_CUTOFF}" ] || [ ${MIN_DEPTH_CUTOFF} -lt 0 ]
then
  MIN_DEPTH_CUTOFF=0 # auto detect
fi
KMER_LENGTH=${KMER_LENGTH:=$(get_config mer_size)}
MIN_CONTIG_LENGTH=${MIN_CONTIG_LENGTH:=${KMER_LENGTH}}
NUM_ALIGNMENTS=${NUM_ALIGNMENTS:=192000}

EXCLUDE_REPEATS=$(get_config gap_close_rpt_depth_ratio)
if [ "$EXCLUDE_REPEATS" == "" ]; then
  echo "# No value set for EXCLUDE_REPEATS, using 2.0"    
  EXCLUDE_REPEATS="2.0"
else
  echo "# using ${EXCLUDE_REPEATS} for EXCLUDE_REPEATS (a.k.a. gap_close_rpt_depth_ratio or merauderRepeatCopyCount)"
fi

is_diploid=$(get_config is_diploid)
if [ "$is_diploid" == "" ]; then
  is_diploid=0
fi

if [ $is_diploid -eq 0 ]; then
  echo "# Diploid NOT set"
  POLYMODE=""
else
  echo "# Diploid SET"
  POLYMODE="-P"
fi

#ILLUMINA_VERSION=${ILLUMINA_VERSION:=18}
HIPMER_SW_CACHE=${HIPMER_SW_CACHE:=6144}
RUNDIR=${RUNDIR:=$(pwd)/$(hostname)-$$-$(date '+%Y%m%d-%H%M%S')}


trap ctrl_c INT

function ctrl_c() {
    trap "" INT
    [ -n "KEEP_RUNDIR" ] || echo "deleting $RUNDIR... hit ctrl-c again to abort purge"
    sleep 5
    echo "aborting "
    pkill -2 -P $$
    [ -n "KEEP_RUNDIR" ] || rm -rf $RUNDIR
}


guess_qual_offset()
{
  fastq=$1
  if [ ! -f "$fastq" ]
  then
    echo "Error: $fastq does not exist"
    exit 1
  fi

  check0=$(head -n50000 $fastq | grep -A1 "^\+" | grep -v "^\+" | grep "[\"#$%&'()*,./0123456789:;<=>?]" |wc -l)
  if [ $check0 -eq 0 ]
  then
    echo 64
  else
    echo 33
  fi
}

if [ -z "${CORES_PER_NODE}" ] || [ $CORES_PER_NODE -eq 0 ]
then
  # auto-detect the cores per node
  if [ -x $(which lscpu) ]
  then
    export CORES_PER_NODE=$(($(lscpu -p | tail -1 | awk -F, '{print $2}')+1))
    echo "# Autodetected CORES_PER_NODE=${CORES_PER_NODE} from lscpu"
  else 
    SOCKETS=$((grep '^physical id' /proc/cpuinfo | sort | uniq | wc -l))
    CORES_PER_SOCKET=$((grep '^core id' /proc/cpuinfo | sort | uniq | wc -l))
    export CORES_PER_NODE=$((SOCKETS*CORES_PER_SOCKET))
    echo "# Autodetected CORES_PER_NODE=${CORES_PER_NODE} from /proc/cpuinfo"
  fi
fi

HYPERTHREADS=${HYPERTHREADS:=1}
CORES_PER_NODE=$((CORES_PER_NODE*HYPERTHREADS))
UFX_HLL=${UFX_HLL:=0}
UFX_OPT=${UFX_OPTS:=}
if [ "${UFX_HLL}" != "0" ]
then
  UFX_OPT="${UFX_OPT} -H"
fi
BUBBLE_MIN_DEPTH_CUTOFF=${BUBBLE_MIN_DEPTH_CUTOFF:=$(get_config bubble_min_depth_cutoff)}

POLITE_SYNC=
if [ ${HYPERTHREADS} -gt 1 ]
then
  POLITE_SYNC="-polite-sync"
fi

THREADS=${THREADS:=${CORES_PER_NODE}}
if [ ${THREADS} -lt ${CORES_PER_NODE} ]
then
  # support 1 thread or less than the cores in one node on a single node
  echo "Overriding CORES_PER_NODE from ${CORES_PER_NODE} to ${THREADS}"
  export CORES_PER_NODE=${THREADS}
elif [ $((THREADS % CORES_PER_NODE)) -ne 0 ]
then
  echo "Overriding THREADS from ${THREADS} to be a mulitple of ${CORES_PER_NODE}"
  THREADS=$((THREADS - (THREADS % CORES_PER_NODE)))
fi
echo "# Executing on a total of ${THREADS} threads"

NODES=$((${THREADS}/${CORES_PER_NODE}))
ONO_ACTIVE_THREADS=${ONO_ACTIVE_THREADS:=${CORES_PER_NODE}}

echo "# Using $NODES nodes, with $CORES_PER_NODE cores per node and ${HYPERTHREADS} hyperthreads per core"
echo "# Using $ONO_ACTIVE_THREADS active threads for oNo"


if [ -z "${CAN_SPLIT_JOB}" ]
then
  if [ -n "${SLURM_JOB_ID}" ]
  then
    CAN_SPLIT_JOB=1
  else
    CAN_SPLIT_JOB=0
  fi
fi
if [ -n "${CAN_SPLIT_JOB}" ]; then
    echo "# Can split job: ${CAN_SPLIT_JOB}"
fi


MAX_UPC_SHARED_HEAP_MB=$((PHYS_MEM_MB/CORES_PER_NODE))
if [ -z "${UPC_SHARED_HEAP_MB}" ]
then
  UPC_SHARED_HEAP_MB=${MAX_UPC_SHARED_HEAP_MB}
fi
MAX_ALLOWED_SHARED_HEAP_MB=${MAX_ALLOWED_SHARED_HEAP_MB:=15500}
# maximum allowed shared heap in bupc...
if [ ${UPC_SHARED_HEAP_MB} -gt ${MAX_ALLOWED_SHARED_HEAP_MB} ]
then
  UPC_SHARED_HEAP_MB=${MAX_ALLOWED_SHARED_HEAP_MB}
fi
LFS_STRIPE=${LFS_STRIPE:=72}

echo "# Using -shared-heap=${UPC_SHARED_HEAP_MB}M"

MPIRUN=${MPIRUN:="mpirun -n"}
MPIRUN="${MPIRUN} $((THREADS/HYPERTHREADS))"

UPCRUN=${UPCRUN:="upcrun -n"}
# parse into _UPCRUN="upcrun" and the arguments
_UPCRUN=${UPCRUN%% *}
if [ -n "${DEBUG}" ]
then
  ulimit -c unlimited || /bin/true
  DEBUG_UPC="-backtrace -abort -verbose"
fi
UPCRUN_NO_PTHREADS="${_UPCRUN} ${DEBUG_UPC} ${POLITE_SYNC} -shared-heap=${UPC_SHARED_HEAP_MB}M ${UPCRUN#${_UPCRUN}}"
UPCRUN="${_UPCRUN} ${UPC_PTHREADS_OPT} ${DEBUG_UPC} ${POLITE_SYNC} -shared-heap=${UPC_SHARED_HEAP_MB}M ${UPCRUN#${_UPCRUN}}"

ENVS="
environmental variables with some effects:

Physical memory properties:
  PHYS_MEM_MB=${PHYS_MEM_MB}

The number of UPC threads in the job (MPI will not use hyperthreads):
  THREADS=${THREADS} HYPERTHREADS=${HYPERTHREADS}

The rundirectory to place all the files (will be created if necessary):
  RUNDIR=${RUNDIR}
  DATADIR=${DATADIR}
  LFS_STRIPE=${LFS_STRIPE}

Gasnet properties (by default 80% of physmem):
  GASNET_PHYSMEM_MAX=${GASNET_PHYSMEM_MAX}
  GASNET_PHYSMEM_NOPROBE=${GASNET_PHYSMEM_NOPROBE}

UPC properties:
  UPC_SHARED_HEAP_MB=${UPC_SHARED_HEAP_MB} (Do not set to use 80% of the node memory)
  UPC_PTHREADS=${UPC_PTHREADS}
  UPC_PTHREADS_OPT=${UPC_PTHREADS_OPT}

HipMer options (will override config file defaults):
  UFX_HLL=${UFX_HLL}
  MIN_DEPTH_CUTOFF=${MIN_DEPTH_CUTOFF} # use 0 for auto-detect after UFX generation
  BUBBLE_MIN_DEPTH_CUTOFF=${BUBBLE_MIN_DEPTH_CUTOFF}
  KMER_LENGTH=${KMER_LENGTH}
  MIN_CONTIG_LENGTH=${MIN_CONTIG_LENGTH}
  NUM_ALIGNMENTS=${NUM_ALIGNMENTS}
  ILLUMINA_VERSION=${ILLUMINA_VERSION}
  HIPMER_SW_CACHE=${HIPMER_SW_CACHE}
  NO_CONCAT=${NO_CONCAT}
  CAN_SPLIT_JOB=${CAN_SPLIT_JOB}
  ONO_ACTIVE_THREADS=${ONO_ACTIVE_THREADS}

MPI/UPC environment:
  MPIRUN=${MPIRUN}
  UPCRUN=${UPCRUN}
"
USAGE="${USAGE}
environmental variables with some effects:
${ENVS}

"

set -e

MEMTIME=memtime
which ${MEMTIME} > /dev/null || MEMTIME=
which ${UPCRUN%% *} > /dev/null || bail "${USAGE} upcrun (${UPCRUN}) must be installed and in your PATH!"
which ${MPIRUN%% *} > /dev/null || bail "${USAGE} mpirun (${MPIRUN}) must be installed and in your PATH!"
which meraculous-${KMER_LENGTH} > /dev/null || bail "${USAGE}Could not find meraculous-${KMER_LENGTH}!  Have you compiled with the proper option: '-DHIPMER_KMER_LENGTHS=${KMER_LENGTH}'?"
which ufx-${KMER_LENGTH} > /dev/null || bail "could not find ufx-${KMER_LENGTH}! Have you compiled with the proper option '-D$KMER_LENGTH'?"

echo "Using:
   upcrun: ${UPCRUN}
   mpirun: ${MPIRUN}
   $(which ufx-${KMER_LENGTH})
   $(which meraculous-${KMER_LENGTH})
${ENVS}
"


[ -n "${RUNDIR}" ] && cd ${RUNDIR} || bail "Could not find ${RUNDIR} from $(pwd) on $(uname -n)"
RUNDIR=$(pwd -P)
PER_THREAD_DIR=per_thread
rankpath=$(get_rankpath 0)

mkdir -p ${PER_THREAD_DIR}
# if on luster, set stripe to be 1 for per_thread files
[ ! -x "$(which lfs)" ] || lfs setstripe -c 1 ${PER_THREAD_DIR} 2>/dev/null || /bin/true

mkdir -p ${rankpath}

# set all other new files defaulted to be on a large number the OSTs
[ ! -x "$(which lfs)" ] || lfs setstripe -c ${LFS_STRIPE} ${RUNDIR}/. 2>/dev/null || /bin/true

CONFIG_PATH=${CONFIG%/*}
# split each lib_seq entry in the config file into its own separate file
if [ ! -f config.LIBSEQaa ]
then
  grep '^lib_seq' ${CONFIG} | split -l 1 - ${PER_THREAD_DIR}/config.LIBSEQ
  [ -f ${PER_THREAD_DIR}/config.LIBSEQaa ] || bail "Could not generate config.LIBSEQaa!"
  for i in ${PER_THREAD_DIR}/config.LIBSEQ* ; do ln $i ; done
fi

# now each config.LIBSEQ* file will have 1 entry of the form:
# lib_seq [ wildcard ][ prefix ][ insAvg ][ insSdev ][ avgReadLen ][ hasInnieArtifact ][ isRevComped ][ useForContigging ][ onoSetId ][ useForGapClosing ][ 5pWiggleRoom ][3pWiggleRoom] [FilesPerPair] [ useForSplinting ]
#is_diploid=$(get_config is_diploid)
if [ $is_diploid -eq 0 ] && [ $MIN_CONTIG_LENGTH -lt $((KMER_LENGTH*2)) ]
then
  echo "# Setting MIN_CONTIG_LENGTH to 2*${KMER_LENGTH} because this is not a diploid assembly"
  MIN_CONTIG_LENGTH=$((KMER_LENGTH*2))
fi

GAPLIBS=
ALLLIBS=
MAX_INSERT_SIZE=0
MAX_INSERT_SIGMA=0
MAX_ONO_SET_ID=0
IN=ALL_INPUTS.fofn
rm -f ${IN} ${PER_THREAD_DIR}/${IN}
# Check for existance of libname.fofn for each library
for x in config.LIBSEQ*
do
  get_libseq_vars $x
  if [ -n "$ALLLIBS" ]
  then
    ALLLIBS="${ALLLIBS},${libname}"
  else
    ALLLIBS="${libname}"
  fi
  if [ $libgapclosing -eq 1 ]
  then
    if [ -n "$GAPLIBS" ]
    then
      GAPLIBS="${GAPLIBS},${libname}"
    else
      GAPLIBS="${libname}"
    fi
  fi
  if [ ${MAX_ONO_SET_ID} -lt ${libonosetid} ]
  then
    MAX_ONO_SET_ID=${libonosetid}
  fi
  if [ -f "${libname}.fofn" ]
  then
    echo "# For library ${libfiles[@]} ${libname} using ${libname}.fofn:"
    echo "$(cat ${libname}.fofn)" 
  else
    for file in ${libfiles[@]}
    do
      if [ -f ${file} ]
      then
        echo $(readlink -f ${file})
      elif [ -d "${DATADIR}" ] && [ -f ${DATADIR}/${file} ]
      then
        echo $(readlink -f ${DATADIR}/${file})
      else
        for file1 in ${CONFIG_PATH}/${file}
        do
          if [ -f ${file1} ]
          then
            echo $(readlink -f ${file1})
          else
            bail "Could not find $file or $file1 specified as ${libfile} for Library ${libname}"
         fi
       done
      fi
    done > ${PER_THREAD_DIR}/${libname}.fofn
    ln -f ${PER_THREAD_DIR}/${libname}.fofn
  fi
  if [ ${libcontigging} -eq 1 ]
  then
    cat ${libname}.fofn >> ${PER_THREAD_DIR}/${IN} 
  fi
done 
ln -f ${PER_THREAD_DIR}/${IN}

INUFX=${IN}.ufx.bin
OUT=${OUT:=contigs}

TIMINGS=${TIMINGS:=${PER_THREAD_DIR}/timings.log}

# needed in the rerun_stage.sh script 
echo "# Starting '$0 $@' at $(date) in $(pwd) on $(hostname)
# Install path: $INSTALL_PATH
# RUNDIR: ${RUNDIR}" >> $TIMINGS
ln -f ${TIMINGS} 2>/dev/null || /bin/true

log_and_run()
{
  if [ -z "${RETRY}" ]
  then
    RETRY=0
  else
    RETRY=$((RETRY-1))
  fi
  log=${PER_THREAD_DIR}/$1
  loglink=$1
  LAST_LOG=$(pwd)/${log}.$(date '+%Y%m%d_%H%M%S')
  shift
  if [ -f ${log} ] && [  -f ${loglink} ]
  then
    echo "# Skipping, already executed: '$@'"
  else
    echo "# Starting '$@' at $(date) log in ${LAST_LOG}"
    echo "# Starting $@ at $(date)" >> ${TIMINGS}
    secs=$SECONDS
    if [ "$DRYRUN" != "0" ]
    then
      echo "# DRY run would execute:
   $@
"
      return
    else
      echo "$@" >> ${TIMINGS}
      echo "Starting '$@' at $(date)" > ${LAST_LOG}
      $@ >> ${LAST_LOG} 2>&1 \
        || ( [ ${RETRY} -gt 0 ] && RETRY=${RETRY} log_and_run ${loglink} $@ ) \
        || bail "Could not execute $@:
    (exit: $? after $((SECONDS-$secs)) s) See the log:
    ${LAST_LOG}"
    fi
    runsecs=$((SECONDS-$secs))
    echo "# Finished (${runsecs} s) '$@' at $(date)"
    echo "# ${log} ${runsecs} walltime-secs $((runsecs*THREADS)) core-seconds at $(date)" >> ${TIMINGS}
    echo "# ${log} ${runsecs} walltime-secs $((runsecs*THREADS)) core-seconds at $(date)" >> ${LAST_LOG}
    ln -f ${LAST_LOG} $log
    ln -f ${log} ${loglink}
    LAST_LOG=

    sleep 1 # insert a small delay between invocations of srun / upcrun to allow shared memory to clean up
  fi
}

if [ -z "$ILLUMINA_VERSION" ]; then
  a_lib_fname=$(head -n1 ${IN})
  [ ! -f "$a_lib_fame" ] || bail "Could not find $a_lib_fname to guess ILLUMINA format"
  # guess the illumina version
  qual_offset=$(guess_qual_offset $a_lib_fname)
  echo "# Guessed phred offset=${qual_offset}. Set ILLUMINA_VERSION to override."
else
  if [ $ILLUMINA_VERSION -eq 13 ] || [ $ILLUMINA_VERSION -eq 15 ]; then
    qual_offset=64
  elif [ $ILLUMINA_VERSION -eq 18 ]; then
    qual_offset=33
  fi
  echo "# Phred offset $qual_offset"
fi

if [ -z "$qual_offset" ]; then
  echo "Could not guess illumina version... assuming phred offset is 33"
  qual_offset=33
fi

if [ "$USE_SHM" ]
then
    SHM_FLAG=-s
    log_and_run loadfq.log $MEMTIME $UPCRUN ${THREADS} loadfq -f ${ALLLIBS} 
fi

UFX_ERR_COUNT=${MIN_DEPTH_CUTOFF}
if [ ${UFX_ERR_COUNT} -le 1 ]
then
  UFX_ERR_COUNT=2
fi
log_and_run ufx.log $MEMTIME $MPIRUN ufx-${KMER_LENGTH} -Q $qual_offset -f ${IN} -e ${UFX_ERR_COUNT} ${UFX_OPT} $SHM_FLAG
# This should create ALL_INPUTS.fofn.ufx.bin
if [ "$USE_SHM" ]
then
    EXTUFX=${BASE_DIR}/${INUFX}0
else
    EXTUFX=${INUFX}
    [ -f ${EXTUFX} ] || bail "No ${EXTUFX} file!"
fi
HIST_FILE=${PER_THREAD_DIR}/histogram_k${KMER_LENGTH}.txt
[ -f ${HIST_FILE} ] || bail "ufx-${KMER_LENGTH} did not produce ${HIST_FILE}"
ln -f ${HIST_FILE}

if [ -z "${MIN_DEPTH_CUTOFF}" ] || [ ${MIN_DEPTH_CUTOFF} -eq 0 ]
then
  echo "Calculating min_depth_cutoff..."
  if [ ! -f findDMin2.txt ]
  then
    MIN_DEPTH_CUTOFF=$(findDMin2.pl ${HIST_FILE})
    [ -z "${MIN_DEPTH_CUTOFF}" ] || (echo ${MIN_DEPTH_CUTOFF} > ${PER_THREAD_DIR}/findDmin2.txt && ln -f ${PER_THREAD_DIR}/findDmin2.txt)
  else
    MIN_DEPTH_CUTOFF=$(cat findDmin2.txt)
  fi
  if [ -z "${MIN_DEPTH_CUTOFF}" ] || [ ${MIN_DEPTH_CUTOFF} -lt 2 ]
  then
      bail "Could not calculate a proper MIN_DEPTH_CUTOFF (${MIN_DEPTH_CUTOFF}).
	Please evaluate the ${HIST_FILE} and provide either:
		MIN_DEPTH_CUTOFF= (ENVironmental variable) or
		min_depth_cutoff paramter in the ${HIPMER_CONFIG}
"
  fi
  echo "Using auto-calculated value for MIN_DEPTH_CUTOFF=${MIN_DEPTH_CUTOFF}"
fi

if [ -n "$DEBUG" ] && [ ! ${USE_SHM} ] && [ -z "{NO_CONCAT}" ]
then
   log_and_run b2tufx.log $MEMTIME $MPIRUN b2tufx-${KMER_LENGTH} ${INUFX}
   [ -f ${INUFX}.txt.full ] || bail "No ${INUFX}.txt.full file!"
fi

TIGS=UUtigs

log_and_run contigs.log $MEMTIME $UPCRUN ${THREADS} meraculous-${KMER_LENGTH} -N ${CORES_PER_NODE} -i ${INUFX} -m ${MIN_CONTIG_LENGTH} -o ${TIGS}_${OUT} -c 100 -l 1 -d ${MIN_DEPTH_CUTOFF} -B ${BASE_DIR}

if [ ! "$USE_SHM" ]; then
    [ -f ${BASE_DIR}/${rankpath}/${TIGS}_${OUT}_0.fasta ] || bail "meraculous-${KMER_LENGTH} failed to produce ${TIGS}_${OUT}_0.fasta"
    [ -f ${BASE_DIR}/${rankpath}/my${TIGS}_${OUT}_0.txt ] || bail "meraculous-${KMER_LENGTH} failed to produce my${TIGS}_${OUT}_0.txt"
#    [ -f ${BASE_DIR}/${rankpath}/FASTAcontigs_0 ] || bail "meraculous-${KMER_LENGTH} failed to produce FASTAcontigs_0"
fi
[ -f ${PER_THREAD_DIR}/n${TIGS}_${OUT}.txt ] || bail "meraculous-${KMER_LENGTH} failed to produce ${PER_THREAD_DIR}/n${TIGS}_${OUT}.txt"

merDepthPrefix=merDepth_${TIGS}_${OUT}
log_and_run contigMerDepth.log $MEMTIME $UPCRUN ${THREADS} contigMerDepth-${KMER_LENGTH} -i ${INUFX} -k ${KMER_LENGTH} -c ${TIGS}_${OUT} -d ${MIN_DEPTH_CUTOFF} -s 100 -B ${BASE_DIR}

if [ ! "$USE_SHM" ]; then
    [ -f ${BASE_DIR}/${rankpath}/${merDepthPrefix}_0.txt ] || bail "contigMerDepth-${KMER_LENGTH} failed to produce ${merDepthPrefix}_0.txt"
fi

if [ $is_diploid -eq 1 ]
then

  log_and_run contigEndAnalyzer.log $MEMTIME $UPCRUN ${THREADS} contigEndAnalyzer-${KMER_LENGTH} -i ${INUFX} -k ${KMER_LENGTH} -c ${TIGS}_${OUT} -d ${MIN_DEPTH_CUTOFF} -s 100 -B ${BASE_DIR}

  if [ ! "$USE_SHM" ]; then
    [ -f ${BASE_DIR}/${rankpath}/${TIGS}_${OUT}_0.cea ] || bail "contigEndAnalyzer-${KMER_LENGTH} failed to produce ${TIGS}_${OUT}_0.cea"
  fi

  nContigs=$(cat ${PER_THREAD_DIR}/n${TIGS}_${OUT}.txt)
  bubble_input=${TIGS}_${OUT}
  OUT=diplotigs
  TIGS=Bubbletigs
  log_and_run bubbleFinder.log $MEMTIME $UPCRUN ${THREADS} bubbleFinder-${KMER_LENGTH} -N ${nContigs} -d ${bubble_input} -c ${bubble_input} -o ${TIGS}_${OUT} -B ${BASE_DIR} -D ${BUBBLE_MIN_DEPTH_CUTOFF:=0}

  if [ ! "$USE_SHM" ]; then
    [ -f ${BASE_DIR}/${rankpath}/${TIGS}_${OUT}_0.fasta ] || bail "bubbleFinder-${KMER_LENGTH} failed to produce bubble_${TIGS}_${OUT} "
    [ -f ${BASE_DIR}/${rankpath}/mapping_${TIGS}_${OUT}_0.txt ] || bail "bubbleFinder-${KMER_LENGTH} failed to produce mapping_${TIGS}_${OUT}_0.txt"
    [ -f ${BASE_DIR}/${rankpath}/unfiltered_${TIGS}_${OUT}_0.fasta ] || bail "bubbleFinder-${KMER_LENGTH} failed to produce unfiltered_${TIGS}_${OUT}_0.fasta"
  fi
  [ -f ${PER_THREAD_DIR}/n${TIGS}_${OUT}.txt ] || bail "bubbleFinder-${KMER_LENGTH} failed to produce ${PER_THREAD_DIR}/n${TIGS}_${OUT}.txt"
  ln -f ${PER_THREAD_DIR}/contig_histogram.txt && ln -f ${PER_THREAD_DIR}/contig_histogram_diplotigs.txt && ln -f ${PER_THREAD_DIR}/contig_histogram_isotigs.txt  || bail "bubbleFinder failed to produce the 3 histograms!"

  
  merDepthPrefix=merDepth_${TIGS}_${OUT}
  # merDepth is now calculated in bubbleFinder
  if [ ! "${USE_SHM}" ]; then
    [ -f ${BASE_DIR}/${rankpath}/${merDepthPrefix}_0.txt ] || bail "bubbleFinder-${KMER_LENGTH} failed to produce the new merdepth file(s): ${merDepthPrefix}_0.txt"
  fi

else
  echo "# Skipping Bubble finding as this is not diploid"
fi

nContigs=$(cat ${PER_THREAD_DIR}/n${TIGS}_${OUT}.txt)
if [ -z "${NO_CONCAT}" ]
then
    if [ ${NODES} -gt 16 ]
    then
        echo "Skipping canonical contigs as it is not efficient to run at scale"
    else
        log_and_run upc-canonical-contigs.log $UPCRUN ${THREADS} upc-canonical-assembly -o $BASE_DIR -s $nContigs -n $CORES_PER_NODE -f ${TIGS}_${OUT}_ -F .fasta -O canonical_contigs.fa
        [ -f canonical_contigs.fa ] || bail "upc-canonical-assembly failed to create canonical_contigs.fa"
    fi
fi

# Now run meraligner, histogrammer and splinter on all the libraries

bmaMeta=""
nTotalAlignments=0
list_revcomp=
list_5pwiggle=
list_3pwiggle=
max_readlen=0
for x in config.LIBSEQ*
do
    get_libseq_vars $x

    which merAligner-${KMER_LENGTH} > /dev/null || bail "Could not find merAligner-${KMER_LENGTH} have you compiled with -DKMER_LENGTH=${KMER_LENGTH}?"

    cache_contig_length=400
    # make sure the -e parameter is at least 3 * the overlap
    min_cache_contig_length=$((3 * 2 * (${libreadlen} + 20 - ${KMER_LENGTH})))
    if [ ${cache_contig_length} -lt ${min_cache_contig_length} ] 
    then
        cache_contig_length=${min_cache_contig_length}
    fi

    log_and_run merAligner-${libname}.log $MEMTIME $UPCRUN ${THREADS} merAligner-${KMER_LENGTH} -N ${CORES_PER_NODE} -j 1.2 -r ${libname}.fofn -l ${libname} -P ${libfilesperpair} -c ${TIGS}_${OUT} -C ${HIPMER_SW_CACHE} -k 16228 -x 200000 -m ${KMER_LENGTH} -e ${cache_contig_length} -B ${BASE_DIR} -L ${libreadlen}

    [ -f ${PER_THREAD_DIR}/${libname}-nTotalAlignments.txt ] || bail "merAligner-${KMER_LENGTH} failed to produce ${libname}-nTotalAlignments.txt"

    if [ ! "$USE_SHM" ]; then
        [ -f ${BASE_DIR}/${rankpath}/${libname}-merAlignerOutput_0_Read1 ] || bail "merAligner-${KMER_LENGTH} failed to produce ${BASE_DIR}/${libname}-merAlignerOutput_0_Read1"
    fi

    nTotalAlignments=$((nTotalAlignments+$(cat ${PER_THREAD_DIR}/${libname}-nTotalAlignments.txt)))

    if [ $libgapclosing -eq 1 ]
    then
        list_revcomp="$list_revcomp$librevcomp,"
        if [ $libreadlen -eq 0 ]; then
            libreadlen=$(sort -n ${PER_THREAD_DIR}/${libname}-readlen.txt | tail -1)
            echo "# Found libreadlen $libreadlen"
        fi
        if [ $max_readlen -lt $libreadlen ]
        then
             max_readlen=$libreadlen
        fi
        list_5pwiggle="$list_5pwiggle$lib5pwiggle,"
        list_3pwiggle="$list_3pwiggle$lib3pwiggle,"
    fi

    if [ $libsplinting -eq 1 ]
    then
        opts="-l ${libname} -m ${KMER_LENGTH} -r ${libreadlen}"
        [ ${lib3pwiggle} -eq 0 ] || opts="$opts -T ${lib3pwiggle}" # 3 prime wiggle
        [ ${lib5pwiggle} -eq 0 ] || opts="$opts -F ${lib5pwiggle}" # 5 prime wiggle
        [ $THREADS -gt 960 ] && opts="$opts -N $CORES_PER_NODE"
        log_and_run splinter-${libname}.log $MEMTIME $UPCRUN ${THREADS} splinter ${opts} -B ${BASE_DIR}
        if [ ! "$USE_SHM" ]; then
            [ -f ${BASE_DIR}/${rankpath}/${libname}-splints_0 ] || bail "splinter did not create its files"
        fi
        [ -f "${PER_THREAD_DIR}/$libname-bmaMeta-splints" ] || bail "splinter did not create $libname-bmaMeta-splints"
        [ "$bmaMeta" == "" ] || bmaMeta+=","
        bmaMeta+="$libname-bmaMeta-splints"
        # now fix the splints file to have the insert values other than 0, although that shouldn't 
        # matter because these values should never get used
        echo -e "${libinsavg}\t${libinsstddev}\t0\t${libreadlen}\t${libname}-splints" > ${PER_THREAD_DIR}/${libname}-bmaMeta-splints
    fi
done

[ "$bmaMeta" == "" ] && echo "# There were no libraries used for splinting!"


#
# loop spanner (optional in order of onoSetIds)
#      bma2Links (splint + spanner)
#      ono
#
list_insertSize=
list_insertSigma=
lastSRF=
bestP=
onoSetId=1
prevRoundMaxInsert=0
while [ 1 ]
do

    nLinks=0
    nContigs=$(cat ${PER_THREAD_DIR}/n${TIGS}_${OUT}.txt)
    maxInsertSize=0
    for x in config.LIBSEQ*
    do
        get_libseq_vars $x
        if [ ${onoSetId} -eq ${libonosetid} ] && [ ${libfilesperpair} -gt 0 ]
        then
            [ -z "${lastSRF}" ] || nScaffolds=$(cat ${PER_THREAD_DIR}/nScaffolds_$lastSRF.txt)
            # use histogrammer to get the insert and stddev sizes
            minThreads=${CORES_PER_NODE}
            [ "$NODES" -gt 1 ] && minThreads=$((2*$CORES_PER_NODE))
            threadsToUsePerNode=$(($minThreads/$NODES))
            [ "$threadsToUsePerNode" -lt 1 ] && threadsToUsePerNode=1
            totThreadsToUse=$(($threadsToUsePerNode*$NODES))
            numAlignsPerThread=$(($NUM_ALIGNMENTS/$totThreadsToUse))
            opts="-l ${libname} -m ${KMER_LENGTH} -i ${libinsavg} -s ${libinsstddev} -c ${numAlignsPerThread} "
            [ ${librevcomp} -eq 0 ] || opts="$opts -R" # reverse complement
            [ ${libinnie} -eq 0 ] || opts="$opts -A" # innie removal
            [ ${lib3pwiggle} -eq 0 ] || opts="$opts -T ${lib3pwiggle}" # 3 prime wiggle
            [ ${lib5pwiggle} -eq 0 ] || opts="$opts -F ${lib5pwiggle}" # 5 prime wiggle
            [ -z "${lastSRF}" ] || opts="${opts} -Z ${lastSRF} -C $nContigs -G $nScaffolds"
            log_and_run histogrammer-${libname}.log $MEMTIME $UPCRUN ${THREADS} merAlignerAnalyzer ${opts} -B ${BASE_DIR} -N $CORES_PER_NODE -p ${threadsToUsePerNode}
            ( [ -f ${PER_THREAD_DIR}/${libname}-insert.txt ] && [ -f ${PER_THREAD_DIR}/${libname}-std.txt ] ) || bail "merAlignerAnalyzer did not make its files"

            insertSize=`cat ${PER_THREAD_DIR}/${libname}-insert.txt`
            insertSigma=`cat ${PER_THREAD_DIR}/${libname}-std.txt`
            [ "$insertSize" -lt "$maxInsertSize" ] || maxInsertSize=$insertSize
            
            if [ $libgapclosing -eq 1 ]; then
                list_insertSize="$list_insertSize$insertSize,"
                list_insertSigma="$list_insertSigma$insertSigma,"
            fi

            # run spanner on this libary now
            opts="-l ${libname} -m ${KMER_LENGTH}"
            [ ${librevcomp} -eq 0 ] || opts="$opts -R" # reverse complement
            [ ${libinnie} -eq 0 ] || opts="$opts -A" # innie removal
            [ ${lib3pwiggle} -eq 0 ] || opts="$opts -T ${lib3pwiggle}" # 3 prime wiggle
            [ ${lib5pwiggle} -eq 0 ] || opts="$opts -F ${lib5pwiggle}" # 5 prime wiggle
            opts="${opts} -i ${insertSize} -s ${insertSigma} -r ${libreadlen}"
            [ -n "${lastSRF}" ] && opts="${opts} -S ${lastSRF} -C $nContigs -Z $nScaffolds"
            log_and_run spanner-${libname}-${onoSetId}.log $MEMTIME $UPCRUN ${THREADS} spanner ${opts} -B ${BASE_DIR}
            if [ ! "$USE_SHM" ]; then
                [ -f ${BASE_DIR}/${rankpath}/${libname}-spans_0 ] || bail "spanner did not write spans files!"
            fi
            [ -f ${PER_THREAD_DIR}/${libname}-bmaMeta-spans ] || bail "spanner did not write spans bmaMeta file!"

            [ "$bmaMeta" != "" ] && bmaMeta+=","
            bmaMeta+="${libname}-bmaMeta-spans"
        fi
    done

    
    # bmaToLinks

    log_and_run bmaToLinks-${onoSetId}.log $MEMTIME $UPCRUN ${THREADS}  bmaToLinks -b ${bmaMeta} -m ${KMER_LENGTH} -L ${onoSetId} -n ${nTotalAlignments} -B ${BASE_DIR} -P ${prevRoundMaxInsert}
    if [ ! "$USE_SHM" ]; then
        [ -f ${BASE_DIR}/${rankpath}/LINKS_OUTPUT_${onoSetId}_0 ] || bail "bmaToLinks did not create LINKS_OUTPUT_${onoSetId}_0!"
    fi
    [ -f ${PER_THREAD_DIR}/nLinks-${onoSetId}.txt ] || bail "bmaToLinks did not create nLinks-${onoSetId}.txt!"
    [ -f ${PER_THREAD_DIR}/linksMeta-${onoSetId} ] || bail "bmaToLinks did not create linksMeta-${onoSetId}!"

    # unset bmaMeta for the next round
    bmaMeta=""
    nLinks=$((nLinks+$(cat ${PER_THREAD_DIR}/nLinks-${onoSetId}.txt)))
    thisSRF=SsrfFile-${onoSetId}
    
    # oNo
    # create merged files to enable multiple simultaneous runs of oNo
    nOnoRuns=8
    # nOnoRuns must be a factor of THREADS, check for this
    if [ $((THREADS%nOnoRuns)) -ne 0 ]
    then
      echo "Overriding CAN_SPLIT_JOB, as THREADS(${THREADS}) is not a factor of nOnoRuns (${nOnoRuns})"
      CAN_SPLIT_JOB=0
    fi
    minP=3
    threadsPerRun=${THREADS}
    if [ "${CAN_SPLIT_JOB}" == "1" ]
    then
      # threadsPerRun must be a multiple of CORES_PER_NODE
      threadsPerRun=$((THREADS/nOnoRuns))
      threadsPerRun=$((CORES_PER_NODE * (threadsPerRun / CORES_PER_NODE)))
      if [ ${threadsPerRun} -lt ${CORES_PER_NODE} ]
      then
        echo "Overriding threadsPerRun (from ${threadsPerRun}) to ${CORES_PER_NODE}"
        threadsPerRun=${CORES_PER_NODE}
        ONO_ACTIVE_THREADS=${CORES_PER_NODE}
      fi
      # is this still true?? [ $(($THREADS%$nOnoRuns)) -eq 0 ] || bail "nOnoRuns (${nOnoRuns}) must be a factor of THREADS (${THREADS})"
    fi
    if [ -n "$USE_SHM" ]; then
        # if using shm,  the sets of threads for each oNo run must fit on non-overlapping nodes 
        if [ $threadsPerRun -lt $CORES_PER_NODE ]; then
            [ $(($CORES_PER_NODE%$threadsPerRun)) -eq 0 ] || bail "threadsPerRun ($threadsPerRun) must be a factor of cores per node (${CORES_PER_NODE})"
        else 
            [ $(($threadsPerRun%$CORES_PER_NODE)) -eq 0 ] || bail "cores per node (${CORES_PER_NODE}) must be a factor of threadsPerRun ($threadsPerRun)"
        fi
    fi
    if [ ${ONO_ACTIVE_THREADS} -gt ${threadsPerRun} ]
    then
       echo "Overriding ONO_ACTIVE_THREADS to ${threadsPerRun}"
       ONO_ACTIVE_THREADS=${threadsPerRun}
    elif [ $((threadsPerRun % ONO_ACTIVE_THREADS)) -ne 0 ]
    then
       y=1
       for x in $(seq 1 ${ONO_ACTIVE_THREADS})
       do
         if [ $((threadsPerRun % x)) -eq 0 ]
         then 
           y=${x}
         fi
       done
       echo "Overriding ONO_ACTIVE_THREADS from ${ONO_ACTIVE_THREADS} to ${y}"
       ONO_ACTIVE_THREADS=${y}
    fi
    [ $((threadsPerRun % ONO_ACTIVE_THREADS)) -eq 0 ] || bail "ONO_ACTIVE_THREADS (${ONO_ACTIVE_THREADS}) must be a factor of threadsPerRun (${threadsPerRun})"

    echo "# Launching $nOnoRuns oNo runs with $threadsPerRun threads each with ${ONO_ACTIVE_THREADS} active threads"

    opts="-l merged_linksMeta-${onoSetId}  -m ${KMER_LENGTH} -L ${nLinks} -C ${nContigs} "
    if [ -n "$lastSRF" ]
    then
        nScaffolds=$(cat ${PER_THREAD_DIR}/nScaffolds_${lastSRF}.txt)
        # the lastSRF is with the best p value from the previous round
        opts="$opts -n ${nScaffolds} -s ${lastSRF}-${bestP} "
    else
        log_and_run merger-merDepth-${onoSetId}.log $MEMTIME $UPCRUN ${threadsPerRun} merger -P ${THREADS} -f ${merDepthPrefix} -s .txt
        opts="$opts -c merged_${merDepthPrefix} "
    fi

    linksInsert=`awk '{print $1}' ${PER_THREAD_DIR}/linksMeta-${onoSetId}`
    linksStd=`awk '{print $2}' ${PER_THREAD_DIR}/linksMeta-${onoSetId}`
    linksPrefix=`awk '{print $3}' ${PER_THREAD_DIR}/linksMeta-${onoSetId}`
    echo -e "$linksInsert\t$linksStd\tmerged_${linksPrefix}" > ${PER_THREAD_DIR}/merged_linksMeta-${onoSetId}

    log_and_run merger-links-${onoSetId}.log $MEMTIME $UPCRUN ${threadsPerRun} merger -P ${THREADS} -f ${linksPrefix}

    secs=$SECONDS
    echo "# Starting parallel ($CAN_SPLIT_JOB) oNo-${onoSetId} at $(date)"
    echo "# Starting parallel ($CAN_SPLIT_JOB) oNo-${onoSetId} at $(date)" >> ${TIMINGS}

    nodesPerRun=$((${NODES}/${nOnoRuns}))
    if [ "$nodesPerRun" -eq 0 ]; then
        nodesPerRun=1
    fi
    upcrunNodesOpt=""
    if [ "${CAN_SPLIT_JOB}" == "1" ]
    then
      upcrunNodesOpt="-cpus-per-node=${CORES_PER_NODE} -nodes=0"
    fi
    # parallel launch of ono to find best value of p 
    for p in `seq ${minP} $((nOnoRuns+minP-1))`
    do 
        # launch as background process - note the ampersand!
        # adding a second -n will override the first one which is part of $UPCRUN
        RETRY=3 log_and_run oNo-${onoSetId}-${p}.log $MEMTIME $UPCRUN ${threadsPerRun} ${upcrunNodesOpt} oNo ${opts} -N ${ONO_ACTIVE_THREADS} -B ${BASE_DIR} -p ${p} -o ${thisSRF}-${p} &
        [ "${CAN_SPLIT_JOB}" == "1" ] || wait || bail "Could not execute oNo for p ${p}"
        if [ "${CAN_SPLIT_JOB}" != "1" ] && [ $((nodesPerRun * nOnoRuns)) -le ${NODES} ]
        then
            sleep 1 # attempt to rate limit the parallel submission
        fi
    done

    # wait for all background processes to complete
    [ "${CAN_SPLIT_JOB}" != "1" ] || wait || bail "Could not execute all of the oNo jobs"
    runsecs=$((SECONDS-$secs))
    echo "# Finished parallel oNo-${onoSetId} at $(date) in ${runsecs} s"
    echo "# parallel_oNo-${onoSetId} ${runsecs} walltime-secs $((runsecs*THREADS)) core-seconds at $(date)" >> ${TIMINGS}

    # now check outputs
    for p in `seq ${minP} $((nOnoRuns+minP-1))`
    do 
        if [ ! "$USE_SHM" ]; then
            [ -f ${BASE_DIR}/${rankpath}/${thisSRF}-${p}_0 ] || bail "oNo failed to produce ${thisSRF}-${p}_0"
        fi
        [ -f ${PER_THREAD_DIR}/nScaffolds_${thisSRF}-${p}.txt ] || bail "ono failed to produce nScaffolds_${thisSRF}-${p}.txt"
        [ -f ${PER_THREAD_DIR}/nGaps_${thisSRF}-${p}.txt ] || bail "ono failed to produce nGaps_${thisSRF}-${p}.txt "
        [ -f ${PER_THREAD_DIR}/n50_${thisSRF}-${p}.txt ] || bail "ono failed to produce n50_${thisSRF}-${p}.txt "
    done

    # get the best N50
    bestP=`cat ${PER_THREAD_DIR}/n50_${thisSRF}-* |sort -nr -k2 |head -n1|awk '{print $1}'`
    bestN50=`cat ${PER_THREAD_DIR}/n50_${thisSRF}-${bestP}.txt|awk '{print $2}'`

    echo "# Selected p = ${bestP} for round ${onoSetId}, N50 ${bestN50}"

    # split the best ssrf files produced by oNo back to one per core - they are used by scaffolding stages other than oNo
    log_and_run splitter-${onoSetId}.log $MEMTIME $UPCRUN ${threadsPerRun} splitter -P ${THREADS} -s split -f ${thisSRF}-${bestP}
    # rename the files so they can be used by the other stages 
    for i in `seq 0 $((THREADS-1))`
    do
        rd=$(get_rankpath ${i})
        ln -f ${rd}/${thisSRF}-${bestP}_split_${i} ${rd}/${thisSRF}_${i}
    done
    ln -f ${PER_THREAD_DIR}/nScaffolds_${thisSRF}-${bestP}.txt ${PER_THREAD_DIR}/nScaffolds_${thisSRF}.txt 
    ln -f ${PER_THREAD_DIR}/nGaps_${thisSRF}-${bestP}.txt ${PER_THREAD_DIR}/nGaps_${thisSRF}.txt 

    lastSRF=${thisSRF}

    onoSetId=$((onoSetId+1))
    [ ${onoSetId} -gt ${MAX_ONO_SET_ID} ] && break

    prevRoundMaxInsert=${maxInsertSize}

done  # onoSetId loop



if [ ! "$USE_SHM" ]; then
    [ -f "${BASE_DIR}/${rankpath}/${lastSRF}_0" ] || bail "Could not find the srf!"
fi

if [ -n "${GAPLIBS}" ]
then
  opts="-Q ${qual_offset} -s ${lastSRF} -c ${TIGS}_${OUT} -a merAlignerOutput -b ${GAPLIBS} -i ${list_insertSize} -I ${list_insertSigma} -m ${KMER_LENGTH} -D ${MIN_DEPTH_CUTOFF} -R $EXCLUDE_REPEATS -A $POLYMODE -B $BASE_DIR -r $list_revcomp -F $list_5pwiggle -T $list_3pwiggle -l $max_readlen "
  log_and_run gapclosing.log $MEMTIME $UPCRUN ${THREADS} gapclosing ${opts} 
  if [ ! "$USE_SHM" ]; then
    [ -f $BASE_DIR/${rankpath}/assembly.0.fa ] || bail "gapclosing failed to produce assembly"
    [ -f $BASE_DIR/${rankpath}/assembly-1000.0.fa ] || bail "gapclosing failed to produce assembly-1000"
  fi
else
  bail "TODO: generate scaffolds without gap closing"
fi

nScaffolds=$(cat ${PER_THREAD_DIR}/nScaffolds_SsrfFile-1.txt)
log_and_run upc-canonical-assembly.log $UPCRUN ${THREADS} upc-canonical-assembly -o $BASE_DIR -s $nScaffolds -n $CORES_PER_NODE -O final_assembly.fa
[ -f final_assembly.fa ] || bail "upc-canonical-assembly failed to produce final_assembly.fa"

log_and_run upc-canonical-assembly-1000.log $UPCRUN ${THREADS} upc-canonical-assembly -o $BASE_DIR -s $nScaffolds -n $CORES_PER_NODE -f assembly-1000. -F .fa -O final_assembly-1000.fa
[ -f final_assembly-1000.fa ] || bail "upc-canonical-assembly failed to produce final_assembly-1000.fa"

#if [ "$USE_SHM" ]
#then
#    find /dev/shm/ -maxdepth 1 -type f -delete
#fi

echo "# Total ${SECONDS} walltime-sec $((SECONDS*THREADS)) core-secs at $(date)" >> $TIMINGS

if [ "$RUNDIR" != "." ]; then
    cd ..
    rm -f latest_output
    echo "adding link latest_output->${RUNDIR}"
    ln -s ${RUNDIR} latest_output || /bin/true
fi

if [ -n "${POST_RUN}" ]
then
  echo "Running POST_RUN: ${POST_RUN}"
  ( source ${POST_RUN} )
  echo "# After POST_RUN ${SECONDS} walltime-sec $((SECONDS*THREADS)) core-secs at $(date)" >> $RUNDIR/$TIMINGS
fi

echo "# Finished at $(date) in $SECONDS s: $0 $@ (${SLURM_JOB_ID}${PBS_JOBID}${JOB_ID} pid:$$)"
