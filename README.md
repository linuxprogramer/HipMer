______________________________________________________________________________

    HipMer v 1.0, Copyright (c) 2016, The Regents of the University of California,
    through Lawrence Berkeley National Laboratory (subject to receipt of any
    required approvals from the U.S. Dept. of Energy).  All rights reserved.
 
    If you have questions about your rights to use or distribute this software,
    please contact Berkeley Lab's Innovation & Partnerships Office at  IPO@lbl.gov.
 
    NOTICE.  This Software was developed under funding from the U.S. Department
    of Energy and the U.S. Government consequently retains certain rights. As such,
    the U.S. Government has been granted for itself and others acting on its behalf
    a paid-up, nonexclusive, irrevocable, worldwide license in the Software to
    reproduce, distribute copies to the public, prepare derivative works, and
    perform publicly and display publicly, and to permit other to do so. 

______________________________________________________________________________

    This is HipMer, the high performance distributed memory scalable version of Meraculous.

    HipMer is a high performance parallelization and port of Meraculous 
    (http://jgi.doe.gov/data-and-tools/meraculous/).  It is largely written in UPC,
    with the exception of the UFX generation, which is written in C++/MPI.

    This project is a joint collaboration between JGI (http://jgi.doe.gov), 
    NERSC (http://www.nersc.gov/) and CRD (http://crd.lbl.gov/)

    Primary authors are:
    Evangelos Georganas, Aydin Buluc, Steven Hofmeyr, Leonid Oliker and Rob Egan, 
    with direction and advice from Kathy Yelick.
    The original Meraculous was developed by Jarrod Chapman, Isaac Ho, Eugene Goltsman,
    and Daniel Rokhsar.


Building and installing
-----------------------

You can download the source from sourceforge: https://sourceforge.net/projects/hipmer/

To build, install and run test cases, use scripts from the appropriate 
.platform_deploy, where 'platform' is one of several different platforms, e.g. 
'.edison_deploy' for NERSC's Edison system, '.cori_deploy' for NERSC's cori 
system, and '.generic' for a generic Linux system.

You should specify the environmental variable SCRATCH for default placement 
of the build and install paths. You can also change the default build and install 
paths by overriding the environmental variables: BUILD and PREFIX respectively

To build:

    .platform_deploy/build.sh


To install:

    .platform_deploy/install.sh

By default, the build will be in $SCRATCH/build-platform and the install will 
be in $SCRATCH/install-platform.

There are environmental variables that are automatically set for a release 
(non-debug) build (.platform_deploy/env.sh). To build a debug version, first 
execute:

    source .platform_deploy/env-debug.sh

Then run .platform_deploy/build.sh && .platform_deploy/install.sh as before.

By default, the debug build will be in $SCRATCH/build-platform-debug, and the 
install in $SCRATCH/install-platform-debug.

To force a complete rebuild:

    CLEAN=1 .platform_deploy/build.sh

To force a rebuild with all the environment checks:

    DIST_CLEAN=1 .platform_deploy/build.sh

Note that running .platform_deploy/install.sh should do partial rebuilds for 
changed files.

WARNING: the build process does not detect header file dependencies for UPC 
automatically, so changes to header files will not necessarily trigger 
rebuilds. The dependencies need to be manually added. This has been done for 
some, but not all, stages.

Some features of the cmake build process:
* Builds multiple binaries based on the build parameters:
    * export HIPMER_BUILD_OPTS="-DHIPMER_KMER_LENGTHS='21;51'"
* Properly builds UPC source (if you name the source .upc or set the LANGUAGE 
and LINKER_LANGUAGE property to UPC)
    * Sets the -D definition flags consistently 
    * Supports -DCMAKE_BUILD_TYPE=Release or Debug


Running
-------

To run, use the src/hipmer/run_hipmer.sh script which requires the install 
directory and meraculous.config file:

    ${PREFIX}/bin/run_hipmer.sh ${PREFIX} meraculous.config

There are several configuration files in test/pipeline/*.config:

| Config File                     | Description                         |
|:--------------------------------|:------------------------------------|
| meraculous-validation.config    | a small validation test             |
| meraculous-ecoli.config         | ecoli dataset, easy to run on single node systems with limited cores & memory |
| meraculous-chr14.config         | human chromosome 14 (diploid). Can be run on single node systems, but will be slower than ecoli. |
| meraculous-human.config         | full human dataset, requires around 1TB of memory |


For convenience, there are run scripts in .platform_deploy that make it easier 
to run jobs. For systems like Edison, these scripts can be submitted directly 
to the job queue (with overridden queue, mppwidth and wall time options). For 
.generic_deploy, the scripts can be executed directly. The scripts expect the 
data for the tests to be in hipmer_name_data, where name is 'ecoli', 
'validation', 'human', 'chr14'. If a dataset doesn't exist, the script will 
download and install it (with a stripe of 72 on Edison).


The run scripts automatically set the CORES_PER_NODE and THREADS variables. In 
the case of Edison, the number of threads is determined from the number of 
processors found at runtime, and the CORES_PER_NODE is fixed to 24. In the case 
of generic, the number of threads by default is all those available on the 
single node, and the CORES_PER_NODE is the same value. You can override these 
values, but make sure to set the CORES_PER_NODE appropriately if you change 
anything.


The pipeline can also be run with all the intermediate per-thread files in 
shared memory (/dev/shm), plus the FASTQ inputs. Set the environment variable 
USE_SHM=1 to achieve this. There are some scripts that have the shared memory 
option, e.g. .edison_deploy/test_hipmer_human-edison-shm.sh. Using shared 
memory will be faster, especially at larger concurrencies, but require a lot 
more memory (around 1TB for human).


Before launching run_hipmer.sh, some settings can be changed through 
environmental variables:

    Physical memory properties:
      PHYS_MEM_MB=413135
    
    The number of UPC threads in the job (MPI will not use hyperthreads):
      THREADS=20 HYPERTHREADS=1
    
    The rundirectory to place all the files (will be created if necessary):
      RUNDIR=/global/homes/r/regan/workspace/hipmeraculous/cori09-116095-20160309-160941
      DATADIR=
      LFS_STRIPE=72
    
    Gasnet properties (by default 80% of physmem):
      GASNET_PHYSMEM_MAX=413135MB
      GASNET_PHYSMEM_NOPROBE=1
     
    UPC properties:
      UPC_SHARED_HEAP_MB=15500 (Do not set to use 80% of the node memory)
      UPC_PTHREADS=
      UPC_PTHREADS_OPT=
    
    HipMer options (will override config file defaults):
      UFX_HLL=0
      MIN_DEPTH_CUTOFF=0 # use 0 for auto-detect after UFX generation
      BUBBLE_MIN_DEPTH_CUTOFF=
      KMER_LENGTH=
      MIN_CONTIG_LENGTH=
      NUM_ALIGNMENTS=192000
      ILLUMINA_VERSION=
      HIPMER_SW_CACHE=6144
      NO_CONCAT=
      CAN_SPLIT_JOB=0
      ONO_ACTIVE_THREADS=20
    
    MPI/UPC environment:
      MPIRUN=mpirun -n 20
      UPCRUN=upcrun    -shared-heap=15500M  -n


Note: the Illumina Version is automatically detected, and will be reported in 
the output for the run. To override, set the ILLUMINA_VERSION environment 
variable.

Some features of run_hipmer.sh:

  * Should detect and run the alternate diploid workflow if set in the config 
file
  * Runs the proper set and ordering of splinter, spanner, bmaToLinks and oNo 
for all libraries based on the config file
  * Organizes inputs by library, finds the files specified in the config file
  * Calls the proper version of binaries (KMER_LENGTH & READ_LENGTH)
  * Logs all commands and timings in timings.log
  * Aborts on error, continues on first failed step
  * Validates that outputs are generated

To rerun specific stages in the pipeline, first delete the .log file for the 
stage or stages, and then execute:

    export RUNDIR=<name of output dir>
    ${PREFIX}/bin/run_hipmer.sh <install_path> <output_dir>/<config_file>

Other helper scripts:

    ${PREFIX}/bin/rerun_stage.sh <stage>

Execute from within the output dir for a run, and it will scan the 
timings.log to determine what stages have run. Without any arguments, it will 
show a list of stages, one or more of which can then be passed in a comma 
separated list to rerun those stages. The old files will be overwritten, except 
for timings.log, which will not be affected. Instead, the results are appended 
to a new file, timings-rerun.log.

    ${PREFIX}/bin/rerun_single_stage.sh <stage>

Similar to the rerun_stage.sh, execute this from within an output dir. 
However, when run with a stage name, it will simple bring up the command line 
to execute that stage without actually running it, so you can edit the command 
lin and then hit enter to execute it. Also, it will not update the 
<command>.log file or the timings.log, but it will change any files that 
running that stage will normally change. This script needs to be run from 
within an interactive session.

    ${PREFIX}/bin/compare_results.sh <dir1> <dir2>

Pass in two different directories and it will compare the outputs using a 
number of statistics. Because of non-determinism in the scaffolding process, 
this is the best we can do for checking for similarities. 

    ${PREFIX}/bin/get-stage-timings.sh

Execute from within the output dir for a run and it will extract the time 
taken by each stage, both the internal and the overall time (including the job 
launch time).


Workflow
--------

The HipMer workflow is controlled within the configuration file, when the 
libraries are specified. For each library, you can specify what round of oNo to 
use it in, and you can specify whether or not to use it for splinting. The 
workflow is as follows (see run_hipmer.sh for details):

1. (prepare input fastq files)
  1. They must be uncompressed
  2. They ought to be striped for efficient parallel access
2. prepare meraculous.config
3. ufx
4. contigs
5. contigMerDepth
6. if diploid:
    1. contigEndAnalyzer
    2. bubbleFinder
7. (optionally upc_canonical_assembly: canonical_contigs.fa)
8. for each library:
    1. merAligner
    2. splinter (if specified in config file)
9. for each oNoSetID:
    1. for each library in that oNoSetId:
        1. merAlignerAnalyzer (histogrammer)
        2. spanner
    2. bmaToLinks
    3. merger
    4 for each oNoRuns choose a -p
        1. oNo
    5. splitter
10. gapclosing
11. upc_canonical_assembly: final_assembly.fa

This means that the first round of bmaToLinks could end up processing the 
outputs from multiple iterations of splinter plus multiple ones of spanner. The 
subsequent calls to bmaToLinks will only process outputs from spanner.