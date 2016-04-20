# Plasmodium genome assembler and evaluator

The assembler project is designed to perform all the logical steps of a genome
assembly project, from preprocessing to comparative analysis. The pipeline was
built to document the process of genome assembly for viral genomes.
The assembly pipeline defaults to the best parameters we found through extensively
testing and optimizing each tool, to achieve the best assembly for our sequencing
runs. The evaluation pipeline, in vision, would allow you to automate the process
of evaluating assemblies and fine tuning the parameters to arrive at the best assembly.
The evaluation pipeline, hence would be a more generic pipeline.
This module will eventually be added to the NEST toolkit.
This project is funded by CDC and the pipeline is built from the combined work
at the Vannberg lab at Georgia Institute of Technology and the malaria research group
at CDC.


# Setup:
This sections will go through the setup steps for each of the software used
within this pipeline. This setup procedure will try to install all the software
to the local library, but some system modules may need super user privileges.

1. AbySS:

  1. OpenMPI:
  ```{sh}
  wget https://www.open-mpi.org/software/ompi/v1.10/downloads/openmpi-1.10.2.tar.gz
  tar -zxvf openmpi-1.10.2.tar.gz
  cd openmpi-1.10.2
  ./configure --prefix=~/local
  make all install
  ```
  2. Boost:
  ```{sh}
  wget "http://downloads.sourceforge.net/project/boost/boost/1.60.0/boost_1_60_0.zip?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fboost%2F%3Fsource%3Dtyp_redirect&ts=1458155252&use_mirror=iweb" -O boost_1_60_0.zip
  unzip boost_1_60_0.zip
  cd boost_1_60_0
  ./bootstrap.sh --prefix=~/local
  ./b2 install
  ```

  3. Sparsehash:
  ```{sh}
  git clone https://github.com/sparsehash/sparsehash.git
  cd sparsehash
  ./configure --prefix=~/local
  make
  make install
  ```
  4. AbySS:
  ```{sh}
  wget https://github.com/bcgsc/abyss/releases/download/1.9.0/abyss-1.9.0.tar.gz
  tar -zxvf abyss-1.9.0.tar.gz
  cd abyss-1.9.0
  ./configure --prefix=~/local --with-boost=~/local/include/ --with-mpi=~/local/lib/openmpi
  make
  make install
  ```

2. Spades:

  1. Gcc:
  ```{sh}
  svn checkout svn://gcc.gnu.org/svn/gcc/trunk gcc
  cd gcc
  ./contrib/download_prerequisites
  cd ..
  mkdir objdir
  cd objdir
  ./configure --prefix=~/local --enable-languages=c,c++,fortran,go --disable-multilib
  make
  make install
  ```
  2. Cmake:
  ```{sh}
  wget https://cmake.org/files/v3.5/cmake-3.5.0.tar.gz
  tar -zxvf cmake-3.5.0.tar.gz
  cd cmake-3.5.0
  ./bootstrap --prefix=~/local
  make
  make install
  ```
  2. Spades:
  ```{sh}
   wget http://spades.bioinf.spbau.ru/release3.7.1/SPAdes-3.7.1.tar.gz
   tar -zxvf SPAdes-3.7.1.tar.gz
   cd SPAdes-3.7.1
   PREFIX=~/local ./spades_compile.sh
  ```
3. NGOPT:
  1. NGOPT:
  ```{sh}
  wget https://sourceforge.net/projects/ngopt/files/a5_miseq_linux_20150522.tar.gz/download
  tar -zxvf download
  export PATH=$PATH:~/a5_miseq_linux_20150522/bin
  ```
