#	Lattice Blind Signature
This repository contains the code for the paper **Improved Lattice Blind Signatures from Recycled Entropy**. Part of the code was taken from the implementation of [lattice anonymous credentials](https://github.com/Chair-for-Security-Engineering/lattice-anonymous-credentials).

The repository contains
- the parameter selection and security estimation scripts, and
- the implementation in C

## Parameter Scripts

In the `scripts` subdirectory, we provide a parameter estimation script which depends on the lattice-estimator.
The lattice-estimator is included as submodule, so please do a recursive pull.
```shell
git submodule update --init --recursive
```

You can then run the python script `parameters.py`

## C Implementation

Requirements:
- FLINT, which depends on GMP and MPFR
- AES-NI instructions
- cmake

FLINT is included as submodule to this repository, so please do a recursive pull.
To install FLINT, run the following commands (from the root directory of this repository). It only needs to be built once.
```shell
cd code/libs/flint
mkdir _build
cd _build
cmake .. -DBUILD_SHARED_LIBS=ON
cmake --build . --target install
cd ..
./bootstrap.sh 
./configure --enable-avx2 --disable-pthread CFLAGS="-O3 -Wall -march=native" CC="gcc-11"
make -j
```

To build our implementation, run the following commands (from the root directory of this repository).
```shell
cd code
mkdir build
cd build
cmake ..
make -j
```

Then, running `./test` or `./bench` runs the tests or the benchmarks, respectively. Each change of the source code requires to run the last command (`make` or `make -j`).


## Building with Docker
For information on how to install docker on your system visit the [docker docs](https://docs.docker.com/).
Depending on your setup, you may need to prefix the following commands with `sudo`.

### 1. Build the docker image
Navigate to the top-level directory which contains this README file. Then run (notice the trailing dot!)
```shell
docker build -t lbsig-docker .
```
This takes some time as docker needs to pull the base-image for Ubuntu and then builds all dependencies
(i.e. GMP, MPFR and FLINT) as well as the test and benchmark executables from scratch.
In particular the `make check` instructions are *very* time consuming.

The libraries are built with standard options and fine-tuning is possible by adjusting the corresponding
`make` and `configure` commands in the Dockerfile. However, installation paths *must not* be changed.

Note that this step needs to be repeated for any changes in the code base or Dockerfile but not
for subsequent runs of the code.

### 2. Run the docker image
Issue the following command to start the docker container, where `XXX` is either `./test` or `./bench` to select
the tests or benchmarks respectively.
```shell
docker run -t lbsig-docker XXX
```
You should now see the output of either the test or the benchmarking executables.

### 3. Stopping the running docker image
To get a list of running docker containers use `docker ps`. Then you can use `docker stop <container_id>` to stop
the execution. For more information see the [docker docs](https://docs.docker.com/engine/reference/builder/).

### 4. Cleaning up
To remove all stopped docker containers use `docker container prune`.
To remove docker images use `docker images` to get a list of all images and then `docker rm <image_id>` to remove
the image in question.