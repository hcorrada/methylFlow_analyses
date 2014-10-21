# How to make things work

## Getting started

Create a directory for the methylFlow project. Once there do

```bash
git co https://github.com/hcorrada/methylFlow_analyses.git exps
mkdir src && cd src
git co https://github.com/hcorrada/methylFlow.git
```

## Setting up

Set the location of the parent directory and install
directories variables (starting from project directory)

```bash
cd exps
source setup.sh
```

## Compiling and installing methylFlow

```bash
pushd ${MF_PROJECT_ROOT}/src/methylFlow
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=${MF_INSTALL_DIR} ..
make
make install
popd
```

## Compiling and installing evaluation code

```bash
pushd ${MF_PROJECT_ROOT}/exps/simulations
mkdir build && cd build
cmake -DMF_INSTALL_DIR=${MF_INSTALL_DIR} -DCMAKE_INSTALL_PREFIX=${MF_INSTALL_DIR} ..
make
make install
popd
```

## Directory layout

At this point, the layout should be

<pre>
parent_dir
| exps (this will contain this repository)
| src (this will contain the source repo)
| install (this is where all binaries, libraries and headers were installed)
</pre>

## Running installed programs

```bash
# the methylFlow executable
${MF_INSTALL_DIR}/bin/methylFlow

# the simulator
${MF_INSTALL_DIR}/bin/mfSimulate

## the evaluation programs
${MF_INSTALL_DIR}/bin/mfEvaluation
${MF_INSTALL_DIR}/bin/avgEvaluation
```

## Running on bscapture data

```bash
pushd ${MF_PROJECT_ROOT}/exps/bscapture
sh run_cpgmf.sh 2>&1 >/dev/null
popd
```
