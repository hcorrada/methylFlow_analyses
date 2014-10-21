# How to make things work

## Directory layout

Setup the methylflow project directory as follows:

<pre>
parent_dir
| exps (this will contain this repository)
| src (this will contain the source repo)
| install (this is where all binaries, libraries and headers are installed)
</pre>

To set this up, go to your project directory
and run

```bash
git co https://github.com/hcorrada/methylFlow_analyses.git exps  
mkdir src && cd src  
git co https://github.com/hcorrada/methylFlow.git  
```

## Setting up

First, set the location of the parent directory and install
directories variables (starting from project directory)

```bash
cd exps
source setup.sh
```

## Compiling and installing methyl flow

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


