export MF_PROJECT_ROOT=${PWD}/..
export MF_INSTALL_DIR=${MF_PROJECT_ROOT}/install

if [[ $PS1 =~ "mfset" ]]; then
    echo 'already set'
else
    export PS1="(mfset)$PS1"
fi

