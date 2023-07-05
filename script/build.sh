
# TODO: change the libPath to your specific path
libPath=""

cmakePath="$libPath/cmake-3.17.2"
gccPath="$libPath/gcc-9.1.0"
CLI11Path="$libPath/CLI11-1.9.0"
spdlogPath="$libPath/spdlog-1.5.0"
libdeflatePath="$libPath/libdeflate-1.5"
htslibPath="$libPath/htslib-1.9"
yggPath="$libPath/ygg-master"
taskflowPath="$libPath/taskflow-2.5.0"
sparseppPath="$libPath/sparsepp"
py="$libPath/python3.6/bin/python"

binPath="/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/bin"
export PATH="$gccPath/bin:$cmakePath/bin:$binPath"

export LD_LIBRARY_PATH="$libdeflatePath/lib:$htslibPath/lib:$gccPath/lib64:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/lib64"
export LIBRARY_PATH=$LD_LIBRARY_PATH

export C_INCLUDE_PATH="$sparseppPath/include:$C_INCLUDE_PATH"
export C_INCLUDE_PATH="$libdeflatePath/include:$yggPath/include:$C_INCLUDE_PATH"
export C_INCLUDE_PATH="$htslibPath/include:$gccPath/include:$C_INCLUDE_PATH"
export C_INCLUDE_PATH="$CLI11Path/include:$spdlogPath/include:$C_INCLUDE_PATH"
export C_INCLUDE_PATH="$taskflowPath/include:$C_INCLUDE_PATH"
export CPLUS_INCLUDE_PATH=$C_INCLUDE_PATH

export CC="$gccPath/bin/gcc"
export CXX="$gccPath/bin/g++"

absPath(){
    relativePath=$1
    mkdir -p $relativePath
    if [ ${relativePath:0:1} == "/" ]
    then
        echo $relativePath
        return 0
    fi
    ap="$(cd $relativePath; pwd)"
    echo "$ap"
}

srcPath="$(cd $(dirname $(dirname $0)); pwd)"
buildPath="$(mktemp -d -p /dev/shm $(basename $srcPath).XXXXXXXXXX)"

buildPath="$(absPath $buildPath)"
if [ -e "$buildPath" ]
then
    rm -rf $buildPath
fi
mkdir -p $buildPath

cd $buildPath

installPath="$buildPath/install"
mkdir -p $installPath/bin $installPath/lib

timeStart=$(date +%s)
test="OFF"
if [ $# == 1 ]
then
    if [ $1 == "ON" ]
    then
        test="ON"
    fi
fi

cmake $srcPath -DINSTALL_PATH=$installPath -DCMAKE_INSTALL_PREFIX=$installPath

thread=$(grep -c ^processor /proc/cpuinfo)
make -j $thread install #VERBOSE=1
if [ $? != 0 ]
then
    exit 1
fi

install(){
    libFile="$1"
    tlib="$installPath/lib"
    if [ -d "$tlib" ]
    then
        if [ -e "$libFile" ]
        then
            echo "Installing: $libFile into $tlib"
            cp $libFile $tlib
        else
            echo "[WARN] File not exists: $libFile"
        fi
    fi
}

install $htslibPath/lib/libhts.so.2
install $libdeflatePath/lib/libdeflate.so.0


# Copy reference data
if [ -e "$installPath/bin/anno" ]
then
    echo "Reference data exists"
else
    cp -R $srcPath/anno $installPath/bin/
    echo "Installing: reference data"
fi

# Copy README
cp -R $srcPath/README.md $installPath/bin/
cp -R $srcPath/CHANGELOG.md $installPath/bin/


# compile py code
$py -m py_compile $srcPath/d2c/plot.py
cp $srcPath/d2c/__pycache__/plot*.pyc $installPath/bin/plot.pyc
cp $srcPath/d2c/rank.R $installPath/bin/

cd $srcPath
cp -R $installPath $srcPath
#ln -fs $installPath .
if [ -e "$buildPath" ]
then
    rm -rf $buildPath
fi

timeEnd=$(date +%s)
secs=$(($timeEnd - $timeStart))
mins=$(($secs/60))
secs=$(($secs%60))
echo "Cost: $mins Mins $secs Secs"
