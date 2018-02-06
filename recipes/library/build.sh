echo $CIL_VERSION

mkdir ${SRC_DIR}/build
cp -rv ${RECIPE_DIR}/../../Core/ ${SRC_DIR}/build
mkdir ${SRC_DIR}/build/build
cd ${SRC_DIR}/build/build
cmake -G "Unix Makefiles" -DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Release -DINSTALL_LIB_DIR="${CONDA_PREFIX}/lib" -DLIBRARY_LIB="${CONDA_PREFIX}/lib" -DLIBRARY_INC="${CONDA_PREFIX}" -DCMAKE_INSTALL_PREFIX="${PREFIX}" ../Core

make -j2 VERBOSE=1
make install
