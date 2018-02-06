IF NOT DEFINED CIL_VERSION (
ECHO CIL_VERSION Not Defined.
exit 1
)

mkdir "%SRC_DIR%\build"
ROBOCOPY /E "%RECIPE_DIR%\..\..\Core" "%SRC_DIR%\build"
cd "%SRC_DIR%\build"

echo "we should be in %SRC_DIR%\build"

cmake -G "NMake Makefiles" -DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Release -DLIBRARY_LIB="%CONDA_PREFIX%\lib" -DLIBRARY_INC="%CONDA_PREFIX%" -DCMAKE_INSTALL_PREFIX="%PREFIX%\Library" -DINSTALL_LIB_DIR="%PREFIX%\Library\lib" "%SRC_DIR%\build"

:: Build C library
nmake install
if errorlevel 1 exit 1

:: Install step