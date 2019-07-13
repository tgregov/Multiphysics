:: This file is intended to generate debug and release Code::Blocks project under windows.
:: Please replace the windows compiler by the latest 64 bits compiler from here
:: -> https://sourceforge.net/projects/mingw-w64/files/latest/download
:: How to create a new compile profile in Code::Blocks
:: -> https://medium.com/@yzhong.cs/code-blocks-compile-64-bit-under-windows-with-mingw-w64-79101f5bbc02
:: Do not forget to change the CB project compiler when you launch the .cbp !

@ECHO OFF

set GMSHSDK=C:\Program Files (x86)\CodeBlocks\gmsh-4.4.0-Windows64-sdk
set EIGENSDK=C:\Program Files (x86)\CodeBlocks\eigen-eigen-323c052e1731

:: where is gmsh.exe and gmsh-**.dll ? (HINT: copy gmsh-**.dll to the bin folder)
set PATH=%GMSHSDK%\bin;%GMSHSDK%\lib;%PATH%
:: where is gmsh.h ? (rename gmsh.h_cwrap => gmsh.h)
set INCLUDE=%EIGENSDK%;%GMSHSDK%\include;%INCLUDE%
:: where is gmsh.lib ?
set LIB=%GMSHSDK%\lib;%LIB%
:: where is gmsh.py ? (required only if you want to use the python API)
set PYTHONPATH=%GMSHSDK%\lib;%PYTHONPATH%

cd ../../

rd /s /q build
mkdir build
cd build

mkdir Release
cd Release
cmake ../../ -DCMAKE_BUILD_TYPE=Release  -G "CodeBlocks - MinGW Makefiles"
copy "%GMSHSDK%\bin\gmsh-4.4.dll" "%cd%\bin"
xcopy /E /I "../../Geometry" "%cd%\bin\Geometry"
xcopy /E /I "../../Params" "%cd%\bin\Params"

cd ../

mkdir Debug
cd Debug
cmake ../../ -DCMAKE_BUILD_TYPE=Debug  -G "CodeBlocks - MinGW Makefiles"
copy "%GMSHSDK%\bin\gmsh-4.4.dll" "%cd%\bin"
xcopy /E /I "../../Geometry" "%cd%\bin\Geometry"
xcopy /E /I "../../Params" "%cd%\bin\Params"

cd ../../
PAUSE
