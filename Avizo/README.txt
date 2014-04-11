Compilation on Windows
1) Install prebuilt binary version of VTK 5.8 into Dependencies folder.
   http://pointclouds.org/downloads/windows.html
2) Open Avizo XPand Development Wizard. Set Local Avizo Directory to <SVN>/Avizo folder
3) Create Build System -> Select Single Package -> CCPi and Press OK
4) Open file <SVN>/Avizo/src/CCPi/CCPi.VC100.proj this will prompt you to enter the solution space.
set the solution space to <SVN>/Avizo/Avizo.sln
5) Currently on release version can only be built.