# Basic layout for a c++ program compiled with cmake  
The first step is to make the run time aware of the created libraries  
This only needs to be done once per cmake project, but absolutely required! This code, as well as any code linked to this code will not work, if the LD_LIBRARY_PATH is not set correctly.  
```
echo -e "\n\nexport LD_LIBRARY_PATH=$(pwd -P)/install/lib:\${LD_LIBRARY_PATH}" >> $HOME/.bashrc
source $HOME/.bashrc
```

## How to compile  
```
cd build && cmake .. -DCMAKE_INSTALL_PREFIX=../install  
--with cmake3 
make -j $(($(nproc)-1))  
make -j $(($(nproc)-1)) install  
```

## How to execute  
```
cd install/bin  
./Program
```

## Troubleshooting  
If your program does not compile, e.g. if you added something and it is not picked up correctly, try first  
```
make clean && make -j $(($(nproc)-1)) install
```

If this does not work try to run cmake again  
```
cd build  
rm -rf *
cmake .. -DCMAKE_INSTALL_TARGET=../install 
make -j $(($(nproc)-1))  
make -j $(($(nproc)-1)) install  
```

## How to add new structures  
If you want to add a new subdirectory with code, just execute ./InitCMAKESubDirectory.sh <Your_Folder_Name>  
If you want to have additional cpp, h files, just add them in your favorite subdirectory.  
Have a look into the CMakeLists.txt file for more explanations.  
Follow all steps under `How to compile` again, whenever you add new files.  

## How to link against other libraries  
1) Locate the folder where the shared object (.so) files of your external library lie. This is usually in some kind of lib folder and make sure, that this folder is known to runtime, i.e. look if it is in your LD_LIBRARY_PATH. If not, repeat the first steps described in this readme for your external library  
2) Locate the folder where the header files (.h) of your external library lie. This is usually in some kind of include folder  
3) In your CMakeLists.txt, add the following line before you add subdirectories! There is also some space reserved for it.  
```
include_directories(<Path_To_Header_File_Directory>)
link_directories(<Path_To_.so_File_Directory>)
```
4) Given you external library is called libTest.so, then add Test as target to target_link_libraries() in the CMakeLists.txt  
4a) Note, that the lib and the .so is dropped!
5) Congratulations! Now you can use whatever code is defined in the external library, everywhere in your code.  

## How to create multiple Programs with different mains?  
1) Add a new .cpp file in your top level directory  
2) Look into the CMakeLists.txt file and locate the code lines  
```
add_executable(Program main.cpp)
target_link_libraries(Program ${ROOT_LIBRARIES} CodeDir1)
install(TARGETS Program DESTINATION bin)
```
3) Copy and past them and substitute Program for your new program name and extend the included libraries if needed  
4) Repeat the `How to compile` steps

