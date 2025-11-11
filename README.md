# minimesh

## Related software
**Student Name:** Ziyu Sun  
**Student Number:** 98696719  

---

## Guide (On Windows)

```bash
# Unzip
unzip mesh.zip
unzip third-party.zip

# Create build directory
mkdir build 
cd build

# Configure with CMake (Debug mode)
cmake .. -G "Visual Studio 17 2022" -A x64 -DCMAKE_BUILD_TYPE=Debug -DFREEGLUT_DIR="%cd%/../third-party/freeglut/bin-win64-msvc2017/debug" -DGLUI_DIR="%cd%/../third-party/glui/bin-win64-msvc2017/debug" -DEIGEN3_DIR="%cd%/../third-party/eigen"

# Build and copy DLLs
cmake --build . --config debug 
cmake --build . --config Debug --target COPY_DLLS

# select loop / butterfly subdivision by command line parameter
./bin/minimeshgui.exe OBJ_PATH 

# Configure with CMake (Release mode)
cmake .. -G "Visual Studio 17 2022" -A x64 -DCMAKE_BUILD_TYPE=Release -DFREEGLUT_DIR="%cd%/../third-party/freeglut/bin-win64-msvc2017/release" -DGLUI_DIR="%cd%/../third-party/glui/bin-win64-msvc2017/release" -DEIGEN3_DIR="%cd%/../third-party/eigen"

# Build and copy DLLs
cmake --build . --config Release
cmake --build . --config Release --target COPY_DLLS 
