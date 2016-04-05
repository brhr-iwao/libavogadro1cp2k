## libavogadro1cp2k
CP2K input generator and output visualizer for Avogadro 1.
It also generates CHARMM-Style Parameter (pot) / Protein Structure File (psf) .

###  Description
A generator of [CP2K](http://cp2k.org/) input, [CHARMM](http://www.charmm.org)-style Force Field Parameter (pot) and Protein Structure File (psf)  and a [CP2K](http://cp2k.org/)  output visualizer for [Avogadro](http://avogadro.cc/) 1.1 (or later).

###  Installation
#### Using with precompiled binary

An easy way to use libavogadro1cp2k is to use with precompiled [Avogadro](http://avogadro.cc/) for Windows.

1. Download Avogadro Installer for Windows (Avogadro-1.1.1-win32.exe) from [Sourceforge] (http://sourceforge.net/projects/avogadro/) and install the Avogadro binary on your Windows PC.

2. Copy [released cp2kextension.dll] (https://github.com/brhr-iwao/libavogadro1cp2k/releases) to (Avogadro install directory, maybe C:\Program Files\Avogadro)\lib\Avogadro\1_1\ extensions directory.

3. In order to use the pot generator, it is necessary to put gaff.dat, which is a plain text file and found in http://archive.ambermd.org/201111/att-0689/gaff.dat, in (Avogadro install directory, maybe C:\Program Files\Avogadro)/bin directory.

4. Launch Avogadro. If you find 'Psf/Pot...' and 'CP2K' entries in 'Extensions' menu, the installation was successful.

#### Compiling with Visual Studio
Here I described only a way to intall on Windows from sources since I worked only on Windows...

1. Build and install Avogadro 1.1.1 (the latest version) from [github master](https://github.com/cryos/avogadro) source codes
by following the instruction [Compiling on Windows with MSVC2008](http://avogadro.cc/wiki/Compiling_on_Windows_with_MSVC_2008)
or [Compiling on Windows with MSVC 2010](http://avogadro.cc/wiki/Compiling_on_Windows_with_MSVC_2010) accroding as your environment.

2. Copy the cp2k folder and all contents to (avogadro root)/libavogadro/src/extensions directory. Add a sentence "add_subdirectory(cp2k)" onto the end of (avogadro root)/libavogadro/src/extensions/CMakeLists.txt.

3. Launch Visual Studio 20XX command prompt, change the current working diretory to (avogadro root)/scripts and execute the command "nmake" or "nmake cp2kextension".

4. If nmake completes successfully, cp2kextension.dll is created in (avogadro root)/scripts/lib (or (avogadro root)/scripts/lib/Release).
 Put this dll in (Avogadro install dir, maybe C:/Program Files/Avogadro)/lib/avogadro/1_1/extensions/.
 
5. If you get "Psf/Pot..." and "CP2K"entries in Avogadro "Extensions" menu, the installation was successful.

### Usage
Load the molecule you want to work with on Avogadro and then select "Psf/Pot..." or "CP2K > Generate Input..." in "Extensions" menu.
This opens a dialog where correspondent parameters can be configured.
To Visualize a CP2K calculation output, select "CP2K > Analyze Output..." in "Extensions" menu and open a "Analyze CP2K output "dialog. 
Select proper file type via the "Data Type" combo box and load a CP2K output file with the "Load a file..." button.


**Note:** This software is in early developement stage and should be used for
experimentation only. There is absolutely no guarantee of correctness of
generated input files and there are probably many bugs.

###  Contact
[Aoyama Iwao](https://github.com/brhr-iwao)
