## libavogadro1cp2k
cp2k input and psf/pot generator for avogadro

### CP2K input and CHARMM-Style Parameter (pot) / Protein Structure File (psf) generator for Avogadro 1.x

###  Description
A generator of [CP2K](http://cp2k.org/) input, CHARMM-Style Force Field Parameter (pot) and Protein Structure File (psf)
for [Avogadro](http://avogadro.cc/) 1.1 (or later).

###  Installation
Here I described only the way to intall on Windows since I worked on only Windows...

Build and install Avogadro 1.1.1 (the latest version) from [github master](https://github.com/cryos/avogadro) source codes
according to the instruction 
[Compiling on Windows with MSVC 2008](http://avogadro.cc/wiki/Compiling_on_Windows_with_MSVC_2008)
or
[Compiling on Windows with MSVC 2010](http://avogadro.cc/wiki/Compiling_on_Windows_with_MSVC_2010).

Copy the cp2k folder and all contents to (avogadro root)/libavogadro/src/extensions.
Add a sentence "add_subdirectory(cp2k)" at the end of (avogadro root)/libavogadro/src/extensions/CMakeLists.txt.

Launch Visual Studio 20XX command prompt, change the current working diretory to (avogadro root)/scripts,
and run the command "nmake cp2kinputextension".

If nmake completes successfully, cp2kinputextension.dll is created in (avogadro root)/scripts/lib.
Put this dll in (Avogadro install dir, maybe C:/Program Files/Avogadro)/lib/avogadro/1_1/extensions/.

If you get "Psf/Pot..." and "CP2K Input..."entry in Avogadro "Extensions" menu, the installation was successful.

###  Usage

Load the molecule you want to work with on Avogadro and then select "Psf/Pot..." or "CP2K Input..." in "Extensions" menu.
This opens a dialog where correspondent parameters can be configured.

**Note:** This software is in early developement stage and should be used for
experimentation only. There is absolutely no guarantee of correctness of
generated input files; they might produce unrealistic results, crash randomly,
or scare your cat. Or any combination of those.
(This note was shamelessly copied from a leading project [CP2K input generator for 
Avogadro 2](https://github.com/infuniri/avogadrolibs-cp2k) )

###  Contact
[Aoyama Iwao](https://github.com/brhr-iwao)