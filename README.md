Community Contribution and AI Policy
====================================

ExoPlaSim is currently maintained by a single developer in their spare time. As such the pace of development can be slow, and as an open-source project community contributions are greatly welcomed. Contributions developed with the assistance of LLMs/AI are also welcome, as ExoPlaSim has users outside of the scientific community for whom physics and Fortran may not be fluent skills, and AI tools may allow a larger subset of the ExoPlaSim user community to contribute to both bugfixes and feature development. 

However, the use of such tools does not change expectations for contributor responsibility, engagement, or code quality. As many peer-reviewed publications now require detailed AI contribution statements, there are additional expectations which apply specifically to AI-assisted contributions, for the sake of transparency and to expedite the PR approval process.

Please note that contributions for which the AI was the main contributor in all 3 of discovery, ideation, and implementation may be rejected. All contributions must be fully-understood by a human contributor; therefore potential contributors are encouraged to use AI as a tool to accelerate work (like a very advanced debugger or autocomplete), as opposed to allowing broad, self-directed work by the AI. Science is a human endeavour and ExoPlaSim is first and foremost a scientific tool. All contributions must therefore have ultimately human provenance.

The following requirements apply to all contributions.

1. **Pull request descriptions:** Although code contributions written by LLM/AI are acceptable (after detailed review and understanding by a contributor, see 2), all pull request descriptions must be written without LLM/AI assistance. The aim of this policy is to help ensure that the descriptions are concise and interpretable for maintainers, and that contributors are able to explain and justify their changes.
2. **Human ownership and accountability:**
    1. Contributors are responsible for all submitted content, regardless of whether LLMs/AI tools were used.
    2. Contributors must be able to explain all changes during review.
    3. Contributors are responsible for ensuring all submitted code complies with ExoPlaSim’s licensing.
4. **Disclosure:** If LLMs/AI tools were used for substantive portions of the contribution, contributors must describe their use in the pull request description. These descriptions must be detailed and itemized, with a breakdown of both provenance of ideas and provenance of work, as enumerated below. If a PR contains multiple work items, the following should be described for each work item:
    1. Whether the problem being addressed was discovered and characterized by a human or AI.
    2. Whether the solution was primarily determined by a human or AI.
    3. Whether implementation in the code was led by a human or AI.
    4. Estimated token usage for the task, and which AI/LLM was used for the task. AI-assisted PRs for which the work was done by an unknown AI will be rejected immediately without further review. Similarly, due to the high environmental cost of AI usage and ExoPlaSim's domain as a climate tool, AI-assisted contributions for which there is no estimate of AI inference cost (token usage is an acceptable proxy) will also be rejected immediately without further review.
5. **Authentic engagement:** The pull request process is collaborative and iterative. Contributors are expected to actively engage with reviewer feedback themselves. Copying and pasting replies to / from an LLM/AI tool does not count as engaging with the reviewer; contributors are expected to understand the changes they are proposing, not act as message courier for a third party. As such, PRs opened by agents without human engagement may be closed.
6. **Consistency with conventions and existing style:** This particularly means conveying information (comments, documentation, etc) in a style similar to that already used in ExoPlaSim and the underlying PlaSim, and that is appropriate and to the point. Commentary within the code that is deemed to be excessively verbose or unclear may need to be removed or revised.


This policy is adapted from the [Astropy](<https://github.com/astropy/astropy-project/blob/main/policies/ai-policy.md>) and [COSMIC](<https://github.com/COSMIC-PopSynth/COSMIC>) "Usage of Generative AI" policies.

Version information
==========================

Final maintenance release for 3.x. Includes numerous quality-of-life improvements, modernizations, and bugfixes.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2533357.svg)](https://doi.org/10.5281/zenodo.2533357)

EXOPLASIM
=========

![Two rows of planets, progressing from yellow to blue from top left to bottom right. The top row appears to represent tidally-locked planets, while the bottom row appears to represent Earth-like planets.](mixplanets.png "Two rows of planets, progressing from yellow to blue from top left to bottom right. The top row appears to represent tidally-locked planets, while the bottom row appears to represent Earth-like planets.")
*A range of planets modeled by ExoPlaSim, and postprocessed with SBDART. The top row consists of tidally-locked aquaplanets at T21 orbiting stars ranging from 2500 K to 4000 K, with orbital periods increasing with stellar mass. The bottom row consists of aquaplanets with 24-hour rotation at T42, orbiting stars ranging from 4000 K to 8000 K.*

The PlaSim 3D general climate model, extended for terrestrial planets. This model contains the PlaSim GCM, as well as all necessary modifications to run tidally-locked planets, planets with substantially different surface pressures than Earth, planets orbiting stars with different effective temperatures, super-Earths, and more. This model also includes the ability to compute carbon-silicate weathering, dynamic orography through the glacier module (though only accumulation and ablation/evaporation/melting are included; glacial flow and spreading are not), and storm climatology. Future features will include support for multiple celestial light sources (e.g. for a habitable moon orbiting a Jovian planet, or circumbinary planets), coupling with N-body integrators such as REBOUND, and CO2 condensation.

This model ships with a Python API, described below. Full documentation of the API is available at <http://exoplasim.readthedocs.io>.

Documentation of the original PlaSim model is found in exoplasim/plasim/doc.

ExoPlaSim Python API README
===========================

Created by Adiv Paradise

Copyright 2020, Distributed under the General Public License

This API was written with Python 3 in mind, but should work with Python
2 and outdated versions of NumPy.

Requirements
------------

-   numpy
-   scipy
-   matplotlib (only needed for additional utilities)
-   GNU C (gcc/g++) and Fortran (gfortran) compilers (for Python utilities)
-   (optionally) Other compilers whose use you prefer for the model itself
-   (optionally) MPI libraries for those compilers
-   netCDF4 (optional)
-   h5py (optional)

Compatibility
*************

-   Linux (tested on Ubuntu 18.04, CentOS 6.10): **Yes**
-   Google Colaboratory: Yes (note that OpenMPI support on Colaboratory is limited due to automatic root privileges; look up how to run OpenMPI executables with root permissions and note that this is not recommended)
-   Windows 10: Yes, via Windows Subsystem for Linux
-   Mac OS X: Yes, requires Xcode and developer tools, and [OpenMPI support requires that Fortran-compatible libraries be built.](<https://www.open-mpi.org/faq/?category=osx#not-using-osx-bundled-ompi>) Tested on Mac OS X Catalina and Big Sur (with MacPorts, GCC10, OpenMPI, and Anaconda3), Apple M1 compatibility has not been tested.

Installation
------------

    pip install exoplasim

OR:

    python setup.py install

The first time you import the module and try to create a model
after either installing or updating, ExoPlaSim will run a 
configuration script.

Multiple output formats are supported by the built-in `pyburn`
postprocessor. If you wish to use HDF5 or NetCDF output formats, you
will need the netCDF4-python and h5py libraries, respectively. You
can ensure these are included at install-time by specifying them:

    pip install exoplasim[netCDF4]
    
OR:

    pip install exoplasim[HDF5]
    
OR:

    pip install exoplasim[netCDF4,HDF5]

You may also configure and compile the model manually if you wish to not
use the Python API, by entering the exoplasim/ directory and running
first configure.sh, then compile.sh (compilation flags are shown by
running `./compile.sh -h`).

PlaSim Documentation
--------------------

Original PlaSim documentation is available in the exoplasim/docs/
folder.

Usage
-----

To use the ExoPlaSim Python API, you must import the module, create a
Model or one of its subclasses, call its configure method and/or modify
method, and then run it.

An IPython notebook is included with ExoPlaSim; which demonstrates
basic usage. It can be found in the ExoPlaSim installation directory,
or [downloaded directly here.](https://raw.githubusercontent.com/alphaparrot/ExoPlaSim/master/exoplasim/exoplasim_tutorial.ipynb)

Basic example::

    import exoplasim as exo
    mymodel = exo.Model(workdir="mymodel_testrun",modelname="mymodel",resolution="T21",layers=10,ncpus=8)
    mymodel.configure()
    mymodel.exportcfg()
    mymodel.run(years=100,crashifbroken=True)
    mymodel.finalize("mymodel_output")

In this example, we initialize a model that will run in the directory
"mymodel\_testrun", and has the name "mymodel", which will be used to
label output and error logs. The model has T21 resolution, or 32x64, 10
layers, and will run on 8 CPUs. By default, the compiler will use 8-byte
precision. 4-byte may run slightly faster, but possibly at the cost of
reduced stability. If there are machine-specific optimization flags you
would like to use when compiling, you may specify them as a string to
the optimization argument, e.g. `optimization='mavx'`. ExoPlaSim will
check to see if an appropriate executable has already been created, and
if not (or if flags indicating special compiler behavior such as
debug=True or an optimization flag are set) it will compile one. We then
configure the model with all the default parameter choices, which means
we will get a model of Earth. We then export the model configurations to
a `.cfg` file (named automatically after the model), which will allow
the model configuration to be recreated exactly by other users. We run
the model for 100 years, with error-handling enabled. Finally, we tell
the model to clean up after itself. It will take the most recent output
files and rename them after the model name we chose, and delete all the
intermediate output and configuration files.
