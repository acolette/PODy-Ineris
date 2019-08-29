PODySPEC
================

## Overview

These scripts are designed to evaluate the impact of ozone on
agriculture by computing PODySPEC (*Phytotoxic Ozone Dose*). The
methodology is described in the
[manual](https://icpvegetation.ceh.ac.uk/chapter-3-mapping-critical-levels-vegetation)
by the **IPC Vegetation**.

The scripts provided allow to compute PODy for various species of crops
(Bread Wheat, Potato, Tomato), trees (Spruce, Temperate Oak, Beech) and
for Grass. They have been used to study ***Long-term evolution of the
impacts of ozone air pollution on agricultural yields in Europe*** in
the [Eionet
report](https://acm.eionet.europa.eu/reports/EIONET_Rep_ETCACM_2018_15_O3impactTrends).

The following section describes the scripts, the input variables needed
& the documentation.

## Description

### Tool & documentation

In the folder, we provide some documentations explaining the methodology
used in the scripts: the manual and the annexes and the eionet report
(application to Europe of the script for wheat). We also provide a small
input files (*pody\_input\_EDT\_CHIF\_2010\_light.nc*) regrouping all
the variables needed as well as some parameters for the ozone deposition
derived from the chemistry and transport model
[CHIMERE](http://www.lmd.polytechnique.fr/chimere/) (*DEPO\_PARS*) and
species related parameter derived from the mapping manual
(*SPECIES\_PAR.csv*). These files are regrouped in the folder *data* &
*input* respectively. Finally two scripts are provided: *phenology.R* &
*pody.R*.

First *phenology.R* must be launched to preprocess the phenology, then
*pody.R* must be launched to compute PODy for one specie. Both script
are based on hourly input variables described in the last section.

### Input variables

The following table regroups the variables used as input in the code
with the dimention and unit of all the
variables.

| Variables | Names                                    | Dimension | Timestep | Units                              |
| --------- | ---------------------------------------- | --------- | -------- | ---------------------------------- |
| alti      | Altitude                                 | 2D        |          | m                                  |
| FC        | Water Content at Field Capacity          | 2D        |          | cm<sup>3</sup>.cm<sup>-3</sup>     |
| PWP       | Water Content at Permanent Wilting Point | 2D        |          | cm<sup>3</sup>.cm<sup>-3</sup>     |
| lon       | Longitude                                | 2D        |          | \[-180,180\]                       |
| lat       | Latitude                                 | 2D        |          | \[-180,180\]                       |
| sreh      | Relative humidity                        | 3D        | hourly   | m<sup>3</sup>.m<sup>-3</sup>       |
| swvl      | Soil moisture                            | 3D        | hourly   | m<sup>3</sup>.m<sup>-3</sup>       |
| surf temp | Surface temperature                      | 3D        | hourly   | °C                                 |
| swrd      | Short Wave Radiation                     | 3D        | hourly   | W.m<sup>-2</sup>                   |
| ppfd      | Photosynthetic Photon Flux Density       | 3D        | hourly   | µmol.m<sup>-2</sup>.s<sup>-1</sup> |
| ustar     | Friction velocity                        | 3D        | hourly   | m.s<sup>-1</sup>                   |
| obuk      | Monin-Obhukov Length                     | 3D        | hourly   | m                                  |
| airm      | Air molecule                             | 3D        | hourly   | molec.cm<sup>-3</sup>              |
| O3        | Ozone                                    | 3D        | hourly   | ppb                                |

### Output variables

The following table regroups the output variables (lon, lat & PODy) of
the light
output.

| Variables | Names                 | Dimension | Timestep | Units                  |
| --------- | --------------------- | --------- | -------- | ---------------------- |
| lon       | Longitude             | 2D        |          | \[-180,180\]           |
| lat       | Latitude              | 2D        |          | \[-180,180\]           |
| PODy      | Phytotoxic Ozone Dose | 2D        | yearly   | mmol.m<sup>2</sup>.PLA |

If you set the option *heavy output* to *TRUE* then a **hourly netcdf
file** of the following variables is written (could be a heavy output
depending on your grid). It is usefull to understand how evolve the PODy
by looking at the limitation functions, the ozone (downscalled depending
on the height of the first model
layer)…

| Variables | Names                             | Dimension | Timestep | Units                  |
| --------- | --------------------------------- | --------- | -------- | ---------------------- |
| lon       | Longitude                         | 1D        |          | \[-180,180\]           |
| lat       | Latitude                          | 1D        |          | \[-180,180\]           |
| PODy      | Phytotoxic Ozone Dose             | 3D        | hourly   | mmol.m<sup>2</sup>.PLA |
| O3        | Ozone canopy                      | 3D        | hourly   | ppb                    |
| accper    | Accumulation period               | 3D        | hourly   | \[0,1\]                |
| flight    | Phytotoxic Ozone Dose             | 3D        | hourly   | \[0,1\]                |
| fO3       | Ozone limitation function         | 3D        | hourly   | \[0,1\]                |
| fphen     | Phenology limitation function     | 3D        | hourly   | \[0,1\]                |
| fsmi      | Soil moisture limitation function | 3D        | hourly   | \[0,1\]                |
| ftemp     | Temperature limitation function   | 3D        | hourly   | \[0,1\]                |
| fvpd      | Air humidity limitation function  | 3D        | hourly   | \[0,1\]                |

## Future development

Here we list some future developments:

1.  Add biogeographical region & spatialize some values
2.  Add PODyIAM (generic PODy)
3.  Add sumVPD is missing to correct the PODy
4.  Enhance the code (simplify, OOB programing) & documentation
5.  Find all the parametrization (all values for each types of trees…).
    Here we have used Temperate Oak Spain Italy for example.
