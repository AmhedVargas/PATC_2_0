# PATC shiny app

This project represents the first stage of the integration of bioinformatic tools for the enhancement of *C. elegans* transgenesis. Here, we implement a user interface for the deployment of an algorithm written in Pascal as well as visualization tools.

## General description

Periodic A/T Clusters, or PATCs, are common DNA structures seen in the genome of *C. elegans* particularly in the non-coding regions of germline genes. 

Recently, PATCs have been used to allow the germline expression of foreign genes (trans-genes) in *C. elegans* demonstrating this way the importance and usability of PATCs as tools for streamlined transgenesis.

PATCs are difficult to detect by human eye; fortunately, there is a plethora of computational algorithms that can be used to detect these structures. In this project, we make use of an old (but surprisingly fast) program that has been written in the programming language Pascal. In order to facilitate the use of this program, we have deployed a shinny app that transforms user inputs into the appropiate format for the program and viceversa. 

This shiny app was orginially deployed [within KAUST](https://wormbuilder.kaust.edu.sa/) thanks to the aid of the IT Linux and Advanced Platforms team. Nowadays, the app is hosted and powered via the Amazon Web Services (AWS) and can be accesed at [https://wormbuilder.dev/PATC/](https://wormbuilder.dev/PATC/).


## Structure of the application

The core of the program recides in a R shinny app that is divided into a server.R and ui.R code. While the ui handles the user queries, the server computes them and executes system calls for the pascal code.

Basic processing diagram:

User -> ui.R -> server.R -> NewPATC052115Balanced (compiled Pascal program) -> server.R -> ui.R -> User

The communication between the server.R code and the Pascal program is mediated by a single shell script acting in temporary directory created by the shiny app. Upon disconection, the folder is removed automatically.  In addition, there is extra files that helps the user to browse pre-computed PATC values.

Additionally, the PATC app show a genome browser that relies on a custom modified version of the [igv javascript](https://github.com/igvteam/igv.js/) 2.5.5 .

**Directory structure**

**root**
*   server.R
*   ui.R
*   **www** `Folder containing web sources accesible to the shiny app`
*   *   elegans3.gif `Small gif that is displayed while the server is busy`
*   *   igv.html `HTML script that calls and defines tracks for the genome browser. Please modify if different tracks are needed`
*   *   igv.min_all-features.js `Javascript modified for adhoc display of labelings. For specific igv pop-overs, the igv.html file was modified instead`
*  **PATC**     `Directory for computation and system calls. It **must** have rw privileges for the working user`
*  *    PATCsGenesChristianData.tsv `File that contains gene info as well as their PATC values. Christian Froek-jaer Jensen kindly provided the data. See his 2016 paper for more info.`
*  *    patc.sh    `Shell script to handle Pascal code. Receives input from console within R via a "system" command`
*  *    user folder    `Directory that contains temporary folders for each connected user; it won't be there till the code runs sucessfully at least once`
*  *    NewPATC052115Balanced    `Unix x86 compiled Pascal script. In case of running in other operative system, the script has to be recompiled (see below)
*  **source** `Only needed if running not in an Unix based system.`
*   *   NewPATC052115Balanced.pas   `Pascal script to be compiled`

## Dependencies
**R**

While the program has been succesfully tested in R 3.6.1 (2019-07-05, and later versions) on a x86_64-pc-linux-gnu platform, the application should work fine in any other version as long as the following libraries are succesfully installed:
*   install.packages("shiny")
*   install.packages("shinythemes")
*   install.packages("ggvis")
*   install.packages("ggplot2")
*   install.packages("DT")
*   install.packages("shinyWidgets")
*   install.packages("Cairo")
*   install.packages("Biostrings")

**Shell**

`Programs already installed in most common unix distributions`
*  awk
*  perl


`Any program to compile Pascal code such as fpc`

*   [fpc](https://www.freepascal.org/)

## Installation procedure

Please first make sure to have installed R along with all its dependencies in your working enviroment. Similarly, make sure to have installed a compiler of pascal code. In the following lines, *fpc* has been used but any other compiler should do just fine.

1. Copy or clone the contents of this git in to a new directory

`git clone https://gitlab.kaust.edu.sa/velazqam/wormbuilder-patc.git`

`cd wormbuilder-patc/`

2. Compile Pascal code and move to main PATC directory

`cd source/`

`fpc NewPATC052115Balanced.pas`

`mv NewPATC052115Balanced ../PATC/`

`cd ..`

2. (Optional) Change the rw privilegies of PATC directory if needed

`chmod 777 PATC/`

3. Deploy shiny server

The final app should be looking something like [this](https://wormbuilder.dev/PATC/)


## Troubleshout

Please feel free to [e-mail me](mailto:amhed.velazquez@kaust.edu.sa) for any question, doubt or error in the app.

## Acknowledgements

Different public available shiny apps from the [R-studio shiny gallery](https://shiny.rstudio.com/gallery/) were source of inspiration in the development of this app. Similarly, the deployment of the app used a similar schema seen in Jacques Serizay's [RegAtlas](https://github.com/js2264/RegAtlas); the deployment and construction of the app greatly benefited from its dcoumentaion. 
All the templates, libraries, and programs used to produce this site are under the MIT and GNU licenses.

## The PATC algorithm
Andrew Fire et al. (2006) developed the original PATC algorithm. This server uses a modified "balanced" PATC algorithm (Frøkjær-Jensen et al., 2016). This app acts only as a front-end for the PATC algorihtm which runs in the background. For a proper description of how it works, see Fire et al. 2006 paper.


