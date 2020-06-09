########PATC User Interace####
#Run a long with server.R in an interactive version or deploy shiny server
###Amhed Vargas
###amhed.velazquez@kaust.edu.sa
###UI

##Required packages
#install.packages("shiny")
#install.packages("shinythemes")
#install.packages("ggvis")
#install.packages("ggplot2")
#install.packages("DT")
#install.packages("shinyWidgets")
#install.packages("Cairo")

#Load libraries
library(shiny)
library(shinythemes)
library(ggvis)
library(ggplot2)
library(DT)
library(shinyWidgets)
library(Cairo)


####Busy indicator function based from shinybusy function
###Its only a wrapper for a javascript timeout/setinterval function
busyIndicator <- function(text = "Server is worming up... Please wait", wait=2000) {
  shiny::tagList(
    shiny::div(class="coso",id="loadmessage",text,img(src="elegans3.gif"))
      ,shiny::tags$script(sprintf(
      " setInterval(function(){
         if ($('html').hasClass('shiny-busy')) {
          setTimeout(function() {
            if ($('html').hasClass('shiny-busy')) {
              $('div.coso').show()
            }
          }, %d)          
        } else {
          $('div.coso').hide()
        }
      },100)
      ",wait)
        )
  ) 
}


#########################################
####Start here: Definition User interface
shinyUI(
    fluidPage(
      ###HTML header that describes the object "loadmessage" used in shiny busy indicator.
      tags$head(tags$style(type="text/css", "
             #loadmessage {
               position: fixed;
               top: 60px;
               left: 0px;
               width: 100%;
               padding: 5px 0px 5px 0px;
               text-align: center;
               font-weight: bold;
               font-size: 100%;
               color: #000000;
               background-color: #D3D3D3;
               z-index: 105;
             }
          ")),
    ##Costum extra css styles: Remove sliders background and change color of navbar brand (title of page)  
    tags$style(type = 'text/css', 
               ".js-irs-none .irs-single, .js-irs-none .irs-bar-edge, .js-irs-none .irs-bar {
                          background: transparent;
                          border-top-color: transparent;
                          border-bottom-color: transparent;
                          border-left-color: transparent;
                          border-right-color: transparent}
               .navbar-default .navbar-brand:hover {color: #ffffff;}
               "),
    #Call busy indicator on top of all tabs.
    busyIndicator(),
    #script to renew igv browser everytime a tab is changed. This solves issue regarding Nan positioning when browser is not in main tab.
    tags$script("
                     Shiny.addCustomMessageHandler('igvstat-change', function(panel) {
                     igv.visibilityChange()
                     });
                                 "),
    ######Main tab pages coded in a navbarPage######
    navbarPage(
      title=actionLink("link_to_tabpanel_title", HTML("<b>PATC app</b>")),
      windowTitle="PATCs app",
      ###Theme of shiny
        theme = shinytheme("flatly"),
        id = "panels",
        
        ###Analysis panel
        tabPanel("Analysis",
                 mainPanel(h1("Analyze your DNA sequence"),
                           tabsetPanel(
                               tabPanel("Sequence Input",
                                        ##Sequence Input
                                        textAreaInput("text", label = HTML("<h3>Paste any DNA sequence</h3><h5>(limited to 50 Kb)</h5>"), value = "", cols= 100, rows=5, width = "600px"),
                                        fluidRow(
                                          column(width = 8,
                                                 plotOutput("plot2", height = 300,
                                                            brush = brushOpts(
                                                                id = "plot2_brush",
                                                                direction=c("x"),
                                                                resetOnNew = TRUE
                                                            )
                                                 )
                                          ),
                                          column(width = 4,
                                                 plotOutput("plot3", height = 300)
                                          ),
                                          htmlOutput("computedfft")
                           ),
                           fluidRow(
                             column(4,
                                    h3("Options:"),
                                    #custom slider sliderInput("thres", "Phasing treshold", min = 0, max = 95, value = 55, step = 5),
                                    tags$div(class = "js-irs-none", sliderInput("thres", HTML("Phasing threshold<sup>1</sup>"), value = 60, min = 5, max = 95, step = 5)),
                                    checkboxInput("flabal",HTML("<b>Balanced algorithm<sup>2</sup></b>"), value=TRUE),
                                    checkboxInput("flag_hist",HTML("<b>Analyze PATC periodicity<sup>3</sup></b>"), value=TRUE),
                                    actionButton("action2", label = "Calculate PATCs")
                             )
                             ),
                                        hr(),
                                        HTML(
                                            "<font size=\"2\">",
                                            "<b><sup>1</sup>: Threshold</b> =<br> 60 <b>-></b> ~ 1% Phasing in random DNA<br>95 <b>-></b> ~ .001% Phasing in random DNA<br>
                                       Fire <i>et al</i>. (2006)<br>"),
                                        HTML("<p align=\"justify\"><b><sup>2</sup>:</b> The balanced algorithm (Frøkjær-Jensen <i>et al</i>., 2016) substracts
                                       \"off-helical\" A<sub>n</sub>/T<sub>n</sub> signals in an effort to reduce false positive PATC signals in repeat regions and 
                                       A/T rich genomes.</p>"),
                                        HTML("<p align=\"justify\"><b><sup>3</sup>:</b> Display histogram of frequencies between A/T<sub>4</sub> motifs and highest 
                                        signals periodicty (for the first 20bp) calculated via a Fourier transform.</p>",
                                             "</font>"),
                                        hr(),
                                        DT::dataTableOutput("table2"),
                                        htmlOutput("inc")
                               ),
                               tabPanel("Fasta Input",
                                        ###Fasta Input
                                        fileInput("file", label = h3("Fasta(s) file input"), accept = c(".fasta",".fa",".fna"), multiple = FALSE),
                                        fluidRow(
                                            column(6,
                                        h3("Options:"),
                                        #custom slide
                                        tags$div(class = "js-irs-none", sliderInput("thres", HTML("Phasing threshold<sup>1</sup>"), value = 60, min = 5, max = 95, step = 5)),
                                        checkboxInput("flabal2",HTML("<b>Balanced algorithm<sup>2</sup></b>"), value=TRUE),
                                        checkboxInput("flag_hist2",HTML("<b>Analyze PATC periodicity<sup>3</sup></b>"), value=TRUE),
                                        ##Download demo trhough handler
                                        downloadButton("downloadData", "Demo File"),
                                        actionButton("action1", label = "Calculate PATCs")
                                            )),
                                        fluidRow(
                                          column(width = 8,
                                                 plotOutput("plot4", height = 300,
                                                            brush = brushOpts(
                                                              id = "plot2_brush",
                                                              direction=c("x"),
                                                              resetOnNew = TRUE
                                                            )
                                                 )
                                          ),
                                          column(width = 4,
                                                 plotOutput("plot5", height = 300)
                                          )),
                                          htmlOutput("computedfft2"),
                                        hr(),
                                        HTML(
                                            "<font size=\"2\">",
                                            "<b><sup>1</sup>: Threshold</b> =<br> 60 <b>-></b> ~ 1% Phasing in random DNA<br>95 <b>-></b> ~ .001% Phasing in random DNA<br>
                                       Fire <i>et al</i>. (2006)<br>"),
                                        HTML("<p align=\"justify\"><b><sup>2</sup>:</b> The balanced algorithm (Frøkjær-Jensen <i>et al</i>., 2016) substracts
                                       \"off-helical\" A<sub>n</sub>/T<sub>n</sub> signals in an effort to reduce false positive PATC signals in repeat regions and 
                                       A/T rich genomes.</p>"),
                                        HTML("<p align=\"justify\"><b><sup>3</sup>:</b> Display histogram of frequencies between A/T<sub>4</sub> motifs and highest 
                                        signals periodicty (for the first 20bp) calculated via a Fourier transform.</p>",
                                             "</font>"),
                                        HTML("<b>Please note that the results are reported in bins of 1 megabase in length. For fasta files of length greater than 1MB, consecutive windows will be reported.</b>"),
                                        hr(),
                                        DT::dataTableOutput("table1") 
                               )
                           )
                 )),
                ####Interactive panel
        tabPanel("Interactive",
                 mainPanel(
                     HTML("<h1>PATCs across the <i>C. elegans</i> genome</h1>"), ###introduction
                                  div(style="text-align:justify",
                                  
                                  h4("")
                                  ),
                                  #####interactive plot
                                  fluidRow(
                                      column(4, ###sliders
                                             selectInput("chromo","Chromosome",c("All","I","II","III", "IV", "V", "X", "MtDNA")),
                                             sliderInput("cpos", "Relative position in chromosome (%)",0, 100, c(0, 100), step = 1),
                                             sliderInput("papos", "Display top genes ranked by PATC content (%)",0, 100, c(0, 100), step = 1),
                                             checkboxInput("flascal","Log scale Y axis", value=FALSE),
                                             textAreaInput("genetext", label = "Gene search", value = "", placeholder= "Wormbase ID or common name"),
                                             actionButton("actionsearch", label = "Go!"),
                                             hr(),
                                             htmlOutput("geneid")
                                             ),
                                      column(8, ggvisOutput("plot1"))
                                  )
                         )
                 ),
        ####Genome browser tab
        tabPanel("Genome Browser",
                     HTML("<h1>PATC signal across the <i>C. elegans</i> genome<br></h1> "),
                     h4(""),
                     tags$head(tags$script(src ="igv.min_all-features.js")),
                     includeHTML("www/igv.html")
        ),
        
        
		tabPanel("Background",
		         mainPanel(
		             HTML("<h1>Periodic A<sub>n</sub>/T<sub>n</sub> Clusters (PATCs)</h1><br>"),
		             HTML("<p align=\"justify\">Fire <i>et al.</i> (2006) first identified an unusual non-coding DNA structure that was associated with genes expressed in the germline of
                      <i>C. elegans</i>. The unusual DNA structure consists of periodic clusters of A or T nucleotides spaced by 10 bp and extends over relatively long distances (hundreds to thousands of base pairs).
                      For this reason, the structure was named Periodic A<sub>n</sub>/T<sub>n</sub> Clusters (PATCs) (Fire <i>et al</i>., 2006).
                      PATCs are strongly associated with genes that are expressed in the germline and, when present, are strongly enriched in non-coding regions  (5\', intronic, and 3\').
                      At a genome-scale, PATCs are strongly enriched on autosome arms, which are also enriched for repressive chromatin modifications and transposable elements (Liu <i>et al</i>., 2010). 
                      However, on a local scale, PATCs are anti-correlated with repressive histone marks (Gu and Fire, 2010) which suggests that PATCs may confer resistance to repressive chromatin.
                      </p><br>"),
		             HTML("<p align=\"justify\">
                      Functionally, PATCs can partially prevent (trans)gene silencing in the germline of <i>C. elegans</i>.
                      We showed that single-copy fluorescent transgenes containing PATC-rich sequences in introns were resistant to positional silencing in repressive
                      genomic domains (\"arms\") and stochastic silencing in euchromatic domains (\"centers\") (Frøkjær-Jensen <i>et al</i>., 2016).
                      The effect appeared to be specific to germline expression; we could not detect enhanced somatic expression.  Other laboratories have shown that transgenes engineered to contain PATCs are resistant to silencing in different contexts:
                      Zhang <i>et al</i>. (2018) showed that PATC-rich transgenes were resistant to small RNA-mediated silencing via the piRNA pathway.
                      Fielmich <i>et al</i>. (2018) demonstrated that PATCs alleviated the silencing of endogenous genes tagged by CRISPR. 
                      </p><br>"),
		             HTML("<p align=\"justify\">
                      At present, we have very little understanding of how PATCs influence (trans)gene silencing but a working model is that PATCs constrain DNA 
                      and nucleosome interactions to resist the assembly of higher-order heterochromatic structures (Fire <i>et al</i>., 2006). 
                      <br><br>Here we have developed a set of tools based on the original PATC algorithm (Fire <i>et al</i>., 2006) 
                      to enhance the study of PATCs in <i>C. elegans</i> and other organisms:<br>An "),
		             actionLink("link_to_tabpanel_interactive", "interactive"),
		             HTML(" overview of pre-calculated PATC values for all <i>C. elegans</i> protein-coding genes.
                      <br>A "),
		             actionLink("link_to_tabpanel_genome_browser", "genome browser"),
		             HTML(" with a continuous track of PATCs values.
                      <br>An online "),
		             actionLink("link_to_tabpanel_analysis", "analysis"),
		             HTML(" web service that calculates the PATC value of arbitrary DNA sequences. 
                      </p><br>"),
		             HTML("<p align=\"justify\">
                      <b>References</b>
                      <font size=\"1\">
                      <br>Fire, A., Alcazar, R., and Tan, F. (2006). Unusual DNA structures associated with germline genetic activity in <i>Caenorhabditis elegans</i>. <i>Genetics 173</i>, 1259–1273.
                      <br>Gu, S.G., and Fire, A. (2010). Partitioning the <i>C. elegans</i> genome by nucleosome modification, occupancy, and positioning. <i>Chromosoma 119</i>, 73–87.
                      <br>Liu, T., Rechtsteiner, A., Egelhofer, T.A., Vielle, A., Latorre, I., Cheung, M.-S., Ercan, S., Ikegami, K., Jensen, M., Kolasinska-Zwierz, P., <i>et al</i>. (2010). Broad chromosomal domains of histone modification patterns in <i>C. elegans</i>. <i>Genome Res</i>.
                      <br>Frøkjær-Jensen, C., Jain, N., Hansen, L., Davis, M.W., Li, Y., Zhao, D., Rebora, K., Millet, J.R.M., Liu, X., Kim, S.K., <i>et al</i>. (2016). An Abundant Class of Non-coding DNA Can Prevent Stochastic Gene Silencing in the <i>C. elegans</i> Germline. <i>Cell 166</i>, 343–357.
                      <br>Fielmich, L.-E., Schmidt, R., Dickinson, D.J., Goldstein, B., Akhmanova, A., and van den Heuvel, S. (2018). Optogenetic dissection of mitotic spindle positioning in vivo. <i>ELife 7</i>, e38198.
                      <br>Zhang, D., Tu, S., Stubna, M., Wu, W.-S., Huang, W.-C., Weng, Z., and Lee, H.-C. (2018). The piRNA targeting rules and the resistance to piRNA silencing in endogenous genes. <i>Science 359</i>, 587–592.
                      </font>
                      </p><br>")
		         )),
		
		
		
		
		###Software panel: please do update links if any of them are down.
        tabPanel("Software",
                 mainPanel(
                 h1("Download the PATC algorithm"),
                 HTML("If required, use \"worm\" as password."),
                 h3("Software"),
                 HTML("<h4><a href=\"https://exrcsdrive.kaust.edu.sa/exrcsdrive/index.php/s/COn5idUy2PgUnSp\">PATC 2006 for mac</a></h4>"),
                 HTML("<h4><a href=\"https://exrcsdrive.kaust.edu.sa/exrcsdrive/index.php/s/2Hq0PNjmcAVajQZ\">PATC 2006 for linux</a></h4>"),
                 HTML("<h4><a href=\"https://exrcsdrive.kaust.edu.sa/exrcsdrive/index.php/s/SCfVrGqhjKRyXi0\">PATC 2016 for mac (balanced option)</a></h4>"),
                 HTML("<h4><a href=\"https://exrcsdrive.kaust.edu.sa/exrcsdrive/index.php/s/dfoNulv57Luu1kA\">PATC 2016 for linux (balanced option)</a></h4>"),
                 HTML("<h4><a href=\"https://exrcsdrive.kaust.edu.sa/exrcsdrive/index.php/s/dFODDwTEPrKdEJq\">PATC Source files</a></h4>"),
                 h3("Data Files"),
                 HTML("<h4><a href=\"https://exrcsdrive.kaust.edu.sa/exrcsdrive/index.php/s/Nc2CTOAk96xO0zi\">PATC algorithm documentation (Fire <i>et al</i>., 2006)</a></h4>"),
                 HTML("<h4><a href=\"https://exrcsdrive.kaust.edu.sa/exrcsdrive/index.php/s/8ufIeDwuVgulry6\">Pre-calculated PATC Genomic trace for 
                      <i>C. elegans</i> ce11/WS245 (Frøkjær-Jensen <i>et al</i>., 2016)</a></h4>")
                 )
                 ),
        ###About Panel
        tabPanel("About",
                 mainPanel(
                 h3("The app"),
                 HTML("<p align=\"justify\">This website is generated via custom modified css/html code running in R via the shiny library.
                 Different public available shiny apps from the <a href=\"https://shiny.rstudio.com/gallery/\">R-studio shiny gallery</a> were source of inspiration in the development of this app. 
                 <br>The genome browser is produced via a custom modified script based on <a href=\"https://github.com/igvteam/igv.js/\">igv.js</a>.
                 <br>All the templates, libraries, and programs used to produce this site are under the MIT and GNU licenses.</p>"),
                 h3("The PATC algorithm"),
                 HTML("<p align=\"justify\">
                      Andrew Fire <i>et al.</i> (2006) developed the original PATC algorithm. This server uses a modified \"balanced\" PATC algorithm (Frøkjær-Jensen <i>et al</i>., 2016). This app acts as a front-end for the balanced alogirhtm which runs in the background.
                      </p>"),
                 h3("The Laboratory of Synthetic Genome Biology"),
                 HTML("<p align=\"justify\">
                 The Laboratory of Synthetic Genome Biology is located in building 2 - level 3 (Ibn Al-Haytham – Above Spine) at King Abdullah University of Science and Technology (KAUST).
                 <br><i>Contact info</i>:<br>Christian-Froekjaer Jensen, Ph.D. 
                 <br>Assistant Professor of Bioscience
                 <br>Laboratory of Synthetic Genome Biology
                 <br>Email: <a href=\"mailto:cfjensen@kaust.edu.sa\">cfjensen@kaust.edu.sa</a>
                 
                      </p>"),
                 h3("The people behind the app"),
                 HTML("<p align=\"justify\">
                      The app was originally conceived by Christian Frøkjær-Jensen and implemented by <a href=\"https://www.researchgate.net/profile/Amhed_Vargas_Velazquez\">Amhed Missael Vargas Velazquez</a>.
                      Furthermore, the PATC app is on the web thanks to Amazon Web Services (AWS) and the original efforts of the Linux and Advanced Platforms team in KAUST. 
                      </p>")
        )
        )
    ),
    hr(),
    HTML("<a href=\"https://syngenbio.kaust.edu.sa\">Syntetic genome biology laboratory</a> @ <a href=\"https://kaust.edu.sa/en\"> KAUST</a><br>"),
    HTML("<a href=\"http://www.wormbuilder.org/\">Wormbuilder</a><br>"),
    HTML("<a href=\"mailto:amhed.velazquez@kaust.edu.sa\">Contact us!</a>")
    
)
)



