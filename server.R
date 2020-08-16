########PATC server####
###Amhed Vargas
###amhed.velazquez@kaust.edu.sa
###Server

##Required packages
#install.packages("shiny")
#install.packages("shinythemes")
#install.packages("ggvis")
#install.packages("ggplot2")
#install.packages("DT")
#install.packages("shinyWidgets")
#install.packages("Cairo")
#install.packages("Biostrings")

#Load libraries
library(shiny)
library(shinythemes)
library(ggvis)
library(ggplot2)
library(DT)
library(shinyWidgets)
library(Cairo)
library(Biostrings)

  shinyServer(function(input, output, session) {
    ####Initial sequence#######
    ##Retrieve unique ID for the session
    session_id <- session$token
    
    ##Create temporary folder for unique user
    system(paste("mkdir -p PATC/users/",session_id,sep=""))
    
    options(shiny.maxRequestSize=300*1024^2) ###Limits max upload to 300 Megabytes. #Copy and paste 50Kb limit comes from default html maxlength parameter (Default value is 524288)


    #####Read data for functions##########################
        ###Data explorer
        ChrisPAT=read.table("PATC/PATCsGenesChristianData.tsv",sep="\t", header= TRUE)
        rownames(ChrisPAT)=as.character(ChrisPAT[,2])
        ##Assign Yes, No, unknown insted of TRUE, FALSE, NA
        ChrisPAT$GermlineExpression[which(!(ChrisPAT$GermlineExpression))] <- "No"
        ChrisPAT$GermlineExpression[which(ChrisPAT$GermlineExpression == TRUE)] <- "Yes"
        ChrisPAT$GermlineExpression[which(is.na(ChrisPAT$GermlineExpression))] <- "unknown"
    ###################################################
    
    #####Functions#####################################
    ###Exit function#
    ###On exit, force the remove of directory
    ##Attention: Interactive sessions create problems because the listening of the server stops in main directory and not in sub directory
    session$onSessionEnded(function(){
      system(paste("rm -rf PATC/users/",session_id,sep=""))
    }
    )
    ################################################################


    ###Control panels functions##########################################################################################
    ##FUnctions needed to generate links between panels
    observeEvent(input$link_to_tabpanel_interactive, {
    newvalue <- "Interactive"
    updateTabsetPanel(session, "panels", newvalue)
  })

    observeEvent(input$link_to_tabpanel_genome_browser, {
    newvalue <- "Genome Browser"
    updateTabsetPanel(session, "panels", newvalue)
  })

    observeEvent(input$link_to_tabpanel_analysis, {
    newvalue <- "Analysis"
    updateTabsetPanel(session, "panels", newvalue)
  })
    
    observeEvent(input$link_to_tabpanel_title, {
      newvalue <- "Analysis"
      updateTabsetPanel(session, "panels", newvalue)
    })
    #########################################################################################################


    ###Gene explorer functions#######################################################################################
    
    ##Hover action in graphic
    gene_tooltip <- function(x) {
      if (is.null(x)) return(NULL)
      if (is.null(x$Wormbase.ID)) return(NULL)
      # Pick out the gene with this ID
      gene <- ChrisPAT[as.character(x$Wormbase.ID),]
      
      paste0("<b>Gene: ", gene$Gene.name, "</b><br>Wormbase ID: ",
             gene$Wormbase.ID, "<br>",
             "Total PATC value: ", gene$Total.value.of.PATC.algorithm
      )
    }
    
    ##Click action - inputs to reactive variable geneval that is used in a function to generate table
    geneweb_tooltip <- function(x) {
      if (is.null(x)) return(NULL)
      if (is.null(x$Wormbase.ID)) return(NULL)
      # Pick out the gene with this ID
      geneval$id <- as.character(x$Wormbase.ID)
      return(NULL)
    }
    
    geneval <- reactiveValues(id=c(""))
    
    ##Function that returns the gene info from main table and formats it into html
    geneinfo = function(x){
      if (is.null(x)) return(NULL)
      gidi=""
      if(x %in% as.character(ChrisPAT$Wormbase.ID)){gidi = x}
      if(x %in% as.character(ChrisPAT$Gene.name)){gidi = rownames(ChrisPAT[which(as.character(ChrisPAT$Gene.name) == x),])}
      if(gidi ==""){return(paste("Gene not found! try with Wormbase ID"))}
      
      gene <- ChrisPAT[as.character(gidi),]
      
      paste0("<b>Gene: ", gene$Gene.name, "</b><br>Wormbase ID: ",
             "<a href=\"https://www.wormbase.org/species/c_elegans/gene/", gene$Wormbase.ID,"\">", gene$Wormbase.ID, "</a><br>",
             "Transcript analyzed: ", gene$Transcript.name, "<br>",
             "Chromosome: ", gene$Chromosome, "<br>",
             "Relative position: ", gene$Relative.loc.chromosome, "<br>",
             "Gene size: ", gene$Bin.size, "<br>",
             "Number of phased bases: ", gene$Bases.with.PATC.55, "<br>",
             "Phased bases in gene: ", as.integer(gene$Frequency.of.PATCs.55 * 100), "%<br>",
             "Total PATC score: ", gene$Total.value.of.PATC.algorithm, "<br>",
             "PATC density: ", gene$PATC.density, "<br>",
             "Germline expression: ", gene$GermlineExpression, "<sup>1</sup><br>",
             "RPKM Oocyte expression<sup>2</sup>: ", gene$stoekius_oocyte_rpkm, "<br>",
             "RPKM Sperm expression<sup>3</sup>: ", gene$spermatogenic_gonad_fem.3_RPKM_Ortiz_2014,"<br>","<br>",
             "<p align=\"justify\">",
             "<font size=\"2\">",
             "<sup>1</sup>: Oocyte Reads Per Kilobase of transcript per Million mapped reads (RPKM) > 2<br>",
             "<sup>2</sup>: Stoeckius, <i>et al.</i> (2019). Large-scale sorting of <i>C. elegans</i> embryos reveals the dynamics of small RNA expression.
                                  <i>Nat. Methods 6(10)</i>: 745-751<br>",
             "<sup>3</sup>: Ortiz, <i>et al.</i> (2014). A New Dataset of Spermatogenic <i>vs.</i> Oogenic Transcriptomes in the Nematode <i>Caenorhabditis elegans</i>.
             <i>G3: GENES, GENOMES, GENETICS 4(9)</i>: 1765-1772<br>",
             "</font>",
             "</p>"
      )
    }
    
    ##Observers for action button search
    observeEvent(input$actionsearch, {
      mygene =input$genetext
      output$geneid <- renderUI({
        HTML(
        geneinfo(mygene)
        )
        })
    }, ignoreInit = T)
    
    #Observer reactive to changes in geneval to create table
    observeEvent(geneval$id, {
      mygene = geneval$id
      output$geneid <- renderUI({
        HTML(
        geneinfo(mygene)
        )
        })
    }, ignoreInit = T)
    
    ##Reactive subtable used for plotting
    ChrisPAps <- reactive({
      chromos = input$chromo
      xone=input$cpos[1]/100
      xtwo=input$cpos[2]/100
      yone=input$papos[1]/100
      ytwo=input$papos[2]/100
      flagscalelog=input$flascal
      crisp = ChrisPAT
      
      if(chromos != "All"){ crisp = subset(crisp, Chromosome == chromos)}
      if(input$flascal){crisp[which(crisp[,9] > .001),9] = log(crisp[which(crisp[,9] > .001),9])}
      
      crisp = subset(crisp, Relative.loc.chromosome >= xone & Relative.loc.chromosome <= xtwo)
      maxpat=quantile(crisp$PATC.density, c(ytwo), na.rm=TRUE)
      minpat=quantile(crisp$PATC.density, c(yone), na.rm=TRUE)
      crisp = subset(crisp, PATC.density >= (minpat) & PATC.density <= (maxpat))
      
      ##Highlight data if genee names present
      if(input$genetext %in% rownames(crisp)){crisp[which(input$genetext == rownames(crisp)),7] = "No" }
      if(input$genetext %in% as.character(crisp$Gene.name)){crisp[which(input$genetext == as.character(crisp$Gene.name)),7] = "highlit" }
      
      ##Return table to plot
      crisp
    })
    
    # A reactive expression with the ggvis plot that takes in changes in subtable
    vis <- reactive({
      # Lables for axes
      xvar_name <- "Normalized position in Chromosome"
      yvar_name <- "PATC density"
      
      ChrisPAps %>%
        ggvis(x = ~Relative.loc.chromosome, y = ~PATC.density) %>%
        layer_points(size := 40, size.hover := 200,
                     fillOpacity := 0.2, fillOpacity.hover := 0.5,
                     fill = ~GermlineExpression, 
                     stroke = ~GermlineExpression,
                     shape = ~GermlineExpression,
                     key := ~Wormbase.ID) %>%
        add_tooltip(gene_tooltip, "hover") %>%
        add_tooltip(geneweb_tooltip, "click") %>%
        add_axis("x", title = xvar_name) %>%
        add_axis("y", title = yvar_name) %>%
        add_legend(c("fill","stroke", "shape"), title = "Germline (oogenic) expression", values = c("Yes", "No", "unknown")) %>%
        ###Nominal scales that modifies the shapes of the points
        scale_nominal("shape", domain = c("Yes", "No", "unknown","highlit"), range = c("circle", "diamond", "triangle-up","cross")) %>%
        scale_nominal("fill", domain = c("Yes", "No","unknown","highlit"),
                      range = c("lightblue", "lightgrey","gold","tomato")) %>%
        scale_nominal("stroke", domain = c("Yes", "No", "unknown","highlit"),
                      range = c("blue", "grey", "orange","red")) %>%
        set_options(width = 1000, height = 600)
    })
    ##########################

    ##Send to output
    vis %>% bind_shiny("plot1")
    
    ######################################################################################################

    ###PATC analysis tab#################################################################################
    ####
        ######Makett################### table on top of DNA sequence
    makett = function(){
        tt=read.table(paste("PATC/users/",session_id,"/tab.file", sep=""),sep="\t",header=TRUE)
        colnames(tt)=c("Sequence ID", "Number of bases","Bases in phase", "Phasing frequency", "Total PATC value", "Phasing density")
        ##Change name
        if(as.character(tt[1,1])=="file.fasta"){tt[1,1]="Sequence input"}
        ##Rround stuff
        tt[,4]=round(tt[,4], 2)
        tt[,5]=round(tt[,5], 2)
        tt[,6]=round(tt[,6], 2)
        return(tt)
    }

    ##########higfft function
    ##Function to get highest bp periodicity
        higfft=function(x){
          ##If A/T clusters appears  less than5 times, do not perform analysis
          if(length(x) < 5){return(NULL)}else{
            ##Calculate spectrum based on fft
            spect = spectrum(x, plot = FALSE)
            #Create period to analyze
            period = seq(2, 20, .5)
            ##Convert spectral signal to frequency by ordering indexes based on 1/period
            consig = spect$spec[unlist(lapply(c(1/(period)), function (freq) {
              idx = which.min(abs(freq - spect$freq))
            }))]
            ##Return period with highest signal
            return(period[which.max(consig)])
          }}
        
        ############
   ##################
        ###Fasta input#####################
    ###Observer of submit button in fasta input
    ##The whole structure of this function is to A) retrieve user inputs, B)System call to PATC algorithm, C) get back output and place them into ui
    observeEvent(input$action1, { 
      strfile=input$file
      path=strfile$datapath
      if(input$flabal2){ ###Blanced flag
        balflag="B"}else{balflag="0"
          }
  
      ##Move temporary file into usage path
      system(paste("mv ", path, " PATC/users/",session_id,"/file.fasta", sep=""))
      
      ####New format of sh
      ##patc.sh directory_id treshold balanced_flag type_flag
      system(paste("sh PATC/patc.sh",session_id,input$thres, balflag,0)) 
      
      output$table1 = DT::renderDataTable(DT::datatable({
       makett()
      }, options = list(dom = 't')))
      
      ###Add visualization 4 & 5 if checkbox
      if(input$flag_hist2){
        ###Read file to analyze
        inseq=readLines(paste("PATC/users/",session_id,"/file.fasta", sep=""))
        if(length(grep(">",inseq)) > 0) {inseq=inseq[-c(grep(">",inseq))]}
        
        inseq=paste(inseq,collapse="")
        
        bseq=DNAString(gsub(pattern="\n",replacement="",x=toupper(inseq)))
        
        dists=dist(start(matchPattern("WWWW",bseq, fixed=FALSE)))
        vals=table(dists)
        freqAT= as.integer(vals)
        xas=as.integer(names(vals))
        freqPAT = data.frame(BPseparation=xas[-c(1:5)], Ocurrences=freqAT[-c(1:5)])
        
        valu=higfft(vals)
        if(!is.null(valu)){output$computedfft2 <- renderUI({ HTML(paste0("Highest periodicity signal found at ","<b>",valu,"bp </b>",sep=""))})}
        
        
        output$plot4 <- renderPlot({
          ggplot(freqPAT, aes(BPseparation, Ocurrences)) +
            ggtitle("A/T spacing (click and drag to see details)") +
            geom_line() +
            coord_cartesian(xlim=c(0,max(freqPAT$BPseparation)), ylim=c(0,max(freqPAT$Ocurrences)),expand = FALSE)
        })
        
        output$plot5 <- renderPlot({
          ggplot(freqPAT, aes(BPseparation, Ocurrences)) +
            ggtitle("Detailed view") +
            geom_line() +
            coord_cartesian(xlim = ranges2$x, ylim = ranges2$y, expand = FALSE)
        })
      }else{
        output$plot4 <- renderPlot({ggplot()})
        output$plot5 <- renderPlot({ggplot()})
        output$computedfft2 <- renderUI({ HTML(paste0(""))})
      }
    })
    
    ###############Observer DemoSeq button
    #observe({
    #  input$DemoSeq
    #  Demova=readLines("PATC/Demo.fasta")
    #  updateTextAreaInput(session, "text", value = paste(Demova[1],"\n",Demova[2],"\n",sep=""))
    #})
    #######################
    
    #################Observer Reset button
    #observe({
    #  input$ResetSeq
    #  updateTextAreaInput(session, "text", value = c(""))
    #  })
    ############

    ##################Download handler
    output$downloadData <- downloadHandler(
      filename <- function() {
        paste("Smu-1_unspliced", "fasta", sep=".")
      },
      
      content <- function(file) {
        file.copy("PATC/Demo2.fasta", file)
      },
    )
    
    output$DownSeqOut <- downloadHandler(
      filename <- function() {
        paste("Genebuild", "fasta", sep=".")
      },
      
      content <- function(file) {
        file.copy(paste("PATC/users/",session_id,"/SeqOpop.fasta", sep=""), file)
      },
    )
    
    ######################################################

    ###Observer of submit button in DNA sequence input#################################
    observeEvent(input$action2, {   
      #write(">DNA_input",paste("PATC/users/",session_id,"/file.fasta", sep=""))
      #write(input$text,paste("PATC/users/",session_id,"/file.fasta", sep=""),append= TRUE)
      
      write(input$text,paste("PATC/users/",session_id,"/file.fasta", sep=""))
      
      if(input$flabal){ ###Blanced flag
        balflag="B"}else{balflag="0"
        }
      
      ####New format of sh
      ##patc.sh directory_id treshold balanced_flag type_flag
      system(paste("sh PATC/patc.sh",session_id,input$thres, balflag,1)) 
      
      output$table2 = DT::renderDataTable(DT::datatable({
       makett()
      }, options = list(dom = 't')))
      
      ##External function to include the HTML output
      #It has to be now relocated within temporary directory
      getPage<-function(dir) {
        return(includeHTML(paste(paste("PATC/users/",dir,"/DNA-patc.html", sep=""))))
      }
      output$inc<-renderUI({getPage(session_id)})
      
      ###Add visualization 4 & 5 if checkbox
      if(input$flag_hist){
        inseq=readLines(paste("PATC/users/",session_id,"/file.fasta", sep=""))
        if(length(grep(">",inseq)) > 0) {inseq=inseq[-c(grep(">",inseq))]}
           
           inseq=paste(inseq,collapse="")
           
           bseq=DNAString(gsub(pattern="\n",replacement="",x=toupper(inseq)))
           
           dists=dist(start(matchPattern("WWWW",bseq, fixed=FALSE)))
           vals=table(dists)
           freqAT= as.integer(vals)
           xas=as.integer(names(vals))
           freqPAT = data.frame(BPseparation=xas[-c(1:5)], Ocurrences=freqAT[-c(1:5)])
           
           valu=higfft(vals)
           if(!is.null(valu)){output$computedfft <- renderUI({ HTML(paste0("Highest periodicity signal found at ","<b>",valu,"bp </b>",sep=""))})}
           
        
        output$plot2 <- renderPlot({
          ggplot(freqPAT, aes(BPseparation, Ocurrences)) +
            ggtitle("A/T spacing (click and drag to see details)") +
            geom_line() +
            coord_cartesian(xlim=c(0,max(freqPAT$BPseparation)), ylim=c(0,max(freqPAT$Ocurrences)),expand = FALSE)
        })
        
        output$plot3 <- renderPlot({
          ggplot(freqPAT, aes(BPseparation, Ocurrences)) +
            ggtitle("Detailed view") +
            geom_line() +
            coord_cartesian(xlim = ranges2$x, ylim = ranges2$y, expand = FALSE)
        })
      }else{
        output$plot2 <- renderPlot({ggplot()})
        output$plot3 <- renderPlot({ggplot()})
        output$computedfft <- renderUI({ HTML(paste0(""))})
      }
      
    })
    ######################################################


    #######Brush for ggplot###########
    ranges2 <- reactiveValues(x = NULL, y = NULL)
    
    # When a double-click happens, check if there's a brush on the plot.
    # If so, zoom to the brush bounds; if not, reset the zoom.
    observe({
      brush <- input$plot2_brush
      if (!is.null(brush)) {
        ranges2$x <- c(brush$xmin, brush$xmax)
        if((brush$xmax - brush$xmin) > 1000){ranges2$x <- c(brush$xmin, brush$xmin + 1000)}
        ranges2$y <- c(brush$ymin, brush$ymax)
        
      } else {
        ranges2$x <- c(0, 100)
        ranges2$y <- NULL
      }
    })
    #############
    ###Initialize webpage##########################
    init=function(){
      load(file="PATC/init/Initial.RData")
      updateTextAreaInput(session, "text", value = paste(Demova[1],"\n",Demova[2],"\n",sep=""))
      
      output$table2 = DT::renderDataTable(DT::datatable({
        tt
      }, options = list(dom = 't')))
      
      ##External function to include the HTML output
      #It has to be now relocated within temporary directory
      output$inc<-renderUI({
        includeHTML(paste(paste("PATC/init/DNA-patc.html", sep="")))
        })
      
      ###Add visualization 4 & 5 if checkbox
      if(input$flag_hist){        
        if(!is.null(valu)){output$computedfft <- renderUI({ HTML(paste0("Highest periodicity signal found at ","<b>",valu,"bp </b>",sep=""))})}
        
        output$plot2 <- renderPlot({
          ggplot(freqPAT, aes(BPseparation, Ocurrences)) +
            ggtitle("A/T spacing (click and drag to see details)") +
            geom_line() +
            coord_cartesian(xlim=c(0,max(freqPAT$BPseparation)), ylim=c(0,max(freqPAT$Ocurrences)),expand = FALSE)
        })
        ranges2$x = c(0,100)
        output$plot3 <- renderPlot({
          ggplot(freqPAT, aes(BPseparation, Ocurrences)) +
            ggtitle("Detailed view") +
            geom_line() +
            coord_cartesian(xlim = ranges2$x, ylim = ranges2$y, expand = FALSE)
        })
      }else{
        output$plot2 <- renderPlot({ggplot()})
        output$plot3 <- renderPlot({ggplot()})
        output$computedfft <- renderUI({ HTML(paste0(""))})
      }
      
    }
    ###########################################################################################################
    
    
    ################Starting function#############
    
    #On starting, put Smu-1##########
    observeEvent(session$clientData,{init()})
    
    # prints actuall tab
    observeEvent(input$panels,{
        sendM(input$panels);
    })
    
    ##function to send check tab status and change visibility of browser
    sendM = function(x){
      session$sendCustomMessage("igvstat-change", x)
    }
    
  })
