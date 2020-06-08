########PATC-v server####
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
#install.packages

#Load libraries
library(shiny)
library(shinythemes)
library(ggvis)
library(ggplot2)
library(DT)
library(shinyWidgets)
library(Biostrings)
library(Cairo)
library(stringdist)

  shinyServer(function(input, output, session) {
    ####Initial sequence#######
    ##Retrieve unique ID for the session
    session_id <- session$token
    
    ##Create temporary folder for unique user
    system(paste("mkdir -p PATC/users/",session_id,sep=""))
    
    options(shiny.maxRequestSize=300*1024^2) ###Limits max upload to 300 Megabytes. #Copy and paste 50Kb limit comes from default html maxlength parameter (Default value is 524288)


    #####Data##########################
        ###Data explorer
    ###Format input table
        ChrisPAT=read.table("PATC/PATCsGenesChristianData.tsv",sep="\t", header= TRUE)
        rownames(ChrisPAT)=as.character(ChrisPAT[,2])
        ##Assign Yes, No, unknown insted of TRUE, FALSE, NA
        ChrisPAT$GermlineExpression[which(!(ChrisPAT$GermlineExpression))] <- "No"
        ChrisPAT$GermlineExpression[which(ChrisPAT$GermlineExpression == TRUE)] <- "Yes"
        ChrisPAT$GermlineExpression[which(is.na(ChrisPAT$GermlineExpression))] <- "unknown"
    ###################################################
        
        ##Codon Adapt
        CAIS=read.table("PATC/Ultimate_aminos.txt",sep="\t",header=T)
        rownames(CAIS)= toupper(as.character(CAIS$Codon))
        codons=unique(CAIS$Amino)
        
        AAtoCodF=list()
        
        for(i in 1:length(codons)){
          AAtoCodF=append(AAtoCodF, list(CAIS[which(CAIS$Amino == codons[i]),]))
          names(AAtoCodF)[i]=as.character(codons[i])
        }
        
        IntronSeqs=read.table("PATC/Introns.csv",sep=",",header=F,row.names=1)
        
        ##piRNA
        pi2UGs=read.table("PATC/Ultimate_seed-remainder_iterations_2GUs.txt",sep="\t",header=F)
        
        Piseqs=list()
        for(i in 1:nrow(pi2UGs)){
          Piseqs=append(Piseqs,strsplit(as.character(pi2UGs[i,2]),";"))
          names(Piseqs)[i]=as.character(pi2UGs[i,1])
        }
        
        Pies=readLines("PATC/HengPies.txt")
        
    ####Functions#######
        ##Codon
        sampcod=function(aa,list,cai){
          newcod=sample((list[[aa]])[,6],1,prob=(list[[aa]])[,cai])
          return(toupper(as.character(newcod)))
        }
        
        repcds=function(x,tabibi,list,cai){
          x=toupper(x)
          if(!is.character(x)){return(c())}
          vecseq=unlist(strsplit(x,""))
          if((length(vecseq) %% 3) != 0){return(c())}
          nnseq=c()
          for(i in seq(1,length(vecseq),by=3)){
            nncod=sampcod(as.character(tabibi[paste(c(vecseq[i:(i+2)]),sep="",collapse=""),"Amino"]),list,cai)
            nnseq=append(nnseq,nncod)
          }
          return(paste(nnseq,sep="",collapse=""))
        }
        
        sampnewcod=function(aa,oldcodon,list,cai){
          oldcodon=toupper(oldcodon)
          if(nrow(list[[aa]]) < 2 ){return(oldcodon)}
          oldcodon =which(as.character(rownames(list[[aa]])) == oldcodon)
          newcod=sample((list[[aa]])[-c(oldcodon),6],1,prob=(list[[aa]])[-c(oldcodon),cai])
          return(toupper(as.character(newcod)))
        }
        
        repnewcds=function(x,tabibi,list,cai){
          x=toupper(x)
          if(!is.character(x)){return(c())}
          vecseq=unlist(strsplit(x,""))
          if((length(vecseq) %% 3) != 0){return(c())}
          nnseq=c()
          for(i in seq(1,length(vecseq),by=3)){
            nncod=sampnewcod(as.character(tabibi[paste(c(vecseq[i:(i+2)]),sep="",collapse=""),"Amino"]),paste(c(vecseq[i:(i+2)]),sep="",collapse=""),list,cai)
            nnseq=append(nnseq,nncod)
          }
          return(paste(nnseq,sep="",collapse=""))
        }
        
        modbyposiz=function(x,starts,tabibi,list,cai){
          x=toupper(x)
          if(!is.character(x)){return(c())}
          vecseq=unlist(strsplit(x,""))
          if((length(vecseq) %% 3) != 0){return(c())}
          if(length(starts)>0){starts=unique(starts - (starts %% 3) +1 )}
          for(pos in c(starts)){
            nncod=sampnewcod(as.character(tabibi[paste(c(vecseq[pos:(pos+2)]),sep="",collapse=""),"Amino"]),paste(c(vecseq[pos:(pos+2)]),sep="",collapse=""),list,cai)
            vecseq[pos:(pos+2)]=unlist(strsplit(nncod,""))
          }
          return(paste(vecseq,sep="",collapse=""))
        }
        ##Pis
        countpies=function(x,y){
          if(!is.character(x)){return(c())}
          vecseq=unlist(strsplit(x,""))
          if(length(vecseq) != 21){return(c())}
          return((stringdist(paste(vecseq[1:14],collapse=""),unlist(y[paste(vecseq[15:20],collapse="")]), method="h")))
        }
        
        Strcountpies=function(x,y){
          if(!is.character(x)){return(c())}
          vecseq=unlist(strsplit(x,""))
          if(length(vecseq) < 21){return(c())}
          con=0
          for(i in 1:(length(vecseq)-20)){
            con=con+countpies(paste(vecseq[i:(i+20)],collapse=""),y)
          }
          return(con)
        }
        
        countmatpies=function(x,y){
          if(!is.character(x)){return(c())}
          vecseq=unlist(strsplit(x,""))
          if(length(vecseq) != 21){return(c())}
          return(length(which(stringdist(paste(vecseq[1:14],collapse=""),unlist(y[paste(vecseq[15:20],collapse="")]), method="h")<6)))
        }
        
        Strcountmatpies=function(x,y){
          if(!is.character(x)){return(c())}
          vecseq=unlist(strsplit(x,""))
          if(length(vecseq) < 21){return(c())}
          con=0
          for(i in 1:(length(vecseq)-20)){
            con=con+countmatpies(paste(vecseq[i:(i+20)],collapse=""),y)
          }
          return(con)
        }
        
        condepies=function(x,pies,mm){
          if(!is.character(x)){return(c())}
          vecseq=unlist(strsplit(x,""))
          if(length(vecseq) != 20){return(c())}
          return(sum(stringdist(x,pies,method="hamming") <= mm))
        }
        
        Strcondepies=function(x,y,mm){
          if(!is.character(x)){return(c())}
          vecseq=unlist(strsplit(x,""))
          if(length(vecseq) < 20){return(c())}
          con=0
          for(i in 1:(length(vecseq)-19)){
            con=con+condepies(paste(vecseq[i:(i+19)],collapse=""),y,mm)
          }
          return(con)
        }
        
        findpies=function(x,pies,mm){
          if(!is.character(x)){return(c())}
          vecseq=unlist(strsplit(x,""))
          if(length(vecseq) != 20){return(c())}
          idx=c(stringdist(x,pies,method="hamming") <= mm)
          if(sum(idx)>0){return(pies[idx])}else{return()}
        }
        
        Strfindpies=function(x,y,mm){
          if(!is.character(x)){return(c())}
          vecseq=unlist(strsplit(x,""))
          if(length(vecseq) < 20){return(c())}
          con=c()
          for(i in 1:(length(vecseq)-19)){
            con=append(con,findpies(paste(vecseq[i:(i+19)],collapse=""),y,mm))
          }
          return(con)
        }
        ###Analysis
        
        ##Function fft
        higfft=function(x){
          if(length(x) < 10){return(NULL)}else{
            spect = spectrum(x, plot = FALSE)
            
            period = seq(2, 20, .5)
            
            mtm <- spect$spec[unlist(lapply(c(1/(period)), function (freq) {
              idx <- which.min(abs(freq - spect$freq))
            }))]
            
            return(period[which.max(mtm)])
          }}
        
        ############
    ##############Exit function############################
    ###On exit, force the remove of directory
    ##Attention: Interactive sessions create problems because the listening of the server stops in main directory and not in sub directory
    session$onSessionEnded(function(){
      system(paste("rm -rf PATC/users/",session_id,sep=""))
    }
    )
    ################################################################


    ###Control panels##########################################################################################
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


    ###Gene explorer#######################################################################################
    
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
      
      ##Code to test highlighting
      #if(geneval$id %in% rownames(crisp)){crisp = subset(crisp,Wormbase.ID == geneval$id)}
      #if(input$genetext %in% rownames(crisp)){crisp = subset(crisp,Wormbase.ID == input$genetext)}
      #if(input$genetext %in% as.character(crisp$Gene.name)){crisp = subset(crisp,Gene.name == input$genetext)}
      #if(input$genetext %in% rownames(crisp)){crisp[which(input$genetext == rownames(crisp)),9] = 9000 }
      #if(input$genetext %in% as.character(crisp$Gene.name)){crisp[which(input$genetext == as.character(crisp$Gene.name)),9] = 9000 }
      #if(input$genetext %in% rownames(crisp)){crisp$GermlineExpressio[which(input$genetext == rownames(crisp))] <- "unknown" }
      #if(input$genetext %in% as.character(crisp$Gene.name)){crisp$GermlineExpressio[which(input$genetext == as.character(crisp$Gene.name))] <- "unknown" }
      if(input$genetext %in% rownames(crisp)){crisp[which(input$genetext == rownames(crisp)),7] = "No" }
      if(input$genetext %in% as.character(crisp$Gene.name)){crisp[which(input$genetext == as.character(crisp$Gene.name)),7] = "highlit" }
      
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
        #add_legend("stroke", title = "", values = c("", "")) %>%
        #hide_legend("fill") %>%
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
      Demova=readLines("PATC/Demo.fasta")
      updateTextAreaInput(session, "text", value = paste(Demova[1],"\n",Demova[2],"\n",sep=""))
      
      write(paste(Demova[1],"\n",Demova[2],"\n",sep=""),paste("PATC/users/",session_id,"/file.fasta", sep=""))
      
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
    
    #####Trangene generation#####
    observeEvent(input$actionSeq, {  
      
      ErrorFlag=0
      
      CodonAl=as.integer(input$selectCAI)
      FlaPi=input$checkPirna
      FlaIn=input$checkIntron
      FlaBs=input$checkBsaI
      FlaRi=input$checkboxRibo
      FlaSp=input$checkSapI
      FlaEp=input$checkEsp
      
      output$ErrorMessage <- renderText({})
      seqDNA=unlist(strsplit(toupper(input$seqDNA),""))
      
      if((ErrorFlag == 0) & ((length(seqDNA) %% 3) != 0)){ ##Check for errors in size
        output$ErrorMessage <- renderText({
          paste("Error: Sequence is not multiple of three")
        })
        ErrorFlag=1
      }
      
      if((ErrorFlag == 0) & ((length(seqDNA) < 50) & (FlaRi))){ ##Check for errors in size
        output$ErrorMessage <- renderText({
          paste("Error: Sequence is too short to optimize for ribosome binding")
        })
        ErrorFlag=1
      }
      
      if((ErrorFlag == 0) & ((length(seqDNA) < 450) & (FlaIn))){ ##Check for errors in size
        output$ErrorMessage <- renderText({
          paste("Error: Sequence is too short to add ~150bp spaced introns. The option will be deactivated")
        })
        FlaIn=FALSE
        updateCheckboxInput(session, "checkIntron", value = FALSE)
      }
      
      if((ErrorFlag == 0) & (nchar(gsub("A|T|C|G","",toupper(input$seqDNA))) != 0)){ ##Check for strange non ATCG characters
        output$ErrorMessage <- renderText({
          paste("Error: Unrecognized characters in sequence:",input$seqDNA)
        })
        ErrorFlag=1
      }
      
      if((ErrorFlag == 0) & (paste(seqDNA[1:3],sep="",collapse="") != "ATG")){ ##Check for errors in size
        output$ErrorMessage <- renderText({
          paste("Error: Sequence does not start with ATG")
        })
        ErrorFlag=1
      } 
      
      if((ErrorFlag == 0) & ((FlaPi)&(CodonAl == 5))){ ##Check for errors in size
        output$ErrorMessage <- renderText({
          paste("Error: Sequence cannot remove piRNAs and mantain a Codon Adaptation Index equal to 1. piRNAs will not be removed")
        })
        FlaPi=FALSE
        updateCheckboxInput(session, "checkPirna", value = FALSE)
      } 
      
      if(ErrorFlag == 0){##Check stop codon
        seqiqi=paste(seqDNA,sep="",collapse="")
        stpos=c()
        stpos=start(matchPattern("*",translate(DNAString(seqiqi))))
        if(length(stpos) == 0){
          output$ErrorMessage <- renderText({
            paste("Error: Sequence does not have stop codon")
          })
          ErrorFlag=1
        }
        if(length(stpos) > 1){
          output$ErrorMessage <- renderText({
            paste("Error: Sequence have multiple stop codons")
          })
          ErrorFlag=1
        }
        if(length(stpos) == 1){
          if( (stpos) != (length(translate(DNAString(seqiqi)))) ){
          output$ErrorMessage <- renderText({
            paste("Error: Sequence has a stop codon but not at its end.")
          })
          ErrorFlag=1
          }
        }
        }
      
      if(ErrorFlag == 0){ ##Main routine
        
        withProgress(message = 'Generating transgene', style = "notification", detail = "(~4 min per kb for piRNA optimization)", value = 0, {
        ###Internal parameters
        RetrieveTop=100
        
        if(FlaRi){ ###Ribosomal binding optimization
        testSeq=paste(c(seqDNA[1:39]),sep="",collapse="")
        write(testSeq,paste("PATC/users/",session_id,"/gene.fasta", sep=""))
        system(paste("sh bin/rnafold.sh",paste("PATC/users/",session_id,sep=""),RetrieveTop))
        inseqs=readLines(paste("PATC/users/",session_id,"/seqswithoutaaaas.txt", sep=""))
        fwd=grep("GGTCTC",x=inseqs)
        rev=grep("GAGACC",x=inseqs)
        all=unique(c(fwd,rev))
        if((length(all) != 0)&(length(all) != length(inseqs))){inseqs=inseqs[-c(all)]}
        id=order(sapply(inseqs,function(x){Strcondepies(x,Pies,4)}))[1]
        SeqStart=inseqs[id]
        }else{ ###Not optimization
          SeqStart=repcds(paste(c(seqDNA[1:39]),sep="",collapse=""),CAIS,AAtoCodF,CodonAl)
        }
        
        incProgress(3/10)
        
        if((CodonAl == 5)|(!(FlaPi))){ ###Use CAI equal to one
          SeqEnd=repcds(paste(c(seqDNA[40:length(seqDNA)]),sep="",collapse=""),CAIS,AAtoCodF,CodonAl)
          
        }else{ ###Use Codon Al as column of frequencies to mimic, and sample 100 times to reduce piRNAs
            setrep=c()
            for(j in 1:100){
              setrep=append(setrep,repcds(paste(c(seqDNA[40:length(seqDNA)]),sep="",collapse=""),CAIS,AAtoCodF,CodonAl))
            }
            inseqs=setrep
            fwd=grep("GGTCTC",x=inseqs)
            rev=grep("GAGACC",x=inseqs)
            all=unique(c(fwd,rev))
            if((length(all) != 0)&(length(all) != length(inseqs))){inseqs=inseqs[-c(all)]}
            id=order(sapply(inseqs,function(x){Strcondepies(x,Pies,4)}))[1]
            SeqEnd=inseqs[id]
        }
        
        SeqtoOpt=paste(c(SeqStart,SeqEnd),sep="",collapse="")
        #############
        
        incProgress(3/10)
        
        ##If PiRNA removal
        if(FlaPi){
          pipipis=Strfindpies(SeqtoOpt,Pies,4)
          if(length(pipipis) > 0 ){
            stpos=c()
            for(pipi in pipipis){
              stpos=append(stpos,start(matchPattern(DNAString(pipi),DNAString(SeqtoOpt),max.mismatch=4,fixed=T)))
            }
            stpos=unique(c(stpos,stpos+3,stpos+6,stpos+9,stpos+12,stpos+15))
            SeqtoOpt=modbyposiz(SeqtoOpt,stpos,CAIS,AAtoCodF,1)
            
            pipipis=Strfindpies(SeqtoOpt,Pies,4)
            if(length(pipipis) > 0 ){
              Iter=1
              nflag=TRUE
              while((Iter < 100)&(nflag)){
                
                stpos=c()
                for(pipi in pipipis){
                  stpos=append(stpos,start(matchPattern(DNAString(pipi),DNAString(SeqtoOpt),max.mismatch=4,fixed=T)))
                }
                stpos=unique(c(stpos,stpos+3,stpos+6,stpos+9,stpos+12,stpos+15))
                SeqtoOpt=modbyposiz(SeqtoOpt,stpos,CAIS,AAtoCodF,1)
                
                pipipis=Strfindpies(SeqtoOpt,Pies,4)
                
                if(length(pipipis) > 0 ){Iter=1+Iter}else{nflag=FALSE}
              }
              ##Error Message for iterations
              if(Iter==100){output$ErrorMessage <- renderText({paste("Error: PiRNA removal did not work even after 100 iterations. Final number of piRNA sites found was: ",length(pipipis))})}
            }
          }
          
        }
        
        incProgress(3/10)
        
        #If Bsa Sites
        if(FlaBs){
          ##Rough match
          fwd=grep("GGTCTC",x=SeqtoOpt)
          rev=grep("GAGACC",x=SeqtoOpt)
          all=unique(c(fwd,rev))
          
          if(length(all) > 0 ){ ##Do proper biostrings match
            stpos=c()
            stpos=append(stpos,start(matchPattern(DNAString("GGTCTC"),DNAString(SeqtoOpt),fixed=T)))
            stpos=append(stpos,start(matchPattern(DNAString("GAGACC"),DNAString(SeqtoOpt),fixed=T)))
            stpos=unique(c(stpos+3))
            SeqtoOpt=modbyposiz(SeqtoOpt,stpos,CAIS,AAtoCodF,1)
            }
        }
        
        #IfSapI
        if(FlaSp){
          ##Rough match
          fwd=grep("GCTCTTC",x=SeqtoOpt)
          rev=grep("GAAGAGC",x=SeqtoOpt)
          all=unique(c(fwd,rev))
          
          if(length(all) > 0 ){ ##Do proper biostrings match
            stpos=c()
            stpos=append(stpos,start(matchPattern(DNAString("GCTCTTC"),DNAString(SeqtoOpt),fixed=T)))
            stpos=append(stpos,start(matchPattern(DNAString("GAAGAGC"),DNAString(SeqtoOpt),fixed=T)))
            stpos=unique(c(stpos+3))
            SeqtoOpt=modbyposiz(SeqtoOpt,stpos,CAIS,AAtoCodF,1)
          }
        }
        
        if(FlaEp){
          ##Rough match
          fwd=grep("CGTCTC",x=SeqtoOpt)
          rev=grep("GAGACG",x=SeqtoOpt)
          all=unique(c(fwd,rev))
          
          if(length(all) > 0 ){ ##Do proper biostrings match
            stpos=c()
            stpos=append(stpos,start(matchPattern(DNAString("CGTCTC"),DNAString(SeqtoOpt),fixed=T)))
            stpos=append(stpos,start(matchPattern(DNAString("GAGACG"),DNAString(SeqtoOpt),fixed=T)))
            stpos=unique(c(stpos+3))
            SeqtoOpt=modbyposiz(SeqtoOpt,stpos,CAIS,AAtoCodF,1)
          }
        }
        
        incProgress(1/20)
        
        #If introns
        if(FlaIn){
          typeIn=as.integer(input$intropt)
          finalvec=unlist(strsplit(toupper(SeqtoOpt),""))
          stpos=c()
          inpos=c()
          stpos=append(stpos,start(matchPattern(DNAString("AGR"),DNAString(SeqtoOpt),fixed=F)))
          stpos=stpos + 1
          if(length(stpos)>3){
          if(sum((stpos > 50)&(stpos<150))>1){inpos=append(inpos,c(stpos[(stpos > 50)&(stpos<150)])[1])}else{inpos=c(50)}
          if(sum(stpos > (inpos[1]+150))>1){inpos=append(inpos,c(stpos[(stpos > (inpos[1]+150))])[1])}else{inpos=c(inpos[1],inpos[1]+150)}
          if(sum(stpos > (inpos[2]+150))>1){inpos=append(inpos,c(stpos[(stpos > (inpos[2]+150))])[1])}else{inpos=append(inpos,sample((inpos[2]+50):(length(finalvec)-1),1))}
          }else{
            inpos=c(100,250,400)
            }
          inposis=inpos[order(inpos)]
          SeqtoOpt=paste(c(finalvec[1:inposis[1]],as.character(IntronSeqs[typeIn,1]),finalvec[(inposis[1]+1):inposis[2]],as.character(IntronSeqs[typeIn,2]),finalvec[(inposis[2]+1):inposis[3]], as.character(IntronSeqs[typeIn,3]),finalvec[(inposis[3]+1):length(finalvec)]),sep="",collapse="")
        }
        
        ########
        
        incProgress(1/20)
        
        aaaads=""
        if(FlaRi){aaaads="aaaa"}
        #PartialResult
        output$PartialResult <- renderText({
          paste(c("Result:\n",aaaads,SeqtoOpt),sep="",collapse="")
        })
        
        output$downloadoptseq <- renderUI({
          if((ErrorFlag == 0) & !is.null(SeqtoOpt)) {
            optsin=colnames(CAIS)[CodonAl]
            if(FlaRi){optsin=paste(optsin,"_OptimalRibosomalBinding",sep="",collapse="")}
            if(FlaPi){optsin=paste(optsin,"_RemovepiRNAHomology",sep="",collapse="")}
            if(FlaBs){optsin=paste(optsin,"_noBsaI",sep="",collapse="")}
            if(FlaSp){optsin=paste(optsin,"_noSpaI",sep="",collapse="")}
            if(FlaEp){optsin=paste(optsin,"_noEsp3I",sep="",collapse="")}
            if(FlaIn){optsin=paste(optsin,"_withSyntheticIntrons",sep="",collapse="")}
            write(paste(">Optimized_cDNA:Codon-",optsin,"\n",aaaads,SeqtoOpt,"\n",sep="",collapse=""),paste("PATC/users/",session_id,"/SeqOpop.fasta", sep=""))
            downloadButton('DownSeqOut', 'Download Optimized DNA sequence')
          }
        })
        ###End Main routine
        })
      }
      
    }) 
    
    
    ################
    
    #On starting, put Smu-1##########
    observeEvent(session$clientData,{init()})
    
    # prints actuall tab
    observeEvent(input$panels,{
        sendM(input$panels);
    })
    
    ##function to send message
    sendM = function(x){
      session$sendCustomMessage("igvstat-change", x)
    }
    
  })
