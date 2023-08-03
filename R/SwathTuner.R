#      
#      The SwathTUNER software:
#           Building optimized precursor isolation windows for
#           data-independent acquisition methods in Mass Spectrometry.
#      
#      Author: Aivett Bilbao 
#           aivett.bilbao@isb-sib.ch, bilbao.aivett@gmail.com
#      Copyright (C) 2015  
#           Life Sciences Mass Spectrometry, University of Geneva, Geneva, Switzerland.
#           Proteome Informatics Group, SIB Swiss Institute of Bioinformatics, Geneva, Switzerland.
#      
#      This program is free software: you can redistribute it and/or modify
#      it under the terms of the GNU General Public License as published by
#      the Free Software Foundation, either version 3 of the License, or
#      (at your option) any later version.
#      
#      This program is distributed in the hope that it will be useful,
#      but WITHOUT ANY WARRANTY; without even the implied warranty of
#      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#      GNU General Public License for more details.
#      
#      You should have received a copy of the GNU General Public License
#      along with this program.  If not, see <http://www.gnu.org/licenses/>
#      


# --------------------------------------------------------
# --- Graphical user interface using Shiny package
# --------------------------------------------------------
launchGUI = function()
{
	if(exists("myPath"))
        addResourcePath(prefix="www", directoryPath=paste(myPath, "/inst/extdata/www/", sep=""))
    
    shinyApp(options=list(port=5863),
             
        ui = fluidPage(
            fluidRow(
                column(6,
                       h2("SwathTUNER", align="center")),
                column(6,
                       br(),
                       p("If you use SwathTUNER, please cite it: DOI: 10.1021/acs.jproteome.5b00543", align="left"),
                       a(h6("The use of variable Q1 isolation windows improves selectivity in LC-SWATH-MS acquisition"), href="http://pubs.acs.org/doi/abs/10.1021/acs.jproteome.5b00543", target="_blank", align="left"))
            ),
            hr(),
            # File upload manager
            fluidRow(
                column(6,
                       h4("1. Input file (LC-MS peaks)", tags$button(id="btnHelp", onClick="{var e = document.getElementById('helpPanel'); if (e.style.display == 'block') e.style.display = 'none'; else e.style.display = 'block';}", icon("question-circle", class = NULL, lib="font-awesome"))),
                       actionButton("btnFileName", "Browse file..."),
                       textInput("txtFileName", "",""),
                       tags$div(style="color:red; text-align:center", textOutput(outputId="errorMessage"))),
                column(6, id="helpPanel", style="display:none",
                       helpText("A tab-separated text file containing a list of LC-MS peaks with the columns: mz, rt, intensity."),
                       tags$img(src="www/example.png"))
            ),
            hr(),
            h4("2. Visualize your input-data"),

    
            fluidRow(
                column(6, plotOutput("lcmsScatter")),
                column(6, plotOutput("lcHist"))
            ),
            hr(),
            br(),
            
            fluidRow(
                column(6, wellPanel(sliderInput("n", "Bin width (m/z):", min=5, max=100, value=20, width="100%")),
                       plotOutput("precMzHist")),
                column(6, wellPanel(uiOutput("nRT")),
                       plotOutput("precRTHist"))
            ),
            hr(),
            
            h4("3. Build precursor isolation windows"),
            wellPanel(
                fluidRow(
                    column(10, uiOutput("mzRange")),
                    column(2, numericInput("overlap", "Overlap (m/z units):", min=0, max=10, value=1, step=1))),
                fluidRow(
                    column(5, uiOutput("rtRange")),
                    column(5, sliderInput("nWindows", "Number of windows:", min=1, max=200, value=20, width="100%")),
                    column(2, numericInput("decPlaces", "Decimal places:", min=0, max=4, value=1, step=1))),
                uiOutput("copyTable1"),
                uiOutput("copyTable2"),
                uiOutput("copyTable3")
                ),
            actionButton("btnSavePDF", "Save images to PDF file", icon=icon("floppy-o", class = NULL, lib="font-awesome")),
            uiOutput("savePDF"),
            fluidRow(align="center", 
                     column(4, h5("Fixed width"),
                            plotOutput("fixedWidthHist"), 
                            plotOutput("fixedWidthHistInt"), br(),
                            actionButton("btnCopyFixed", "", icon=icon("clipboard", class = NULL, lib="font-awesome")),
                            tableOutput("fixedWidthTable")),
                     column(4, h5("Equalized frequency"),
                            plotOutput("equalizedFreqHist"),
                            plotOutput("equalizedFreqHistInt"), br(),
                            actionButton("btnCopyEqFreq", "", icon=icon("clipboard", class = NULL, lib="font-awesome")),
                            tableOutput("equalizedFreqTable")),
                     column(4, h5("Equalized TIC"),
                            plotOutput("equalizedTICHist"),
                            plotOutput("equalizedTICHistInt"), br(),
                            actionButton("btnCopyEqTIC", "", icon=icon("clipboard", class = NULL, lib="font-awesome")),
                            tableOutput("equalizedTICTable"))
                     )
        ), 
        
        server = function(input, output, session) 
        {
            
            getMzPeaks = reactive(
            {
                output$errorMessage = renderText("")
                if (input$btnFileName == 0) 
                    return(NULL)
                    
                mzPeaks = tryCatch({
                    inFile = file.choose()
                    read.table(inFile, header=TRUE, sep="\t")
                }, warning = function(w) {
                    message(w)
                    return(w)
                }, error = function(e) {
                    message(e)
                    return(e)
                })
                
                isolate(updateTextInput(session, "txtFileName",  value = inFile))
                
                if(is.character(mzPeaks))
                {
                    output$errorMessage = renderText(paste("Invalid file, please check the help (?). Error message:",mzPeaks))
                    return(NULL)
                }
            
                summary(mzPeaks)
                
                if(!is.numeric(mzPeaks$mz) || !is.numeric(mzPeaks$intensity) || !is.numeric(mzPeaks$rt))
                {
                    cols = colnames(mzPeaks)
                    if(is.null(cols))
                    {   
                        output$errorMessage = renderText("Column names are required, please check the help (?).")
                    }else
                    {    
                        cols = paste(as.character(unlist(cols)), collapse=", ")
                        output$errorMessage = renderText(paste("Invalid column names: ", cols, ". Please check the help (?).", sep=""))
                    }
                    return(NULL)
                }
                
                mzPeaks$log10_intensity = cut(log10(mzPeaks$intensity), breaks=9, include.lowest=TRUE)
                return(mzPeaks)
            })
            
            getRangeRT = reactive(
            {
                mzPeaks = getMzPeaks()
                if(is.null(mzPeaks))
                    return(NULL)
                
                rtMin = min(mzPeaks$rt)
                rtMax = max(mzPeaks$rt)
                rtStep = (rtMax - rtMin) / 100
                return(c(rtMin,rtMax,rtStep))
            })
            
            output$nRT = renderUI(
            {
                rangeRT = getRangeRT()
                if(is.null(rangeRT))
                    return(NULL)                

                rangeRT = ceiling(rangeRT)
                sliderInput("nRT", "Bin width (RT):", min=rangeRT[1], max=rangeRT[2], value=(rangeRT[3]*5), step=rangeRT[3], width="100%")
            })
    
            output$rtRange = renderUI(
            {
                rangeRT = getRangeRT()
                if(is.null(rangeRT))
                    return(NULL)
                
                rtMin = floor(rangeRT[1])
                rtMax = ceiling(rangeRT[2])
                rtStep = ceiling(rangeRT[3])
                sliderInput("rtRange", "RT range:", min=rtMin, max=rtMax, value=c(rtMin,rtMax), step=rtStep, width="100%")
            })
            
            output$mzRange = renderUI(
            {
                mzPeaks = getMzPeaks()
                if(is.null(mzPeaks))
                    return(NULL)
                
                mzMax = ceiling(max(mzPeaks$mz))
                
                sliderInput("mzRange", "m/z range:", min=50, max=(mzMax+50), value=c(50,mzMax - (mzMax %% 10)), step=5, width="100%")
            })
            
            getMzPeaksFiltered = reactive(
            {
                mzPeaks = getMzPeaks()
                if(is.null(mzPeaks) || is.null(input$rtRange) || is.null(input$mzRange))
                    return(NULL)
            
                mzPeaks = mzPeaks[which(mzPeaks$rt >= input$rtRange[1] & mzPeaks$rt <= input$rtRange[2] & mzPeaks$mz >= input$mzRange[1] & mzPeaks$mz <= input$mzRange[2]),]
                return(mzPeaks)
            })
            
            output$lcmsScatter = renderPlot(
            {
                mzPeaks = getMzPeaks()
                if(is.null(mzPeaks))
                    return(NULL)
                
                p = build_lcmsScatter(mzPeaks)
                plot(p)
            })
            
            output$lcHist = renderPlot(
            {       
                mzPeaks = getMzPeaks()
                if(is.null(mzPeaks))
                    return(NULL)
                
                binwd = (max(mzPeaks$rt) - min(mzPeaks$rt)) / 500                
                p = build_lcHist(mzPeaks, binwd)
                plot(p)
            })
            
            output$precMzHist = renderPlot(
            {       
                mzPeaks = getMzPeaks()
                if(is.null(mzPeaks) || is.null(input$n))
                    return(NULL)
                
                p = build_precMzHist(mzPeaks, input$n)
                plot(p)
            })
            
            output$precRTHist = renderPlot(
            {       
                mzPeaks = getMzPeaks()
                if(is.null(mzPeaks) || is.null(input$nRT))
                    return(NULL)
                
                p = build_precRTHist(mzPeaks, input$nRT)
                plot(p)
            })
            
            getFixedWidthWindows = reactive(
            {                
                mzPeaks = getMzPeaksFiltered()
                if(is.null(mzPeaks))
                    return(NULL)
                
                win2 = buildFixedWidthWindows(input$nWindows, input$mzRange[1], input$mzRange[2])
                return(win2)
            })

            output$fixedWidthHist = renderPlot(
            {
                mzPeaks = getMzPeaksFiltered()
                win2 = getFixedWidthWindows()
                if(is.null(win2))
                    return(NULL)
                
                brks = c(win2$start, win2$end[length(win2$end)])
                plot(build_HistMz(mzPeaks, brks))
            })   
            
            output$fixedWidthHistInt = renderPlot(
            {
                mzPeaks = getMzPeaksFiltered()
                win2 = getFixedWidthWindows()
                if(is.null(win2))
                    return(NULL)
                
                brks = c(win2$start, win2$end[length(win2$end)])
                plot(build_HistInt(mzPeaks, brks))
            })
            
            output$fixedWidthTable = renderTable(digits=4,
            {
                win2 = getFixedWidthWindows()
                if(is.null(win2))
                    return(NULL)

                win2 = addOverlap(win2, input$overlap)
                win2$width = win2$end - win2$start
                win2 = round(win2, digits=input$decPlaces)
                return(win2)
            })

            getWindowsEqualizedPrecFreq = reactive(
            {                
                mzPeaks = getMzPeaksFiltered()
                if(is.null(mzPeaks))
                    return(NULL)
                
                win2 = buildWindowsEqualizedPrecFreq(mzPeaks, input$nWindows, input$mzRange[1], input$mzRange[2])
                return(win2)
            })

            output$equalizedFreqHist = renderPlot(
            {
                mzPeaks = getMzPeaksFiltered()
                win2 = getWindowsEqualizedPrecFreq()
                if(is.null(win2))
                    return(NULL)
                
                brks = c(win2$start, win2$end[length(win2$end)])
                plot(build_HistMz(mzPeaks, brks))                
            })   
            
            output$equalizedFreqHistInt = renderPlot(
            {
                mzPeaks = getMzPeaksFiltered()
                win2 = getWindowsEqualizedPrecFreq()
                if(is.null(win2))
                    return(NULL)
                
                brks = c(win2$start, win2$end[length(win2$end)])                
                plot(build_HistInt(mzPeaks, brks))
            })   

            output$equalizedFreqTable = renderTable(digits=4,
            {
                win2 = getWindowsEqualizedPrecFreq()
                if(is.null(win2))
                    return(NULL)

                win2 = addOverlap(win2, input$overlap)
                win2$width = win2$end - win2$start
                win2 = round(win2, digits=input$decPlaces)
                return(win2)
            })
            
            getWindowsEqualizedTIC = reactive(
            {                
                mzPeaks = getMzPeaksFiltered()
                if(is.null(mzPeaks))
                    return(NULL)
                
                win2 = buildWindowsEqualizedTIC(mzPeaks, input$nWindows, input$mzRange[1], input$mzRange[2])
                return(win2)
            })

            output$equalizedTICHist = renderPlot(
            {
                mzPeaks = getMzPeaksFiltered()
                win2 = getWindowsEqualizedTIC()
                if(is.null(win2))
                    return(NULL)
                
                brks = c(win2$start, win2$end[length(win2$end)])
                plot(build_HistMz(mzPeaks, brks))
            })   

            output$equalizedTICHistInt = renderPlot(
            {
                mzPeaks = getMzPeaksFiltered()
                win2 = getWindowsEqualizedTIC()
                if(is.null(win2))
                    return(NULL)
                
                brks = c(win2$start, win2$end[length(win2$end)])
                plot(build_HistInt(mzPeaks, brks))
            })
            
            output$equalizedTICTable = renderTable(digits=4,
            {
                win2 = getWindowsEqualizedTIC()
                if(is.null(win2))
                    return(NULL)
                
                win2 = addOverlap(win2, input$overlap)
                win2$width = win2$end - win2$start
                win2 = round(win2, digits=input$decPlaces)
                return(win2)
            })
            
            output$savePDF = renderUI(
            {
                input$btnSavePDF
                
                mzPeaks = isolate(getMzPeaks())
                mzPeaksFilt = isolate(getMzPeaksFiltered())
                win2Fix = isolate(getFixedWidthWindows())
                win2Mz = isolate(getWindowsEqualizedPrecFreq())
                win2TIC = isolate(getWindowsEqualizedTIC())
                inFile = isolate(input$txtFileName)
                if(is.null(mzPeaks) || is.null(mzPeaksFilt) || is.null(win2Fix) || is.null(win2Mz) || is.null(win2TIC) || inFile == "" || is.null(input$n) || is.null(input$nRT))
                    return(NULL)
                
                # create .pdf file
                figuresFileName = sub("\\..*", ".pdf", inFile)
                
                pdf(figuresFileName, paper="a4r")
                
                plot(build_lcmsScatter(mzPeaks))
                
                binwd = (max(mzPeaks$rt) - min(mzPeaks$rt)) / 500                
                plot(build_lcHist(mzPeaks, binwd))
                
                plot(build_precMzHist(mzPeaks, input$n))
                
                plot(build_precRTHist(mzPeaks, input$nRT))
                
                brks = c(win2Fix$start, win2Fix$end[length(win2Fix$end)])
                plot(build_HistMz(mzPeaksFilt, brks))
                plot(build_HistInt(mzPeaksFilt, brks))
                
                brks = c(win2Mz$start, win2Mz$end[length(win2Mz$end)])
                plot(build_HistMz(mzPeaksFilt, brks))
                plot(build_HistInt(mzPeaksFilt, brks))

                brks = c(win2TIC$start, win2TIC$end[length(win2TIC$end)])
                plot(build_HistMz(mzPeaksFilt, brks))
                plot(build_HistInt(mzPeaksFilt, brks))
                
                # close .pdf file
                dev.off()
                
                return(NULL)
            })
            
            output$copyTable1 = renderUI(
            {
                input$btnCopyFixed
                win2 = isolate(getFixedWidthWindows())
                if(is.null(win2))
                    return(NULL)
                
                win2 = addOverlap(win2, isolate(input$overlap))
                win2 = round(win2, digits=isolate(input$decPlaces))
                write.table(format(win2, nsmall=isolate(input$decPlaces), trim=TRUE), "clipboard", sep="\t", row.names=FALSE, quote=FALSE)
                return(NULL)
            })
            
            output$copyTable2 = renderUI(
            {
                input$btnCopyEqFreq
                win2 = isolate(getWindowsEqualizedPrecFreq())
                if(is.null(win2))
                    return(NULL)
                
                win2 = addOverlap(win2, isolate(input$overlap))
                win2 = round(win2, digits=isolate(input$decPlaces))
                write.table(format(win2, nsmall=isolate(input$decPlaces), trim=TRUE), "clipboard", sep="\t", row.names=FALSE, quote=FALSE)
                return(NULL)
            })
            
            output$copyTable3 = renderUI(
            {
                input$btnCopyEqTIC
                win2 = isolate(getWindowsEqualizedTIC())
                if(is.null(win2))
                    return(NULL)
                
                win2 = addOverlap(win2, isolate(input$overlap))
                win2 = round(win2, digits=isolate(input$decPlaces))
                write.table(format(win2, nsmall=isolate(input$decPlaces), trim=TRUE), "clipboard", sep="\t", row.names=FALSE, quote=FALSE)
                return(NULL)
            })

        }
    )
}

# --------------------------------------------------------
# --- BUILD PLOTS
# --------------------------------------------------------
build_lcmsScatter = function(mzPeaks)
{
    p = ggplot(mzPeaks, aes(x=rt, y=mz))
    p = p + labs(x="Retention time (RT)", y="m/z", title="LC-MS")
    p = p + theme_bw()
    p = p + theme(text=element_text(size = 15), axis.title = element_text(size = 20))
    return(p + geom_point(alpha=1/5))
}

build_lcHist = function(mzPeaks, binwd)
{
    p = ggplot(mzPeaks, aes(x=rt, weight=intensity))
    p = p + labs(x="Retention time (RT)", y="Intensity", title="Pseudo-TIC")
    p = p + theme_bw()
    p = p + theme(text=element_text(size = 15), axis.title = element_text(size = 20))
    return(p + geom_histogram(binwidth=binwd, colour="darkblue", fill = "white"))
}

build_precMzHist = function(mzPeaks, n)
{
    p = ggplot(mzPeaks, aes(x=mz, fill=log10_intensity))
    p = p + labs(x="m/z", y="Number of ions", title="m/z distribution of ions")
    p = p + theme_bw()
    p = p + theme(text=element_text(size = 15), axis.title = element_text(size = 20))
    return(p + geom_histogram(binwidth=n, colour="darkblue") + scale_fill_brewer(name = "log10(Intensity)"))
}

build_precRTHist = function(mzPeaks, nRT)
{
    p = ggplot(mzPeaks, aes(x=rt, fill=log10_intensity))
    p = p + labs(x="Retention time (RT)", y="Number of ions", title="RT distribution of ions")
    p = p + theme_bw()
    p = p + theme(text=element_text(size = 15), axis.title = element_text(size = 20))
    return(p + geom_histogram(binwidth=nRT, colour="darkblue") + scale_fill_brewer(guide=FALSE))
}

build_HistMz = function(mzPeaks, brks)
{
    p = ggplot(mzPeaks, aes(x=mz, fill=log10_intensity))
    p = p + labs(x="m/z", y="Number of ions", title="m/z distribution of ions")
    p = p + theme_bw()
    p = p + theme(text=element_text(size = 15), axis.title = element_text(size = 20))
    return(p + geom_histogram(breaks=brks, colour="darkblue") + scale_fill_brewer(guide=FALSE))
}

build_HistInt = function(mzPeaks, brks)
{
    p = ggplot(mzPeaks, aes(x=mz, weight=intensity))
    p = p + labs(x="m/z", y="TIC", title="TIC per window")
    p = p + theme_bw()
    p = p + theme(text=element_text(size = 15), axis.title = element_text(size = 20))
    return(p + geom_histogram(breaks=brks, colour="darkblue", fill = "gray"))
}

# --------------------------------------------------------
# --- BUILD FIXED WIDTH WINDOWS
# --------------------------------------------------------
buildFixedWidthWindows = function(nWindows, startMz, endMz)
{
    # validate parameters
    if (missing(nWindows) || missing(startMz) || missing(endMz))
        stop("All parameters are required.")
    if (is.null(nWindows) || is.null(startMz) || is.null(endMz) || nWindows<1 || startMz<0 || endMz<0)
        stop("All parameters must be positive numbers.")
    
    # create data frame to store windows
    win2 = data.frame('start'=rep(0,1,nWindows),'end'=rep(0,1,nWindows))
    # compute window width
    wwidth = (endMz - startMz) / nWindows
    # build mz range
    start1 = startMz
    k = 1
    while(k <= nWindows)
    {
        win2$start[k] = start1
        start1 = start1 + wwidth
        win2$end[k] = start1
        k = k + 1
    }
    win2$end[k-1] = endMz
    
    return(win2)
}

# --------------------------------------------------------
# --- BUILD VARIABLE WINDOWS equalize precursor frequency per window
# --------------------------------------------------------
buildWindowsEqualizedPrecFreq = function(ions, nWindows, startMz, endMz)
{
    # validate parameters
    if (missing(ions) || missing(nWindows) || missing(startMz) || missing(endMz))
        stop("All parameters are required.")
    if (is.null(ions) || is.null(nWindows) || is.null(startMz) || is.null(endMz) || nWindows<1 || startMz<0 || endMz<0)
        stop("All parameters must be positive numbers.")

    ions = ions[which(ions$mz >= startMz & ions$mz <= endMz),]
    # create data frame to store windows
    win2 = data.frame('start'=rep(0,1,nWindows),'end'=rep(0,1,nWindows))
    # compute constant precursor frequency
    precs = ions$mz
    precFreq = floor(length(precs) / nWindows)
    residue = length(precs) %% nWindows
    # build mz range
    precs = sort(precs)
	win2$start[1] = startMz
	win2$end[1] = sum(precs[precFreq] + precs[precFreq + 1]) / 2
	k = 2
	index = precFreq * 2
    while(k <= nWindows)
    {
        win2$start[k] = win2[(k-1),'end']
        if((nWindows-k+1) == residue)
            precFreq = precFreq + 1
        win2$end[k] = sum(precs[index] + precs[index + 1]) / 2
        k = k + 1
        index = index + precFreq
    }
    win2$end[k-1] = endMz

    return(win2)
}

# --------------------------------------------------------
# --- BUILD VARIABLE WINDOWS equalize TIC per window
# --------------------------------------------------------
buildWindowsEqualizedTIC = function(ions, nWindows, startMz, endMz)
{
    # validate parameters
    if (missing(ions) || missing(nWindows) || missing(startMz) || missing(endMz))
        stop("All parameters are required.")
    if (is.null(ions) || is.null(nWindows) || is.null(startMz) || is.null(endMz) || nWindows<1 || startMz<0 || endMz<0)
        stop("All parameters must be positive numbers.")
    
    ions = ions[which(ions$mz >= startMz & ions$mz <= endMz),]
    ions$intensity = ions$intensity / max(ions$intensity) # normalize intensities to 0-1
    # create data frame to store windows
    win2 = data.frame('start'=rep(0,1,nWindows),'end'=rep(0,1,nWindows))
    # compute constant TIC
    maxtic = sum(ions$intensity) / nWindows
    # build mz range
    ions = ions[with(ions, order(mz)),]
    start1 = startMz
    tic = 0
    k = 1
    for(i in (1:(length(ions$mz) - 1)))
    {
        tic = tic + ions$intensity[i]
        if(tic >= maxtic)
        {
            if(k < nWindows && i < length(ions$intensity))
                maxtic = sum(ions$intensity[(i+1):length(ions$intensity)]) / (nWindows - k )
            tic = 0
            win2$start[k] = start1
            start1 = sum(ions$mz[i] + ions$mz[i+1]) / 2
            win2$end[k] = start1
            k = k + 1
        }
    }
    
	win2$start[nWindows] = start1
    win2$end[nWindows] = endMz
    win2 = win2[1:nWindows,]
    
    return(win2)
}


addOverlap = function(win2, overlap)
{
    # validate parameters
    if (missing(win2) || missing(overlap))
        stop("All parameters are required.")
    if (is.null(win2) || is.null(overlap) || overlap<0)
        stop("All parameters must be positive numbers.")
    
    overlap = overlap / 2
    win2$start = win2$start - overlap
    win2$end = win2$end + overlap    
    
    return(win2)
}
