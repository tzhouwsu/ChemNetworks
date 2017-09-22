#########################################
#                                       #
#    POLYHEDRA Output Analysis Tools    #
#                                       #
#           By Marissa Masden           #
#       marissa.masden@gmail.com        #
#      Washington State University      #
#        Department of Chemistry        #
#                                       #
#########################################



#analyze PD can take a ".pdfs.csv" file or a data frame from this file. 
#It saves two histograms to your working directory.
#It WILL overwrite previous files! 

analyzePD <- function(pdffile)
{
    #I don't anticipate these numbers changing, 
    #but if they do, you can change them here.
     
    MINPOLYHEDRASIZE<-4
    MAXPOLYHEDRASIZE<-12

    pdfs <-pdffile
    if(mode(pdffile)=="character")
    {
    pdfs <- read.csv(pdffile, header=FALSE)
    }
    else if(mode(pdffile)!="list")
    {
        cat("Input must be file name or existing data frame.\n")
        return;
    }
    colnames(pdfs) = c("shsz","dist")
    
    if(mode(pdffile)=="character")
    {
    pdffile <- gsub(".csv","",pdffile,fixed=TRUE)
    pdffile <- gsub(".input","",pdffile,fixed=TRUE)
    }
    else
    {
    pdffile <- "consolePD"
    }
    

    basehist <- paste(pdffile, ".hist.png", sep="")
    tempfn <-basehist
    
    count=0
    
    while(file.exists(basehist))
    {
        basehist <- paste(tempfn, count)
        count=count+1
        if(count>100)
        {
            cat("Too many files already have the same name. Quitting.", "\n")
            return;
        }
    }

    png(basehist)
    hist(pdfs$dist,breaks=30, freq=FALSE, main="Histogram of Distances Between \n Solvation Shell Atoms", xlab="Distance (units)")
    dev.off()
    
    histwithshsz <- paste(pdffile, ".hist.shsz.png",sep="")
       
    tempfn <-histwithshsz
    
    count=0
    
    while(file.exists(histwithshsz))
    {
        histwithshsz <- paste(tempfn, count)
        count=count+1
        if(count>100)
        {
            cat("Too many files already have the same name. Quitting.", "\n")
            return;
        }
    }
    

    png(histwithshsz)
    hist(pdfs$dist,breaks=30,freq=FALSE)
    newymax=par("usr")[4]*1.25

    hist(pdfs$dist,breaks=30, freq=FALSE, main="Histogram of Distances Between \n Solvation Shell Atoms\nBy Shell Size", xlab="Distance (units)", ylim=c(0,newymax))
    
    shellsizes <- vector()
    for(i in MINPOLYHEDRASIZE:MAXPOLYHEDRASIZE)
    {
        if(sum(pdfs$shsz==i)>2){
        shellsizes=c(shellsizes,i)
        }
    }
    
    for (nums in shellsizes)
    {
        lines(density(pdfs$dist[pdfs$shsz==nums]), col=(nums-min(shellsizes)+1), lwd=2)
    }
    
    params <- par("usr")
    clrs = shellsizes-min(shellsizes)+1

    legend(params[1]+.5, (floor(params[4]*10))/10, shellsizes, fill=clrs)

    dev.off()
}

#makePolyFreqs will analyze frequency data for the polyhedra which
#appear in your simulation. In order fot the "persistence" calcuation
#to work, ensure that the file names are in chronological order in
#the .Polys file. This outputs a .txt file instead of console output.
#There is also a .csv file output. Input should be .Polys filename.

makePolyFreqs<-function(polyfile)
{
    ##### READ FILE #####

    polys <- read.csv(polyfile, header=FALSE, strip.white=TRUE)

    colnames(polys) <- c("fname","molnum","atmnum","shsz","polynm","ID","time")
    
    ##### CREATE FILE #####
    fname <- polyfile
    fname <- gsub(".input","",fname,fixed=TRUE)
    
    tempfn <- fname

    fname <- paste(tempfn, ".txt", sep="")
    csvname <- paste(tempfn, ".csv", sep="")
    
    count=1
    
     while(file.exists(fname)||file.exists(csvname))
    {
        fname <- paste(tempfn, count, ".txt", sep="")
        csvname <- paste(tempfn, count, ".csv", sep="")
        count=count+1
        if(count>100)
        {
            cat("Too many files already have the same name. Quitting.", "\n")
            return;
        }
    }
    
    fileout <- file(fname, open="a")
    csvout <- file(csvname, open="a")
    
    
    ##### OUTPUT DATE/TIME HEADER #####
    
    curtm = date()

    cat("This analysis of ", polyfile, " was performed ", curtm, ".\n\n", file=fileout, sep="")

    #SHELL SIZE DATA
        
    cat("******** SHELL SIZE DISTRIBUTION DATA ********\n\n", file=fileout)
    
    

    ##### GET DISTRIBUTION OF SHELL SIZES #####
    
    shellfreqs = as.data.frame(table(polys$shsz))
    
    shellfreqs=cbind(shellfreqs, round(shellfreqs$Freq/length(polys$shsz)*100, digits=2))

    colnames(shellfreqs) = c("Shell Size", "Frequency", "Percent")
    
    writeformatted(shellfreqs, fileout)
    
    
    ##### CALCULATE PERSISTENCE VALUES ##### 
    
    steplist=list(length=length(unique(polys$fname)))
    
    stepnum=1
    
    totallength=length(unique(polys$molnum))*length(unique(polys$atmnum))
    
    countvec=vector(length=totallength)
    IDvec=vector(length=totallength)
    
    listallIDs = sort(unique(polys$ID))
    nrunbyID = vector(length=length(listallIDs))
    runtimebyID = vector(length=length(listallIDs))
    

    for (step in sort(unique(polys$fname)))
    {
        steplist[[stepnum]]=polys[polys$fname == step, ]

        if(stepnum==1)
        {
        sameaslast=1
        }else{
        sameaslast=(steplist[[stepnum]][,"ID"]==steplist[[stepnum-1]][,"ID"] | steplist[[stepnum]][,"ID"]==-1 )
        }
        
        #Make list of which IDs were changed from, and what the run length was for them. 
        liststoaddto=IDvec[!sameaslast]
        numbertoadd=countvec[!sameaslast]

        if(length(liststoaddto)>0)
        {
            for(i in 1:length(liststoaddto))
            {
                nrunbyID[listallIDs==liststoaddto[i]] = nrunbyID[listallIDs==liststoaddto[i]]+1
                runtimebyID[listallIDs==liststoaddto[i]] = runtimebyID[listallIDs==liststoaddto[i]]+numbertoadd[i]

            }
        }
        
        countvec = countvec*sameaslast +1 #Drop count to 1 if new poly

        IDvec[!sameaslast] = steplist[[stepnum]][!sameaslast, "ID"]
        
        stepnum=stepnum+1
    }
    persistencebyID=runtimebyID/nrunbyID
    

    ##### GET POLYHEDRAL RECOGNITION PROPORTIONS FOR EACH SHELL SIZE #####
    
    cat("\n\n******** POLYHEDRAL RECOGNITION DATA ********", file=fileout)
    
    write.table(rbind(c("ID", "Polyhedron Name", "Frequency", "Tot. %", "Rel. %", "Persistence (# snaps)")), sep=",", csvout, row.names=FALSE, col.names=FALSE)

    for(i in sort(unique(polys$shsz[polys$shsz>3 & polys$shsz<12])))
    {
        cat("\n\nSHELL SIZE:", i, "\n", file=fileout)

        polysseen=as.data.frame(table(droplevels(polys$polynm[polys$shsz==i])))
        
        matchprcnt=sum(polysseen$Freq[!is.na(polysseen$Var1)])/sum(polys$shsz==i)*100
        
        cat("MATCHED ", round(matchprcnt,2), "% IN THIS SHELL SIZE\n\n", file=fileout)

        ##Get IDs##
        polyIDs=as.data.frame(table((polys$ID[polys$shsz==i])))
        polyIDs=polyIDs$Var1[polyIDs$Var1!=-1]
        

        ##Arrange table##
        polysseen=cbind(polyIDs, polysseen, round(polysseen$Freq/length(polys$shsz), 4),round(polysseen$Freq/sum(polys$shsz==i)*100, digits=2), round(persistencebyID[match(polyIDs, listallIDs)], 3))
        
        ##Output in 2 formats##
        colnames(polysseen)=c("ID", "Polyhedron Name", "Frequency", "Tot. %", "Rel. %", "Persistence (# snaps)")
        
        polysseen <- polysseen[order(-polysseen[,4]),]

        write.table(polysseen, csvout, row.names=FALSE, col.names=FALSE, sep=",")
        writeformatted(polysseen, fileout)
    
    }

    ##### CLOSE FILES #####
    
    close(fileout)
    close(csvout)
   
}

writeformatted <- function(dtafrm, fileout)
{
    linerow=vector(length=ncol(dtafrm),mode="character")

    for (i in 1:ncol(dtafrm))
    {
        dtafrm[,i]=as.character(dtafrm[,i])
        newwidth=max(c(nchar(dtafrm[,i]),nchar(colnames(dtafrm)[i])))+2
        linerow[i]=paste(rep("-",newwidth),collapse="")
        dtafrm[,i] = format(dtafrm[,i], width=newwidth, justify="left")
        colnames(dtafrm)[i] = format(colnames(dtafrm)[i], width=newwidth, justify="centre")
    }
    dtafrm=rbind(linerow,dtafrm)

    suppressWarnings(write.table(dtafrm, file=fileout, quote=FALSE, row.names=FALSE, append=TRUE, sep=" "))


}













