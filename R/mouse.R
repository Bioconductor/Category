
cb_parse_band_mmus <- function(x) {
    ## Given a chromosome band annotation (see examples below),
    ## return a vector giving the path from the root:
    ##
    ## 20p12.2 => 20, 20p, 20p1, 20p12, 20p12.2
    ##
    ## x1 <- "2q32"
    ## x2 <- "4q21.22"
    ## x3 <- "20p12.2"
    ## x4 <- "2qter"
    ##

    ##remove cM ou centrome info for the moment
    x <- gsub("[[:space:]].*cM$","" , x)
     x <- gsub("[[:space:]].*centromere$","", x)
    
    if(regexpr("[[:space:]]", x)==-1){
        return(x)
    }else{
        chr.pos <- regexpr("[[:space:]]", x)
        chr <-  substr(x, 1, chr.pos-1)
       
        x <- gsub("[[:space:]]", "", x) 
        chrBand <- substr(x, chr.pos, chr.pos)
       

        sbs <- strsplit(substr(x, chr.pos+1, nchar(x)), "")[[1]]
        if(length(sbs)==0) {
            ##if (nchar(x) == nchar(chrBand))
              return(c(chr, paste(chr, chrBand, sep=" ")))
        }else{          
            bands <- character(length(sbs)+1)
            bands[1] <- chrBand
            i <- 1
            j <- 2
            prev <- bands[1]
            
            while (TRUE) {
                if (sbs[i] == ".") {
                    prev <- paste(prev, ".", sep="")
                    i <- i + 1
                    next
                }
                bands[j] <- paste(prev, sbs[i], sep="")
                prev <- bands[j]
                i <- i + 1
                j <- j + 1
                if (i > length(sbs))
                  break
            }
        }
            if (i == j)                         # there was a '.'
              return(c(chr, paste(chr, bands[1:(i-1)], sep=" ")))
            else
              return(c(chr, paste(chr, bands, sep=" ")))
        
    }
}

    
