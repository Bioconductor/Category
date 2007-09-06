
cb_parse_band_Mm <- function(x) {
    ## Given a chromosome band annotation (see examples below),
    ## return a vector giving the path from the root:
    ##
    ## 10 B3.5 => 10, 10 B, 10 B3, 10B3.5
    ##

    ##remove cM ou centrome info for the moment
    ## not sure how to able it otherwise
    x <- gsub("[[:space:]].*cM$","" , x)
    x <- gsub("[[:space:]].*centromere$","", x)

    if(regexpr("[[:space:]]", x)==-1){
        return(x)
    }else{
        chr.pos <- regexpr("[[:space:]]", x)
        chr <-  substr(x, 1, chr.pos-1)

        chrBand <- substr(x, chr.pos+1, chr.pos+1)
        sbs <- strsplit(substr(x, chr.pos+2, nchar(x)), "")[[1]]

        if(length(sbs)==0) {
            return(c(chr, paste(chr, chrBand, sep=" ")))
        }else{

            bands <- character(length(sbs)+1)
            bands[1] <- paste(chr, chrBand, sep=" ")
            i <- 1
            j <- 2
            prev <- bands[1]
            print(prev)
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
        if (i == j)
          return(bands[1:(i-1)])
        else
          return(bands)
    }
}


