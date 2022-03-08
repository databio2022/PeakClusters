#!/user/bin/Rscript
# super enhancer identify

#' Identify super enhancers
#'
#' @param d A dataframe from output of chipseeker
#' @param distance The distance for merging two adjacent peaks
#' @param min_peak_nums min peak numbers
#'
#' @return A dataframe
#' @export
#'
#' @examples
#' scanse(d)
scanse <- function(d,distance=125000,min_peak_nums=3){
    de=d[!grepl("Promoter",d$annotation),] # remove Promoter peaks
    de=de[!grepl("G|K",de$seqnames),] # remove Promoter peaks

    chrs=unique(de$seqnames)
    enhanceralone=NULL
    process=0
    for(x in chrs){
        dch=de[de$seqnames == x,]
        eachchrloc=NULL
        for(i in 1:length(dch$abs_summit)){
            process=process+1
            print(process)
            n=0
            for( j in 1:length(dch$abs_summit)){
                if(j != i){
                    if(abs(dch$abs_summit[j]-dch$abs_summit[i]) < 25000){
                        n=n+1
                        break
                    }
                }
            }
            if(n == 0){
                eachchrloc=append(eachchrloc,i)
            }
        }
        enhanceralone=rbind(enhanceralone,dch[eachchrloc,])
    }


    de$ppois <- ppois(q=de$pileup,
                              lambda=mean(enhanceralone$pileup),
                              lower.tail=FALSE,log=FALSE)

    de$fdr=p.adjust(de$ppois,method="fdr")

    dec=de[de$fdr <= 0.05,]

    chrs=unique(dec$seqnames)
    #peak_res=NULL
    bed_se=NULL
    for(x in chrs){
        chrj=dec[dec$seqnames == x,]
        xx=as.matrix(data.frame(pos=chrj$abs_summit))
        #kNNdistplot(xx, k =  2)
        db <- dbscan::dbscan(xx, eps = distance, minPts = min_peak_nums )
        chrj$cluster=db$cluster
        chrj=chrj[order(chrj$abs_summit),]
        #peak_res=rbind(peak_res,chrj)
        cluster_num=unique(chrj$cluster)
        for(a in cluster_num){
            if(a > 0){
                peak_locs=chrj[chrj$cluster == a,]$abs_summit
                start=peak_locs[1]-200
                end=peak_locs[length(peak_locs)]+200
                se_loc=data.frame(chr=x,start=start,end=end)
                bed_se=rbind(bed_se,se_loc)
            }
        }
    }
   return(bed_se)
}
