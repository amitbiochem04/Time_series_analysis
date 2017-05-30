library( hwriter )

########## here the output of html file name is open. images is the folder where all the png files are.
page <- openPage( "open.html",
                  head = paste( sep="\n",
                                "<script>",
                                "   function set_image( name ) {",
                                "      document.getElementById( 'plot' ).setAttribute( 'src', 'images/' + name + '.png' );",
                                "   }",
                                "</script>" ) )
cat(file=page,
    '<table><tr><td style="vertical-align:top"><div style="height:500px; overflow-y:scroll">' )
#####
hwrite(ressig, border=NULL, page=page,
       onmouseover = sprintf( "set_image( '%s' );", ressig$ens.id  ) )
cat( file=page,
     '</div></td><td style="vertical-align:top"><img id="plot" width="200px"></td></tr></table>' )
closePage(page)
browseURL( "open.html" )

#########another way of writeing html file.
#load( "final.rda" )
final$log2FoldChange <- sprintf( "%.2f", final$log2FoldChange )
plotfiles <- paste0( rownames(final), ".png" )
plotfiles2<-paste0( rownames(final), ".jpeg" )
final$plot <-hwriteImage(plotfiles, height="70px", table=FALSE, link = plotfiles )
final$plot2<-hwriteImage(plotfiles2, height="70px", table=FALSE, link = plotfiles2 )
maintable <- hwrite( final, 
                     table.style = "border-collapse: collapse;",
                     cellpadding = 5,row.bgcolor='#ffdc98',
                     center=TRUE,onmouseover="this.bgColor='#ffaaaa'", 
                     onmouseout="this.bgColor='white'", bgcolor='white')

hwrite(maintable, page="table.html")

#####hetamp2
heatmap.2(as.matrix( dd[,c(1:12)]), col=greenred(30), trace="none",
                    Colv=FALSE, Rowv=FALSE, dendrogram = "none",scale="row",
          labRow = FALSE,cexRow=0.075,main="Significant rythmic genes in dark phase",keysize=1)






########ressig plot 
ressig$log2FoldChange <- sprintf( "%.2f", ressig$log2FoldChange )
plotfiles <- paste0( rownames(ressig), ".png" )
ressig$plot <-hwriteImage(plotfiles, height="70px", table=FALSE, link = plotfiles )
#final$plot2<-hwriteImage(plotfiles2, height="70px", table=FALSE, link = plotfiles2 )
maintable <- hwrite( ressig, 
                     table.style = "border-collapse: collapse;",
                     cellpadding = 5,row.bgcolor='#ffdc98',
                     center=TRUE,onmouseover="this.bgColor='#ffaaaa'", 
                     onmouseout="this.bgColor='white'", bgcolor='white')
hwrite(maintable, page="differential regulated gene table.html")
