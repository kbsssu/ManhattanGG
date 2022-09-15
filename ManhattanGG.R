#' ManhattanGG v.0.0.6 alpha
#' @author Beomsu Kim
#' kbsssu@gmail.com
#' https://github.com/kbsssu/ManhattanGG

suppressPackageStartupMessages({
    require(data.table)
    require(stringr)
    require(ggplot2)
    require(scales)
    require(argparse)
    require(ggrepel)
})

### Parse the command line
parser <- argparse::ArgumentParser(description=":: ManhattanGG v.0.0.6 alpha; by Beomsu Kim. ::", formatter_class="argparse.ArgumentDefaultsHelpFormatter")
parser$add_argument("--assoc", required=TRUE,
                    help="GWAS results file")
parser$add_argument("--metal", nargs="?", const=TRUE, default=FALSE, metavar="{blank}",
                    help="ManhattanGG set column names automatically and generate CHR and POS column if necessary; MarkerName format sould be CHR:POS(:...)")
parser$add_argument("--out", default=FALSE,
                    help="output prefix")
parser$add_argument("--paper-material", nargs="?", const=TRUE, default=FALSE, metavar="{blank}",
                    help="output for manhattan plot: 'manhattan.point.png', 'manhattan.label.pdf'")
parser$add_argument("--qq", nargs="?", const=TRUE, default=FALSE, metavar="{blank}",
                    help="Plotting Manhattan plot and Q-Q plot")
parser$add_argument("--qq-only", nargs="?", const=TRUE, default=FALSE, metavar="{blank}",
                    help="Plotting only Q-Q plot")
parser$add_argument("--qq-lambda", default=FALSE,
                    help="represent lambda value on Q-Q plot")
parser$add_argument("--pdf", nargs="?", const=TRUE, default=FALSE, metavar="{blank}",
                    help="output format: *.pdf")
parser$add_argument("--Rdata", nargs="?", const=TRUE, default=FALSE, metavar="{blank}",
                    help="Write *.Rdata file for converted data (contents: df, df_deco, df_color, manplot, qqplot)")
parser$add_argument("--chr", default="CHR",
                    help="colname: Chromosome")
parser$add_argument("--pos", default="BP",
                    help="colname: Base position")
parser$add_argument("--snp", default="SNP",
                    help="colname: Marker name")
parser$add_argument("--pval", default="P",
                    help="colname: P-value")
parser$add_argument("--y-lim", default=FALSE,
                    help="max value of ylim on manhattan plot")
parser$add_argument("--y-breaks", default=FALSE,
                    help="ybreaks on manhattan plot (delimeter: ',')")
parser$add_argument("--y-ceiling", default=FALSE,
                    help="similar to --y-lim but draw excess points on ymax (yvale[yvale > --y-ceiling] = --y-ceiling)")
parser$add_argument("--point-size", default=1.2,
                    help="point size on plots")
parser$add_argument("--text-size", default=3.5,
                    help="text size on plots (deco)")
parser$add_argument("--width", default=30,
                    help="plot width (cm)")
parser$add_argument("--height", default=15,
                    help="plot height (cm)")
parser$add_argument("--deco", default=FALSE,
                    help="decoration file name")
parser$add_argument("--deco-color", default=FALSE,
                    help="colname in --deco: color for points")
parser$add_argument("--deco-label", default=FALSE,
                    help="colname in --deco: text label")
parser$add_argument("--deco-label-color", default=FALSE,
                    help="colname in --deco: color for text label")
parser$add_argument("--highlight-color", default=FALSE,
                    help="highlighting color for points (yes: default color(#3399FF))")
parser$add_argument("--highlight-flank", default=1000,
                    help="highlighting flank for points (no: coloring only variants in --deco file)")
parser$add_argument("--remain-only", default=FALSE,
                    help="remain only points in this file")
parser$add_argument("--trait", nargs="+", default=FALSE,
                    help="colnames for trait-marking bottom plot")
parser$add_argument("--trait-color", nargs="+", default=FALSE,
                    help="colors for each trait in --trait plot")
parser$add_argument("--suggestiveline", default=1e-05,
                    help="P-value threshold for suggestive significance (no: no line)")
parser$add_argument("--genomewideline", default=5e-08,
                    help="P-value threshold for genome-wide significance (no: no line)")
parser$add_argument("--line-type", default="solid",
                    help="linetype: solid(default)/dashed")
parser$add_argument("--text-angle", default=0,
                    help="text angle")
parser$add_argument("--dark", nargs="?", const=TRUE, default=FALSE, metavar="{blank}",
                    help="dark colors for default dots.")
parser$add_argument("--title", default="",
                    help="title on plots")

argnames_ordered <- c("assoc", "metal", "out", "paper_material", "qq", "qq_only", "qq_lambda", "pdf", "Rdata", "chr", "pos", "snp", "pval", "y_lim", "y_breaks", "y_ceiling", "point_size", "text_size", "width", "height", "deco", "deco_color", "deco_label", "deco_label_color", "highlight_color", "highlight_flank", "remain_only", "trait", "trait_color", "suggestiveline", "genomewideline", "line_type", "text_angle", "dark", "title")
argvals_default <- c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, "CHR", "BP", "SNP", "P", FALSE, FALSE, FALSE, 1.2, 3.5, 30, 15, FALSE, FALSE, FALSE, FALSE, FALSE, 1000, FALSE, FALSE, FALSE, 1e-05, 5e-08, "solid", 0, FALSE, "")


#################
##=============##
##  Functions  ##
##=============##
#################

#######################
## Utility Functions ##
#######################

##### calculating -log10(P) - can calculate with extremely small P-value
log10c <- function(x){
    name.ori = names(x)
    x = toupper(x); names(x) = 1:length(x)
    idx.exp = grep("E-", x)
    if(length(idx.exp)==0) return(log10(as.numeric(x)))
    x.exp = x[idx.exp]; name.exp = names(x.exp)
    x.nonexp = x[-idx.exp]; name.nonexp = names(x.nonexp)
    # exp
    x.exp = strsplit(x.exp, "E-")
    x.exp = matrix(as.numeric(unlist(x.exp)), ncol=2, byrow=TRUE)
    x.exp1 = log10(x.exp[, 1])
    x.exp2 = x.exp[, 2]
    x.exp = x.exp1 - x.exp2
    names(x.exp) = name.exp
    # nonexp
    x.nonexp = log10(as.numeric(x.nonexp))
    names(x.nonexp) = name.nonexp
    # result
    rslt = c(x.exp, x.nonexp)
    rslt = rslt[order(as.numeric(names(rslt)))]
    names(rslt) = name.ori
    return(rslt)
}

##### exit()
exit <- function(...){
    args <- list(...)
    if(all(sapply(args, is.null))){
        rc <- 0
    }else if((length(args) > 1) || !is.numeric(args[[1]])){
        args = c(args, list("\n"))
        args$file = stderr()
        args$sep = ""
        do.call(cat, args)
        rc <- 1
    }else{
        rc <- as.integer(args[[1]])
    }
    q(save="no", status=rc)
}

##### running time
runtime <- function(starttime, endtime){
    tdiff = difftime(endtime, starttime, units="secs")
    tdiff = format(.POSIXct(tdiff, tz="GMT"), "%H:%M:%S")
    return(tdiff)
}

##### convert args (from argparse package) to text for message
args_to_txt <- function(args, argnames_ordered, argvals_default){
    argnames = names(args)
    idx = argnames_ordered %in% argnames
    argnames = argnames_ordered[idx]
    argvals_default = argvals_default[idx]
    txtout = "Input Arguments:\n\n"
    for(i in 1:length(argnames)){
        argname = argnames[i]
        argval = args[[argname]]
        if(length(argval)!=1 | (argval[1]!=argvals_default[i])){
            txtout = paste0(txtout, "--", paste(strsplit(argname, "_")[[1]], collapse="-"), " ", paste(argval, collapse=" "), "\n")
        }
    }
    paste0(txtout, "\n")
    return(txtout)
}

##### ggcolors functions
ggcolors <- function(n){
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}


#############################
## Formatting - Input data ##
#############################
Formatting_inputs <- function(assocf, decof=FALSE, metal=FALSE, chr="CHR", bp="BP", snp="SNP", p="P", deco.color=FALSE, deco.label=FALSE, deco.label.color=FALSE, trait=FALSE, qq_only=FALSE){
    ## read association file
    df = data.table::fread(assocf, na.string=c("NA", ""), colClasses="character", data.table=FALSE)
    if(metal){
        snp = "MarkerName"; p = "P-value"
        # generating chr and bp
        if(!all(c(chr, bp) %in% colnames(df)) & !(qq_only)){
            df = df[, -which(colnames(df) %in% c(chr, bp))]
            df[, c(chr, bp)] = as.data.frame(stringr::str_split_fixed(df[, "MarkerName"], ":", 4))[, 1:2]
            df = data.frame(df)
        }
    }
    if(!(p %in% colnames(df)) & ("P.value" %in% colnames(df))) p = "P.value"
    if(qq_only==FALSE){
        df = df[, c(snp, chr, bp, p)]
        colnames(df) = c("snp", "chr", "bp", "p")
        # convert chromosom to numeric
        XYMmap = list(X=23, Y=24, XY=25, M=26)
        XYM = intersect(c("X", "Y", "XY", "M"), toupper(unique(df$chr)))
        if(length(XYM)!=0){
            for(xym in XYM){
                df[toupper(df$chr)==xym, "chr"] = XYMmap[[xym]]
            }
        }
        df$chr = as.numeric(df$chr)
        df$bp = as.numeric(df$bp)
    }else{
        df = df[, c(snp, p)]
        colnames(df) = c("snp", "p")
    }
    df = df[!is.na(df$p), ]
    df[as.numeric(df$p) < 1e-300, "p"] = "1e-300"
    df$logP = -log10c(df$p)
    ## read decoration file
    if(decof!="FALSE"){
        deco = data.table::fread(decof, na.string=c("NA", ""), colClasses="character", data.table=FALSE)
        cols = c(deco.color, deco.label, deco.label.color, trait)
        if(all(c("chr", "bp") %in% colnames(deco))) cols = c(cols, "chr", "bp")
        cols = cols[cols!="FALSE"]
        deco = deco[, c(snp, cols)]
        colnames(deco) = c("snp", cols)
    }else{
        deco = FALSE
    }

    # return
    return(list(df=df, deco=deco))
}


########################################
## Plotting Function - Manhattan plot ##
########################################
manhattan <- function(df, metal=FALSE, chr="CHR", bp="BP", snp="SNP", p="P", y.lim=FALSE, y.breaks=FALSE, y.ceiling=FALSE, p.size=2, text.size=4,
                      deco=FALSE, deco.color=FALSE, deco.label=FALSE, deco.label.color=FALSE,
                      remain.only=FALSE,
                      trait=FALSE, trait.color=FALSE,
                      highlight.color=FALSE, highlight.flank=1000, suggestiveline=1e-05, genomewideline=5e-08,
                      l.type="solid", title="", text_angle=0, dark=FALSE, ManhattanGG=FALSE){
    # metal, ManhattanGG :  options for using this function independently
    #####
    if(!ManhattanGG){
        if(metal){
            snp = "MarkerName"; p = "P-value"
            # generating chr and bp
            if(!all(c(chr, bp) %in% colnames(df)) & !(qq_only)){
                df = df[, -which(colnames(df) %in% c(chr, bp))]
                df[, c(chr, bp)] = as.data.table(stringr::str_split_fixed(df[, "MarkerName"], ":", 4))[, 1:2]
                df = data.frame(df)
            }
        }
        if(!(p %in% colnames(df)) & ("P.value" %in% colnames(df))) p = "P.value"
        df = df[, c(snp, chr, bp, p)]
        colnames(df) = c("snp", "chr", "bp", "p")
        # convert chromosom to numeric
        XYMmap = list(X=23, Y=24, XY=25, M=26)
        XYM = intersect(c("X", "Y", "XY", "M"), toupper(unique(df$chr)))
        if(length(XYM)!=0){
            for(xym in XYM){
                df[toupper(df$chr)==xym, "chr"] = XYMmap[[xym]]
            }
        }
        df$chr = as.numeric(df$chr)
        df$bp = as.numeric(df$bp)
        df = df[(!is.na(df$p)) & (df$p != "0"), ]
        df$logP = -log10c(df$p)
    }
    #####

    df$pos = NA        # The position for each SNPs on the x axis
    ticks = NULL       # The center of each chromosome
    df$maxlogP = NA    # for decoration
    df = df[order(df$chr, df$bp), ]

    # y ceiling (limit ymax and draw excess points on ymax)
    if(is.numeric(y.ceiling)){
        df$ceiling = ifelse(df$logP>y.ceiling, 1, 0)
        df[df$ceiling==1, "logP"] = y.ceiling
        point_shape = c(16, 17)
        point_size = c(p.size, p.size*2.5)
        if(y.lim==FALSE) y.lim = y.ceiling
    }else{
        df$ceiling = 0
        point_shape = 16
        point_size = p.size
    }

    # max value of y axis (ylim)
    if(is.numeric(y.lim)){
        yMax = y.lim
    }else{
        yMax = max(ceiling((max(df$logP)+1)), genomewideline)    # default
    }

    # interval of y axis
    if(is.numeric(y.breaks)){
        ybreaks = seq(0, yMax, y.breaks)
    }else{
        ybreaks = scales::pretty_breaks()
    }

    ### mapping
    chrs = sort(unique(df$chr))
    n_chr = length(chrs)
    if(n_chr==1){
        df$pos = df$bp
        ticks = floor(length(df$pos))/2 + 1
        df$maxp = max(df$p)
        xlabel = paste("Chromosome", unique(df$chr))
        labs = ""
        snpcolors = ifelse(dark==TRUE, "gray20", "gray60")
    }else{
        lastbase = 0
        for(i in 1:n_chr){
            if(i==1){
                df[df$chr==chrs[i], "pos"] = df[df$chr==chrs[i], "bp"]
            }else{
                lastbase = lastbase + tail(df[df$chr==chrs[i-1], "bp"], 1)
                df[df$chr==chrs[i], "pos"] = df[df$chr==chrs[i], "bp"] + lastbase
            }
            ticks = c(ticks, (min(df[df$chr==chrs[i], "pos"]) + max(df[df$chr==chrs[i], "pos"]))/2 + 1)
            df[df$chr==chrs[i], "maxlogP"] = max(df[df$chr==chrs[i], "logP"])
        }
        xlabel = "Chromosome"
        labs = unique(df$chr)
        snpcolors = rep(c("gray60", "gray80"), length.out=n_chr)
        if(dark==TRUE) snpcolors = rep(c("gray20", "gray50"), length.out=n_chr)
    }

    xMax = ceiling(max(df$pos) * 1.03)
    xMin = floor(max(df$pos) * -0.03)


    ### decoration data formatting and coordinating
    if(is.data.frame(deco)){
        ## basic formatting
        # read input columns
        cols = c(deco.color, deco.label, deco.label.color)
        cols_format = c("color", "label", "labelcolor")
        if(all(c("chr", "bp") %in% colnames(deco))){
            cols = c(cols, c("chr", "bp"))
            cols_format = c(cols_format, "chr", "bp")
        }
        cols_format = cols_format[cols!="FALSE"]
        cols = cols[cols!="FALSE"]
        df_deco = deco[, c(snp, cols)]
        colnames(df_deco) = c("snp", cols_format)
        # color
        if(deco.color!=FALSE) df_deco[is.na(df_deco$color) | (df_deco$color==""), "color"] = "transparent"
        if(deco.label.color!=FALSE) df_deco[(!is.na(df_deco$label)) & (is.na(df_deco$labelcolor) | (df_deco$labelcolor=="")), "labelcolor"] = "transparent"
        ## coordinating
        # get coordinates from df
        if(!all(c("chr", "bp") %in% colnames(df_deco))){ 
            df_deco = merge(df_deco, df[, c("snp", "chr", "bp", "p", "pos", "logP")], key="snp", all.x=TRUE)
            df_deco = na.omit(df_deco)
            df_deco = df_deco[order(df_deco$color), ]
            rownames(df_deco) <- NULL
        }else{
            df_deco$bp = as.numeric(df_deco$bp)
        }
        # highlight variants
        if(highlight.flank){
            flank = highlight.flank*1000
            df_color = data.frame(matrix(nrow=0, ncol=ncol(df)+1))
            colnames(df_color) = c(colnames(df), "color")
            for(i in 1:nrow(df_deco)){
                if(!all(c("chr", "bp") %in% colnames(df_deco))){
                    if(any(is.na(df_deco$pos), is.na(df_deco$logP))) next
                    df_filtered = df[(df$pos>=df_deco$pos[i]-flank) & (df$pos<=df_deco$pos[i]+flank), ]
                }else{
                    df_filtered = df[(df$chr==df_deco$chr[i]) & (df$bp>=df_deco$bp[i]-flank) & (df$bp<=df_deco$bp[i]+flank), ]
                }
                df_filtered$color = df_deco$color[i]
                df_color = rbind(df_color, df_filtered)
            }
            rm(df_filtered)
            df_color = unique(df_color)
        }else{
            df_color = df_deco
        }
        # set point colors
        n_color = length(unique(df_color$color))
        if(deco.color!=FALSE){
            colorset = c(snpcolors, unique(df_color$color))
            names(colorset) = c(as.character(chrs), as.character(unique(df_color$color)))
        }else if(highlight.color!=FALSE){
            colorset = c(snpcolors, rep(highlight.color, length.out=n_color))
            names(colorset) = c(as.character(chrs), as.character(rep(highlight.color, length.out=n_color)))
        }else{
            colorset = c(snpcolors, ggcolors(n_color))
            names(colorset) = c(as.character(chrs), as.character(ggcolors(n_color)))
        }
        # set text colors
        if(!any(deco.label==FALSE, deco.label=="FALSE")){
            df_label = df_deco[order(df_deco$label, df_deco$p), ]
            df_label = na.omit(df_label)
            df_label = df_label[df_label$label!="", ]
            df_label = df_label[!duplicated(df_label[, c("snp", "label")]), ]
            rownames(df_label) <- NULL
            # label color
            if(deco.label.color!=FALSE){
                n_labelcolor = length(unique(df_label$labelcolor))
                colorset_label = unique(df_label$labelcolor)
                names(colorset_label) = as.character(unique(df_label$labelcolor))
            }else{
                n_label = length(unique(df_label$label))
                df_label$labelcolor = "black"
                colorset_label = "black"
                names(colorset_label) = "black"
            }
        }else{
            colorset_label = NULL
        }
    }else{
        df_deco = NULL
        df_color = NULL
        df_label = NULL
        colorset = snpcolors
        names(colorset) = c(as.character(chrs))
	colorset_label = NULL
    }
    colorset = c(colorset, colorset_label)

    ### --remain-only: color setting
    if(remain.only!=FALSE){
        df[!(df$snp %in% remain.only), "chr"] = NA
        if(!is.null(df_color)){
            df_color[!(df_color$snp %in% remain.only), "color"] = NA
        }
        if(!is.null(df_label)){
            df_label[!(df_label$snp %in% remain.only), "labelcolor"] = NA
        }
    }

    ### formatting and coordinating for trait-marking plot
    if(is.data.frame(deco)){
        if(all(trait!=FALSE)){
            df_trait = deco[, c(snp, trait)]
            colnames(df_trait) = c("snp", trait)
            df_trait = merge(df_trait, df_deco[!duplicated(df_deco[, c("snp", "pos")]), c("snp", "pos")], key="snp", all.x=TRUE)
            df_trait = df_trait[order(df_trait$pos), ]
            rownames(df_trait) <- NULL
            # length for detecting point width with a unit as 'pos'
            maxpos = max(df$pos)
            xLength = maxpos - min(df$pos)
            p.width = xLength / 55    # 0.5cm (width) : 27.5cm (total x length)
            ## x-coordinating
            # left to right
            df_trait$diff_pos = c(0, df_trait$pos[2:nrow(df_trait)] - df_trait$pos[1:(nrow(df_trait)-1)])
            i = 2
            while(TRUE){
                if(i>nrow(df_trait)) break
                if(df_trait$diff_pos[i]<=p.width){
                    # detecting bundle of close points
                    j = 1
                    while(TRUE){
                        # stop if not near point
                        if(i+j>nrow(df_trait)) break
                        if(df_trait$diff_pos[i+j]>p.width) break
                        j = j+1
                    }
                    idx_bundle = (i-1):(i+j-1)
                    x_bundle = mean(df_trait$pos[idx_bundle]) + c(-p.width*0.5*(length(idx_bundle)-1), p.width*0.5*(length(idx_bundle)-1))
                    if(i>2){
                        if(x_bundle[1]<df_trait$pos[i-2]+p.width){
                            # if bundle area invade the previous point, ...
                            x_bundle = x_bundle + (df_trait$pos[i-2] - x_bundle[1]) + 2*p.width
                        }
                    }
                    # generating x coordinates
                    x_bundle = seq(x_bundle[1], x_bundle[2], length.out=length(idx_bundle))    # generated x coordinates
                    # injection
                    df_trait[idx_bundle, "pos"] = x_bundle
                    # recalculating pos-diff
                    df_trait$diff_pos = c(0, df_trait$pos[2:nrow(df_trait)] - df_trait$pos[1:(nrow(df_trait)-1)])
                    # i jumping
                    i = i + j
                }else{
                    i = i + 1
                }
            }
            # right to left for matching end point to manhattan ---> sliding to left
            sliding_size = df_trait$pos[nrow(df_trait)] + p.width - maxpos
            if(sliding_size>0){
                # sliding
                df_trait$diff_pos = c(0, df_trait$pos[2:nrow(df_trait)] - df_trait$pos[1:(nrow(df_trait)-1)])
                i = nrow(df_trait)
                while(TRUE){
                    # detecting sliding start point
                    x_aftersliding = df_trait$pos[i] - sliding_size
                    crt = ifelse(df_trait$diff[i]>p.width, df_trait$pos[i-1] + 2*p.width, df_trait$pos[i-1] + p.width)
                    if(x_aftersliding>crt) break
                    i = i - 1
                }
                # sliding injection
                df_trait[i:nrow(df_trait), "pos"] = df_trait$pos[i:nrow(df_trait)] - sliding_size
            }

            ## y-coordinating and shape(color) defining
            # indicators for color and fill
            sign_color = c()
            sign_fill = c()
            i = 0
            for(t in trait){
                i = i + 1
                sign_color = c(sign_color, ifelse(df_trait[, t]==0, length(trait)*2, i))
                sign_fill = c(sign_fill, ifelse(df_trait[, t]==1, i, length(trait)*2))
            }
            colorset_trait = c(trait.color, "grey80")
            fillset_trait = c(trait.color, "white")
            # wide to long format
            df_trait = reshape(df_trait, direction="long", varying=trait, v.names="sign", idvar=c("snp", "pos"), timevar="trait", times=trait)
            rownames(df_trait) <- NULL
            # indicators injection
            df_trait$sign_color = sign_color
            df_trait$sign_fill = sign_fill
            # y coords
            df_trait$trait_y = rep(seq(-yMax/7, by=-yMax/20, length.out=length(trait)), each=nrow(df_trait))

            colorset = c(colorset, colorset_trait)
            fillset = fillset_trait
        }else{
            colorset = colorset
            fillset = NULL
        }
    }


    ### Plotting
    plot <- ggplot() +
      geom_point(data=df, aes(x=pos, y=logP, colour=factor(chr), shape=factor(ceiling), size=factor(ceiling)), show.legend=FALSE) +
      scale_x_continuous(name=xlabel, breaks=ticks, labels=labs, expand=c(0,0)) +
#      scale_y_continuous(name=expression(paste(-log[10], " (", italic(P), "-value)", collapse="")), limits=c(0, yMax), breaks=ybreaks, expand=c(0,0)) +
      scale_y_continuous(name=expression(paste(-log[10], " (", italic(P), ")", collapse="")), breaks=ybreaks, expand=c(0,0)) +
      scale_colour_manual(values=colorset, na.value="transparent") +
      scale_shape_manual(values=point_shape) +
      scale_size_manual(values=point_size) +
      coord_cartesian(xlim=c(xMin, xMax), ylim=c(0, yMax), clip="off") +
      ggtitle(title) + theme_classic() +
      theme(axis.line.x = element_blank(),
            axis.title = element_text(size=14),
            axis.text = element_text(size=11))
    # highlight deco data
    if(is.data.frame(deco)){
        plot <- plot + geom_point(data=df_color, aes(x=pos, y=logP, colour=factor(color), shape=factor(ceiling), size=factor(ceiling)), show.legend=FALSE)
        if(all(trait!=FALSE)){
            plot <- plot + geom_point(data=df_trait, aes(x=pos, y=trait_y, colour=factor(sign_color), fill=factor(sign_fill)), alpha=0.6, shape=22, size=5, stroke=2, show.legend=FALSE) +
                           scale_fill_manual(values=fillset)
        }
    }
    # draw line
    if(l.type=="solid") l.type=1
    if(l.type=="dashed") l.type=2
    if(is.numeric(suggestiveline))  plot <- plot + geom_hline(yintercept=-log10(suggestiveline), size=0.5, colour="#3399FF", lty=l.type)
    if(is.numeric(genomewideline))  plot <- plot + geom_hline(yintercept=-log10(genomewideline), size=0.5, colour="#FF3399", lty=l.type)
    # deco label
    if(is.data.frame(deco)){
        if(is.character(deco.label)&(deco.label!="FALSE")){
            label.parse = ifelse(any(grepl("italic", df_deco$label)), TRUE, FALSE)
            if(text_angle==0){
                plot <- plot + ggrepel::geom_text_repel(data=df_label, aes(x=pos, y=logP, label=label, colour=factor(labelcolor)), size=text.size, point.padding=unit(0, "lines"), nudge_y=yMax/50, vjust=1, parse=label.parse, angle=text_angle, show.legend=FALSE)
            }else{
                plot <- plot + ggrepel::geom_text_repel(data=df_label, aes(x=pos, y=logP, label=label, colour=factor(labelcolor)), size=text.size, point.padding=unit(0, "lines"), nudge_y=yMax/50, vjust=1, parse=label.parse, angle=text_angle, show.legend=FALSE)
            }
        }
    }

    ### return
    return(list(plot=plot, df=df, df_deco=df_deco, df_color=df_color))
}


##################################
## Plotting Function - Q-Q plot ##
##################################
qq <- function(p, CI=TRUE, y.breaks=FALSE, gc.lambda=FALSE, p.size=1.2, title=""){
  
#  p = as.numeric(p)
  p = p[(!is.na(p)) & (p != "0")]
  logP = -log10c(p)
#  o = -log10(sort(p, decreasing=F))
  o = sort(logP, decreasing=T)
  e = -log10(1:length(o)/(length(o)+1))
  
  df = as.data.frame(cbind(e, o))
  
  if(CI){
    c975 = rep(0, length(o))
    c025 = rep(0, length(o))
    for(i in 1:length(o)){
      c975[i] = -log10(qbeta(0.975, i, length(o)-i+1))
      c025[i] = -log10(qbeta(0.025, i, length(o)-i+1))
    }
    df = as.data.frame(cbind(df, c025, c975))
  }
  
  if(y.breaks){
    ybreaks = seq(0, max(o*1.05), y.breaks)
  }else{
    ybreaks = scales::pretty_breaks()
  }
  
  ymax = max(o)
  xmax = max(e)
  plot <- ggplot(data=df) +
    geom_point(aes(x=e, y=o), color="gray20", size=p.size) +
    geom_abline(slope=1, intercept=0, size=0.5, color="red") +
    coord_cartesian(xlim=c(0, max(e)*1.05), ylim=c(0, max(o)*1.05)) +
    # scale_x_continuous(name=expression(Expected~~paste(-log[10], " (", italic(P), "-value)", collapse="")), expand=c(0,0)) +
    # scale_y_continuous(name=expression(Observed~~paste(-log[10], " (", italic(P), "-value)", collapse="")), breaks=ybreaks, expand=c(0,0)) +
    scale_x_continuous(name=expression(Expected~~paste(-log[10], " (", italic(P), ")", collapse="")), expand=c(0,0)) +
    scale_y_continuous(name=expression(Observed~~paste(-log[10], " (", italic(P), ")", collapse="")), breaks=ybreaks, expand=c(0,0)) +
    ggtitle(title) + theme_classic() +
    theme(panel.border = element_rect(fill=NA, colour="grey40"),
          axis.title = element_text(size=16),
          axis.text = element_text(size=14),
          panel.grid.major.y = element_line(color="grey80"))
  if(CI==TRUE)  plot <- plot + geom_ribbon(aes(x=e, ymin=c025, ymax=c975), alpha=0.3, fill="gray50")
  if(is.numeric(gc.lambda)) plot <- plot + geom_text(x=xmax*0.05, y=ymax*0.95, label=as.expression(bquote(λ[GC] == .(gc.lambda))), size=7, hjust=0)
  
  return(plot)
}




#======================================================================================================================#
#======================================================================================================================#

###########################
##=======================##
##  Running ManhattanGG  ##
##=======================##
###########################

starttime = Sys.time()
### arg options - startwith "opt." (just convenience for me)
args <- parser$parse_args()
opt.assoc = args$assoc
opt.metal = args$metal
opt.out = args$out
opt.paper_material = args$paper_material
opt.qq = args$qq
opt.qq_only = args$qq_only
opt.qq_lambda = args$qq_lambda
opt.pdf = args$pdf
opt.Rdata = args$Rdata
opt.chr = args$chr
opt.pos = args$pos
opt.snp = args$snp
opt.pval = args$pval
opt.y_lim = args$y_lim
opt.y_breaks = args$y_breaks
opt.y_ceiling = args$y_ceiling
opt.point_size = args$point_size
opt.text_size = args$text_size
opt.width = args$width
opt.height = args$height
opt.deco = args$deco
opt.deco_color = args$deco_color
opt.deco_label = args$deco_label
opt.deco_label_color = args$deco_label_color
opt.highlight_color = args$highlight_color
opt.highlight_flank = args$highlight_flank
opt.remain_only = args$remain_only
opt.trait = args$trait
opt.trait_color = args$trait_color
opt.suggestiveline = args$suggestiveline
opt.genomewideline = args$genomewideline
opt.line_type = args$line_type
opt.text_angle = args$text_angle
opt.dark = args$dark
opt.title = args$title


### starting message
writeLines(c("",
             "=============================================",
             "      :: ManhattanGG v.0.0.6 :: (alpha)      ",
             "=============================================",
             ": Manhattan/Q-Q plot for GWAS results       :",
             ":                             using ggplot2 :",
             "  -. Beomsu Kim",
             "  -. SAIHST, SKKU",
             "  -. kbsssu@gmail.com",
             "  -. https://github.com/kbsssu/ManhattanGG",
             "---------------------------------------------",
             ""))
message(args_to_txt(args, argnames_ordered, argvals_default), "\n")


### Check and convert inputs and options
# y
if(opt.y_lim!=FALSE) opt.y_lim = as.numeric(opt.y_lim)
if(opt.y_breaks!=FALSE) opt.y_breaks = as.numeric(strsplit(opt.y_breaks, ",")[[1]])
if(opt.y_ceiling!=FALSE) opt.y_ceiling = as.numeric(opt.y_ceiling)
# size and angle
opt.point_size = as.numeric(opt.point_size)
opt.text_size = as.numeric(opt.text_size)
opt.text_angle = as.numeric(opt.text_angle)
opt.width = as.numeric(opt.width)
opt.height = as.numeric(opt.height)
# highlight
if(tolower(opt.highlight_color)=="yes") opt.highlight_color = "#3399FF"
opt.highlight_flank = ifelse(tolower(opt.highlight_flank)=="no", FALSE, as.numeric(opt.highlight_flank))
# lines
opt.suggestiveline = ifelse(tolower(opt.suggestiveline)=="no", FALSE, as.numeric(opt.suggestiveline))
opt.genomewideline = ifelse(tolower(opt.genomewideline)=="no", FALSE, as.numeric(opt.genomewideline))
# qq
opt.qq_ci = ifelse((tolower(opt.qq)=="ci")|(tolower(opt.qq_only)=="ci"), TRUE, FALSE)
opt.qq = ifelse(opt.qq==FALSE, FALSE, TRUE)
opt.qq_only = ifelse(opt.qq_only==FALSE, FALSE, TRUE)
if(opt.qq_lambda!=FALSE) opt.qq_lambda = as.numeric(opt.qq_lambda)
# remain only
if(opt.remain_only!=FALSE){
    opt.remain_only = data.table::fread(opt.remain_only, na.string=c("NA", ""), colClasses="character", data.table=FALSE)
    opt.remain_only = opt.remain_only[, opt.snp]
}
# output file name
manpref = ifelse(opt.out!=FALSE, paste(opt.out, "manhattan", sep="."), "manhattan")
qqpref = ifelse(opt.out!=FALSE, paste(opt.out, "qq", sep="."), "qq")
Rdatapref = ifelse(opt.out!=FALSE, opt.out, "ManhatanGG")
# output format
opt.pdf = ifelse(opt.pdf==TRUE, "pdf", "png")



### Read inputs
message("Read data..")
dfs <- suppressWarnings(Formatting_inputs(assocf=opt.assoc, decof=opt.deco, metal=opt.metal, chr=opt.chr, bp=opt.pos, snp=opt.snp, p=opt.pval, deco.color=opt.deco_color, deco.label=opt.deco_label, deco.label.color=opt.deco_label_color, trait=opt.trait, qq_only=opt.qq_only))
df = dfs$df
deco = dfs$deco
rm(dfs)




### plotting
if(opt.qq_only){
    message("::Procedure:: Q-Q plot")
    suppressWarnings(qqplot <- qq(df$p, CI=opt.qq_ci, y.breaks=opt.y_breaks, gc.lambda=opt.qq_lambda, p.size=opt.point_size, title=opt.title))
    ggsave(paste(qqpref, opt.pdf, sep="."), plot=qqplot, dpi=300, width=20, height=20, units="cm")
    df = NULL
    df_deco = NULL
    df_color = NULL
    if(opt.Rdata=="plot") save(qqplot, file=paste(Rdatapref, "Rdata", sep="."))
}else{
    if(opt.paper_material==FALSE){
        message("::Procedure:: Manhattan plot")
        suppressWarnings(manplot <- manhattan(df=df, meta=opt.metal, chr="chr", bp="bp", snp="snp", p="p", y.lim=opt.y_lim, y.breaks=opt.y_breaks, y.ceiling=opt.y_ceiling, p.size=opt.point_size, text.size=opt.text_size,
                                              deco=deco, deco.color=opt.deco_color, deco.label=opt.deco_label, deco.label.color=opt.deco_label_color,
                                              remain.only=opt.remain_only,
                                              trait=opt.trait, trait.color=opt.trait_color,
                                              highlight.color=opt.highlight_color, highlight.flank=opt.highlight_flank, suggestiveline=opt.suggestiveline, genomewideline=opt.genomewideline,
                                              l.type=opt.line_type, title=opt.title, text_angle=opt.text_angle, dark=opt.dark, ManhattanGG=TRUE))
        ggsave(paste(manpref, opt.pdf, sep="."), plot=manplot$plot, dpi=300, width=opt.width, height=opt.height, units="cm")
        if(opt.Rdata=="plot") save(manplot, file=paste(Rdatapref, "Rdata", sep="."))
    }else{
        ## point plot
        message("::Procedure:: Manhattan plot (materials for additional decoration)")
        message("                -. plotting: *.point.png")
        suppressWarnings(manplot <- manhattan(df=df, meta=opt.metal, chr="chr", bp="bp", snp="snp", p="p", y.lim=opt.y_lim, y.breaks=opt.y_breaks, y.ceiling=opt.y_ceiling, p.size=opt.point_size, text.size=opt.text_size,
                                              deco=deco, deco.color=opt.deco_color, deco.label=FALSE, deco.label.color=FALSE,
                                              remain.only=opt.remain_only,
                                              trait=opt.trait, trait.color=opt.trait_color,
                                              highlight.color=opt.highlight_color, highlight.flank=opt.highlight_flank, suggestiveline=opt.suggestiveline, genomewideline=opt.genomewideline,
                                              l.type=opt.line_type, title=opt.title, text_angle=opt.text_angle, dark=opt.dark, ManhattanGG=TRUE))
        plotpoint <- manplot$plot +
                        theme(title = element_text(colour="white"),
                              axis.title = element_text(size=14, colour="white"),
                              axis.text = element_text(size=11, colour="white"),
                              axis.line = element_line(colour="white"),
                              axis.ticks = element_line(colour="white"))
        ggsave(paste(manpref, "point", opt.pdf, sep="."), plot=plotpoint, dpi=300, width=opt.width, height=opt.height, units="cm")

        ## label plot
        # dataframe for label plot
        message("                -. plotting: *.label.pdf")
        df_labelplot = df[df$snp %in% deco$snp, ]
        rownames(df_labelplot) <- NULL
        chrs = sort(unique(df$chr))
        for(i in 1:length(chrs)){
            chr = chrs[i]
            minbp = min(df[df$chr==chr, "bp"])
            maxbp = max(df[df$chr==chr, "bp"])
            minmax = df[df$chr==chr & (df$bp==minbp | df$bp==maxbp), ]
            df_labelplot = rbind(df_labelplot, minmax)
        }
        suppressWarnings(manplot.label <- manhattan(df=df_labelplot, meta=opt.metal, chr="chr", bp="bp", snp="snp", p="p", y.lim=opt.y_lim, y.breaks=opt.y_breaks, y.ceiling=opt.y_ceiling, p.size=opt.point_size, text.size=opt.text_size,
                                                    deco=deco, deco.color=opt.deco_color, deco.label=opt.deco_label, deco.label.color=opt.deco_label_color,
                                                    remain.only=opt.remain_only,
                                                    trait=opt.trait, trait.color=opt.trait_color,
                                                    highlight.color=opt.highlight_color, highlight.flank=opt.highlight_flank, suggestiveline=opt.suggestiveline, genomewideline=opt.genomewideline,
                                                    l.type=opt.line_type, title=opt.title, text_angle=opt.text_angle, dark=opt.dark, ManhattanGG=TRUE))
        ggsave(paste(manpref, "label.pdf", sep="."), plot=manplot.label$plot, dpi=300, width=opt.width, height=opt.height, units="cm")
        if(opt.Rdata=="plot") save(manplot, manplot.label, file=paste(Rdatapref, "Rdata", sep="."))
    }
    if(opt.qq){
        message("::Procedure:: Q-Q plot")
        suppressWarnings(qqplot <- qq(df$p, CI=opt.qq_ci, y.breaks=opt.y_breaks, gc.lambda=opt.qq_lambda, p.size=opt.point_size, title=opt.title))
        ggsave(paste(qqpref, opt.pdf, sep="."), plot=qqplot, dpi=300, width=20, height=20, units="cm")
        if(opt.Rdata=="plot"){
            if(opt.paper_material==FALSE){
                save(qqplot, manplot, manplot.label, file=paste(Rdatapref, "Rdata", sep="."))
            }else{
                save(qqplot, manplot, file=paste(Rdatapref, "Rdata", sep="."))
            }
        }
    }
}

# save *.Rdata
if(opt.Rdata=="plot"){
    message(sprintf("Write %s", paste(Rdatapref, "Rdata", sep=".")))
}else if(opt.Rdata==TRUE){
    message(sprintf("Write %s", paste(Rdatapref, "Rdata", sep=".")))
    save.image(file=paste(Rdatapref, "Rdata", sep="."))
}

message(sprintf("Done:)  Running time: %s\n", runtime(starttime, endtime=Sys.time())))


