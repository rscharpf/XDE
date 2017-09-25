setMethod(".studySplit", "ExpressionSet",
          function(object, 
                   phenotypeLabel,
                   nsplit=2,
                   balanced=TRUE,
                   classSize=100,
                   studyName="study",
                   SCALE.SD=1, ...){

		  if(missing(phenotypeLabel)) stop("must indicate phenotypeLabel as an element in varLabel")
		  phenotypeLabels <- pData(object)[, grep(phenotypeLabel, colnames(pData(object)))]
            group0 <- grep(0, phenotypeLabels)
            group1 <- grep(1, phenotypeLabels)
            if(balanced){
              smallestGroup <- min(length(group0), length(group1))
              splitSize <- floor(smallestGroup/nsplit)
              if(splitSize > classSize)    splitSize <- classSize
              splitSize0 <- splitSize1 <- splitSize
            } else{##do not require balanced study
              splitSize0 <- floor(length(group0)/nsplit)
              splitSize1 <- floor(length(group1)/nsplit)
              if(splitSize0 > classSize) splitSize0 <- classSize
              if(splitSize1 > classSize) splitSize1 <- classSize
            }
            "%w/o%" <- function(x,y) x[!x %in% y] #--  x without y
            eset <- list()
            for(i in 1:nsplit){
              samples0 <- sample(group0, size=splitSize0, replace=FALSE)
              group0 <- group0 %w/o% samples0
              samples1 <- sample(group1, size=splitSize1, replace=FALSE)
              group1 <- group1 %w/o% samples1
              tmp <- object[, c(samples0, samples1)]
              ##Add independent noise
              s <- rowSds(exprs(tmp))
              e <- matrix(rnorm(length(s), 0, SCALE.SD*s), nrow(tmp), ncol(tmp), byrow=FALSE)
              exprs(tmp) <- exprs(tmp) + e
              ##    eset[[i]] <- dat[, c(samples0, samplnes1)]
              eset[[i]] <- tmp
            }
            if(!identical(featureNames(eset[[1]]), featureNames(eset[[2]])))
              stop("feature names not identical in the study split")
            gs <- matrix(1, nrow=nrow(eset[[1]]), ncol=nsplit)
            rownames(gs) <- featureNames(eset[[1]])
            xset <- new("ExpressionSetList",
                        list(split1=eset[[1]],
                             split2=eset[[2]],
                             split3=eset[[3]]))
            xset
          })
