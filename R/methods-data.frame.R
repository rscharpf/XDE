setMethod("pairs", "data.frame",
          function(x, ...){
            x <- as.matrix(x)
            XDE:::pairs(x, ...)
          })
