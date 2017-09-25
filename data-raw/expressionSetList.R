##Example with adenocarcinoma versus squamous
##Exclude samples known to have problems
##integrative correlation
library(lungExpression)
#library(XDE)
data(stanford)
library(XDE)
##ignore warnings
source("~/projects/sandbox/myXde/R/functions.R") 
set.seed(50)
xstan <- avExample(object=stanford,
                   nsplit=3,
                   classSize=4,
                   balanced=FALSE,
                   SCALE.SD=1.5)
xstan <- as(xstan, "ExpressionSetList")
expressionSetList <- xstan[sample((1:nrow(xstan[[1]])), 500, replace=FALSE), ]
save(expressionSetList, file="data/expressionSetList.RData", compress=TRUE)


