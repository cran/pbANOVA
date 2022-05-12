#' potato data
#'
#'This dataset potato is from an experiment on how plants adapt to
#'cold climates. The investigators decided to study this problem after
#'observing that plants that have been conditioned to cold previously
#'appear to suffer less damage from the cold. Two species of potato
#'were studied (species 1 and 2). Each plant was exposed to one of two
#'acclimatization regimes (1= plant was kept in cold room; 0= plant
#'                         was kept at room temperature) for several days. Later, plants were
#'subjected to one of two cold temperatures (-4 degrees C is coded as
#'                                           1; -8 degrees C is coded as 2). Two responses were measured: damage
#'score for photosynthesis (photo), and damage score for ion leakage
#'(leak).
#'
#'Use ion leakage to be the response variable. Some of the 80 plants
#'originally assigned to the treatment combinations were lost during
#'the experiment. Analyze the data from the plants that made it
#'through, and assess the effects of the three experimental factors
#'species, regime, and temperature on the response leakage.
#'
#'
#'
#' @format This data frame contains the following columns:
#'
#'\describe{
#'\item{variety:}{ Two species of potato
#'were studied (species 1 and 2)}
#'
#'\item{regime:}{Each plant was exposed to one of two
#'acclimatization regimes (1= plant was kept in cold room; 0= plant
#'                         was kept at room temperature) for several days.}
#'
#'\item{temp:}{plants were
#'subjected to one of two cold temperatures (-4 degrees C is coded as
#'                                           1; -8 degrees C is coded as 2)}
#'
#'\item{photo:}{ damage
#'score for photosynthesis}
#'
#'\item{leak:}{damage score for ion leakage}
#'}
#'
#' @docType data
#'
#' @usage data(potato)
#'
#'
#'
#' @keywords datasets
#'
#'@references Alver and Zhang (2022), \href{https://math.unm.edu/~gzhang12/paper/multiANOVA_PBcondensed_Feb2022.pdf}{Parametric Bootstrap Procedures for
#'Three-Factor ANOVA and Multiple Comparison Procedures with Unequal Group Variances}
#'@references Alver and Zhang (2022), \href{https://math.unm.edu/~gzhang12/paper/PB_Dunnett.pdf}{Multiple Comparisons of Treatment vs Control
#'Under Unequal Variances Using Parametric Bootstrap}
#'
#'
#'
#'
#'
"potato"





