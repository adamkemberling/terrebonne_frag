# Community Ordination Practice
# Package : vegan
# vignette 1 : https://cran.r-project.org/web/packages/vegan/vignettes/intro-vegan.pdf
# vignette 2 : http://cc.oulu.fi/~jarioksa/opetus/metodi/vegantutor.pdf

library(vegan)
library(MASS)

# Ordination basic method

'
Non-metric multidimensional scaling can be performed using isoMDS function
in the MASS package. This function needs dissimilarities as input.
Function vegdist in vegan contains dissimilarities which are found good
in community ecology. The default is Bray-Curtis dissimilarity, nowadays
often known as Steinhaus dissimilarity, or in Finland as Sørensen index.
The basic steps are:
'
#load practice dataset, colnames = species
data(varespec)

#calculate dissimilarity indexes
vare.dis <- vegdist(varespec)
marsh.dis <- vegdist(pcov.mat)
bin.dis <- vegdist(bin.mat)



#non-metric multidimensional scaling
vare.mds0 <- isoMDS(vare.dis)
marsh.mds0 <- isoMDS(marsh.dis)
bin.mds0 <- isoMDS(bin.dis)

stressplot(vare.mds0, vare.dis)

'
Function stressplot draws a Shepard plot where ordination distances
are plotted against community dissimilarities, and the fit is shown as a
monotone step line. In addition, stressplot shows two correlation like
statistics of goodness of fit. The correlation based on stress is R2 = 1−S
2
.
The “fit-based R2” is the correlation between the fitted values θ(d) and
ordination distances ˜d, or between the step line and the points. This
should be linear even when the fit is strongly curved and is often known
as the “linear fit”.
'

#null models
'
In linear fit, the
null model is that all ordination distances are equal, and the fit is a flat
horizontal line. This sounds sensible, but you need N − 1 dimensions for
the null model of N points, and this null model is geometrically impossible
in the ordination space. The basic stress uses the null model where all
observations are put in the same point, which is geometrically possible.
Finally a word of warning: you sometimes see that people use correlation
between community dissimilarities and ordination distances. This is dangerous
and misleading since nmds is a nonlinear method: an improved ordination with 
more nonlinear relationship would appear worse with this criterion.
Functions scores and ordiplot in vegan can be used to handle the
results of nmds:
'

ordiplot(vare.mds0, type = "t")

'
Only site scores were shown, because dissimilarities did not have information
about species.
The iterative search is very difficult in nmds, because of nonlinear relationship
between ordination and original dissimilarities. The iteration
easily gets trapped into local optimum instead of finding the global optimum.
Therefore it is recommended to use several random starts, and
select among similar solutions with smallest stresses. This may be tedious,
but vegan has function metaMDS which does this, and many more
things. The tracing output is long, and we suppress it with trace = 0,
but normally we want to see that something happens, since the analysis
can take a long time:
'

#demo data
vare.mds <- metaMDS(varespec, trace = FALSE)
vare.mds

pcov.mds <- metaMDS(pcov.mat, trace = FALSE) #yay it worked

'
We did not calculate dissimilarities in a separate step, but we gave the
original data matrix as input. The result is more complicated than previously,
and has quite a few components in addition to those in isoMDS results:
nobj, nfix, ndim, ndis, ngrp, diss, iidx, jidx, xinit, istart,
isform, ities, iregn, iscal, maxits, sratmx, strmin, sfgrmn,
dist, dhat, points, stress, grstress, iters, icause, call,
model, distmethod, distcall, data, distance, converged, tries,
engine, species. The function wraps recommended procedures into one
command. So what happened here?

1.  The range of data values was so large that the data were square root
transformed, and then submitted to Wisconsin double standardization,
or species divided by their maxima, and stands standardized
to equal totals. These two standardizations often improve the quality
of ordinations, but we forgot to think about them in the initial
analysis.

2. Function used Bray–Curtis dissimilarities.

3. Function run isoMDS with several random starts, and stopped either
after a certain number of tries, or after finding two similar
configurations with minimum stress. In any case, it returned the
best solution.

4. Function rotated the solution so that the largest variance of site
scores will be on the first axis.

5. Function scaled the solution so that one unit corresponds to halving
of community similarity from the replicate similarity.

6. Function found species scores as weighted averages of site scores,
but expanded them so that species and site scores have equal variances.
This expansion can be undone using shrink = TRUE in display
commands.
'


'
Function metaMDS used Bray-Curtis dissimilarity as default, which
usually is a good choice. Jaccard (Ruˇziˇcka) index has identical rank
order, but has better metric properties, and probably should be preferred.
Function rankindex in vegan can be used to study which of the indices
best separates communities along known gradients using rank correlation
as default. The following example uses all environmental variables in data
set varechem, but standardizes these to unit variance:

'
data(varechem)
rankindex(scale(varechem), varespec, c("euc","man","bray","jac","kul"))

#need corresponding environmental data for pcov matrix i.e. not including empty rows
#rows that aren't bare, columns of environmental variables
pcov.env <- marsh[rowSums(marsh[,12:34] != 0), c("bogaert","month",  "depth_cm", "dist_from_edge_m")]
pcov.env$month <- ifelse(pcov.env$month == "April", 0, 1)

rankindex(scale(pcov.env), pcov.mat, c("euc","man","bray","jac","kul"))
rankindex(scale(pcov.env), bin.mat, c("euc","man","bray","jac","kul"))

'
I took a very practical approach on indices emphasizing their ability
to recover underlying environmental gradients. Many textbooks emphasize
metric properties of indices. These are important in some methods,
but not in nmds which only uses rank order information. The metric
properties simply say that
1. if two sites are identical, their distance is zero,
2. if two sites are different, their distance is larger than zero,
3. distances are symmetric, and
4. the shortest distance between two sites is a line, and you cannot
improve by going through other sites.

These all sound very natural conditions, but they are not fulfilled by all
dissimilarities. Actually, only Euclidean distances – and probably Jaccard
index – fulfill all conditions among the dissimilarities discussed here, and
are metrics. Many other dissimilarities fulfill three first conditions and
are semimetrics.

There is a school that says that we should use metric indices, and
most naturally, Euclidean distances. One of their drawbacks was that
they have no fixed limit, but two sites with no shared species can vary
in dissimilarities, and even look more similar than two sites sharing some
species. This can be cured by standardizing data. Since Euclidean distances
are based on squared differences, a natural transformation is to
standardize sites to equal sum of squares, or to their vector norm using
function decostand:

'
dis <- vegdist(decostand(varespec, "norm"), "euclid")
pcov.dis <- vegdist(decostand(pcov.mat, "norm"), "euclid")
bin.dis <- vegdist(decostand(bin.mat, "norm"), "euclid")

'
This gives chord distances which reach a maximum limit of √2 when
there are no shared species between two sites. Another recommended
alternative is Hellinger distance which is based on square roots of sites
standardized to unit total:

'
dis <- vegdist(decostand(varespec, "hell"), "euclidean")

'
Despite standardization, these still are Euclidean distances with all their
good properties, but for transformed data. Actually, it is often useful to
transform or standardize data even with other indices. If there is a large
difference between smallest non-zero abundance and largest abundance,
we want to reduce this difference. Usually square root transformation is
sufficient to balance the data. Wisconsin double standardization often
improves the gradient detection ability of dissimilarity indices; this can
be performed using command wisconsin in vegan. Here we first divide
all species by their maxima, and then standardize sites to unit totals.
After this standardization, many dissimilarity indices become identical in
rank ordering and should give equal results in nmds.

'

#----------- Comparing Ordinations, Procrustes rotation --------------------
'
Two ordinations can be very similar, but this may be difficult to see,
because axes have slightly different orientation and scaling. Actually, in
nmds the sign, orientation, scale and location of the axes are not de-
fined, although metaMDS uses simple method to fix the last three components.
The best way to compare ordinations is to use Procrustes rotation.
Procrustes rotation uses uniform scaling (expansion or contraction) and
rotation to minimize the squared differences between two ordinations.
Package vegan has function procrustes to perform Procrustes analysis.
How much did we gain with using metaMDS instead of default isoMDS?
'

tmp <- wisconsin(sqrt(varespec))
dis <- vegdist(tmp)
vare.mds0 <- isoMDS(dis, trace = 0)
pro <- procrustes(vare.mds, vare.mds0)
pro
plot(pro)
#identify()
plot(pro, kind = 2) #plot residuals


# ----------------- PCA, and CA - eigenvector analyses ----------------------
'
Non-metric multidimensional scaling was a hard task, because any kind
of dissimilarity measure could be used and dissimilarities were nonlinearly
mapped into ordination. If we accept only certain types of dissimilarities
and make a linear mapping, the ordination becomes a simple task of
rotation and projection. In that case we can use eigenvector methods.
Principal components analysis (pca) and correspondence analysis (ca)
are the most important eigenvector methods in community ordination.
In addition, principal coordinates analysis a.k.a. metric scaling (mds) is
used occasionally. Pca is based on Euclidean distances, ca is based on
Chi-square distances, and principal coordinates can use any dissimilarities
(but with Euclidean distances it is equal to pca).
Pca is a standard statistical method, and can be performed with base
R functions prcomp or princomp. Correspondence analysis is not as ubiquitous,
but there are several alternative implementations for that also. In
this tutorial I show how to run these analyses with vegan functions rda
and cca which actually were designed for constrained analysis.
'

#                             Principal components analysis can be run as:
vare.pca <- rda(varespec)
vare.pca
plot(vare.pca)
biplot(vare.pca, scaling = -1)
#negative scaling to emphasize less abunant species
vare.pca <- rda(varespec, scale = TRUE)
vare.pca
plot(vare.pca, scaling = 3)
biplot(vare.pca, scaling = 3)
dim(varespec)

pcov.df <- as.data.frame(pcov.mat)
marsh.pca <- rda(pcov.df)
marsh.pca
plot(marsh.pca)
biplot(marsh.pca, scaling = -1)
#negative scaling to emphasize less abunant species
marsh.pca <- rda(pcov.df, scale = TRUE)
marsh.pca
plot(marsh.pca, scaling = 3)
biplot(marsh.pca, scaling = 3)



#                                      correspondance analysis
vare.ca <- cca(varespec)
vare.ca
plot(vare.ca)
chisq.test(varespec/sum(varespec))

plot(vare.ca, scaling = 1)
'
We already saw an example of scaling = 3 or symmetric scaling in pca.
The other two integers mean that either species are weighted averages of
sites (2) or sites are weighted averages of species (1). When we take
weighted averages, the range of averages shrinks from the original values.
The shrinkage factor is equal to the eigenvalue of ca, which has a
theoretical maximum of 1.
'
marsh.ca <- cca(pcov.df)
marsh.ca
plot(marsh.ca)
chisq.test(pcov.df/sum(pcov.df))









#                                      Fitting with env variables
'
Fitting environmental vectors is easy using function envfit. The
example uses the previous nmds result and environmental variables in
the data set varechem:

'

data(varechem)
ef <- envfit(vare.mds, varechem, permu = 999)
ef

ef.bin <- envfit(bin.bc.mds, pcov.env, permu = 999)
ef

'
The first two columns give direction cosines of the vectors, and r2 gives
the squared correlation coefficient. For plotting, the axes should be scaled
by the square root of r2. The plot function does this automatically, and
you can extract the scaled values with scores(ef, "vectors"). The
significances (Pr>r), or P-values are based on random permutations of
the data: if you often get as good or better R2 with randomly permuted
data, your values are insignificant.
You can add the fitted vectors to an ordination using plot command.
You can limit plotting to most significant variables with argument p.max.
As usual, more options can be found in the help pages.

'
plot(vare.mds, display = "sites")
plot(ef, p.max = 0.1)

plot(bin.bc.mds, display = "sites")
plot(ef.bin, p.max = 0.1)


#                                  Factors
'
Class centroids are a natural choice for factor variables, and R2
can be used as a goodness-of-fit statistic. The “significance” can be tested with
permutations just like in vector fitting. Variables can be defined as factors
in R, and they will be treated accordingly without any special tricks.
As an example, we shall inspect dune meadow data which has several
class variables. Function envfit also works with factors:

'

data(dune)
data(dune.env)
str(dune.env)
dune.ca <- cca(dune)

ef <- envfit(dune.ca, dune.env, permutations = 999)
ef

#marsh data
pcov.env <- pcov.env %>% mutate(month = factor(month),
                                dist_from_edge_m = factor(dist_from_edge_m))
str(pcov.env)

ef.bin <- envfit(bin.bc.mds, pcov.env, permu = 999)
plot(bin.bc.mds, display = "sites")
plot(ef.bin)









