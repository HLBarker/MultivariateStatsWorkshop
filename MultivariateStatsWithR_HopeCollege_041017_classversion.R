# -----------------------------------------------------------------------------
# Hope College Multivariate Stats + R Workshop ################################  
# -----------------------------------------------------------------------------
  # Author = Hilary Barker, 4/10/17
  # Script adapted from "Numerical Ecology with R"
    # Borcard, Gillet, Legendre 2011

# Purpose of this script: 
  # An introduction to R
  # Prepare community data (transform, relativize)
  # Analyze community data (NMDS, dbRDA, PerMANOVA, indicator species analysis)


# -----------------------------------------------------------------------------
# Load packages ###############################################################
# -----------------------------------------------------------------------------
library(vegan) # community statistics package
library(labdsv) # calculates indicator values for species in different groups
  # (i.e., treatments)


# -----------------------------------------------------------------------------
# Load Doubs fish data ########################################################  
# -----------------------------------------------------------------------------
spe <- read.csv('http://www.davidzeleny.net/anadat-r/data-download/DoubsSpe.csv', row.names = 1)

  # contains the abundance of 27 fish species surveyed from 30 sites along the 
    # Doubs river

env <- read.csv('http://www.davidzeleny.net/anadat-r/data-download/DoubsEnv.csv', row.names = 1)

  # contains 11 environmental variables collected from the 30 Doubs river sites
    # distance from source, altitude, slope, mean minimum discharge, water pH,
    # water hardness, [phosphate], [nitrate], [ammonium], dissolved O2, biological 
    # oxygen demand

spa <- read.csv('http://www.davidzeleny.net/anadat-r/data-download/DoubsSpa.csv', row.names = 1)

  # contains the geographical (Cartesian, X and Y) coordinates of the 30 Doubs 
    # river sites


# -----------------------------------------------------------------------------
# Explore site locations ######################################################  
# -----------------------------------------------------------------------------
  # We will create a 'map' of the site locations along the Doubs river

# create an empty plot with proportional axes 1:1 and titles
plot(spa, asp = 1, type = "n", main = "Site Locations",
     xlab = "X coordinate (km)", ylab = "Y coordinate (km)")

# add a blue line that connects the sites (Doubs river)


# add the site numbers along the river



# -----------------------------------------------------------------------------
# Explore fish community data ################################################# 
# -----------------------------------------------------------------------------
  # Community data is zero-inflated! 
# number of 0s in the dataset 
# proportion of 0s


  # Many species are rare

  # count the # of times a species is absent from a site
  # it is common in community analyses to delete any species that occur in less
  # than 5% of your samples (this reduces noise). Since these fish species all 
  # occur in more than 5% of the river sites, we will keep all of the fish data 
  # for further analyses


# try out various kinds of data transformations, standardizations
  # compare distributions of fish species with boxplots and select
  # a transformation that best suites the data and your underlying
  # question

spe.scal <- decostand(spe, "max") # scale abundances by dividing them by the 
  # maximum values for each species
spe.relsp <- decostand(spe, "total", MARGIN = 2) # scale abundances by dividing
  # them by the species totals (relative abundance per species)
spe.rel <- decostand(spe, "total", MARGIN = 1) # scale abundances by dividing 
  # them by the site totals (relative abundance, or relative frequencies,
  # per site)
spe.norm <- decostand(spe, "normalize") # scale abundances to mean of 0 and SD of 1
spe.hel <- decostand(spe, "hellinger") # compute relative frequencies by rows 
  # (site profiles), then square root them
spe.chi <- decostand(spe, "chi.square") # Chi-square transformation (by both species
  # and sites)
spe.wi <- wisconsin(spe) # abundances are first ranged by species maxima and
  # then by site totals

# start plotting! (we'll just look at the "LOC" species for now)
par(mfrow = c(2,2)) # change the figure output to display 2x2 figures

boxplot(spe$LOC, sqrt(spe$LOC), log1p(spe$LOC),
        las = 1, main = "Simple transformation", 
        names = c("raw", "sqrt", "log"), col = "bisque")
boxplot(spe.scal$LOC, spe.relsp$LOC,
        las = 1, main = "Standardization by species", 
        names = c("max", "total"), col = "lightgreen")
boxplot(spe.hel$LOC, spe.rel$LOC, spe.norm$LOC,
        las = 1, main = "Standardization by sites", 
        names = c("Hellinger", "total", "norm"), col = "lightblue")
boxplot(spe.chi$LOC, spe.wi$LOC,
        las = 1, main = "Double standardization", 
        names = c("Chi-square", "WI"), col = "orange")

# How might these different transformations/relativizations change your statistical 
  # results?


# We'll continue the analyses with sqrt() transformed fish community data


# -----------------------------------------------------------------------------
# Explore environmental data ##################################################  
# -----------------------------------------------------------------------------
  # investigate relationships/correlations among environmental variables
  # load this "panelutils.R" script that makes effective correlation plots
    # panelutils.R is from "Numerical Ecology with R" and can be downloaded from
    # https://github.com/JoeyBernhardt/NumericalEcology/blob/master/panelutils.R

source("~/path_to_/panelutils.R") # change the "~/path_to_" to the directory path
  # that leads to the panelutils.R script saved on your computer

op <- par(mfrow = c(1,1), pty = "s") # save the default figure output formatting

  # Since the environmental data are all in various units (which is like comparing
    # apples to oranges), we need to first standardize the data (mean = 0, SD = 1)


pairs(env.norm, lower.panel = panel.smooth, upper.panel = panel.cor, 
      method = "pearson", diag.panel = panel.hist, 
      main = "Bivariate plots with histograms and smooth curves")

 # put the figure output format back to the default conditions


# -----------------------------------------------------------------------------
# Do the fish communities vary with environmental conditions? #################
# -----------------------------------------------------------------------------
  # We'll use distance based Redundancy Analysis (dbRDA) to identify which
  # environmental variables are correlated with fish community composition.
  # Then, we'll visualize this with nonmetric multidimensional scaling (NMDS)
  # and overlaid environmental vectors


# make a bray-curtis dissimilarity
  # matrix for community data. Bray-curtis dissimilarity is well suited to community
  # data since it takes into account both differences in species richness and 
  # abundance between communities. Dissimilarity values vary from 0 to 1. 0 indicates
  # that two communities share all of the same species with the same  
  # abundances and 1 indicates that two communities are completely different

# we have a study site where no fish were found!

# since we can't compute distances with sites that contain no community data 
# (i.e., fish), we must delete the empty site (row 8) from all dataframes 
# for analyses


# we'll analyze these data with distance-based RDA 
  # dbRDA transforms the many-dimensioned distance matrix into a series of 
  # coordinates using Principal Coordinates Analysis (PCoA). This step reduces the
  # number of dimensions in the y-variable while still accounting for the variation
  # in the dissimiliarity matrix (similar to a PCA), and then redundancy analysis
  # summarizes linear relationships between the PCoA coordinates and the x variables.
  # Typical RDA can only use euclidean distance matrices, but distance-based RDA
  # can use many different distance metrics, such as Bray-Curtis, and thus why it
  # is more often used with community data.

# function "capscale" in vegan performs distance-based RDA and function "ordiR2step"
# conducts stepwise model selection in vegan (for RDA models)
  # thus, we'll first create two dbRDA models, a null model and model with all of the
  # environmental variables included, and then we'll use model selection to determine
  # which of the environmental variables are significant in explaining variation in 
  # fish communities and these variables will be the ones that are included in the 
  # final dbRDA model (we'll call this rda.select, below)

 # null model with intercept only
rda1 <- capscale(spe.dist ~ das + alt + pen + deb + pH + dur + pho + nit + amm + oxy + dbo, 
         data = env.norm, dist = "bray") # full model with all environmental variables



# we'll then visualize these data with NMDS and overlaid environmental vectors 
  # NMDS is similar to dbRDA, in that it reduces the dimemnsions in the y-variable 
  # (the dissimilarity matrix) into something more manageable, (e.g., just 2 axes)
  # while at the same time, keeping as much of the original variation as possible.
  # NMDS is non-metric and thus is well suited for community data that is non-normal
  # and rarely linear in nature. 

# the "metaMDS" function in vegan performs NMDS

  # the NMDS stress is a measure of how well the data are modeled with this method
  # the lower the stress, the better the fit
     # typically we want our stress levels to be below 15-20% 

 # these are the NMDS points that we'll use for graphing

# the "envfit" function in vegan identifies the correlation between various
  # environmental variables with the community NMDS, here we will just look at the 
  # environmental variables that were selected using dbRDA



# make an empty NMDS plot

# Add the NMDS points for each site
  # and label them as the site number (1-30)

# overlay the environmental vectors


# In any ordination plot (for instance this NMDS plot), points that are close
  # together are more similar. Thus, sites 22, 28, and 27 all have very similar 
  # fish communities. Points that are far away from each other are more different.
  # Thus sites 1 and 23 have very different fish communities. 

# What do you think the environmental vectors show in this plot?


# -----------------------------------------------------------------------------
# Do the fish communities vary with geographic distance? ######################
# -----------------------------------------------------------------------------
  # We'll compare the distance among fish communities with the geographic
  # distance among site locations using Procrustes correlation.

# Procrustes correlation requires two distance matrices, in this case our bray-curtis
  # fish community distance matrix (our y-variable), and a geophraphic distance matrix 
  # (our x variable). It then rotates and scales the distance matrices to
  # maximize the similarity between the two matrices (minimizing the sum of squared
  # differences) and then calculates the correlation between the two matrices

# We need to make the geographic distance matrix, using euclidean distance


# function "protest" in vegan performs this correlation, tests for the significance
  # of the correlation, and gives the r correlation coefficient value

 # Procrustes p-value
 # Procrustes r value

# Plot fish abundance on the river for different species on the site 'map'
par(mfrow = c(2,2)) # divide up plot into 4 quadrants

plot(spa[,1:2], asp = 1, col = "brown", cex = spe$TRU, main = "Brown trout",
     xlab = "X coordinate (km)", ylab = "Y coordinate (km)")
lines(spa, col = "light blue")

plot(spa[,1:2], asp = 1, col = "gray", cex = spe$OMB, main = "Grayling",
     xlab = "X coordinate (km)", ylab = "Y coordinate (km)")
lines(spa, col = "light blue")

plot(spa[,1:2], asp = 1, col = "blue", cex = spe$BAR, main = "Barbel",
     xlab = "X coordinate (km)", ylab = "Y coordinate (km)")
lines(spa, col = "light blue")

plot(spa[,1:2], asp = 1, col = "green", cex = spe$BCO, main = "Common bream",
     xlab = "X coordinate (km)", ylab = "Y coordinate (km)")
lines(spa, col = "light blue")


# -----------------------------------------------------------------------------
# Do the fish communities vary by river zone? #################################
# -----------------------------------------------------------------------------
  # We'll use PerMANOVA to test whether fish communities vary by upstream vs.
  # downstream groups and nested zone groups. Then, we'll use indicator species
  # analysis to identify which fish species explain the differences among these 
  # groups.


# add grouping variables that split up the river into 4 different zones

# Up vs. downstream variable


# Within the "updown" variable, we'll split up the upstream and downstream
  # into 4 different zones


# This is a "nested design" = one variable (zone) is nested within another
  # variable (up vs. downstream)


# adonis in vegan performs perMANOVA analyses
# like ANOVA/MANOVA, perMANOVA partitions variation in the y-variable according
  # to various x-variables or factors. Yet ANOVA/MANOVA assume normally-distributed 
  # data + euclidean (linear) distances in the y-variable, whereas perMANOVA works 
  # with any distance metric and uses permutations and thus any kind of data 
  # distribution is allowed with this model (i.e., data do not need to be normal)
    # use "A/B" notation to nest variable "A" within variable "B"


# Let's look back out our NMDS plot and see whether upstream vs. downstream sites
  # cluster together
# make an empty NMDS plot

 # upstream sites will be light blue

 # overlay the environmental vectors

# What does the upstream vs. downstream divide show in the NMDS plot?


# Indicator species analysis: 

# which fish species are indicative
  # of the upstream part of the river compared to the downstream section 
  # (and vice versa)?
    # look at the indval information in the help window in Rstudio to better 
    # understand the output (relfrq, relabu, indval, etc.)

