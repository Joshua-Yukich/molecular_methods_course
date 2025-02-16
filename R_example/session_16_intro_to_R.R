# clear console, global environment  ----
rm(list = ls())
cat("\014")  #sends 'CTRL + L' to console optional


#-----------------------------------------------------------------------------#
#                                                                             #
#   Introduction to R                                                         #
#     Session 16 Intro to Statistical Computing                               #
#         J. Yukich October 17th 2019  jyukich@tulane.edu                     #
#                                                                             #
#   Description: Script is intended to show best practices for developing a   #
#   basic R script by showing good basic R programming conventions            #
#   as well as most standard datatypes handled in R.                          #
#                                                                             #
#   Dependencies: Script requires use of the tidyverse and devtools           #
#    packages but only for illustrative purposes. All other function calls    #
#    are contained in base R installs as of version 3.5.3                     #
#                                                                             #
#                                                                             #
#-----------------------------------------------------------------------------#



# loading required libraries ----
library(tidyverse)
library(devtools)

#set any other options ----

# print and get runtime ----
date()
current_date <- date()

# session info ----
devtools::session_info()
sessionInfo()

# above is an example of a sensible preamble default when working in R


# User defined functions for specific script ----



# Introduction to data types in R ----

# R has six basic data types, character, numeric, integer,
# logical, complex and raw. we will mainly use the first 4 of these,
# raw is used to store information as raw bytes and will not be discussed
# further. Some basics about the other types are below.

# character

x <- "a string which can be any abitrary length, like 100,000"

x
typeof(x)
attributes(x)
class(x)
str(x)

# numeric

y <- 1.5

y
typeof(y)
attributes(y)
class(y)
str(y)


# integer

z <- 2L

z
typeof(z)
attributes(z)
class(z)
str(z)

# logical

a <- TRUE

a
typeof(a)
attributes(a)
class(a)
str(a)


# complex
b <- 1+4i

b
typeof(b)
attributes(b)
class(b)
str(b)


# Introduction to Data structures in R----

# R has many, many data structures available but today we'll
# focus on the main basic types atomic vectors, lists, matrices,
# data.frames, and factors.

# atomic vectors ----

# these data structures are one-dimensional objects which contain
# an arbitrary number of elements all of which are the same data type.
# They can contain data of a number of types, but within a single
# vector all elements will be of the same data type.
# The c() function is the most simple way to make a vector.
# understanding vectors is one of the two fundamental building blocks
# to understanding R,the other is the function.

# a character vector

char_vec <- c("Sade", "Sarah Vaughn", "Selena", "Siouxsie Sioux")

char_vec
typeof(char_vec)
length(char_vec)
class(char_vec)
str(char_vec)


# a numeric vector
num_vec <- c(seq(from = 1, to = 10, by = 0.1))

num_vec
typeof(num_vec)
length(num_vec)
class(num_vec)
str(num_vec)

# an integer vector

int_vec <- c(1L, 4L, 100L)

int_vec
typeof(int_vec)
length(int_vec)
class(int_vec)
str(int_vec)

# a logical vector

log_vec <- c(TRUE, FALSE, FALSE, TRUE, TRUE)

log_vec
typeof(log_vec)
length(log_vec)
class(log_vec)
str(log_vec)


### vectors can have missing data


mis_vec <- c(1, 3.4, 2.5, NA, 7)

# can test if elements of a vector are missing with is.na() function.

is.na(mis_vec)

logical_na_vec <- is.na(mis_vec)

# produces a result identical in length to original
# vector but of logical type.

str(logical_na_vec)
length(logical_na_vec) == length(mis_vec)

sum(logical_na_vec)
# R interpretes true as 1 and FALSE as zero when used in numerical opeartions
# most R functions are designed to work on a vector.

# getting information out of a vector

# element wise using the '[' function.


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! R is [ROW, COL] !!!!!!!!!!!!!!!!

char_vec[3]
num_vec[12:27]
char_vec[c(1, 4, 3, 2)]
char_vec[1, 4]
# needs the c() function

#using a named vector

names(char_vec) <- c("top", "second", "middle", "best")

attributes(char_vec)
names(char_vec)


char_vec["best"]
char_vec[c("best", "second")]

# THE FACTOR VECTOR ----

# the factor vector is a special type of vector with
# attributes that is designed to hold a
# categorical or class variable data.
# It is restricted to containing a limited range of entries.
# these are stored with an integer index and factor labels
# for character descriptions. Similar approaches
# to labeled variables in STATA or class variables in SAS.

# well create a character vector with categorical
# data in it to show how this works.

new_char_vec <- c(rep("high", times = 3),
                  rep("low", times = 5),
                  rep("middle", times = 10))

new_char_vec
class(new_char_vec)

fac_vec <- factor(new_char_vec)

fac_vec
class(fac_vec)
attributes(fac_vec)
str(fac_vec)
levels(fac_vec)
nlevels(fac_vec)

as.numeric(fac_vec)

# note here that we get high as 1, low as 2, and middle as 3, why?
# This is because that was the order the levels were created in
# there are some obvious problems with this approach so never use it.

# one common challenge is that you want to add data to a factor
# that does not already have a level for it

# for example replace an obs in fac_vec with "lowest"

fac_vec[4] <- "lowest"

#error for invalid factor level.

# must first add the level to the factor.

levels(fac_vec) <- c(levels(fac_vec), "lowest")

fac_vec[4] <- "lowest"

fac_vec

# Matrices ----

# a matrix is just like a matrix in Math...
# in R it is just an atomic vector with a dimension
# dim() attribute. Like vectors they can contain only one datatype.
# They are two dimensional. # matrices with more than two dimensions
# are called arrays in R.

# making a matrix using the matrix function.


mat_num <- matrix(data = c(1:100), nrow = 10, ncol = 10)

mat_num
class(mat_num)
attributes(mat_num)
length(mat_num)
# note that this isn't rows or columns,
# but the same as the vector that filled it
str(mat_num)
typeof(mat_num)
dim(mat_num)

# character data might also be common.

mat_char <- matrix(state.name, nrow = 25, ncol = 2)

mat_char
class(mat_char)
attributes(mat_char)
length(mat_char)
str(mat_char)
typeof(mat_char)
dim(mat_char)

# getting things out of a matrix with the '[' function.


mat_char[1, 2]
mat_char[23, 1]
mat_char[4, 24]
# get errors when you ask for an element that does not exist.
# "subscript out of bounds"

# can use vector like indexing
mat_char[48]


# matrices can have row name attributes.

rownames(mat_char) <- paste0("row", c(1:25))
colnames(mat_char) <- paste0("col", c(1, 2))

attributes(mat_char)

# these can be used to access elements of a matrix as well.

mat_char["row1", "col2"]
mat_char["row23", "col1"]
mat_char["row4", "col24"]

mat_char["row3"]
# can't use vector subset this way

mat_char[1:10, 1:2]
# can subset to get smaller matrices.

mat_char[c(10, 5, 11, 13), c(2, 1)]
# in any rearrangement using vector inputs!

mat_char[, 2]
# all rows of column 2
mat_char[2, ]
# all columns for row 2

# other ways to make matrices, using rbind or cbind, or dim.

#row bind.

mat1 <- rbind(char_vec, char_vec)

mat1
# note how the names came along!
class(mat1)
attributes(mat1)
length(mat1)
str(mat1)
typeof(mat1)
dim(mat1)

# column bind

mat2 <- cbind(char_vec, char_vec)

mat2
# note how the names came along!
class(mat2)
attributes(mat2)
length(mat2)
str(mat2)
typeof(mat2)
dim(mat2)

#turn a numeric vector into a matrix with dim
num_vec <- 1:100

dim(num_vec) <- c(10, 10)
# what happened to num_vec in the environment pane # moved from values to data

dim(num_vec) <- c(10, 20)
# dims must match length of vector making up matrix

# Lists ----

# Lists are like vectors! one dimensional, but each element can be any valid R
# datatype and data structure.

#lets look at an example by using all the junk in our environment to make one.

test_list <- list(char_vec, num_vec, int_vec,
                  logical_na_vec, mis_vec, current_date,
                  b, mat_char)

# R studio gives you some handy viewing options in the environment pane.

typeof(test_list)
length(test_list)
# length of what comes out here?
# The list (does not consider the potential lengths)
# of the elements in the list!).

class(test_list)
str(test_list)
attributes(test_list)

#  why do we want this?  sometimes we want a number of
# different data structures and types
# packaged together for any number of purposes
# especially when it comes to applying functions to a list

# how do we get things back out of a list?

# main way is with function '[['

test_list[[2]]
test_list[[8]]

#can get information about the elements
attributes(test_list[[1]])
class(test_list[[8]])

#lists can have names

names(test_list) <- c("char_vec", "num_vec", "int_vec",
                      "logical_na_vec", "mis_vec",
                        "current_date", "b", "mat_char")

# therefore you can index a list with these.

test_list[["b"]]
typeof(test_list[["b"]])

# if you more than one index at a time returned as a list use the '[' function.

test_list[1:3]
class(test_list[1:3])

# works in any order and with names too (if the list is named).

test_list[c("b", "current_date")]


#getting sub-elements depends on what is inside the list itself.
# Look at the matrix.

test_list[["mat_char"]][23, 1:2]
# gives the 1st and 2nd columns of the 23rd row of the matrix stored
# in test_list as "mat_char"

test_list[[1]]["top"]
# gives the element named "top" of the vector stored in test_list place 1.

test_list[[4]][3]
# gives the third element of a vector stored in the 4 element of test_list


# DATA FRAMES ----

# data frames are the main data type that contains information
# for statistical analyses in R  these are more or less equivalent to a
# dataset in STATA and SAS. For lay purposes.
# They are 2 dimensional like matrices, N rows by M columns.
# Each column must be same length, each row, also must be the same length,
# each column must have the same datatype contained and each column, must
# have a name.

# by convention and standard in other applications,
# variables are stored in columns and 'observations'
# or records are in the rows.
# So it looks just like any dataset in any statistical software

# In R data frames are actually a special type of list,
# one which has to be rectangular. meaning that every element
# of the list is actually a vector of the same length.


#we are going to first construct a data.frame
# using the data.frame() function.


# make a few vectors to start with.

names <- factor(c(rep("Blanche", times = 5),
                  rep("Heavenly", times = 5),
                  rep("Maribelle", times = 5),
                  rep("Alma", times = 5)))

numbers <- floor(runif(n = 20, min = 1, max = 10))

author <- rep("Tennessee Williams", times = 20)

# vectors of the same length can easily be
# compiled into a dataframe using data.frame

data_df <- data.frame(names, numbers, author)

data_df
head(data_df)
tail(data_df)
dim(data_df)
length(data_df)
names(data_df)
nrow(data_df)
ncol(data_df)
sapply(data_df, class)
summary(data_df)

# data frames are also made by a number of importing functions.

# well load a data.frame that comes with base R for the next demonstrations
# the data() function can be deployed to load a preloaded R dataset.

data("USArrests")

head(USArrests)

# well assign this dataset to a name explicity
# so that you can see explicitly how this works

data_df <- USArrests

summary(data_df)

# getting information out of a data.frame

# they can be indexed/subset in very similar
# ways to matrices. use the '[' function.

data_df[2, ]
# second observation
data_df[, 2]
# second variable
data_df[3, 4]
# third observation of fourth variable
data_df[3, "Assault"]
# third observation of column named "assault"
data_df[, "Assault"]
# all observations of assault
str(data_df[, "Assault"])
# all observations of assault
# note this is a vector

# the '$' function can also be used to subset
# data.frames by specific variables.

data_df$Murder
# returns all observations of a variable as a vector

# other useful functions

class(data_df)
attributes(data_df)

# since this one has row names you can also subset with row names.

data_df["Mississippi", ]

is.list(data_df)
# data frames are lists

data_df[[2]]
# so you can also subset them like lists
data_df[["Murder"]]
data_df[["Assault"]][4]

# This is the end of the script.

