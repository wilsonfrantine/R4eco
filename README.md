# R4eco
A set of R scripts to help ecologists with no fancy tasks

At this very moment the package only includes functions to handle Renkonen's similarity index (S) calculation, its bootstrap implementation and some of conveniences for plotting a UPGMA with bootstrap values.

In this documentation you'll find a very simple example of how to use it

## Installation

to install this package you will must have the `devtools` package in your machine. To get that run the code bellow


```{r
#First install devtools
install.packages("devtools")

#then load it
library(devtools)

#to download R4eco
install_github("wilsonfrantine/R4eco")

#now you can load it
library(R4bio)
```

## Using some functions

We must have some data to start with. If you have any, you may generates it from the package with `simulate_data()`.

```{r}
#Simple like this:
simulate_data()
```
To calculate renkonen similarity you can just:

```{r}
df <- simulate_data(5, 9)

S <- renkonen(df)
head(S)
```
With bootstrap:
```{r}
df<-simulate_data()
S<-renkonen(df, boot=1000)
head(S)
```

To plot it as a dendrogram:

```{r}
df<-simulate_data()
S<-renkonen(df, boot=1000)
trees<-upgma(S)
plot_upgma(trees)
```

You can also compute the consensus dendrogram to save as newick or use in another software.

```{r}
df <- simulate_data()
S <- renkonen(df, boot=1000)
trees <- upgma(S)
tree_consensus <- tree_consensus(trees, type="mcc")
```

## Contatct

If you got any crash or need some help, please mail-me:
Wilson Frantine-Silva
wilsonfrantine@gmail.com
