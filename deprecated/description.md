# Part 1
####Working to <i>standardize flanks</i> in order to achieve better accuracy

Results:
* <b>Trends are more easily noticeable</b> when points are directly inspected
* Linear regression trend lines are <b>less consistent</b>
* Flank count sums were found to be approximately normally distributed; thus, their means
were used in the above-mentioned graphs
* Standard deviation lines around the trend lines were computed by taking evenly sized
groups of adjacent points and computing their standard deviation; then, the standard
deviation lines were further smoothed/linearized through linear regression

These tests were run under the following conditions:
* <i>Only copy numbers at/below 8 were allowed</i>
* <i>Array lengths above 5000 were discarded</i>
* <i>Indistinguishable sequences were discarded</i>
* <i>Rows with flank/inside-array counts of 0 were discarded</i>
* <i>Flanks were held constant depending on read length</i>

```
{read length of 100 -> flank sum of 102}
{read length of 148 -> flank sum of 103}
{read length of 250 -> flank sum of 194}
```

An approach was made to discard the lowest percentiles of pattern sizes; however, this
did not yield any useful results as it only gave more power to points of high influence.

<b>Finding note 1:</b> the mean copy number specific to the dataset, when excluding all
copy numbers of 8, was found to be ```3.8```. Thus, theoretical gain/loss trend lines
were computed by scaling the 0/0 homozygous genotype trend line.

The factors used in scaling were computed as the copy number, ```C```, in the equation
* ```(C + G) / C``` where ```G``` is a gain/loss value between -1 and 1, inclusive

# Part 2
####Working to <i>classify datasets</i> and <i>individual points</i> within those datasets

<i>In the previous meeting, it was discovered that as array length increased past 1000, in
order to detect statistical differences in array length between two different arrays, the
difference had to be greater than 25.3%. What is described below is a different approach.</i>

Basic summary: the overall objective was to be able to determine whether a dataset itself,
in addition to individual points from the dataset, could be classified as gains/losses
relative to the 0/0 homozygous simulated dataset.

These tests were run under the following conditions:
* <i>Only copy numbers at/below 8 were allowed</i>
* <i>Array lengths above 5000 were discarded</i>
* <i>Indistinguishable sequences were discarded</i>
* <i>Rows with flank/inside-array counts of 0 were discarded</i>

Results -- emulating a "nearest neighbors" approach:
* The most successful algorithm for dataset <b>genotype</b> classification functioned by
performing a residual-sum-of-squares comparison of the input dataset to the simulated
data/regression lines; another alternative to this algorithm featured a difference
calculation between the parameters of the standard linear regression equations for each
dataset (this yielded more varied results that tended to be more accurate at face value) 
* The most successful algorithm for dataset <b>read length</b> classification functioned by
performing a residual-sum-of-squares comparison of the input dataset to the simulated
data/regression lines
* An algorithm was developed for classification of individual points as belonging to a
given genotype; it was found that the <i>correct functionality</i> of this algorithm was
largely impacted by the pattern size/array length of the data point in question
(specifically, with higher ratios of array-length-to-pattern-length, the classification
was more accurate; with lower ratios, it was more prone to errors)
* Each genotypes line was calculated using the method described in <b>finding note 1.</b>

<i>The algorithms were only tested on HG001, HG002, and HG007. More data would be needed
in order to validate the algorithm's functionality.</i>

# Part 3: misc. notes
####Predicting array length based on <i>inside hit/flank ratio</i> and <i>pattern size</i>

Results:
* In order to predict array length of a given data point, a quadratic curve was fit to
the data; it was able to model the long run trend relatively accurately
* The short run trend, however, was better modeled by a linear model; thus, it was
concluded that using a piece-wise function for predicting array length would be optimal 

Github repository: https://github.com/mkorovkin2/tandem-rep-mkorovkin